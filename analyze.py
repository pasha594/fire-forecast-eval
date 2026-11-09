#!/usr/bin/env python3
"""Growth-only precision/recall of archived ELMFIRE forecasts vs Cornea perimeters.

Reads the archive in raw/ (populate with `collect.py --local`, or sync a bucket
with --sync-url) and writes data/metrics.csv. One row per (model, fire, forecast run,
percentile, perimeter): each run is scored against perimeters at each perimeter's own
time (its exact lead, no tolerance), as scored_perimeters() picks them, and every row is
kept (no growth, offset or off-grid filters). score_grid() is the scorer every model
uses (the FARSITE hindcast calls it too).

  analyze.py                          compute metrics
  analyze.py --sync-url https://pub-xxx.r2.dev   sync bucket -> raw/ first
  analyze.py --refresh-perims         re-fetch Cornea perimeters before computing
  analyze.py --inspect PATH.tif       print raster calibration info
  analyze.py --overlay SLUG RUN PCT H write an eyeball GeoJSON to raw/
  analyze.py --skill-only             rebuild + upload per-fire skill files from the
                                      published metrics.csv (needs S3_* env)
"""
import argparse
import concurrent.futures
import csv
import json
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

# --- calibration constants (inspected ca-bug 20260813_092300/50.tif, 2026-08-13) ---
# float32, EPSG:326xx (per-fire UTM), 30m pixels, nodata=0.0.
# Values are HOURS since forecast start (max ~279 ≈ 11.6d on a ≤14d horizon).
# Already-burned-at-start area is NOT encoded (no 0/negative values; it is nodata),
# so the predicted mask is simply valid & (toa <= H). Growth is what burns outside
# what had burned when the run began: the Cornea baseline perimeter AND Pyrecast's own
# start (see START_VAR), which it seeds with satellite detections on top of the
# perimeter (in 41 of 42 ca-plaskett runs run_ts is a VIIRS overpass to the minute).
# run_ts is UTC (tif Last-Modified ≈ run_ts + ~3h processing lag).
TOA_TO_HOURS = 1.0

# Each run is scored against perimeters over its first SCORE_DAYS days (a day = 24 h from
# the run's start), each at its own time: every perimeter with a mapped time; on a day with
# none, the latest perimeter published that day (timed by its last edit, not a re-publish of
# an earlier shape); and, so no analysis shown before Nov 2026 is lost, the perimeter the
# previous method scored at each of its fixed lead times (OLD_HORIZONS_H: the one nearest
# run + H within ACTUAL_TOL_H, shown when within OLD_MAX_OFFSET_H), rescored at its own time.
SCORE_DAYS = 7
OLD_HORIZONS_H = (12, 24, 48, 72, 96, 120, 168, 240, 336)
OLD_MAX_OFFSET_H = 3
MODEL = "pyrecast"
PUBLIC_BUCKET_URL = "https://f005.backblazeb2.com/file/fire-forecast-archive"
# baseline = latest perimeter mapped <= run + this. run_ts is cut to the minute, and is
# the time of whatever the run started from: sometimes a perimeter, usually a satellite pass
BASELINE_TOL_H = 59.999 / 3600
# Pyrecast's start: the first granule of its hours-since-burned output (named for run_ts)
# marks the cells it had burned (400) or had burning (401) when the forecast began. Its
# time of arrival is never set there, so growth into those cells can't be predicted:
# scored as growth, they'd all be misses. The start is the same at every percentile.
# Without that granule (the tar wasn't archived, or the hourly purge took the first
# hour first), the satellite detections the run could have seen stand in for it: every
# detection footprint up to run_ts. Each row records which it used (start_basis).
START_VAR = "hours-since-burned"
START_CODES = (400, 401)
START_FILE = "start.tif"           # raw/forecast_archive/<slug>/<run_ts>/, the granule as published
START_NONE = "start.none.json"     # beside it when no archived tar holds that granule: {"checked": [pct, ...]}
START_PCT_ORDER = ("50", "30", "70", "10", "90")
ACTUAL_TOL_H = 12.0      # the previous method: perimeter mapped (else edited/published) within +/- this of run + H
# Perimeter times (see PerimeterSet): a snapshot's file key / Cornea date is its last
# EDIT (poly_DateCurrent), a median ~17h after the polygon was mapped. The mapped time is
# poly_PolygonDateTime, used only if it is plausible and the geometry still matches the
# first edit that carried it (editors often redraw a perimeter without updating it);
# otherwise the edit time stands in, and each row records which one it used.
CAPTURE_REL_TOL = 0.02   # later edits within this area change of the first edit,
CAPTURE_TOL_AC = 25.0    # and within this many acres, are corrections of the same capture
CAPTURE_LATE_S = 300     # a capture time may sit at most this far after its own edit
SQM_PER_ACRE = 4046.8564224
ISO_Z = "%Y-%m-%dT%H:%M:%SZ"
PERCENTILES = ("10", "30", "50", "70", "90")

REPO = Path(__file__).resolve().parent
RAW = REPO / "raw"
SKILL_DIR = REPO / "data" / "skill"
START_DIR = SKILL_DIR / "start"    # <slug>/<run_ts>.tif, published as skill/start/...

SKILL_UPLOAD_WORKERS = 16
SYNC_WORKERS = 16          # parallel downloads in sync_bucket
SKILL_INDEX_EVERY_S = 300   # --publish: the merged skill index goes up at most this often mid-run

# lead_h: hours from the run's start to the perimeter's time (mapped, else last edited);
# day: 1-7 (24 h periods from the run's start; later only for "previous" picks);
# chosen: why that perimeter is scored (see scored_perimeters): mapped | daily | previous
CSV_COLUMNS = ["model", "slug", "cornea_id", "run_ts", "run_dt_utc", "percentile",
               "lead_h", "day", "chosen", "baseline_perim_ms", "baseline_time_utc",
               "baseline_time_basis", "baseline_offset_h", "actual_perim_ms", "actual_time_utc",
               "actual_time_basis", "pixel_m", "baseline_acres",
               "start_basis", "start_acres", "burned_before_acres", "pred_total_acres", "pred_new_acres", "act_new_acres", "inter_acres",
               "precision", "recall", "iou", "act_outside_acres", "act_offgrid_frac",
               "hs_n", "hs_new_acres", "hs_hit_acres", "hs_capture",
               "hs_artifact_n", "hs_transient_n", "expired_run", "notes"]

# satellite active-fire footprint radii (m). A detection is a fire somewhere in its pixel,
# often far smaller than the pixel, so a footprint is only as big as it needs to be:
#  VIIRS    375 m px: half the pixel, the least at which neighbouring detections (375 m
#           apart at nadir) still touch; off-nadir pixels are bigger, and theirs don't
#  MODIS    1 km px: 250 m, half the pixel's half-width. Its detections are ~1 km apart and
#           don't touch either way; the full 500 m claimed ~194 acres burned per detection
#  Landsat  30 m px: half the pixel (it used to get VIIRS's 187.5 m)
# Every detection also takes the grid cell it falls in, which a disk smaller than a cell
# can miss (HotspotSet.mask).
HS_RADIUS_VIIRS = 187.5
HS_RADIUS_MODIS = 250.0
HS_RADIUS_LANDSAT = 15.0
# Detections that don't look like the fire are left out of the hotspot metric instead
# of counting as misses against the forecast (see HotspotSet._screen):
#  artifact:  farther than HS_ARTIFACT_M from the perimeter current when it was detected,
#             and its footprint never touches a perimeter mapped after it. Without a
#             perimeter before it, or until one is mapped after it, it can't be judged
#             and is kept.
#  transient: no other detection within HS_PERSIST_M on any satellite's nearest pass
#             before or after its own (its own satellite's included), so a source that
#             revisits rarely, like Landsat every 8-16 days, is confirmed by the daily
#             ones. A radius, not the same spot: the detected location shifts between
#             passes with the view angle and as the fire moves.
HS_ARTIFACT_M = 1000.0
HS_PERSIST_M = 1000.0
HS_PASS_GAP_S = 20 * 60   # one source's detections closer in time than this are one overpass


def log(msg):
    print(f"[analyze] {msg}", flush=True)


def parse_iso(s):
    s = s.rstrip("Z")
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
        try:
            return datetime.strptime(s, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"bad timestamp: {s}")


def load_json(path):
    p = Path(path)
    return json.loads(p.read_text()) if p.exists() else None


# ---------------------------------------------------------------- sync

# B2's public URL throws occasional 500/503s and the odd dropped connection; a
# single one used to kill the whole analyze run mid-sync (Oct 29 + 30)
RETRY_STATUS = {429, 500, 502, 503, 504}
RETRY_DELAYS_S = (2, 5, 15, 30, 60)


def http_get(sess, url, timeout=(10, 120), headers=None):
    """GET with backoff on transient errors (5xx, 429, connection/timeout). A 404
    or other 4xx raises at once: a manifest-listed object should exist."""
    import requests
    for delay in RETRY_DELAYS_S + (None,):
        try:
            r = sess.get(url, timeout=timeout, headers=headers)
            if r.status_code not in RETRY_STATUS or delay is None:
                r.raise_for_status()
                # urllib3 1.x returns a cut-off body silently; sync would cache it
                want = r.headers.get("Content-Length")
                if want and not r.headers.get("Content-Encoding") and len(r.content) != int(want):
                    raise requests.exceptions.ChunkedEncodingError(
                        f"truncated body: {len(r.content)}/{want} bytes")
                return r
            err = f"HTTP {r.status_code}"
        except (requests.ConnectionError, requests.Timeout,
                requests.exceptions.ChunkedEncodingError) as e:
            if delay is None:
                raise
            err = type(e).__name__
        log(f"retrying in {delay}s ({err}): {url}")
        time.sleep(delay)


# A fire's perimeter index.json and latest hotspot days are MUTABLE while collect
# still writes them, and a cached copy then hides newer data (the rows-frozen-at-463
# bug). But collect only writes fires that are listed upstream or had a run in the
# last 15 days (collect.main's scoring set); after that their files never change.
# Re-fetching every finished fire anyway cost ~2,500 downloads (~13 min) per run by
# Nov. A finished fire is re-fetched once more after it freezes, then left alone
# until its manifest entries change. 3 days of slack covers clock skew and a collect
# run overlapping this sync.
FROZEN_AFTER = timedelta(days=15 + 3)
# collect rewrites hotspot days since (now - 3d) by date: up to 4 calendar days
HOTSPOT_REWRITE_DAYS = 4
# in raw/, so it travels in the actions cache together with the files it vouches for
FROZEN_STATE = "sync_frozen.json"


def frozen_fingerprints(manifest, now):
    """{slug: fingerprint} for the fires collect can no longer write. A fire with no
    datable run counts as active, so it is always re-fetched."""
    last = {}
    for e in manifest["runs"].values():
        try:
            t = datetime.strptime(e["run_ts"], "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
            if e.get("last_seen"):
                t = max(t, parse_iso(e["last_seen"]))
        except (KeyError, ValueError):
            continue
        last[e["slug"]] = max(last.get(e["slug"], t), t)
    perims = manifest.get("perimeters") or {}
    hotspots = manifest.get("hotspots") or {}
    out = {}
    for slug in set(perims) | set(hotspots):
        if slug not in last or now - last[slug] < FROZEN_AFTER:
            continue
        ep, days = perims.get(slug) or [], hotspots.get(slug) or []
        # any manifest change for the fire (it came back, a late backfill) un-freezes it
        out[slug] = {"last_activity": last[slug].strftime(ISO_Z),
                     "perimeters": [len(ep), ep[-1] if ep else None],
                     "hotspots": [len(days), days[-1] if days else None]}
    return out


def tar_first_member(get_range, key):
    """(name, bytes) of a tar's first member by two ranged reads (its 512-byte header,
    then the member), not the whole tar"""
    head = get_range(key, 0, 511)
    name = head[:100].rstrip(b"\0").decode("utf-8", "replace")
    size = int(head[124:136].rstrip(b"\0 ") or b"0", 8)
    return name, (get_range(key, 512, 511 + size) if size > 0 else b"")


def fetch_start_frame(get_range, run_key, entry):
    """Pyrecast's run-start granule (see START_VAR) -> raw/forecast_archive/<run_key>/start.tif.
    collect tars the granules sorted by name, i.e. by time, so it's the first member of a
    {pct}_hours-since-burned.tar, if that tar still caught it. Each archived tar is looked
    at once: when none has it, START_NONE lists the ones checked, and only a tar archived
    later is tried. -> "new", "none" or None (nothing to look at)"""
    d = RAW / "forecast_archive" / run_key
    if (d / START_FILE).exists():
        return None
    vs = entry.get("vars") or {}
    pcts = [p for p in START_PCT_ORDER if ((vs.get(p) or {}).get(START_VAR) or {}).get("ok")]
    checked = (load_json(d / START_NONE) or {}).get("checked", [])
    todo = [p for p in pcts if p not in checked]
    if not todo:
        return None
    want = f"{START_VAR}_{entry['run_ts']}.tif"
    d.mkdir(parents=True, exist_ok=True)
    for p in todo:
        name, body = tar_first_member(get_range, f"forecast_archive/{run_key}/{p}_{START_VAR}.tar")
        if name == want and body[:4] in (b"II*\x00", b"MM\x00*"):
            tmp = d / (START_FILE + ".part")
            tmp.write_bytes(body)
            tmp.replace(d / START_FILE)
            (d / START_NONE).unlink(missing_ok=True)
            return "new"
        checked.append(p)
    (d / START_NONE).write_text(json.dumps({"checked": checked}))
    return "none"


def sync_start_frames(url, manifest, only_slug=None):
    """fetch_start_frame for every run with a time-of-arrival tif. A failure only logs:
    the run is scored with the stand-in start (start_basis says so) and retried next sync."""
    import threading
    import requests
    local = threading.local()

    def get_range(key, a, b):
        if not hasattr(local, "sess"):
            local.sess = requests.Session()
        r = http_get(local.sess, f"{url}/{key}", headers={"Range": f"bytes={a}-{b}"})
        return r.content if r.status_code == 206 else r.content[a:b + 1]   # 200: range ignored

    def one(item):
        run_key, entry = item
        try:
            return fetch_start_frame(get_range, run_key, entry)
        except Exception as e:   # noqa: BLE001 -- one run's start must not stop the sync
            log(f"start frame for {run_key} failed, retried next sync: {e}")
            return "failed"

    runs = [(k, e) for k, e in manifest["runs"].items()
            if (not only_slug or e["slug"] == only_slug)
            and any(i.get("ok") for i in e["files"].values())]
    with concurrent.futures.ThreadPoolExecutor(8) as ex:
        got = list(ex.map(one, runs))
    log(f"start frames: {got.count('new')} new, {got.count('none')} not archived, "
        f"{got.count('failed')} failed; {sum((RAW / 'forecast_archive' / k / START_FILE).exists() for k, _ in runs)}"
        f" of {len(runs)} runs have one")


def sync_bucket(url, only_slug=None):
    """Mirror the R2 bucket (public r2.dev URL) into raw/, manifest-driven; with
    only_slug, just that fire's files (a quick local copy to test with)."""
    import requests
    url = url.rstrip("/")
    sess = requests.Session()

    def get(key):
        return http_get(sess, f"{url}/{key}")

    n = 0
    for key in ("manifest.json", "fire_matches.json"):
        p = RAW / key
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(get(key).content)
    manifest = load_json(RAW / "manifest.json")
    frozen = frozen_fingerprints(manifest, datetime.now(timezone.utc))
    state_path = RAW / FROZEN_STATE
    settled = load_json(state_path) or {}
    mine = (lambda slug: slug == only_slug) if only_slug else (lambda slug: True)

    def is_settled(slug):
        """local copies were fetched after collect stopped writing this fire"""
        return slug in frozen and settled.get(slug) == frozen[slug]

    wanted = []
    for run_key, entry in manifest["runs"].items():
        if not mine(entry["slug"]):
            continue
        for pct, info in entry["files"].items():
            if info.get("ok"):
                wanted.append(f"forecast_archive/{run_key}/{pct}.tif")
    refreshing = set()
    for slug, epochs in (manifest.get("perimeters") or {}).items():
        if not mine(slug):
            continue
        # the perimeter timeline grows: refresh unless the fire is settled
        idx = RAW / "perimeter_archive" / slug / "index.json"
        if not is_settled(slug):
            refreshing.add(slug)
            if idx.exists():
                idx.unlink()
        wanted.append(f"perimeter_archive/{slug}/index.json")
        wanted += [f"perimeter_archive/{slug}/{e}.geojson" for e in epochs]
    for slug, days in (manifest.get("hotspots") or {}).items():
        if not mine(slug):
            continue
        # collect rewrites the most recent days; refresh them unless the fire is settled
        recent = set() if is_settled(slug) else set(days[-HOTSPOT_REWRITE_DAYS:])
        if recent:
            refreshing.add(slug)
        for d in days:
            p = RAW / "hotspot_archive" / slug / f"{d}.geojson"
            if d in recent and p.exists():
                p.unlink()
            wanted.append(f"hotspot_archive/{slug}/{d}.geojson")
    # in parallel: one at a time, a cold start (no Actions cache: ~22k objects in Nov) took
    # over 2 h and, with the scoring, ran past the job's time limit
    import threading
    missing = [k for k in wanted if not (RAW / k).exists()]
    local = threading.local()      # one HTTP session per worker thread

    def fetch(key):
        if not hasattr(local, "sess"):
            local.sess = requests.Session()
        p = RAW / key
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(p.suffix + ".part")
        tmp.write_bytes(http_get(local.sess, f"{url}/{key}").content)
        tmp.replace(p)

    with concurrent.futures.ThreadPoolExecutor(SYNC_WORKERS) as ex:
        for _ in ex.map(fetch, missing):
            n += 1
    sync_start_frames(url, manifest, only_slug)
    if only_slug:
        log(f"sync ({only_slug} only): {n} new objects ({len(wanted)} tracked)")
        return
    # only once every download landed: each frozen fire is now settled, either from
    # an earlier run or because this run just fetched its final copies
    tmp = state_path.with_suffix(".json.part")
    tmp.write_text(json.dumps(frozen, sort_keys=True))
    tmp.replace(state_path)
    log(f"sync: {n} new objects ({len(wanted)} tracked); refreshed {len(refreshing)} fires, "
        f"{len(frozen) - len(refreshing & set(frozen))} finished fires left as cached")


# ---------------------------------------------------------------- inspect

def inspect(path):
    import rasterio
    with rasterio.open(path) as ds:
        a = ds.read(1)
        nod = ds.nodata
        print(f"crs {ds.crs} | res {ds.res} | dtype {a.dtype} | nodata {nod} | shape {a.shape}")
        valid = (a != nod) if nod is not None else np.isfinite(a)
        v = a[valid]
        print(f"valid px {v.size}/{a.size} ({100 * v.size / a.size:.1f}%)")
        if v.size:
            print(f"min {v.min():.3f} max {v.max():.3f} | zeros {(v == 0).sum()} | neg {(v < 0).sum()}")
            print("pct[1,5,25,50,75,95,99]:",
                  [round(float(np.percentile(v, p)), 2) for p in (1, 5, 25, 50, 75, 95, 99)])
            print("px<=24h:", int((v <= 24).sum()), "| <=48h:", int((v <= 48).sum()))


# ---------------------------------------------------------------- geometry

def _parts(gj):
    """[(geometry, properties)] of a snapshot, however Cornea wrapped it."""
    if gj.get("type") == "FeatureCollection":
        g = gj
    else:
        g = gj.get("geometry") if "geometry" in gj else gj
    g = g or {}
    if g.get("type") == "FeatureCollection":
        return [((f or {}).get("geometry") or {}, (f or {}).get("properties") or {})
                for f in g.get("features") or []]
    return [(g, gj.get("properties") or {})]


def _ring_acres(ring):
    """Approximate area of one ring: shoelace on degrees, scaled at the ring's latitude."""
    a = np.asarray([p[:2] for p in ring], dtype=float)
    if len(a) < 3:
        return 0.0
    deg2 = abs(np.dot(a[:, 0], np.roll(a[:, 1], -1)) - np.dot(a[:, 1], np.roll(a[:, 0], -1))) / 2
    return deg2 * 111320.0 ** 2 * np.cos(np.radians(a[:, 1].mean())) / SQM_PER_ACRE


def _acres(parts):
    """Polygon area in acres, measured from the geometry: poly_GISAcres can't be used, as
    redraws leave it stale too. Only compared between edits of one perimeter."""
    total = 0.0
    for g, _ in parts:
        polys = {"Polygon": [g.get("coordinates")], "MultiPolygon": g.get("coordinates")}.get(g.get("type")) or []
        for rings in polys:
            if rings:
                total += max(0.0, _ring_acres(rings[0]) - sum(_ring_acres(r) for r in rings[1:]))
    return total


class PerimeterSet:
    """Perimeter timeline for one fire + lazy per-grid rasterization cache.

    Each Cornea snapshot is one edit of a WFIGS polygon, keyed by its edit time. Every
    snapshot gets the most reliable time it has: when it was mapped (poly_PolygonDateTime),
    if that is plausible and the polygon still matches the first edit that carried it
    (later matching edits are corrections of that capture); otherwise when it was last
    edited (published), which bounds when it was mapped. poly_CreateDate adds nothing: it
    is either the upload time or the polygon record's original creation, often days before."""

    def __init__(self, slug):
        self.slug = slug
        self.timeline = []      # sorted [(time, epochms, basis)], basis "mapped" or "edited";
                                # a mapped capture appears once, as its latest correction
        self.corrections = {}   # mapped time -> sorted [(edit dt, epochms)] of that capture
        self.published = []     # sorted [(edit dt, epochms)]: every snapshot
        self.time_of = {}       # epochms -> (time, basis) for every snapshot
        self.republished = set()   # edited snapshots that re-publish an earlier snapshot's shape
        self.unreadable = 0
        by_capture, shape_of = {}, {}
        idx = load_json(RAW / "perimeter_archive" / slug / "index.json")
        for item in (idx or {}).get("index", []):
            path, date = item.get("path"), item.get("date")
            if not path or not date:
                continue
            epochms = path.rsplit("/", 1)[-1]
            try:
                gj = load_json(RAW / "perimeter_archive" / slug / f"{epochms}.geojson")
                if gj is None:
                    continue
                parts = _parts(gj)
                acres = _acres(parts)
            except (ValueError, TypeError, AttributeError, IndexError) as e:
                # one bad snapshot must not stop the whole run: leave it out of the timeline
                log(f"{slug}/{epochms}.geojson unreadable, skipped: {e}")
                self.unreadable += 1
                continue
            pr = parts[0][1] if parts else {}
            edited = parse_iso(date)
            self.published.append((edited, epochms))
            shape_of[epochms] = hash(json.dumps([g for g, _ in parts], sort_keys=True))
            cap, disc = pr.get("poly_PolygonDateTime"), pr.get("attr_FireDiscoveryDateTime")
            # plausible: at most CAPTURE_LATE_S after its own edit, and no earlier than a day
            # before discovery (61 days before the edit if discovery is unknown): catches
            # typos like year 0926
            lo = (disc / 1000 if isinstance(disc, (int, float))
                  else edited.timestamp() - 60 * 86400) - 86400
            if (isinstance(cap, (int, float)) and np.isfinite(acres)
                    and lo <= cap / 1000 <= edited.timestamp() + CAPTURE_LATE_S):
                by_capture.setdefault(int(cap), []).append((edited, epochms, acres))
            else:
                self._add(edited, epochms, "edited")
        for cap, edits in by_capture.items():
            edits.sort()
            mapped, a0 = datetime.fromtimestamp(cap / 1000, timezone.utc), edits[0][2]
            tol = min(CAPTURE_REL_TOL * a0, CAPTURE_TOL_AC)
            same = [(edited, e) for edited, e, area in edits if a0 > 0 and abs(area - a0) <= tol]
            for edited, e, area in edits:
                if (edited, e) in same:
                    self.time_of[e] = (mapped, "mapped")
                else:
                    self._add(edited, e, "edited")   # redrawn under a stale capture time
            if same:
                self.corrections[mapped] = same
                self.timeline.append((mapped, same[-1][1], "mapped"))
        self.timeline.sort()
        self.published.sort()
        seen = set()
        for _, e in self.published:
            if shape_of.get(e) in seen and self.time_of.get(e, (None, ""))[1] == "edited":
                self.republished.add(e)
            seen.add(shape_of.get(e))
        self._geom_utm = {}   # (epochms, crs_str) -> shapely geom in raster CRS
        self._mask = {}       # (epochms, grid_sig) -> bool ndarray

    def _add(self, t, epochms, basis):
        self.timeline.append((t, epochms, basis))
        self.time_of[epochms] = (t, basis)

    def latest_at_or_before(self, dt):
        """The most recent perimeter known to be mapped by dt -> (time, epochms, basis), as
        published by dt: of a mapped capture, the latest correction out by then (else its
        first edit), so later growth can't creep in. An edit time is safe here too: it is
        no earlier than the mapping it stands in for."""
        best = None
        for t, e, basis in self.timeline:
            if t <= dt and (best is None or t > best[0]):
                best = (t, e, basis)
        if best and best[2] == "mapped":
            edits = self.corrections[best[0]]
            out = [e for edited, e in edits if edited <= dt]
            best = (best[0], out[-1] if out else edits[0][1], "mapped")
        return best

    def nearest(self, dt, tol_h):
        """The perimeter showing the fire at dt -> (time, epochms, offset_h, basis), by
        reliability: the one mapped closest to dt within tol_h; else the one last edited
        closest to it (its mapped time unknown); else the one published closest to dt, at
        its mapped time even though that is further off (a late re-edit of an older map).
        A mapped capture is always its latest correction."""
        for want in ("mapped", "edited"):
            best = None
            for t, e, basis in self.timeline:
                off = (t - dt).total_seconds() / 3600.0
                if basis == want and abs(off) <= tol_h and (best is None or abs(off) < abs(best[2])):
                    best = (t, e, off, basis)
            if best:
                return best
        pub = [(abs((t - dt).total_seconds()), e) for t, e in self.published
               if abs((t - dt).total_seconds()) <= tol_h * 3600]
        if not pub:
            return None
        e = min(pub)[1]
        t, basis = self.time_of[e]
        if basis == "mapped":
            e = self.corrections[t][-1][1]
        return (t, e, (t - dt).total_seconds() / 3600.0, basis)

    def geom_utm(self, epochms, crs):
        from rasterio.warp import transform_geom
        from shapely.geometry import shape
        key = (epochms, str(crs))
        if key not in self._geom_utm:
            gj = load_json(RAW / "perimeter_archive" / self.slug / f"{epochms}.geojson")
            geom = gj.get("geometry", gj)
            if geom.get("type") == "FeatureCollection":
                from shapely.ops import unary_union
                parts = [shape(transform_geom("EPSG:4326", crs, f["geometry"]))
                         for f in geom["features"]]
                g = unary_union([p.buffer(0) for p in parts])
            else:
                g = shape(transform_geom("EPSG:4326", crs, geom)).buffer(0)
            self._geom_utm[key] = g
        return self._geom_utm[key]

    def mask(self, epochms, crs, transform, shape_):
        from rasterio import features
        from shapely.geometry import mapping
        key = (epochms, str(crs), tuple(transform)[:6], shape_)
        if key not in self._mask:
            g = self.geom_utm(epochms, crs)
            if g.is_empty:
                self._mask[key] = np.zeros(shape_, dtype=bool)
            else:
                self._mask[key] = features.rasterize(
                    [(mapping(g), 1)], out_shape=shape_, transform=transform,
                    fill=0, dtype="uint8").astype(bool)
        return self._mask[key]


class HotspotSet:
    """Timestamped VIIRS/MODIS detections for one fire. Unlike perimeters these
    carry exact acquisition times, so windows need no tolerance: a detection
    either falls inside (run, run+H] or it doesn't. Footprints are disks sized by
    instrument (see HS_RADIUS_VIIRS: VIIRS 187.5 m, MODIS 250 m, Landsat 15 m).

    Every detection is screened once against the fire's whole record (see
    HS_ARTIFACT_M): only those that pass are scored; the rest are left out of
    the metric, neither captured nor missed."""

    def __init__(self, slug, ps=None):
        self.pts = []  # (dt, lon, lat, radius_m, source)
        for p in sorted((RAW / "hotspot_archive" / slug).glob("*.geojson")):
            try:
                gj = json.loads(p.read_text())
            except ValueError:
                continue
            for f in gj.get("features", []):
                pr = f.get("properties") or {}
                d, c = pr.get("acq_date"), (f.get("geometry") or {}).get("coordinates")
                if not d or not c:
                    continue
                try:
                    hm = int(float(pr.get("acq_time") or 0))
                except ValueError:
                    hm = 0
                dt = datetime.strptime(d, "%Y-%m-%d").replace(
                    hour=hm // 100 % 24, minute=hm % 100 % 60, tzinfo=timezone.utc)
                src = (pr.get("source") or "").lower()
                r = (HS_RADIUS_MODIS if ("terra" in src or "aqua" in src or "modis" in src)
                     else HS_RADIUS_LANDSAT if "landsat" in src
                     else HS_RADIUS_VIIRS)
                self.pts.append((dt, c[0], c[1], r, src))
        self.pts.sort(key=lambda x: (x[0], x[4]))
        # "ok", "artifact" or "transient", per detection
        self.status = (self._screen(ps) if ps is not None and ps.timeline and self.pts
                       else ["ok"] * len(self.pts))
        self._mask = {}

    def _screen(self, ps):
        import shapely
        from rasterio.warp import transform as warp_transform
        n = len(self.pts)
        lons = np.array([p[1] for p in self.pts])
        lats = np.array([p[2] for p in self.pts])
        zone = int((np.median(lons) + 180) // 6) % 60 + 1
        crs = f"EPSG:{(32600 if np.median(lats) >= 0 else 32700) + zone}"
        xs, ys = (np.asarray(v) for v in warp_transform("EPSG:4326", crs, lons.tolist(), lats.tolist()))
        status = np.full(n, "ok", dtype=object)

        # 1. artifacts: far from the perimeter current at detection, never inside a later one
        current, groups = {}, {}
        for i, p in enumerate(self.pts):
            if p[0] not in current:
                c = ps.latest_at_or_before(p[0])
                current[p[0]] = c[1] if c else None
            if current[p[0]] is not None:
                groups.setdefault(current[p[0]], []).append(i)

        def near(epochms, idx, dist):
            """which detections idx lie within dist (m, scalar or per detection) of a perimeter"""
            # prepared containment first, distance only for the few outside: plain
            # distance() on a 47k-vertex perimeter cost ~17s for one fire
            g = ps.geom_utm(epochms, crs)
            shapely.prepare(g)
            pts = shapely.points(xs[idx], ys[idx])
            out = shapely.intersects(g, pts)
            rest = ~out
            if rest.any():
                out[rest] = shapely.dwithin(g, pts[rest], np.broadcast_to(dist, out.shape)[rest])
            return out

        far = np.zeros(n, dtype=bool)
        for epochms, idx in groups.items():
            idx = np.asarray(idx)
            far[idx] = ~near(epochms, idx, HS_ARTIFACT_M)
        far_idx = np.flatnonzero(far)
        if far_idx.size:
            t = np.array([self.pts[i][0].timestamp() for i in far_idx])
            r = np.array([self.pts[i][3] for i in far_idx])
            judged = np.zeros(far_idx.size, dtype=bool)    # some perimeter was mapped after it
            touched = np.zeros(far_idx.size, dtype=bool)   # and its footprint touches one of those
            for pt, epochms, _ in ps.timeline:
                later = t < pt.timestamp()
                judged |= later
                chk = later & ~touched
                if chk.any():
                    touched[chk] = near(epochms, far_idx[chk], r[chk])
            status[far_idx[judged & ~touched]] = "artifact"

        # 2. transients, among the rest: a pass is one source's detections less than
        # HS_PASS_GAP_S apart; a detection needs one within HS_PERSIST_M on some
        # source's nearest pass before or after its own (its own source's included)
        by_src = {}
        for i in np.flatnonzero(status == "ok"):
            by_src.setdefault(self.pts[i][4], []).append(i)
        passes = {}   # source -> its passes in time order
        for s, idx in by_src.items():
            runs, run = [], [idx[0]]
            for a, b in zip(idx, idx[1:]):
                if (self.pts[b][0] - self.pts[a][0]).total_seconds() > HS_PASS_GAP_S:
                    runs.append(run)
                    run = []
                run.append(b)
            runs.append(run)
            passes[s] = [np.asarray(r) for r in runs]
        starts = {s: np.array([self.pts[r[0]][0].timestamp() for r in pl])
                  for s, pl in passes.items()}
        r2 = HS_PERSIST_M ** 2

        def confirmed(a, b):
            """which of pass a's detections have one of pass b's within HS_PERSIST_M"""
            out = np.zeros(a.size, dtype=bool)
            step = max(1, 2_000_000 // b.size)   # caps the pairwise block at ~16 MB
            for s in range(0, a.size, step):
                aa = a[s:s + step]
                out[s:s + step] = (((xs[aa, None] - xs[None, b]) ** 2 +
                                    (ys[aa, None] - ys[None, b]) ** 2) <= r2).any(axis=1)
            return out

        for s0, pl0 in passes.items():
            for k, ix in enumerate(pl0):
                nb = []
                for s, pl in passes.items():
                    if s == s0:
                        j = k   # its own source: the passes either side of this one
                        before, after = j - 1, j + 1
                    else:       # another: its last pass starting before, first after
                        j = int(np.searchsorted(starts[s], starts[s0][k]))
                        before, after = j - 1, j
                    nb += [pl[x] for x in (before, after) if 0 <= x < len(pl)]
                ok = (confirmed(ix, np.concatenate(nb)) if nb
                      else np.zeros(ix.size, dtype=bool))
                status[ix[~ok]] = "transient"
        return status.tolist()

    def mask(self, t0, t1, crs, transform, shape_, screened=True):
        """-> (footprint mask of the screened-in detections in (t0, t1], how many there
        are, {"artifact": n, "transient": n} of those in the window left out).
        screened=False takes every detection (RunStart: what a run was seeded from)."""
        from rasterio import features
        from rasterio.warp import transform as warp_transform
        from shapely.geometry import Point, mapping
        key = (t0, t1, str(crs), tuple(transform)[:6], shape_, screened)
        if key not in self._mask:
            win = [i for i, p in enumerate(self.pts) if t0 < p[0] <= t1]
            sel = [self.pts[i] for i in win if self.status[i] == "ok" or not screened]
            left_out = {"artifact": 0, "transient": 0}
            for i in win:
                if self.status[i] != "ok" and screened:
                    left_out[self.status[i]] += 1
            if not sel:
                self._mask[key] = (np.zeros(shape_, dtype=bool), 0, left_out)
            else:
                xs, ys = warp_transform("EPSG:4326", crs,
                                        [p[1] for p in sel], [p[2] for p in sel])
                geoms = [(mapping(Point(x, y).buffer(p[3], resolution=8)), 1)
                         for x, y, p in zip(xs, ys, sel)]
                # and the cell each falls in: a 15 m Landsat disk misses every cell
                # centre of a 30 m grid about 1 time in 5
                geoms += [(mapping(Point(x, y)), 1) for x, y in zip(xs, ys)]
                m = features.rasterize(geoms, out_shape=shape_, transform=transform,
                                       fill=0, dtype="uint8").astype(bool)
                self._mask[key] = (m, len(sel), left_out)
        return self._mask[key]


# ---------------------------------------------------------------- metrics

def acres(npix, px_area_sqm):
    return round(npix * px_area_sqm / SQM_PER_ACRE, 1)


def on_grid(arr, src_crs, src_tfm, crs, tfm, shape_):
    """a bool mask on another grid (nearest cell); as is when the grids match"""
    if arr.shape == shape_ and tuple(src_tfm)[:6] == tuple(tfm)[:6] and str(src_crs) == str(crs):
        return arr
    from rasterio.warp import Resampling, reproject
    out = np.zeros(shape_, dtype=np.uint8)
    reproject(arr.astype(np.uint8), out, src_transform=src_tfm, src_crs=src_crs,
              dst_transform=tfm, dst_crs=crs, resampling=Resampling.nearest)
    return out.astype(bool)


class RunStart:
    """What had burned when one run began, beyond the baseline perimeter (see START_VAR):
    Pyrecast's own start if its granule was archived, else every satellite detection
    footprint up to run_ts, else nothing. -> mask(crs, transform, shape) on a forecast grid."""

    def __init__(self, run_key, run_dt, hs):
        import rasterio
        self.basis, self._src = "", None
        p = RAW / "forecast_archive" / run_key / START_FILE
        if p.exists():
            with rasterio.open(p) as ds:
                self._src = (np.isin(ds.read(1), START_CODES), ds.crs, ds.transform)
            self.basis = "pyrecast"
        elif any(t <= run_dt for t, *_ in hs.pts):
            self.basis = "hotspots"
        self._run_dt, self._hs, self._grids = run_dt, hs, {}

    def mask(self, crs, tfm, shape_):
        key = (str(crs), tuple(tfm)[:6], shape_)
        if key not in self._grids:
            if self.basis == "pyrecast":
                m = on_grid(*self._src, crs, tfm, shape_)
            elif self.basis == "hotspots":
                m = self._hs.mask(datetime(1970, 1, 1, tzinfo=timezone.utc), self._run_dt,
                                  crs, tfm, shape_, screened=False)[0]
            else:
                m = np.zeros(shape_, dtype=bool)
            self._grids[key] = m
        return self._grids[key]


def write_start_mask(path, m, crs, tfm):
    """the start footprint a run was scored with, for map.html and report.html to redraw
    the calculation: uint8 0/1 on the forecast's grid"""
    import rasterio
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=m.shape[0], width=m.shape[1], count=1,
                       dtype="uint8", crs=crs, transform=tfm, compress="deflate") as ds:
        ds.write(m.astype(np.uint8), 1)


def scored_perimeters(ps, run_dt):
    """The perimeters a run is scored against, each at its own time -> sorted
    [(time, epochms, basis, chosen, day)], day = 24 h periods from the run's start:
      mapped    every perimeter with a mapped time in the run's first SCORE_DAYS days
      daily     on a day with none, the latest perimeter published that day: timed by its
                last edit (its mapped time unknown), not a re-publish of an earlier shape
      previous  the perimeter the previous method scored at each of its lead times (nearest
                run + H within ACTUAL_TOL_H, shown when within OLD_MAX_OFFSET_H), so no analysis
                shown before is lost; it may lie beyond SCORE_DAYS (the old 10 and 14 days)"""
    end = run_dt + timedelta(days=SCORE_DAYS)
    day = lambda t: int(((t - run_dt).total_seconds() - 1e-6) // 86400) + 1
    out = {}
    for t, e, b in ps.timeline:
        if b == "mapped" and run_dt < t <= end:
            out[e] = (t, e, b, "mapped")
    timed_days = {day(t) for t, *_ in out.values()}
    for k in range(1, SCORE_DAYS + 1):
        if k in timed_days:
            continue
        d0, d1 = run_dt + timedelta(days=k - 1), run_dt + timedelta(days=k)
        cands = [(t, e, b) for t, e, b in ps.timeline
                 if b == "edited" and d0 < t <= d1 and e not in ps.republished]
        if cands:
            t, e, b = max(cands)
            out[e] = (t, e, b, "daily")
    for H in OLD_HORIZONS_H:
        a = ps.nearest(run_dt + timedelta(hours=H), ACTUAL_TOL_H)
        if a and abs(a[2]) <= OLD_MAX_OFFSET_H and a[0] > run_dt and a[1] not in out:
            out[a[1]] = (a[0], a[1], a[3], "previous")
    return sorted((t, e, b, c, day(t)) for t, e, b, c in out.values())


def score_grid(toa_h, valid, crs, tfm, grid_box, ps, hs, run_dt, perims, before, b_geom, start_extra):
    """Growth-only scores of one forecast (time of arrival in hours since run_dt on a grid)
    against each perimeter in perims (scored_perimeters), at the perimeter's own time: the
    predicted area is where the forecast has burned by then. The same for every model.
      before       bool grid mask of what had burned when the run began: the baseline
                   perimeter plus every model's start for the run (one common start)
      b_geom       the baseline perimeter in the grid's CRS, or None
      start_extra  bool grid mask of the part of `before` beyond the baseline
    No growth: precision 0 and IoU 0 if the forecast predicted growth, recall blank; all
    blank if neither grew. -> [dict] of the per-perimeter columns"""
    px_area = abs(tfm.a * tfm.e)
    out = []
    for a_dt, a_ems, a_basis, chosen, day in perims:
        lead = (a_dt - run_dt).total_seconds() / 3600.0
        A = ps.mask(a_ems, crs, tfm, toa_h.shape)
        pred = valid & (toa_h <= lead)
        pred_new = pred & ~before
        act_new = A & ~before
        inter = pred_new & act_new
        n_pred, n_act, n_int = int(pred_new.sum()), int(act_new.sum()), int(inter.sum())
        n_union = n_pred + n_act - n_int
        hs_mask, hs_n, hs_left_out = hs.mask(run_dt, a_dt, crs, tfm, toa_h.shape)
        hs_new = hs_mask & ~before
        n_hs_new = int(hs_new.sum())
        n_hs_hit = int((pred_new & hs_new).sum())
        # the share of the growth beyond the forecast's grid (recorded, not filtered on): the
        # start footprint lies on the grid, so only the in-grid part loses its cells
        a_geom = ps.geom_utm(a_ems, crs)
        act_new_geom = a_geom.difference(b_geom) if b_geom else a_geom
        outside = act_new_geom.difference(grid_box)
        new_area = max(outside.area, act_new_geom.area - int((A & start_extra).sum()) * px_area)
        out.append({
            "lead_h": round(lead, 2), "day": day, "chosen": chosen,
            "actual_perim_ms": a_ems, "actual_time_utc": a_dt.strftime(ISO_Z),
            "actual_time_basis": a_basis,
            "pred_total_acres": acres(int(pred.sum()), px_area),
            "pred_new_acres": acres(n_pred, px_area),
            "act_new_acres": acres(n_act, px_area),
            "inter_acres": acres(n_int, px_area),
            "precision": round(n_int / n_pred, 4) if n_pred else "",
            "recall": round(n_int / n_act, 4) if n_act else "",
            "iou": round(n_int / n_union, 4) if n_union else "",
            "act_outside_acres": round(outside.area / SQM_PER_ACRE, 1),
            "act_offgrid_frac": round(outside.area / new_area, 4) if new_area > 0 else 0.0,
            "hs_n": hs_n,
            "hs_new_acres": acres(n_hs_new, px_area),
            "hs_hit_acres": acres(n_hs_hit, px_area),
            "hs_capture": round(n_hs_hit / n_hs_new, 4) if n_hs_new else "",
            "hs_artifact_n": hs_left_out["artifact"],
            "hs_transient_n": hs_left_out["transient"],
        })
    return out


def compute(only_slug=None, on_fire=None):
    """Score every archived run, newest fires first (a fire's newest run decides). on_fire,
    if given, is called with (slug, rows) as each fire finishes, to publish it right away."""
    import rasterio
    from shapely.geometry import box

    manifest = load_json(RAW / "manifest.json")
    matches = load_json(RAW / "fire_matches.json")
    if not manifest or not matches:
        log("FATAL raw/manifest.json or raw/fire_matches.json missing — run collect.py --local or --sync-url first")
        sys.exit(1)

    perimsets = {}
    hssets = {}
    rows = []
    skipped = {"no_match": 0, "no_perims": 0, "no_perimeter_to_score": 0, "missing_tif": 0,
               "unreadable_perims": 0}
    starts = {"pyrecast": 0, "hotspots": 0, "": 0}
    # each run's start footprint, as scored, for the pages (published by publish_skill)
    import shutil
    shutil.rmtree(START_DIR / only_slug if only_slug else START_DIR, ignore_errors=True)

    # newest fires first (by each fire's newest run), a fire's runs together
    newest = {}
    for e in manifest["runs"].values():
        newest[e["slug"]] = max(newest.get(e["slug"], ""), e["run_ts"])
    order = sorted(manifest["runs"], key=lambda k: (newest[manifest["runs"][k]["slug"]],
                                                    manifest["runs"][k]["slug"], k), reverse=True)
    fire_rows, cur_slug = [], None

    def fire_done():
        if on_fire and cur_slug and fire_rows:
            on_fire(cur_slug, list(fire_rows))
        fire_rows.clear()

    for run_key in order:
        entry = manifest["runs"][run_key]
        slug = entry["slug"]
        if only_slug and slug != only_slug:
            continue
        if slug != cur_slug:
            fire_done()
            cur_slug = slug
        m = matches["matches"].get(slug)
        if not m:
            skipped["no_match"] += 1
            continue
        if slug not in perimsets:
            # sorted run_keys iterate grouped by slug and no cache is reused
            # across fires, so evicting here caps peak memory at one fire's
            # masks instead of the whole fleet's (OOM-killed CI runners Sep 7+)
            perimsets.clear()
            hssets.clear()
            perimsets[slug] = PerimeterSet(slug)
            skipped["unreadable_perims"] += perimsets[slug].unreadable
        ps = perimsets[slug]
        if not ps.timeline:
            skipped["no_perims"] += 1
            continue
        if slug not in hssets:
            hssets[slug] = HotspotSet(slug, ps)
        hs = hssets[slug]
        hs._mask.clear()  # keyed per grid+window; new run = new grid, drop old

        run_dt = datetime.strptime(entry["run_ts"], "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
        baseline = ps.latest_at_or_before(run_dt + timedelta(hours=BASELINE_TOL_H))
        perims = scored_perimeters(ps, run_dt)
        if not perims:
            skipped["no_perimeter_to_score"] += 1
            continue
        start = RunStart(run_key, run_dt, hs)
        starts[start.basis] += 1
        start_written = False

        for pct in PERCENTILES:
            if not entry["files"].get(pct, {}).get("ok"):
                continue
            tif = RAW / "forecast_archive" / run_key / f"{pct}.tif"
            if not tif.exists():
                skipped["missing_tif"] += 1
                continue
            with rasterio.open(tif) as ds:
                arr = ds.read(1)
                crs, tfm = ds.crs, ds.transform
                nod = ds.nodata
                grid_box = box(*ds.bounds)
            px_area = abs(tfm.a * tfm.e)
            valid = (arr != nod) if nod is not None else np.isfinite(arr)
            toa_h = arr * TOA_TO_HOURS

            B = (ps.mask(baseline[1], crs, tfm, arr.shape)
                 if baseline else np.zeros(arr.shape, dtype=bool))
            b_geom = ps.geom_utm(baseline[1], crs) if baseline else None
            S = start.mask(crs, tfm, arr.shape)
            BS = B | S   # burned before: what nothing after the run can count as growth
            if S.any() and not start_written:
                write_start_mask(START_DIR / slug / f"{entry['run_ts']}.tif", S, crs, tfm)
                start_written = True

            notes = [] if baseline else ["baseline=none"]
            for sc in score_grid(toa_h, valid, crs, tfm, grid_box, ps, hs, run_dt, perims,
                                 BS, b_geom, S & ~B):
                row = {
                    "model": MODEL, "slug": slug, "cornea_id": m["cornea_id"],
                    "run_ts": entry["run_ts"],
                    "run_dt_utc": run_dt.strftime(ISO_Z),
                    "percentile": pct,
                    "baseline_perim_ms": baseline[1] if baseline else "",
                    "baseline_time_utc": baseline[0].strftime(ISO_Z) if baseline else "",
                    "baseline_time_basis": baseline[2] if baseline else "",
                    "baseline_offset_h": (round((baseline[0] - run_dt).total_seconds() / 3600, 2)
                                          if baseline else ""),
                    "pixel_m": round(abs(tfm.a), 1),
                    "baseline_acres": acres(int(B.sum()), px_area),
                    "start_basis": start.basis,
                    "start_acres": acres(int(S.sum()), px_area),
                    "burned_before_acres": acres(int(BS.sum()), px_area),
                    "expired_run": entry.get("expired", False),
                    "notes": ";".join(notes),
                    **sc,
                }
                rows.append(row)
                fire_rows.append(row)
        log(f"{run_key}: {len(perims)} perimeters x pcts done")
    fire_done()

    out = REPO / "data" / "metrics.csv"
    out.parent.mkdir(exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {len(rows)} rows -> {out} | skipped runs: {skipped} | start footprint per "
        f"scored run: {starts['pyrecast']} Pyrecast's own, {starts['hotspots']} from detections, "
        f"{starts['']} none")


# ---------------------------------------------------------------- skill files

def _num(s):
    return float(s) if s else None


def _r3(s):
    v = _num(s)
    return None if v is None else round(v, 3)


def _perim_id(s):
    return int(float(s)) if s else None


def skill_entry(r):
    """A metrics row as a skill-file row: [lead_h, precision, recall, iou, baseline_perim_ms,
    actual_perim_ms, 0, actual_mapped]. The scores, then what the map needs to redraw the
    calculation: the perimeter archive ids (null baseline when none was mapped before the
    run), the scored perimeter's hours from the window end (always 0 now: each perimeter is
    scored at its own time; v3 files carry the old offset), and 1 if that time is when it
    was mapped, 0 if it's its last edit."""
    basis = r.get("actual_time_basis")
    return [round(float(r["lead_h"]), 2),
            _r3(r["precision"]), _r3(r["recall"]), _r3(r["iou"]),
            _perim_id(r["baseline_perim_ms"]), _perim_id(r["actual_perim_ms"]),
            0, None if not basis else int(basis == "mapped")]


def skill_doc(slug, rows, generated):
    """map.html's skill file for one fire from its metrics rows (every row: no filters).
    v4: rows at each perimeter's exact lead (skill_entry); "starts": {run_ts: "pyrecast" |
    "hotspots"} for the runs scored with a start footprint beyond the baseline (see
    START_VAR), published as skill/start/<slug>/<run_ts>.tif, which the pages add to the
    baseline perimeter."""
    runs, starts = {}, {}
    for r in rows:
        if r.get("model", MODEL) != MODEL:
            continue
        if r.get("start_basis") and (_num(str(r.get("start_acres") or "")) or 0) > 0:
            starts[r["run_ts"]] = r["start_basis"]
        runs.setdefault(r["run_ts"], {}).setdefault(str(r["percentile"]), []).append(skill_entry(r))
    runs = {run_ts: {pct: sorted(by_pct[pct], key=lambda e: e[0]) for pct in sorted(by_pct, key=int)}
            for run_ts, by_pct in sorted(runs.items())}
    n = sum(len(v) for by_pct in runs.values() for v in by_pct.values())
    return {"v": 4, "slug": slug, "generated": generated, "method": "exact-lead",
            "score_days": SCORE_DAYS, "starts": dict(sorted(starts.items())), "runs": runs}, n


def build_skill(csv_path, out_dir):
    """Per-fire skill files for map.html from a metrics CSV: out_dir/<slug>.json (skill_doc)
    + out_dir/index.json {slug: rows}. Returns (fires, rows)."""
    by_slug = {}
    with open(csv_path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("lead_h") in (None, ""):
                continue      # an older CSV (fixed horizons): rebuild from a full compute
            by_slug.setdefault(r["slug"], []).append(r)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for p in out_dir.glob("*.json"):
        p.unlink()
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    index = {}
    for slug in sorted(by_slug):
        doc, n = skill_doc(slug, by_slug[slug], generated)
        if not n:
            continue
        index[slug] = n
        (out_dir / f"{slug}.json").write_text(json.dumps(doc, separators=(",", ":")))
    (out_dir / "index.json").write_text(json.dumps(
        {"v": 1, "generated": generated, "slugs": index}, separators=(",", ":")))
    return len(index), sum(index.values())


class FirePublisher:
    """--publish: each fire's skill file (and its runs' start footprints) goes up as soon
    as compute() finishes it, newest fires first, so recent fires update first. The index
    goes up merged with the bucket's (fires not reached yet keep their last files) at most
    every SKILL_INDEX_EVERY_S; the full publish after compute() replaces everything."""

    def __init__(self, client, bucket):
        import requests
        self.client, self.bucket = client, bucket
        try:
            self.index = http_get(requests.Session(), f"{PUBLIC_BUCKET_URL}/skill/index.json").json()["slugs"]
        except Exception as e:
            log(f"no published skill index to merge into ({e}); starting empty")
            self.index = {}
        self.generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        self.last_index, self.n = 0.0, 0

    def __call__(self, slug, rows):
        try:
            for p in sorted((START_DIR / slug).glob("*.tif")):
                self.client.upload_file(str(p), self.bucket, f"skill/start/{slug}/{p.name}",
                                        ExtraArgs={"ContentType": "image/tiff"})
            doc, n = skill_doc(slug, rows, self.generated)
            if not n:
                return
            self.client.put_object(Bucket=self.bucket, Key=f"skill/{slug}.json",
                                   Body=json.dumps(doc, separators=(",", ":")).encode(),
                                   ContentType="application/json")
            self.index[slug] = n
            self.n += 1
            if time.time() - self.last_index >= SKILL_INDEX_EVERY_S:
                self.put_index()
        except Exception as e:   # the full publish at the end still covers this fire
            log(f"early publish of {slug} failed: {e}")

    def put_index(self):
        self.client.put_object(Bucket=self.bucket, Key="skill/index.json",
                               Body=json.dumps({"v": 1, "generated": self.generated, "slugs": self.index},
                                               separators=(",", ":")).encode(),
                               ContentType="application/json")
        self.last_index = time.time()
        log(f"early publish: {self.n} fires up so far")


def publish_starts(client, bucket):
    """Upload the start footprint of every run compute() scored with one, to
    skill/start/<slug>/<run_ts>.tif. Every run, not just the ones the skill files name:
    the report redraws whichever row it's asked to, default filters or not. Called before
    metrics.csv goes up, so a page never meets a row whose footprint isn't there yet; a
    --skill-only rebuild has none here, and the bucket keeps the last full run's."""
    masks = sorted(START_DIR.glob("*/*.tif"))

    def put_mask(p):
        client.upload_file(str(p), bucket, f"skill/start/{p.parent.name}/{p.name}",
                           ExtraArgs={"ContentType": "image/tiff"})

    with concurrent.futures.ThreadPoolExecutor(SKILL_UPLOAD_WORKERS) as ex:
        list(ex.map(put_mask, masks))
    log(f"published {len(masks)} start footprints")
    return len(masks)


def publish_skill(client, bucket, csv_path):
    # the start footprints the skill files name are already up (publish_starts)
    n_fires, n_rows = build_skill(csv_path, SKILL_DIR)
    index = SKILL_DIR / "index.json"
    files = sorted(SKILL_DIR.glob("*.json"))

    def put(p):
        client.upload_file(str(p), bucket, f"skill/{p.name}",
                           ExtraArgs={"ContentType": "application/json"})

    with concurrent.futures.ThreadPoolExecutor(SKILL_UPLOAD_WORKERS) as ex:
        list(ex.map(put, [p for p in files if p != index]))
    # index last, so it never lists a slug whose file is not uploaded yet
    put(index)
    log(f"published skill files: {n_fires} fires, {n_rows} rows, {len(files)} files")


def skill_only():
    import collect
    import requests
    # before the download, so missing S3_* env fails fast
    client, bucket = collect.make_r2()
    r = http_get(requests.Session(), f"{PUBLIC_BUCKET_URL}/metrics.csv")
    dst = RAW / "metrics_published.csv"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(r.content)
    log(f"downloaded published metrics.csv ({len(r.content)} bytes) -> {dst}")
    publish_skill(client, bucket, dst)


# ---------------------------------------------------------------- overlay

def overlay(slug, run_ts, pct, H):
    import rasterio
    from rasterio import features
    from rasterio.warp import transform_geom
    from shapely.geometry import mapping, shape

    manifest = load_json(RAW / "manifest.json")
    entry = manifest["runs"][f"{slug}/{run_ts}"]
    ps = PerimeterSet(slug)
    run_dt = datetime.strptime(run_ts, "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
    baseline = ps.latest_at_or_before(run_dt + timedelta(hours=BASELINE_TOL_H))
    actual = ps.nearest(run_dt + timedelta(hours=H), ACTUAL_TOL_H)
    tif = RAW / "forecast_archive" / slug / run_ts / f"{pct}.tif"
    with rasterio.open(tif) as ds:
        arr = ds.read(1)
        crs, tfm = ds.crs, ds.transform
        nod = ds.nodata
    valid = (arr != nod) if nod is not None else np.isfinite(arr)
    B = ps.mask(baseline[1], crs, tfm, arr.shape) if baseline else np.zeros(arr.shape, bool)
    B |= RunStart(f"{slug}/{run_ts}", run_dt, HotspotSet(slug)).mask(crs, tfm, arr.shape)
    pred_new = (valid & (arr * TOA_TO_HOURS <= H)) & ~B

    feats = []
    for geom, val in features.shapes(pred_new.astype("uint8"), transform=tfm):
        if val == 1:
            feats.append({"type": "Feature",
                          "properties": {"role": "pred_new", "H": H, "pct": pct},
                          "geometry": transform_geom(crs, "EPSG:4326", geom)})
    for role, sel in (("baseline", baseline), ("actual", actual)):
        if sel:
            gj = load_json(RAW / "perimeter_archive" / slug / f"{sel[1]}.geojson")
            feats.append({"type": "Feature",
                          "properties": {"role": role, "date": sel[0].isoformat()},
                          "geometry": gj.get("geometry", gj)})
    out = RAW / f"overlay_{slug}_{run_ts}_{pct}_{H}h.geojson"
    out.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    log(f"wrote {out} ({len(feats)} features)")


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sync-only", action="store_true",
                    help="mirror the bucket into raw/ (as --sync) and stop: the CI workflow syncs, "
                         "saves its cache, then analyzes")
    ap.add_argument("--sync", action="store_true",
                    help=f"mirror the R2 bucket ({PUBLIC_BUCKET_URL}) into raw/ first")
    ap.add_argument("--sync-url", help="override the public bucket URL for --sync")
    ap.add_argument("--refresh-perims", action="store_true",
                    help="re-fetch Cornea perimeter snapshots before computing")
    ap.add_argument("--inspect", metavar="TIF")
    ap.add_argument("--overlay", nargs=4, metavar=("SLUG", "RUN", "PCT", "H"))
    ap.add_argument("--only", help="restrict metrics (and --sync) to one slug")
    ap.add_argument("--publish", action="store_true",
                    help="after computing, upload data/metrics.csv + skill/*.json to the bucket "
                         "(needs S3_* env; used by the CI analyze workflow)")
    ap.add_argument("--skill-only", action="store_true",
                    help="only rebuild + upload skill/*.json from the published "
                         "metrics.csv, without computing (needs S3_* env)")
    args = ap.parse_args()

    if args.inspect:
        inspect(args.inspect)
        return
    if args.skill_only:
        skill_only()
        return
    if args.sync or args.sync_url or args.sync_only:
        sync_bucket(args.sync_url or PUBLIC_BUCKET_URL, only_slug=args.only)
        if args.sync_only:
            return
    if args.refresh_perims:
        import collect
        archive = collect.Archive(RAW)
        manifest = archive.load_json("manifest.json") or {"runs": {}, "perimeters": {}}
        matches = archive.load_json("fire_matches.json") or {"matches": {}, "unmatched": {}}
        slugs = sorted({e["slug"] for e in manifest["runs"].values()})
        collect.snapshot_perimeters(archive, manifest, matches, slugs)
        archive.save_json("manifest.json", manifest)
        log(f"refreshed perimeters: {collect.counts}")
    if args.overlay:
        slug, run_ts, pct, h = args.overlay
        overlay(slug, run_ts, pct, int(h))
        return
    client = bucket = pub = None
    if args.publish:
        import collect
        client, bucket = collect.make_r2()   # before computing, so missing S3_* env fails fast
        pub = FirePublisher(client, bucket)
    compute(only_slug=args.only, on_fire=pub)
    if not args.publish:   # local / dry run: the skill files the page would get, in data/skill
        log("built skill files (not published): %d fires, %d rows" % build_skill(REPO / "data" / "metrics.csv", SKILL_DIR))
        return
    if args.only:
        # one fire: its skill file and start footprints went up as it finished; the index goes
        # up merged with the bucket's. metrics.csv and every other fire stay as published (a
        # one-fire CSV or index would wipe the rest)
        pub.put_index()
        log(f"published {args.only} only (metrics.csv and other fires unchanged)")
        return
    if args.publish:
        publish_starts(client, bucket)   # before metrics.csv, which names them
        client.upload_file(str(REPO / "data" / "metrics.csv"), bucket, "metrics.csv",
                           ExtraArgs={"ContentType": "text/csv"})
        log("published metrics.csv to bucket")
        publish_skill(client, bucket, REPO / "data" / "metrics.csv")


if __name__ == "__main__":
    main()
