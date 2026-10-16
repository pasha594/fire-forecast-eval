#!/usr/bin/env python3
"""Growth-only precision/recall of archived ELMFIRE forecasts vs Cornea perimeters.

Reads the archive in raw/ (populate with `collect.py --local`, or sync a bucket
with --sync-url) and writes data/metrics.csv. One row per
(fire, forecast run, percentile, horizon) that has ground truth close enough.

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
# so the predicted mask is simply valid & (toa <= H); growth is isolated purely by
# subtracting the Cornea baseline perimeter.
# run_ts is UTC (tif Last-Modified ≈ run_ts + ~3h processing lag).
TOA_TO_HOURS = 1.0

HORIZONS_H = (12, 24, 48, 72, 96, 120, 168, 240, 336)
PUBLIC_BUCKET_URL = "https://f005.backblazeb2.com/file/fire-forecast-archive"
# baseline = latest perimeter mapped <= run + this. run_ts is cut to the minute: Pyrecast
# names a run after the time of the perimeter it started from, seconds dropped
BASELINE_TOL_H = 59.999 / 3600
ACTUAL_TOL_H = 12.0      # actual = perimeter mapped (else edited/published) within +/- this of run + H
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

# must stay equal to report.html's default filter inputs
SKILL_FILTERS = {"min_act_new_acres": 10, "max_abs_offset_h": 3, "max_offgrid_frac": 0.10}
SKILL_UPLOAD_WORKERS = 16

CSV_COLUMNS = ["slug", "cornea_id", "run_ts", "run_dt_utc", "percentile",
               "horizon_h", "valid_dt_utc", "baseline_perim_ms", "baseline_time_utc",
               "baseline_time_basis", "baseline_offset_h", "actual_perim_ms", "actual_time_utc",
               "actual_time_basis", "actual_offset_h", "pixel_m", "baseline_acres",
               "pred_total_acres", "pred_new_acres", "act_new_acres", "inter_acres",
               "precision", "recall", "iou", "act_outside_acres", "act_offgrid_frac",
               "hs_n", "hs_new_acres", "hs_hit_acres", "hs_capture",
               "expired_run", "notes"]

# satellite active-fire footprint radii (m): VIIRS 375m px, MODIS 1km px
HS_RADIUS_VIIRS = 187.5
HS_RADIUS_MODIS = 500.0


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

def sync_bucket(url):
    """Mirror the R2 bucket (public r2.dev URL) into raw/, manifest-driven."""
    import requests
    url = url.rstrip("/")
    sess = requests.Session()

    def get(key):
        r = sess.get(f"{url}/{key}", timeout=(10, 120))
        r.raise_for_status()
        return r

    n = 0
    for key in ("manifest.json", "fire_matches.json"):
        p = RAW / key
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(get(key).content)
    manifest = load_json(RAW / "manifest.json")
    wanted = []
    for run_key, entry in manifest["runs"].items():
        for pct, info in entry["files"].items():
            if info.get("ok"):
                wanted.append(f"forecast_archive/{run_key}/{pct}.tif")
    for slug, epochs in (manifest.get("perimeters") or {}).items():
        # index.json is MUTABLE (the perimeter timeline grows) — always refresh,
        # or a cached copy hides newer snapshots from PerimeterSet
        idx = RAW / "perimeter_archive" / slug / "index.json"
        if idx.exists():
            idx.unlink()
        wanted.append(f"perimeter_archive/{slug}/index.json")
        wanted += [f"perimeter_archive/{slug}/{e}.geojson" for e in epochs]
    for slug, days in (manifest.get("hotspots") or {}).items():
        # recent day files are rewritten upstream; re-fetch the last 3 always
        recent = set(days[-3:])
        for d in days:
            p = RAW / "hotspot_archive" / slug / f"{d}.geojson"
            if d in recent and p.exists():
                p.unlink()
            wanted.append(f"hotspot_archive/{slug}/{d}.geojson")
    for key in wanted:
        p = RAW / key
        if p.exists():
            continue
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(p.suffix + ".part")
        tmp.write_bytes(get(key).content)
        tmp.replace(p)
        n += 1
        time.sleep(0.02)
    log(f"sync: {n} new objects ({len(wanted)} tracked)")


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
        self.unreadable = 0
        by_capture = {}
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
    either falls inside (run, run+H] or it doesn't. Footprints are rasterized
    at instrument resolution (VIIRS 375m px -> r=187.5m, MODIS 1km -> r=500m)."""

    def __init__(self, slug):
        self.pts = []  # (dt, lon, lat, radius_m)
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
                     else HS_RADIUS_VIIRS)
                self.pts.append((dt, c[0], c[1], r))
        self.pts.sort(key=lambda x: x[0])
        self._mask = {}

    def mask(self, t0, t1, crs, transform, shape_):
        from rasterio import features
        from rasterio.warp import transform as warp_transform
        from shapely.geometry import Point, mapping
        key = (t0, t1, str(crs), tuple(transform)[:6], shape_)
        if key not in self._mask:
            sel = [p for p in self.pts if t0 < p[0] <= t1]
            if not sel:
                self._mask[key] = (np.zeros(shape_, dtype=bool), 0)
            else:
                xs, ys = warp_transform("EPSG:4326", crs,
                                        [p[1] for p in sel], [p[2] for p in sel])
                geoms = [(mapping(Point(x, y).buffer(p[3], resolution=8)), 1)
                         for x, y, p in zip(xs, ys, sel)]
                m = features.rasterize(geoms, out_shape=shape_, transform=transform,
                                       fill=0, dtype="uint8").astype(bool)
                self._mask[key] = (m, len(sel))
        return self._mask[key]


# ---------------------------------------------------------------- metrics

def acres(npix, px_area_sqm):
    return round(npix * px_area_sqm / SQM_PER_ACRE, 1)


def compute(only_slug=None):
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
    skipped = {"no_match": 0, "no_perims": 0, "no_horizon": 0, "missing_tif": 0,
               "unreadable_perims": 0}

    for run_key in sorted(manifest["runs"]):
        entry = manifest["runs"][run_key]
        slug = entry["slug"]
        if only_slug and slug != only_slug:
            continue
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
            hssets[slug] = HotspotSet(slug)
        hs = hssets[slug]
        hs._mask.clear()  # keyed per grid+window; new run = new grid, drop old

        run_dt = datetime.strptime(entry["run_ts"], "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
        baseline = ps.latest_at_or_before(run_dt + timedelta(hours=BASELINE_TOL_H))
        # horizons that have ground truth close enough — cheap pre-check
        horizons = []
        for H in HORIZONS_H:
            actual = ps.nearest(run_dt + timedelta(hours=H), ACTUAL_TOL_H)
            if actual:
                horizons.append((H, actual))
        if not horizons:
            skipped["no_horizon"] += 1
            continue

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

            for H, (a_dt, a_ems, a_off, a_basis) in horizons:
                A = ps.mask(a_ems, crs, tfm, arr.shape)
                pred = valid & (toa_h <= H)
                pred_new = pred & ~B
                act_new = A & ~B
                inter = pred_new & act_new
                n_pred, n_act, n_int = int(pred_new.sum()), int(act_new.sum()), int(inter.sum())
                n_union = n_pred + n_act - n_int
                hs_mask, hs_n = hs.mask(run_dt, run_dt + timedelta(hours=H), crs, tfm, arr.shape)
                hs_new = hs_mask & ~B
                n_hs_new = int(hs_new.sum())
                n_hs_hit = int((pred_new & hs_new).sum())

                a_geom = ps.geom_utm(a_ems, crs)
                act_new_geom = a_geom.difference(b_geom) if b_geom else a_geom
                outside = act_new_geom.difference(grid_box)
                offgrid_frac = (outside.area / act_new_geom.area
                                if act_new_geom.area > 0 else 0.0)

                notes = [] if baseline else ["baseline=none"]
                rows.append({
                    "slug": slug, "cornea_id": m["cornea_id"],
                    "run_ts": entry["run_ts"],
                    "run_dt_utc": run_dt.strftime(ISO_Z),
                    "percentile": pct, "horizon_h": H,
                    "valid_dt_utc": (run_dt + timedelta(hours=H)).strftime(ISO_Z),
                    "baseline_perim_ms": baseline[1] if baseline else "",
                    "baseline_time_utc": baseline[0].strftime(ISO_Z) if baseline else "",
                    "baseline_time_basis": baseline[2] if baseline else "",
                    "baseline_offset_h": (round((baseline[0] - run_dt).total_seconds() / 3600, 2)
                                          if baseline else ""),
                    "actual_perim_ms": a_ems,
                    "actual_time_utc": a_dt.strftime(ISO_Z),
                    "actual_time_basis": a_basis,
                    "actual_offset_h": round(a_off, 2),
                    "pixel_m": round(abs(tfm.a), 1),
                    "baseline_acres": acres(int(B.sum()), px_area),
                    "pred_total_acres": acres(int(pred.sum()), px_area),
                    "pred_new_acres": acres(n_pred, px_area),
                    "act_new_acres": acres(n_act, px_area),
                    "inter_acres": acres(n_int, px_area),
                    "precision": round(n_int / n_pred, 4) if n_pred else "",
                    "recall": round(n_int / n_act, 4) if n_act else "",
                    "iou": round(n_int / n_union, 4) if n_union else "",
                    "act_outside_acres": round(outside.area / SQM_PER_ACRE, 1),
                    "act_offgrid_frac": round(offgrid_frac, 4),
                    "hs_n": hs_n,
                    "hs_new_acres": acres(n_hs_new, px_area),
                    "hs_hit_acres": acres(n_hs_hit, px_area),
                    "hs_capture": round(n_hs_hit / n_hs_new, 4) if n_hs_new else "",
                    "expired_run": entry.get("expired", False),
                    "notes": ";".join(notes),
                })
        log(f"{run_key}: {len(horizons)} horizons x pcts done")

    out = REPO / "data" / "metrics.csv"
    out.parent.mkdir(exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {len(rows)} rows -> {out} | skipped runs: {skipped}")


# ---------------------------------------------------------------- skill files

def _num(s):
    return float(s) if s else None


def _r3(s):
    v = _num(s)
    return None if v is None else round(v, 3)


def _passes_skill_filters(r):
    # blank offset / off-grid count as 0, as report.html's null coercion does
    grow = _num(r["act_new_acres"])
    off = _num(r["actual_offset_h"]) or 0.0
    offgrid = _num(r["act_offgrid_frac"]) or 0.0
    return (grow is not None and grow >= SKILL_FILTERS["min_act_new_acres"]
            and abs(off) <= SKILL_FILTERS["max_abs_offset_h"]
            and offgrid <= SKILL_FILTERS["max_offgrid_frac"])


def build_skill(csv_path, out_dir):
    """Per-fire skill files for map.html from a metrics CSV: out_dir/<slug>.json
    + out_dir/index.json, rows filtered like the report's default view.
    Returns (fires, rows)."""
    fires = {}
    with open(csv_path, newline="") as fh:
        for r in csv.DictReader(fh):
            if not _passes_skill_filters(r):
                continue
            entry = [int(float(r["horizon_h"])),
                     _r3(r["precision"]), _r3(r["recall"]), _r3(r["iou"])]
            (fires.setdefault(r["slug"], {}).setdefault(r["run_ts"], {})
             .setdefault(r["percentile"], []).append(entry))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for p in out_dir.glob("*.json"):
        p.unlink()
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    index = {}
    for slug in sorted(fires):
        runs = {run_ts: {pct: sorted(by_pct[pct], key=lambda e: e[0])
                         for pct in sorted(by_pct, key=int)}
                for run_ts, by_pct in sorted(fires[slug].items())}
        index[slug] = sum(len(v) for by_pct in runs.values() for v in by_pct.values())
        doc = {"v": 1, "slug": slug, "generated": generated,
               "filters": SKILL_FILTERS, "runs": runs}
        (out_dir / f"{slug}.json").write_text(json.dumps(doc, separators=(",", ":")))
    (out_dir / "index.json").write_text(json.dumps(
        {"v": 1, "generated": generated, "slugs": index}, separators=(",", ":")))
    return len(index), sum(index.values())


def publish_skill(client, bucket, csv_path):
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
    r = requests.get(f"{PUBLIC_BUCKET_URL}/metrics.csv", timeout=(10, 120))
    r.raise_for_status()
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
    ap.add_argument("--sync", action="store_true",
                    help=f"mirror the R2 bucket ({PUBLIC_BUCKET_URL}) into raw/ first")
    ap.add_argument("--sync-url", help="override the public bucket URL for --sync")
    ap.add_argument("--refresh-perims", action="store_true",
                    help="re-fetch Cornea perimeter snapshots before computing")
    ap.add_argument("--inspect", metavar="TIF")
    ap.add_argument("--overlay", nargs=4, metavar=("SLUG", "RUN", "PCT", "H"))
    ap.add_argument("--only", help="restrict metrics to one slug")
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
    if args.sync or args.sync_url:
        sync_bucket(args.sync_url or PUBLIC_BUCKET_URL)
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
    compute(only_slug=args.only)
    if args.publish:
        import collect
        client, bucket = collect.make_r2()
        client.upload_file(str(REPO / "data" / "metrics.csv"), bucket, "metrics.csv",
                           ExtraArgs={"ContentType": "text/csv"})
        log("published metrics.csv to bucket")
        publish_skill(client, bucket, REPO / "data" / "metrics.csv")


if __name__ == "__main__":
    main()
