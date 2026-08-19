#!/usr/bin/env python3
"""Archive pyrecast ELMFIRE fire-spread forecasts + Cornea perimeter snapshots.

Pyrecast deletes forecast runs after ~1 day, so this must run every ~6h.
Default mode uploads to Cloudflare R2 (GitHub Actions); --local archives to
raw/ instead. --push uploads a local raw/ archive into R2 (do NOT run while
the CI workflow is active — the manifest assumes a single writer).

Layout mirrored into the archive:
  forecast_archive/{slug}/{run_ts}/{pct}.tif     pct in 10,30,50,70,90
  perimeter_archive/{slug}/index.json, {epochms}.geojson
  manifest.json                                  collector state
  fire_matches.json                              slug -> cornea fire match
"""
import argparse
import concurrent.futures
import io
import json
import math
import os
import re
import sys
import tarfile
import tempfile
import threading
import time
import urllib.parse
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests

PYRECAST_BASE = "https://data.pyrecast.org/fire_spread_forecast/"
GEOSERVER_OWS = "https://geoserver-usw1.pyrecast.org/geoserver02/ows"
CORNEA_BASE = "https://fire-api-prod.web.app"
PERCENTILES = ("10", "30", "50", "70", "90")
# hourly-mosaic variables; each granule dir is tarred into ONE archive object per
# (run, pct, var) so R2 write-ops stay in the free tier (169 files -> 1 PUT)
VAR_DIRS = ("crown-fire", "flame-length", "hours-since-burned", "spread-rate")
ISO_EXTS = (".shp", ".dbf", ".prj", ".shx", ".qix")
FETCH_WORKERS = 12
VAR_COMPLETE_H = 167  # granules reaching run+167h == variable fully published
RUN_RE = re.compile(r"^\d{8}_\d{6}/$")
HREF_RE = re.compile(r'href="([^"]+)"')
TIFF_MAGIC = (b"II*\x00", b"MM\x00*")
MATCH_NAME_DIST_KM = 50.0
MATCH_SPATIAL_ONLY_KM = 15.0

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "fire-forecast-eval/0.1 (research; contact: repo owner)"

counts = {}


def log(msg):
    print(f"[collect] {msg}", flush=True)


def bump(key, n=1):
    counts[key] = counts.get(key, 0) + n


def utcnow():
    return datetime.now(timezone.utc)


def iso(dt):
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def fetch(url, tries=2, timeout=(10, 120), sleep_s=1.0):
    """GET with retries; returns Response or None. Retries 429 with backoff."""
    last = None
    for attempt in range(tries):
        try:
            r = SESSION.get(url, timeout=timeout)
            if r.status_code == 429:
                last = "HTTP 429"
                time.sleep(max(sleep_s * (attempt + 1), 2.0))
                continue
            if r.status_code == 200:
                return r
            last = f"HTTP {r.status_code}"
            if r.status_code == 404:
                break  # not-yet-published; no point hammering
        except requests.RequestException as e:
            last = f"{type(e).__name__}: {e}"
        time.sleep(sleep_s)
    fetch.last_error = last
    return None


fetch.last_error = None


class Archive:
    """Local mirror at root; optionally backed by an R2 (S3-compatible) bucket.

    All keys are relative POSIX paths. The manifest is the source of truth for
    what is already archived — safe because CI serializes runs (concurrency
    group) and --push must not run concurrently with CI.
    """

    def __init__(self, root, r2=None):
        self.root = Path(root)
        self.r2 = r2  # (client, bucket) or None

    def local_path(self, key):
        p = self.root / key
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def load_json(self, key):
        if self.r2:
            client, bucket = self.r2
            try:
                body = client.get_object(Bucket=bucket, Key=key)["Body"].read()
                return json.loads(body)
            except client.exceptions.NoSuchKey:
                return None
        p = self.root / key
        if not p.exists():
            return None
        return json.loads(p.read_text())

    def save_json(self, key, obj):
        data = json.dumps(obj, separators=(",", ":"))
        p = self.local_path(key)
        tmp = p.with_suffix(p.suffix + ".part")
        tmp.write_text(data)
        os.replace(tmp, p)
        if self.r2:
            client, bucket = self.r2
            client.put_object(Bucket=bucket, Key=key, Body=data.encode(),
                              ContentType="application/json")

    def commit_file(self, key):
        """Upload an already-downloaded local file when R2-backed."""
        if self.r2:
            client, bucket = self.r2
            client.upload_file(str(self.root / key), bucket, key)


def make_r2():
    """Build the S3-compatible client. Prefers the provider-agnostic S3_* env
    (any S3-compatible store: Backblaze B2, R2, ...); falls back to the legacy
    R2_* names so old secrets keep working during migration."""
    import boto3  # only needed in bucket mode; local mode stays boto3-free
    env = os.environ.get
    endpoint = env("S3_ENDPOINT") or (
        f"https://{env('R2_ACCOUNT_ID')}.r2.cloudflarestorage.com"
        if env("R2_ACCOUNT_ID") else None)
    if endpoint and "://" not in endpoint:  # boto3 requires a scheme
        endpoint = f"https://{endpoint}"
    access = env("S3_ACCESS_KEY_ID") or env("R2_ACCESS_KEY_ID")
    secret = env("S3_SECRET_ACCESS_KEY") or env("R2_SECRET_ACCESS_KEY")
    bucket = env("S3_BUCKET") or env("R2_BUCKET") or "fire-forecast-archive"
    if not (endpoint and access and secret):
        log("FATAL missing env: need S3_ENDPOINT + S3_ACCESS_KEY_ID + "
            "S3_SECRET_ACCESS_KEY (or legacy R2_* equivalents)")
        sys.exit(1)
    client = boto3.client(
        "s3", endpoint_url=endpoint, aws_access_key_id=access,
        aws_secret_access_key=secret, region_name="auto",
    )
    return client, bucket


# ---------------------------------------------------------------- scraping

def list_hrefs(url):
    r = fetch(url)
    if r is None:
        return None
    return HREF_RE.findall(r.text)


def scrape_runs(only=None):
    """Return {slug: [run_ts, ...]} from the pyrecast autoindex."""
    hrefs = list_hrefs(PYRECAST_BASE)
    if hrefs is None:
        log(f"FATAL root listing failed: {fetch.last_error}")
        sys.exit(1)
    slugs = sorted(h.rstrip("/") for h in hrefs
                   if h.endswith("/") and not h.startswith((".", "/", "?")))
    if not slugs:
        # tripwire: autoindex format drift must fail loudly, not archive nothing
        log("FATAL root listing parsed to zero fire directories — format change?")
        sys.exit(1)
    if only:
        slugs = [s for s in slugs if s in only]
    out = {}
    for slug in slugs:
        hrefs = list_hrefs(f"{PYRECAST_BASE}{slug}/")
        if hrefs is None:
            log(f"{slug}: listing failed ({fetch.last_error}) — SKIPPING")
            bump("slug_list_errors")
            continue
        out[slug] = sorted(h.rstrip("/") for h in hrefs if RUN_RE.match(h))
        time.sleep(0.05)
    return out


def tif_centroid_lonlat(path):
    """(lon, lat) of the raster bounds center; None if rasterio unavailable."""
    try:
        import rasterio
        from rasterio.warp import transform as warp_transform
    except ImportError:
        return None
    try:
        with rasterio.open(path) as ds:
            cx = (ds.bounds.left + ds.bounds.right) / 2
            cy = (ds.bounds.bottom + ds.bounds.top) / 2
            xs, ys = warp_transform(ds.crs, "EPSG:4326", [cx], [cy])
            return (xs[0], ys[0])
    except Exception as e:
        log(f"centroid failed for {path}: {e}")
        return None


def fetch_wcs_toa(slug, run_ts, pct, dest):
    """Fallback for ToA tifs that 403 on the static host: pyrecast's files stage
    behind a public-ACL window (~5-6h, sometimes never opening), but their
    geoserver reads the same store internally. WCS output verified bit-identical
    to the static tif (2026-08-19). ToA only — every GetCoverage is server-side
    work for them, so this must never be used for the hourly granules."""
    layer = f"fire-spread-forecast_{slug}_{run_ts}__elmfire_landfire_{pct}_time-of-arrival"
    r = fetch(f"{GEOSERVER_OWS}?service=WCS&version=2.0.1&request=GetCoverage"
              f"&coverageId={layer}&format=image/geotiff", tries=1, timeout=(10, 180))
    if r is None or r.content[:4] not in TIFF_MAGIC:
        return None
    raw = dest.with_suffix(".wcs.part")
    raw.write_bytes(r.content)
    tmp = dest.with_suffix(".tif.part")
    try:
        # geoserver emits uncompressed (~25MB vs ~2MB) — recompress for the archive
        import rasterio
        with rasterio.open(raw) as ds:
            profile = ds.profile
            data = ds.read()
        profile.update(compress="deflate")
        with rasterio.open(tmp, "w", **profile) as out:
            out.write(data)
        os.replace(tmp, dest)
    except Exception as e:
        log(f"{slug}/{run_ts} pct{pct}: WCS recompress failed ({e}); archiving raw")
        os.replace(raw, dest)
    finally:
        raw.unlink(missing_ok=True)
        tmp.unlink(missing_ok=True)
    return dest.stat().st_size


def download_run(archive, manifest, slug, run_ts):
    """Fetch missing percentile tifs for one run. Returns True if entry changed."""
    run_key = f"{slug}/{run_ts}"
    entry = manifest["runs"].setdefault(run_key, {
        "slug": slug, "run_ts": run_ts, "first_seen": iso(utcnow()),
        "complete": False, "expired": False, "files": {}, "errors": {},
    })
    entry["last_seen"] = iso(utcnow())  # still listed upstream (geoserver still serves it)
    if entry["complete"]:
        return False
    changed = False
    for pct in PERCENTILES:
        if entry["files"].get(pct, {}).get("ok"):
            continue
        url = f"{PYRECAST_BASE}{slug}/{run_ts}/elmfire/landfire/{pct}/time-of-arrival.tif"
        key = f"forecast_archive/{slug}/{run_ts}/{pct}.tif"
        dest = archive.local_path(key)
        part = dest.with_suffix(".tif.part")
        r = fetch(url)
        if r is None:
            err = fetch.last_error or ""
            if "403" in err:  # staging-ACL window: same data is servable via WCS
                size = fetch_wcs_toa(slug, run_ts, pct, dest)
                if size:
                    archive.commit_file(key)
                    entry["files"][pct] = {"bytes": size, "etag": "", "ok": True, "via": "wcs"}
                    entry["errors"].pop(pct, None)
                    bump("tif_wcs_rescued")
                    changed = True
                    if "centroid" not in entry:
                        ll = tif_centroid_lonlat(dest)
                        if ll:
                            entry["centroid"] = [round(ll[0], 5), round(ll[1], 5)]
                    time.sleep(0.5)  # GetCoverage is server-side work; be gentle
                    continue
            entry["errors"][pct] = f"{err} at {iso(utcnow())}"
            bump("tif_errors")
            changed = True
            continue
        part.write_bytes(r.content)
        clen = r.headers.get("Content-Length")
        if not r.content[:4] in TIFF_MAGIC or (clen and int(clen) != len(r.content)):
            part.unlink(missing_ok=True)
            entry["errors"][pct] = f"bad content ({len(r.content)}b) at {iso(utcnow())}"
            bump("tif_errors")
            changed = True
            continue
        os.replace(part, dest)
        archive.commit_file(key)
        entry["files"][pct] = {"bytes": len(r.content),
                               "etag": r.headers.get("ETag", ""), "ok": True}
        entry["errors"].pop(pct, None)
        bump("tif_downloaded")
        changed = True
        if "centroid" not in entry:
            ll = tif_centroid_lonlat(dest)
            if ll:
                entry["centroid"] = [round(ll[0], 5), round(ll[1], 5)]
        time.sleep(0.1)
    entry["complete"] = all(entry["files"].get(p, {}).get("ok") for p in PERCENTILES)
    if entry["complete"]:
        bump("runs_completed")
    elif not entry["files"]:
        # zero successes: check the layout assumption before writing it off
        hrefs = list_hrefs(f"{PYRECAST_BASE}{slug}/{run_ts}/elmfire/") or []
        dirs = [h for h in hrefs if h.endswith("/") and h != "../"]
        if dirs and "landfire/" not in dirs:
            log(f"WARNING {run_key}: unexpected layout under elmfire/: {dirs}")
    return changed


_tls = threading.local()


def thread_session():
    if not hasattr(_tls, "s"):
        _tls.s = requests.Session()
        _tls.s.headers.update(SESSION.headers)
    return _tls.s


def fetch_content(url, tries=2):
    """Thread-safe GET returning bytes or None (for parallel granule fetches)."""
    for attempt in range(tries):
        try:
            r = thread_session().get(url, timeout=(10, 60))
            if r.status_code == 200:
                return r.content
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(0.5)
    return None


def write_tar(path, files):
    """Atomically write {name: bytes} as an uncompressed tar (tifs are already
    deflate-compressed internally, so gzip would buy little)."""
    tmp = path.with_suffix(path.suffix + ".part")
    with tarfile.open(tmp, "w") as tf:
        for name in sorted(files):
            info = tarfile.TarInfo(name)
            info.size = len(files[name])
            info.mtime = 0  # deterministic output
            tf.addfile(info, io.BytesIO(files[name]))
    os.replace(tmp, path)


def archive_run_variables(archive, manifest, slug, run_ts):
    """Tar + upload the hourly-mosaic variables and isochrones for one listed run.
    Retries while the run is listed; a variable is final once its granules reach
    the forecast horizon or the granule set is stable across two invocations."""
    entry = manifest["runs"][f"{slug}/{run_ts}"]
    vars_state = entry.setdefault("vars", {})
    run_dt = datetime.strptime(run_ts, "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
    base = f"{PYRECAST_BASE}{slug}/{run_ts}/elmfire/landfire"
    changed = False
    for pct in PERCENTILES:
        pstate = vars_state.setdefault(pct, {})
        if not pstate.get("isochrones", {}).get("ok"):
            files = {}
            for ext in ISO_EXTS:
                data = fetch_content(f"{base}/{pct}/isochrones{ext}")
                if data:
                    files[f"isochrones{ext}"] = data
                elif ext != ".qix":  # qix spatial index is optional
                    files = None
                    break
            if files:
                key = f"forecast_archive/{slug}/{run_ts}/{pct}_isochrones.tar"
                write_tar(archive.local_path(key), files)
                archive.commit_file(key)
                pstate["isochrones"] = {"ok": True, "n": len(files)}
                bump("iso_tars")
                changed = True
        for var in VAR_DIRS:
            vstate = pstate.get(var, {})
            if vstate.get("complete"):
                continue
            hrefs = list_hrefs(f"{base}/{pct}/{var}/")
            if hrefs is None:
                continue
            gran = sorted(h for h in hrefs
                          if re.fullmatch(rf"{re.escape(var)}_\d{{8}}_\d{{6}}\.tif", h))
            if not gran:
                continue
            last_ts = gran[-1][len(var) + 1:-4]
            last_dt = datetime.strptime(last_ts, "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
            complete = last_dt >= run_dt + timedelta(hours=VAR_COMPLETE_H)
            stable = (vstate.get("n") == len(gran) and vstate.get("last") == last_ts)
            if vstate.get("ok") and (complete or stable):
                vstate["complete"] = True  # existing tar already holds the final set
                pstate[var] = vstate
                changed = True
                continue
            contents = {}
            with concurrent.futures.ThreadPoolExecutor(FETCH_WORKERS) as ex:
                futs = {ex.submit(fetch_content, f"{base}/{pct}/{var}/{g}"): g for g in gran}
                for fut in concurrent.futures.as_completed(futs):
                    data = fut.result()
                    if data and data[:4] in TIFF_MAGIC:
                        contents[futs[fut]] = data
            if len(contents) < len(gran) * 0.9:
                # upstream mid-publish or flaky; try again next invocation
                bump("var_fetch_deferred")
                continue
            key = f"forecast_archive/{slug}/{run_ts}/{pct}_{var}.tar"
            write_tar(archive.local_path(key), contents)
            archive.commit_file(key)
            pstate[var] = {"ok": True, "n": len(gran), "got": len(contents),
                           "last": last_ts, "complete": complete}
            bump("var_tars")
            bump("var_granules", len(contents))
            changed = True
    return changed


def mark_expired(manifest, listed):
    listed_keys = {f"{s}/{r}" for s, runs in listed.items() for r in runs}
    for run_key, entry in manifest["runs"].items():
        if run_key not in listed_keys and not entry["complete"] and not entry["expired"]:
            entry["expired"] = True
            missing = [p for p in PERCENTILES if not entry["files"].get(p, {}).get("ok")]
            log(f"{run_key}: expired incomplete (missing pct {','.join(missing)})")
            bump("runs_expired")


# ---------------------------------------------------------------- matching

def norm_name(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = math.radians(lat2 - lat1), math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def parse_coords(s):
    try:
        lat, lon = (float(x) for x in s.split(","))
        return lat, lon
    except (AttributeError, ValueError):
        return None


def cornea_search(name):
    q = urllib.parse.urlencode({"search": name, "fire_type": "wildfire",
                                "active": "all", "limit": 250})
    r = fetch(f"{CORNEA_BASE}/fires?{q}", tries=3)
    if r is None:
        return None
    try:
        return r.json().get("fires", [])
    except ValueError:
        return None


def attempt_match(slug, centroid_lonlat):
    """Return (match_dict, None) or (None, unmatched_reason_dict)."""
    state, _, name = slug.partition("-")
    if not name:
        return None, {"reason": "unparseable slug"}
    query = name.replace("-", " ")
    tokens = name.split("-")
    # multi-word slugs may concatenate names (complexes) or add qualifiers, so
    # fall back to sub-token queries when the full name finds nothing in-state
    queries = [query]
    if len(tokens) > 1:
        queries += [" ".join(tokens[i:]) for i in range(1, len(tokens))]
        queries += [t for t in tokens if len(t) > 3]
    cands, seen = [], set()
    for q in queries:
        fires = cornea_search(q)
        if fires is None:
            if not cands:
                return None, {"reason": f"search failed: {fetch.last_error}"}
            continue
        for f in fires:
            if ((f.get("state") or "").upper() == state.upper()
                    and f.get("cornea_id") not in seen):
                seen.add(f["cornea_id"])
                cands.append(f)
        if cands and q == query:
            break  # full-name hits are enough; no need for fallback queries
        time.sleep(0.1)
    if not cands:
        return None, {"reason": "no candidates in state"}
    want = norm_name(query)
    scored = []
    for f in cands:
        title = norm_name(f.get("post_title"))
        name_ok = bool(title) and (title == want or want in title or title in want)
        dist = None
        coords = parse_coords(f.get("fire_coordinates"))
        if coords and centroid_lonlat:
            dist = haversine_km(coords[0], coords[1], centroid_lonlat[1], centroid_lonlat[0])
        scored.append((f, name_ok, dist))
    passing = [s for s in scored
               if s[1] and (s[2] is None or s[2] <= MATCH_NAME_DIST_KM)]
    # slug that concatenates >=2 distinct passing incident names is likely a
    # multi-fire complex — one fire's perimeter is the wrong ground truth
    titles = {norm_name(f.get("post_title")) for f, _, _ in passing}
    proper_subs = [t for t in titles if t and t != want and t in want]
    if len(proper_subs) >= 2:
        return None, {"reason": f"probable multi-fire complex: {sorted(proper_subs)}; "
                                "set overrides.json to force or skip"}
    if passing:
        if len(passing) > 1:
            passing.sort(key=lambda s: s[2] if s[2] is not None else 1e9)
            log(f"{slug}: {len(passing)} name matches, taking nearest")
        f, _, dist = passing[0]
        return {"cornea_id": f["cornea_id"], "post_title": f.get("post_title"),
                "method": "name+spatial" if dist is not None else "name_only",
                "dist_km": round(dist, 1) if dist is not None else None,
                "matched_at": iso(utcnow())}, None
    spatial = [s for s in scored if s[2] is not None and s[2] <= MATCH_SPATIAL_ONLY_KM]
    if spatial:
        spatial.sort(key=lambda s: s[2])
        f, _, dist = spatial[0]
        log(f"{slug}: WARNING spatial-only match -> {f.get('post_title')} ({dist:.1f}km)")
        return {"cornea_id": f["cornea_id"], "post_title": f.get("post_title"),
                "method": "spatial_only", "dist_km": round(dist, 1),
                "matched_at": iso(utcnow())}, None
    return None, {"reason": "no name or spatial pass",
                  "candidates": [{"post_title": f.get("post_title"),
                                  "dist_km": round(d, 1) if d is not None else None}
                                 for f, _, d in scored[:5]]}


def run_matching(archive, manifest, slugs, overrides):
    matches = archive.load_json("fire_matches.json") or {"matches": {}, "unmatched": {}}
    for slug in slugs:
        if slug in overrides:
            forced = overrides[slug]
            if forced is None:
                matches["matches"].pop(slug, None)
                matches["unmatched"][slug] = {"reason": "override: skip"}
            elif matches["matches"].get(slug, {}).get("cornea_id") != forced:
                matches["matches"][slug] = {"cornea_id": forced, "method": "override",
                                            "matched_at": iso(utcnow())}
                matches["unmatched"].pop(slug, None)
            continue
        if slug in matches["matches"]:
            continue  # sticky
        centroid = None
        for run_key, entry in manifest["runs"].items():
            if entry["slug"] == slug and entry.get("centroid"):
                centroid = entry["centroid"]
                break
        m, why = attempt_match(slug, centroid)
        if m:
            matches["matches"][slug] = m
            matches["unmatched"].pop(slug, None)
            log(f"{slug}: matched -> {m.get('post_title')} ({m['method']})")
            bump("fires_matched")
        else:
            why["last_tried"] = iso(utcnow())
            matches["unmatched"][slug] = why
            bump("fires_unmatched")
        time.sleep(0.1)
    matches["generated"] = iso(utcnow())
    archive.save_json("fire_matches.json", matches)
    return matches


# ---------------------------------------------------------------- perimeters

def snapshot_perimeters(archive, manifest, matches, slugs):
    archived = manifest.setdefault("perimeters", {})
    for slug in slugs:
        m = matches["matches"].get(slug)
        if not m:
            continue
        cid = urllib.parse.quote(m["cornea_id"], safe="")
        r = fetch(f"{CORNEA_BASE}/fires/{cid}/perimeters", tries=3)
        if r is None:
            log(f"{slug}: perimeter index failed ({fetch.last_error}) — SKIPPING")
            bump("perim_index_errors")
            continue
        try:
            idx = r.json()
        except ValueError:
            log(f"{slug}: perimeter index not JSON — SKIPPING")
            bump("perim_index_errors")
            continue
        if not isinstance(idx, list):
            idx = []
        archive.save_json(f"perimeter_archive/{slug}/index.json",
                          {"generated": iso(utcnow()), "cornea_id": m["cornea_id"],
                           "index": idx})
        have = set(archived.get(slug, []))
        for item in idx:
            path, date = item.get("path"), item.get("date")
            if not path:
                continue
            epochms = path.rsplit("/", 1)[-1]
            if epochms in have:
                continue
            rp = fetch(CORNEA_BASE + path, tries=3)
            if rp is None:
                bump("perim_errors")
                continue
            try:
                gj = rp.json()
            except ValueError:
                bump("perim_errors")
                continue
            key = f"perimeter_archive/{slug}/{epochms}.geojson"
            archive.save_json(key, gj)
            have.add(epochms)
            bump("perims_downloaded")
            time.sleep(0.1)
        archived[slug] = sorted(have)


# ---------------------------------------------------------------- hotspots

def _geom_bounds(gj):
    xs, ys = [], []

    def walk(c):
        if isinstance(c[0], (int, float)):
            xs.append(c[0])
            ys.append(c[1])
        else:
            for cc in c:
                walk(cc)
    geom = gj.get("geometry", gj)
    if geom.get("type") == "FeatureCollection":
        for f in geom["features"]:
            walk(f["geometry"]["coordinates"])
    else:
        walk(geom["coordinates"])
    return min(xs), min(ys), max(xs), max(ys)


def snapshot_hotspots(archive, manifest, matches, slugs):
    """Archive VIIRS/MODIS detections per (fire, acq_date). Detections are
    immutable once a day is past; we rewrite the last ~3 day files each run to
    catch late-arriving detections."""
    hs_state = manifest.setdefault("hotspots", {})
    since = (utcnow() - timedelta(days=3)).strftime("%Y-%m-%d")
    for slug in slugs:
        if slug not in matches["matches"]:
            continue
        bounds = None
        idx = archive.load_json(f"perimeter_archive/{slug}/index.json")
        if idx and idx.get("index"):
            latest = max(idx["index"], key=lambda p: p.get("date", ""))
            ems = latest["path"].rsplit("/", 1)[-1]
            gj = archive.load_json(f"perimeter_archive/{slug}/{ems}.geojson")
            if gj:
                w, s, e, n = _geom_bounds(gj)
                bounds = (s - 0.2, w - 0.2, n + 0.2, e + 0.2)  # lat-first for the API
        if bounds is None:
            cent = next((en.get("centroid") for en in manifest["runs"].values()
                         if en["slug"] == slug and en.get("centroid")), None)
            if not cent:
                continue
            bounds = (cent[1] - 0.4, cent[0] - 0.4, cent[1] + 0.4, cent[0] + 0.4)
        r = fetch(f"{CORNEA_BASE}/hotspots?bbox={bounds[0]},{bounds[1]},{bounds[2]},{bounds[3]}"
                  f"&since={since}&limit=50000", tries=3)
        if r is None:
            bump("hs_errors")
            continue
        try:
            gj = r.json()
        except ValueError:
            bump("hs_errors")
            continue
        byday = {}
        for f in gj.get("features", []):
            d = (f.get("properties") or {}).get("acq_date")
            if d:
                byday.setdefault(d, []).append(f)
        for d, feats in byday.items():
            archive.save_json(f"hotspot_archive/{slug}/{d}.geojson",
                              {"type": "FeatureCollection", "features": feats})
            bump("hs_day_files")
        bump("hs_detections", sum(len(v) for v in byday.values()))
        hs_state[slug] = sorted(set(hs_state.get(slug, [])) | set(byday))
        time.sleep(0.1)


# ---------------------------------------------------------------- push mode

def push_local_to_r2(local_root):
    """Upload a --local archive into R2 and merge manifests. Single writer only:
    do not run while the CI workflow is active."""
    r2 = make_r2()
    local = Archive(local_root)
    remote = Archive(local_root, r2=r2)  # same paths; uploads read local files
    lman = local.load_json("manifest.json")
    if not lman:
        log("FATAL no local manifest to push")
        sys.exit(1)
    rman = remote.load_json("manifest.json") or {"runs": {}, "perimeters": {}}
    for run_key, lentry in lman["runs"].items():
        rentry = rman["runs"].setdefault(run_key, dict(lentry, files={}, errors={}))
        for pct, info in lentry["files"].items():
            if info.get("ok") and not rentry["files"].get(pct, {}).get("ok"):
                key = f"forecast_archive/{run_key}/{pct}.tif"
                remote.commit_file(key)
                rentry["files"][pct] = info
                bump("pushed_tifs")
                if counts["pushed_tifs"] % 25 == 0:
                    log(f"pushed {counts['pushed_tifs']} tifs...")
        rentry["complete"] = all(rentry["files"].get(p, {}).get("ok") for p in PERCENTILES)
        if lentry.get("centroid") and not rentry.get("centroid"):
            rentry["centroid"] = lentry["centroid"]
        # variable/isochrone tars produced by --local runs
        for pct, lp in (lentry.get("vars") or {}).items():
            rp = rentry.setdefault("vars", {}).setdefault(pct, {})
            for name, linfo in lp.items():
                if not linfo.get("ok") or rp.get(name, {}).get("ok"):
                    continue
                suffix = "isochrones" if name == "isochrones" else name
                key = f"forecast_archive/{run_key}/{pct}_{suffix}.tar"
                if (local.root / key).exists():
                    remote.commit_file(key)
                    rp[name] = linfo
                    bump("pushed_tars")
    for slug, epochs in (lman.get("perimeters") or {}).items():
        have = set(rman["perimeters"].get(slug, []))
        for ems in epochs:
            if ems not in have:
                remote.commit_file(f"perimeter_archive/{slug}/{ems}.geojson")
                have.add(ems)
                bump("pushed_perims")
        rman["perimeters"][slug] = sorted(have)
        idx = local.load_json(f"perimeter_archive/{slug}/index.json")
        if idx:
            remote.save_json(f"perimeter_archive/{slug}/index.json", idx)
    for slug, days in (lman.get("hotspots") or {}).items():
        have = set(rman.setdefault("hotspots", {}).get(slug, []))
        for d in days:
            key = f"hotspot_archive/{slug}/{d}.geojson"
            if (local.root / key).exists():
                remote.commit_file(key)
                have.add(d)
                bump("pushed_hs_days")
        rman["hotspots"][slug] = sorted(have)
    lmatch = local.load_json("fire_matches.json")
    rmatch = remote.load_json("fire_matches.json")
    if lmatch and not rmatch:
        remote.save_json("fire_matches.json", lmatch)
    elif lmatch and rmatch:
        for slug, m in lmatch["matches"].items():
            rmatch["matches"].setdefault(slug, m)
        remote.save_json("fire_matches.json", rmatch)
    rman["generated"] = iso(utcnow())
    remote.save_json("manifest.json", rman)
    log(f"push done: {counts}")


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--local", action="store_true",
                    help="archive into ./raw instead of R2")
    ap.add_argument("--push", action="store_true",
                    help="upload an existing --local archive into R2, then exit")
    ap.add_argument("--only", action="append",
                    help="restrict to slug (repeatable, for debugging)")
    ap.add_argument("--skip-perimeters", action="store_true")
    ap.add_argument("--skip-vars", action="store_true",
                    help="skip the hourly-mosaic variable tars (time-of-arrival only)")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent
    if args.push:
        push_local_to_r2(repo / "raw")
        return

    if args.local:
        archive = Archive(repo / "raw")
    else:
        archive = Archive(Path(tempfile.mkdtemp(prefix="fce-")), r2=make_r2())

    started = time.time()
    manifest = archive.load_json("manifest.json") or {"runs": {}, "perimeters": {}}
    listed = scrape_runs(only=set(args.only) if args.only else None)
    bump("slugs_listed", len(listed))
    bump("runs_listed", sum(len(v) for v in listed.values()))

    for slug, runs in sorted(listed.items()):
        changed = False
        for run_ts in runs:
            changed = download_run(archive, manifest, slug, run_ts) or changed
            if not args.skip_vars:
                changed = archive_run_variables(archive, manifest, slug, run_ts) or changed
        if changed:
            manifest["generated"] = iso(utcnow())
            archive.save_json("manifest.json", manifest)  # checkpoint per slug

    if not args.only:
        mark_expired(manifest, listed)

    try:
        overrides = json.loads((repo / "overrides.json").read_text()).get("overrides", {})
    except (OSError, ValueError):
        overrides = {}
    # ground truth must keep accumulating for up to 14 days AFTER a fire's last
    # run (long scoring horizons), even once pyrecast delists it
    cutoff = utcnow() - timedelta(days=15)
    scoring = set(listed)
    for e in manifest["runs"].values():
        try:
            if datetime.strptime(e["run_ts"], "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc) >= cutoff:
                scoring.add(e["slug"])
        except ValueError:
            pass
    if args.only:
        scoring &= set(args.only)
    matches = run_matching(archive, manifest, sorted(scoring), overrides)

    if not args.skip_perimeters:
        snapshot_perimeters(archive, manifest, matches, sorted(scoring))
        snapshot_hotspots(archive, manifest, matches, sorted(scoring))

    manifest["generated"] = iso(utcnow())
    archive.save_json("manifest.json", manifest)
    n_matched = sum(1 for s in listed if s in matches["matches"])
    log(f"done in {time.time() - started:.0f}s: {len(listed)} fires "
        f"({n_matched} matched), {counts}")


if __name__ == "__main__":
    main()
