#!/usr/bin/env python3
"""Publish a FARSITE hindcast version to the archive bucket, in the dashboard's format.

The FARSITE pilot (branch farsite-hindcast, farsite/) writes each run's ensemble as
data/fires/<slug>/runs-<variant>/<run_ts>/toa_p{10..90}.tif: arrival in hours since the run,
0 = no data, 1e-3 = burning at the start, in EPSG:5070 on the fire's whole landscape grid.
The dashboard reads forecasts like Pyrecast's archive copies: per-fire UTM GeoTIFFs. This
reprojects each grid onto the same run's Pyrecast grid (its 50th-percentile tif in the
archive: CRS, transform and size; nearest neighbour, so arrival hours are kept as they
were), so both models are scored pixel for pixel on the same ground and the map redraws
exactly what was scored. Without a Pyrecast tif, the run is cropped to where it burned
(plus a margin) on a 30 m grid in the fire's UTM zone. Writes

  farsite_archive/<slug>/<run_ts>/<pct>.tif   float32 hours since the run, nodata 0
  farsite_archive/<slug>/<run_ts>/ign.geojson  the lit edge (EPSG:4326): part of every model's
                                               common start in analyze.py's scoring

What burned before the run is left out of the grids, as in Pyrecast's (its start is nodata):
FARSITE's own start (the cells it lit at the start, its grid's 1e-3 hours) becomes nodata,
and ign.geojson becomes the lit edge plus those cells (the ignition as FARSITE rasterized it,
a cell or so wider), so all of it is in the common start. Arrival inside the baseline
perimeter (outside its unburned islands) is nodata too: the run burns the scar as
non-burnable except a 30 m strip around the lit edge (hindcast.write_inputs), which it can
re-burn, and nothing there can burn again. Neither changes growth beyond the common start.
  farsite_archive/manifest.json               {"model", "variant", "label", "horizon_h",
                                               "fires": {slug: {"runs": {run_ts: {...}}}}}

The manifest is merged with the published one: publishing one fire keeps the others, and a
fire's runs are added to (or replace) the ones already published, so a re-run of a few runs
publishes on its own; --replace drops the fire's published runs first.

  farsite_publish.py RESULTS_DIR --slug ca-dome [--variant cal] [--replace] (--out DIR | --upload)

RESULTS_DIR is searched for */fires/<slug>/runs-<variant>/<run_ts>/meta.json, so a
merged set of the workflow's shard-* artifacts (or its farsite-results) works as is.
"""
import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PUBLIC_BUCKET_URL = "https://f005.backblazeb2.com/file/fire-forecast-archive"
PREFIX = "farsite_archive"
PCTS = ("10", "30", "50", "70", "90")
LABELS = {"cal": "FARSITE (calibrated)", "calr": "FARSITE (calibrated, easing)",
          "calbar": "FARSITE (calibrated + control lines)", "base": "FARSITE (uncalibrated)"}
HORIZON_H = 96          # the pilot's forecast length (farsite/hindcast.py HORIZON_H)
MARGIN_M = 2000.0
RES_M = 30.0
START_H = 0.01          # FARSITE's start in its grids (1e-3 h, hindcast.py)
FORMAT = 2              # 2: start and re-burned scar left out (see above); in each run's entry


def utm_epsg(lon, lat):
    return (32600 if lat >= 0 else 32700) + int(math.floor((lon + 180) / 6)) + 1


def pyrecast_grid(slug, run_ts, cache):
    """(crs, transform, width, height) of the run's Pyrecast 50th-percentile tif, or None"""
    import requests
    import rasterio
    p = cache / slug / f"{run_ts}_50.tif"
    if not p.exists():
        r = requests.get(f"{PUBLIC_BUCKET_URL}/forecast_archive/{slug}/{run_ts}/50.tif", timeout=120)
        if r.status_code != 200:
            return None
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(r.content)
    with rasterio.open(p) as ds:
        return ds.crs, ds.transform, ds.width, ds.height


_PERIMS = {}


def baseline_perimeter(slug, ms):
    """the run's baseline perimeter (GeoJSON geometry, EPSG:4326) from the archive, or None"""
    import requests
    key = (slug, str(ms))
    if key not in _PERIMS:
        r = requests.get(f"{PUBLIC_BUCKET_URL}/perimeter_archive/{slug}/{ms}.geojson", timeout=120)
        gj = r.json() if r.status_code == 200 else None
        feats = (gj.get("features") if gj and "features" in gj else [gj]) if gj else []
        geoms = [f["geometry"] for f in feats if f and f.get("geometry")]
        _PERIMS[key] = ({"type": "GeometryCollection", "geometries": geoms} if len(geoms) > 1
                        else geoms[0] if geoms else None)
    return _PERIMS[key]


def convert_run(rdir, out_dir, grid=None, slug=None):
    """one run's five grids + lit edge -> out_dir (on `grid` when given, see pyrecast_grid),
    less what burned before the run (see the module notes); -> manifest entry, or None if
    unusable"""
    import rasterio
    from rasterio import features
    from rasterio.warp import Resampling, reproject, transform_geom
    from shapely.geometry import mapping, shape
    from shapely.ops import unary_union
    meta = json.loads((rdir / "meta.json").read_text())
    if meta.get("status") != "ok":
        return None
    ign = json.loads((rdir / "ign.geojson").read_text()) if (rdir / "ign.geojson").exists() else None
    with rasterio.open(rdir / "toa_p50.tif") as ds:
        src_crs, src_tfm, shape_ = ds.crs, ds.transform, ds.shape
    # one window for all five: where any percentile burned, plus the lit edge, plus a margin
    any_ = np.zeros(shape_, bool)
    grids = {}
    for p in PCTS:
        with rasterio.open(rdir / f"toa_p{p}.tif") as ds:
            a = ds.read(1).astype("float32")
        grids[p] = a
        any_ |= a > 0
    if not any_.any():
        return None
    rows, cols = np.nonzero(any_)
    x0, y1 = src_tfm * (cols.min(), rows.min())
    x1, y0 = src_tfm * (cols.max() + 1, rows.max() + 1)
    if ign:
        g = unary_union([shape(f["geometry"]) for f in ign["features"]])
        bx = g.bounds
        x0, y0, x1, y1 = min(x0, bx[0]), min(y0, bx[1]), max(x1, bx[2]), max(y1, bx[3])
    x0, y0, x1, y1 = x0 - MARGIN_M, y0 - MARGIN_M, x1 + MARGIN_M, y1 + MARGIN_M
    from rasterio.warp import transform as warp_transform
    lon, lat = warp_transform(src_crs, "EPSG:4326", [(x0 + x1) / 2], [(y0 + y1) / 2])
    if grid:
        dst_crs, dst_tfm, w, h = grid
    else:
        dst_crs = f"EPSG:{utm_epsg(lon[0], lat[0])}"
        dst_tfm, w, h = _target(src_crs, dst_crs, x0, y0, x1, y1)
    out_dir.mkdir(parents=True, exist_ok=True)
    prof = dict(driver="GTiff", width=w, height=h, count=1, dtype="float32", crs=dst_crs,
                transform=dst_tfm, nodata=0.0, compress="deflate", predictor=3, tiled=True,
                blockxsize=256, blockysize=256)
    dsts = {}
    for p in PCTS:
        dst = np.zeros((h, w), "float32")
        reproject(grids[p], dst, src_transform=src_tfm, src_crs=src_crs, src_nodata=0.0,
                  dst_transform=dst_tfm, dst_crs=dst_crs, dst_nodata=0.0, resampling=Resampling.nearest)
        dsts[p] = dst
    # what burned before the run: FARSITE's start (every member lights the same cells), and
    # the baseline perimeter's burned ground (rasterized with its holes, the unburned islands)
    start = np.zeros((h, w), bool)
    for dst in dsts.values():
        start |= (dst > 0) & (dst <= START_H)
    perim = baseline_perimeter(slug, meta["baseline_ms"]) if slug and meta.get("baseline_ms") else None
    scar = (features.rasterize([(transform_geom("EPSG:4326", dst_crs, perim), 1)], out_shape=(h, w),
                               transform=dst_tfm, fill=0, dtype="uint8").astype(bool)
            if perim else np.zeros((h, w), bool))
    px_ac = abs(dst_tfm.a * dst_tfm.e) / 4046.8564224
    reburn = 0
    for p in PCTS:
        dst = dsts[p]
        if p == "50":
            reburn = int((scar & ~start & (dst > 0)).sum())
        dst[start | scar] = 0.0
        with rasterio.open(out_dir / f"{p}.tif", "w", **prof) as ds:
            ds.write(dst, 1)
    # the lit edge, plus the cells FARSITE lit (in the grid's CRS, then EPSG:4326)
    polys = [transform_geom(src_crs, dst_crs, f["geometry"]) for f in (ign or {}).get("features", [])]
    polys += [g for g, v in features.shapes(start.astype("uint8"), mask=start, transform=dst_tfm) if v == 1]
    if polys:
        lit = unary_union([shape(g) for g in polys]).buffer(0)
        (out_dir / "ign.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": [
            {"type": "Feature", "properties": {}, "geometry": transform_geom(dst_crs, "EPSG:4326", mapping(lit))}]}))
    meta["start_ac"] = round(float(start.sum()) * px_ac)
    meta["reburn_ac_p50"] = round(reburn * px_ac)
    meta["scar_masked"] = bool(perim)
    keep = ("status", "members", "res_m", "ignition_ac", "baseline_time", "baseline_ms", "cal_status",
            "cal_factor", "cal_obs_ac", "cal_pred_ac", "cal_window_h", "central_new_ac", "new_ac_min",
            "new_ac_max", "barrier_map", "barrier_km", "member_timeout_s", "coarse_retries",
            "start_ac", "reburn_ac_p50", "scar_masked")
    (out_dir / "meta.json").write_text(json.dumps({k: meta[k] for k in keep if k in meta}, separators=(",", ":")))
    return {"pcts": list(PCTS), "ign": bool(polys), "crs": str(dst_crs), "pyrecast_grid": bool(grid),
            "fmt": FORMAT, **{k: meta[k] for k in ("cal_factor", "cal_status", "central_new_ac") if k in meta}}


def _target(src_crs, dst_crs, x0, y0, x1, y1):
    """a north-up 30 m grid in dst_crs covering the source window"""
    from rasterio.transform import from_origin
    from rasterio.warp import transform as warp_transform
    xs, ys = warp_transform(src_crs, dst_crs, [x0, x1, x0, x1], [y0, y0, y1, y1])
    left, right = math.floor(min(xs) / RES_M) * RES_M, math.ceil(max(xs) / RES_M) * RES_M
    bottom, top = math.floor(min(ys) / RES_M) * RES_M, math.ceil(max(ys) / RES_M) * RES_M
    return from_origin(left, top, RES_M, RES_M), int(round((right - left) / RES_M)), int(round((top - bottom) / RES_M))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", type=Path)
    ap.add_argument("--slug", required=True, action="append", help="fire(s) to publish")
    ap.add_argument("--variant", default="cal", choices=sorted(LABELS))
    ap.add_argument("--replace", action="store_true",
                    help="publish only these results' runs for each fire (default: add to the published ones)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out", type=Path, help="write the bucket layout under this directory instead")
    g.add_argument("--upload", action="store_true", help="upload to the bucket (needs S3_* env)")
    a = ap.parse_args()

    import requests
    stage = a.out or Path("farsite_stage")
    r = requests.get(f"{PUBLIC_BUCKET_URL}/{PREFIX}/manifest.json", timeout=60)
    if r.status_code == 404:
        manifest = {}     # the first publish
    else:
        r.raise_for_status()  # never merge into (and overwrite) a manifest we couldn't read
        manifest = r.json()
    # only what this script writes (a first publish once merged the bucket's 404 body in)
    manifest = {k: manifest[k] for k in ("variant", "fires") if k in manifest}
    if manifest.get("variant") not in (None, a.variant):
        sys.exit(f"the published FARSITE version is {manifest.get('variant')}, not {a.variant}: refusing to mix them")
    manifest.update({"v": 1, "model": "farsite", "variant": a.variant, "label": LABELS[a.variant],
                     "horizon_h": HORIZON_H, "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")})
    fires = manifest.setdefault("fires", {})
    written = []
    for slug in a.slug:
        runs = {}
        for m in sorted(a.results.rglob(f"fires/{slug}/runs-{a.variant}/*/meta.json")):
            rdir, run_ts = m.parent, m.parent.name
            if run_ts in runs:
                continue
            entry = convert_run(rdir, stage / PREFIX / slug / run_ts,
                                pyrecast_grid(slug, run_ts, stage / "_pyrecast"), slug)
            if entry:
                runs[run_ts] = entry
                written += [f"{PREFIX}/{slug}/{run_ts}/{n}" for n in [f"{p}.tif" for p in PCTS] + ["meta.json"] +
                            (["ign.geojson"] if entry["ign"] else [])]
        kept = {} if a.replace else {k: v for k, v in ((fires.get(slug) or {}).get("runs") or {}).items()
                                     if k not in runs}
        fires[slug] = {"runs": dict(sorted({**kept, **runs}.items()))}
        print(f"{slug}: {len(runs)} runs converted ({a.variant}), {len(kept)} already published kept, "
              f"{len(fires[slug]['runs'])} in all", flush=True)
    (stage / PREFIX).mkdir(parents=True, exist_ok=True)
    (stage / PREFIX / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))
    if a.upload:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import collect
        client, bucket = collect.make_r2()
        types = {".tif": "image/tiff", ".json": "application/json", ".geojson": "application/geo+json"}
        import concurrent.futures
        def put(key):
            client.upload_file(str(stage / key), bucket, key,
                               ExtraArgs={"ContentType": types[Path(key).suffix]})
        with concurrent.futures.ThreadPoolExecutor(16) as ex:
            list(ex.map(put, written))
        put(f"{PREFIX}/manifest.json")    # last: it never names a file that isn't up yet
        print(f"uploaded {len(written) + 1} files", flush=True)


if __name__ == "__main__":
    main()
