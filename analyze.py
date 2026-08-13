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
"""
import argparse
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
BASELINE_TOL_H = 3.0     # baseline = latest perimeter <= run + this
ACTUAL_TOL_H = 12.0      # actual must be within +/- this of run + H
SQM_PER_ACRE = 4046.8564224
PERCENTILES = ("10", "30", "50", "70", "90")

REPO = Path(__file__).resolve().parent
RAW = REPO / "raw"

CSV_COLUMNS = ["slug", "cornea_id", "run_ts", "run_dt_utc", "percentile",
               "horizon_h", "valid_dt_utc", "baseline_perim_ms", "baseline_offset_h",
               "actual_perim_ms", "actual_offset_h", "pixel_m", "baseline_acres",
               "pred_total_acres", "pred_new_acres", "act_new_acres", "inter_acres",
               "precision", "recall", "iou", "act_outside_acres", "act_offgrid_frac",
               "expired_run", "notes"]


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
        wanted.append(f"perimeter_archive/{slug}/index.json")
        wanted += [f"perimeter_archive/{slug}/{e}.geojson" for e in epochs]
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

class PerimeterSet:
    """Perimeter timeline for one fire + lazy per-grid rasterization cache."""

    def __init__(self, slug):
        self.slug = slug
        self.timeline = []  # sorted [(dt, epochms)]
        idx = load_json(RAW / "perimeter_archive" / slug / "index.json")
        for item in (idx or {}).get("index", []):
            path, date = item.get("path"), item.get("date")
            if not path or not date:
                continue
            epochms = path.rsplit("/", 1)[-1]
            if (RAW / "perimeter_archive" / slug / f"{epochms}.geojson").exists():
                self.timeline.append((parse_iso(date), epochms))
        self.timeline.sort()
        self._geom_utm = {}   # (epochms, crs_str) -> shapely geom in raster CRS
        self._mask = {}       # (epochms, grid_sig) -> bool ndarray

    def latest_at_or_before(self, dt):
        best = None
        for t, e in self.timeline:
            if t <= dt:
                best = (t, e)
        return best

    def nearest(self, dt, tol_h):
        best = None
        for t, e in self.timeline:
            off = (t - dt).total_seconds() / 3600.0
            if abs(off) <= tol_h and (best is None or abs(off) < abs(best[2])):
                best = (t, e, off)
        return best

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
    rows = []
    skipped = {"no_match": 0, "no_perims": 0, "no_horizon": 0, "missing_tif": 0}

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
            perimsets[slug] = PerimeterSet(slug)
        ps = perimsets[slug]
        if not ps.timeline:
            skipped["no_perims"] += 1
            continue

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

            for H, (a_dt, a_ems, a_off) in horizons:
                A = ps.mask(a_ems, crs, tfm, arr.shape)
                pred = valid & (toa_h <= H)
                pred_new = pred & ~B
                act_new = A & ~B
                inter = pred_new & act_new
                n_pred, n_act, n_int = int(pred_new.sum()), int(act_new.sum()), int(inter.sum())
                n_union = n_pred + n_act - n_int

                a_geom = ps.geom_utm(a_ems, crs)
                act_new_geom = a_geom.difference(b_geom) if b_geom else a_geom
                outside = act_new_geom.difference(grid_box)
                offgrid_frac = (outside.area / act_new_geom.area
                                if act_new_geom.area > 0 else 0.0)

                notes = [] if baseline else ["baseline=none"]
                rows.append({
                    "slug": slug, "cornea_id": m["cornea_id"],
                    "run_ts": entry["run_ts"],
                    "run_dt_utc": run_dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "percentile": pct, "horizon_h": H,
                    "valid_dt_utc": (run_dt + timedelta(hours=H)).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "baseline_perim_ms": baseline[1] if baseline else "",
                    "baseline_offset_h": (round((baseline[0] - run_dt).total_seconds() / 3600, 2)
                                          if baseline else ""),
                    "actual_perim_ms": a_ems,
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
    args = ap.parse_args()

    if args.inspect:
        inspect(args.inspect)
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


if __name__ == "__main__":
    main()
