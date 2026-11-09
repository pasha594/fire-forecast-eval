"""Score the FARSITE hindcast versions and PyreCast on the same runs, with the dashboard's own
scorer, so every model is scored the same way (analyze.scored_perimeters + analyze.score_grid):

  - each run against every perimeter in its first 7 days, at the perimeter's own time (its
    exact lead; no window or tolerance): every perimeter with a mapped time, the latest one
    published on a day without one, and the ones the previous fixed-horizon method showed;
    the same perimeters for every model
  - growth beyond ONE common start for all models: the baseline perimeter B, plus PyreCast's
    start footprint S (analyze.RunStart), plus FARSITE's lit edge I (ign.geojson)
  - no filters: every perimeter is scored (no minimum growth, timing or off-grid rule)
  - hotspot capture over the same window (run start to the perimeter's time)

Writes runs/hindcast_scores.csv (one row per model, run, percentile, perimeter) and prints
summaries over rows matched across every model.  usage: score2.py [SLUG ...]
"""
import csv
import json
import sys
from collections import defaultdict

import numpy as np
import rasterio
from rasterio import features
from rasterio.warp import transform_geom
from shapely.geometry import box

from common import DATA, FIRES, RUNS, analyze, baseline, perimset, run_dt, runs_for

# FARSITE versions by output dir (see hindcast.py variants); only those present are scored
FARSITE_DIRS = {"farsite": "runs", "farsite-cal": "runs-cal", "farsite-calr": "runs-calr",
                "farsite-calbar": "runs-calbar", "farsite-calbarr": "runs-calbarr"}
SHORT = {"pyrecast": "P", "farsite": "F", "farsite-cal": "Fcal", "farsite-calr": "Fcalr",
         "farsite-calbar": "Fbar", "farsite-calbarr": "Fbarr"}


def models_for(slugs):
    return ["pyrecast"] + [m for m, d in FARSITE_DIRS.items() if any((FIRES / s / d).is_dir() for s in slugs)]


def tif_for(model, slug, run, pct):
    if model == "pyrecast":
        return DATA / "raw" / "forecast_archive" / slug / run / f"{pct}.tif"
    return FIRES / slug / FARSITE_DIRS[model] / run / f"toa_p{pct}.tif"


def farsite_run_dirs(slug, run):
    return [FIRES / slug / d / run for d in FARSITE_DIRS.values() if (FIRES / slug / d / run / "meta.json").exists()]


def ignition_mask(slug, run, crs, tfm, shape_):
    """FARSITE's lit edge (EPSG:5070 polygons, the same for every version) on a forecast grid"""
    p = next((d / "ign.geojson" for d in farsite_run_dirs(slug, run) if (d / "ign.geojson").exists()), None)
    if p is None:
        return np.zeros(shape_, bool)
    gs = [transform_geom("EPSG:5070", crs, f["geometry"]) for f in json.load(open(p))["features"]]
    if not gs:
        return np.zeros(shape_, bool)
    return features.rasterize([(g, 1) for g in gs], out_shape=shape_, transform=tfm, fill=0,
                              dtype="uint8").astype(bool)


def score_fire(slug, models):
    ps = perimset(slug)
    hs = analyze.HotspotSet(slug, ps)
    rows = []
    for run in runs_for(slug):
        if not farsite_run_dirs(slug, run):
            continue
        dt = run_dt(run)
        b = baseline(ps, dt)
        perims = analyze.scored_perimeters(ps, dt)
        if not perims:
            continue
        start = analyze.RunStart(f"{slug}/{run}", dt, hs)
        for model in models:
            for pct in analyze.PERCENTILES:
                p = tif_for(model, slug, run, pct)
                if not p.exists():
                    continue
                with rasterio.open(p) as ds:
                    arr = ds.read(1)
                    crs, tfm, nod, grid_box = ds.crs, ds.transform, ds.nodata, box(*ds.bounds)
                valid = (arr != nod) if nod is not None else np.isfinite(arr)
                B = ps.mask(b[1], crs, tfm, arr.shape) if b else np.zeros(arr.shape, bool)
                before = B | start.mask(crs, tfm, arr.shape) | ignition_mask(slug, run, crs, tfm, arr.shape)
                b_geom = ps.geom_utm(b[1], crs) if b else None
                for sc in analyze.score_grid(arr * analyze.TOA_TO_HOURS, valid, crs, tfm, grid_box, ps, hs,
                                             dt, perims, before, b_geom, before & ~B):
                    rows.append(dict(model=model, slug=slug, run_ts=run, percentile=pct,
                                     start_basis=start.basis, **sc))
        ps._mask.clear()
        hs._mask.clear()
    return rows


def matched(rows, models):
    """rows for the same (run, pct, perimeter) from every model in `models`"""
    by = defaultdict(dict)
    for r in rows:
        if r["model"] in models:
            by[(r["slug"], r["run_ts"], r["percentile"], r["actual_perim_ms"])][r["model"]] = r
    return [tuple(v[m] for m in models) for v in by.values() if len(v) == len(models)]


def summarize(rows, title, models):
    f = lambda x: float(x) if x not in ("", None) else None
    sh = [SHORT.get(m, m) for m in models]
    mm = matched(rows, models)
    print(f"\n=== {title}  ({', '.join(f'{a} = {m}' for a, m in zip(sh, models))})")
    print(f"{len(mm)} matched rows (every percentile), {len({x[0]['run_ts'] for x in mm})} runs; "
          f"by day after the run (24 h periods), means over rows with a score:")
    print(f"{'day':>4} {'n':>5} " + " ".join(f"{'IoU ' + a:>9}" for a in sh) + " " +
          " ".join(f"{'over ' + a:>9}" for a in sh) + " " + " ".join(f"{'capt ' + a:>9}" for a in sh))
    for d in sorted({int(x[0]["day"]) for x in mm}):
        g = [x for x in mm if int(x[0]["day"]) == d]
        iou = lambda i: np.mean([f(x[i]["iou"]) for x in g if f(x[i]["iou"]) is not None] or [np.nan])
        over = lambda i: sum(f(x[i]["pred_new_acres"]) for x in g) / max(sum(f(x[i]["act_new_acres"]) for x in g), 1e-9)
        cap = lambda i: np.mean([f(x[i]["hs_capture"]) for x in g if f(x[i]["hs_capture"]) is not None] or [np.nan])
        k = range(len(models))
        print(f"{d:4d} {len(g):5d} " + " ".join(f"{iou(i):9.3f}" for i in k) + " " +
              " ".join(f"{over(i):9.1f}" for i in k) + " " + " ".join(f"{cap(i):9.2f}" for i in k))


def main():
    slugs = sys.argv[1:] or ["ca-dome", "ca-plaskett"]
    models = models_for(slugs)
    rows = []
    for s in slugs:
        rows += score_fire(s, models)
    RUNS.mkdir(exist_ok=True)
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(RUNS / "hindcast_scores.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {RUNS / 'hindcast_scores.csv'}")
    for s in slugs:
        present = [m for m in models if any(r["model"] == m and r["slug"] == s for r in rows)]
        for ms in ([present] if len(present) <= 2 else [present[:2], present]):
            summarize([r for r in rows if r["slug"] == s], s, ms)


if __name__ == "__main__":
    main()
