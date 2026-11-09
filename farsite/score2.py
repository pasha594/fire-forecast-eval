"""Score the FARSITE hindcast vs PyreCast on the same runs, at the resolution the truth
allows. Both forecasts are continuous in time (arrival time in hours), so the limit is
the observations:

  perim   every perimeter in (run, run + 96 h], scored at its exact lead time
          (predicted = arrival <= lead). Growth-only precision / recall / IoU, as the
          dashboard computes them. truth_basis records whether that time is when the
          perimeter was mapped or only when it was published (an upper bound).
  hotspot every satellite overpass in (run, run + 96 h]: share of that pass's new
          detection footprints (outside the run-time perimeter) that the forecast had
          burned by the overpass time. Recall-like only, but several passes a day.
  dash    the dashboard's own fixed horizons (12/24/48/72/96 h, nearest perimeter within
          +/-12 h, default filters), for continuity with report.html.

Growth is measured beyond a COMMON start, the same for both models: the official baseline
perimeter B, plus PyreCast's own start footprint S (its hours-since-burned granule, else every
detection footprint up to run time, as analyze.py's RunStart does), plus FARSITE's ignition I
(the active edge it was lit from). Each model was told some of that area was already burning;
scoring against B alone hands each model free hits for what it was given at the start.
perim rows are also scored against B alone (baseline="perimeter") to show the effect.

Writes runs/hindcast_scores.csv (every row) and runs/hindcast_bins.csv (1-hour lead bins),
prints paired summaries.  usage: score2.py [SLUG ...]
"""
import csv
import sys
from collections import defaultdict
from datetime import timedelta

import numpy as np
import rasterio
from rasterio import features
from rasterio.warp import transform_geom
from shapely.geometry import Point, box, mapping, shape
from shapely.ops import unary_union

from common import DATA, FIRES, RUNS, analyze, baseline, perimset, run_dt, runs_for

MAX_LEAD_H = 96
DASH_H = (12, 24, 48, 72, 96)
F = analyze.SKILL_FILTERS
PASS_GAP_S = 20 * 60
MODELS = ("pyrecast", "farsite")


def tif_for(model, slug, run, pct):
    if model == "pyrecast":
        return DATA / "raw" / "forecast_archive" / slug / run / f"{pct}.tif"
    return FIRES / slug / "runs" / run / f"toa_p{pct}.tif"


def passes(hs, t0, t1):
    """Group detections in (t0, t1] into overpasses: same radius class, < 20 min apart."""
    pts = [p for p in hs.pts if t0 < p[0] <= t1]
    groups = []
    for p in pts:
        for g in groups:
            if g[0][3] == p[3] and abs((p[0] - g[-1][0]).total_seconds()) < PASS_GAP_S:
                g.append(p)
                break
        else:
            groups.append([p])
    return [(max(x[0] for x in g), g) for g in groups]


def ignition_mask(slug, run, crs, tfm, shape_):
    """FARSITE's ignition polygons (EPSG:5070) on a forecast grid"""
    p = FIRES / slug / "runs" / run / "ign.geojson"
    if not p.exists():
        return np.zeros(shape_, bool)
    import json
    gs = [transform_geom("EPSG:5070", crs, f["geometry"]) for f in json.load(open(p))["features"]]
    if not gs:
        return np.zeros(shape_, bool)
    return features.rasterize([(g, 1) for g in gs], out_shape=shape_, transform=tfm, fill=0,
                              dtype="uint8").astype(bool)


def growth_row(arr, valid, B, bg, A, ag, gb, lead, px):
    pn = valid & (arr <= lead) & ~B
    an = A & ~B
    n_p, n_a, n_i = int(pn.sum()), int(an.sum()), int((pn & an).sum())
    ang = ag.difference(bg) if bg is not None else ag
    off = ang.difference(gb).area / ang.area if ang.area else 0.0
    return dict(act_new_ac=analyze.acres(n_a, px), pred_new_ac=analyze.acres(n_p, px),
                inter_ac=analyze.acres(n_i, px),
                precision=round(n_i / n_p, 4) if n_p else "",
                recall=round(n_i / n_a, 4) if n_a else "",
                iou=round(n_i / (n_p + n_a - n_i), 4) if (n_p + n_a - n_i) else "",
                offgrid=round(off, 4))


def score_fire(slug):
    ps = perimset(slug)
    hs = analyze.HotspotSet(slug)
    rows = []
    for run in runs_for(slug):
        if not (FIRES / slug / "runs" / run / "meta.json").exists():
            continue
        dt = run_dt(run)
        b = baseline(ps, dt)
        perims = [(t, e, basis) for t, e, basis in ps.timeline
                  if dt < t <= dt + timedelta(hours=MAX_LEAD_H) and e != b[1]]
        dash = [(H, ps.nearest(dt + timedelta(hours=H), analyze.ACTUAL_TOL_H)) for H in DASH_H]
        hpasses = passes(hs, dt, dt + timedelta(hours=MAX_LEAD_H))
        start = analyze.RunStart(f"{slug}/{run}", dt, hs)
        for model in MODELS:
            fp_cache = {}
            for pct in analyze.PERCENTILES:
                p = tif_for(model, slug, run, pct)
                if not p.exists():
                    continue
                with rasterio.open(p) as ds:
                    arr = ds.read(1)
                    crs, tfm, nod, gb = ds.crs, ds.transform, ds.nodata, box(*ds.bounds)
                px = abs(tfm.a * tfm.e)
                valid = (arr != nod) if nod is not None else np.isfinite(arr)
                Bp = ps.mask(b[1], crs, tfm, arr.shape)
                gkey = (str(crs), tuple(tfm)[:6], arr.shape)
                if gkey not in fp_cache:
                    S = start.mask(crs, tfm, arr.shape)
                    I = ignition_mask(slug, run, crs, tfm, arr.shape)
                    fp_cache[gkey] = (Bp | S | I, S & ~Bp, I & ~Bp & ~S)
                B, s_extra, i_extra = fp_cache[gkey]
                bg = ps.geom_utm(b[1], crs)
                base = dict(slug=slug, run_ts=run, model=model, pct=pct, baseline="common",
                            start_basis=start.basis,
                            start_extra_ac=analyze.acres(int(s_extra.sum()), px),
                            ign_extra_ac=analyze.acres(int(i_extra.sum()), px))
                for t, e, basis in perims:
                    lead = (t - dt).total_seconds() / 3600
                    A, ag = ps.mask(e, crs, tfm, arr.shape), ps.geom_utm(e, crs)
                    for bl, BB in (("common", B), ("perimeter", Bp)):
                        r = growth_row(arr, valid, BB, bg, A, ag, gb, lead, px)
                        rows.append(dict(base, baseline=bl, truth="perim", truth_id=e, truth_basis=basis,
                                         lead_h=round(lead, 3), **r,
                                         passes=int(r["act_new_ac"] >= F["min_act_new_acres"]
                                                    and r["offgrid"] <= F["max_offgrid_frac"])))
                for H, a in dash:
                    if not a:
                        continue
                    a_dt, a_ms, a_off, a_basis = a
                    r = growth_row(arr, valid, B, bg, ps.mask(a_ms, crs, tfm, arr.shape),
                                   ps.geom_utm(a_ms, crs), gb, H, px)
                    rows.append(dict(base, truth="dash", truth_id=a_ms, truth_basis=a_basis,
                                     lead_h=H, actual_offset_h=round(a_off, 2), **r,
                                     passes=int(r["act_new_ac"] >= F["min_act_new_acres"]
                                                and abs(a_off) <= F["max_abs_offset_h"]
                                                and r["offgrid"] <= F["max_offgrid_frac"])))
                for t, grp in hpasses:
                    key = (t, len(grp))
                    if key not in fp_cache:
                        geoms = [mapping(Point(lon, lat)) for _, lon, lat, rad in grp]
                        polys = [shape(transform_geom("EPSG:4326", crs, g)).buffer(rad)
                                 for g, (_, _, _, rad) in zip(geoms, grp)]
                        fp = features.rasterize([(mapping(unary_union(polys)), 1)], out_shape=arr.shape,
                                                transform=tfm, fill=0, dtype="uint8").astype(bool)
                        fp_cache[key] = fp
                    new = fp_cache[key] & ~B
                    n_new = int(new.sum())
                    if not n_new:
                        continue
                    lead = (t - dt).total_seconds() / 3600
                    hit = int((new & valid & (arr <= lead)).sum())
                    rows.append(dict(base, truth="hotspot", truth_id=t.strftime("%Y%m%dT%H%MZ"),
                                     truth_basis=f"{len(grp)} det", lead_h=round(lead, 3),
                                     hs_new_ac=analyze.acres(n_new, px), hs_hit_ac=analyze.acres(hit, px),
                                     capture=round(hit / n_new, 4),
                                     pred_new_ac=analyze.acres(int((valid & (arr <= lead) & ~B).sum()), px),
                                     passes=int(n_new * px / analyze.SQM_PER_ACRE >= 10)))
        ps._mask.clear()
    return rows


def paired(rows, truth, baseline="common"):
    by = defaultdict(dict)
    for r in rows:
        if r["truth"] == truth and r["baseline"] == baseline:
            by[(r["slug"], r["run_ts"], r["pct"], r["truth_id"], r["lead_h"])][r["model"]] = r
    return [(v["pyrecast"], v["farsite"]) for v in by.values()
            if len(v) == 2 and v["pyrecast"]["passes"] and v["farsite"]["passes"]]


def binned(pairs, width):
    out = defaultdict(list)
    for p, f in pairs:
        out[int(float(p["lead_h"]) // width) * width].append((p, f))
    return dict(sorted(out.items()))


def summarize(rows, title):
    iou = lambda r: float(r["iou"] or 0)
    print(f"\n=== {title}")
    pp = paired(rows, "perim")
    print(f"perimeters at exact lead time: {len(pp)} paired rows (all percentiles; mapped-time "
          f"{sum(p['truth_basis'] == 'mapped' for p, _ in pp)})")
    print(f"{'lead':>9} {'n':>4} {'IoU P':>6} {'IoU F':>6} {'F wins':>6} {'over P':>7} {'over F':>7}")
    for lo, g in binned(pp, 12).items():
        ov = lambda s: sum(float(x[s]["pred_new_ac"]) for x in g) / max(sum(float(x[s]["act_new_ac"]) for x in g), 1e-9)
        print(f"{lo:3d}-{lo + 12:3d}h {len(g):4d} {np.mean([iou(p) for p, _ in g]):6.3f} "
              f"{np.mean([iou(f) for _, f in g]):6.3f} {np.mean([iou(f) > iou(p) for p, f in g]):6.0%} "
              f"{ov(0):7.1f} {ov(1):7.1f}")
    hp = paired(rows, "hotspot")
    print(f"satellite passes: {len(hp)} paired rows")
    print(f"{'lead':>9} {'n':>4} {'capt P':>6} {'capt F':>6} {'predac P':>9} {'predac F':>9}")
    for lo, g in binned(hp, 12).items():
        print(f"{lo:3d}-{lo + 12:3d}h {len(g):4d} {np.mean([float(p['capture']) for p, _ in g]):6.2f} "
              f"{np.mean([float(f['capture']) for _, f in g]):6.2f} "
              f"{np.median([float(p['pred_new_ac']) for p, _ in g]):9.0f} {np.median([float(f['pred_new_ac']) for _, f in g]):9.0f}")
    bp = paired(rows, "perim", "perimeter")
    if bp:
        print(f"same perimeters scored against the official perimeter only ({len(bp)} rows): mean IoU "
              f"P {np.mean([iou(p) for p, _ in bp]):.3f}  F {np.mean([iou(f) for _, f in bp]):.3f}")
    dp = paired(rows, "dash")
    print(f"dashboard horizons (default filters): {len(dp)} paired rows")
    for H in DASH_H:
        g = [(p, f) for p, f in dp if int(float(p["lead_h"])) == H]
        if g:
            print(f"  {H:3d}h n={len(g):3d} IoU P {np.mean([iou(p) for p, _ in g]):.3f}  "
                  f"IoU F {np.mean([iou(f) for _, f in g]):.3f}  F wins {np.mean([iou(f) > iou(p) for p, f in g]):.0%}")


def bins_table(rows):
    out = []
    for slug in sorted({r["slug"] for r in rows}) + ["all"]:
        sub = [r for r in rows if slug == "all" or r["slug"] == slug]
        for truth in ("perim", "hotspot"):
            for pct in ("all",) + analyze.PERCENTILES:
                pr = [(p, f) for p, f in paired(sub, truth) if pct == "all" or p["pct"] == pct]
                for lo, g in binned(pr, 1).items():
                    rec = dict(fire=slug, truth=truth, pct=pct, lead_bin_h=lo, n=len(g))
                    for side, m in ((0, "pyrecast"), (1, "farsite")):
                        k = "iou" if truth == "perim" else "capture"
                        rec[f"{k}_{m}"] = round(float(np.mean([float(x[side][k] or 0) for x in g])), 4)
                        rec[f"pred_new_ac_{m}"] = round(float(np.mean([float(x[side]["pred_new_ac"]) for x in g])), 1)
                    out.append(rec)
    return out


def main():
    slugs = sys.argv[1:] or ["ca-dome", "ca-plaskett"]
    rows = []
    for s in slugs:
        rows += score_fire(s)
    RUNS.mkdir(exist_ok=True)
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(RUNS / "hindcast_scores.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    tb = bins_table(rows)
    cols = []
    for r in tb:
        cols += [k for k in r if k not in cols]
    with open(RUNS / "hindcast_bins.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(tb)
    for s in slugs:
        summarize([r for r in rows if r["slug"] == s], s)
    summarize(rows, "both fires")


if __name__ == "__main__":
    main()
