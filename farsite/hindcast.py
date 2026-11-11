"""FARSITE hindcast for archived PyreCast runs (pilot on Plaskett + Dome).

For every complete PyreCast run of a fire: start from the same perimeter at the same time,
weather = NBM as published by then, analyst-style rules (active-edge ignition on a
burn-scar landscape, burn periods centred on 16:00, no spotting), a 25-member ensemble
(wind speed x, wind direction +, dead fuel moisture +, burn period length), published as
the members at the 10/30/50/70/90th percentile of burned area (PyreCast's five levels)
plus a burn-probability grid.

Variants (--variant):
  base  the v1 pilot (runs/): gust speed as the wind during burn periods, no calibration.
  cal   sustained wind all day, plus calibration to the previous day's growth. Before the
        forecast, CAL_MEMBERS members replay the previous day between two official perimeters:
        from the latest perimeter P0 20-54 h before the run's baseline perimeter B, with the
        NBM published by P0's time, up to B's time. Their median new area is compared with the
        observed growth, B beyond P0. (Detection footprints are not used: 375 m VIIRS pixels
        along a slow edge read as ~5x the perimeter growth on Dome.) The ratio is a factor on
        burn-period length, i.e. hours of active spread per day: the lever analysts tune to
        observed activity, and unlike FARSITE's ROS adjustment file it slows crown fire too.
        A second replay at the first factor refines it (secant in log space). Two flavours
        share the replay:
          runs-cal/   the factor holds for the whole forecast
          runs-calr/  the factor's effect (in log space) halves every CAL_HALF_LIFE burn periods,
                      easing back toward uncalibrated at longer leads
  calbar  cal, plus completed control lines as FARSITE barriers (runs-calbar/, runs-calbarr/):
        dozer, hand, road-as-line and mixed construction lines from the latest vector IAP ops
        map published (FTP time) by the run's start, within BARRIER_NEAR_M of the fire
        (interventions/<slug>.geojson, from interventions.py). The replay uses the lines
        published by its own start. Fires with no vector maps run as cal.
  Every run records cal_* fields in meta.json (factor, observed and replayed acres, window),
  and calbar runs barrier_* fields (map, km of line).

usage: hindcast.py SLUG [--variant base|cal|calbar] [--members N] [--workers N] [--concurrent-runs N]
                         [--limit N] [--only RUN_TS] [--shard I/N] [--prep-only]
--shard I/N runs the I-th of N contiguous slices of the run list (GitHub Actions matrix);
--prep-only builds the landscape and stops.
outputs: data/fires/SLUG/runs[-FLAVOUR]/RUN_TS/{toa_p10..p90.tif, burnprob.tif, members.csv, meta.json}
"""
import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import rasterio
from rasterio import features
from rasterio.shutil import copy as rio_copy
from rasterio.warp import transform as warp_transform, transform_geom
from shapely.geometry import MultiPolygon, Point, Polygon, mapping, shape
from shapely.ops import unary_union

from common import AC, FARSITE_BIN, FIRES, analyze, baseline, geom_5070, perimset, run_dt, runs_for
import wx

FARSITE = FARSITE_BIN
TZ = ZoneInfo("America/Los_Angeles")           # pilot fires are all Pacific time
HORIZON_H = 96
MEMBER_TIMEOUT_S = int(os.environ.get("FARSITE_MEMBER_TIMEOUT_S", 600))  # p99 of 500 Dome members: 307 s
# FARSITE's cost grows with the active perimeter. Big active edges (Plaskett: 2,500-7,500 ac lit)
# start at half resolution, as an analyst would for a big fire; a timeout halves it again.
LARGE_IGNITION_AC, RES_M, RES_LARGE_M = 3000, (30, 60), (60, 120)   # (distance, perimeter) m
BURN_CENTER = 16.0                              # burn periods centred on 16:00 local
BASE_MOIST = (5, 6, 8, 30, 90)                  # 1h 10h 100h liveH liveW
EDGE_BAND_M, HS_PAD_M, MIN_PART_AC = 150.0, 300.0, 1.0
MARGIN_KM = 15.0
PCTS = ("10", "30", "50", "70", "90")
LCP_OPTS = dict(ELEVATION_UNIT="METERS", SLOPE_UNIT="DEGREES", ASPECT_UNIT="AZIMUTH_DEGREES",
                FUEL_MODEL_OPTION="NO_CUSTOM_AND_NO_FILE", CANOPY_COV_UNIT="PERCENT",
                CANOPY_HT_UNIT="METERS_X_10", CBH_UNIT="METERS_X_10",
                CBD_UNIT="KG_PER_CUBIC_METER_X_100", LINEAR_UNIT="METER")

# flavours: output name -> calibration half-life in burn periods (None = constant factor)
VARIANTS = {
    "base": dict(gust=True, calib=False, barriers=False, flavors={"base": None}),
    "cal": dict(gust=False, calib=True, barriers=False, flavors={"cal": None, "calr": 4}),
    "calbar": dict(gust=False, calib=True, barriers=True, flavors={"calbar": None, "calbarr": 4}),
}
INTERVENTIONS = Path(__file__).resolve().parent / "interventions"
BARRIER_NEAR_M = 5000.0       # control lines farther than this from the fire are left out
CAL_WINDOW_H = (20, 54)       # replay from the latest perimeter 20-54 h before the baseline perimeter
CAL_MEMBERS = 9               # replay members (m00-m08); the median sets the factor
CAL_RANGE = (0.05, 2.5)       # factor clamp
DET_EDGE_M = 200.0            # (diagnostic only) detections this close to P0 are not growth
BP_RANGE_H = (0.25, 16.0)     # calibrated burn period length clamp (16 h = 08:00-24:00)


def run_dir(slug, flavor, run_ts):
    return FIRES / slug / ("runs" if flavor == "base" else f"runs-{flavor}") / run_ts


# ---------------------------------------------------------------- landscape (per fire)
def landscape(slug, run_list):
    import landscape as L
    d = FIRES / slug
    if (d / "lcp.tif").exists():
        return json.load(open(d / "grid.json"))
    d.mkdir(parents=True, exist_ok=True)
    ps = perimset(slug)
    boxes = np.array([geom_5070(slug, baseline(ps, run_dt(r))[1]).bounds for r in run_list])
    m = MARGIN_KM * 1000
    xmin = L.X0 + np.floor((boxes[:, 0].min() - m - L.X0) / L.RES) * L.RES
    ymax = L.Y0 - np.floor((L.Y0 - (boxes[:, 3].max() + m)) / L.RES) * L.RES
    w = int(np.ceil((boxes[:, 2].max() + m - xmin) / L.RES))
    h = int(np.ceil((ymax - (boxes[:, 1].min() - m)) / L.RES))
    bbox = (float(xmin), float(ymax - h * L.RES), float(xmin + w * L.RES), float(ymax))
    stack = np.stack([L.fetch(path, bbox, w, h) for _, path in L.LAYERS]).astype(np.int16)
    nod = stack == -9999
    stack[3][nod[3]] = 98 if nod[0].any() else 99          # ocean/no-data -> water
    for i in (0, 1, 4, 5, 6, 7):
        stack[i][nod[i]] = 0
    stack[2][nod[2]] = -1
    tfm = rasterio.transform.from_origin(bbox[0], bbox[3], L.RES, L.RES)
    prof = dict(driver="GTiff", width=w, height=h, count=8, dtype="int16", crs="EPSG:5070",
                transform=tfm, compress="deflate")
    with rasterio.open(d / "lcp.tif", "w", **prof) as ds:
        ds.write(stack)
    lat = warp_transform("EPSG:5070", "EPSG:4326", [(bbox[0] + bbox[2]) / 2], [(bbox[1] + bbox[3]) / 2])[1][0]
    g = {"bbox": bbox, "w": w, "h": h, "lat": lat}
    json.dump(g, open(d / "grid.json", "w"))
    print(f"{slug}: landscape {w}x{h} ({w * 0.03:.0f}x{h * 0.03:.0f} km), "
          f"nodata cells {int(nod[3].sum())}", flush=True)
    return g


# ---------------------------------------------------------------- ignition (per run)
_HS = {}


def hotspots(slug):
    if slug not in _HS:
        _HS[slug] = analyze.HotspotSet(slug)
    return _HS[slug]


def _polys(g):
    if isinstance(g, Polygon):
        return [g]
    return [p for p in getattr(g, "geoms", []) if isinstance(p, Polygon)]


def _footprints(pts, pad=0.0):
    return unary_union([shape(transform_geom("EPSG:4326", "EPSG:5070", mapping(Point(lon, lat)))).buffer(r + pad)
                        for _, lon, lat, r, *_ in pts]) if pts else Polygon()   # (dt, lon, lat, radius_m[, source])


def active_edge(slug, base, base_ms, dt):
    """Active edge of perimeter `base` at time dt: growth since the perimeter 18 h earlier,
    plus the edge band near detections from the last 24 h. -> (parts, growth, detections)"""
    ps = perimset(slug)
    prev = ps.latest_at_or_before(dt - timedelta(hours=18))
    growth = Polygon()
    if prev and prev[1] != base_ms:
        growth = base.difference(geom_5070(slug, prev[1])).buffer(-15).buffer(15)
    pts = [p for p in hotspots(slug).pts if dt - timedelta(hours=24) < p[0] <= dt]
    foot = _footprints(pts, HS_PAD_M)
    band = base.difference(base.buffer(-EDGE_BAND_M))
    active = unary_union([growth, band.intersection(foot)]).intersection(base)
    parts = [Polygon(p.exterior) for p in _polys(active) if p.area / AC >= MIN_PART_AC]
    return parts, growth, pts


def write_inputs(slug, d, base, parts):
    """ign.shp (+ ign.geojson) and scar.lcp in d: the active edge is lit, the rest of the
    burn scar is non-burnable"""
    d.mkdir(parents=True, exist_ok=True)
    gj = d / "ign.geojson"
    json.dump({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"id": i}, "geometry": mapping(p)} for i, p in enumerate(parts)]},
        open(gj, "w"))
    subprocess.run(["ogr2ogr", "-f", "ESRI Shapefile", str(d / "ign.shp"), str(gj), "-a_srs",
                    "EPSG:5070", "-nlt", "POLYGON"], check=True, capture_output=True)
    with rasterio.open(FIRES / slug / "lcp.tif") as ds:
        stack, prof = ds.read(), ds.profile
        scar = base.difference(unary_union(parts).buffer(30))
        m = features.rasterize([(mapping(scar), 1)], out_shape=(ds.height, ds.width),
                               transform=ds.transform, fill=0, dtype="uint8").astype(bool) \
            if not scar.is_empty else np.zeros((ds.height, ds.width), bool)
    stack[3] = np.where(m, 99, stack[3])
    for i in (4, 5, 6, 7):
        stack[i] = np.where(m, 0, stack[i])
    with rasterio.open(d / "scar.tif", "w", **prof) as ds:
        ds.write(stack)
    lat = json.load(open(FIRES / slug / "grid.json"))["lat"]
    rio_copy(d / "scar.tif", d / "scar.lcp", driver="LCP", LATITUDE=str(int(round(lat))), **LCP_OPTS)
    (d / "scar.tif").unlink()


def ignition(slug, run_ts, rdir):
    """Active-edge ignition polygons + burn-scar LCP for one run. Returns info dict."""
    dt = run_dt(run_ts)
    b = baseline(perimset(slug), dt)
    base = geom_5070(slug, b[1])
    parts, growth, pts = active_edge(slug, base, b[1], dt)
    info = dict(baseline_ms=b[1], baseline_time=f"{b[0]:%Y-%m-%dT%H:%MZ}", baseline_ac=round(base.area / AC),
                growth_ac=round(growth.area / AC), detections_24h=len(pts),
                ignition_ac=round(sum(p.area for p in parts) / AC), ignition_parts=len(parts),
                weather_point=list(warp_transform("EPSG:5070", "EPSG:4326",
                                                  [base.centroid.x], [base.centroid.y])))
    info["weather_point"] = [info["weather_point"][0][0], info["weather_point"][1][0]]
    if parts:
        write_inputs(slug, rdir, base, parts)
    return info


def prior_perimeter(slug, run_ts):
    """The calibration replay's start: the latest perimeter 20-54 h before the run's baseline
    perimeter, or None"""
    ps = perimset(slug)
    tb = baseline(ps, run_dt(run_ts))[0]
    p0 = ps.latest_at_or_before(tb - timedelta(hours=CAL_WINDOW_H[0]))
    if p0 and tb - p0[0] <= timedelta(hours=CAL_WINDOW_H[1]):
        return p0
    return None


_IV = {}


def barriers_at(slug, when, near):
    """Completed control lines from the latest vector ops map published by `when`, within
    BARRIER_NEAR_M of geometry `near`. -> (lines, map filename or None)"""
    p = INTERVENTIONS / f"{slug}.geojson"
    if not p.exists():
        return [], None
    if slug not in _IV:
        _IV[slug] = json.load(open(p))["features"]
    pubs = sorted({f["properties"]["published"] for f in _IV[slug]})
    known = [x for x in pubs if datetime.fromisoformat(x.replace("Z", "+00:00")) <= when]
    if not known:
        return [], None
    zone = near.buffer(BARRIER_NEAR_M)
    fs = [f for f in _IV[slug] if f["properties"]["published"] == known[-1]]
    lines = [(f["properties"]["kind"], shape(f["geometry"])) for f in fs]
    return [(k, g) for k, g in lines if g.intersects(zone)], ",".join(sorted({f["properties"]["map"] for f in fs}))


def write_barriers(d, lines):
    """barrier.shp (+ barrier.geojson) in d for FARSITE's barrier slot; none if no lines"""
    for f in d.glob("barrier.*"):
        f.unlink()
    if not lines:
        return
    gj = d / "barrier.geojson"
    json.dump({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"kind": k}, "geometry": mapping(g)} for k, g in lines]}, open(gj, "w"))
    subprocess.run(["ogr2ogr", "-f", "ESRI Shapefile", str(d / "barrier.shp"), str(gj), "-a_srs",
                    "EPSG:5070", "-nlt", "LINESTRING", "-explodecollections"], check=True, capture_output=True)


# ---------------------------------------------------------------- members
def members(slug, run_ts, n):
    rng = random.Random(zlib.crc32(f"{slug}/{run_ts}".encode()))
    out = [dict(m=0, ws=1.0, wd=0.0, dfm=0.0, bp=8)]
    for i in range(1, n):
        out.append(dict(m=i, ws=round(min(max(np.exp(rng.gauss(0, 0.25)), 0.6), 1.6), 3),
                        wd=round(rng.gauss(0, 15), 1), dfm=round(rng.uniform(-2, 2), 2),
                        bp=rng.choice((6, 7, 8, 9, 10))))
    return out


def input_text(start, end, rows, elev_ft, mem, res=RES_M, gust=True, factors=None):
    """FARSITE inputs file. factors(k) scales the k-th burn period's length (calibration);
    gust=True uses gust speed as the wind inside the member's nominal burn period (v1)."""
    s, e = start.astimezone(TZ), end.astimezone(TZ)
    hhmm = lambda h: int(h) * 100 + int(round((h % 1) * 60))
    periods, d = [], s.date()
    while d <= e.date():
        L = mem["bp"] if factors is None else min(max(mem["bp"] * factors(len(periods)), BP_RANGE_H[0]), BP_RANGE_H[1])
        bp0, bp1 = BURN_CENTER - L / 2, BURN_CENTER + L / 2
        p0 = datetime.combine(d, datetime.min.time(), TZ) + timedelta(hours=bp0)
        p1 = datetime.combine(d, datetime.min.time(), TZ) + timedelta(hours=bp1) - timedelta(minutes=1)
        if p1 > s:
            p0 = max(p0, s)
            if p0 < min(p1, e):
                periods.append(f"{d:%m %d} {p0.hour * 100 + p0.minute:04d} {p1.hour * 100 + p1.minute:04d}")
        d += timedelta(days=1)
    g0, g1 = hhmm(BURN_CENTER - mem["bp"] / 2), hhmm(BURN_CENTER + mem["bp"] / 2)
    lines = []
    for r in rows:
        lt = r["t"].astimezone(TZ)
        in_bp = gust and g0 <= lt.hour * 100 < g1
        ms = (r["G"] if in_bp else r["WS"]) * mem["ws"]
        lines.append(f"{lt.year} {lt.month:02d} {lt.day:02d} {lt.hour * 100:04d} {r['T'] * 9 / 5 + 32:.1f} "
                     f"{min(max(r['RH'], 1), 99):.0f} {r['P'] * 0.03937:.2f} {ms * 2.23694 / 1.15:.1f} "
                     f"{(r['WD'] + mem['wd']) % 360:.0f} 0")
    m1, m10, m100, lh, lw = BASE_MOIST
    f = lambda v: max(2, round(v + mem["dfm"]))
    txt = f"""FARSITE INPUTS FILE VERSION 1.0
FARSITE_START_TIME: {s:%m %d %H%M}
FARSITE_END_TIME: {e:%m %d %H%M}
FARSITE_TIMESTEP: 30
FARSITE_DISTANCE_RES: {res[0]:.1f}
FARSITE_PERIMETER_RES: {res[1]:.1f}
FARSITE_MIN_IGNITION_VERTEX_DISTANCE: 15.0
FARSITE_SPOT_GRID_RESOLUTION: 30.0
FARSITE_SPOT_PROBABILITY: 0.0
FARSITE_SPOT_IGNITION_DELAY: 0
FARSITE_MINIMUM_SPOT_DISTANCE: 30
FARSITE_ACCELERATION_ON: 1
FARSITE_FILL_BARRIERS: 1
SPOTTING_SEED: 1000
FUEL_MOISTURES_DATA: 1
0 {f(m1)} {f(m10)} {f(m100)} {lh} {lw}
RAWS_ELEVATION: {elev_ft:.0f}
RAWS_UNITS: English
RAWS: {len(lines)}
""" + "\n".join(lines) + """
FOLIAR_MOISTURE_CONTENT: 100
CROWN_FIRE_METHOD: Finney
NUMBER_PROCESSORS: 1
"""
    if periods:
        txt += f"FARSITE_BURN_PERIODS: {len(periods)}\n" + "\n".join(periods) + "\n"
    return txt, periods


def run_member(md, wdir, mem, txt, res=RES_M):
    """One FARSITE member in md, lit from wdir/ign.shp on wdir/scar.lcp. -> (mem, arrival min, secs, err)"""
    md.mkdir(parents=True, exist_ok=True)
    (md / "in.input").write_text(txt)
    bar = wdir / "barrier.shp" if (wdir / "barrier.shp").exists() else 0
    (md / "cmd.txt").write_text(f"{wdir / 'scar.lcp'} {md / 'in.input'} {wdir / 'ign.shp'} {bar} {md / 'o'} 0\n")
    asc = md / "o_ArrivalTime.asc"
    t0, err = time.time(), ""
    # FARSITE occasionally blows up on perimeter vertices (40+ CPU-min vs a 22 s median).
    # Kill at MEMBER_TIMEOUT_S and retry once at half the perimeter/distance resolution.
    coarse = (txt.replace(f"FARSITE_DISTANCE_RES: {res[0]:.1f}", f"FARSITE_DISTANCE_RES: {2 * res[0]:.1f}")
                 .replace(f"FARSITE_PERIMETER_RES: {res[1]:.1f}", f"FARSITE_PERIMETER_RES: {2 * res[1]:.1f}"))
    for attempt, t in enumerate((txt, coarse)):
        (md / "in.input").write_text(t)
        for old in md.glob("o_*"):
            old.unlink()
        try:
            p = subprocess.run([str(FARSITE), str(md / "cmd.txt")], capture_output=True, text=True,
                               timeout=MEMBER_TIMEOUT_S)
            err = (p.stdout + p.stderr)[-400:]
            break
        except subprocess.TimeoutExpired:
            err = f"timeout {MEMBER_TIMEOUT_S}s attempt {attempt}"
    mem = dict(mem, coarse=attempt)
    secs = time.time() - t0
    if not asc.exists() or err.startswith("timeout"):
        return mem, None, secs, err
    with rasterio.open(asc) as ds:
        a = ds.read(1)
    return mem, a, secs, ""


# ---------------------------------------------------------------- calibration (per run)
def calibrate(slug, run_ts, wdir, mems, pool, res, gust, wxrows, elev_ft, barriers=False):
    """Replay the previous day (perimeter P0 -> the baseline perimeter B) and return the
    burn-period factor that matches its observed growth, with the evidence (cal_* fields)."""
    out = dict(cal_factor=1.0)
    dt = run_dt(run_ts)
    p0 = prior_perimeter(slug, run_ts)
    if not p0:
        return dict(out, cal_status="no_prior_perimeter")
    if not wxrows:
        return dict(out, cal_status="no_weather")
    t0, P0 = p0[0], geom_5070(slug, p0[1])
    b = baseline(perimset(slug), dt)
    t1, B = b[0], geom_5070(slug, b[1])
    win = (t1 - t0).total_seconds() / 3600
    out.update(cal_p0_ms=p0[1], cal_p0_time=f"{t0:%Y-%m-%dT%H:%MZ}", cal_p1_time=f"{t1:%Y-%m-%dT%H:%MZ}",
               cal_window_h=round(win, 1))
    # observed: the baseline perimeter beyond P0 (detection footprints recorded for reference)
    obs_ac = B.difference(P0).area / AC
    det = [p for p in hotspots(slug).pts if t0 < p[0] <= t1]
    det_ac = _footprints(det).difference(P0.buffer(DET_EDGE_M)).difference(B).area / AC if det else 0.0
    out.update(cal_obs_ac=round(obs_ac), cal_det_outside_ac=round(det_ac), cal_obs_det=len(det))
    parts, _, _ = active_edge(slug, P0, p0[1], t0)
    if not parts:
        return dict(out, cal_status="no_replay_edge")
    cdir = wdir / "cal"
    write_inputs(slug, cdir, P0, parts)
    if barriers:
        lines, src = barriers_at(slug, t0, P0)
        write_barriers(cdir, lines)
        out.update(cal_barrier_map=src, cal_barrier_km=round(sum(g.length for _, g in lines) / 1000, 1))
    out["cal_ignition_ac"] = round(sum(p.area for p in parts) / AC)
    with rasterio.open(FIRES / slug / "lcp.tif") as ds:
        p0m = features.rasterize([(mapping(P0), 1)], out_shape=(ds.height, ds.width),
                                 transform=ds.transform, fill=0, dtype="uint8").astype(bool)

    def replay(f, tag):
        futs = []
        for mem in mems[:CAL_MEMBERS]:
            txt, _ = input_text(t0, t1, wxrows, elev_ft, mem, res, gust=gust, factors=lambda k: f)
            futs.append(pool.submit(run_member, cdir / f"{tag}{mem['m']:02d}", cdir, mem, txt, res))
        g = []
        for fu in futs:
            _, a, _, _ = fu.result()
            if a is not None:
                h = np.where(a != -9999, a / 60.0, np.nan)
                g.append(float((np.isfinite(h) & (h > 0) & (h <= win) & ~p0m).sum() * 900 / AC))
        return float(np.median(g)) if g else None

    lo, hi = CAL_RANGE
    clamp = lambda v: min(max(v, lo), hi)
    if obs_ac < 1:                  # the perimeter did not grow: no replay needed, the floor it is
        return dict(out, cal_status="ok_no_growth", cal_factor=lo)
    g1 = replay(1.0, "a")
    out["cal_pred_ac"] = None if g1 is None else round(g1)
    if not g1:
        return dict(out, cal_status="replay_failed" if g1 is None else "replay_no_growth")
    r1 = obs_ac / g1
    f = clamp(r1)
    if lo < f < hi and abs(np.log(r1)) > np.log(1.25):
        g2 = replay(f, "b")
        if g2:
            s = min(max(np.log(g2 / g1) / np.log(f), 0.3), 3.0)
            out.update(cal_pred2_ac=round(g2), cal_pred2_factor=round(f, 4), cal_slope=round(float(s), 2))
            f = clamp(f * (obs_ac / g2) ** (1 / s))
    return dict(out, cal_status="ok", cal_factor=round(float(f), 4))


# ---------------------------------------------------------------- one run
def do_run(slug, run_ts, wxrows, calrows, n_members, pool, lcp_profile, variant="base"):
    cfg = VARIANTS[variant]
    dirs = {fl: run_dir(slug, fl, run_ts) for fl in cfg["flavors"]}
    wdir = next(iter(dirs.values()))       # holds the ignition + scar landscape for every flavour
    if all((d / "meta.json").exists() for d in dirs.values()):
        return json.load(open(wdir / "meta.json"))
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    dt = run_dt(run_ts)
    t0 = time.time()
    info = ignition(slug, run_ts, wdir)
    meta = dict(slug=slug, run_ts=run_ts, variant=variant, **info, members=n_members)
    start, end = dt, dt + timedelta(hours=HORIZON_H + 1)
    shape_ = (lcp_profile["height"], lcp_profile["width"])
    if not info["ignition_parts"] or not wxrows:
        meta["status"] = "no_active_edge" if not info["ignition_parts"] else "no_weather"
        empty = np.zeros(shape_, "float32")
        for fl, d in dirs.items():
            for p in PCTS:
                _write(d / f"toa_p{p}.tif", empty, lcp_profile)
            json.dump(dict(meta, flavor=fl), open(d / "meta.json", "w"), indent=1)
        return meta
    with rasterio.open(FIRES / slug / "lcp.tif") as ds:
        x, y = warp_transform("EPSG:4326", "EPSG:5070", [info["weather_point"][0]], [info["weather_point"][1]])
        r, c = ds.index(x[0], y[0])
        el = ds.read(1)[max(r - 40, 0):r + 40, max(c - 40, 0):c + 40]
    elev_ft = float(el[el > 0].mean() if (el > 0).any() else 0) * 3.28084
    mems = members(slug, run_ts, n_members)
    res = RES_LARGE_M if info["ignition_ac"] > LARGE_IGNITION_AC else RES_M
    meta["res_m"] = list(res)
    if cfg["barriers"]:
        lines, src = barriers_at(slug, dt, geom_5070(slug, info["baseline_ms"]))
        write_barriers(wdir, lines)
        meta.update(barrier_map=src, barrier_km=round(sum(g.length for _, g in lines) / 1000, 1),
                    barrier_km_by_kind={k: round(sum(g.length for kk, g in lines if kk == k) / 1000, 1)
                                        for k in sorted({k for k, _ in lines})})
    if cfg["calib"]:
        meta.update(calibrate(slug, run_ts, wdir, mems, pool, res, cfg["gust"], calrows, elev_ft,
                              barriers=cfg["barriers"]))
    fac = meta.get("cal_factor", 1.0)
    jobs = {}
    for fl, half in cfg["flavors"].items():
        if not cfg["calib"]:
            factors = None
        elif half is None:
            factors = lambda k, f=fac: f
        else:
            factors = lambda k, f=fac, h=half: f ** (0.5 ** (k / h))
        for mem in mems:
            txt, periods = input_text(start, end, wxrows, elev_ft, mem, res, gust=cfg["gust"], factors=factors)
            jobs[pool.submit(run_member, dirs[fl] / f"m{mem['m']:02d}", wdir, mem, txt, res)] = (fl, periods)
    toas = {fl: {} for fl in dirs}
    rows = {fl: [] for fl in dirs}
    for fu in as_completed(jobs):
        fl, periods = jobs[fu]
        mem, a, secs, err = fu.result()
        if a is None:
            rows[fl].append(dict(mem, secs=round(secs, 1), new_ac="", stalled="", error=err.replace("\n", " ")))
            continue
        h = np.where(a != -9999, a / 60.0, np.nan).astype("float32")   # hours since run (start = run)
        toas[fl][mem["m"]] = h
        new = np.isfinite(h) & (h > 0) & (h <= HORIZON_H)
        rows[fl].append(dict(mem, secs=round(secs, 1), new_ac=round(new.sum() * 900 / AC),
                             stalled=int(_stalled(h, dt, periods)), error=""))
    wall = round(time.time() - t0, 1)
    out = {}
    for fl, d in dirs.items():
        m = dict(meta, flavor=fl)
        if cfg["calib"]:
            m["cal_half_life_bp"] = cfg["flavors"][fl]
        _finish(d, rows[fl], toas[fl], m, lcp_profile)
        m["wall_secs"] = wall
        if d != wdir:
            for f in list(wdir.glob("ign.*")) + list(wdir.glob("barrier.geojson")):
                shutil.copy(f, d / f.name)
        if (wdir / "cal" / "ign.geojson").exists():
            shutil.copy(wdir / "cal" / "ign.geojson", d / "cal_ign.geojson")
        for md in d.glob("m[0-9][0-9]"):
            shutil.rmtree(md)
        json.dump(m, open(d / "meta.json", "w"), indent=1)
        out[fl] = m
    shutil.rmtree(wdir / "cal", ignore_errors=True)
    for ext in ("lcp", "lcp.aux.xml"):
        (wdir / f"scar.{ext}").unlink(missing_ok=True)
    for f in wdir.glob("barrier.*"):
        if f.suffix != ".geojson":
            f.unlink()
    first = out[next(iter(dirs))]
    return dict(first, flavors={fl: m.get("central_new_ac") for fl, m in out.items()})


def _finish(rdir, rows, toas, meta, lcp_profile):
    """members.csv, the five percentile members' arrival grids, burnprob; fills meta"""
    rows.sort(key=lambda r: r["m"])
    with open(rdir / "members.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    fails = sum(r["new_ac"] == "" for r in rows)
    ok = [r for r in rows if r["new_ac"] != ""]
    if not ok:
        meta.update(status="all_failed", failed=fails)
        return
    order = sorted(ok, key=lambda r: (r["new_ac"], r["m"]))
    pick = {}
    for p in PCTS:
        i = min(int(round(int(p) / 100 * (len(order) - 1))), len(order) - 1)
        pick[p] = order[i]["m"]
        _write(rdir / f"toa_p{p}.tif", np.nan_to_num(np.where(toas[pick[p]] <= 0, 1e-3, toas[pick[p]]), nan=0.0), lcp_profile)
    bp = np.mean([np.isfinite(t) & (t <= HORIZON_H) for t in toas.values()], axis=0).astype("float32")
    _write(rdir / "burnprob.tif", bp, lcp_profile)
    meta.update(status="ok", pick=pick, failed=fails, stalled=sum(r["stalled"] for r in ok),
                coarse_retries=sum(r["coarse"] for r in rows), member_timeout_s=MEMBER_TIMEOUT_S,
                member_secs=round(sum(r["secs"] for r in rows), 1),
                new_ac_min=order[0]["new_ac"], new_ac_max=order[-1]["new_ac"],
                central_new_ac=next((r["new_ac"] for r in ok if r["m"] == 0), ""))


def _stalled(h, dt, periods):
    """The day-1-then-stop failure: still spreading at the end of the first burn period,
    later burn periods exist, but nothing burns after the first one."""
    if len(periods) < 2 or not np.isfinite(h).any():
        return False
    mo, dy, a, b = periods[0].split()
    y = dt.astimezone(TZ).year
    p1_end = datetime(y, int(mo), int(dy), int(b) // 100, int(b) % 100, tzinfo=TZ)
    e = (p1_end - dt).total_seconds() / 3600
    late = np.isfinite(h) & (h > e - 1) & (h <= e + 0.1)
    after = np.isfinite(h) & (h > e + 1)
    return late.sum() > 50 and after.sum() == 0


def _write(path, arr, prof):
    p = dict(prof, count=1, dtype="float32", nodata=0, compress="deflate")
    p.pop("blockxsize", None); p.pop("blockysize", None); p.pop("tiled", None)
    with rasterio.open(path, "w", **p) as ds:
        ds.write(arr.astype("float32"), 1)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slug")
    ap.add_argument("--variant", default="base", choices=sorted(VARIANTS))
    ap.add_argument("--members", type=int, default=25)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--concurrent-runs", type=int, default=8)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--only")
    ap.add_argument("--shard", help="I/N: the I-th (0-based) of N contiguous slices")
    ap.add_argument("--prep-only", action="store_true")
    a = ap.parse_args()
    calib = VARIANTS[a.variant]["calib"]
    run_list = runs_for(a.slug)
    if a.only:
        run_list = [r for r in run_list if r in a.only.split(",")]
    if a.limit:
        run_list = run_list[:a.limit]
    if a.shard:        # contiguous, so neighbouring runs share NBM cycle downloads
        i, n = map(int, a.shard.split("/"))
        run_list = run_list[i * len(run_list) // n:(i + 1) * len(run_list) // n]
    landscape(a.slug, runs_for(a.slug))
    if a.prep_only or not run_list:
        print(f"{a.slug}: landscape ready, {len(run_list)} runs in this slice", flush=True)
        return
    d = FIRES / a.slug
    lat = json.load(open(d / "grid.json"))["lat"]
    if not (d / "lcp.lcp").exists():
        rio_copy(d / "lcp.tif", d / "lcp.lcp", driver="LCP", LATITUDE=str(int(round(lat))), **LCP_OPTS)
    with rasterio.open(d / "lcp.tif") as ds:
        prof = dict(ds.profile)
    # weather for every run at its baseline centroid (and, for calibration, for every replay
    # from its earlier perimeter's time at that perimeter's centroid), cached per run
    wxd = d / "wx"
    wxd.mkdir(exist_ok=True)
    todo = []
    ps = perimset(a.slug)

    def job(key, ms, when):
        if not (wxd / f"{key}.json").exists():
            g = geom_5070(a.slug, ms).centroid
            x, y = warp_transform("EPSG:5070", "EPSG:4326", [g.x], [g.y])
            todo.append((key, when, x[0], y[0]))

    for r in run_list:
        job(r, baseline(ps, run_dt(r))[1], run_dt(r))
        p0 = prior_perimeter(a.slug, r) if calib else None
        if p0:
            job(f"cal_{r}", p0[1], p0[0])
    if todo:
        t0 = time.time()
        for r, rows in wx.build_nbm(todo).items():
            json.dump([dict(x, t=x["t"].isoformat()) for x in rows], open(wxd / f"{r}.json", "w"))
        print(f"{a.slug}: weather for {len(todo)} runs/replays in {time.time() - t0:.0f}s", flush=True)
    t_all = time.time()

    def load(key):
        p = wxd / f"{key}.json"
        return [dict(x, t=datetime.fromisoformat(x["t"])) for x in json.load(open(p))] if p.exists() else []

    def one(r):
        rows = load(r)
        if len(rows) < HORIZON_H:
            print(f"  {r}: weather short ({len(rows)} h)", flush=True)
        return r, do_run(a.slug, r, rows, load(f"cal_{r}") if calib else [], a.members, pool, prof, a.variant)

    # several runs in flight at once: one run's slow members no longer idle the pool
    hotspots(a.slug)                          # warm caches before threads share them
    with ThreadPoolExecutor(a.workers) as pool, ThreadPoolExecutor(a.concurrent_runs) as outer:
        for i, f in enumerate(as_completed([outer.submit(one, r) for r in run_list])):
            r, m = f.result()
            cal = (f" cal {m.get('cal_status')} x{m.get('cal_factor')} (obs {m.get('cal_obs_ac')} ac vs replay "
                   f"{m.get('cal_pred_ac')} ac over {m.get('cal_window_h')} h), central by flavour {m.get('flavors')}"
                   + (f", barriers {m.get('barrier_km')} km ({m.get('barrier_map')})" if m.get("barrier_map") else "")
                   if calib else "")
            print(f"  [{i + 1}/{len(run_list)}] {r} {m.get('status')} ign {m.get('ignition_ac')} ac/"
                  f"{m.get('ignition_parts')} parts, det24 {m.get('detections_24h')}, new ac "
                  f"{m.get('new_ac_min')}-{m.get('new_ac_max')} (central {m.get('central_new_ac')}), "
                  f"stalled {m.get('stalled')}, failed {m.get('failed')}, coarse {m.get('coarse_retries')}, "
                  f"{m.get('wall_secs')}s{cal}", flush=True)
    print(f"{a.slug}: done {len(run_list)} runs in {time.time() - t_all:.0f}s", flush=True)


if __name__ == "__main__":
    main()
