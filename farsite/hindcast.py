"""FARSITE hindcast for archived PyreCast runs (v1 spec, pilot on Plaskett + Dome).

For every complete PyreCast run of a fire: start from the same perimeter at the same time,
weather = NBM as published by then, analyst-style rules (active-edge ignition on a
burn-scar landscape, 12:00-20:00 burn periods, gusts during burn periods, no spotting),
a 25-member ensemble (wind speed x, wind direction +, dead fuel moisture +, burn period
length), published as the members at the 10/30/50/70/90th percentile of burned area
(PyreCast's five levels) plus a burn-probability grid.

usage: hindcast.py SLUG [--members N] [--workers N] [--concurrent-runs N] [--limit N]
                         [--only RUN_TS] [--shard I/N] [--prep-only]
--shard I/N runs the I-th of N contiguous slices of the run list (GitHub Actions matrix);
--prep-only builds the landscape and stops.
outputs: data/fires/SLUG/runs/RUN_TS/{toa_p10..p90.tif, burnprob.tif, members.csv, meta.json}
"""
import argparse
import csv
import json
import random
import shutil
import subprocess
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
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
MEMBER_TIMEOUT_S = 600            # p99 of 500 normal members was 307 s
BURN_CENTER = 16.0                              # burn periods centred on 16:00 local
BASE_MOIST = (5, 6, 8, 30, 90)                  # 1h 10h 100h liveH liveW
EDGE_BAND_M, HS_PAD_M, MIN_PART_AC = 150.0, 300.0, 1.0
MARGIN_KM = 15.0
PCTS = ("10", "30", "50", "70", "90")
LCP_OPTS = dict(ELEVATION_UNIT="METERS", SLOPE_UNIT="DEGREES", ASPECT_UNIT="AZIMUTH_DEGREES",
                FUEL_MODEL_OPTION="NO_CUSTOM_AND_NO_FILE", CANOPY_COV_UNIT="PERCENT",
                CANOPY_HT_UNIT="METERS_X_10", CBH_UNIT="METERS_X_10",
                CBD_UNIT="KG_PER_CUBIC_METER_X_100", LINEAR_UNIT="METER")


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


def _polys(g):
    if isinstance(g, Polygon):
        return [g]
    return [p for p in getattr(g, "geoms", []) if isinstance(p, Polygon)]


def ignition(slug, run_ts, rdir):
    """Active-edge ignition polygons + burn-scar LCP for one run. Returns info dict."""
    ps = perimset(slug)
    dt = run_dt(run_ts)
    b = baseline(ps, dt)
    base = geom_5070(slug, b[1])
    prev = ps.latest_at_or_before(dt - timedelta(hours=18))
    growth = Polygon()
    if prev and prev[1] != b[1]:
        growth = base.difference(geom_5070(slug, prev[1])).buffer(-15).buffer(15)
    if slug not in _HS:
        _HS[slug] = analyze.HotspotSet(slug)
    pts = [(t, lon, lat, r) for t, lon, lat, r in _HS[slug].pts if dt - timedelta(hours=24) < t <= dt]
    foot = unary_union([shape(transform_geom("EPSG:4326", "EPSG:5070", mapping(Point(lon, lat))))
                        .buffer(r + HS_PAD_M) for _, lon, lat, r in pts]) if pts else Polygon()
    band = base.difference(base.buffer(-EDGE_BAND_M))
    active = unary_union([growth, band.intersection(foot)]).intersection(base)
    parts = [Polygon(p.exterior) for p in _polys(active) if p.area / AC >= MIN_PART_AC]
    info = dict(baseline_ms=b[1], baseline_time=f"{b[0]:%Y-%m-%dT%H:%MZ}", baseline_ac=round(base.area / AC),
                growth_ac=round(growth.area / AC), detections_24h=len(pts),
                ignition_ac=round(sum(p.area for p in parts) / AC), ignition_parts=len(parts),
                weather_point=list(warp_transform("EPSG:5070", "EPSG:4326",
                                                  [base.centroid.x], [base.centroid.y])))
    info["weather_point"] = [info["weather_point"][0][0], info["weather_point"][1][0]]
    if not parts:
        return info
    gj = rdir / "ign.geojson"
    json.dump({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"id": i}, "geometry": mapping(p)} for i, p in enumerate(parts)]},
        open(gj, "w"))
    subprocess.run(["ogr2ogr", "-f", "ESRI Shapefile", str(rdir / "ign.shp"), str(gj), "-a_srs",
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
    with rasterio.open(rdir / "scar.tif", "w", **prof) as ds:
        ds.write(stack)
    lat = json.load(open(FIRES / slug / "grid.json"))["lat"]
    rio_copy(rdir / "scar.tif", rdir / "scar.lcp", driver="LCP", LATITUDE=str(int(round(lat))), **LCP_OPTS)
    (rdir / "scar.tif").unlink()
    return info


# ---------------------------------------------------------------- members
def members(slug, run_ts, n):
    rng = random.Random(zlib.crc32(f"{slug}/{run_ts}".encode()))
    out = [dict(m=0, ws=1.0, wd=0.0, dfm=0.0, bp=8)]
    for i in range(1, n):
        out.append(dict(m=i, ws=round(min(max(np.exp(rng.gauss(0, 0.25)), 0.6), 1.6), 3),
                        wd=round(rng.gauss(0, 15), 1), dfm=round(rng.uniform(-2, 2), 2),
                        bp=rng.choice((6, 7, 8, 9, 10))))
    return out


def input_text(start, end, rows, elev_ft, mem):
    s, e = start.astimezone(TZ), end.astimezone(TZ)
    bp0, bp1 = BURN_CENTER - mem["bp"] / 2, BURN_CENTER + mem["bp"] / 2
    hhmm = lambda h: int(h) * 100 + int(round((h % 1) * 60))
    periods, d = [], s.date()
    while d <= e.date():
        p0 = datetime.combine(d, datetime.min.time(), TZ) + timedelta(hours=bp0)
        p1 = datetime.combine(d, datetime.min.time(), TZ) + timedelta(hours=bp1) - timedelta(minutes=1)
        if p1 > s:
            p0 = max(p0, s)
            if p0 < min(p1, e):
                periods.append(f"{d:%m %d} {p0.hour * 100 + p0.minute:04d} {p1.hour * 100 + p1.minute:04d}")
        d += timedelta(days=1)
    lines = []
    for r in rows:
        lt = r["t"].astimezone(TZ)
        in_bp = hhmm(bp0) <= lt.hour * 100 < hhmm(bp1)
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
FARSITE_DISTANCE_RES: 30.0
FARSITE_PERIMETER_RES: 60.0
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


def run_member(rdir, mem, txt, start, end):
    md = rdir / f"m{mem['m']:02d}"
    md.mkdir(exist_ok=True)
    (md / "in.input").write_text(txt)
    (md / "cmd.txt").write_text(f"{rdir / 'scar.lcp'} {md / 'in.input'} {rdir / 'ign.shp'} 0 {md / 'o'} 0\n")
    asc = md / "o_ArrivalTime.asc"
    t0, err = time.time(), ""
    # FARSITE occasionally blows up on perimeter vertices (40+ CPU-min vs a 22 s median).
    # Kill at MEMBER_TIMEOUT_S and retry once at half the perimeter/distance resolution.
    for attempt, t in enumerate((txt, txt.replace("FARSITE_DISTANCE_RES: 30.0", "FARSITE_DISTANCE_RES: 60.0")
                                       .replace("FARSITE_PERIMETER_RES: 60.0", "FARSITE_PERIMETER_RES: 120.0"))):
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


# ---------------------------------------------------------------- one run
def do_run(slug, run_ts, wxrows, n_members, pool, lcp_profile):
    rdir = FIRES / slug / "runs" / run_ts
    if (rdir / "meta.json").exists():
        return json.load(open(rdir / "meta.json"))
    rdir.mkdir(parents=True, exist_ok=True)
    dt = run_dt(run_ts)
    t0 = time.time()
    info = ignition(slug, run_ts, rdir)
    meta = dict(slug=slug, run_ts=run_ts, **info, members=n_members)
    start, end = dt, dt + timedelta(hours=HORIZON_H + 1)
    shape_ = (lcp_profile["height"], lcp_profile["width"])
    if not info["ignition_parts"] or not wxrows:
        meta["status"] = "no_active_edge" if not info["ignition_parts"] else "no_weather"
        empty = np.zeros(shape_, "float32")
        for p in PCTS:
            _write(rdir / f"toa_p{p}.tif", empty, lcp_profile)
        json.dump(meta, open(rdir / "meta.json", "w"), indent=1)
        return meta
    with rasterio.open(FIRES / slug / "lcp.tif") as ds:
        x, y = warp_transform("EPSG:4326", "EPSG:5070", [info["weather_point"][0]], [info["weather_point"][1]])
        r, c = ds.index(x[0], y[0])
        el = ds.read(1)[max(r - 40, 0):r + 40, max(c - 40, 0):c + 40]
    elev_ft = float(el[el > 0].mean() if (el > 0).any() else 0) * 3.28084
    mems = members(slug, run_ts, n_members)
    futs, periods_by_m = [], {}
    for mem in mems:
        txt, periods = input_text(start, end, wxrows, elev_ft, mem)
        periods_by_m[mem["m"]] = periods
        futs.append(pool.submit(run_member, rdir, mem, txt, start, end))
    toas, rows, fails = {}, [], 0
    base_m = None
    for f in as_completed(futs):
        mem, a, secs, err = f.result()
        if a is None:
            fails += 1
            rows.append(dict(mem, secs=round(secs, 1), new_ac="", stalled="", error=err.replace("\n", " ")))
            continue
        h = np.where(a != -9999, a / 60.0, np.nan).astype("float32")   # hours since run (start = run)
        toas[mem["m"]] = h
        new = np.isfinite(h) & (h > 0) & (h <= HORIZON_H)
        stalled = _stalled(h, dt, periods_by_m[mem["m"]])
        rows.append(dict(mem, secs=round(secs, 1), new_ac=round(new.sum() * 900 / AC), stalled=int(stalled), error=""))
    rows.sort(key=lambda r: r["m"])
    with open(rdir / "members.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    ok = [r for r in rows if r["new_ac"] != ""]
    if ok:
        order = sorted(ok, key=lambda r: (r["new_ac"], r["m"]))
        pick = {}
        for p in PCTS:
            i = min(int(round(int(p) / 100 * (len(order) - 1))), len(order) - 1)
            pick[p] = order[i]["m"]
            _write(rdir / f"toa_p{p}.tif", np.nan_to_num(np.where(toas[pick[p]] <= 0, 1e-3, toas[pick[p]]), nan=0.0), lcp_profile)
        bp = np.mean([np.isfinite(t) & (t <= HORIZON_H) for t in toas.values()], axis=0).astype("float32")
        _write(rdir / "burnprob.tif", bp, lcp_profile)
        meta.update(status="ok", pick=pick, failed=fails, stalled=sum(r["stalled"] for r in ok),
                    coarse_retries=sum(r["coarse"] for r in rows),
                    member_secs=round(sum(r["secs"] for r in rows), 1),
                    new_ac_min=order[0]["new_ac"], new_ac_max=order[-1]["new_ac"],
                    central_new_ac=next((r["new_ac"] for r in ok if r["m"] == 0), ""))
    else:
        meta.update(status="all_failed", failed=fails)
    meta["wall_secs"] = round(time.time() - t0, 1)
    for md in rdir.glob("m[0-9][0-9]"):
        shutil.rmtree(md)
    for ext in ("lcp", "lcp.aux.xml"):
        (rdir / f"scar.{ext}").unlink(missing_ok=True)
    json.dump(meta, open(rdir / "meta.json", "w"), indent=1)
    return meta


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
    ap.add_argument("--members", type=int, default=25)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--concurrent-runs", type=int, default=8)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--only")
    ap.add_argument("--shard", help="I/N: the I-th (0-based) of N contiguous slices")
    ap.add_argument("--prep-only", action="store_true")
    a = ap.parse_args()
    run_list = runs_for(a.slug)
    if a.only:
        run_list = [r for r in run_list if r == a.only]
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
    # weather for every run at its baseline centroid, cached per run
    wxd = d / "wx"
    wxd.mkdir(exist_ok=True)
    todo = []
    ps = perimset(a.slug)
    for r in run_list:
        if not (wxd / f"{r}.json").exists():
            g = geom_5070(a.slug, baseline(ps, run_dt(r))[1]).centroid
            x, y = warp_transform("EPSG:5070", "EPSG:4326", [g.x], [g.y])
            todo.append((r, run_dt(r), x[0], y[0]))
    if todo:
        t0 = time.time()
        for r, rows in wx.build_nbm(todo).items():
            json.dump([dict(x, t=x["t"].isoformat()) for x in rows], open(wxd / f"{r}.json", "w"))
        print(f"{a.slug}: weather for {len(todo)} runs in {time.time() - t0:.0f}s", flush=True)
    t_all = time.time()

    def one(r):
        rows = [dict(x, t=datetime.fromisoformat(x["t"])) for x in json.load(open(wxd / f"{r}.json"))]
        if len(rows) < HORIZON_H:
            print(f"  {r}: weather short ({len(rows)} h)", flush=True)
        return r, do_run(a.slug, r, rows, a.members, pool, prof)

    # several runs in flight at once: one run's slow members no longer idle the pool
    if a.slug not in _HS:                     # warm caches before threads share them
        _HS[a.slug] = analyze.HotspotSet(a.slug)
    with ThreadPoolExecutor(a.workers) as pool, ThreadPoolExecutor(a.concurrent_runs) as outer:
        for i, f in enumerate(as_completed([outer.submit(one, r) for r in run_list])):
            r, m = f.result()
            print(f"  [{i + 1}/{len(run_list)}] {r} {m.get('status')} ign {m.get('ignition_ac')} ac/"
                  f"{m.get('ignition_parts')} parts, det24 {m.get('detections_24h')}, new ac "
                  f"{m.get('new_ac_min')}-{m.get('new_ac_max')} (central {m.get('central_new_ac')}), "
                  f"stalled {m.get('stalled')}, failed {m.get('failed')}, coarse {m.get('coarse_retries')}, {m.get('wall_secs')}s", flush=True)
    print(f"{a.slug}: done {len(run_list)} runs in {time.time() - t_all:.0f}s", flush=True)


if __name__ == "__main__":
    main()
