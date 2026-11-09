"""Hourly point weather for FARSITE runs, as it was forecast at each run's start time.

Hours 0..~44: HRRR, latest extended cycle (00/06/12/18Z) published by the run time
(cycle + 2 h), read from the HRRR-Zarr archive (one 150x150 chunk per variable holds
all 48 forecast hours, ~2 MB). Winds are rotated from grid-relative to earth-relative.
Remainder to +97 h: NBM, latest 00Z/12Z cycle published by the run time (+2 h),
3-hourly GRIB messages fetched by byte range, linear in time.

usage (library): build([(key, run_dt, lon, lat), ...]) -> {key: rows}
"""
import concurrent.futures as cf
import math
import re
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np

HRRR = "https://hrrrzarr.s3.amazonaws.com/sfc"
NBM = "https://noaa-nbm-grib2-pds.s3.amazonaws.com"
HVARS = {"T": ("2m_above_ground", "TMP"), "RH": ("2m_above_ground", "RH"),
         "U": ("10m_above_ground", "UGRD"), "V": ("10m_above_ground", "VGRD"),
         "G": ("surface", "GUST"), "P": ("surface", "APCP_1hr_acc_fcst")}
NVARS = {"T": "TMP:2 m above ground", "RH": "RH:2 m above ground", "WS": "WIND:10 m above ground",
         "WD": "WDIR:10 m above ground", "G": "GUST:10 m above ground"}
HRRR_PROJ = ("+proj=lcc +lat_1=38.5 +lat_2=38.5 +lat_0=38.5 +lon_0=-97.5 +a=6371229 "
             "+b=6371229 +units=m +no_defs")
HORIZON_H = 97
LAG = timedelta(hours=2)
BLEND_H = 3


def _get(url, rng=None, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                raise
            err = e
        except Exception as e:  # noqa: BLE001  transient network errors
            err = e
        import time
        time.sleep(2 * (i + 1))
    raise err


def hrrr_cycle(dt):
    c = (dt - LAG).replace(minute=0, second=0, microsecond=0)
    return c - timedelta(hours=c.hour % 6)


def nbm_cycle(dt):
    c = (dt - LAG).replace(minute=0, second=0, microsecond=0)
    return c - timedelta(hours=c.hour % 12)


# ------------------------------------------------------------------ HRRR-Zarr
_hgrid = None


def hrrr_grid():
    """x/y coordinate vectors of the HRRR grid (projection metres)."""
    global _hgrid
    if _hgrid is None:
        import json
        import numcodecs
        base = f"{HRRR}/20260915/20260915_00z_fcst.zarr/2m_above_ground/TMP"
        out = []
        for ax in ("projection_x_coordinate", "projection_y_coordinate"):
            meta = json.loads(_get(f"{base}/{ax}/.zarray"))
            raw = _get(f"{base}/{ax}/0")
            if meta.get("compressor"):
                raw = numcodecs.get_codec(meta["compressor"]).decode(raw)
            out.append(np.frombuffer(raw, dtype=meta["dtype"]))
        _hgrid = tuple(out)
    return _hgrid


def hrrr_index(lon, lat):
    from rasterio.warp import transform
    xs, ys = transform("EPSG:4326", HRRR_PROJ, [lon], [lat])
    gx, gy = hrrr_grid()
    return int(np.abs(gy - ys[0]).argmin()), int(np.abs(gx - xs[0]).argmin())


def hrrr_series(cycle, pts):
    """{key: {valid_dt: dict(T,RH,WS,WD,G,P)}} for points [(key, lon, lat)]."""
    import numcodecs
    blosc = numcodecs.Blosc()
    d, hh = cycle.strftime("%Y%m%d"), cycle.strftime("%H")
    idx = {k: hrrr_index(lon, lat) for k, lon, lat in pts}
    chunks = {(iy // 150, ix // 150) for iy, ix in idx.values()}
    data = {}

    def fetch(job):
        v, (cy, cx) = job
        lvl, var = HVARS[v]
        raw = _get(f"{HRRR}/{d}/{d}_{hh}z_fcst.zarr/{lvl}/{var}/{lvl}/{var}/0.{cy}.{cx}")
        return job, np.frombuffer(blosc.decode(raw), dtype="<f4").reshape(-1, 150, 150)

    with cf.ThreadPoolExecutor(12) as ex:
        for job, arr in ex.map(fetch, [(v, c) for v in HVARS for c in chunks]):
            data[job] = arr
    out = {}
    for k, lon, lat in pts:
        iy, ix = idx[k]
        c = (iy // 150, ix // 150)
        a = {v: data[(v, c)][:, iy % 150, ix % 150] for v in HVARS}
        alpha = math.radians(math.sin(math.radians(38.5)) * (lon + 97.5))
        ue = math.cos(alpha) * a["U"] + math.sin(alpha) * a["V"]
        ve = -math.sin(alpha) * a["U"] + math.cos(alpha) * a["V"]
        rows = {}
        for i in range(a["T"].shape[0]):
            ws = float(np.hypot(ue[i], ve[i]))
            wd = float((np.degrees(np.arctan2(-ue[i], -ve[i])) + 360) % 360)  # direction from
            rows[cycle + timedelta(hours=i + 1)] = dict(
                T=float(a["T"][i]) - 273.15, RH=float(a["RH"][i]), WS=ws, WD=wd,
                G=float(a["G"][i]), P=max(float(a["P"][i]), 0.0), src="hrrr")
        out[k] = rows
    return out


# ------------------------------------------------------------------ NBM GRIB
def nbm_series(cycle, pts, fhours):
    import rasterio
    from rasterio.warp import transform as wt
    d, hh = cycle.strftime("%Y%m%d"), cycle.strftime("%H")
    base = f"{NBM}/blend.{d}/{hh}/core/blend.t{hh}z.core"

    def one(fh):
        url = f"{base}.f{fh:03d}.co.grib2"
        try:
            idx = _get(url + ".idx").decode().strip().splitlines()
        except urllib.error.HTTPError:
            return fh, None
        offs = [int(l.split(":")[1]) for l in idx]
        picks = {}
        for i, l in enumerate(idx):
            for k, pat in NVARS.items():
                if l.endswith(f":{pat}:{fh} hour fcst:"):
                    picks.setdefault(k, i)
            m = re.search(r":APCP:surface:(\d+)-(\d+) hour acc fcst:$", l)
            if m and int(m.group(2)) == fh:
                w = fh - int(m.group(1))
                if "P" not in picks or w < picks["P"][1]:
                    picks["P"] = (i, w)
        vals = {k: {} for k, _, _ in pts}
        for v, p in picks.items():
            i = p[0] if v == "P" else p
            blob = _get(url, (offs[i], offs[i + 1] - 1))
            with rasterio.MemoryFile(blob) as mf, mf.open() as ds:
                a = ds.read(1)
                xs, ys = wt("EPSG:4326", ds.crs, [p_[1] for p_ in pts], [p_[2] for p_ in pts])
                for (k, _, _), x, y in zip(pts, xs, ys):
                    r, c = ds.index(x, y)
                    val = float(a[r, c])
                    vals[k][v] = val / p[1] if v == "P" else val
        return fh, vals

    out = {k: {} for k, _, _ in pts}
    with cf.ThreadPoolExecutor(12) as ex:
        for fh, vals in ex.map(one, sorted(fhours)):
            if vals is None:
                continue
            for k in vals:
                r = vals[k]
                if {"T", "RH", "WS", "WD", "G"} <= set(r):
                    r.setdefault("P", 0.0)
                    r["src"] = "nbm"
                    out[k][cycle + timedelta(hours=fh)] = r
    return out


# ------------------------------------------------------------------ merge
def _interp(series, t):
    """Linear interpolation of a {dt: row} series at t (direction via vector mean)."""
    ts = sorted(series)
    if t in series:
        return dict(series[t])
    lo = max((s for s in ts if s < t), default=None)
    hi = min((s for s in ts if s > t), default=None)
    if lo is None or hi is None:
        return None
    w = (t - lo) / (hi - lo)
    a, b = series[lo], series[hi]
    r = {k: (1 - w) * a[k] + w * b[k] for k in ("T", "RH", "WS", "G", "P")}
    u = (1 - w) * math.sin(math.radians(a["WD"])) + w * math.sin(math.radians(b["WD"]))
    v = (1 - w) * math.cos(math.radians(a["WD"])) + w * math.cos(math.radians(b["WD"]))
    r["WD"] = (math.degrees(math.atan2(u, v)) + 360) % 360
    r["src"] = a["src"]
    return r


def build_nbm(jobs, step=3):
    """NBM only: the latest 00Z/12Z cycle published by each run time, every `step` hours,
    linear to hourly. One download per cycle serves every point that uses it.
    jobs: [(key, run_dt, lon, lat)] -> {key: [rows with 't' (UTC hour)]}"""
    by_n = {}
    for key, dt, lon, lat in jobs:
        by_n.setdefault(nbm_cycle(dt), []).append((key, dt, lon, lat))
    nser = {}
    for cyc, js in sorted(by_n.items()):
        need = set()
        for key, dt, lon, lat in js:
            f0 = int((dt - cyc).total_seconds() // 3600)
            f1 = int((dt + timedelta(hours=HORIZON_H + 1) - cyc).total_seconds() // 3600) + 1
            need |= {f for f in range(f0 - f0 % step, f1 + step, step) if f >= 1}
            need.add(max(f0, 1))      # run hour itself: f0 < step left the first hour unbracketed (0 rows)
        nser.update(nbm_series(cyc, [(k, lon, lat) for k, _, lon, lat in js], need))
    out = {}
    for key, dt, lon, lat in jobs:
        n = nser.get(key, {})
        t = dt.replace(minute=0, second=0, microsecond=0)
        rows = []
        while t <= dt + timedelta(hours=HORIZON_H + 1):
            r = _interp(n, t) if n else None
            if r is None:
                break
            rows.append(dict(r, t=t))
            t += timedelta(hours=1)
        out[key] = rows
    return out


def build(jobs):
    """HRRR (0-48 h) + NBM: jobs: [(key, run_dt, lon, lat)] -> {key: [rows with 't' (UTC hour)]}."""
    by_h, by_n = {}, {}
    for key, dt, lon, lat in jobs:
        by_h.setdefault(hrrr_cycle(dt), []).append((key, lon, lat))
        by_n.setdefault(nbm_cycle(dt), []).append((key, dt, lon, lat))
    hser = {}
    for cyc, pts in sorted(by_h.items()):
        hser.update(hrrr_series(cyc, pts))
    nser = {}
    for cyc, js in sorted(by_n.items()):
        need = set()
        for key, dt, lon, lat in js:
            hend = hrrr_cycle(dt) + timedelta(hours=48)
            t0 = hend - timedelta(hours=BLEND_H + 3)
            t1 = dt + timedelta(hours=HORIZON_H + 3)
            f0 = int((t0 - cyc).total_seconds() // 3600)
            f1 = int((t1 - cyc).total_seconds() // 3600) + 1
            need |= {f for f in range(f0, f1 + 1) if f % 3 == 0 or f <= 36}
        nser.update(nbm_series(cyc, [(k, lon, lat) for k, _, lon, lat in js], need))
    out = {}
    for key, dt, lon, lat in jobs:
        h, n = hser.get(key, {}), nser.get(key, {})
        hend = max(h) if h else dt
        t = dt.replace(minute=0, second=0, microsecond=0)
        rows = []
        while t <= dt + timedelta(hours=HORIZON_H + 1):
            rh, rn = (h.get(t) if t <= hend else None), _interp(n, t) if n else None
            if rh and rn and t > hend - timedelta(hours=BLEND_H):
                w = (hend - t).total_seconds() / 3600 / BLEND_H     # 1 -> 0 across the seam
                r = {k: w * rh[k] + (1 - w) * rn[k] for k in ("T", "RH", "WS", "G", "P")}
                r["WD"] = rh["WD"] if w >= 0.5 else rn["WD"]
                r["src"] = "blend"
            else:
                r = rh or rn
            if r is None:
                break
            r = dict(r, t=t)
            rows.append(r)
            t += timedelta(hours=1)
        out[key] = rows
    return out
