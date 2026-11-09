"""Shared helpers for the FARSITE hindcast: paths, run lists, and the dashboard's own
perimeter logic (the repo's analyze.py, so baselines match report.html exactly).

Paths default to farsite/work/ next to this file; override with
  FARSITE_WORK      working dir (data/raw mirror, data/fires/<slug>/..., runs/)
  FARSITE_MANIFEST  bucket manifest.json (default: <work>/manifest.json)
  FARSITE_BIN       TestFARSITE binary (default: <work>/farsite/src/TestFARSITE)
"""
import json
import os
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
FT = Path(os.environ.get("FARSITE_WORK", HERE / "work")).resolve()
DATA = FT / "data"
RUNS = FT / "runs"
FIRES = DATA / "fires"
MANIFEST = Path(os.environ.get("FARSITE_MANIFEST", FT / "manifest.json"))
FARSITE_BIN = Path(os.environ.get("FARSITE_BIN", FT / "farsite" / "src" / "TestFARSITE"))
BUCKET = "https://f005.backblazeb2.com/file/fire-forecast-archive"

sys.path.insert(0, str(HERE.parent))           # the repo root's analyze.py
import analyze  # noqa: E402

analyze.RAW = DATA / "raw"
AC = analyze.SQM_PER_ACRE


def run_dt(run_ts):
    return datetime.strptime(run_ts, "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)


def runs_for(slug, complete_only=True):
    """PyreCast run times for a fire, as listed in the manifest (all five ToA percentiles)"""
    m = json.load(open(MANIFEST))
    out = []
    for k, v in m["runs"].items():
        if v["slug"] != slug:
            continue
        if complete_only and not all(v["files"].get(p, {}).get("ok") for p in analyze.PERCENTILES):
            continue
        out.append(v["run_ts"])
    return sorted(out)


_ps = {}


def perimset(slug):
    if slug not in _ps:
        _ps[slug] = analyze.PerimeterSet(slug)
    return _ps[slug]


def baseline(ps, dt):
    return ps.latest_at_or_before(dt + timedelta(hours=analyze.BASELINE_TOL_H))


def geom_5070(slug, epochms):
    from rasterio.warp import transform_geom
    from shapely.geometry import shape
    from shapely.ops import unary_union
    gj = analyze.load_json(analyze.RAW / "perimeter_archive" / slug / f"{epochms}.geojson")
    return unary_union([shape(transform_geom("EPSG:4326", "EPSG:5070", g)).buffer(0)
                        for g, _ in analyze._parts(gj) if g])
