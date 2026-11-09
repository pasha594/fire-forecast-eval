"""Mirror one or more fires' archive objects (ToA tifs, perimeters, hotspots) into data/raw.
usage: fetch_archive.py [--no-toa] SLUG [SLUG ...]
Downloads the bucket's manifest.json to FARSITE_MANIFEST first if it isn't there.
The hindcast needs only perimeters + hotspots (--no-toa); scoring also needs the ToA tifs."""
import concurrent.futures as cf
import json
import sys
import time
import urllib.request
from collections import Counter

from common import BUCKET, DATA, MANIFEST


def jobs(manifest, slugs, toa=True):
    m = json.load(open(manifest))
    for k, v in m["runs"].items():
        if toa and v["slug"] in slugs:
            for p in ("10", "30", "50", "70", "90"):
                if v["files"].get(p, {}).get("ok"):
                    yield f"forecast_archive/{k}/{p}.tif"
    for s in slugs:
        yield f"perimeter_archive/{s}/index.json"
        for e in m["perimeters"].get(s, []):
            yield f"perimeter_archive/{s}/{e}.geojson"
        for d in m["hotspots"].get(s, []):
            yield f"hotspot_archive/{s}/{d}.geojson"


def get(rel):
    dst = DATA / "raw" / rel
    if dst.exists() and dst.stat().st_size > 0 and not rel.endswith("index.json"):
        return "cached"
    dst.parent.mkdir(parents=True, exist_ok=True)
    for i in range(4):                    # B2 throws the odd 500/503; 404 is final
        try:
            with urllib.request.urlopen(f"{BUCKET}/{rel}", timeout=120) as r:
                dst.write_bytes(r.read())
            return "ok"
        except Exception as e:
            if getattr(e, "code", None) == 404 or i == 3:
                return f"err {getattr(e, 'code', e)}"
            time.sleep(5 * (i + 1))


if __name__ == "__main__":
    args = sys.argv[1:]
    toa = "--no-toa" not in args
    slugs = {s for s in args if not s.startswith("--")}
    if not MANIFEST.exists():
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(f"{BUCKET}/manifest.json", timeout=120) as r:
            MANIFEST.write_bytes(r.read())
    js = list(jobs(MANIFEST, slugs, toa))
    with cf.ThreadPoolExecutor(16) as ex:
        res = Counter(ex.map(get, js))
    print(len(js), res)
    if any(k.startswith("err") for k in res):
        sys.exit(1)
