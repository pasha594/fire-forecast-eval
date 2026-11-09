"""Completed control lines from IAP operations maps, read from the GeoPDF's vector structure
(never from the rendered image).

Source: Incibrief's mirror of the NIFC incident-map FTP (the responder-debrief-data bucket).
state/state.json indexes every incident folder's files with the FTP publish time ("lm"), which
is when a map became public, so a run only uses maps published by its start.

For each ops-type map of a fire (ops, ops-contingency, stratops, branch):
  - ArcGIS exports keep map symbols as vectors (PyMuPDF get_drawings, rawdict text). PDFs
    flattened to a single image have no drawings; they are listed as raster_only and skipped.
  - The map's own legend decodes the symbology. For each "Completed ... Line" row, the symbol
    left of the label gives marker glyphs (font + character + colour) and filled glyph outlines
    (matched by path shape, so H / R / M letters are told apart at any scale).
  - Completed hand line / road as line / mixed construction: the thin black strokes those
    letter glyphs sit on. Completed dozer line: the chain of hatch glyphs (no stroke under them),
    joined by a minimum spanning tree with an edge cap.
  - Page points -> map CRS with GDAL's geotransform at 72 dpi (1 px = 1 pt) -> EPSG:5070.
  - QA: the map's own fire-edge strokes are compared with the nearest official perimeter.

Writes farsite/interventions/<slug>.geojson (EPSG:5070 LineStrings: kind, map, product, op_date,
published) with a "maps" list (status and counts per PDF).
usage: interventions.py SLUG [SLUG ...]      (needs pymupdf, scipy, gdalinfo)
"""
import json
import re
import subprocess
import sys
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

import numpy as np
import pymupdf
from rasterio.crs import CRS
from rasterio.warp import transform as warp_transform
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString, Point, mapping, shape
from shapely.ops import linemerge, unary_union
from shapely import wkt as shp_wkt

from common import DATA, analyze, geom_5070, perimset

RDD = "https://f005.backblazeb2.com/file/responder-debrief-data"
OUT = Path(__file__).resolve().parent / "interventions"
PRODUCTS = re.compile(r"^(ops|ops-contingency|opsir|stratops|branch)", re.I)
KINDS = {"dozer": "dozer", "hand": "hand", "mixed construction": "mixed", "road as": "road"}
LEGEND = re.compile(r"(?i)^\s*completed (dozer|hand|mixed construction|road as) line\s*$")
MARKER_FONT = re.compile(r"(?i)^esri")    # ArcGIS marker-symbol fonts
THIN_MAX_PT = 2.0          # line strokes under letter glyphs (0.78 pt on 9/5 Plaskett, 1.56 on 9/1)
ON_LINE_PT = 6.0           # a letter glyph's centre within this of a stroke labels it
DEDUPE_PT = 3.0            # hatch glyph drawn as outline + text: one mark


def _merge(segs):
    segs = [g for g in segs if g.length > 0]
    if not segs:
        return None
    u = unary_union(segs)
    return u if isinstance(u, LineString) else linemerge(u)


def _get(url, path):
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=120) as r:
            path.write_bytes(r.read())
    return path


def maps_for(slug, state):
    """ops-type PDFs of a fire: [(filename, url, product, op_date, published UTC)]"""
    name = slug.split("-", 1)[1]
    out = []
    for key, inc in state["incidents"].items():
        if inc.get("fire_slug") != name or "/2026_" not in key:
            continue
        for fk, f in inc["files"].items():
            fn = fk.split("/")[-1]
            if not fk.startswith("products/") or not fn.lower().endswith(".pdf") or not PRODUCTS.match(fn):
                continue
            prefix = (state.get("tiled", {}).get(f.get("sha16"), {}) or {}).get("prefix") or name
            m = re.search(r"_(\d{4})(day|night)?(?:[-_]|\.pdf)", fn)
            out.append(dict(filename=fn, url=f"{RDD}/raw/incidents/{prefix}/{fk.replace(' ', '_')}",
                            product=fn.split("_")[0].lower(), op_date=fk.split("/")[1],
                            published=parsedate_to_datetime(f["lm"]).astimezone(timezone.utc)))
    return sorted(out, key=lambda x: x["published"])


def georef(pdf):
    """GDAL's 72-dpi geotransform (page points -> map coords) and CRS"""
    j = json.loads(subprocess.run(["gdalinfo", "-json", "--config", "GDAL_PDF_DPI", "72", str(pdf)],
                                  capture_output=True, text=True, check=True).stdout)
    gt = j.get("geoTransform")
    wkt = (j.get("coordinateSystem") or {}).get("wkt")
    if not gt or not wkt or gt == [0, 1, 0, 0, 0, 1]:
        return None
    neat = (j.get("metadata", {}).get("", {}) or {}).get("NEATLINE")
    return dict(gt=gt, crs=CRS.from_wkt(wkt), neatline=neat, size=j["size"])


def _col(c):
    return tuple(round(x, 2) for x in c) if c else None


def _sig(d):
    return "".join(it[0] for it in d["items"])


def _center(r):
    return ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)


def _chars(page):
    for b in page.get_text("rawdict")["blocks"]:
        for ln in b.get("lines", []):
            for s in ln["spans"]:
                for ch in s["chars"]:
                    yield s, ch


def decode_legend(page, drawings):
    """{kind: {"text": {(font, char, color)}, "outline": {(type, fill, sig)}}} from the legend rows"""
    labels = []
    for b in page.get_text("dict")["blocks"]:
        for ln in b.get("lines", []):
            for sp in ln["spans"]:
                m = LEGEND.match(sp["text"])
                if m:
                    labels.append((KINDS[m.group(1).lower()], pymupdf.Rect(sp["bbox"])))
    styles = {}
    for kind, r in labels:
        yc = (r.y0 + r.y1) / 2
        zone = pymupdf.Rect(r.x0 - 75, yc - 8, r.x0 - 2, yc + 8)
        st = styles.setdefault(kind, {"text": set(), "outline": set()})
        for s, ch in _chars(page):      # marker fonts only (not a neighbouring column's label text)
            if ch["c"].strip() and MARKER_FONT.match(s["font"]) \
                    and zone.contains(pymupdf.Point(*_center(pymupdf.Rect(ch["bbox"])))):
                st["text"].add((s["font"], ch["c"], s["color"]))
        for d in drawings:
            if d["type"] in ("f", "fs") and _col(d.get("fill")) == (0.0, 0.0, 0.0) and zone.contains(d["rect"]):
                st["outline"].add((d["type"], _sig(d)))
    return styles, [r for _, r in labels]


def extract(pdf, meta):
    page = pymupdf.open(pdf)[0]
    drawings = page.get_drawings()
    if len(drawings) < 50:
        return dict(status="raster_only"), []
    g = georef(pdf)
    if not g:
        return dict(status="no_georef"), []
    styles, label_rects = decode_legend(page, drawings)
    if not styles:
        return dict(status="no_legend"), []
    gt = g["gt"]
    inv = np.linalg.inv(np.array([[gt[1], gt[2]], [gt[4], gt[5]]]))
    to_page = lambda X, Y: tuple(inv @ np.array([X - gt[0], Y - gt[3]]))
    frame = None
    if g["neatline"]:
        nl = shp_wkt.loads(g["neatline"])
        frame = LineString([to_page(x, y) for x, y in nl.exterior.coords]).convex_hull
    inside = lambda pt: frame is None or frame.contains(Point(pt))
    legend_box = unary_union([LineString([(r.x0 - 80, r.y0 - 10), (r.x1 + 10, r.y1 + 10)]).envelope
                              for r in label_rects])
    in_map = lambda pt: inside(pt) and not legend_box.contains(Point(pt))

    # glyph marks per kind: marker-font characters and filled outlines matching the legend
    marks = defaultdict(list)
    text_kind = {t: k for k, st in styles.items() for t in st["text"]}
    for s, ch in _chars(page):
        k = text_kind.get((s["font"], ch["c"], s["color"]))
        c = _center(pymupdf.Rect(ch["bbox"]))
        if k and in_map(c):
            marks[k].append(c)
    out_kind = {o: k for k, st in styles.items() for o in st["outline"]}
    for d in drawings:
        if d["type"] in ("f", "fs") and _col(d.get("fill")) == (0.0, 0.0, 0.0):
            k = out_kind.get((d["type"], _sig(d)))
            c = _center(d["rect"])
            if k and in_map(c) and d["rect"].width < 40:
                marks[k].append(c)

    lines = []
    # hand / road / mixed: thin black solid strokes in the map layer, labelled by the letters on them
    letters = [(k, p) for k in ("hand", "road", "mixed") for p in marks.get(k, [])]
    if letters:
        tree = cKDTree([p for _, p in letters])
        for d in drawings:
            if (d["type"] != "s" or _col(d.get("color")) != (0.0, 0.0, 0.0) or (d.get("width") or 9) > THIN_MAX_PT
                    or not d.get("layer") or (d.get("dashes") or "[] 0").split("]")[0].strip("[ ")):
                continue
            geom = _merge([LineString([it[1], it[-1]]) for it in d["items"] if it[0] in ("l", "c")])
            if geom is None:
                continue
            for part in getattr(geom, "geoms", [geom]):   # one path can hold many lines: label each
                pts = np.array(part.coords)
                near = Counter()
                for i in set(j for js in tree.query_ball_point(pts, part.length + ON_LINE_PT) for j in js) \
                        if len(pts) < 3 else set(j for js in tree.query_ball_point(pts, ON_LINE_PT * 3) for j in js):
                    k, p = letters[i]
                    if part.distance(Point(p)) <= ON_LINE_PT:
                        near[k] += 1
                if near and in_map(pts[0]) and in_map(pts[-1]):
                    lines.append((near.most_common(1)[0][0], part, sum(near.values())))
    # dozer: chain of hatch marks
    pts = marks.get("dozer", [])
    if len(pts) >= 2:
        P = np.array(pts)
        keep = []
        for i in np.argsort(P[:, 0]):
            if all(np.hypot(*(P[i] - P[j])) > DEDUPE_PT for j in keep[-50:]):
                keep.append(i)
        P = P[keep]
        if len(P) >= 2:
            tree = cKDTree(P)
            nn = tree.query(P, k=2)[0][:, 1]
            cap = 2.0 * float(np.median(nn))
            pairs = np.array(sorted(tree.query_pairs(cap)))
            if len(pairs):
                from scipy.sparse import coo_matrix
                w = np.hypot(*(P[pairs[:, 0]] - P[pairs[:, 1]]).T)
                mst = minimum_spanning_tree(coo_matrix((w, (pairs[:, 0], pairs[:, 1])), shape=(len(P), len(P))))
                a, b = mst.nonzero()
                geom = _merge([LineString([P[i], P[j]]) for i, j in zip(a, b)])
                for gpart in getattr(geom, "geoms", [geom]) if geom is not None else []:
                    lines.append(("dozer", gpart, len(gpart.coords)))

    # page -> map CRS -> EPSG:5070
    feats = []
    for kind, geom, n in lines:
        for part in getattr(geom, "geoms", [geom]):
            xy = np.array(part.coords)
            X = gt[0] + xy[:, 0] * gt[1] + xy[:, 1] * gt[2]
            Y = gt[3] + xy[:, 0] * gt[4] + xy[:, 1] * gt[5]
            x5, y5 = warp_transform(g["crs"], "EPSG:5070", X.tolist(), Y.tolist())
            ls = LineString(zip(x5, y5))
            if ls.length > 0:
                feats.append(dict(kind=kind, geom=ls, glyphs=n))
    # QA: the map's fire-edge strokes vs the official perimeter nearest its publish time
    edge = [d for d in drawings if d["type"] == "s" and (d.get("width") or 0) >= 4
            and _col(d.get("color")) in ((0.9, 0.0, 0.0), (0.0, 0.0, 0.0)) and d.get("layer")]
    qa = None
    if edge:
        xy = np.array([it[1] for d in edge for it in d["items"] if it[0] == "l"])[::5]
        xy = np.array([p for p in xy if in_map(p)])
        if len(xy):
            X = gt[0] + xy[:, 0] * gt[1] + xy[:, 1] * gt[2]
            Y = gt[3] + xy[:, 0] * gt[4] + xy[:, 1] * gt[5]
            x5, y5 = warp_transform(g["crs"], "EPSG:5070", X.tolist(), Y.tolist())
            qa = list(zip(x5, y5))
    info = dict(status="vector", legend={k: dict(text=len(v["text"]), outline=len(v["outline"])) for k, v in styles.items()},
                marks={k: len(v) for k, v in marks.items()},
                lines=dict(Counter(f["kind"] for f in feats)),
                km=dict((k, round(sum(f["geom"].length for f in feats if f["kind"] == k) / 1000, 1))
                        for k in sorted({f["kind"] for f in feats})))
    return dict(info, _edge=qa), feats


def main():
    st = json.load(urllib.request.urlopen(f"{RDD}/state/state.json", timeout=120))
    OUT.mkdir(exist_ok=True)
    for slug in sys.argv[1:]:
        ps = perimset(slug)
        feats, maps = [], []
        for m in maps_for(slug, st):
            pdf = _get(m["url"], DATA / "raw" / "iap" / slug / m["filename"])
            info, fs = extract(pdf, m)
            edge = info.pop("_edge", None)
            if edge:
                near = ps.latest_at_or_before(m["published"])
                if near:
                    per = geom_5070(slug, near[1]).boundary
                    d = np.array([per.distance(Point(p)) for p in edge])
                    info["qa_edge_to_perimeter_m"] = dict(median=round(float(np.median(d))), p90=round(float(np.percentile(d, 90))),
                                                          perimeter_time=f"{near[0]:%Y-%m-%dT%H:%MZ}")
            row = dict(filename=m["filename"], product=m["product"], op_date=m["op_date"],
                       published=m["published"].strftime("%Y-%m-%dT%H:%M:%SZ"), **info)
            maps.append(row)
            print(f"{slug} {m['filename'][:60]:60s} {row['published']} {info['status']} "
                  f"{info.get('lines', '')} {info.get('km', '')} {info.get('qa_edge_to_perimeter_m', '')}", flush=True)
            for f in fs:
                feats.append({"type": "Feature", "geometry": mapping(f["geom"]),
                              "properties": dict(kind=f["kind"], status="completed", glyphs=f["glyphs"],
                                                 map=m["filename"], product=m["product"], op_date=m["op_date"],
                                                 published=row["published"])})
        json.dump({"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "EPSG:5070"}},
                   "maps": maps, "features": feats}, open(OUT / f"{slug}.geojson", "w"))
        print(f"{slug}: {len(feats)} line features from {sum(x['status'] == 'vector' for x in maps)} vector maps "
              f"({sum(x['status'] == 'raster_only' for x in maps)} raster-only) -> {OUT / (slug + '.geojson')}")


if __name__ == "__main__":
    main()
