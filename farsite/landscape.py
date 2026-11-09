"""LANDFIRE layers for a FARSITE .lcp, from the public ArcGIS image services (no account,
unlike LFPS). Fuels/canopy LF2024, topography LF2020, native EPSG:5070 snapped to LANDFIRE's
30 m lattice, so fuels are never resampled."""
import time
import urllib.parse
import urllib.request

import rasterio

SVC = "https://lfps.usgs.gov/arcgis/rest/services"
LAYERS = [  # (band name, service path) in LCP band order
    ("elev", "Landfire_Topo/LF2020_Elev_CONUS"),      # m
    ("slope", "Landfire_Topo/LF2020_SlpD_CONUS"),     # degrees
    ("aspect", "Landfire_Topo/LF2020_Asp_CONUS"),     # degrees azimuth, -1 flat
    ("fbfm40", "Landfire_LF2024/LF2024_FBFM40_CONUS"),
    ("cc", "Landfire_LF2024/LF2024_CC_CONUS"),        # percent
    ("ch", "Landfire_LF2024/LF2024_CH_CONUS"),        # m * 10
    ("cbh", "Landfire_LF2024/LF2024_CBH_CONUS"),      # m * 10
    ("cbd", "Landfire_LF2024/LF2024_CBD_CONUS"),      # kg/m3 * 100
]
X0, Y0, RES = -2362425, 3267405, 30   # LANDFIRE CONUS lattice origin


def fetch(path, bbox, w, h):
    q = urllib.parse.urlencode({
        "bbox": ",".join(str(v) for v in bbox), "bboxSR": 5070, "imageSR": 5070,
        "size": f"{w},{h}", "format": "tiff", "pixelType": "S16", "noData": -9999,
        "interpolation": "RSP_NearestNeighbor", "f": "image"})
    url = f"{SVC}/{path}/ImageServer/exportImage?{q}"
    for i in range(5):                     # a network drop killed one launch with no retry
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                blob = r.read()
            break
        except OSError:
            if i == 4:
                raise
            time.sleep(15 * (i + 1))
    with rasterio.MemoryFile(blob) as mf, mf.open() as ds:
        a = ds.read(1)
        assert ds.width == w and ds.height == h, (ds.width, ds.height)
        assert abs(ds.transform.a - RES) < 1e-6 and abs(ds.transform.c - bbox[0]) < 1e-3, ds.transform
    return a
