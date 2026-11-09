# FARSITE hindcast (pilot)

A second forecast model next to PyreCast: for every archived PyreCast run of a fire, FARSITE
starts from the same perimeter at the same time, with only inputs that were public at run time.
It runs as a 25-member ensemble, published at PyreCast's five percentile levels.

- **Engine:** FARSITE from [mbedward/farsite](https://github.com/mbedward/farsite) PR #1 (public domain), built on Linux.
- **Landscape:** LANDFIRE LF2024 fuels and canopy plus LF2020 topography, from the public ArcGIS image services, 30 m in EPSG:5070.
- **Weather:** NOAA NBM as published by run time, from the `noaa-nbm-grib2-pds` bucket.
- **Analyst-style setup:**
  - Ignition is the active edge only: growth since the perimeter 18 h earlier, plus the edge band near detections from the last 24 h.
  - Fuel is set to non-burnable inside the rest of the burn scar.
  - Burn periods are centred on 16:00 local.
  - Spotting is off.
- **Ensemble:**
  - Member 0 uses the forecast as published.
  - Each other member scales wind speed by a factor of ×0.6–1.6, shifts wind direction by about ±15°, shifts dead fuel moisture by up to ±2 points, and uses burn periods of 6–10 h.
  - Members are seeded by fire and run, so reruns are reproducible.
- **Resolution:** 30 m distance and 60 m perimeter resolution. A run whose active-edge ignition is over 3,000 ac (FARSITE's cost grows with the active perimeter) starts at 60/120 m instead (`res_m` in `meta.json`).
- **Hung members:** a member that runs past its timeout (600 s locally, 1,500 s on Actions) is retried once at half that resolution (`coarse_retries` in `meta.json`).

## Versions

`hindcast.py --variant` picks the version. Each one writes its own output directory.

| Variant | Output directories | What it adds |
|---|---|---|
| `base` | `runs/` | The v1 pilot: gust speed is the wind during burn periods, no calibration. |
| `cal` | `runs-cal/`, `runs-calr/` | Sustained wind all day, plus calibration to the previous day's growth. |
| `calbar` | `runs-calbar/`, `runs-calbarr/` | `cal`, plus completed control lines from the IAP ops maps as FARSITE barriers. |

**How calibration works**

1. Find two official perimeters about a day apart: the run's baseline perimeter B, and the latest perimeter P0 from 20–54 h before it.
2. Light P0's active edge with 9 members and run them to B's time, using the NBM published by P0's time.
3. Compare their median new area with B beyond P0.
4. The ratio becomes a factor on burn-period length, meaning hours of active spread per day.
   - Analysts tune this lever to observed fire activity.
   - Unlike FARSITE's ROS adjustment file, it also slows crown fire.
5. A second replay refines the factor. It is clamped to 0.05–2.5.
6. If the perimeter didn't grow, the replay is skipped and the factor is the floor.

The two flavours share the same replay:

- `runs-cal/` holds the factor for the whole 96 h.
- `runs-calr/` halves the factor's effect, in log space, every 4 burn periods. This eases the forecast back toward uncalibrated at longer leads.

Satellite detections are not used to measure growth. Footprints of 375 m VIIRS pixels along a slow edge read as about 5 times the perimeter growth on Dome.

**Control lines** come from `interventions.py`, which writes `interventions/<slug>.geojson`.

- It reads the IAP operations maps that Incibrief mirrors from the NIFC incident-map FTP, in the `responder-debrief-data` bucket.
- It reads them from the GeoPDF's vector structure, not the image:
  - The map's own legend gives the symbols.
  - Completed dozer lines are chains of marker-font hatch glyphs.
  - Hand, road-as-line and mixed lines are the strokes their H, R and M glyphs sit on.
  - GDAL's 72 dpi geotransform turns page points into coordinates.
- On Plaskett, the map's fire edge matches the official perimeter to a median of 1 m.
- Each run uses lines from the latest vector map published (FTP time) by its start, within 5 km of the fire. The replay uses lines published by its own start.
- Maps flattened to a single image carry no vectors and are skipped. That covers every Dome 2026 map and Plaskett's 9/24 and 9/26 maps.

## Running it

GitHub Actions (`.github/workflows/farsite-hindcast.yml`) runs `targets.json` on every push to
the `farsite-hindcast` branch that touches the hindcast code, the targets, or `interventions/`:

```json
{"runs": [{"slug": "ca-dome", "variant": "cal"}, {"slug": "ca-plaskett", "variant": "calbar"}],
 "members": 25, "shards": 70, "concurrent_runs": 1}
```

- **Entries:** an entry can take `"only": [run_ts, ...]` to run just those runs.
- **Shards:** each shard job runs the shard's slice of every entry.
- **Concurrency:** at most 16 jobs run at once, which leaves runners for other projects. One workflow run goes at a time, and a new push waits for the run in progress.

The result is the run's `farsite-results` artifact, laid out as:

```
data/fires/<slug>/runs[-flavour]/<run_ts>/toa_p{10,30,50,70,90}.tif   arrival, hours since run_ts (0 = no data)
                                burnprob.tif                share of members burned by 96 h
                                members.csv, meta.json      per-member inputs, timings, status
```

To run locally, use the same scripts. The working directory defaults to `farsite/work/`:

```bash
python farsite/fetch_archive.py ca-plaskett            # manifest, ToA tifs, perimeters, hotspots
python farsite/interventions.py ca-plaskett             # control lines (pymupdf, scipy, gdalinfo)
python farsite/hindcast.py ca-plaskett --variant calbar --members 25   # needs TestFARSITE at work/farsite/src/
python farsite/score2.py ca-plaskett && python farsite/plot2.py
```

`score2.py` scores both models at the time resolution the observations allow:

- every perimeter at its exact lead time,
- every satellite pass, measured as capture,
- 1-hour lead bins,
- plus the dashboard's fixed horizons.

Growth is measured beyond a common start for both models: the baseline perimeter, plus
PyreCast's start footprint (`analyze.RunStart`), plus FARSITE's ignition. Every FARSITE
version present is scored side by side. Summaries use only rows where every version has a
forecast.
