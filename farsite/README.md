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
  - Burn periods are centred on 16:00 local, with gust speed during burn periods.
  - Spotting is off.
- **Ensemble:**
  - Member 0 uses the forecast as published.
  - Each other member scales wind speed by a factor of ×0.6–1.6, shifts wind direction by about ±15°, shifts dead fuel moisture by up to ±2 points, and uses burn periods of 6–10 h.
  - Members are seeded by fire and run, so reruns are reproducible.
- **Resolution:** 30 m distance and 60 m perimeter resolution. A run whose active-edge ignition is over 3,000 ac (FARSITE's cost grows with the active perimeter) starts at 60/120 m instead (`res_m` in `meta.json`).
- **Hung members:** a member that runs past its timeout (600 s locally, 1,500 s on Actions) is retried once at half that resolution (`coarse_retries` in `meta.json`).

## Running it

GitHub Actions (`.github/workflows/farsite-hindcast.yml`) runs `targets.json` on every push to
the `farsite-hindcast` branch:

```json
{"slugs": ["ca-plaskett"], "members": 25, "shards": 52, "concurrent_runs": 1, "limit": null}
```

The result is the run's `farsite-results` artifact, laid out as:

```
data/fires/<slug>/runs/<run_ts>/toa_p{10,30,50,70,90}.tif   arrival, hours since run_ts (0 = no data)
                                burnprob.tif                share of members burned by 96 h
                                members.csv, meta.json      per-member inputs, timings, status
```

To run locally, use the same scripts. The working directory defaults to `farsite/work/`:

```bash
python farsite/fetch_archive.py ca-plaskett            # manifest, ToA tifs, perimeters, hotspots
python farsite/hindcast.py ca-plaskett --members 25    # needs TestFARSITE at work/farsite/src/
python farsite/score2.py ca-plaskett && python farsite/plot2.py
```

`score2.py` scores both models at the time resolution the observations allow:

- every perimeter at its exact lead time,
- every satellite pass, measured as capture,
- 1-hour lead bins,
- plus the dashboard's fixed horizons.

Growth is measured beyond a common start for both models: the baseline perimeter, plus
PyreCast's start footprint (`analyze.RunStart`), plus FARSITE's ignition.
