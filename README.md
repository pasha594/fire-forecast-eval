# fire-forecast-eval

Survey of pyrecast/ELMFIRE fire-spread forecast metrics vs actual perimeters:
area-based precision & recall per (fire, forecast run, weather percentile,
forecast window), plus how skill improves as forecasts refresh.

Pyrecast deletes forecast runs after ~1 day (verified: both
`data.pyrecast.org/fire_spread_forecast/` and their geoserver retain only the
latest ~2 runs per fire), so a GitHub Actions cron archives them every 2h into
Backblaze B2. Ground truth is Cornea's fire API (`fire-api-prod.web.app`),
which retains dated perimeter snapshots per fire (often several edits a day, but
typically about one newly mapped perimeter per fire-day).

## Pieces

| file | what |
|---|---|
| `collect.py` | scraper/archiver + fire matcher. CI runs it against the bucket; `--local` archives to `raw/` |
| `.github/workflows/collect.yml` | cron `17 * * * *` + manual dispatch; concurrency group = single writer |
| `.github/workflows/analyze.yml` | cron 2x daily: syncs ToA+perimeters on the runner, recomputes and publishes `metrics.csv` + per-fire `skill/*.json` to the bucket — no laptop needed; `report.html` falls back to the bucket copy |
| `analyze.py` | local: syncs the archive, computes growth-only P/R → `data/metrics.csv` |
| `report.html` | static Plotly report over `data/metrics.csv` |
| `map.html` | interactive MapLibre map, flat or over 3D terrain, over Cornea perimeter history + hotspots. Forecasts are decoded in the browser from the archive (all variables/percentiles); runs still on pyrecast are read live through `relay/` (its WMS tiles for the hourly variables, its raw time-of-arrival coverage until the archive has it). Incident maps (IAP / ops GeoPDFs) are Incibrief's: the page reads the fire's manifest from Incibrief's public bucket by Cornea id (`catalogs/incidents/id/<cornea_id>.json`) and overlays its tiles; nothing is processed here. Fires Incibrief wrote no manifest for show none |
| `relay/` | Cloudflare Worker (free plan) that passes pyrecast's geoserver through with CORS headers, which its geoserver refuses to send for any site but pyrecast.org; GetMap and GetCoverage only. It also fetches a clicked hotspot's full NASA FIRMS record (`/firms/area`, `/firms/availability`) with the `FIRMS_KEY` secret, so the key stays out of the public page. Deploy with `cd relay && npx wrangler deploy` |
| `overrides.json` | manual slug→cornea_id match overrides (string forces, null skips) |

Archive layout (bucket or local `raw/`):

```
forecast_archive/{slug}/{run_ts}/{pct}.tif              # ELMFIRE time-of-arrival rasters
forecast_archive/{slug}/{run_ts}/{pct}_{var}.tar        # hourly-mosaic granules, one tar per
                                                        #   (run, pct, var): crown-fire,
                                                        #   flame-length, hours-since-burned,
                                                        #   spread-rate + {pct}_isochrones.tar
perimeter_archive/{slug}/index.json, {epochms}.geojson  # Cornea perimeter snapshots
hotspot_archive/{slug}/{YYYY-MM-DD}.geojson             # VIIRS/MODIS detections by acq day
manifest.json                                           # collector state (source of truth)
fire_matches.json                                       # slug -> cornea fire
metrics.csv                                             # analyze.yml output (bucket only)
skill/index.json, skill/{slug}.json                     # per-fire skill for map.html (bucket only)
```

Since ~Aug 15 the static host 403s ToA tifs during a ~5-6h staging window
(occasionally forever); the collector falls back to a WCS `GetCoverage` from
their geoserver (verified bit-identical, recompressed before archiving,
`via: wcs` in the manifest). The fallback is ToA-only by design — GetCoverage
is server-side work for them, and the hourly granules would need thousands of
calls per run.

Tarring the ~169 hourly granules per variable into one object is deliberate:
it keeps bucket write operations ~100x below free-tier caps. Note two
upstream quirks: percentile tifs appear over ~2h after a run is listed, and
pyrecast **purges the hourly granule dirs within ~2–4h** — hence the 2h cron,
whose interval bounds worst-case granule loss; a variable tar holds whatever
survived at capture time (`got`/`n` in the manifest records the gap).

### Finished fires are carried over

A fire that has been finished for 18 days (`frozen_fingerprints`: no new runs, and its
perimeter and hotspot record unchanged) is not rescored every run. `metrics_state.json`,
published beside `metrics.csv`, records the scoring version (`SCORING_VERSION` in
`analyze.py`) and each finished fire's fingerprint and Cornea match when it was scored. The
next sync carries over every finished fire whose fingerprint and match are unchanged: its
published rows are copied into the new `metrics.csv`, and its archive files leave `raw/`, so
the Actions cache holds only active fires. Bump `SCORING_VERSION` whenever scoring changes:
finished fires are then rescored once, and their files are downloaded again for that run.
A one-fire run (`--only`) always rescores its fire.

### Per-fire skill files

`analyze.py --publish` also turns `metrics.csv` into small per-fire JSON files
for `map.html` (built in git-ignored `data/skill/`, uploaded to `skill/` in the
bucket, e.g. `https://f005.backblazeb2.com/file/fire-forecast-archive/skill/index.json`):

```
skill/index.json   {"v": 1, "generated": "<UTC ISO>", "slugs": {"<slug>": <scored rows>, ...}}
skill/{slug}.json  {"v": 4, "slug", "generated", "method": "exact-lead", "score_days": 7,
                    "starts": {"<run_ts>": "pyrecast" | "hotspots", ...},
                    "runs": {"<run_ts>": {"<pct>": [[lead_h, precision, recall, iou,
                      baseline_perim_ms, actual_perim_ms, 0, actual_mapped], ...]}}}
skill/start/{slug}/{run_ts}.tif   the start footprint a run was scored with, for every run
                    scored with one (uint8, 1 = burned when the run began, on the forecast's grid)
```

v4 (Nov 2026) holds every metrics row (no filters), one per scored perimeter, sorted by
`lead_h`: the perimeter's hours after the run (2 decimals). Scores are rounded to 3
decimals (`null` when undefined, e.g. recall when the fire didn't grow), and
`index.json` lists only fires with ≥ 1 scored row. The last four fields let the map
redraw a row's calculation as hits and misses: the two perimeter archive ids
(`baseline_perim_ms` is `null` when no perimeter was mapped before the run), 0 (each
perimeter is scored at its own time; v3 files had the old offset from a fixed horizon
here), and `actual_mapped` 1 when that time is when it was mapped, 0 when it's its last
edit (the map shows those in italics). v3 files (fixed horizons, filtered) still load. v3 adds `starts`: the runs whose burned-before
area takes in a start footprint beyond the baseline perimeter (see Metric
definitions), and where it came from; the map and report add
`skill/start/{slug}/{run_ts}.tif` to the baseline when they redraw a row. Footprints
go up for every run scored with one, not just the runs in `starts`, since the report
redraws rows outside the default filters too, and they go up before `metrics.csv`,
so no published row names a footprint that isn't there yet. If one still won't load,
the report draws the row without it and says so. A `--skill-only` rebuild uploads no
footprints: the last full run's stay in the bucket. On localhost the map reads `data/skill/` first,
so `analyze.py` output can be checked before it's published. Fires that
drop out keep a stale `skill/{slug}.json` in the bucket, but are no longer
listed. To rebuild just these files from the already-published `metrics.csv`
(no sync or recompute), run the analyze workflow manually with **skill_only**
checked (`analyze.py --skill-only`, needs the `S3_*` env).

## One-time setup

1. **GitHub**: push this repo to GitHub (`gh repo create ... --push`).
2. **Backblaze B2** (any S3-compatible store works; B2 chosen for its **hard spending
   caps** — set them under Account → Caps & Alerts):
   - Create bucket `fire-forecast-archive`, files **Public** (read-only public URL;
     the data is public anyway). Note the S3 endpoint (e.g.
     `s3.us-east-005.backblazeb2.com`) and friendly URL
     (`https://f005.backblazeb2.com/file/<bucket>`).
   - App Keys → create a key scoped to the bucket, read+write.
   - A bucket CORS rule allowing cross-origin GET: `map.html` draws every
     forecast layer from the archive, straight from the bucket.
3. **GitHub repo secrets** (repo → Settings → Secrets and variables → Actions), set these yourself:
   - `S3_ENDPOINT` (full https URL), `S3_ACCESS_KEY_ID`, `S3_SECRET_ACCESS_KEY`
   - Optional variable `S3_BUCKET` (defaults to `fire-forecast-archive`).
   - Legacy `R2_*` secrets still work as a fallback (collect.py checks `S3_*` first).
4. Trigger the workflow once manually (Actions → collect → Run workflow) and check the summary line + objects in the bucket.
5. **Pyrecast relay** for `map.html`'s live runs: `cd relay && npx wrangler deploy` (wrangler logs in
   to your Cloudflare account in a browser the first time), then set `RELAY` in `map.html` to the
   printed `*.workers.dev` URL. The free plan allows 100k requests a day; the map makes one per tile
   (only for runs still on pyrecast) and one per live time-of-arrival raster.
   The hotspot card's FIRMS lookups need the key as a Worker secret: `cd relay && npx wrangler secret put FIRMS_KEY`
   (a free MAP_KEY from firms.modaps.eosdis.nasa.gov/api/map_key: 5,000 transactions per 10 minutes; a click
   costs one, plus one per page load for the products' date ranges). For `npx wrangler dev`, put
   `FIRMS_KEY=...` in `relay/.dev.vars` (gitignored).

## Local usage

```
python3 -m venv .venv && ./.venv/bin/pip install requests "rasterio>=1.4,<1.5" "shapely>=2,<2.1"
./.venv/bin/python collect.py --local        # capture now, into raw/ (no cloud needed)
./.venv/bin/python collect.py --push         # upload a raw/ archive into the bucket  ⚠ see below
./.venv/bin/python analyze.py --sync         # mirror the bucket + write data/metrics.csv
./.venv/bin/python analyze.py --inspect raw/forecast_archive/ca-bug/<run>/50.tif
python3 serve.py                             # then http://localhost:8090/report.html
                                             #  and http://localhost:8090/map.html
```

⚠ `--push` and the CI workflow must not run at the same time — the manifest
assumes a single writer. Disable the workflow (or wait for it to finish) first.

## Metric definitions

For a forecast issued at `run_ts`, percentile `p`, scored against perimeter A at its
own time, `lead_h` hours after the run (a worked example is on the site at `/methodology`):

- **baseline B** = latest Cornea perimeter mapped (or, lacking a reliable mapped time, last edited) by the end of `run_ts`'s minute, as published by then, rasterized onto the forecast's own UTM grid — never one from after the run. (`run_ts` is cut to the minute.)
- **start S** = where the forecast itself started burning. Pyrecast seeds its runs with satellite detections on top of the perimeter (`run_ts` is usually a VIIRS overpass to the minute), so its start reaches well beyond B — a median ~2,400 ac, even when B was mapped within a minute of the run. Its time of arrival is never set there, so growth into S can't be predicted; counted as growth, it was all misses (median recall was ~0.19 instead of ~0.68 before S was subtracted, Nov 2026). S is the first granule of the run's `hours-since-burned` output (values 400 burned / 401 burning at the start; the same at every percentile), which `analyze.py --sync` fetches with two ranged reads of an archived `{pct}_hours-since-burned.tar` into `raw/forecast_archive/<slug>/<run_ts>/start.tif`. When no archived tar caught that granule (~8% of scored runs), every satellite detection footprint up to `run_ts` stands in for it. Scored with the stand-in instead of the real start, the 167 runs on four test fires that have both get recall within a median 0.10 of the real-start recall, and precision and IoU within 0.004 (Nov 2026, with the footprint sizes under Hotspot capture; the old, larger MODIS and Landsat footprints gave 0.11 and 0.007). `start_basis` records which (`pyrecast`, `hotspots`, or blank for neither); `start_acres` is S's area and `burned_before_acres` B ∪ S's.
- **actual A**: each run is scored against every perimeter it picks (`analyze.scored_perimeters`, recorded in `chosen`), over its first 7 days (`day` = 24 h periods from the run's start): every perimeter with a mapped time (`mapped`); on a day with none, the latest perimeter published that day, timed by its last edit and not a re-publish of an earlier shape (`daily`); and the perimeter the previous fixed-horizon method showed at 12h–14d (nearest `run_ts + H` within ±12h, shown when within ±3h), rescored at its own time so no earlier analysis is lost (`previous`)
- **pred_new** = pixels with time-of-arrival ≤ `lead_h`, minus B ∪ S (every model gets the same B ∪ S: one common start, `analyze.score_grid`)
- **precision** = area(pred_new ∩ act_new) / area(pred_new), **recall** = area(pred_new ∩ act_new) / area(act_new), **IoU** = area(pred_new ∩ act_new) / area(pred_new ∪ act_new), where act_new = A − (B ∪ S). Every row is kept (no minimum growth, timing or off-grid filter; `act_offgrid_frac` is recorded). No growth: precision and IoU 0 if the forecast predicted growth, recall blank; all blank if neither grew

Treating all of S as burned is safe even though S overstates the burned area
(only ~42% of S beyond B ends up inside the fire's final perimeter, as 375m
pixels overhang the fire's edge): the part that never burns is neither predicted
nor actual growth, so it can't move a score. Padding the start doesn't pay
either — the forecast's front then starts too far out, which costs precision.

**Perimeter times.** Every Cornea snapshot is one *edit* of a WFIGS polygon, and
its date (and archive file name) is the edit time, `poly_DateCurrent` — a median
~17h after the polygon was mapped (as of Oct 2026; p90 ~5 days, since old
captures get re-edited). Each snapshot is scored at the most reliable time it
has:

1. **mapped** — `poly_PolygonDateTime`, if it is at most 5 minutes after its own
   edit, no earlier than a day before discovery (61 days before the edit when
   discovery is unknown; typos like year 0926 exist), and the polygon's area is
   within 2% and 25 acres of the first edit that carried that time (editors often
   redraw a growing perimeter without updating it, and leave `poly_GISAcres`
   stale too, so the area is measured from the geometry). Matching edits are
   corrections of that capture: an actual uses the latest one, a baseline the
   latest one published by the run (else the first), so no growth published
   after the run can enter it.
2. **edited** — otherwise `poly_DateCurrent`, the publish time, which is no
   earlier than the mapping it stands in for. (`poly_CreateDate` was checked
   and adds nothing: it is either the upload time or the polygon record's
   original creation, often days before the perimeter.)

No snapshot is dropped. `baseline_time_basis` / `actual_time_basis` say which
time each row used, and report.html can restrict to mapped times only.
`*_perim_ms` are snapshot ids (edit times); `*_time_utc` are the times scored.

**Hotspot capture** (`hs_*` columns) is the deliberately *separate* satellite
companion metric: detections within exactly `(run, run+H]` (their timestamps
are precise, so no tolerance is needed), each rasterized as a disk around it,
B ∪ S subtracted; `hs_capture` = predicted ∩ hotspot area / hotspot area. A
detection is a fire somewhere in its pixel, often far smaller than the pixel,
so each disk is only as big as it needs to be. VIIRS (375 m px) gets r = 187.5 m,
the least at which neighbouring detections, 375 m apart at nadir, still touch;
off-nadir pixels are bigger, and those don't touch. MODIS (1 km px) gets
r = 250 m, not the full 500 m, which claimed ~194 acres burned per detection;
its detections are ~1 km apart and don't touch either way. Landsat (30 m px)
gets r = 15 m; it used to borrow VIIRS's 187.5 m. Every detection also takes the
grid cell it falls in. The same disks make the `hotspots` stand-in for S and
decide whether a far-off detection ever touches a later perimeter (Nov 2026;
before then MODIS was 500 m and Landsat 187.5 m). It is
recall-flavored only — hotspots undersample burned area (clouds, smoke,
overpass gaps), so a hotspot-based "precision" would punish correct
predictions; perimeters remain the precision source. The two truths are kept
separate rather than merged: merging 375m–1km footprints into perimeter truth
would inflate areas and destroy error attribution.

Detections are screened first, once per fire against its whole record, and the
ones that don't look like the fire are left out (neither captured nor missed):
- **artifact** (`hs_artifact_n`): more than 1 km from the perimeter current
  when it was detected, and its footprint never touches a perimeter mapped
  after it. With no later perimeter yet it can't be judged and is kept, so a
  live fire's newest rows can change as perimeters arrive.
- **transient** (`hs_transient_n`): no other detection within 1 km on any
  satellite's nearest pass before or after its own, its own satellite
  included. A pass is one source's detections less than 20 minutes apart; each
  source's sequence is over passes that detected this fire (non-artifacts
  only), so cloudy or empty overpasses don't break it. Every satellite's
  nearest pass counts, not just the one pass right before and after in time:
  that one often saw the fire somewhere else, and checking only it left out
  23% of Landsat and 18% of MODIS detections on four test fires (Nov 2026),
  nearly all of which another satellite's pass confirms. Now 1–2%.

`hs_n` counts only the detections scored; the two `_n` columns count those
left out in the same window. The start footprint's fallback (detections up to
the run) takes every detection, screened or not, since it stands in for what
the run was seeded from.

Growth-only by design: comparing full footprints would be dominated by
already-burned area. All raw areas are in `metrics.csv` so other metrics can
be derived. Caveats worth knowing: recall is measured against *on-grid* actual
growth (`act_offgrid_frac` quantifies clipping); and predicted area is *not* strictly
monotone in percentile (observed inversions in ~1/3 of groups, median 11%
relative) — the percentile runs are independent weather draws, not nested
contours, so this is a property of the source data.

## Costs (full ingestion)

~10–14GB/day when all variables are captured (fire-count dependent). On B2:
storage **$6/TB-month ≈ $2.50/month per month of accumulated archive** (a full
season ~2.5TB ≈ $15/month); egress free up to 3× stored volume per month
(incremental `--sync` uses ~1×); transaction ops negligible thanks to the
tarred layout. **Hard caps in Account → Caps & Alerts make overruns
impossible — B2 halts service instead of billing.** Retention is the other
lever: prune old `*_​{var}.tar` objects (keep the ToA tifs — the skill metrics
need them) to cut storage ~10x. GitHub Actions: free (public repo).

## Teardown after the survey

- Disable the workflow (Actions → collect → ⋯ → Disable) or delete the cron block.
- The relay can go too (`cd relay && npx wrangler delete`); the map then shows only archived runs.
- Keep `data/metrics.csv` + `report.html` committed.
- Bucket can sit free under 10GB, or archive a final `tar` of it and delete —
  **metrics are unrecomputable without the tifs**.
