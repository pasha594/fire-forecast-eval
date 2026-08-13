# fire-forecast-eval

Survey of pyrecast/ELMFIRE fire-spread forecast skill vs actual perimeters:
area-based precision & recall per (fire, forecast run, weather percentile,
forecast window), plus how skill improves as forecasts refresh.

Pyrecast deletes forecast runs after ~1 day (verified: both
`data.pyrecast.org/fire_spread_forecast/` and their geoserver retain only the
latest ~2 runs per fire), so a GitHub Actions cron archives them every 6h into
Cloudflare R2. Ground truth is Cornea's fire API (`fire-api-prod.web.app`),
which retains dated perimeter snapshots per fire (~2–3/day for active fires).

## Pieces

| file | what |
|---|---|
| `collect.py` | scraper/archiver + fire matcher. CI runs it against R2; `--local` archives to `raw/` |
| `.github/workflows/collect.yml` | cron `17 */6 * * *` + manual dispatch; concurrency group = single writer |
| `analyze.py` | local: syncs the archive, computes growth-only P/R → `data/metrics.csv` |
| `report.html` | static Plotly report over `data/metrics.csv` |
| `map.html` | interactive map: pyrecast forecast layers (live WMS, all variables/percentiles) over Cornea perimeter history + hotspots |
| `overrides.json` | manual slug→cornea_id match overrides (string forces, null skips) |

Archive layout (R2 bucket or local `raw/`):

```
forecast_archive/{slug}/{run_ts}/{pct}.tif              # ELMFIRE time-of-arrival rasters
forecast_archive/{slug}/{run_ts}/{pct}_{var}.tar        # hourly-mosaic granules, one tar per
                                                        #   (run, pct, var): crown-fire,
                                                        #   flame-length, hours-since-burned,
                                                        #   spread-rate + {pct}_isochrones.tar
perimeter_archive/{slug}/index.json, {epochms}.geojson  # Cornea perimeter snapshots
manifest.json                                           # collector state (source of truth)
fire_matches.json                                       # slug -> cornea fire
```

Tarring the ~169 hourly granules per variable into one object is deliberate:
it keeps R2 Class A (write) operations ~100x below the free-tier cap. Note two
upstream quirks: percentile tifs appear over ~2h after a run is listed, and
pyrecast **purges the hourly granule dirs within hours** — hence the 3h cron;
a variable tar holds whatever survived at capture time (`got`/`n` in the
manifest records the gap).

## One-time setup

1. **GitHub**: push this repo to GitHub (`gh repo create ... --push`).
2. **Backblaze B2** (any S3-compatible store works; B2 chosen for its **hard spending
   caps** — set them under Account → Caps & Alerts):
   - Create bucket `fire-forecast-archive`, files **Public** (read-only public URL;
     the data is public anyway). Note the S3 endpoint (e.g.
     `s3.us-east-005.backblazeb2.com`) and friendly URL
     (`https://f005.backblazeb2.com/file/<bucket>`).
   - App Keys → create a key scoped to the bucket, read+write.
   - Optional: a bucket CORS rule allowing cross-origin GET lets `map.html`
     render archived rasters straight from the bucket on non-localhost origins.
3. **GitHub repo secrets** (repo → Settings → Secrets and variables → Actions), set these yourself:
   - `S3_ENDPOINT` (full https URL), `S3_ACCESS_KEY_ID`, `S3_SECRET_ACCESS_KEY`
   - Optional variable `S3_BUCKET` (defaults to `fire-forecast-archive`).
   - Legacy `R2_*` secrets still work as a fallback (collect.py checks `S3_*` first).
4. Trigger the workflow once manually (Actions → collect → Run workflow) and check the summary line + objects in the bucket.

## Local usage

```
python3 -m venv .venv && ./.venv/bin/pip install requests "rasterio>=1.4,<1.5" "shapely>=2,<2.1"
./.venv/bin/python collect.py --local        # capture now, into raw/ (no cloud needed)
./.venv/bin/python collect.py --push         # upload a raw/ archive into R2  ⚠ see below
./.venv/bin/python analyze.py --sync         # mirror the bucket + write data/metrics.csv
./.venv/bin/python analyze.py --inspect raw/forecast_archive/ca-bug/<run>/50.tif
python3 serve.py                             # then http://localhost:8090/report.html
                                             #  and http://localhost:8090/map.html
```

⚠ `--push` and the CI workflow must not run at the same time — the manifest
assumes a single writer. Disable the workflow (or wait for it to finish) first.

## Metric definitions

For a forecast issued at `run_ts`, percentile `p`, window `H` hours:

- **pred_new** = pixels with time-of-arrival ≤ H, minus baseline
- **baseline B** = latest Cornea perimeter at/before `run_ts`(+3h), rasterized onto the forecast's own UTM grid
- **actual A** = Cornea perimeter nearest `run_ts + H` (required within ±12h; signed offset recorded)
- **precision** = area(pred_new ∩ act_new) / area(pred_new), **recall** = area(pred_new ∩ act_new) / area(act_new), where act_new = A − B

Growth-only by design: comparing full footprints would be dominated by
already-burned area. All raw areas are in `metrics.csv` so other metrics can
be derived. Caveats worth knowing: recall is measured against *on-grid* actual
growth (`act_offgrid_frac` quantifies clipping); ELMFIRE's internal ignition
perimeter differs slightly from Cornea's baseline perimeter, which penalizes
precision near the perimeter edge; and predicted area is *not* strictly
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
- Keep `data/metrics.csv` + `report.html` committed.
- Bucket can sit free under 10GB, or archive a final `tar` of it and delete —
  **metrics are unrecomputable without the tifs**.
