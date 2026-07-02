# TODO / Notes

## ENTSO-E bidding zone list

`data/entsoe_cache/entsoe_bidding_zones.csv` is currently hand-crafted (60 zones).
The `entsoe` Python library already enumerates all bidding zones via its `entsoe.Area`
enum. Replace the CSV with a file generated from that enum, or validate directly against
it in `retrieve_entsoe.py`, so the list stays up to date automatically without manual edits.

## Cutout cache (like the grid data cache)

**Level 1 — keyed exact-match cache: DONE** (`common/_cutout_cache.py`). Each download is
keyed on its *actual* request parameters (module, rounded bbox, dx/dy, time range) and stored
under `cutouts/cache/<area>_<start>_<end>_<key>.nc` (+ a `.json` params sidecar). `02_make_cutouts`
checks the cache before hitting CDS and hardlinks the entry into the rule output on a hit
(copy fallback). This closes the stale-bounds hole of the old backup hack: a `mainland_bbox` /
`offshore_max_distance_km` edit changes the key, so it re-downloads instead of silently reusing
a differently-bounded cutout. The `_backup.nc` sibling is still honoured as a fallback (and gets
promoted into the cache on use), so previously pinned backups keep working — it can be retired
once nothing relies on it.

**Level 2 — coverage-aware reuse + partial fills: DEFERRED.** Slice a sub-request out of a larger
cached cutout when it spatially/temporally covers the request, and download only missing months
before concatenating. Worth it only if we start requesting overlapping sub-regions or extending
time ranges — today we request whole country-years, which Level 1 already serves. The added
slicing / grid-alignment / concat logic is where the subtle bugs live, hence deferred.

## Cutout bounds = land ∪ offshore (RESOLVED)

`download_cutout` now reads the pre-built `{cf_area}_geo.parquet` +
`{cf_area}_offshore_geo.parquet` (outputs of `01`/`01b`), unions them, and takes the
bounding box padded by `cutout.bbox_pad_deg` — matching the reference `02_make_cutouts`.
This fixes the earlier land-only + 1° pad, which clipped offshore-wind cells for wide-EEZ
areas (the offshore reach is `offshore_max_distance_km: 200` km ≈ 1.8°; e.g. AUS's offshore
zone extends ~2° past the land bbox). The geometry is no longer re-derived from Natural Earth
inside `02`, so the `geo → cutout → timeseries` chain is a single clean dependency line.

Re-trigger insulation now rests entirely on the `_backup.nc` hack (below), not on a
`download_cutout` → NE decoupling: if geometry regenerates, `download_cutout` re-runs but
copies the backup instead of hitting CDS. **Caveat:** that backup was built with whatever
bounds were current when it was cached, so a `mainland_bbox` / `offshore_max_distance_km`
edit won't actually re-bound a cutout that has a stale backup until the backup is refreshed —
a known property of the stopgap, to be retired with the proper cutout cache below.

## Wire in best-site CF (`07`) — drop the quarterly-cutout requirement

`07_make_bestsite_cf_timeseries` (and `_helpers.cutout_path()` / `QUARTERS`) still build the
CF year by concatenating four quarterly cutouts `{cc}_{year}_{q}.nc`. Stage 2 now emits a
single **annual** cutout `{cf_area}_{start_date}_{end_date}.nc`, so `07` cannot run against the
live pipeline (this is the WIP gap flagged in its `build_cf_year` docstring).

TODO:
1. Make `build_cf_year` read the single annual cutout instead of looping `QUARTERS`
   (drop `cutout_path()` / `QUARTERS` from `_helpers`).
2. Add a Snakemake rule wiring `07` the way `build_res_cf_profile` is wired
   (inputs: cutout + regions + offshore_regions; params from `config.res_cf`; output a
   bestsite parquet under `resources/res_cf/`).
3. Keep `07` driving its wind conversion off `config.res_cf.wind_cf` (`smooth` +
   `add_cutout_windspeed`), so best-site and national (`03`) CFs use identical settings.
   The reference pipeline was inconsistent here — national `03` was unsmoothed while
   best-site `07` used `smooth=True` — which quietly broke the "best-site ≥ national,
   same CF fields" claim in the res_cf README.

This also retires the last reason the reference `04_concat_quaters` / `05_combine_techs` logic
existed.

The diagnostic `06_resource_spread` shares the same quarterly-cutout migration and is also
unwired. If it gets wired, it needs one extra fix: `national_mean_from_csv()` reads
`{cc}_cf_{year}.parquet` (the combined per-country file that `05_combine_techs` used to write),
which the current pipeline no longer produces — so `nat_mean` silently falls back to the
spatial mean and the `uplift_*_vs_national` columns become meaningless. Repoint it at the
per-tech `{cf_area}_{tech}_country-average_{start}_{end}.parquet` files (from `03`), or drop
the national-mean column.

### Decision: keep best-site wind single-cell (do NOT restore the 3×3 block average)

The reference `07` averaged a 3×3 cell block for wind to "approximate wind farm-scale
variability and avoid artificial saturation at CF ≈ 1." Dropping it (single-cell sampling) was
correct — do not reintroduce it:

- **Cell-scale averaging is far too coarse.** ERA5 cells are tens of km across, so a 3×3 block
  spans ~100 km — vastly larger than a wind farm. It smears together genuinely different
  weather, not farm-internal spread.
- **Farm-scale spread is already covered by the `smooth` flag.** atlite's `smooth=True` is
  exactly the wind-farm-aggregation correction (it spreads the single-turbine power curve to
  represent a distribution of turbines), so the thing the 3×3 was meant to approximate is
  already handled at the right scale.
- **Saturation at CF ≈ 1 is physically real, not an artifact.** Above rated wind speed a
  turbine deliberately holds output flat at rated power by pitching its blades (pitch
  control/regulation) to shed the excess — so flat CF≈1 stretches in high wind are the turbine
  behaving correctly, and should be preserved, not smoothed away.

## Wire in the complementarity screen (`08`) + fix its latent traps

`08_complementarity_screen` chooses a RES-mix by *firmness* — it screens (onshore, offshore,
solar) cell triplets within `max_radius_km` and ranks them by
`score = w_coincidence·coincidence − w_correlation·mean_corr` (reward combined uptime, penalise
correlated profiles). It's the only thing in the pipeline that picks a mix by complementarity
rather than peak resource, and it's the closest stranded script to wireable (it already reads
`wildcards.country` / `output.top` / `output.avg`). It feeds the best-site RES-mix scenarios in
`h2_dri`. To wire it:

1. Shares `07`'s quarterly→annual migration (it imports `build_cf_year` etc. from `07`), plus a
   rule in `res_cf.smk` and moving the tuning constants out of `load_cfg`'s hardcoded block back
   into `config.res_cf`.
2. **Verify/fix the greedy bootstrap.** In `greedy_screen._find_best_triplet`, the step that
   picks the offshore complement (before solar is chosen) passes the onshore series into the
   solar slot: `score_triplet(ts_on[:, best_i], ts_off[:, j], ts_on[:, best_i], ...)`. So it
   scores `onshore + offshore + onshore` — double-counting onshore, ignoring solar. Either an
   intentional anchor-as-stand-in bootstrap or a copy-paste bug; only affects the greedy path
   (large candidate spaces — brute force is unaffected). Confirm before trusting greedy output.
3. **`quality_floor` is loaded but never applied** — the filter line is commented out
   (`qualified = valid #& (...)`) so it screens all valid cells, and the docstring's
   "pre-filter to quality_floor percentile" overstates what runs. Either re-enable it (shrinks
   the candidate space, speeds the screen) or drop the param + docstring line.
4. `save_average_profiles` reads the same nonexistent `{cc}_cf_{year}.parquet` as `06` (see the
   `07` section) — but here it `FileNotFoundError`s rather than falling back. Repoint it when
   fixing the `06` case.

## Surface resource-spread diagnostics in viz (in spirit)

The intra-country resource heterogeneity that `06_resource_spread` captures (best-site uplift,
spatial P90/P95/max vs national mean) is a genuinely useful story for reports — "how much does
siting matter in this country" — and currently lives only as a standalone CSV/parquet nobody
looks at. At some point, surface it in the `viz` pipeline: e.g. an uplift bar/whisker per
country × tech, or annotate the siting map with the P95-vs-national gap. Doesn't need to reuse
`06` as-is — the point is to get the *concept* into the report, computed from whatever CF grids
the live pipeline produces.

## CDS download monitoring

Atlite ERA5 cutout downloads go through the CDS API. With `monthly_requests=True`
(set in `download_cutout.py`), atlite submits one CDS job per month — up to 12 jobs
for a full year. Each job queues independently; global queue depths of 5000+ are
common and waits of several hours are normal.

### Checking status from the terminal

```python
import cdsapi
c = cdsapi.Client(quiet=True)
for j in c.client.get_jobs()._json_dict['jobs']:
    qos  = j['metadata']['qos']['status']
    user = qos.get('user', [{}])[0]
    print(
        j['status'], j['processID'],
        '| user queue:', user.get('queued', 0),
        'running:', user.get('running', 0),
    )
```

Key fields:
- `status` — `accepted` (queued), `running`, `successful`, `failed`
- `user.queued` / `user.running` — your own per-user position (limit: 1 running at a time)
- Global queue (`qos.limit`) gives context on overall wait times

### Automatic status updates: DONE

`common/_cds_monitor.py` (`cds_progress_logger`) wraps `cutout.prepare()` in `02_make_cutouts`
with a daemon thread that polls the CDS jobs endpoint every `res_cf.cutout.cds_poll_interval_s`
seconds (default 30) and logs one status line — running / queued / successful / failed plus your
per-user queue position — so long cutout runs are observable instead of appearing hung. All CDS
access is best-effort (poll errors are swallowed at DEBUG).

## Plotly PNG export needs Chrome (kaleido v1)

Every viz script writes a static PNG via `fig.write_image(...)`. With the current
pin (`python-kaleido` v1), kaleido no longer bundles a renderer and requires a
**Chrome/Chromium binary** to be present, or `write_image` raises
`ChromeNotFoundError`. The `.html` outputs (`fig.write_html`) need nothing and
always work.

One-time fix per environment:

```bash
# inside the env
plotly_get_chrome            # or, from Python:
python -c "import kaleido; kaleido.get_chrome_sync()"
```

This drops a Chrome-for-Testing build into a local cache (e.g.
`~/Library/Application Support/choreographer/`); it is per-machine and not shared,
so each dev box / HPC node / CI runner needs it (or must run HTML-only).

Options to make this less of a footgun: pin `python-kaleido<1` (the old bundled
renderer) in `environment.yaml`; or add a Snakemake `onstart`/setup check that
fetches Chrome if missing; or make PNG export optional (HTML always, PNG only when
a renderer is available) in the viz scripts.
