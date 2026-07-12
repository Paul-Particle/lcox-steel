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

**Update:** the multi-year need is better served by stitching the *derived CF*, not the weather
cutout — see "CF stitching over time (and space)" below. That keeps Level 1's exact-match cutout
cache untouched and moves temporal (and later spatial) assembly onto the smaller, safer CF-grid
layer. Coverage-aware *cutout* reuse (this Level 2) stays deferred and may not be needed at all.

## CF stitching over time (and space): tile → assemble → select refactor

**Status: DESIGN ONLY — needs sign-off before implementation.** This is invasive: it changes the
res_cf *artifact model* (introduces a per-cell CF grid intermediate) and refactors `03`/`07`/`08`
into a *compute* half and a *select* half. It is well beyond Snakemake rewiring — do not start it
without agreeing the shape below.

### Why this exists

The optimisation must run across multiple years (decades, in principle), but the CF pipeline today
produces one artifact per `(area, tech, variant, start, end)` from a single monolithic script per
variant, each reading one cutout. A multi-year run therefore either rebuilds/re-downloads
everything or cannot reuse the per-year cutouts already on disk.

Key realisation: **stitch time on the derived CF, not on the weather cutout.** Per-cell CF is
*timestep-local* — CF at hour *t* depends only on that cell's weather at *t* (`smooth` is
power-curve smoothing in wind-speed space, not temporal; solar geometry is per-timestamp) — so
computing CF per year and concatenating is bit-identical to computing it on a multi-year cutout.
Concatenating CF grids is also strictly simpler and safer than concatenating weather cutouts: one
value per cell/time/tech, no ERA5/ERA5T `expver` structure and no multi-feature merge — i.e. none
of the machinery where the historical download-corruption lived.

By contrast, every **selection** step is period-dependent and *cannot* be tiled: P95 ranking,
complementarity coincidence/correlation, top-N — all read the whole-period distribution. Tiling
those per year and concatenating gives per-year winners glued together (wrong). So selection must
move *above* the assembly layer and run once on the assembled multi-year CF.

### Layering

```
weather cutout(area, year)                 [cached; existing Level 1 cache]
  └─► CF-grid tile(area, tech, year)        [NEW compute unit: build_cf_year → per-cell grid]
        └─► assemble+slice → CF grid(area, tech, [S,E])   [NEW: concat along time, slice to request]
              ├─► country-average(area, tech, [S,E])        [03 select half: spatial mean]
              ├─► bestsite-p95 / anchored(...)              [07 select half]
              └─► complementarity(area, [S,E])              [08 select half]
```

The per-cell CF grid is effectively a *derived, lightweight cutout* (one CF value per
cell/time/tech). It is the canonical tile; `07` and `08` already share its producer — `08`
imports `build_cf_year` straight from `07`.

### Per-script compute/select split

| Script | Today (monolithic) | Compute → **tile** (per year) | **Assembly** | Select → **post-assembly** |
|---|---|---|---|---|
| `03` | cutout → region-weighted spatial mean → country-average ts | country-average ts / yr (via atlite `indicatormatrix`; timestep-local) | concat + slice | *none* — the mean is the output |
| `07` | cutout → per-cell CF grid → period-mean rank → P95 cell + anchor/counterpart matching → extract cells' ts | per-cell CF grid / yr (`build_cf_year`) | concat + slice grid | `find_p95_cell`, `find_nearest_valid_cell` / land-score matching, cell extraction |
| `08` | cutout → per-cell CF grid → all-candidate ts matrix → triplet screen | per-cell CF grid / yr (`build_cf_year`) | concat + slice grid | candidate masking, ts-matrix build, brute/greedy triplet screen, top-N + top-1 ts |

`07` and `08` consume the **same** per-cell CF-grid tile. `03` is the odd one out: its
country-average uses atlite's `indicatormatrix` (area-weighted by cell/region overlap), whereas
`07`/`08` use an unweighted in-region cell mask. Either keep a separate (small) country-average
tile to preserve `03`'s weighting, or derive country-average from the per-cell grid and accept
unweighted means — **decision needed** (a pre-existing inconsistency; don't paper over it).

### Snakemake wiring (fits repo conventions)

- **Universal assembly, no annual special-case.** pypsa consumes *assembled* outputs for every
  request. A single full year is a degenerate one-tile assembly (passthrough); a sub-year (month)
  request is a one-tile slice; multi-year is concat-of-N-tiles + edge slice. One rule, one code
  path — this is what dissolves the "extra rule for everything except annual" awkwardness.
- **Declarative fan-in, not imperative commissioning.** The assembly rule *declares* its
  covering-year tiles as `input:`; Snakemake builds whatever is missing (→ missing cutouts →
  cache/CDS). No script that pokes at disk and "commissions" downloads.
- The covering-year set is wildcard-dependent, which would normally need an input function
  (forbidden by CLAUDE.md). Avoidable: `projects.csv` is a **closed** request set, so expand each
  `(area, variant, start, end)` into its covering years once at Snakefile module scope (header
  Python is allowed), then fan tiles in with `collect(...)` / `lookup(...)`.
- **Granularity = calendar year.** Matches the per-year weather cutout (the actually-expensive
  artifact). Sub-year requests slice; month-tiles would add file sprawl for zero reuse benefit,
  since the cutout is per-year regardless.

### Open decisions for the coworker discussion

1. **Tile storage:** persist per-cell CF grids (netCDF like a cutout, or wide parquet) vs
   assemble-on-the-fly from cutouts. Large native-res countries (BRA/AUS) over decades make this
   non-trivial — needs sizing. Ties into #6.
2. **Country-average weighting:** unify `03` (`indicatormatrix`) with `07`/`08` (unweighted mask),
   or keep a separate country-average tile.
3. **Filename/namespace change:** `tiles/…` vs assembled `resources/res_cf/…`; confirm
   `h2_dri/solve_network.py` consumes the assembled files unchanged.
4. **Selection semantics over multi-year:** confirm P95 / complementarity should rank on the
   *full-period* distribution (assumed yes). Leap / short years are harmless — per-cell means are
   over all timesteps.
5. **Fold in the already-logged latent fixes** (the `07`/`08` quarterly→annual wiring, the greedy
   bootstrap `score_triplet` bug, the unused `quality_floor`, the stale `{cc}_cf_{year}.parquet`
   reads) as part of this refactor, or keep separate. **This refactor *is* the "broader
   cutout-cache refactor" those WIP notes said to wait for.**
6. **Drop `coarse`/resolution** (see the cutout-cache discussion): removes `dx/dy` from the cutout
   key and puts every area/year on the identical native ERA5 grid — which is what makes time-concat
   (and any future spatial slicing) clean and uniform. Best done here.

### Spatial (lower priority)

The same tile model extends to space: slice a sub-region out of a cached larger grid (pure `sel`) —
safe once resolution is uniform (#6). *Stitching* adjacent tiles into a bigger domain is only
needed if one analysis spans multiple areas' domains at once; defer until then.

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

- script07 Colocation heuristic varification  @PP
- write tests to make sure ERA5/CF/Atlite downloads dont fail silently for all 5 countries, year 2023 and hourly timesteps @PP
- grid pipeline compare to Energycharts.info API
- (low prio) azimuth sweep 270-0-90 degrees only for Southern hemisphere right now

- VIC test run finding: offshore-anchored scenario finds 0 candidates for both wind_onshore and solar within max_radius_km=100 — the offshore anchor's file has no co-located land options at all. Likely the offshore P95 cell sits >100km from any valid land cell for this region. Investigate: increase radius for offshore-anchor scenarios specifically, or revisit offshore P95 cell selection. Related to the broader anchor-selection review (BRA P95 issue).

- Dead spatial_matching_res_mix config, but still wired into Snakemake rules: Script 07's code no longer reads spatial_matching_res_mix (the old anchor co-location logic that used it — MAX_RADIUS_KM/QUALITY_FLOOR_FRAC — was deleted when we moved anchor co-location to 07b). However, res_cf.smk still passes spatial_matching_res_mix as a params: lookup to both build_bestsite_p95 and build_anchored_res_mix_cfs rules (lines ~124-125, ~148-149), and the config value itself still exists in config/config.yaml (~line 76). Since Snakemake params: values don't need to be consumed by the script to avoid erroring, this isn't currently broken — just stale plumbing. Needs someone to remove the spatial_matching_res_mix params from both rules in res_cf.smk, then remove the config section, in one coordinated edit (removing config first would break the rules' lookup). Low priority, not urgent.