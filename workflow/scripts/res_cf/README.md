# Renewable capacity factors (res_cf)

These scripts turn ERA5 reanalysis into hourly, per-unit (0–1) capacity-factor
series for three technologies — `wind-onshore`, `wind-offshore` and `solar` —
using [atlite](https://atlite.readthedocs.io). The rules are in
`workflow/rules/res_cf.smk`; the settings are the `res_cf:` block and the
per-area entries under `areas:` in `config/config.yaml`. The results land in
`resources/timeseries/`, where `solve_network` picks them up.

The series are purely climate-driven: apart from the land-sea cut described
below, no land-use, grid or permitting exclusions are applied.

## Stages

| Stage | Rule | Output |
|---|---|---|
| `a_make_area_geometry.py` | `make_area_geometry` | `resources/shapes/{area}_geo.parquet` |
| `b_make_offshore_geometry.py` | `make_offshore_geometry` | `resources/shapes/{area}_offshore_geo.parquet` |
| `c_retrieve_area_cutout.py` | `retrieve_area_cutout` | `cutouts/{area}_{start_date}_{end_date}.nc` |
| `d1_area_average.py` | `area_average` | `resources/timeseries/{area}_{tech}_area-average_{start_date}_{end_date}.parquet` |
| `d2_bestsite_p95.py` | `bestsite_p95` | `resources/timeseries/{area}_{tech}_bestsite-p95_{start_date}_{end_date}.parquet` |
| `d3_anchor_colo.py` | `anchor_colo` | `resources/timeseries/{area}_{tech}_anchor-colo-n{N}_{start_date}_{end_date}.parquet` |
| `d4_tilt_mix.py` | `tilt_mix` | `resources/timeseries/{area}_solar_tilt-mix-n{N}_{start_date}_{end_date}.parquet` |
| `d5_multi.py` | `multi` | `resources/timeseries/{area}_{tech}_multi-n{N}_{start_date}_{end_date}.parquet` |

`{area}` is a key of `areas:` (`DEU`, `VIC1`, …), `{tech}` is one of the three
technologies, and the dates are `YYYYMMDD`. Every d-stage reads the cutout and the
onshore geometry; all except `tilt_mix` also read the offshore geometry. Shared
code sits in `_helpers_res_cf.py` and `workflow/common/`.

## Geometry (a, b)

**a — onshore.** Dissolves the Natural Earth 1:110m Admin 0 country whose
`ADM0_A3`, `SOV_A3` or `ISO_A3` is `areas.<area>.iso3`. With
`areas.<area>.mainland_bbox` (`[lon_min, lon_max, lat_min, lat_max]`), only
polygon parts whose centroid lies inside the box are kept (this is how `FRA`
drops its overseas territories). The file holds one row, `(region, geometry)`;
`region` is `areas.<area>.region`, the tag every later stage selects on.

The rule builds top-level areas (those in no one's `zones:`). Zone geometries are
committed under `resources/shapes/`.

**b — offshore.** Takes the World EEZ v12 polygons with `ISO_TER1 == iso3` and
`POL_TYPE == "200NM"` (joint-regime areas are excluded), dissolves them, and keeps

    offshore = (EEZ − land) ∩ (land buffered by res_cf.offshore_max_distance_km)

with the buffer built in the equal-area EPSG:6933. An area with no sea within that
distance ends up with an empty geometry.

Where to download both source files is under
[External data files](../../../README.md#2-external-data-files).

## Cutout (c)

One ERA5 cutout per area and window. Its extent is the bounding box of the onshore
and offshore geometries together, padded by `res_cf.cutout.bbox_pad_deg`, so the
offshore zone is inside it. `areas.<area>.coarse: true` requests a 0.5° grid
instead of ERA5's native 0.25°. The time axis runs from `start_date` 00:00 to
`end_date` 23:00.

The script takes the first of these that exists:

1. an entry in the keyed cache `cutouts/cache/`, matched on the actual request
   (module, bounding box, resolution, time range) and hardlinked into place;
2. a pinned `cutouts/{area}_{start_date}_{end_date}_backup.nc`, copied into place
   and then added to the cache;
3. a fresh CDS download, which is then added to the cache.

Whichever source is used, the result must pass `common/_cutout_qc.py` (exact
hourly coverage, no gaps or duplicates, no NaN, no ERA5/ERA5T mix) before the
rule succeeds. See [Cutout caching](../../../README.md#cutout-caching) and
[Watching CDS / ERA5 cutout downloads](../../../README.md#watching-cds--era5-cutout-downloads)
in the root README.

| `res_cf.cutout.` key | Effect |
|---|---|
| `bbox_pad_deg` | padding around the geometry bounds, degrees |
| `monthly_requests` | split each CDS request by month (smaller requests, more queue time) |
| `cds_poll_interval_s` | how often the CDS queue status is polled while waiting |
| `cache_warn_size_gb` | warn when `cutouts/cache/` grows past this (0 disables) |
| `min_free_disk_gb` | warn before a download when free disk is below this (0 disables) |

With `LCOX_STORE` set, `cutouts/cache/` and atlite's scratch `.atlite-cache/`
live under that path, so worktrees share downloads (`common/_paths.py`).

## Capacity factors (d1–d5)

### Technology settings

| Tech | Config keys (`res_cf.`) |
|---|---|
| `wind-onshore` | `wind_onshore_turbine` |
| `wind-offshore` | `wind_offshore_turbine` |
| `solar` | `pv_panel`, `pv_orientation` |

Wind uses `wind_cf.smooth` (atlite's smoothed power curve, which approximates a
farm better than a single turbine) and `wind_cf.add_cutout_windspeed` (keeps the
cut-out wind speed when smoothing). `d4` uses `pv_panel` but sets its own
orientations. All series are clipped to [0, 1] and indexed by an hourly
`DatetimeIndex` named `time`.

### Land-sea eligibility and area weighting

Coarse ERA5 coastal cells mix sea and land, which inflates onshore wind CF. For
the land techs (`wind-onshore`, `solar`), a cell is dropped when its land fraction
is below `min_land_fraction` times the largest land fraction in the area (0
disables the cut). The land fraction comes from `eligibility_source`:

- `indicatormatrix`: overlap of the cell with the area polygon;
- `availabilitymatrix`: a finer coastline from `res_cf.availability.land_shapes`
  (Natural Earth 1:10m land), rasterised at `availability.res` metres in
  `availability.crs`.

Offshore wind is never cut. The surviving fractions are multiplied by cos(lat),
where lat is the cell's latitude, so that averages and percentiles are weighted by
physical area rather than by degree area. `d1`, `d2` and `d3` use these weights;
`d4` and `d5` select cells as described in their sections. Background is in
`docs/land_eligibility_design.md`.

### The P95 cell

`d2`, `d3` (for its anchor) and `d4` use the same selection: compute each cell's
annual-mean CF, take the area-weighted 95th percentile of those means over the
weighted cells, and pick the cell whose mean is closest to it. The output is that
one cell's hourly series, so it is the profile of an actual site.

### d1 — `area-average`

atlite's aggregation with `per_unit=True` over the weight matrix, which gives the
area-weighted mean hourly CF over the eligible cells. Offshore uses the offshore
geometry. One column, named by the tech.

### d2 — `bestsite-p95`

The P95 cell's hourly series. One column, named by the tech. A tech with no
eligible cell in the area is skipped.

### d3 — `anchor-colo-n{N}`

`{tech}` names the anchor. The anchor cell is the anchor tech's P95 cell. For
each of the other two techs, the candidate cells are those within
`anchor_colocation.max_radius_km` of the anchor whose centre is inside that tech's
geometry, whose mean CF is positive, and (for land techs) which pass the land-sea
cut. Each candidate is scored against the anchor:

    score = w_coincidence × coincidence − w_correlation × correlation

where `coincidence` is the share of hours in which the mean of the two CFs exceeds
`coincidence_threshold`, `correlation` is the Pearson correlation of the two
hourly series, and the weights and threshold come from `res_cf.anchor_colocation`.
The top N per tech are kept, with N taken from the variant. Distance is not part
of the score, since the solve prices the transmission link by distance. A tech
with no valid cell in range is left out.

Columns are `{tech}@anchor` and `{tech}@c00`, `{tech}@c01`, … (best score first).
The Arrow schema metadata carries `site_coords` (`{column: {lat, lon}}`),
`demand_site` (the anchor column) and `anchor_colocation` (score, coincidence,
correlation and distance per candidate, candidate counts, the settings used, and
a timestamp).

### d4 — `tilt-mix-n{N}`

Solar only. Finds the P95 cell of the `latitude_optimal` CF, weighting each cell
by its overlap with the geometry × cos(lat), then sweeps N
azimuths evenly spaced over the equator-facing direction ± 90° (south, 180°, in
the northern hemisphere; north, 0°, in the southern). With an odd N the
equator-facing azimuth is one of them. For each azimuth, the slope (0–90° in 1°
steps) that maximises annual plane-of-array irradiance is found with the
Hay-Davies model on that cell's irradiance and solar-position data, and the CF is
then computed with atlite at that slope and azimuth. Columns are
`solar_az{azimuth}`.

### d5 — `multi-n{N}`

Candidate cells for siting inside the solve. It takes every cell whose centre lies
in the geometry (the offshore one for `wind-offshore`) and keeps a regular
lattice of them, every k-th cell in both grid directions with
k = round(sqrt(cells inside / N)). If that gives more than N cells it keeps the N
with the highest mean CF, and if fewer it adds the highest-mean remaining cells.
An area with N or fewer cells uses all of them. Columns are `{tech}@c00`,
`{tech}@c01`, … in grid order, with `site_coords` in the schema metadata.

## How the solve uses the series

Each row of `config/scenarios.csv` names a `tech`, `variant`, `area` and window,
and `solve_network` reads
`resources/timeseries/{area}_{tech}_{variant}_{start_date}_{end_date}.parquet`
for every non-grid row of the run. Costs are looked up in
`config/assumptions.yaml` by the column's tech key; `solar_az*` columns use the
`solar` entry.

- **Single site** (`area-average`, `bestsite-p95`, `tilt-mix-n{N}`): the columns
  of all inputs are joined into one frame and each becomes one generator on a
  single electricity bus. Two inputs with the same column name are an error.
- **Multi-site** (`anchor-colo-n{N}`, `multi-n{N}`): triggered when any column
  contains `@`. Each column becomes a site on its own bus, joined by an HVDC link
  to the demand site, which is the `demand_site` column's cell for `anchor-colo`, or
  the area's representative point for `multi`. Every CF input of such a run must
  then use `@` columns, so the two groups cannot be mixed in one run.

## Tools outside the workflow

- `spot_check_cutout.py` checks that CDS serves a year cleanly: it
  downloads a 1° square around a fixed point per country (`de`, `fr`, `es`,
  `aus`, `bra`) for whole years into `cutouts/spotcheck/`, runs the cutout QC and
  exits non-zero on failure:
  `python workflow/scripts/res_cf/spot_check_cutout.py --areas de fr --years 2024`.
- `download_cutouts_batch.sh` builds cutout targets one at a time through
  Snakemake (CDS runs one job per user), skips files that already exist, carries
  on past failures, logs to `logs/overnight_downloads.log` and ends with a QC
  summary. Pass the targets as arguments, e.g.
  `bash workflow/scripts/res_cf/download_cutouts_batch.sh cutouts/DEU_20240101_20241231.nc`.
- Each Python script has hardcoded defaults at the top so it can run without
  Snakemake.
- `reference/hannah-pypsa-lcoe-lcoh/` is a copy of an earlier pipeline kept as
  reference material. It is not part of the workflow.
