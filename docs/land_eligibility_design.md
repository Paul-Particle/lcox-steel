# Land-eligibility weighting for onshore CF (issue #41)

Scoping doc. Branch: `land-eligibility-mask-issue41`.

## Problem (recap)

Coarse ERA5 (~27–31 km) coastal cells blend smooth-sea and rough-land wind, so
onshore wind CF is inflated wherever a cell contains sea. Confirmed across all six
cf_areas (FR/DE/ES/AUS/BRA/VIC) by the decomposition in
`notebooks/res_cf_coastal_cf_exploration.ipynb` and the multi-country re-run posted
to #41:

- Removing coastal border cells collapses the weighted-vs-unweighted P95 gap to ≈0
  (`explained% ≈ 100` everywhere).
- The overperformance is **sea-in-cell** (63–162% of it), not windier coastal land —
  coastal land is often *worse* than inland (FR/ES/AUS).

This is the shared root cause behind the symptom issues #24, #26, #27 (Peter's #41
comment: those close as symptoms once this lands).

## What the pipeline does today

CF weighting uses `cutout.indicatormatrix([geom])` — the **land-coverage fraction**
of each cell against the Natural-Earth region polygon (EPSG:4326). Consumers:

- `03_build_cf_timeseries.py` — country-average CF (per_unit over the matrix).
- `07_make_bestsite_cf_timeseries.py` — P95/best-site selection + resource summary.
- `07b_make_anchor_colocated_cf_timeseries.py` — anchor P95 (reuses 07).
- `03b` / `03c` — solar tilt-mix and candidate lattice (point-in-polygon mask).

Two gaps this leaves:
1. The land fraction correctly down-weights sea, but a coastal cell with, say, 23%
   land still carries its **inflated CF value** and can dominate `max` / the P95 tail.
2. There is no per-cell **capacity cap** (#31): nothing bounds how much capacity a
   mostly-sea cell can host.

## Proposed approach — PyPSA-Eur-style eligibility (approach 1)

atlite in our env already exposes the needed API (verified):
`Cutout.availabilitymatrix(shapes, ExclusionContainer)`.

An `ExclusionContainer` rasterises exclusion layers at a fine resolution (e.g. 100 m)
and `availabilitymatrix` returns, per (region, cell), the **eligible-area fraction**
— the buildable share after excluding sea + non-buildable land. This replaces the
raw land-coverage `indicatormatrix` as the weighting/eligibility source and yields:

- **Selection fix:** majority-sea and non-buildable cells get ~0 eligible area, so
  they drop out of P95/best-site/`max` instead of dominating the tail.
- **Capacity cap (#31):** `p_nom_max = eligible_fraction · cell_area · capacity_per_sqkm`
  (with `cell_area` in an equal-area CRS, i.e. the cos(lat) physical area — ties in
  with #37).

Note on the CF *value*: eligibility does **not** correct the ERA5-limited coastal CF
(that is approach 2 — finer/downscaled data). The multi-country decomposition argues
we do not need approach 2 just to recover coastal onshore resource, because coastal
land ≈ or < inland — there is little genuine coastal advantage being masked away.

## Data sources (must be global — we run EU + AUS + BRA)

PyPSA-Eur uses CORINE (Europe only), which won't cover AUS/BRA. Candidates:

| layer | purpose | global option | notes |
|---|---|---|---|
| land cover | exclude water/urban/forest classes | **ESA WorldCover 10 m** (2020/21) or Copernicus GLC 100 m | WorldCover is the modern global default |
| land–sea | exclude sea (the core #41 fix) | derived from land cover, or a coastline raster | MVP can start here alone |
| protected areas | exclude WDPA | **WDPA** (global) | optional, later phase |
| slope | exclude steep terrain | **Copernicus GLO-30 DEM** → slope | optional, later phase |

MVP needs only a **land–sea / water exclusion** (biggest lever for #41). Land-use,
protected areas and slope are refinements that improve realism and the #31 cap.

## Pipeline integration

New rule + script, mirroring the shapes/cutout pattern:

- `retrieve_land_cover` (or a committed/cached raster) → `data/land_cover/…`
- `build_availability`: inputs `cutout`, `regions`, land-cover raster; output
  `resources/availability/{cf_area}_availability.parquet` (or `.nc`) holding the
  per-cell eligible fraction + eligible area (m²) aligned to `cutout.grid` order.
- Consumers switch from `indicatormatrix` to the availability weights:
  - `03` weighting, `07`/`07b` P95 selection + summary, and the `max` column.
  - `p_nom_max` for #31 flows from the same eligible-area column.

Keep the weighting swap behind a small helper so all consumers share one definition
(this is where the #37 `_helpers.py` refactor — `geom_area_weights` etc. — is the
natural seam; see dependency note).

## Phasing

1. **MVP** — land–sea/water exclusion via `availabilitymatrix`; swap `03`/`07`/`07b`
   weights; re-run the decomposition to confirm coastal cells drop out and the P95
   gap collapses without a hand-picked border cutoff.
2. **Capacity cap (#31)** — emit `p_nom_max` from eligible area × capacity density;
   wire into `solve_network`.
3. **Refinements** — add land-use classes, WDPA, slope exclusions; sensitivity-check
   the P95/summary against MVP.

## Open decisions (for Hannah / Peter)

1. Land-cover source: ESA WorldCover 10 m (finer, larger) vs Copernicus GLC 100 m
   (lighter). WorldCover recommended.
2. Availability resolution (e.g. 100 m raster) and which land-use classes count as
   buildable — needs a defensible default per tech.
3. Capacity density (MW/km²) per tech for the #31 cap.
4. Whether to keep a hard minimum-eligibility threshold for selection in addition to
   weighting.

## Dependencies / sequencing

- The `_helpers.py` consolidation from **PR #55 (issue #37)** is the natural seam for
  the shared weighting swap (`geom_area_weights` → an eligibility-aware version). This
  branch is off `main` and does **not** yet include #55; rebase once #55 merges, or
  cherry-pick the helper module.
- Cross-links: #31 (capacity cap), #24/#26/#27 (symptoms), #37/#55 (cos-lat area,
  which the #31 `cell_area` term also needs).

## Validation

Reuse the notebook decomposition as the acceptance test: after the MVP, the coastal
border cells should carry ≈0 eligible weight, the P95 weighted-vs-unweighted gap
should be ~0 *without* the manual `w < 0.98` border cutoff, and the VIC `max` anomaly
(#24: 47.9% coastal boundary cell, 23% land) should resolve to a realistic buildable
value.
