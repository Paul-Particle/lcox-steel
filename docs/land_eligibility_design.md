# Land-eligibility weighting for onshore CF (issue #41)

Scoping doc. Branch: `land-eligibility-mask-issue41`.

**Scope:** this doc is about **#41 only** — fixing sea-contaminated onshore CF via a
land–sea eligibility mask. The related capacity cap (#31) is deliberately *not* solved
here; it is conceptualised in [Future work](#future-work--capacity-cap-31) so the #41
groundwork is built with it in mind, but the two are kept separate.

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

The gap this leaves (and the subject of this doc): the land fraction correctly
down-weights sea, but a coastal cell with, say, 23% land still carries its **inflated
CF value** and can dominate `max` / the P95 tail. A low weight shrinks its pull in a
*mean*, but does not remove it from **P95 / best-site / `max` selection**.

(A separate gap — no per-cell **capacity cap** (#31) — is out of scope here. It is
not part of the #41 fix, but the eligibility groundwork below happens to be its
natural foundation, so it is sketched in [Future work](#future-work--capacity-cap-31)
rather than solved.)

## Proposed approach — PyPSA-Eur-style eligibility (approach 1)

atlite in our env already exposes the needed API (verified):
`Cutout.availabilitymatrix(shapes, ExclusionContainer)`.

An `ExclusionContainer` rasterises exclusion layers at a fine resolution (e.g. 100 m)
and `availabilitymatrix` returns, per (region, cell), the **eligible-area fraction**
— the buildable share after excluding sea. This replaces the raw land-coverage
`indicatormatrix` as the weighting/eligibility source and delivers the #41 fix:

- **Selection fix:** majority-sea cells get ~0 eligible area, so they drop out of
  P95/best-site/`max` selection instead of dominating the tail — not merely
  down-weighted.

For #41 the only exclusion layer needed is **land–sea**. The same eligible-area
output is also the natural input to a future per-cell capacity cap (#31); that reuse
is sketched in [Future work](#future-work--capacity-cap-31), not built here.

Note on the CF *value*: eligibility does **not** correct the ERA5-limited coastal CF
(that is approach 2 — finer/downscaled data). The multi-country decomposition argues
we do not need approach 2 just to recover coastal onshore resource, because coastal
land ≈ or < inland — there is little genuine coastal advantage being masked away.

## Data sources (must be global — we run EU + AUS + BRA)

PyPSA-Eur uses CORINE (Europe only), which won't cover AUS/BRA. Candidates:

| layer | purpose | global option | notes |
|---|---|---|---|
| land–sea | exclude sea (**the #41 fix**) | derived from land cover, or a coastline raster | **the whole scope of this branch** |
| land cover | exclude urban/forest classes | **ESA WorldCover 10 m** (2020/21) or Copernicus GLC 100 m | out of scope — future #31 refinement |
| protected areas | exclude WDPA | **WDPA** (global) | out of scope — future #31 refinement |
| slope | exclude steep terrain | **Copernicus GLO-30 DEM** → slope | out of scope — future #31 refinement |

The #41 fix needs only a **land–sea / water exclusion**. Land-use, protected areas
and slope do not change the #41 outcome (they alter *buildable land*, not the
coastal-CF mechanism); they matter only for the realism of a future capacity cap, so
they live in [Future work](#future-work--capacity-cap-31).

## Pipeline integration

New rule + script, mirroring the shapes/cutout pattern:

- `retrieve_land_cover` (or a committed/cached raster) → `data/land_cover/…`
- `build_availability`: inputs `cutout`, `regions`, land–sea raster; output
  `resources/availability/{cf_area}_availability.parquet` (or `.nc`) holding the
  per-cell eligible fraction aligned to `cutout.grid` order. (Also emit the eligible
  area in m² — cheap here, and the seam a future #31 cap would read; see Future work.)
- Consumers switch from `indicatormatrix` to the availability weights:
  `03` weighting, `07`/`07b` P95 selection + summary, and the `max` column.

Keep the weighting swap behind a small helper so all consumers share one definition
(this is where the #37 `_helpers.py` refactor — `geom_area_weights` etc. — is the
natural seam; see dependency note).

## Plan (this branch)

Single deliverable — the land–sea exclusion fix:

1. Add `retrieve_land_cover` + `build_availability` (land–sea raster → eligible
   fraction), swap `03`/`07`/`07b` from `indicatormatrix` to the availability weights.
2. Re-run the decomposition to confirm coastal cells drop out and the P95 gap
   collapses without a hand-picked border cutoff (see Validation).

No capacity-cap work here — that is #31, sketched in [Future work](#future-work--capacity-cap-31).

## Open decisions (for Hannah / Peter)

1. Land–sea source: ESA WorldCover 10 m (finer, larger) vs Copernicus GLC 100 m
   (lighter). WorldCover recommended. Only the water/sea class is used for #41.
2. Availability raster resolution (e.g. 100 m).
3. Whether to keep a hard minimum-eligibility threshold for selection in addition to
   weighting.

## Future work — capacity cap (#31)

Out of scope for this branch; recorded here so the #41 groundwork is built with the
follow-on in mind. **This is conceptualisation, not a commitment** — do not mix it
into the #41 change.

The eligibility work above produces a per-cell **eligible-area fraction**. A per-cell
capacity cap falls out of the same output:

```
p_nom_max = eligible_fraction · cell_area · capacity_density
```

- `cell_area` must be a true physical area (equal-area CRS / cos(lat)) — the same
  quantity #37 introduces, so #31 rides on the #37 helper.
- `capacity_density` (MW/km²) is a per-tech assumption we do not yet have a defensible
  default for — the main open question before #31 can be built.
- Wiring: emit `p_nom_max` from the availability output and feed it into
  `solve_network`.

**Refinements that belong to #31, not #41.** With only sea excluded, "eligible" still
counts urban/forest/protected/steep land as buildable, so the cap would be optimistic.
Adding land-use classes, WDPA and slope exclusions (the deferred data-source rows
above) tightens that realism. They do **not** affect the #41 CF fix, so they are only
worth doing once the #31 cap actually **binds** in the solve, or a reviewer wants
defensible buildable-area numbers. Open questions when that time comes: which land-use
classes count as buildable per tech, WDPA handling, and slope thresholds.

## Dependencies / sequencing

- The `_helpers.py` consolidation from **PR #55 (issue #37)** is the natural seam for
  the shared weighting swap (`geom_area_weights` → an eligibility-aware version). This
  branch is off `main` and does **not** yet include #55; rebase once #55 merges, or
  cherry-pick the helper module.
- Cross-links: #24/#26/#27 (symptoms that close once #41 lands), #37/#55 (cos-lat
  area helper), #31 (future capacity cap — see Future work, not built here).

## Validation

Reuse the notebook decomposition as the acceptance test: after the land–sea fix, the
coastal border cells should carry ≈0 eligible weight, the P95 weighted-vs-unweighted gap
should be ~0 *without* the manual `w < 0.98` border cutoff, and the VIC `max` anomaly
(#24: 47.9% coastal boundary cell, 23% land) should resolve to a realistic buildable
value.
