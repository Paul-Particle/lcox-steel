# Land-eligibility weighting for onshore CF (issue #41)

Scoping doc. Branch: `land-eligibility-mask-issue41`.

**Scope:** this doc is about **#41 only** — fixing sea-contaminated onshore CF via a
land–sea eligibility mask. The related capacity cap (#31) is deliberately *not* solved
here; it is conceptualised in [Future work](#future-work--capacity-cap-31) so the #41
groundwork is built with it in mind, but the two are kept separate.

## Problem (recap)

Coarse ERA5 (~27–31 km) coastal cells blend smooth-sea and rough-land wind, so
onshore wind CF is inflated wherever a cell contains sea. Confirmed across all six
areas (FR/DE/ES/AUS/BRA/VIC) by the decomposition in
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

- `d1_area_average.py` — area-average CF (per_unit over the matrix).
- `d2_bestsite_p95.py` — P95/best-site selection + resource summary.
- `d3_anchor_colo.py` — anchor P95 (reuses 07).
- `03b` / `03c` — solar tilt-mix and candidate lattice (point-in-polygon mask).

The gap this leaves (and the subject of this doc): the land fraction correctly
down-weights sea, but a coastal cell with, say, 23% land still carries its **inflated
CF value** and can dominate `max` / the P95 tail. A low weight shrinks its pull in a
*mean*, but does not remove it from **P95 / best-site / `max` selection**.

(A separate gap — no per-cell **capacity cap** (#31) — is out of scope here. It is
not part of the #41 fix, but the eligibility groundwork below happens to be its
natural foundation, so it is sketched in [Future work](#future-work--capacity-cap-31)
rather than solved.)

## Approach — land-sea eligibility behind a shared seam (implemented)

The bug is a **selection** artefact: `indicatormatrix` already gives each cell's
land-coverage fraction and down-weights sea, but a mostly-sea coastal cell keeps its
inflated CF and can still win P95 / best-site / `max`. The notebook decomposition
showed the fix directly: classify cells on that existing land fraction
(`full = w ≥ 0.98·wmax`, `border = w < 0.98·wmax`) and drop the border cells — the
P95 gap collapses to ≈0 (`explained% ≈ 100`) in every area, with no new data.

So the MVP is a **minimum-land-fraction cutoff on the weight we already compute**, not
a new raster:

- **Selection fix:** cells below `min_land_fraction · wmax` (config `min_land_fraction`,
  default 0.98) are zeroed, so majority-sea coastal cells drop out of P95/best-site/
  `max` instead of dominating the tail — not merely down-weighted.
- **Onshore only:** offshore wind passes `min_land_fraction = 0` (it is sited on sea
  by design), so the cutoff is a no-op there.

**Pluggable source.** All of this lives behind one definition of per-cell eligibility
in `res_cf/_helpers.py`, and the land-fraction *source* is a config choice —
`res_cf.eligibility_source`, resolved by `cell_eligible_fraction`:

- `indicatormatrix` (default, shipped on): coarse NE-110m land/sea outline, no extra
  data — what the notebook validated.
- `availabilitymatrix` (**operational**): a finer coastline from atlite's
  `Cutout.availabilitymatrix(shapes, ExclusionContainer)` over a high-res land polygon
  (Natural Earth 10 m land). Excludes sea, rasterised in a global equal-area CRS. Same
  per-cell land fraction, sharper coastline.

The threshold rides on top of whichever source is chosen — the cutoff is identical; only
the fraction underneath changes. On DE the two sources agree to ~0.001 on mean/P95/max
(see Validation), so `indicatormatrix` stays the default and `availabilitymatrix` is the
cross-check / finer-coastline option. (Adding *land-use* / buildable-area exclusions on
top of the same machinery is a separate, later concern — see
[Future work](#future-work--capacity-cap-31).)

Note on the CF *value*: eligibility does **not** correct the ERA5-limited coastal CF
(that is approach 2 — finer/downscaled data). The multi-country decomposition argues
we do not need approach 2 just to recover coastal onshore resource, because coastal
land ≈ or < inland — there is little genuine coastal advantage being masked away.

## Data sources

- **`indicatormatrix` (default): no new data.** The land fraction comes from the
  Natural-Earth 110m outline we already ship — the notebook proved a threshold on that
  coarse coastline closes the #41 gap.
- **`availabilitymatrix`: Natural Earth 10 m land** (`ne_10m_land`, ~3 MB), a finer
  coastline. Manual download (gitignored like the other NE data — see README) to
  `data/shapes/ne_10m_land/ne_10m_land.zip`. Chosen over ESA WorldCover because the
  #41 question is land-vs-sea: a coastline vector is the right, light tool. WorldCover
  is a *land-use* product (open ocean is no-data), which is the wrong tool for the
  coastline and the right one only for the land-use refinements below.

The land-use exclusions that WorldCover / WDPA / DEM would provide alter *buildable
land*, not the coastal-CF mechanism, so they do not change the #41 outcome and are left
to [Future work](#future-work--capacity-cap-31):

| layer | purpose | global option | notes |
|---|---|---|---|
| land cover | exclude urban/forest classes | **ESA WorldCover 10 m** or Copernicus GLC 100 m | future #31 refinement |
| protected areas | exclude WDPA | **WDPA** (global) | future #31 refinement |
| slope | exclude steep terrain | **Copernicus GLO-30 DEM** → slope | future #31 refinement |

## Pipeline integration (as built)

No new Snakemake rule. A shared helper, a config block, and a small manual data file:

- **`res_cf/_helpers.py`** — one definition of per-cell eligibility:
  - `cell_eligible_fraction` dispatches on the source → `cell_land_fraction`
    (indicatormatrix) or `cell_availability_fraction` (availabilitymatrix: NE-10m land
    → `ExclusionContainer` → `availabilitymatrix`, land polygon cached + bbox-clipped).
  - `eligibility_fraction_1d` applies the threshold; `eligibility_weights` (07/07b
    (y,x) weights), `eligibility_matrix` (03's `(1, n_cells)` per_unit matrix), and
    `eligibility_mask_2d` (07b candidate validity) all build on it.
- **`config/config.yaml`** — `res_cf.min_land_fraction` (default 0.98; 0 disables),
  `res_cf.eligibility_source` (default `indicatormatrix`), and an `availability` block
  (land-shapes path, `crs: 6933`, `res: 200`) used only by `availabilitymatrix`.
- **Consumers**, onshore only (offshore passes 0):
  - `d1_area_average` — `eligibility_matrix` feeds the per_unit country mean;
    offshore matrix untouched.
  - `d2_bestsite_p95` — `eligibility_weights` per tech drives the
    national mean, P90/P95, `max`, and the P95-cell pick.
  - `d3_anchor_colo` — same cutoff on the anchor P95 pick and
    on onshore candidate validity.

## Plan (this branch)

Land–sea eligibility cutoff + a selectable finer-coastline source:

1. ✅ Eligibility helpers + `min_land_fraction` config; apply onshore in `03`/`07`/`07b`
   (offshore exempt).
2. ✅ `eligibility_source` selector; `availabilitymatrix` wired to NE-10m land via
   `ExclusionContainer` (`crs 6933`, `res 200`).
3. ✅ Verify across all six areas (see Validation): the cutoff collapses the onshore
   `max` tail everywhere (largest where coastline is largest; VIC #24 resolved), and the
   two sources agree — identical in 5/6 areas, DE differs only in P95 by 0.0014.

No capacity-cap work here — that is #31, sketched in [Future work](#future-work--capacity-cap-31).

## Open decisions (for Hannah / Peter)

1. `min_land_fraction` default: **0.98** (reproduces the notebook cutoff). Lower it if
   you want to keep more partial-land cells.
2. Whether the country **mean** (03) should also carry the cutoff, or only the
   selection paths (07/07b). Currently applied to both for consistency; it barely
   moves the mean.
3. Default `eligibility_source`: shipping `indicatormatrix` (no download, and it agrees
   with the finer source across all six areas — see Validation). `availabilitymatrix`
   stays available as the cross-check / finer-coastline option.

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

- The eligibility seam is self-contained in this branch's `res_cf/_helpers.py` — it
  does **not** depend on PR #55 (issue #37). If #55's `geom_area_weights`
  consolidation later lands on `main`, fold `eligibility_weights` into it on rebase so
  there is still one weighting definition; no blocking dependency either way.
- Cross-links: #24/#26/#27 (symptoms that close once #41 lands), #37/#55 (cos-lat
  area helper — independent), #31 (future capacity cap — see Future work, not built
  here).

## Validation

Acceptance test = the notebook decomposition: with the cutoff on, coastal border cells
carry 0 eligible weight and the P95 weighted-vs-unweighted gap collapses *without* a
hand-picked border cutoff (it is now the configured `min_land_fraction`).

### The fix across all six areas (wind-onshore, `min_land_fraction` off → on = 0 → 0.98)

`indicatormatrix` source; 2023 (VIC 2025):

| area | mean | P95 | max | eligible cells |
|---|---|---|---|---|
| de  | 0.199 → 0.198 | 0.283 → 0.267 | 0.453 → **0.380** | 827 → 673 |
| fr  | 0.183 → 0.179 | 0.292 → 0.277 | 0.477 → **0.363** | 1160 → 953 |
| es  | 0.109 → 0.107 | 0.201 → 0.192 | 0.440 → **0.277** | 945 → 772 |
| aus | 0.193 → 0.191 | 0.292 → 0.287 | 0.550 → **0.335** | 2988 → 2620 |
| bra | 0.056 → 0.054 | 0.155 → 0.141 | 0.411 → **0.278** | 3034 → 2687 |
| vic | 0.177 → 0.168 | 0.243 → 0.228 | 0.479 → **0.253** | 178 → 130 |

The mean barely moves (land-fraction weighting already down-weighted sea); the fix bites
on the **selection tail** — `max` and P95 — and the `max` drop scales with coastline
exposure: largest for AUS (−0.215), VIC (−0.227), ES (−0.163), smallest for inland-ish
DE. **The VIC `max` anomaly (#24) is resolved** — 0.479 → 0.253, 48 coastal cells
removed. Offshore is untouched; full `03`/`07`/`07b` runs are clean under Snakemake for
both sources.

### Source cross-check (`indicatormatrix` vs `availabilitymatrix`, both at 0.98)

| area | indicatormatrix (mean / P95 / max) | availabilitymatrix (mean / P95 / max) |
|---|---|---|
| de  | 0.1978 / 0.2665 / 0.3800 | 0.1977 / 0.2651 / 0.3800 |
| fr  | 0.1787 / 0.2769 / 0.3627 | identical |
| es  | 0.1069 / 0.1916 / 0.2767 | identical |
| aus | 0.1914 / 0.2869 / 0.3348 | identical |
| bra | 0.0537 / 0.1409 / 0.2782 | identical |
| vic | 0.1679 / 0.2284 / 0.2525 | identical |

The finer coastline gives **identical** results to the coarse default in 5 of 6 areas
(DE differs only in P95, by 0.0014) — **including the heavily-coastal ES / AUS / VIC**,
which is exactly where a divergence could have appeared. It doesn't: the 0.98 cutoff
keeps only near-fully-land cells, and those are the same under either coastline (the
fine-vs-coarse difference lives in the coastal cells both sources reject). So the cheap
`indicatormatrix` threshold is validated as the default everywhere, not just inland DE.
