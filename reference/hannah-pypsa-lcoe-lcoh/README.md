# Reference: Hannah's PyPSA / LCOE–LCOH work (orientation only)

Verbatim copy of Hannah Maral's PyPSA-related scripts, taken from the
`lcox-steel-backup-full` repo on branch **`feature/pypsa-lcoe-lcoh`**.

**This is reference material for orientation and memory-jogging only.** The
files are copied *as they were* — they are **not** wired into this repo's
Snakemake pipeline, **not** modified to be runnable here, and **not**
maintained. Paths, imports, and inputs inside them point at the old repo
layout. Do not expect them to run.

## Provenance

- Source branch: `feature/pypsa-lcoe-lcoh`
- Relevant commits (author: Hannah Maral):
  - `139d697` (2026-05-21) — script 09 loop complete: 15 bestsite anchor runs across 5 countries
  - `97bf61e` (2026-05-18) — script 09 complete: single file working with oversizing ratios
  - `529a3d5` (2026-05-11) — stage3b WIP: res_to_h2 logic, storage dispatch, bestsite updates

## What's here

- `scripts/pypsa/` — the actual PyPSA LP prototype
  - `09_pypsa_optimise.py` — single electricity-bus LCOE optimiser; RES + battery
    sized against a constant electrolyser *load*; hardcoded placeholder costs.
  - `09b_battery_cost_sweep.py`, `check_cfs.py`, `08_complementarity_screen.py`
- `scripts/res_to_h2/` — a separate *non-PyPSA* analytical/parametric model
  (ratio sweeps, annuity LCOH, "Stage 3A storage proxy", firm-reliability
  diagnostics) plus the "stage3b" storage-dispatch experiments.
- `scripts/res_cf/` — Hannah's original res_cf pipeline (pre-Snakemake), the
  source the current `workflow/scripts/res_cf/` was re-architected from. Includes
  the quarterly-cutout workflow (`04_concat_quaters.py`, `05_combine_techs.py`)
  that the current pipeline dropped, plus QC/diagnostic helpers (`99_qc_yearly_cf.py`,
  `debug_cf.py`, `monthly_cf.py`). Useful for comparing the old vs current
  `03_build_cf_timeseries.py` and `07_make_bestsite_cf_timeseries.py`.

## Relationship to the current repo

The active PyPSA implementation lives in `workflow/scripts/h2_dri/`
(`build_network.py`, `solve_network.py`). It is a from-scratch reimplementation
that subsumes `09_pypsa_optimise.py` into a full two-bus (electricity + hydrogen)
LCOH model — electrolyser as a Link, H2 buffer, DRI demand, annuitised
config-driven costs, optional grid import, and multi-site HVDC siting. The
`res_to_h2/` parametric approach has no counterpart in the current pipeline.
