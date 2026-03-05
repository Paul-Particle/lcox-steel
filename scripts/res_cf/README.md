# RES capacity factors (Atlite-based)

This folder contains **hourly, per-unit (0–1) renewable energy capacity factor (CF) time series**
derived using **Atlite + ERA5** for use in energy system and techno-economic models
(e.g. LCOX, electrolyser optimisation).

The pipeline has been validated and executed for the following countries (year 2023):

- Germany (DE)  
- France (FR)  
- Spain (ES)  
- Australia (AUS)  
- Brazil (BRA)

All countries use identical methodological and technology assumptions to ensure
cross-country comparability.

---

## Outputs (authoritative)

**Primary modelling inputs (per country):**

- `de_cf_2023.csv`
- `fr_cf_2023.csv`
- `es_cf_2023.csv`
- `aus_cf_2023.csv`
- `bra_cf_2023.csv`

Each file contains:

- `time` (hourly, full calendar year)
- `wind_onshore_cf`
- `solar_cf`

Each file contains **8760 hourly values** (non-leap year) and is intended for
downstream modelling.

---

## Method summary

1. **Climate data**
   - ERA5 reanalysis via Atlite
   - Time chunking for numerical/CDS stability:
     - DE: half-year cutouts (H1/H2)
     - FR/ES/AUS/BRA: quarterly cutouts (Q1–Q4)
   - Large-country grid coarsening for CDS job stability and runtime:
     - AUS/BRA cutouts use `dx=0.5`, `dy=0.5`
     - EU countries use default resolution (no `dx/dy` override)

2. **Spatial aggregation**
   - Country geometry from `regions.geojson`
   - Area-weighted aggregation using Atlite **indicator matrix**
   - Resulting CFs represent **national, climate-driven availability**

3. **Technology assumptions**
   - **Onshore wind:** single representative turbine  
     (`Vestas_V112_3MW`)
   - **Solar PV:** crystalline silicon (`CSi`), latitude-optimal tilt
   - Same technology assumptions used across countries for comparability

4. **Temporal handling**
   - CFs computed hourly per cutout segment (H1/H2 or Q1–Q4)
   - Pure concatenation to full year (no temporal averaging)
   - Aggregation (monthly/annual) performed only for QC or reporting

---

## Validation & QC (Multi-country 2023)

The pipeline passed the following checks for all countries:

- Structural integrity:
  - 8760 hourly steps
  - correct start/end timestamps
  - no gaps or duplicates
- Physical plausibility:
  - CFs bounded in [0,1]
  - realistic annual mean values by country
  - expected seasonal profiles
- Seam checks across concatenated segments showed no discontinuities

---

## Interpretation note

Capacity factors represent **country-level, area-weighted technical availability**.
They are suitable for cross-country comparison and system-level modelling.

Project-level LCOX at high-resource sites (e.g., strong wind or solar corridors)
may yield materially higher CFs than these national averages.
These values should therefore not be interpreted as site-specific project benchmarks.

---

## Scripts (pipeline overview)

- `scripts/res_cf/01_build_regions.py`  
  Builds `data/shapes/regions.geojson` from Natural Earth country geometries.

- `scripts/res_cf/02_make_cutouts.py`  
  Creates ERA5 cutouts (`data/cutouts/*.nc`) for a given country and time segment
  (H1/H2 or Q1–Q4). For large domains (AUS/BRA), uses coarser grid (`dx=0.5`, `dy=0.5`)
  to avoid CDS job failures and reduce runtime/memory.

- `scripts/res_cf/03_build_cf_timeseries.py`  
  Computes per-unit hourly CF time series (wind onshore + solar) from a single cutout
  using indicator-matrix aggregation over the country geometry. Outputs per-segment CSVs
  (e.g. `*_q1.csv`).

- `scripts/res_cf/04_concat_quaters.py`  
  Concatenates quarterly outputs (Q1–Q4) into full-year per-technology series:
  `*_wind_onshore_cf_2023.csv`, `*_solar_cf_2023.csv` (8760 hours).

- `scripts/res_cf/05_combine_techs.py`  
  Merges yearly wind + solar into the final modelling inputs:
  `*_cf_2023.csv` with columns `time`, `wind_onshore_cf`, `solar_cf`.

scripts/res_cf/06_resource_spread.py
Computes intra-country spatial resource distribution metrics from the existing Atlite cutouts.
For each country and technology:

derives annual mean CF per grid cell (no national aggregation)

calculates area-weighted spatial statistics (mean, P90, P95, max)

computes uplift factors relative to the national mean
Outputs:
data/res_cf/resource_spread_2023.csv

scripts/res_cf/07_make_bestsite_cf_timeseries.py
Generates “best-site” hourly CF time series by applying uplift factors (P90/P95) to the national hourly CFs.

multiplies national CFs by uplift factor

clips values to [0,1]
Outputs:
data/res_cf/<cc>_cf_2023_bestsite_p90.csv
data/res_cf/<cc>_cf_2023_bestsite_p95.csv

- `scripts/res_cf/99_qc_yearly_cf.py`  
  QC helper used during validation (structural + physical sanity checks).

---

## Archived intermediate files

For reproducibility/debugging, the following are retained (may be archived):

- ERA5 cutouts (`*.nc`) by segment (H1/H2 or Q1–Q4)
- Per-segment CF outputs (`*_h*.csv` and/or `*_q*.csv`)
- Per-technology yearly CFs (wind / solar separate)

These are **not** intended as direct model inputs.

---

## Notes / TODO (optional)

- Optional sensitivity check: rerun wind CFs with lower and higher
  specific-power turbines and compare annual mean CFs.
- Potential extension: add offshore wind using sea-mask or EEZ-based aggregation.

##Spatial resource spread & “best-site uplift” (2023 extension)

In addition to national area-weighted capacity factors, this pipeline quantifies intra-country spatial resource heterogeneity using the same Atlite cutouts (2023).

Purpose

Estimate the potential uplift achievable when projects are located in the strongest resource areas within a country (“best-site” assumption), while maintaining full methodological consistency with national CF calculations.

Method

For each country and technology (wind onshore, solar):

Hourly CFs are computed per grid cell (no aggregation).

Annual mean CF is calculated for each grid cell.

Area-weighted spatial statistics are derived using the same country geometry and indicator-matrix logic as in national aggregation:

Spatial mean (sanity check)

Spatial P90

Spatial P95

Spatial maximum

Uplift factors are computed relative to the national mean:

factor_p90 = spatial_p90_mean / national_mean

factor_p95 = spatial_p95_mean / national_mean

Outputs

data/res_cf/resource_spread_2023.csv
Contains spatial statistics and uplift factors per country and technology.

data/res_cf/<cc>_cf_2023_bestsite_p90.csv

data/res_cf/<cc>_cf_2023_bestsite_p95.csv

These contain hourly CF time series adjusted using uplift factors and clipped to [0,1].

Interpretation note

These “best-site” series represent a climate-resource upper envelope within national borders.

They:

Preserve hourly temporal structure

Scale national CFs by spatial resource strength

Remain fully consistent with the Atlite methodology

They do not include:

Land-use exclusions

Terrain or slope constraints

Grid access limitations

Permitting or zoning restrictions

They should therefore be interpreted as a resource-optimised scenario, not a fully constrained project benchmark.

Resolution note

Australia and Brazil use coarser cutout resolution (dx=0.5, dy=0.5).
Spatial statistics are computed on each country’s native cutout grid.

Scripts (pipeline extension)

scripts/res_cf/06_resource_spread.py
Computes intra-country spatial resource distribution metrics from existing Atlite cutouts.

Derives annual mean CF per grid cell

Calculates area-weighted spatial statistics (mean, P90, P95, max)

Computes uplift factors relative to the national mean
Output:
data/res_cf/resource_spread_2023.csv

scripts/res_cf/07_make_bestsite_cf_timeseries.py
Generates “best-site” hourly CF time series by applying uplift factors (P90/P95) to the national hourly CFs.

Multiplies national CFs by uplift factor

Clips values to [0,1]
Outputs:
data/res_cf/<cc>_cf_2023_bestsite_p90.csv
data/res_cf/<cc>_cf_2023_bestsite_p95.csv