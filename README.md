# lcox-steel: Levelized Cost of Hydrogen for DRI Steel

This project calculates the levelized cost of hydrogen (LCOH) for green-steel DRI facilities powered by captive renewables and/or a grid connection. It combines an atlite-based capacity factor pipeline, ENTSO-E/NEM market data downloads, and a PyPSA investment optimization model.

## Project structure

```
lcox-steel/
├── workflow/               # Snakemake workflow (standard layout)
│   ├── Snakefile           # Thin index: configfiles + includes + rule all
│   ├── rules/
│   │   ├── common.smk      # Constants, wildcard_constraints, helpers
│   │   ├── grid.smk        # download_entsoe, download_nem, process_entsoe
│   │   ├── res_cf.smk      # extract + atlite CF pipeline
│   │   └── h2_dri.smk      # PyPSA optimisation rule
│   ├── scripts/
│   │   ├── grid/           # Grid data pipeline (ENTSO-E + NEM)
│   │   ├── res_cf/         # Atlite capacity factor pipeline
│   │   │   ├── extract_shapefile.py # Generic zip→shp extractor
│   │   │   ├── build_regions.py
│   │   │   ├── build_offshore_regions.py
│   │   │   ├── make_cutout.py
│   │   │   ├── build_cf_timeseries.py
│   │   │   ├── concat_quarters.py
│   │   │   ├── combine_techs.py
│   │   │   ├── resource_spread.py
│   │   │   ├── make_bestsite_cf.py
│   │   │   ├── complementarity.py
│   │   │   └── diag_*.py               # Diagnostic and QC scripts
│   │   ├── h2_dri/         # PyPSA investment model
│   │   │   ├── run.py      # CLI entry point
│   │   │   ├── network.py  # PyPSA network builder
│   │   │   ├── costs.py    # LCOH post-solve accounting
│   │   │   └── sizing.py   # Electrolyser sizing + annuity factor utilities
│   │   ├── viz/utils.py    # Shared plotting helpers
│   │   └── tests/          # End-to-end smoke tests
│   ├── notebooks/          # API exploration notebooks (entsoe, nem)
│   └── common/             # Shared helpers (_paths.py, _stubs.py)
├── config/
│   ├── config.yaml         # Snakemake pipeline config (dates, countries, CF parameters)
│   ├── assumptions.yaml    # Techno-economic defaults (CAPEX, OPEX, WACC, lifetimes)
│   └── projects.yaml       # Project + scenario definitions for PyPSA runs
├── environment.yaml        # Conda environment (lcox-steel)
├── data/                   # Raw / external / expensive (not produced by this repo)
│   ├── entsoe_cache/       # Internal ENTSO-E monthly cache (area/year-month/data_type) — gitignored
│   ├── nem_cache/          # Internal AEMO NEMOSIS cache — gitignored
│   └── shapes/             # Raw shapefiles: ne_110m_admin_0_countries/, eez/ (gitignored — see below)
├── resources/              # Derived / Snakemake-tracked outputs (reproducible)
│   ├── entsoe/{area}/{data_type}.feather  # Stitched per-area, per-data_type ENTSO-E series
│   ├── entsoe_processed.feather             # Wide hourly grid dataset (MultiIndex columns)
│   ├── nem_processed.feather              # Australian NEM equivalent
│   ├── shapes/             # Derived geojsons (regions, offshore_regions — committed in git)
│   └── res_cf/
│       ├── quarterly/      # Per-tech CF time series by quarter (intermediate)
│       ├── annual/         # Per-tech CF time series full-year (intermediate)
│       ├── <cc>_cf_<year>.csv              # Combined national CF
│       ├── <cc>_cf_<year>_bestsite_p95.csv # P95 best-site profiles
│       ├── resource_spread_<year>.csv      # Spatial resource statistics
│       └── <cc>_complementarity_top<N>_<year>.csv
├── cutouts/                # Atlite ERA5 cutout files (gitignored)
├── .atlite-cache/          # Atlite scratch working dir (gitignored)
└── results/                # PyPSA optimization outputs (.nc + summary CSVs)
```

Snakemake commands are invoked from this repo root — `snakemake` auto-discovers `workflow/Snakefile`.

## Setup

### 1. Conda environment

```bash
conda env create -f environment.yaml
conda activate lcox-steel
```

### 2. External data files

Two large geographic datasets must be downloaded manually. Snakemake has two `extract_*_shapefile` rules that handle unzipping — just put the ZIPs at the canonical paths below (create the directories first) and the pipeline will extract them when it needs them.

**World EEZ v12** — https://www.marineregions.org/downloads.php (free registration). Choose "World EEZ v12 (2023)" → Shapefile. Save (or rename) the download as `data/shapes/eez/eez_v12.zip`. Any v11 or v12 works; needs `ISO_TER1` and `POL_TYPE` columns.

**Natural Earth 1:110m Admin-0 countries** — https://www.naturalearthdata.com/downloads/110m-cultural-vectors/. Save as `data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip`.

```bash
mkdir -p data/shapes/eez data/shapes/ne_110m_admin_0_countries
# then drop the two ZIPs into those directories with the names above
```

(A shapefile is a `.shp` + `.shx` + `.dbf` + `.prj` suite that has to travel together — the extract rule handles flattening if the ZIP nests the components under a subfolder.)

### 3. API keys

**ENTSO-E**: email transparency@entsoe.eu with "Restful API access" in the subject. Add to `.env` in the repo root:
```
ENTSOE_API_KEY=<your-key>
```
The `.env` file is gitignored — never commit it.

### 4. ERA5 access (for atlite cutouts)

Register at https://cds.climate.copernicus.eu and configure `~/.cdsapirc` following the [atlite CDS setup instructions](https://atlite.readthedocs.io/en/latest/installation.html).

## Running the pipelines

Both pipelines are managed by Snakemake. `config/config.yaml` controls which countries and years are processed.

### Full workflow (dry-run first)

```bash
snakemake -n          # preview what would run
snakemake --cores 4   # execute
```

### Grid data only

```bash
snakemake resources/entsoe_processed.feather --cores 4 --resources entsoe_api=4
```

Downloads and stitches ENTSO-E data into `resources/entsoe_processed.feather` (MultiIndex columns `(area, metric)`, hourly UTC DatetimeIndex). The `download_entsoe` rule manages a per-month cache in `data/entsoe_cache/` internally; once a month is cached it is never re-fetched.

`--resources entsoe_api=N` caps concurrent ENTSO-E API calls (their rate limit is ~400/min globally per API key — without throttling, a cold-cache build can exhaust the quota).

To force a refresh, delete the relevant cache files then re-run snakemake:

```bash
rm -rf data/entsoe_cache/DE_LU/2024-12   # one month for one area
rm -rf data/entsoe_cache/*/2024-12       # one month for all areas
rm resources/entsoe/DE_LU/prices.feather # force the rule to re-run too
```

Months that fail (transient ENTSO-E errors, network blips) are retried 3× with exponential backoff inside the rule. Months that still fail are logged and skipped; the rule writes a partial output covering only the successful months and fails only if zero months succeeded. The next run will automatically re-attempt the missing months.

### res_cf pipeline for one country

```bash
snakemake "resources/res_cf/de_cf_2023.csv" --cores 4
```

This chains: `build_regions` → `build_offshore_regions` → `make_cutout` (4 quarters, ERA5) → `build_cf_timeseries` → `concat_quarters` → `combine_techs`.

### PyPSA investment optimization

```bash
python workflow/scripts/h2_dri/run.py --project DE_2023_baseline --scenario dedicated_res
python workflow/scripts/h2_dri/run.py --project DE_2023_baseline   # all scenarios
```

Results are written to `results/<project_name>/` as `.nc` (full PyPSA network) and `_summary.csv` (LCOH + optimal capacities).

Edit `config/projects.yaml` to add projects and scenarios. Edit `config/assumptions.yaml` to change techno-economic defaults.

## Data formats

**Grid data** (`resources/entsoe_processed.feather`): pandas DataFrame with a UTC hourly DatetimeIndex and MultiIndex columns `(area_code, metric)`. Metrics include `price`, `load_actual`, `load_forecast`, `res`, `generation`, `crossborder`.

**Capacity factors** (`resources/res_cf/<cc>_cf_<year>.csv`): hourly CSV with a `time` column and columns `wind_onshore_cf`, `wind_offshore_cf`, `solar_cf`.

**Best-site profiles** (`resources/res_cf/<cc>_cf_<year>_bestsite_p95.csv`): same format; P95 grid cell extracted directly from the Atlite CF grid with 3×3 spatial averaging for wind.

**Complementarity results** (`resources/res_cf/<cc>_complementarity_top<N>_<year>.csv`): ranked triplets of (onshore, offshore, solar) grid cells with score, coincidence, correlation, distance, and coordinates.
