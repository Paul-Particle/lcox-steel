# lcox-steel: Levelized Cost of Hydrogen for DRI Steel

This project calculates the levelized cost of hydrogen (LCOH) for green-steel DRI facilities powered by captive renewables and/or a grid connection. It combines an atlite-based capacity factor pipeline, ENTSO-E/NEM market data downloads, and a PyPSA investment optimization model.

## Project structure

```
lcox-steel/
├── Snakefile               # End-to-end workflow (grid + res_cf pipelines)
├── config.yaml             # Pipeline config (dates, countries, CF parameters)
├── config/
│   ├── assumptions.yaml    # Techno-economic defaults (CAPEX, OPEX, WACC, lifetimes)
│   └── projects.yaml       # Project + scenario definitions for PyPSA runs
├── areas.csv               # ENTSO-E areas to download (enabled flag)
├── environment.yaml        # Conda environment (lcox-steel)
├── data/
│   ├── entsoe_cache/       # Raw ENTSO-E feather files (area/year-month/data_type)
│   ├── integrated/         # Integrated ENTSO-E data (one feather per data type)
│   ├── cutouts/            # Atlite ERA5 cutout files (gitignored)
│   ├── shapes/             # Geographic boundaries (gitignored — see below)
│   └── res_cf/
│       ├── quarterly/      # Per-tech CF time series by quarter (intermediate)
│       ├── annual/         # Per-tech CF time series full-year (intermediate)
│       ├── <cc>_cf_<year>.csv              # Combined national CF
│       ├── <cc>_cf_<year>_bestsite_p95.csv # P95 best-site profiles
│       ├── resource_spread_<year>.csv      # Spatial resource statistics
│       └── <cc>_complementarity_top<N>_<year>.csv
├── results/                # PyPSA optimization outputs (.nc + summary CSVs)
└── scripts/
    ├── grid/               # Grid data pipeline (ENTSO-E + NEM)
    ├── res_cf/             # Atlite capacity factor pipeline
    │   ├── check_external_data.py  # Validate required external files
    │   ├── build_regions.py
    │   ├── build_offshore_regions.py
    │   ├── make_cutouts.py
    │   ├── build_cf_timeseries.py
    │   ├── concat_quarters.py
    │   ├── combine_techs.py
    │   ├── resource_spread.py
    │   ├── make_bestsite_cf_timeseries.py
    │   ├── complementarity.py
    │   └── diag_*.py               # Diagnostic and QC scripts
    ├── h2_dri/             # PyPSA investment model
    │   ├── run.py          # CLI entry point
    │   ├── network.py      # PyPSA network builder
    │   ├── costs.py        # LCOH post-solve accounting
    │   └── sizing.py       # Electrolyser sizing + annuity factor utilities
    └── viz/
        └── utils.py        # Shared plotting helpers
```

## Setup

### 1. Conda environment

```bash
conda env create -f environment.yaml
conda activate lcox-steel
```

### 2. External data files

Two large geographic datasets must be downloaded manually before the res_cf pipeline can run. Run the check script to see what's missing:

```bash
python scripts/res_cf/check_external_data.py
```

**World EEZ v12** (Marine Regions)
- Download from: https://www.marineregions.org/eez.php (free registration required)
- Choose "World EEZ v12 (2023)" → Shapefile format
- Extract so that `data/shapes/eez/eez_v12.shp` exists
- Any v11 or v12 release works; check `ISO_TER1` and `POL_TYPE` columns are present

**Natural Earth 1:110m Admin-0 countries**
- Download from: https://www.naturalearthdata.com/downloads/110m-cultural-vectors/
- Extract to `data/shapes/ne_110m_admin_0_countries/`

### 3. API keys

**ENTSO-E**: email transparency@entsoe.eu with "Restful API access" in the subject. Add to `.env` in the repo root:
```
ENTSOE_API_KEY=<your-key>
```
The `.env` file is gitignored — never commit it.

### 4. ERA5 access (for atlite cutouts)

Register at https://cds.climate.copernicus.eu and configure `~/.cdsapirc` following the [atlite CDS setup instructions](https://atlite.readthedocs.io/en/latest/installation.html).

## Running the pipelines

Both pipelines are managed by Snakemake. `config.yaml` controls which countries and years are processed.

### Full workflow (dry-run first)

```bash
snakemake -n          # preview what would run
snakemake --cores 4   # execute
```

### Grid data only

```bash
snakemake data/processed_data.feather --cores 4
```

Downloads and integrates ENTSO-E data into `data/processed_data.feather` (MultiIndex columns `(area, metric)`, hourly UTC DatetimeIndex).

### res_cf pipeline for one country

```bash
snakemake "data/res_cf/de_cf_2023.csv" --cores 4
```

This chains: `build_regions` → `build_offshore_regions` → `make_cutout` (4 quarters, ERA5) → `build_cf_timeseries` → `concat_quarters` → `combine_techs`.

### PyPSA investment optimization

```bash
python scripts/h2_dri/run.py --project DE_2023_baseline --scenario dedicated_res
python scripts/h2_dri/run.py --project DE_2023_baseline   # all scenarios
```

Results are written to `results/<project_name>/` as `.nc` (full PyPSA network) and `_summary.csv` (LCOH + optimal capacities).

Edit `config/projects.yaml` to add projects and scenarios. Edit `config/assumptions.yaml` to change techno-economic defaults.

## Data formats

**Grid data** (`data/processed_data.feather`): pandas DataFrame with a UTC hourly DatetimeIndex and MultiIndex columns `(area_code, metric)`. Metrics include `price`, `load_actual`, `load_forecast`, `vre`, `generation`, `crossborder`.

**Capacity factors** (`data/res_cf/<cc>_cf_<year>.csv`): hourly CSV with a `time` column and columns `wind_onshore_cf`, `wind_offshore_cf`, `solar_cf`.

**Best-site profiles** (`data/res_cf/<cc>_cf_<year>_bestsite_p95.csv`): same format; P95 grid cell extracted directly from the Atlite CF grid with 3×3 spatial averaging for wind.

**Complementarity results** (`data/res_cf/<cc>_complementarity_top<N>_<year>.csv`): ranked triplets of (onshore, offshore, solar) grid cells with score, coincidence, correlation, distance, and coordinates.
