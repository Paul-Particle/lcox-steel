# lcox-steel: Levelized Cost of Hydrogen for DRI Steel

This project calculates the levelized cost of hydrogen (LCOH) for green-steel DRI facilities powered by captive renewables and/or a grid connection. It combines an atlite-based capacity factor pipeline, ENTSO-E/NEM market data downloads, and a PyPSA investment optimization model.

## Project structure

```
lcox-steel/
├── Snakefile               # Grid data download + processing workflow
├── config.yaml             # Snakemake pipeline config (dates, areas, data types)
├── config/
│   ├── assumptions.yaml    # Techno-economic defaults (CAPEX, OPEX, WACC, lifetimes)
│   └── projects.yaml       # Project + scenario definitions for PyPSA runs
├── areas.csv               # ENTSO-E areas to download (enabled flag)
├── environment.yaml        # Conda environment (lcox-steel)
├── data/
│   ├── entsoe_cache/       # Raw ENTSO-E feather files (area/year-month/data_type)
│   ├── integrated/         # Integrated ENTSO-E data (one feather per data type)
│   ├── res_cf/             # Atlite capacity factor outputs
│   └── processed_data.feather  # Final merged grid dataset
├── results/                # PyPSA optimization outputs (.nc + summary CSVs)
└── scripts/
    ├── grid/               # Snakemake-managed grid data pipeline
    │   ├── download_entsoe.py
    │   ├── download_nem.py
    │   ├── integrate_entsoe.py
    │   └── process_entsoe.py
    ├── res_cf/             # Atlite capacity factor pipeline
    │   ├── build_regions.py
    │   ├── build_offshore_regions.py
    │   ├── make_cutouts.py
    │   ├── build_cf_timeseries.py
    │   ├── concat_quarters.py
    │   ├── combine_techs.py
    │   ├── resource_spread.py
    │   ├── make_bestsite_cf_timeseries.py
    │   ├── complementarity.py      # Spatial site selection (wind/solar triplets)
    │   ├── diag_*.py               # Diagnostic and QC scripts
    │   └── README.md
    ├── h2_dri/             # PyPSA investment model
    │   ├── run.py          # CLI entry point
    │   ├── network.py      # PyPSA network builder
    │   ├── costs.py        # LCOH post-solve accounting
    │   └── sizing.py       # Electrolyser sizing + annuity factor utilities
    └── viz/
        └── utils.py        # Shared plotting helpers
```

## Setup

### 1. API key

Get an ENTSO-E Transparency Platform API key by emailing transparency@entsoe.eu with "Restful API access" in the subject. Create a `.env` file in the repo root:

```
ENTSOE_API_KEY=<your-key>
```

The `.env` file is gitignored — never commit it.

### 2. Conda environment

```bash
conda env create -f environment.yaml
conda activate lcox-steel
```

### 3. ERA5 access (for atlite)

Register at https://cds.climate.copernicus.eu and follow the atlite CDS setup instructions to configure `~/.cdsapirc`.

## Running the pipelines

### Grid data (ENTSO-E + NEM)

Configure download period and areas in `config.yaml` and `areas.csv`, then:

```bash
snakemake --cores 4
```

This downloads, integrates, and processes ENTSO-E data into `data/processed_data.feather` (MultiIndex columns `(area, metric)`, hourly UTC DatetimeIndex).

### Capacity factor pipeline (atlite)

Run scripts `01` through `08` in order from `scripts/res_cf/`. Each script reads config from `config_hannah.yaml`. Scripts `02` and `03` require ERA5 access and take significant time on first run; results are cached in `data/res_cf/`.

### PyPSA investment optimization

```bash
cd scripts/h2_dri
python run.py --project DE_2023_baseline --scenario dedicated_res
python run.py --project DE_2023_baseline   # all scenarios
```

Results are written to `results/<project_name>/` as `.nc` (full network) and `_summary.csv` (LCOH + capacities).

Edit `config/projects.yaml` to add projects and scenarios. Edit `config/assumptions.yaml` to change techno-economic defaults.

## Data formats

**Grid data** (`data/processed_data.feather`): pandas DataFrame with a UTC hourly DatetimeIndex and MultiIndex columns `(area_code, metric)`. Metrics include `price`, `load_actual`, `load_forecast`, `vre`, `generation`, `crossborder`.

**Capacity factors** (`data/res_cf/<cc>_cf_<year>.csv`): hourly CSV with a `time` column and one column per technology (`wind_onshore_cf`, `solar_cf`, etc.).

**Best-site profiles** (`data/res_cf/<cc>_cf_<year>_<variant>.csv`): same format, produced by script `07` for P95 best-site cells.

**Complementarity results** (`data/res_cf/<cc>_complementarity_top<N>_<year>.csv`): ranked triplets of (onshore, offshore, solar) grid cells with score, coincidence, correlation, and coordinates.
