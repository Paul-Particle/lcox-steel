# lcox-steel: Levelized Cost of Hydrogen for DRI Steel

This project calculates the levelized cost of hydrogen (LCOH) for green-steel DRI facilities powered by captive renewables and/or a grid connection. It combines an atlite-based capacity factor pipeline, ENTSO-E/NEM market data downloads, and a PyPSA investment optimization model.

## Project structure

```
lcox-steel/
├── workflow/               # Snakemake workflow (standard layout)
│   ├── Snakefile           # configfiles + sys.path + includes + rule all
│   ├── rules/
│   │   ├── grid.smk        # ENTSO-E + NEM download/process rules
│   │   ├── res_cf.smk      # extract shapefiles + atlite CF pipeline
│   │   └── h2_dri.smk      # PyPSA optimisation rule
│   ├── scripts/
│   │   ├── grid/           # Grid data pipeline (ENTSO-E + NEM)
│   │   ├── res_cf/         # Atlite capacity factor pipeline
│   │   │   ├── extract_shapefile.py # Generic zip→shp extractor
│   │   │   ├── build_regions.py
│   │   │   ├── build_offshore_regions.py
│   │   │   ├── make_cutout.py
│   │   │   ├── build_cf_timeseries.py
│   │   │   ├── resource_spread.py
│   │   │   ├── determine_bestsite_p95.py
│   │   │   ├── determine_complementarity.py
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
│   ├── config.yaml         # Pipeline knobs (countries, turbine specs, CF parameters)
│   ├── assumptions.yaml    # Techno-economic defaults (CAPEX, OPEX, WACC, lifetimes)
│   └── projects.csv        # Project + scenario definitions — one row per (project, scenario, tech)
├── environment.yaml        # Conda environment (lcox-steel)
├── data/                   # Raw / external / expensive (not produced by this repo)
│   ├── entsoe_cache/       # Internal ENTSO-E monthly cache (area/year-month/data_type) — gitignored
│   ├── nem_cache/          # AEMO NEMOSIS cache — CSV/parquet gitignored; one committed XLSX (see §External data files)
│   └── shapes/             # Raw shapefiles: ne_110m_admin_0_countries/, eez/ (gitignored — see below)
├── resources/              # Derived / Snakemake-tracked outputs (reproducible)
│   ├── entsoe/{area}/{data_type}_{start}_{end}.parquet  # Raw per-type ENTSO-E downloads
│   ├── entsoe/{area}_grid_dayahead_{start}_{end}.parquet  # Hourly day-ahead price (EUR/MWh)
│   ├── entsoe/{area}_grid_full_{start}_{end}.parquet      # All 6 data types, native resolution
│   ├── nem/raw/{table}_{start}_{end}.parquet              # Raw AEMO MMSDM dispatch tables
│   ├── nem/{area}_grid_dayahead_{start}_{end}.parquet     # Hourly price (AUD→EUR, 1h mean)
│   ├── nem/{area}_grid_full_{start}_{end}.parquet         # All tables + derived columns, 5-min
│   ├── shapes/{cf_area}_geo.parquet                       # Onshore country geometry (GeoParquet)
│   ├── shapes/{cf_area}_offshore_geo.parquet              # EEZ-clipped offshore geometry
│   └── res_cf/{cf_area}_{tech}_country-average_{start}_{end}.parquet  # Hourly CF time series
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
git config core.hooksPath .githooks
```

### 2. External data files

Two large geographic datasets must be downloaded manually. Snakemake has two `extract_*_shapefile` rules that handle unzipping — just put the ZIPs at the canonical paths below (create the directories first) and the pipeline will extract them when it needs them.

**World EEZ v12** — https://www.marineregions.org/downloads.php (free registration). Choose "World EEZ v12 (2023)" → Shapefile. Save (or rename) the download as `data/shapes/offshore_zones/eez_v12.zip`. Any v11 or v12 works; needs `ISO_TER1` and `POL_TYPE` columns.

**Natural Earth 1:110m Admin-0 countries** — https://www.naturalearthdata.com/downloads/110m-cultural-vectors/. Save as `data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip`.

```bash
mkdir -p data/shapes/offshore_zones data/shapes/ne_110m_admin_0_countries
# then drop the two ZIPs into those directories with the names above
```

(A shapefile is a `.shp` + `.shx` + `.dbf` + `.prj` suite that has to travel together — the extract rule handles flattening if the ZIP nests the components under a subfolder.)

**NEM Registration and Exemption List** — AEMO-maintained generator metadata, refreshed by AEMO every so often. A snapshot is **committed** to the repo at `data/nem_cache/NEM Registration and Exemption List.xlsx`. It's small (~1 MB) and AEMO's hosting is flaky enough — they 403 NEMOSIS's default User-Agent, hence the monkey-patch in `download_nem.py` — that a checked-in copy makes bootstrapping much more reliable than relying on the live download. If the file is missing, `download_nem.py` falls back to `nemosis.static_table("Generators and Scheduled Loads", ...)` (which inherits the same UA patch). If that also fails, fetch manually from https://www.aemo.com.au/-/media/Files/Electricity/NEM/Participant_Information/NEM-Registration-and-Exemption-List.xls (served as XLSX despite the extension) and drop it into `data/nem_cache/` as either `.xls` or `.xlsx`. Caveat: when AEMO blocks NEMOSIS this way, the rest of the NEM downloader (MMSDM dispatch tables) is probably also broken via the same gating — better than nothing, but expect to debug.

### 3. API keys

**ENTSO-E**: email transparency@entsoe.eu with "Restful API access" in the subject. Add to `.env` in the repo root:
```
ENTSOE_API_KEY=<your-key>
```
The `.env` file is gitignored — never commit it.

### 4. ERA5 access (for atlite cutouts)

Register at https://cds.climate.copernicus.eu and configure `~/.cdsapirc` following the [atlite CDS setup instructions](https://atlite.readthedocs.io/en/latest/installation.html).

## Running the pipelines

All pipelines are managed by Snakemake. `config/projects.csv` drives the full DAG — each row is one `(project, scenario, tech)` input, and adding a project means adding rows, not editing the Snakefile.

### Full workflow (dry-run first)

```bash
snakemake -n          # preview what would run
snakemake --cores 4   # execute
```

### Grid data for one area

```bash
# ENTSO-E day-ahead prices (e.g. DE_LU, 2023):
snakemake resources/entsoe/DE_LU_grid_dayahead_20230101_20231231.parquet --cores 4 --resources entsoe_api=4

# NEM day-ahead prices (e.g. VIC1, 2025):
snakemake resources/nem/VIC1_grid_dayahead_20250101_20251231.parquet --cores 4
```

The `download_entsoe` rule manages a per-month cache in `data/entsoe_cache/`; once a month is cached it is never re-fetched. `--resources entsoe_api=N` caps concurrent ENTSO-E API calls (their rate limit is ~400/min per key).

To force a refresh, delete the relevant cache files then re-run:

```bash
rm -rf data/entsoe_cache/DE_LU/2024-12   # one month for one area
rm -rf data/entsoe_cache/*/2024-12       # one month for all areas
```

Months that fail (transient ENTSO-E errors, network blips) are retried 3× with exponential backoff inside the rule. Months that still fail are logged and skipped; the rule writes partial output and fails only if zero months succeeded. The next run re-attempts the missing months automatically.

### res_cf pipeline for one area + tech

```bash
# wind onshore CF for Germany, 2023:
snakemake resources/res_cf/de_wind_onshore_country-average_20230101_20231231.parquet --cores 4
```

This chains: `build_regions` → `build_offshore_regions` → `make_cutout` (ERA5) → `build_cf_timeseries`. The `{tech}` wildcard accepts `wind_onshore`, `wind_offshore`, or `solar`.

### PyPSA investment optimization

Snakemake runs this automatically as part of `rule all`. For ad-hoc runs:

```bash
python workflow/scripts/h2_dri/run.py --project DE_2023_baseline --scenario dedicated_res
python workflow/scripts/h2_dri/run.py --project DE_2023_baseline   # all scenarios
```

Results are written to `results/<project>/` as `.nc` (full PyPSA network) and `_summary.csv` (LCOH + optimal capacities).

Edit `config/projects.csv` to add projects and scenarios. Edit `config/assumptions.yaml` to change techno-economic defaults.

## Data formats

**Grid data** (`resources/entsoe/{area}_grid_dayahead_{start}_{end}.parquet`): pandas DataFrame with a UTC hourly DatetimeIndex and a single `price` column (EUR/MWh). The `_full` variant has MultiIndex columns `(area_code, metric)` covering all six ENTSO-E data types at native resolution. NEM equivalents live under `resources/nem/` with the same naming.

**Capacity factors** (`resources/res_cf/{cf_area}_{tech}_country-average_{start}_{end}.parquet`): hourly parquet, DatetimeIndex named `time`, single column `cf` in [0, 1]. One file per `(area, tech, date range)`; `{tech}` is `wind_onshore`, `wind_offshore`, or `solar`.

**Best-site profiles** (`resources/res_cf/<cc>_cf_<year>_bestsite_p95.parquet`): same format; P95 grid cell extracted directly from the Atlite CF grid with 3×3 spatial averaging for wind.

**Complementarity results** (`resources/res_cf/<cc>_complementarity_top<N>_<year>.parquet`): ranked triplets of (onshore, offshore, solar) grid cells with score, coincidence, correlation, distance, and coordinates.
