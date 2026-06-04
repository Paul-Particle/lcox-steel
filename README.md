# lcox-steel: Levelized Cost of Hydrogen for DRI Steel

This project calculates the levelized cost of hydrogen (LCOH) for green-steel DRI facilities powered by captive renewables and/or a grid connection. It combines an atlite-based capacity factor pipeline, ENTSO-E/NEM market data downloads, and a PyPSA investment optimization model.

## Project structure

```
lcox-steel/
├── workflow/               # Snakemake workflow (standard layout)
│   ├── Snakefile           # one configfile + sys.path + includes + rule all
│   ├── rules/
│   │   ├── _optional_shim.smk  # local stand-in for snakemake's optional() helper (not shipped yet)
│   │   ├── grid.smk        # ENTSO-E + NEM download/process rules
│   │   ├── res_cf.smk      # atlite CF pipeline (shapes → cutout → CF timeseries)
│   │   ├── h2_dri.smk      # PyPSA optimisation rule
│   │   └── viz.smk         # compile_report + plot rules
│   ├── scripts/
│   │   ├── grid/           # Grid data pipeline (ENTSO-E + NEM)
│   │   ├── res_cf/         # Atlite capacity factor pipeline
│   │   │   ├── build_regions.py
│   │   │   ├── build_offshore_regions.py
│   │   │   ├── download_cutout.py             # honours `_backup.nc` cache; see below
│   │   │   ├── build_res_cf_profile.py
│   │   │   ├── build_solar_tilt_mix_p95.py
│   │   │   ├── determine_resource_spread.py
│   │   │   ├── determine_bestsite_p95.py
│   │   │   ├── determine_complementarity.py   # WIP, not in active DAG
│   │   │   └── diag_*.py                      # Diagnostic and QC scripts
│   │   ├── h2_dri/         # PyPSA investment model
│   │   │   ├── build_network.py               # PyPSA network construction (pure)
│   │   │   ├── solve_network.py               # snakemake entrypoint: load → build → solve → write
│   │   │   └── _helpers.py                    # annuity factor + electrolyser sizing
│   │   ├── viz/            # Reporting + Plotly figures
│   │   │   ├── compile_report.py              # Post-solve LCOH accounting → CSV
│   │   │   ├── plot_capacity_bars.py          # Per-project scenario bar chart
│   │   │   ├── plot_cf_map.py                 # Spatial CF heatmap
│   │   │   └── _helpers.py                    # FCA Plotly template + colormap
│   │   └── tests/          # End-to-end smoke tests
│   └── common/             # Shared helpers (_paths.py, _stubs.py)
├── config/
│   ├── config.yaml         # Pipeline knobs (countries, turbine specs, CF parameters, FX rates)
│   ├── assumptions.yaml    # Base techno-economic assumptions (CAPEX, OPEX, WACC, lifetimes)
│   ├── assumptions_{project}_{scenario}.yaml  # Optional per-scenario overlay (file presence = on)
│   └── projects.csv        # Project + scenario definitions — one row per (project, scenario, tech)
├── environment.yaml        # Conda environment (lcox-steel)
├── data/                   # Raw / external / expensive (not produced by this repo)
│   ├── entsoe_cache/       # ENTSO-E monthly raw cache (area/year-month/data_type); all ignored except entsoe_bidding_zones.csv
│   ├── nem_cache/          # AEMO NEMOSIS cache; all ignored except NEM Registration and Exemption List.xlsx
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
│   ├── res_cf/{cf_area}_{tech}_country-average_{start}_{end}.parquet     # Hourly CF time series
│   └── res_cf/{cf_area}_solar_tilt-mix-n{N}_{start}_{end}.parquet    # Orientation-resolved CF sweep
├── cutouts/                # Atlite ERA5 cutout files (gitignored). Existing downloads can be
│                           # preserved across reruns by renaming to `<name>_backup.nc`;
│                           # download_cutout.py copies from the backup when present (caching hack).
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

> [!TIP]
> **Bypassing Git Hooks:** A local git hook is configured above to automatically remove email addresses from commit messages for privacy. If you ever need to bypass this hook for a single commit, you can use the `--no-verify` flag: `git commit --no-verify`

### 2. External data files

Two large geographic datasets must be downloaded manually and placed as ZIPs — the pipeline reads from them directly via geopandas, no extraction step needed.

**World EEZ v12** — https://www.marineregions.org/downloads.php (free registration). Choose "World EEZ v12 (2023)" → Shapefile. Save (or rename) the download as `data/shapes/offshore_zones/eez_v12.zip`. Any v11 or v12 works; needs `ISO_TER1` and `POL_TYPE` columns.

**Natural Earth 1:110m Admin-0 countries** — https://www.naturalearthdata.com/downloads/110m-cultural-vectors/. Save as `data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip`.

```bash
mkdir -p data/shapes/offshore_zones data/shapes/ne_110m_admin_0_countries
# then drop the two ZIPs into those directories with the names above
```

Both cache directories follow the same gitignore strategy: everything is ignored except one committed reference file per pipeline.

**ENTSO-E bidding zone registry** — `data/entsoe_cache/entsoe_bidding_zones.csv` is a hand-crafted list of recognised bidding zone codes (e.g. `DE_LU`, `FR`, `NO_1`). `retrieve_entsoe.py` validates the `area` wildcard against this file before making any API calls — an unrecognised code raises a `ValueError` immediately. The file is small and rarely changes; update it manually when ENTSO-E adds or retires a zone. (See `TODO.md` for a planned migration to derive the list automatically from the `entsoe` library's `Area` enum.)

**NEM Registration and Exemption List** — AEMO-maintained generator metadata, refreshed by AEMO every so often. A snapshot is **committed** to the repo at `data/nem_cache/NEM Registration and Exemption List.xlsx`. It's small (~1 MB) and AEMO's hosting is flaky enough — they 403 NEMOSIS's default User-Agent, hence the monkey-patch in `_nemosis_patches.py` — that a checked-in copy makes bootstrapping much more reliable than relying on the live download. If the file is missing, `download_nem.py` falls back to a NEMOSIS static download (which inherits the same UA patch). If that also fails, fetch manually from https://www.aemo.com.au/-/media/Files/Electricity/NEM/Participant_Information/NEM-Registration-and-Exemption-List.xls (served as XLSX despite the extension) and drop it into `data/nem_cache/` as either `.xls` or `.xlsx`.

### 3. API keys

**ENTSO-E**: email transparency@entsoe.eu with "Restful API access" in the subject. Add to `.env` in the repo root:
```
ENTSOE_API_KEY=<your-key>
```
The `.env` file is gitignored — never commit it.

### 4. ERA5 access (for atlite cutouts)

Register at https://cds.climate.copernicus.eu and configure `~/.cdsapirc` following the [atlite CDS setup instructions](https://atlite.readthedocs.io/en/latest/installation.html).

> [!TIP]
> After a successful cutout download, copy `cutouts/{name}.nc` to `cutouts/{name}_backup.nc` to pin it across code-trigger reruns. `download_cutout.py` copies from the backup when present and skips CDS — handy until the proper cache lands (see `TODO.md`).

## Running the pipelines

All pipelines are managed by Snakemake. `config/projects.csv` drives the full DAG — each row is one `(project, scenario, tech)` input, and adding a project means adding rows, not editing the Snakefile.

### Full workflow (dry-run first)

```bash
snakemake -n                                            # preview what would run
snakemake --profile profiles/default --cores 4          # execute (concise output)
snakemake --profile profiles/default --cores 4 --verbose   # loud
```

The `profiles/default/` profile sets sensible defaults (`keep-going`, `rerun-incomplete`,
quiet rule chatter so per-rule log files carry the detail). Drop `--profile profiles/default`
to fall back to the bare Snakemake invocation.

### Grid data for one area

```bash
# ENTSO-E day-ahead prices (e.g. DE_LU, 2023):
snakemake resources/entsoe/DE_LU_grid_dayahead_20230101_20231231.parquet --cores 4 --resources entsoe_api=4

# NEM day-ahead prices (e.g. VIC1, 2025):
snakemake resources/nem/VIC1_grid_dayahead_20250101_20251231.parquet --cores 4
```

The `retrieve_entsoe` rule maintains a per-month raw cache in `data/entsoe_cache/`; once a month is cached it is never re-fetched. `--resources entsoe_api=N` caps concurrent ENTSO-E API calls (their rate limit is ~400/min per key).

To force a refresh, delete the relevant cache files then re-run:

```bash
rm -rf data/entsoe_cache/DE_LU/2024-12   # one month for one area
rm -rf data/entsoe_cache/*/2024-12       # one month for all areas
```

Months that fail (transient ENTSO-E errors, network blips) are retried 3× with exponential backoff inside the rule. Months that still fail are logged and skipped; the rule writes partial output and fails only if zero months succeeded. The next run re-attempts the missing months automatically.

### res_cf pipeline for one area + tech

```bash
# wind onshore CF for Germany, 2023:
snakemake resources/res_cf/de_wind-onshore_country-average_20230101_20231231.parquet --cores 4
```

This chains: `build_regions` → `build_offshore_regions` → `download_cutout` (ERA5) → `build_res_cf_profile`. The `{tech}` wildcard accepts `wind-onshore`, `wind-offshore`, or `solar`. Tech and variant identifiers use `-` as their internal separator (e.g. `wind-onshore`, `country-average`, `tilt-mix-n7`) so that the underscore-delimited filename pattern `{area}_{tech}_{variant}_{start}_{end}.parquet` stays unambiguous to parse — area can contain `_` (e.g. `DE_LU`), dates are always the last two 8-digit tokens.

> [!NOTE]
> **Current limitations (WIP).** The geometry is computed twice: `build_regions` produces the country's onshore geometry as a parquet for `build_res_cf_profile`, and `download_cutout` independently re-derives the same country boundary from the raw ZIP to compute the ERA5 bounding box. The bounding box is padded (`bbox_pad_deg` in config), so in practice it almost always encompasses the feasible offshore wind distance as well — the explicit offshore geometry from `build_offshore_regions` adds little that the cutout doesn't already cover spatially. The inconsistency is more subtle for offshore: the cutout bbox is a rectangle around the country's land area plus padding, while the offshore region is clipped to the EEZ within `offshore_max_distance_km`. For a narrow EEZ (many neighbours), the two align well. For a wide EEZ, the cutout covers only part of it — but `build_res_cf_profile` uses the full clipped offshore geometry as its spatial mask, so ERA5 cells far from the coast that fall outside the cutout are simply absent. Whether that matters depends on the country. A proper fix requires a cutout cache with explicit spatial and temporal coverage checking (see `TODO.md`).

### PyPSA investment optimization

Snakemake runs this automatically as part of `rule all`. For ad-hoc runs target a single network:

```bash
snakemake results/DE-2023-baseline/dedicated-res.nc --cores 4
```

Results are written to `results/<project>/<scenario>.nc` (full PyPSA network); the per-project
report CSV (`results/report_<project>.csv`, produced by `compile_report`) carries LCOH and
optimal capacities for every scenario in that project.

Edit `config/projects.csv` to add projects and scenarios. Edit `config/assumptions.yaml`
to change techno-economic defaults. To override a single scenario's assumptions, create
`config/assumptions_<project>_<scenario>.yaml` with just the keys you want to bump — the
script deep-merges it on top of the base. The overlay is picked up by file presence (no
column to edit in `projects.csv`).

Name projects and scenarios with dashes between words (`VIC-2025-solar-ew`,
`high-el-cost`), not underscores. Underscore is the field separator in the
filenames these names compose into (e.g. `assumptions_<project>_<scenario>.yaml`),
so reserving it for that keeps the boundaries legible at a glance. The `area`
column is the exception — it keeps underscores, since official bidding-zone codes
use them (`DE_LU`).

## Data formats

**Grid data** (`resources/entsoe/{area}_grid_dayahead_{start}_{end}.parquet`): pandas DataFrame with a UTC hourly DatetimeIndex and a single `price` column (EUR/MWh). The `_full` variant has MultiIndex columns `(area_code, metric)` covering all six ENTSO-E data types at native resolution. NEM equivalents live under `resources/nem/` with the same naming.

**Capacity factors** (`resources/res_cf/{cf_area}_{tech}_country-average_{start}_{end}.parquet`): hourly parquet, DatetimeIndex named `time`, a single column whose name *is* the tech key (e.g. `solar` or `wind-onshore`) with values in [0, 1]. One file per `(area, tech, date range)`; `{tech}` is `wind-onshore`, `wind-offshore`, or `solar`.

**Best-site profiles** (`resources/res_cf/{cf_area}_solar_tilt-mix-n{N}_{start}_{end}.parquet`): same time index, **multiple columns** — one per orientation in the sweep (e.g. `solar_az0`, `solar_az30`, … `solar_az330`). solve_network concatenates columns from all CF inputs into a single multi-tech `cf_timeseries`.

**Complementarity results** (`resources/res_cf/<cc>_complementarity_top<N>_<year>.parquet`): ranked triplets of (onshore, offshore, solar) grid cells with score, coincidence, correlation, distance, and coordinates.

## Logging & live output

Every rule writes a per-(rule, wildcards) log file under `logs/` via the Snakemake `log:`
directive. Scripts use Python's `logging` module wired through a shared helper
(`workflow/common/_logging.py`); there are no stray `print()` calls. Long-running loops
(triplet screening, per-cell extraction) show a `tqdm` progress bar with ETA when stderr
is a terminal, and degrade to a single start/finish log line otherwise.

### Layout

```
logs/
├── snakemake.log                                       # combined engine output (if you pipe through tee)
├── retrieve_entsoe/{area}_{variant}_{start}_{end}.log
├── retrieve_nem/{area}_{variant}_{start}_{end}.log
├── build_regions/{cf_area}.log
├── build_offshore_regions/{cf_area}.log
├── download_cutout/{cf_area}_{start}_{end}.log
├── build_res_cf_profile/{cf_area}_{tech}_{start}_{end}.log
├── h2_dri_optimize/{project}_{scenario}.log
└── compile_report/{project}.log
```

### Verbosity

Default is `INFO`. Crank up via either:

```bash
snakemake --profile profiles/default --cores 4 --config 'logging={level: DEBUG}'
# or just
snakemake --profile profiles/default --cores 4 --verbose
```

The third-party libraries that flood at INFO (atlite, cdsapi, urllib3, entsoe, pypsa,
linopy, fiona, matplotlib) are pinned to WARNING by default; bump
`config.logging.third_party_level` to lower their threshold.

### Watching a run live

The most agent-friendly recipe runs Snakemake in the background and tails everything:

```bash
mkdir -p logs
snakemake --profile profiles/default --cores 4 > logs/snakemake.log 2>&1 &
tail -F logs/snakemake.log logs/*/*.log
```

For a single rule, `tail -F logs/download_cutout/de_20230101_20231231.log` shows just
that one. CDS downloads via `download_cutout` log the start of each request and the
completion line; periodic status polling is on the roadmap (see `TODO.md`).

### HPC / cluster

`profiles/slurm/config.yaml.template` is a placeholder for a SLURM executor profile.
Copy it to `profiles/slurm/config.yaml`, fill in the `cluster:` line for your scheduler,
and run with `snakemake --profile profiles/slurm --jobs 200`. The script-side logging is
unchanged — per-rule files still land under `logs/{rule}/...`.
