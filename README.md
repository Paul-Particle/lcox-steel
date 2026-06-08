# lcox-steel: Levelized Cost of Hydrogen for DRI Steel

`lcox-steel` calculates the levelized cost of hydrogen (LCOH) for green-steel
direct-reduced-iron (DRI) facilities powered by captive renewables, a grid
connection, or both. It couples three data pipelines and an optimisation model
into one Snakemake workflow:

- **`res_cf`** — hourly renewable capacity factors from ERA5 reanalysis (via atlite).
- **`grid`** — hourly electricity market data from ENTSO-E (Europe) and NEM (Australia).
- **`h2_dri`** — a PyPSA investment model that sizes generation, storage, and the
  electrolyser to meet a constant hydrogen demand at least cost.
- **`viz`** — a per-project LCOH report plus Plotly figures.

## Architecture at a glance

Everything is a Snakemake workflow keyed off one table, `config/projects.csv`.
Each row is one `(project, scenario, tech)` input; `rule all` expands over every
row, so **adding a project or scenario is a CSV edit, not a Snakefile edit**.

```
config/projects.csv ─┬─► res_cf  ──►  resources/res_cf/*.parquet   (hourly CF)   ─┐
                     └─► grid    ──►  resources/{entsoe,nem}/*.parquet (€/MWh)   ─┤
                                                                                  ▼
                                                        h2_dri (PyPSA optimise) ──► results/{project}/{scenario}.nc
                                                                                  │
                                                          viz ◄────────────────────┘
                                              results/report_{project}.csv + plots/*.png|html
```

A scenario draws one or more capacity-factor series (one per RES tech) and/or a
single grid price series; `h2_dri` builds and solves the network; `viz` compiles
the LCOH report and charts.

## Project structure

```
lcox-steel/
├── workflow/                       # Snakemake workflow (standard layout)
│   ├── Snakefile                   # configfile + sys.path + includes + rule all
│   ├── rules/
│   │   ├── _optional_shim.smk      # local stand-in for Snakemake's optional() (not shipped yet)
│   │   ├── grid.smk                # ENTSO-E + NEM retrieval rules
│   │   ├── res_cf.smk              # atlite CF pipeline (shapes → cutout → CF series)
│   │   ├── h2_dri.smk              # PyPSA optimisation rule
│   │   └── viz.smk                 # compile_report + plot rules
│   ├── scripts/
│   │   ├── grid/                   # ENTSO-E + NEM download/process
│   │   │   ├── retrieve_entsoe.py  # rule entrypoint: warm-cache slice + on-miss download
│   │   │   ├── download_entsoe.py  # ENTSO-E per-month raw-cache primitives
│   │   │   ├── retrieve_nem.py     # rule entrypoint (NEM)
│   │   │   ├── download_nem.py     # NEMOSIS download primitives
│   │   │   ├── _nemosis_patches.py # AEMO User-Agent / URL-encoding workarounds
│   │   │   └── _helpers.py         # month iteration, UTC-naive coercion, cache checks
│   │   ├── res_cf/                 # atlite capacity-factor pipeline
│   │   │   ├── build_regions.py            # onshore country geometry (GeoParquet)
│   │   │   ├── build_offshore_regions.py   # EEZ-clipped offshore geometry
│   │   │   ├── download_cutout.py          # ERA5 cutout (honours `_backup.nc`; see below)
│   │   │   ├── build_res_cf_profile.py     # country-average hourly CF per tech
│   │   │   ├── build_solar_tilt_mix_p95.py # orientation-resolved solar CF sweep
│   │   │   ├── determine_*.py / diag_*.py  # WIP analysis & QC (NOT in the active DAG)
│   │   │   ├── _helpers.py                 # shared helpers for the WIP scripts
│   │   │   └── README.md                   # author's notes on the CF methodology (WIP)
│   │   ├── h2_dri/                 # PyPSA investment model
│   │   │   ├── build_network.py    # network construction (pure, importable)
│   │   │   ├── solve_network.py    # rule entrypoint: load → build → solve → write
│   │   │   └── _helpers.py         # annuity factor + electrolyser sizing
│   │   └── viz/                    # reporting + Plotly figures
│   │       ├── compile_report.py   # post-solve LCOH accounting → per-project CSV
│   │       ├── plot_capacity_bars.py  # per-project scenario capacity bar chart
│   │       ├── plot_cf_map.py      # spatial mean-CF heatmap with P95 site marked
│   │       └── _helpers.py         # FCA Plotly template + colormap (WIP, see TODO.md)
│   └── common/                     # shared, cross-pipeline Python
│       ├── _constants.py           # physical constants (e.g. H2 LHV)
│       ├── _logging.py             # configure_logging + tqdm progress wrapper
│       ├── _paths.py               # repo-relative path roots
│       └── _stubs.py               # snakemake object stub for linters/IDEs
├── config/
│   ├── config.yaml                 # pipeline knobs (logging, entsoe, nem, res_cf)
│   ├── assumptions.yaml            # base techno-economics (CAPEX, OPEX, WACC, lifetimes)
│   ├── assumptions_{project}_{scenario}.yaml   # optional per-scenario overlay (presence = on)
│   └── projects.csv                # one row per (project, scenario, tech)
├── profiles/
│   ├── default/config.yaml         # local-run defaults (keep-going, quiet, per-rule logs)
│   └── slurm/config.yaml.template  # HPC executor placeholder — copy & fill in
├── data/                           # raw / external / expensive (not produced here)
│   ├── entsoe_cache/               # ENTSO-E monthly raw cache (+ committed bidding-zone CSV)
│   ├── nem_cache/                  # NEMOSIS cache (+ committed AEMO registration list)
│   └── shapes/                     # raw shapefiles: Natural Earth, EEZ (see Setup)
├── resources/                      # derived, Snakemake-tracked outputs (reproducible)
├── cutouts/                        # atlite ERA5 cutouts (gitignored; see caching note)
├── .atlite-cache/                  # atlite scratch dir (gitignored)
├── results/                        # PyPSA networks (.nc), report CSVs, plots
├── environment.yaml                # conda environment (lcox-steel)
├── CLAUDE.md                       # project conventions (logging, Snakefile style)
└── TODO.md                         # roadmap / known WIP
```

Run Snakemake from the repo root — it auto-discovers `workflow/Snakefile`.

## Setup

### 1. Conda environment

```bash
conda env create -f environment.yaml
conda activate lcox-steel
git config core.hooksPath .githooks
```

> [!TIP]
> The `commit-msg` hook strips email addresses from commit messages for privacy.
> Bypass it for a single commit with `git commit --no-verify`.

### 2. External data files

Two geographic datasets must be downloaded manually and dropped in as ZIPs — the
pipeline reads them directly via geopandas, no extraction step needed.

**World EEZ v12** — https://www.marineregions.org/downloads.php (free registration).
Choose "World EEZ v12 (2023)" → Shapefile. Save as
`data/shapes/offshore_zones/eez_v12.zip`. Any v11/v12 works; it needs the
`ISO_TER1` and `POL_TYPE` columns.

**Natural Earth 1:110m Admin-0 countries** —
https://www.naturalearthdata.com/downloads/110m-cultural-vectors/. Save as
`data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip`.

```bash
mkdir -p data/shapes/offshore_zones data/shapes/ne_110m_admin_0_countries
# then drop the two ZIPs in with the names above
```

Two cache directories ship with one committed reference file each (everything
else is gitignored):

- **ENTSO-E bidding-zone registry** — `data/entsoe_cache/entsoe_bidding_zones.csv`
  is a hand-maintained list of recognised zone codes (`DE_LU`, `FR`, `NO_1`, …).
  `retrieve_entsoe` validates the `area` wildcard against it before any API call,
  raising on an unrecognised code. Update it by hand when ENTSO-E adds/retires a
  zone. (A planned migration derives it from the `entsoe` library's `Area` enum —
  see `TODO.md`.)
- **NEM Registration and Exemption List** —
  `data/nem_cache/NEM Registration and Exemption List.xlsx` is a committed AEMO
  snapshot (~1 MB). AEMO's hosting is flaky (it 403s NEMOSIS's default User-Agent,
  hence the patch in `_nemosis_patches.py`), so a checked-in copy makes
  bootstrapping reliable. If absent, `download_nem.py` falls back to a NEMOSIS
  static download; if that also fails, fetch it manually from
  https://www.aemo.com.au/-/media/Files/Electricity/NEM/Participant_Information/NEM-Registration-and-Exemption-List.xls
  (served as XLSX despite the extension) and drop it into `data/nem_cache/` as
  `.xls` or `.xlsx`.

### 3. API keys

**ENTSO-E** — email transparency@entsoe.eu with "Restful API access" in the
subject, then add the key to a gitignored `.env` in the repo root:

```
ENTSOE_API_KEY=<your-key>
```

### 4. ERA5 access (atlite cutouts)

Register at https://cds.climate.copernicus.eu and configure `~/.cdsapirc` per the
[atlite CDS setup instructions](https://atlite.readthedocs.io/en/latest/installation.html).

## Running the pipelines

`config/projects.csv` drives the full DAG. The first invocation should be a
dry-run; the `profiles/default` profile sets `keep-going`, `rerun-incomplete`,
and quiet rule chatter (per-rule log files carry the detail).

```bash
snakemake -n                                              # preview the DAG
snakemake --profile profiles/default --cores 4            # run everything
snakemake --profile profiles/default --cores 4 --verbose  # loud
```

Drop `--profile profiles/default` for the bare invocation.

### Targeting one output

```bash
# Grid — ENTSO-E day-ahead prices (DE_LU, 2023). Cap concurrent API calls:
snakemake resources/entsoe/DE_LU_grid_dayahead_20230101_20231231.parquet --cores 4 --resources entsoe_api=4

# Grid — NEM day-ahead prices (VIC1, 2025):
snakemake resources/nem/VIC1_grid_dayahead_20250101_20251231.parquet --cores 4

# res_cf — wind-onshore CF for Germany, 2023:
snakemake resources/res_cf/de_wind-onshore_country-average_20230101_20231231.parquet --cores 4

# h2_dri — one solved network:
snakemake results/DE-2023-baseline/dedicated-res.nc --cores 4
```

The ENTSO-E rule keeps a per-month raw cache under `data/entsoe_cache/`; once a
month is cached it is never re-fetched. Transient month failures are retried 3×
with backoff, then logged and skipped — the rule writes partial output and fails
only if *zero* months succeeded, and the next run re-attempts the gaps. Force a
refresh by deleting cache files:

```bash
rm -rf data/entsoe_cache/DE_LU/2024-12   # one month, one area
rm -rf data/entsoe_cache/*/2024-12       # one month, all areas
```

The `res_cf` chain is `build_regions` → `build_offshore_regions` →
`download_cutout` → `build_res_cf_profile`. The `{tech}` wildcard is
`wind-onshore`, `wind-offshore`, or `solar`.

> [!NOTE]
> **WIP — geometry computed twice.** `build_regions` produces the onshore
> geometry for `build_res_cf_profile`, while `download_cutout` independently
> re-derives the country boundary from the raw ZIP to set the ERA5 bounding box.
> The box is padded (`res_cf.cutout.bbox_pad_deg`), so it usually covers the
> feasible offshore band too. For a wide EEZ the cutout may not cover the full
> offshore region, but `build_res_cf_profile` masks to the clipped offshore
> geometry, so far-offshore cells outside the cutout are simply absent. A proper
> fix needs a cutout cache with explicit coverage checks (see `TODO.md`).

### Naming convention

Name projects and scenarios with **dashes** between words (`VIC-2025-solar-ew`,
`high-el-cost`), not underscores. The underscore is the field separator in the
filenames these names compose into (`assumptions_<project>_<scenario>.yaml`,
`logs/<project>_<scenario>.log`), so reserving it keeps the boundaries legible at
a glance. Tech and variant follow the same rule (`wind-onshore`,
`country-average`, `tilt-mix-n7`). The **`area`** column is the exception — it
keeps underscores, because official bidding-zone codes use them (`DE_LU`).

## Configuration

| File | Holds |
|------|-------|
| `config/config.yaml` | Pipeline knobs: `logging`, `entsoe` (data types), `nem` (`eur_per_aud` FX), `res_cf` (per-country metadata, turbines, CF flags, cutout settings). |
| `config/assumptions.yaml` | Base techno-economics: CAPEX/OPEX, lifetimes, WACC, electrolyser efficiency, plant sizing. Loaded by `h2_dri_optimize` as an **input file**, not a global `configfile:`. Tech keys (`res.wind-onshore`, `res.solar`, …) match the tech wildcard. |
| `config/assumptions_{project}_{scenario}.yaml` | *Optional* per-scenario overlay. **File presence is the toggle** (no CSV column); the `optional()` shim resolves it at job-evaluation time, and the script deep-merges it onto the base so the overlay carries only the keys it bumps. |
| `config/projects.csv` | Flat table, one row per `(project, scenario, tech)` input. Columns: `project, scenario, tech, variant, pipeline, area, start_date, end_date`. |

## Data formats

**Grid** (`resources/{entsoe,nem}/{area}_grid_dayahead_{start}_{end}.parquet`):
UTC hourly `DatetimeIndex`, single `price` column (EUR/MWh). The `_full` variant
has MultiIndex `(area, metric)` columns covering all data types at native
resolution.

**Capacity factors** (`resources/res_cf/{area}_{tech}_country-average_{start}_{end}.parquet`):
hourly parquet, `DatetimeIndex` named `time`, one column whose name *is* the tech
key (`solar`, `wind-onshore`, …) with values in [0, 1].

**Solar tilt-mix** (`resources/res_cf/{area}_solar_tilt-mix-n{N}_{start}_{end}.parquet`):
same index, **multiple columns** — one per orientation in the sweep (`solar_az0`,
`solar_az30`, …). `solve_network` concatenates columns from all CF inputs into one
multi-tech frame.

**Results**: `results/<project>/<scenario>.nc` is the full solved PyPSA network;
`results/report_<project>.csv` (from `compile_report`) carries LCOH and optimal
capacities for every scenario in the project.

### Cutout caching

ERA5 cutouts are expensive to (re-)download, so `cutouts/{name}.nc` is **not**
marked `protected()`. Instead, `download_cutout.py` copies from a sibling
`cutouts/{name}_backup.nc` when present, skipping CDS. To pin an existing cutout
across code-triggered reruns, rename it: `mv foo.nc foo_backup.nc`. This is a
stopgap until a content-addressed res_cf cache lands (see `TODO.md`).

## Logging & live output

Every rule writes a per-(rule, wildcards) log under `logs/` via the Snakemake
`log:` directive. Scripts log through `workflow/common/_logging.py` — no stray
`print()`s. Long loops show a `tqdm` bar with ETA on a TTY and degrade to a single
start/finish line otherwise.

```
logs/
├── retrieve_entsoe/{area}_{variant}_{start}_{end}.log
├── retrieve_nem/{area}_{variant}_{start}_{end}.log
├── build_regions/{cf_area}.log
├── build_offshore_regions/{cf_area}.log
├── download_cutout/{cf_area}_{start}_{end}.log
├── build_res_cf_profile/{cf_area}_{tech}_{start}_{end}.log
├── h2_dri_optimize/{project}_{scenario}.log
└── compile_report/{project}.log
```

Default level is `INFO`. Crank it up two ways:

```bash
snakemake --profile profiles/default --cores 4 --config 'logging={level: DEBUG}'
snakemake --profile profiles/default --cores 4 --verbose
```

Noisy upstream libraries (atlite, cdsapi, urllib3, entsoe, pypsa, linopy, fiona,
matplotlib) sit at `WARNING` by default; lower the threshold via
`logging.third_party_level` in `config.yaml`.

### Watching a run live

```bash
mkdir -p logs
snakemake --profile profiles/default --cores 4 > logs/snakemake.log 2>&1 &
tail -F logs/snakemake.log logs/*/*.log
```

### HPC / cluster

`profiles/slurm/config.yaml.template` is a SLURM executor placeholder. Copy it to
`profiles/slurm/config.yaml`, fill in the `cluster:` line for your scheduler, and
run `snakemake --profile profiles/slurm --jobs 200`. Script-side logging is
unchanged — per-rule files still land under `logs/{rule}/`.

## Conventions & roadmap

Project conventions (logging style, Snakefile/`.smk` rules, the two script
patterns) live in `CLAUDE.md`. Known WIP and planned work — the ENTSO-E zone-list
migration, a proper cutout cache, CDS download monitoring, and `viz/_helpers.py`
cleanup — are tracked in `TODO.md`.
