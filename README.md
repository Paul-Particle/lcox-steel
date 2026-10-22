# lcox-steel: Levelized Cost of Green Steel

`lcox-steel` calculates levelized costs for green-steel plants powered by
captive renewables, a grid connection, or both: the levelized cost of steel
(LCOS) of a process route. It couples three data pipelines and an
optimisation model into one Snakemake workflow:

- **`res_cf`** — hourly renewable capacity factors from ERA5 reanalysis (via atlite).
- **`grid`** — hourly electricity market data from ENTSO-E (Europe) and NEM (Australia).
- **`solve`** — a PyPSA investment model that sizes generation, storage, and one
  steelmaking process route to meet a constant demand at least cost. One solve
  is one run, building one route named by the `route` column of its scenario
  table:
  - `h2-dri-eaf` — electrolyser → H2 → DRI shaft → sponge iron (storable) → EAF → steel
  - `ng-dri-eaf` — fossil benchmark: natural-gas DRI (flat gas price, optional CO2 price) → EAF → steel
  - `mix-dri-eaf` — H2 and NG DRI shafts side by side; the optimiser picks the split
  - `moe-eaf` — molten oxide electrolysis → liquid iron → EAF → steel
  - `ew-eaf` — low-temperature iron electrowinning → iron plates (storable) → EAF → steel

  Every route also has an `-export` twin that makes its iron at the run's
  origin and melts it at its destination. A route without the suffix makes its
  steel at the origin and, if its run has a destination elsewhere, ships the
  steel there. The route ids live in `workflow/scripts/solve/build_network.py`.
- **`viz`** — an LCOS report per run and per scenario, plus Plotly figures.

## Architecture at a glance

Everything is a Snakemake workflow keyed off a scenario table, one row per input
series of a run. The default table is generated from the input workbook,
`data/assumptions/Assumptions_yaml_inputs.xlsx`, by `python -m
common.generate_inputs`; the demo's is kept by hand (see [Inputs](#inputs)).

A **scenario** is one row of the workbook — one origin, one destination, one
route, one power supply, one renewable setup, one date range — together with its
overlay. It is usually one **run**, one solved network. An area on the grid that
trades through market zones (Australia through its NEM regions) is run once per
zone, and the cheapest zone stands for the scenario. Nothing in the workflow reads
meaning out of a scenario or run name, so **adding a run is a table edit, not a
Snakefile edit**.

```
scenario table ──► res_cf ─┬─► resources/timeseries/*.parquet  (hourly CF) ─┐
               └─► grid   ─┘                        (and grid €/MWh)        │
                                                                            ▼
                         solve (PyPSA optimise) ──► results/{scenario}/.runs/{scenario}_{run}/outputs/*.nc
                                                                            │
                                                    viz ◄───────────────────┘
       results/{scenario}/report_{scenario}.csv  (the run that stands for the scenario)
       results/{scenario}/.runs/                  (every run, complete)
```

A run's rows name its input series: the origin's renewables, and its grid series
if it is on the grid, then the destination market's grid series if it exports.
The first area in a run's rows is its origin and the second its destination.
A run is named after its origin, and an export's after its destination too
(`AUS-DEU`, or `NSW1-DEU` for one of Australia's market zones).

## Project structure

```
lcox-steel/
├── workflow/                       # Snakemake workflow (standard layout)
│   ├── Snakefile                   # configfile + includes + rule all
│   ├── rules/
│   │   ├── grid.smk                # one retrieval rule, dispatching on the area's market
│   │   ├── res_cf.smk              # atlite CF pipeline (shapes → cutout → CF series)
│   │   ├── solve.smk               # PyPSA optimisation rule, one job per run
│   │   └── viz.smk                 # compile_report, pick_best_run + plot rules
│   ├── scripts/
│   │   ├── grid/                   # ENTSO-E + NEM download/process
│   │   │   ├── retrieve_grid_data.py # rule entrypoint: dispatches on the area's market
│   │   │   ├── _entsoe.py          # ENTSO-E: raw months → processed window
│   │   │   ├── download_entsoe.py  # ENTSO-E per-month raw-cache primitives
│   │   │   ├── _nem.py             # NEM: raw months → processed window
│   │   │   ├── download_nem.py     # NEMOSIS download primitives
│   │   │   ├── _nemosis_patches.py # AEMO User-Agent / URL-encoding workarounds
│   │   │   └── _helpers_grid.py    # month iteration, UTC-naive coercion, completeness guard
│   │   ├── res_cf/                 # atlite capacity-factor pipeline (numbered by stage)
│   │   │   ├── a_make_area_geometry.py        # onshore area geometry (GeoParquet)
│   │   │   ├── b_make_offshore_geometry.py    # EEZ-clipped offshore geometry
│   │   │   ├── b2_make_zone_offshore_geometry.py # a sub-national area's share of it, by nearest coast
│   │   │   ├── c_retrieve_area_cutout.py      # ERA5 cutout (honours `_backup.nc`; see below)
│   │   │   ├── d1_area_average.py             # area-average hourly CF per tech
│   │   │   ├── d2_bestsite_p95.py             # best-site P95 + anchored RES-mix CF
│   │   │   ├── d3_anchor_colo.py              # anchor co-located CF
│   │   │   ├── d4_tilt_mix.py                 # orientation-resolved solar CF sweep
│   │   │   ├── d5_multi.py                    # per-cell candidate grid (multi-site siting)
│   │   │   ├── spot_check_cutout.py           # standalone cutout inspection (NOT in the DAG)
│   │   │   ├── _helpers_res_cf.py             # shared helpers for the CF scripts
│   │   │   ├── reference/                     # Hannah's original scripts, verbatim (for side-by-side diffs)
│   │   │   └── README.md                      # author's notes on the CF methodology (WIP)
│   │   ├── solve/                  # PyPSA investment model
│   │   │   ├── build_network.py    # network construction (pure, importable)
│   │   │   ├── solve_network.py    # rule entrypoint: load → build → solve → write
│   │   │   └── _helpers_solve.py   # annuity factor, electrolyser sizing, a run's inputs
│   │   └── viz/                    # reporting + Plotly figures
│   │       ├── compile_report.py   # post-solve cost accounting → one run's report
│   │       ├── pick_best_run.py    # copies the run that stands for a scenario up into it
│   │       ├── outputs/            # result plots: capacity bars, LCOS breakdown, per-site capacity
│   │       ├── diagnostics/        # checks: CF heatmap, siting map, DE price vs residual demand
│   │       ├── dashboard/          # HTML dashboard builders + cost taxonomy
│   │       │   ├── templates/      # the pages' HTML templates
│   │       │   └── hub_pages/      # static hub pages (network schematic, workflow DAG)
│   │       ├── _run_display.py     # run labels, shared by the plots
│   │       └── style.py            # FCA Plotly template + colormap
│   └── common/                     # shared, cross-pipeline Python
│       ├── generate_inputs.py      # writes the scenario table, base and overlays from the input workbook
│       ├── _workbook_translation.py # the workbook's keys and vocabulary → the model's
│       ├── _report_schema.py       # the report's fields + the only reader/writer of its layout
│       ├── _provenance.py          # content digests of a solve's inputs (inputs_hash)
│       ├── _constants.py           # physical constants (e.g. H2 LHV)
│       ├── _logging.py             # configure_logging + tqdm progress wrapper
│       ├── _paths.py               # repo-relative path roots
│       └── _stubs.py               # snakemake object stub for linters/IDEs
├── config/
│   ├── config.yaml                 # pipeline knobs, the inputs to run, area registry
│   ├── demo.yaml                   # runs the template instead: `--configfile config/demo.yaml`
│   ├── .template.yaml              # complete assumptions: the demo's base, and the generator's
│   ├── .template.csv               # the demo's scenario table
│   ├── overlays/{scenario}.yaml    # optional per-scenario overlay (presence = on); generated ones gitignored
│   └── generated/                  # base assumptions + scenario table from the workbook (gitignored)
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
├── pyproject.toml                  # declares workflow/ as the installable package (see Setup)
└── CLAUDE.md                       # project conventions (logging, Snakefile style)
```

Run Snakemake from the repo root — it auto-discovers `workflow/Snakefile`.

## Setup

### 1. Conda environment

```bash
conda env create -f environment.yaml
conda activate lcox-steel
git config core.hooksPath .githooks
plotly_get_chrome -y        # once per machine: Chrome for Testing (~325 MB)
```

Static plot export (`write_image`, kaleido v1) needs a Chrome binary. Without one
every PNG/SVG output fails with "Kaleido requires Google Chrome to be installed";
a Chrome already in `/Applications` works as well. HTML outputs need nothing.

> [!TIP]
> The `commit-msg` hook strips email addresses from commit messages for privacy.
> Bypass it for a single commit with `git commit --no-verify`.

### 2. The package, installed per checkout

The code under `workflow/` is an installable package (`pyproject.toml`), with
two top-level packages: `common` and `scripts`. Every module imports the others
through them (`from common._paths import DATA`,
`from scripts.grid._helpers_grid import iso`), in Snakemake scripts, tests and
standalone runs alike. The dependencies are the conda environment's, so the
package itself is installed with `--no-deps`.

An editable install points at one checkout. **With a single checkout**, install
it straight into the environment:

```bash
pip install -e . --no-deps
```

**With several checkouts on one machine** (git worktrees, copies under
`~/lcox-runs`), give each its own small venv layered on the conda environment,
so each imports its own code. Installing all of them into the shared
environment would leave every checkout importing whichever was installed last.

```bash
conda activate lcox-steel     # the venv needs the env's PROJ_DATA / GDAL_DATA
uv venv --system-site-packages --python "$CONDA_PREFIX/bin/python" .venv
uv pip install --python .venv/bin/python -e . --no-deps
source .venv/bin/activate     # in each new shell, after conda activate
```

Inside such a venv, run Snakemake and pytest through the venv's Python:
`python -m snakemake …` and `python -m pytest`. The bare `snakemake` command
belongs to the conda environment and would run outside the venv. With the
single-checkout install, `snakemake` works as it is.

### 3. External data files

Geographic datasets must be downloaded manually and dropped in as ZIPs — the
pipeline reads them directly via geopandas, no extraction step needed.

**World EEZ v12** — https://www.marineregions.org/downloads.php (free registration).
Choose "World EEZ v12 (2023)" → Shapefile. Save as
`data/shapes/offshore_zones/eez_v12.zip`. Any v11/v12 works; it needs the
`ISO_TER1` and `POL_TYPE` columns.

**Natural Earth 1:110m Admin-0 countries** —
https://www.naturalearthdata.com/downloads/110m-cultural-vectors/. Save as
`data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip`.

**Natural Earth 1:10m Admin-1 states and provinces** *(for sub-national areas —
the NEM regions, Brazil's submarkets, Alberta and Ontario — whose geometry is the
dissolve of the states their `admin1_names` list; Admin 1 is not published at
1:110m)* — https://www.naturalearthdata.com/downloads/10m-cultural-vectors/. Save as
`data/shapes/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.zip`,
or fetch directly:

```bash
curl -sSL --create-dirs -o data/shapes/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.zip \
  https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip
```

Those areas' geometries are also committed under `resources/shapes/`, so without
this zip the workflow uses them as they are.

**Natural Earth 1:10m land** *(only for `res_cf.eligibility_source: availabilitymatrix`
— the finer-coastline land-sea mask; the default `indicatormatrix` does not need it)* —
https://www.naturalearthdata.com/downloads/10m-physical-vectors/. Save as
`data/shapes/ne_10m_land/ne_10m_land.zip`, or fetch directly:

```bash
curl -sSL -o data/shapes/ne_10m_land/ne_10m_land.zip \
  https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip
```

```bash
mkdir -p data/shapes/offshore_zones data/shapes/ne_110m_admin_0_countries \
  data/shapes/ne_10m_admin_1_states_provinces data/shapes/ne_10m_land
# then drop the ZIPs in with the names above
```

Two cache directories ship with one committed reference file each (everything
else is gitignored):

- **ENTSO-E bidding-zone registry** — `data/entsoe_cache/entsoe_bidding_zones.csv`
  is a hand-maintained list of recognised zone codes (`DE_LU`, `FR`, `NO_1`, …).
  `retrieve_entsoe` validates the `area` wildcard against it before any API call,
  raising on an unrecognised code. Update it by hand when ENTSO-E adds/retires a
  zone. (Deriving it from the `entsoe` library's `Area` enum instead is issue #17.)
- **NEM Registration and Exemption List** —
  `data/nem_cache/NEM Registration and Exemption List.xlsx` is a committed AEMO
  snapshot (~1 MB). AEMO's hosting is flaky (it 403s NEMOSIS's default User-Agent,
  hence the patch in `_nemosis_patches.py`), so a checked-in copy makes
  bootstrapping reliable. If absent, `download_nem.py` falls back to a NEMOSIS
  static download; if that also fails, fetch it manually from
  https://www.aemo.com.au/-/media/Files/Electricity/NEM/Participant_Information/NEM-Registration-and-Exemption-List.xls
  (served as XLSX despite the extension) and drop it into `data/nem_cache/` as
  `.xls` or `.xlsx`.

### 4. API keys

**ENTSO-E** — email transparency@entsoe.eu with "Restful API access" in the
subject, then add the key to a gitignored `.env` in the repo root:

```
ENTSOE_API_KEY=<your-key>
```

### 5. ERA5 access (atlite cutouts)

Register at https://cds.climate.copernicus.eu and configure `~/.cdsapirc` per the
[atlite CDS setup instructions](https://atlite.readthedocs.io/en/latest/installation.html).

## Running the pipelines

### Quick start — the demo (no credentials, no downloads)

`config/demo.yaml` runs the hand-kept template instead, a small self-contained
**demo** that runs after a fresh clone with **no CDS download, no
EEZ/Natural-Earth zips, and no API keys**. It ships a pre-sliced Victoria
(Australia) cutout at 0.5° for its own `DEMO` area, plus that area's geometry
parquets, and takes the best-site (`d2`) capacity-factor science through a
`solve` to a `viz` report. Its scenarios build `h2-dri-eaf` and `moe-eaf`, each
at home and as the export twin, which melts its iron in `ES_SYN`, a made-up
Spain-like market, and `moe-eaf` once more shipping its steel there, so their
destination figures are illustrative only:

```bash
snakemake --configfile config/demo.yaml --profile profiles/default --cores 4
```

Outputs land at `results/{scenario}/report_{scenario}.csv`, in a few minutes on a
laptop.

### The workbook scenarios

The default target builds every scenario generated from the input workbook. Write
them first; without `config/generated/`, Snakemake stops on the missing file:

```bash
python -m common.generate_inputs                          # writes config/generated/
snakemake -n                                              # preview the full DAG
snakemake --profile profiles/default --cores 4            # build every scenario
snakemake --profile profiles/default --cores 4 --verbose  # loud
```

Building them needs the external data and credentials described under
[Setup](#setup) — ERA5 cutouts via CDS, the EEZ/NE zips, and ENTSO-E, NEM or ONS
access for grid-connected scenarios and for any scenario that builds an `-export`
route, whose furnace stands in the destination market and prices its power there.
Drop `--profile profiles/default` for the bare invocation.

### Targeting one output

```bash
# Grid — ENTSO-E prices + generation mix for DEU (bidding zone DE_LU). Cap concurrent API calls:
snakemake resources/timeseries/DEU_grid_emissions_20250101_20251231.parquet --cores 4 --resources entsoe_api=4

# Grid — NEM prices + generation mix (VIC1, 2025):
snakemake resources/timeseries/VIC1_grid_emissions_20250101_20251231.parquet --cores 4

# res_cf — wind-onshore CF for Germany, 2023:
snakemake resources/timeseries/DEU_wind-onshore_area-average_20250101_20251231.parquet --cores 4

# solve — one run, results/{scenario}/.runs/{scenario}_{run}/outputs/{scenario}_{run}.nc:
snakemake results/demo-moe-eaf/.runs/demo-moe-eaf_DEMO/outputs/demo-moe-eaf_DEMO.nc --configfile config/demo.yaml --cores 4
```

The grid rule reads a per-month raw cache — `data/entsoe_cache/` for ENTSO-E,
NEMOSIS's own `data/nem_cache/` — and once a month is cached it is never
re-fetched. Processing those months into a window is redone per run, and the
rule output holds the result: it is keyed by `{area}_grid_{variant}_{dates}`, so
Snakemake reuses it on a re-run and a different date range simply writes a
different file. Transient month failures are retried 3×
with backoff, then logged and skipped — the rule writes partial output and fails
only if *zero* months succeeded, and the next run re-attempts the gaps. Force a
refresh by deleting cache files:

```bash
rm -rf data/entsoe_cache/DE_LU/2024-12   # one month, one area
rm -rf data/entsoe_cache/*/2024-12       # one month, all areas
```

The `res_cf` chain is `make_area_geometry` → `make_offshore_geometry` →
`retrieve_area_cutout` → `area_average`. The `{tech}` wildcard is
`wind-onshore`, `wind-offshore`, or `solar`. `make_area_geometry` builds an
area's land as its Natural Earth Admin 0 country (by `iso3`), or, for an area
whose registry entry lists `admin1_names`, as the dissolve of those Admin 1
states or provinces — that is how a NEM region or a Brazilian submarket gets its
outline. Such an area's offshore geometry comes from `make_zone_offshore_geometry`
instead: its country's offshore area split among the country's Admin 1 units by
nearest coast, so neighbouring areas never share a cell of sea.

> [!NOTE]
> **Cutout bounds = land ∪ offshore.** `retrieve_area_cutout` reads the pre-built
> `{area}_geo.parquet` and `{area}_offshore_geo.parquet`, unions them, and
> takes the bounding box padded by `res_cf.cutout.bbox_pad_deg`. Unioning the
> offshore zone in matters: it can reach `res_cf.offshore_max_distance_km`
> (~200 km) from the coast — well beyond a land-only bbox + 1° pad — so without
> it the offshore-wind cells get clipped out of the cutout. The mainland_bbox
> filter and EEZ clip are already applied upstream by `make_area_geometry` /
> `make_offshore_geometry` / `make_zone_offshore_geometry`, so the geometry is
> computed once. An area's `cutout_bbox` replaces the box altogether.

### Naming convention

Five identifiers — **route, tech, variant, area, dates** — have one spelling
each, set in the scenario tables, and every layer below spells them the same:
the assumptions blocks, the Snakemake wildcards and rule names, the script
filenames, and the PyPSA component names.

Name scenarios with **dashes** between words (`VIC-2025-solar-ew`,
`DE-2023-grid-p95`), not underscores. The underscore is the field separator in
the names these compose into (a run's folder is `<scenario>_<run>`), so reserving
it keeps the boundaries legible at a glance. Route, tech and variant follow the
same rule (`moe-eaf`, `wind-onshore`, `area-average`, `tilt-mix-n7`).

The **`area`** column is the exception — it keeps underscores, because official
market-zone codes use them (`DE_LU`). The registry is the top-level `areas:`
block in `config/config.yaml`, and nothing infers anything from the shape of a
code: what an area can be used for follows from its entry.

An area is a country (`DEU`, `FRA`, `ESP`, `AUS`, `BRA`) or a sub-national market
zone (`VIC1`), written as its registry key.

On the grid, an area with no market of its own trades through the `zones:` its
entry lists, so input generation runs a grid scenario there once per zone —
Australia on the grid yields its NEM regions separately, while an islanded
scenario keeps the whole country.

## Inputs

`inputs` in `config/config.yaml` names the base assumptions and the scenario table
`snakemake` runs; a scenario's overlay, if it has one, is
`config/overlays/{scenario}.yaml`. There are two sets of them:

| | base assumptions | scenario table | kept by |
|---|---|---|---|
| default | `config/generated/assumptions.yaml` | `config/generated/scenarios.csv` | `python -m common.generate_inputs`; gitignored |
| `--configfile config/demo.yaml` | `config/.template.yaml` | `config/.template.csv` | hand |

The overlays of both live side by side in `config/overlays/`. The generator's
start with a `# GENERATED` line, and it only ever writes over or removes those;
the hand-kept ones are committed with `git add -f`, as the directory is otherwise
gitignored. `python -m common.generate_inputs --include-template` adds the
template's scenarios to the generated table, so one run builds both, on the
generated base.

Like the numbers in the assumptions, the scenarios are placeholders that exercise
the machinery, not a study.

### From the input workbook

`python -m common.generate_inputs [WORKBOOK] [--include-template] [--verbose]`
reads the workbook's `assumptions_inputs` and `scenarios_clean` tabs and writes:

- the base assumptions are `config/.template.yaml` with every value the workbook
  holds replaced by the workbook's, converted from USD to EUR;
- every row of `scenarios_clean` becomes a scenario named by its `scenario_id`,
  run once per market zone where the row puts such an area on the grid;
- a scenario's overlay carries what differs from the base: the origin's and the
  destination's WACC from the workbook's EU and non-EU rows, and the sensitivity
  overrides the row names (`"finance.default_wacc": "+30%"` scales both WACCs).
  A scenario that differs in nothing gets no overlay.

It stops before writing anything if a row cannot be translated, listing each
problem with the workbook rows that have it. How the workbook names keys and
scenario vocabulary, as against the model, is in
`workflow/common/_workbook_translation.py`; `--verbose` lists the backstop values
the template supplies and the workbook rows left unused.

### The demo

The template's five scenarios, `demo-h2-dri-eaf`, `demo-moe-eaf`, their
`-export` twins and `demo-moe-eaf-ship`, which makes its steel at home and ships
it to `ES_SYN`, sit on the shipped Victoria cutout, so they run after a fresh
clone with no CDS download. Each overlay moves one value its route reads — a
salt-cavern H2 store for the H2 routes, a lower MOE turndown for the MOE ones.

### A scenario with several runs

An area on the grid that trades through market zones is run once per zone.
`report.best_run_by` in `config/config.yaml` names the field that decides which
run stands for the scenario, defaulting to `lcos_eur_per_t`, the levelised cost
of steel. The lowest wins, and a run without a value is never picked over one
with it.

`pick_best_run` copies the winner's report, diagnostic, network and assumptions
up into `results/{scenario}/`, so a scenario's folder looks the same whether it
had one run or several, and nothing about zones is left for a reader of
`report_{scenario}.csv` to interpret. The other zones are not lost, because what
they cost is itself a result: every run keeps its complete folder under
`results/{scenario}/.runs/{scenario}_{run}/`, beside `runs.csv`, which lists the
runs and says which one was picked. A scenario with one run has one entry, so
`results/*/.runs/` is the complete record of every run solved.

## Configuration

| File | Holds |
|------|-------|
| `config/config.yaml` | Pipeline knobs: `logging`, `inputs` (the base assumptions and scenario table to run), `entsoe` (data types), `nem` (`eur_per_aud` FX), `report` (`best_run_by`), `res_cf` (turbines, CF flags, cutout settings), `areas` (the area registry). |
| `config/demo.yaml` | Runs the template instead: `snakemake --configfile config/demo.yaml`. |
| `config/.template.yaml` | A complete assumptions file, every key the model reads: the demo runs on it as it stands, and `python -m common.generate_inputs [WORKBOOK]` writes `config/generated/assumptions.yaml` from it, replacing each value the input workbook, `data/assumptions/Assumptions_yaml_inputs.xlsx`, holds with the workbook's (converted from USD to EUR); the rest are backstop values, copied over as they stand. |
| `config/generated/assumptions.yaml` | Generated from the input workbook and `config/.template.yaml`; edit those instead. Base techno-economics: CAPEX/OPEX, lifetimes, the WACC (`finance.origin_wacc` for everything up to the iron, `finance.destination_wacc` for an export route's EAF, steel store and the furnace's grid connection), electrolyser efficiency, plant sizing, the steel process steps (`dri-h2`, `dri-ng`, `eaf`, `moe`, `ew`, `briquetting`, `iron_store`), natural-gas price/CO2 (`natural_gas`), grid connection charges, and what freight costs (`transport`), with the legs from each origin country to each destination in `transport.distance_km` — a run looks up its own, a market zone shipping as its country, and a pair not listed gets `default`, 3141.593 km (π × 1000), a number nobody would quote, so such a run stands out. Also the emission factors the report levelises (`emissions`) — accounting only, never priced into a solve. Numbers only — the route is chosen by the scenario table, never here. Loaded by `solve_network` as an **input file**, not a global `configfile:`. Tech keys (`res.wind-onshore`, `res.solar`, …) match the scenario table's `tech`. |
| `config/overlays/{scenario}.yaml` | *Optional* per-scenario overlay, deep-merged onto the base so it carries only the keys it moves. **File presence is the toggle**: the `optional()` shim resolves it at job-evaluation time, and a scenario without one runs on the base as it stands. It never selects a route. The generator writes and removes its own (marked `# GENERATED`); the rest are kept by hand. |
| `config/generated/scenarios.csv`, `config/.template.csv` | The scenario tables: one row per input series of a run. Columns: `scenario, run, route, tech, variant, area, start_date, end_date`. A run's rows name the origin's series first; its first area is its origin and its second its destination, and a third is an error. Route and dates are the same on every row of a run. `#` lines are comments. |

## Data formats

**Grid** (`resources/timeseries/{area}_grid_{variant}_{start}_{end}.parquet`):
UTC hourly `DatetimeIndex`, leading `price` column (EUR/MWh). Three variants:

| variant | columns | for |
|---|---|---|
| `dayahead` | `price` | a solve, and nothing else |
| `emissions` | `price` + one per carrier, hourly | the default for grid rows — the solve buys at the price, the report weights the mix by the grid emission factors under `emissions` in the assumptions |
| `full` | all data types, native resolution | analysis; adds load, RES forecast and cross-border flows |

The carrier column names are the shared vocabulary the downloaders emit
(`brown_coal`, `hard_coal`, `gas`, `wind_onshore`, …), which is what lets one
factor table serve every market. `full` is not a mix source despite carrying the
carriers: it is at native resolution, so reading an hour off it would take the
`:00` instant for the hour's mean.

Five real sources and the demo's synthetic one sit behind `retrieve_grid_data`,
and they do not all serve all three variants. A scenario's grid row names one
variant for every area it covers, so a
source that has no generation to publish answers an `emissions` request with the
price alone and says so in its log, rather than refusing and taking the price down
with it:

| source | areas | serves | generation mix? |
|---|---|---|---|
| `entsoe` | DEU, ESP, FRA | all three | yes, per carrier |
| `nem` | the five NEM regions | all three | yes, per carrier |
| `synthetic` | ES_SYN | `dayahead`, `emissions` | yes, made up — the demo's destination |
| `ons` | the four SIN submarkets | `dayahead`, `full` | **no** — the balance dataset aggregates everything thermal into one column, and no single emission factor describes gas, coal, oil, biomass and nuclear together |
| `aeso` | AB | `dayahead` | **no** — per-carrier generation needs AESO's keyed API or its bulk metered-volume files |
| `ieso` | ONT | `dayahead` | **no** — Ontario needs a separate XML report |

A run that imports from one of the last three therefore has no emission intensity.
`compile_report` reports its emission fields blank beside
`emissions_unavailable_reason`, which says which area and why — blank on its own
would mean "undefined for this run", which is the schema's existing meaning and the
wrong one here. Everything that does not need a mix is still reported, including
the MWh each user drew. The same goes for a carrier the factor table has no entry
for, except that one is an error: it means the vocabulary has drifted, not that
the data is missing.

**ONS price caveat.** ONS publishes the Marginal Operating Cost (CMO), the shadow
price from the DESSEM dispatch model. Brazil settles on CCEE's PLD, which is the
hourly CMO clipped to an annual floor and cap set by ANEEL, and `_ons.py`
reconstructs PLD by applying those limits (`ons.pld_limits`). The clip is not
cosmetic — a fifth of 2025 SE hours and nearly half of NE hours fall below the
floor — and the reconstruction has been spot-checked against a single published
weekly PLD and is otherwise unvalidated. Raw CMO is kept as `cmo_brl` in the
`full` variant so nothing is lost to the clip.

**Capacity factors** (`resources/timeseries/{area}_{tech}_area-average_{start}_{end}.parquet`):
hourly parquet, `DatetimeIndex` named `time`, one column whose name *is* the tech
key (`solar`, `wind-onshore`, …) with values in [0, 1].

**Solar tilt-mix** (`resources/timeseries/{area}_solar_tilt-mix-n{N}_{start}_{end}.parquet`):
same index, **multiple columns** — one per orientation in the sweep (`solar_az0`,
`solar_az30`, …). `solve_network` concatenates columns from all CF inputs into one
multi-tech frame.

**Results**: each run's folder, `results/{scenario}/.runs/{scenario}_{run}/`,
holds the full solved PyPSA network (`outputs/{scenario}_{run}.nc`), its report
and diagnostic (from `compile_report`) and the assumptions it ran on. The report
carries the levelised cost of steel (€/t) and optimal capacities. `results/{scenario}/` holds the same files of the run that
stands for the scenario (`outputs/{scenario}.nc`, `report_{scenario}.csv`, …).

### The shape of a report

**One column per run, one row per field.** There are far more fields than runs,
so a report reads down the page, and reports put side by side line up. A report
holds one run, named `{scenario}_{run}` like its folder.

```
field,demo-moe-eaf-export_DEMO-ES_SYN
scenario,demo-moe-eaf-export
origin,DEMO
destination,ES_SYN
route,moe-eaf-export
...
lcos_eur_per_t,914.11
```

The **identity fields lead**: the run (`scenario`, `origin` — where the iron is
made — `destination` — where the steel is made, the origin again unless the
route exports — `route, start_date, end_date`; the column's name says which run
it is, and the hidden diagnostic beside the report also carries `run`), one `{tech}_variant` row per input series
(`solar_variant: bestsite-p95`, `wind_onshore_variant: bestsite-p95`), and the
`inputs_hash`. Everything below them is numbers, so a
reader can skip the block on top and take the rest as a numeric table. The
variant matters because a best-site P95 profile and an area average answer
different questions, and the solved network does not carry it (the solve reads
the parquet, not its name), so `compile_report` takes the scenario table as an
input and reads it back off the run's rows.

`grid_variant_origin` says whether the plant at the origin was connected to the
grid, and which version of that market's data priced its power. Blank means the
plant was islanded and ran on its own renewables and storage. `emissions`, the
default, means it bought power at the hourly market price and the report works
out that power's CO₂ from the market's hourly generation mix. `dayahead` means
it bought power at the hourly price only, so the run's emission fields are blank
and `emissions_unavailable_reason` says why; a market that publishes no mix
(Brazil, Alberta, Ontario) ends up the same way even when asked for
`emissions`. Only the origin's connection is reported: an export route's
furnace also buys power at its destination, but that series is always
`emissions` and the `destination` field already names the market.

`lcoe_eur_per_mwh` divides the run's whole electricity-system cost per year by
all the electricity that entered the system that year, not by the electricity
the processes used. The denominator contains a) what the solar and wind farms
actually delivered, b) what was bought from the grid, and c) on an export
route, the power the steel furnace buys at its destination. It does not contain
curtailed power, which the farms could have generated but were switched off for,
so it never entered the system; its cost is still in the numerator, because the
farms were paid for. It does count the power lost in the battery's charge and
discharge and on the internal lines, so it is larger than what the processes
used; that difference is reported as `el_losses_gwh`. The `lcoe_*_eur_per_mwh`
parts below it share the same denominator, so they add up to
`lcoe_eur_per_mwh`.

`lcoe_storage_eur_per_mwh_discharged` is the battery's own cost per MWh it
delivered, the companion to its system share `lcoe_storage_eur_per_mwh`. It
divides a) the battery's yearly cost, energy capacity and inverter together, by
b) the electricity it gave back to the plant in a year after the discharge
loss, which the report also gives as `battery_discharged_gwh`. It does not
contain c) the cost of the power that charged it, because that power is already
paid for in the renewables and grid rows, so the figure is what the battery
adds per MWh it delivers, not the full cost of a MWh that came out of it. It is
blank when no battery was built or it never discharged.

`el_reduction_mwh_per_t_steel` is the electricity that turned ore into iron
over the year, per tonne of steel produced, and `cost_electricity_reduction_eur_per_t_steel`
is the same electricity priced at the system's average cost per MWh. Both
contain a) running the reduction shaft, whether it reduces with hydrogen, gas
or a blend, b) preheating the hydrogen fed into a shaft, and c) the
electrolysis cells on the MOE and EW routes. They do not contain the
electricity that made the hydrogen (`el_electrolyser_mwh_per_t_steel`,
`cost_electricity_electrolyser_eur_per_t_steel`), melting in the furnace, briquetting, or the
battery's and lines' losses. Being per tonne of steel rather than per tonne of
iron, a DRI route's figure carries the furnace's iron yield: about 1.17 t of
iron goes into each tonne of steel (`eaf_iron_t_per_t_steel`).

`cost_electricity_eur_per_t_iron` is the electricity cost of making the iron,
per tonne of iron. Its numerator is ironmaking's share of the origin's whole
electricity system cost (renewables, battery, grid and transmission). That
share is the electricity drawn by a) the direct reduction shaft furnace or the
MOE or EW cells, b) the hydrogen preheat and the electrolyser, and c) the
briquetting press, divided by all the electricity drawn at the origin. The
losses in the battery and the lines are included: ironmaking carries the same
share of them as of the electricity drawn. Only the electric arc furnace is
left out. The shared system is split by the average-cost rule: every MWh,
whichever step drew it and whenever, is priced at the system's average cost
per MWh, because the model does not track which source supplied which step. On
an export route the furnace melts abroad and draws nothing at the origin, so
the share is the whole origin system. The denominator is the tonnes of iron
the ironmaking step produced over the year: what leaves the shaft or the cells,
or the briquettes on a route that briquettes its iron. The rows below it cut
the same total as the steel's electricity is cut, renewables by technology
into capex and fixed O&M, the battery into power and energy capex, the grid
into connection capex, capacity fee, market purchases and volumetric fee, and
transmission, with subtotals beside their parts; they all end in
`_eur_per_t_iron`. `cost_electricity_battery_fom_eur_per_t_iron` is blank until
the battery's fixed O&M reaches the model. Converted to per tonne of steel by
the iron per tonne of steel, and added to the furnace's own electricity and its
share of the losses, it gives back `cost_electricity_eur_per_t_steel` exactly.

`h2_produced_kg_per_t_steel` is the hydrogen the electrolyser produced over the
year, per tonne of steel produced (`h2_produced_kt` over `steel_produced_mt`).
It is a model result, but on `h2-dri-eaf` it always comes back to the inputs,
the shaft's hydrogen per tonne of iron (`h2_intensity_kg_per_t_dri`) times the
iron per tonne of steel, because the hydrogen buffer is lossless and ends the
year where it began. On `mix-dri-eaf` it is the model's own choice of how much
of the reductant is hydrogen, and on the gas, MOE and EW routes it is 0.

`cost_h2_production_meur` is the yearly cost of producing the hydrogen, in M€:
a) the electrolyser itself (`cost_electrolyser_meur`), b) the hydrogen buffer
(`cost_h2_buffer_meur`), and c) the electricity the electrolyser used
(`cost_electricity_electrolyser_meur`), priced at the plant's average
electricity cost. It equals `lcoh_eur_per_kg` times `h2_produced_kt`, and each
of the three equals its `lcoh_*_eur_per_kg` row times `h2_produced_kt`. The
electricity is not extra money: it is the electrolyser's slice of the
renewables, battery, grid and transmission costs, which the plant already pays
in `cost_res_meur`, `cost_battery_meur`, `cost_grid_meur` and
`cost_transmission_meur`. So `cost_electricity_electrolyser_meur` and
`cost_h2_production_meur` add only to each other and to the electrolyser and
buffer rows. Adding either to the electricity cost groups or to
`total_annual_cost_meur` counts the electrolyser's power twice.
`cost_h2_production_eur_per_t_steel` is the same total per tonne of steel, made
of `cost_electrolyser_eur_per_t_steel`, `cost_h2_buffer_eur_per_t_steel` and
`cost_electricity_electrolyser_eur_per_t_steel`, and the same warning holds:
its electricity part is already inside `cost_electricity_eur_per_t_steel`, so
it is never added to the six group totals or to `lcos_eur_per_t`.
`cost_h2_production_eur_per_t_iron` is the same total per tonne of iron, all of
it, since all the hydrogen goes into ironmaking: `cost_electrolyser_eur_per_t_iron`
(itself capex, fixed O&M and water), `cost_h2_buffer_eur_per_t_iron` and
`cost_electricity_electrolyser_eur_per_t_iron`. Its electricity part is already
inside `cost_electricity_eur_per_t_iron`, so it is never added to that.

**Field names use `_` throughout**, including the part that names a tech or a
link. The network and the scenario tables hyphenate those ids, so the report
puts every one through `field_stem`: the tech `wind-onshore` reports as
`wind_onshore_gw_opt` and `cf_wind_onshore`, the link `dri-h2` as
`plant_dri_h2_eur_per_t_steel` and `dri_h2_utilization`. Only the *values* keep the
hyphens (`route: moe-eaf`, `solar_variant: bestsite-p95`) — those are ids, and
they still name files and config keys.

**Every run writes the same fields** whatever route it took. The set and its
order are declared in `workflow/common/_report_schema.py`, and `compile_report`
puts each run's row on it before writing. Otherwise the shape would be the
union of whatever runs a scenario happened to contain: solving only `ng-dri-eaf`
would leave the electrolyser out, adding one `moe-eaf` run would grow every other
run's row, and a consumer asking for a field nobody wrote would get silence.

A cell a run has no value for says which of two things it means:

- **`0`** — an amount, and this run's is nothing. A route without a MOE cell
  reports `plant_moe_eur_per_t_steel: 0`; it built none, and zero is what that step
  contributed to the cost of its steel. The cost groups stack to
  `total_annual_cost_meur` for every run.
- **blank** — undefined here: a ratio with nothing in its denominator. The
  capacity factor of a turbine that was not built, `lcoh_*` on a route that
  makes no hydrogen, the utilisation of a plant that is not there.

So a blank never means "not looked at", and `0` never means "missing". A
multi-site run additionally names a generator per candidate site; those fields
are not declared in advance and follow the fixed ones.

**Reading one back**: `common._report_schema.read_report(path)` returns a frame
of one row per run with the numbers as numbers. It is the only thing that undoes
the layout — every dashboard, taxonomy and plot goes through it — so how the file
is written stays that module's business alone.

### Tying a result to its inputs

A run's name identifies its inputs but does not pin them: assumptions get
edited and timeseries get rebuilt under the same name, so the same run can
describe two different solves. Each row therefore ends with an **`inputs_hash`** — a digest
over the *contents* of every file the solve read (assumptions base + overlay,
area geometry, each CF and grid parquet), stamped into the network by
`solve_network` at solve time.

Two rows with the same hash were built from byte-identical inputs; two that
differ were not, whatever their keys say. Which file moved is answered by the
per-file map that travels with the network, one content digest per input:

```python
n = pypsa.Network(); n.import_from_netcdf("results/{scenario}/outputs/{scenario}.nc")
n.meta["inputs_hash"]   # '4f2a9c1b8e30' — the value in the report
n.meta["inputs"]        # {'config/generated/assumptions.yaml': 'a91c02...', ...}
```

The chain stops at the solve: the CF parquets are pinned by content, not traced
back to the cutout that produced them.

### Cutout caching

ERA5 cutouts are expensive to (re-)download, so `cutouts/{name}.nc` is **not**
marked `protected()`. Instead, `retrieve_area_cutout` keeps a persistent cache under
`cutouts/cache/`, keyed on the *actual* request parameters (module, bounding box,
`dx`/`dy`, time range) — see `workflow/common/_cutout_cache.py`. On a cache hit it
hardlinks the entry into the rule output and skips CDS; on a miss it downloads,
QC-validates, and stores the result. Keying on the real parameters means a
`mainland_bbox` / `offshore_max_distance_km` edit correctly re-downloads instead
of silently reusing a differently-bounded cutout. The legacy sibling
`cutouts/{name}_backup.nc` still works as a fallback (`mv foo.nc foo_backup.nc`
to pin one) and is promoted into the cache on use. A request no entry matches
exactly is cut out of the smallest cached cutout, of any area and the same window
and resolution, whose box covers it — or linked to it, when the boxes are the
same, as for Brazil's submarkets, which take the national box.

Every finished cutout — freshly downloaded, cached, or from a backup — passes
structural QC (`workflow/common/_cutout_qc.py`) before the rule succeeds:
`expver`/ERA5T-mix marker, exact hourly time coverage, no gaps/duplicates, and no
NaN / fully-NaN timesteps. A partial download or ERA5/ERA5T mix therefore fails
the rule loudly instead of flowing into the CF pipeline.

## Logging & live output

Every rule writes a per-(rule, wildcards) log under `logs/` via the Snakemake
`log:` directive. Scripts log through `workflow/common/_logging.py` — no stray
`print()`s. Long loops show a `tqdm` bar with ETA on a TTY and degrade to a single
start/finish line otherwise.

```
logs/
├── retrieve_grid_data/{area}_{variant}_{start}_{end}.log
├── make_area_geometry/{area}.log
├── make_offshore_geometry/{area}.log
├── make_zone_offshore_geometry/{area}.log
├── retrieve_area_cutout/{area}_{start}_{end}.log
├── area_average/{area}_{tech}_{variant}_{start}_{end}.log
├── solve_network/{scenario}_{area}_{route}_{start}_{end}.log
└── compile_report/{scenario}.log
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

### Watching CDS / ERA5 cutout downloads

ERA5 downloads (`retrieve_area_cutout`) can take **hours**: CDS allows one running job
per user, each `retrieve_area_cutout` job submits ~5 feature requests (one per ERA5
variable group), and a single request commonly spends 20-25 min actually running
plus queue time. So a full country-year cutout is ~2-2.5 h. This is normal.

`workflow/common/_cds_monitor.py` logs progress into
`logs/retrieve_area_cutout/{area}_{start}_{end}.log`. To avoid drowning the log (and
misleading anyone tailing it), it prints **only on a status change**, plus a
periodic "still working" heartbeat:

```
[+0m]  CDS: monitoring queue (status logged on change; heartbeat every 5m)
[+1m]  CDS: running=1 queued=4 done_this_run=0 failed=0
[+24m] CDS: running=1 queued=3 done_this_run=1 failed=0  (a CDS job completed)
[+30m] CDS: still working — running=1 queued=3 done_this_run=1 failed=0 (unchanged 6m; a normal ERA5 job runs ~20-25 min, so this is expected, not a stall)
```

Read `done_this_run` (completions since this download started, counts up to ~5)
for progress. Tune the poll cadence with `res_cf.cutout.cds_poll_interval_s`.

Check the account-wide queue directly at any time:

```bash
python -c "import cdsapi; from collections import Counter; \
print(dict(Counter(j['status'] for j in cdsapi.Client(quiet=True).client.get_jobs()._json_dict['jobs'])))"
# {'accepted': 4, 'running': 1, 'successful': 2}   # accepted = queued
```

**Gotchas (learned the hard way):**
- **Daytime congestion — schedule big runs overnight.** The CDS-MARS archive has a
  fixed global cap (~460 jobs running at once); during EU daytime the queue behind
  it routinely holds several thousand jobs, so a request can wait **an hour or more
  just to *start* running** even though your own per-user limit is free. The exact
  same cutout that finished in ~20 min overnight can stall for hours mid-afternoon.
  You can see this in the accepted job's metadata (`qos.status.limit`:
  `{queued: 3892, running: 460}`). Nothing is wrong — it drains as EU load drops in
  the evening. Kick off multi-cutout batches overnight (CET) for far better
  throughput.
- **Cost/fair-use scheduling (QoS) — heavy usage slows *you* down.** The CDS runs a
  Quality-of-Service scheduler that queues requests and picks what to run next from
  "the user profile, the type of request, and the expected request cost (volume of
  data, CPU usage, etc.)", with very large requests given very low priority and the
  number of simultaneous requests per user capped ([CDS docs](https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+documentation)).
  The exact fair-share formula isn't published, but *observed* behaviour matches it:
  after ~30 jobs across 6 cutouts in one day, later cutouts waited **2-4 h for their
  first run slot even as the global queue shrank** — i.e. the account gets
  deprioritised relative to lighter users once it has consumed a lot recently.
  Practical upshot: it's genuinely per-account, so spreading a big batch over more
  days (or overnight) beats hammering it in one sitting; there's no way to buy back
  priority mid-run.
- **Don't read a single snapshot as a stall.** `running=1 queued=4` holding steady
  for 20+ min is one job legitimately running, not a hang. Look at `done_this_run`
  climbing over time, or the heartbeat's "unchanged Xm" note — not one line.
- **One running job per user.** Launching more downloads doesn't parallelise; they
  queue. Leftover/abandoned jobs from an interrupted run keep occupying the slot —
  cancel them with `cdsapi.Client(quiet=True).client.get_remote(job_id).delete()`.
- **Rate limits (HTTP 429)** appear under rapid polling/many API calls; cdsapi
  retries automatically. Keep `cds_poll_interval_s` at 60s+.
- **Recent years may be incomplete.** A year still partly in ERA5T (preliminary,
  ~5-day-to-3-month latency) can download partial or `expver`-mixed; the cutout QC
  gate catches this and fails the rule rather than passing bad data downstream.

**Batch driver.** `workflow/scripts/res_cf/download_cutouts_batch.sh` downloads a
set of cutouts one at a time through the QC-gated rule, logging to
`logs/overnight_downloads.log` and QC-summarising at the end. It is **resume-safe**:
cutouts already on disk / in `cutouts/cache/` are skipped, so re-running after an
interruption or a VM crash picks up where it left off (the expensive downloads
live in the cache, not the script). Run detached and overnight:

```bash
SMK=path/to/snakemake PY=path/to/python \
  nohup bash workflow/scripts/res_cf/download_cutouts_batch.sh &
```

**For long multi-cutout runs, have an agent babysit.** A batch of country-years is
many hours of mostly-waiting punctuated by rare failures — a good fit for a
Claude Code agent that periodically tails `logs/overnight_downloads.log` +
`logs/retrieve_area_cutout/*.log`, distinguishes "still working" from a real failure
using the signals above, and reports/retries. Run it as a background task and let
it wake on a timer.

### HPC / cluster

`profiles/slurm/config.yaml.template` is a SLURM executor placeholder. Copy it to
`profiles/slurm/config.yaml`, fill in the `cluster:` line for your scheduler, and
run `snakemake --profile profiles/slurm --jobs 200`. Script-side logging is
unchanged — per-rule files still land under `logs/{rule}/`.

## Conventions & roadmap

Project conventions (logging style, Snakefile/`.smk` rules, the two script
patterns) live in `CLAUDE.md`. Planned work and known issues are tracked in the
GitHub issues.
