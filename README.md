# lcox-steel: Levelized Cost of Green Steel

`lcox-steel` calculates levelized costs for green-steel plants powered by
captive renewables, a grid connection, or both — the levelized cost of steel
(LCOS) for full process routes, or the levelized cost of hydrogen (LCOH) for
the pure hydrogen-supply model. It couples three data pipelines and an
optimisation model into one Snakemake workflow:

- **`res_cf`** — hourly renewable capacity factors from ERA5 reanalysis (via atlite).
- **`grid`** — hourly electricity market data from ENTSO-E (Europe) and NEM (Australia).
- **`solve`** — a PyPSA investment model that sizes generation, storage, and one
  steelmaking process route to meet a constant demand at least cost. One solve
  is one route, named by the `route` column of `config/scenarios.csv`:
  - `h2-only` — electrolyser meets a flat H2 demand (the original LCOH model)
  - `h2-dri-eaf` — electrolyser → H2 → DRI shaft → sponge iron (storable) → EAF → steel
  - `ng-dri-eaf` — fossil benchmark: natural-gas DRI (flat gas price, optional CO2 price) → EAF → steel
  - `mix-dri-eaf` — H2 and NG DRI shafts side by side; the optimiser picks the split
  - `moe-eaf` — molten oxide electrolysis → liquid iron → EAF → steel
  - `ew-eaf` — low-temperature iron electrowinning → iron plates (storable) → EAF → steel

  `all-routes` in that column means every route above. The registry lives in
  `workflow/common/_routes.py`; `*-export` ids are reserved for the trade
  scenarios and are not implemented.
- **`viz`** — a per-scenario LCOH/LCOS report plus Plotly figures, one row per run.

## Architecture at a glance

Everything is a Snakemake workflow keyed off one table, `config/scenarios.csv`.
A scenario name is an **umbrella label** and means nothing to the code — any
string works. What identifies one network is the **run key**: scenario, area,
date range, route. Every part of it but the name comes from the rows filed under
that name, so **adding a run is a CSV edit, not a Snakefile edit**.

```
config/scenarios.csv ──► res_cf ─┬─► resources/timeseries/*.parquet  (hourly CF) ─┐
                     └─► grid   ─┘                        (and grid €/MWh)        │
                                                                                  ▼
                                                       solve (PyPSA optimise) ──► results/{scenario}/{area}_{route}_{dates}.nc
                                                                                  │
                                                          viz ◄───────────────────┘
                                              results/report_{scenario}.csv + plots/*.png|html
```

A scenario is a name plus its optional `config/assumptions_{scenario}.yaml`
overlay. Its rows group by `(area, start_date, end_date)`; each group draws one
or more capacity-factor series (one per RES tech) and/or a single grid price
series, and is solved once per route. `viz` compiles the report and charts, one
row per run — so a scenario spanning several areas and years produces the
cross-country comparison table directly.

## Project structure

```
lcox-steel/
├── workflow/                       # Snakemake workflow (standard layout)
│   ├── Snakefile                   # configfile + sys.path + includes + rule all
│   ├── rules/
│   │   ├── _optional_shim.smk      # local stand-in for Snakemake's optional() (not shipped yet)
│   │   ├── grid.smk                # one retrieval rule, dispatching on the area's market
│   │   ├── res_cf.smk              # atlite CF pipeline (shapes → cutout → CF series)
│   │   ├── solve.smk               # PyPSA optimisation rule
│   │   └── viz.smk                 # compile_report + plot rules
│   ├── scripts/
│   │   ├── grid/                   # ENTSO-E + NEM download/process
│   │   │   ├── retrieve_grid_data.py # rule entrypoint: dispatches on the area's market
│   │   │   ├── _entsoe.py          # ENTSO-E: raw months → processed window
│   │   │   ├── download_entsoe.py  # ENTSO-E per-month raw-cache primitives
│   │   │   ├── _nem.py             # NEM: raw months → processed window
│   │   │   ├── download_nem.py     # NEMOSIS download primitives
│   │   │   ├── _nemosis_patches.py # AEMO User-Agent / URL-encoding workarounds
│   │   │   └── _helpers.py         # month iteration, UTC-naive coercion, completeness guard
│   │   ├── res_cf/                 # atlite capacity-factor pipeline (numbered by stage)
│   │   │   ├── a_make_area_geometry.py        # onshore area geometry (GeoParquet)
│   │   │   ├── b_make_offshore_geometry.py    # EEZ-clipped offshore geometry
│   │   │   ├── c_retrieve_area_cutout.py      # ERA5 cutout (honours `_backup.nc`; see below)
│   │   │   ├── d1_area_average.py             # area-average hourly CF per tech
│   │   │   ├── d2_bestsite_p95.py             # best-site P95 + anchored RES-mix CF
│   │   │   ├── d3_anchor_colo.py              # anchor co-located CF
│   │   │   ├── d4_tilt_mix.py                 # orientation-resolved solar CF sweep
│   │   │   ├── d5_multi.py                    # per-cell candidate grid (multi-site siting)
│   │   │   ├── spot_check_cutout.py           # standalone cutout inspection (NOT in the DAG)
│   │   │   ├── _helpers.py                    # shared helpers for the CF scripts
│   │   │   ├── reference/                     # Hannah's original scripts, verbatim (for side-by-side diffs)
│   │   │   └── README.md                      # author's notes on the CF methodology (WIP)
│   │   ├── solve/                  # PyPSA investment model
│   │   │   ├── build_network.py    # network construction (pure, importable)
│   │   │   ├── solve_network.py    # rule entrypoint: load → build → solve → write
│   │   │   └── _helpers.py         # annuity factor + electrolyser sizing
│   │   └── viz/                    # reporting + Plotly figures
│   │       ├── compile_report.py   # post-solve cost accounting → per-scenario CSV
│   │       ├── plot_capacity_bars.py  # per-run capacity bar chart
│   │       ├── plot_lcos_bars.py   # per-run cost-group breakdown in €/t steel
│   │       ├── plot_cf_map.py      # spatial mean-CF heatmap with P95 site marked
│   │       ├── _run_display.py     # run labels + best-zone filter, shared by the plots
│   │       └── style.py            # FCA Plotly template + colormap
│   └── common/                     # shared, cross-pipeline Python
│       ├── _constants.py           # physical constants (e.g. H2 LHV)
│       ├── _logging.py             # configure_logging + tqdm progress wrapper
│       ├── _paths.py               # repo-relative path roots
│       └── _stubs.py               # snakemake object stub for linters/IDEs
├── config/
│   ├── config.yaml                 # pipeline knobs (logging, entsoe, nem, res_cf)
│   ├── assumptions.yaml            # base techno-economics (CAPEX, OPEX, WACC, lifetimes)
│   ├── assumptions_{scenario}.yaml  # optional per-scenario overlay (presence = on)
│   └── scenarios.csv               # one row per (run, tech) input
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

Geographic datasets must be downloaded manually and dropped in as ZIPs — the
pipeline reads them directly via geopandas, no extraction step needed.

**World EEZ v12** — https://www.marineregions.org/downloads.php (free registration).
Choose "World EEZ v12 (2023)" → Shapefile. Save as
`data/shapes/offshore_zones/eez_v12.zip`. Any v11/v12 works; it needs the
`ISO_TER1` and `POL_TYPE` columns.

**Natural Earth 1:110m Admin-0 countries** —
https://www.naturalearthdata.com/downloads/110m-cultural-vectors/. Save as
`data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip`.

**Natural Earth 1:10m land** *(only for `res_cf.eligibility_source: availabilitymatrix`
— the finer-coastline land-sea mask; the default `indicatormatrix` does not need it)* —
https://www.naturalearthdata.com/downloads/10m-physical-vectors/. Save as
`data/shapes/ne_10m_land/ne_10m_land.zip`, or fetch directly:

```bash
curl -sSL -o data/shapes/ne_10m_land/ne_10m_land.zip \
  https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip
```

```bash
mkdir -p data/shapes/offshore_zones data/shapes/ne_110m_admin_0_countries data/shapes/ne_10m_land
# then drop the ZIPs in with the names above
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

### Quick start — the demo (no credentials, no downloads)

The default target is a small, self-contained **demo** that runs after a fresh
clone with **no CDS download, no EEZ/Natural-Earth zips, and no API keys**. It
ships a pre-sliced Victoria (Australia) cutout backup plus the derived geometry
parquets, and exercises the best-site (`07`) and complementarity (`08`)
capacity-factor science through a `solve` to a `viz` report:

```bash
snakemake --profile profiles/default --cores 4        # builds the DEMO-VIC-2025 demo
```

Outputs land at `results/report_DEMO-VIC-2025.csv` and
`results/plots/capacity_bars/DEMO-VIC-2025.png`, in a few minutes on a laptop.

### The real projects

`config/scenarios.csv` drives the full DAG. Building every real scenario needs the
external data and credentials described under [Setup](#setup) — ERA5 cutouts via
CDS, the EEZ/NE zips, and (for grid-connected scenarios) ENTSO-E or NEM access.
Target them explicitly with `snakemake all`:

```bash
snakemake -n all                                              # preview the full DAG
snakemake all --profile profiles/default --cores 4            # build every project
snakemake all --profile profiles/default --cores 4 --verbose  # loud
```

Drop `--profile profiles/default` for the bare invocation.

### Targeting one output

```bash
# Grid — ENTSO-E day-ahead prices for DEU (bidding zone DE_LU). Cap concurrent API calls:
snakemake resources/timeseries/DEU_grid_dayahead_20250101_20251231.parquet --cores 4 --resources entsoe_api=4

# Grid — NEM day-ahead prices (VIC1, 2025):
snakemake resources/timeseries/VIC1_grid_dayahead_20250101_20251231.parquet --cores 4

# res_cf — wind-onshore CF for Germany, 2023:
snakemake resources/timeseries/DEU_wind-onshore_area-average_20250101_20251231.parquet --cores 4

# solve — one network (one route):
snakemake results/standard-islanded/DEU_h2-dri-eaf_20250101_20251231.nc --cores 4
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
`wind-onshore`, `wind-offshore`, or `solar`.

> [!NOTE]
> **Cutout bounds = land ∪ offshore.** `download_cutout` reads the pre-built
> `{area}_geo.parquet` and `{area}_offshore_geo.parquet`, unions them, and
> takes the bounding box padded by `res_cf.cutout.bbox_pad_deg`. Unioning the
> offshore zone in matters: it can reach `res_cf.offshore_max_distance_km`
> (~200 km) from the coast — well beyond a land-only bbox + 1° pad — so without
> it the offshore-wind cells get clipped out of the cutout. The mainland_bbox
> filter and EEZ clip are already applied upstream by `build_regions` /
> `build_offshore_regions`, so the geometry is computed once.

### Naming convention

Five identifiers — **route, tech, variant, area, dates** — have one spelling
each, set in `config/scenarios.csv`, and every layer below spells them the same:
the assumptions blocks, the Snakemake wildcards and rule names, the script
filenames, and the PyPSA component names.

Name scenarios with **dashes** between words (`VIC-2025-solar-ew`,
`DE-2023-grid-p95`), not underscores. The underscore is the field separator in
the filenames these names compose into (`assumptions_<scenario>.yaml`,
`logs/<scenario>_<route>.log`), so reserving it keeps the boundaries legible at a
glance. Route, tech and variant follow the same rule (`moe-eaf`, `wind-onshore`,
`area-average`, `tilt-mix-n7`).

The **`area`** column is the exception — it keeps underscores, because official
market-zone codes use them (`DE_LU`). The registry is the top-level `areas:`
block in `config/config.yaml`, and nothing infers anything from the shape of a
code: what an area can be used for follows from its entry.

An area is a country (`DEU`, `FRA`, `ESP`, `AUS`, `BRA`) or a sub-national market
zone (`VIC1`). Write either — a market's own code is an alias for the area it
prices, so `DE_LU`, `FR` and `ES` resolve to `DEU`, `FRA` and `ESP`. The CSV
prefers the three-letter form.

`all-areas` starts from every area that is not someone's zone. A scenario with a
`tech=grid` row then needs prices, so an area with no market of its own is
replaced by the `zones:` through which it trades — naming Australia in a grid
scenario yields its NEM regions separately, while an islanded scenario keeps the
whole country. An area with neither market nor zones (Brazil, for now) drops out
of grid scenarios and stays in islanded ones.

## Scenarios

The scenarios in `config/scenarios.csv` exist to exercise the workflow, not to
answer anything — like the numbers in `config/assumptions.yaml`, they are
placeholders. There are three:

| scenario | covers |
|---|---|
| `standard-islanded` | every area, all routes, 2025, dedicated renewables only |
| `standard-grid` | every area that has prices (Australia via its NEM zones), all routes, 2025 |
| `sensitivity-test` | VIC1 only, two routes named explicitly — a machinery test |

`sensitivity-test` is the default `snakemake` target: it sits on the shipped
Victoria cutout, so it runs after a fresh clone with no CDS download. Its single
overlay carries a cheap salt-cavern H2 store and a lower MOE turndown, and its
two routes each read one of them — which is why it names both routes rather than
using `all-routes`, since a route that reads neither would just re-solve the base
case under a sensitivity's name.

### Reporting across a country's zones

A country that reaches its market through zones is solved once per zone, so
several report rows describe the same place. Each route then picks the zone where
that route came out cheapest — Australia's `moe-eaf` is reported from whichever
NEM region made steel cheapest with `moe-eaf` — and each date range is ranked on
its own.

`report.best_zone_by` in `config/config.yaml` names the ranking column, defaulting
to `lco_output`: the levelised cost of whatever the run produces, steel or (for
`h2-only`) hydrogen, with `lco_output_unit` alongside so a row is readable on its
own. Set it to null to rank nothing.

The losers are flagged `best_in_country: False` rather than dropped — what the
other zones cost is itself a result, and a dropped row costs a re-solve to
recover. The plots show only the flagged ones.

## Configuration

| File | Holds |
|------|-------|
| `config/config.yaml` | Pipeline knobs: `logging`, `entsoe` (data types), `nem` (`eur_per_aud` FX), `res_cf` (turbines, CF flags, cutout settings), `areas` (the area registry), `demo_scenarios`. |
| `config/assumptions.yaml` | Base techno-economics: CAPEX/OPEX, lifetimes, WACC, electrolyser efficiency, plant sizing, the steel process steps (`dri-h2`, `dri-ng`, `eaf`, `moe`, `ew`, `iron_store`), natural-gas price/CO2 (`natural_gas`) and grid connection charges. Numbers only — the route is chosen by the CSV, never here. Loaded by `solve_network` as an **input file**, not a global `configfile:`. Tech keys (`res.wind-onshore`, `res.solar`, …) match the tech wildcard. |
| `config/assumptions_{scenario}.yaml` | *Optional* per-scenario overlay — a scenario is its name plus this file. It covers every run under that name. **File presence is the toggle** (no CSV column); the `optional()` shim resolves it at job-evaluation time, and the script deep-merges it onto the base so the overlay carries only the keys it bumps. It never selects a route. Every route of the scenario shares it, so an override that only one route reads (a gas price, say) is harmless to the rest. |
| `config/scenarios.csv` | Flat table, one row per `(run, tech)` input. Columns: `scenario, route, tech, variant, area, start_date, end_date`. Rows join into a network by `(scenario, area, start_date, end_date)`. `route` holds one route id, several separated by `|`, or `all-routes`; `area` holds one area or `all-areas`. `#` rows are planned scenarios and are not built. **The scenarios shipped are placeholders that exercise the machinery, not a study.** |

## Data formats

**Grid** (`resources/{entsoe,nem}/{area}_grid_dayahead_{start}_{end}.parquet`):
UTC hourly `DatetimeIndex`, single `price` column (EUR/MWh). The `_full` variant
has MultiIndex `(area, metric)` columns covering all data types at native
resolution.

**Capacity factors** (`resources/timeseries/{area}_{tech}_area-average_{start}_{end}.parquet`):
hourly parquet, `DatetimeIndex` named `time`, one column whose name *is* the tech
key (`solar`, `wind-onshore`, …) with values in [0, 1].

**Solar tilt-mix** (`resources/timeseries/{area}_solar_tilt-mix-n{N}_{start}_{end}.parquet`):
same index, **multiple columns** — one per orientation in the sweep (`solar_az0`,
`solar_az30`, …). `solve_network` concatenates columns from all CF inputs into one
multi-tech frame.

**Results**: `results/<project>/<scenario>.nc` is the full solved PyPSA network;
`results/report_<project>.csv` (from `compile_report`) carries the levelized
cost (LCOH for `h2_only` scenarios, LCOS €/t for steel routes) and optimal
capacities for every scenario in the project.

### Cutout caching

ERA5 cutouts are expensive to (re-)download, so `cutouts/{name}.nc` is **not**
marked `protected()`. Instead, `download_cutout` keeps a persistent cache under
`cutouts/cache/`, keyed on the *actual* request parameters (module, bounding box,
`dx`/`dy`, time range) — see `workflow/common/_cutout_cache.py`. On a cache hit it
hardlinks the entry into the rule output and skips CDS; on a miss it downloads,
QC-validates, and stores the result. Keying on the real parameters means a
`mainland_bbox` / `offshore_max_distance_km` edit correctly re-downloads instead
of silently reusing a differently-bounded cutout. The legacy sibling
`cutouts/{name}_backup.nc` still works as a fallback (`mv foo.nc foo_backup.nc`
to pin one) and is promoted into the cache on use. Coverage-aware reuse (slicing
a sub-request out of a larger cached cutout; partial-month fills) is a deferred
follow-up — see `TODO.md`.

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

ERA5 downloads (`download_cutout`) can take **hours**: CDS allows one running job
per user, each `download_cutout` submits ~5 feature requests (one per ERA5
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
`logs/download_cutout/*.log`, distinguishes "still working" from a real failure
using the signals above, and reports/retries. Run it as a background task and let
it wake on a timer.

### HPC / cluster

`profiles/slurm/config.yaml.template` is a SLURM executor placeholder. Copy it to
`profiles/slurm/config.yaml`, fill in the `cluster:` line for your scheduler, and
run `snakemake --profile profiles/slurm --jobs 200`. Script-side logging is
unchanged — per-rule files still land under `logs/{rule}/`.

## Conventions & roadmap

Project conventions (logging style, Snakefile/`.smk` rules, the two script
patterns) live in `CLAUDE.md`. Known WIP and planned work — the ENTSO-E zone-list
migration, coverage-aware cutout-cache reuse (the keyed cache and CDS download
monitoring already landed), and `viz/_helpers.py` cleanup — are tracked in
`TODO.md`.
