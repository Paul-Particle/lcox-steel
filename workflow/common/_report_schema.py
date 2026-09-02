"""The report's fields: which ones there are, in what order, and what a blank means.

A report file is one column per run and one row per field, because there are far
more fields than runs and a reader wants to compare two runs down a page rather
than across one. The identity fields lead — what the run was, in text — so the
rest of the file is a numeric block that can be read past them.

Every run writes the same fields whatever route it took, so the report can be
read as a table rather than probed field by field. Without this the shape was the
union of whatever the runs in one scenario happened to produce: solving only
`h2-only` left out the steel chain entirely, and adding a `moe-eaf` run grew
every other run's row. Two runs of the same model gave CSVs a reader could not
treat alike, and a consumer asking for a field that was never written got silence
rather than an error.

A blank cell therefore means one thing: the quantity is undefined for this run.
That is a ratio with nothing in its denominator — the capacity factor of a
turbine that was not built, the levelised cost of hydrogen on a route that makes
none. Anything that adds into a total is written as `0`, because zero is what it
contributed. So blank is never "we did not look", and `0` is never "missing".

The steel chain's per-tonne fields follow the totals rather than the ratios:
they are shares of the levelised cost of steel and have to stack up to it, so a
route without a MOE cell reads `0` there. A run that makes no steel at all
(`h2-only`) reads `0` down the whole chain, matching its `steel_produced_mt`.

`REPORT_FIELDS` maps each field to how it is filled when the run did not produce
it. Fields a run produces that are not declared here (a multi-site run names a
generator per candidate site) follow the declared ones.
"""

import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

# The fill for a field a run did not produce.
ZERO = 0.0  # an amount, and this run's is nothing
UNDEFINED = None  # a ratio or a label with no meaning for this run

# Techs that contribute an input series to a run, named as `config/scenarios.csv`
# names them. The variant each was solved with rides along in the report.
INPUT_TECHS = ("solar", "wind-onshore", "wind-offshore", "grid")

# Renewable tech -> the stem the report writes its per-tech fields under. The
# tech is hyphenated like the wildcard; the stem is not, because the fields
# read `lcoe_wind_onshore_eur_per_mwh` and `cf_wind_onshore`.
RES_STEMS = (("solar", "solar"), ("wind-onshore", "wind_onshore"),
             ("wind-offshore", "wind_offshore"))

# The steel chain's links, by the id `build_network` gives them.
PROCESS_LINKS = ("dri-h2", "dri-ng", "eaf", "moe", "ew")

# Annualised cost groups, in the order `compile_report._cost_breakdown` builds
# them. Together they are the total annual cost.
COST_GROUPS = ("res", "battery", "grid", "gas", "electrolyser", "h2_buffer",
               "process", "ore_consumables", "iron_store", "steel_store",
               "transmission")

REPORT_FIELDS = {
    # What this run is a result for.
    "scenario": UNDEFINED,
    "area": UNDEFINED,
    "country": UNDEFINED,
    "route": UNDEFINED,
    "start_date": UNDEFINED,
    "end_date": UNDEFINED,
    **{f"{tech}_variant": UNDEFINED for tech in INPUT_TECHS},
    # Diagnostic only: the report drops it, having already selected on it.
    "best_in_country": UNDEFINED,

    # Annual system cost and the groups it is made of.
    "total_annual_cost_meur": ZERO,
    **{f"cost_{group}_meur": ZERO for group in COST_GROUPS},

    # What the run produced, and what a unit of it cost.
    "steel_produced_mt": ZERO,
    "h2_produced_kt": ZERO,
    "ng_gwh_lhv": ZERO,
    "lcos_eur_per_t": UNDEFINED,
    "lcoh_eur_per_kg": UNDEFINED,
    "lco_output": UNDEFINED,
    "lco_output_unit": UNDEFINED,

    # Electricity: the system LCOE, the contributions that sum to it, then each
    # technology's own cost over its own generation.
    "lcoe_eur_per_mwh": UNDEFINED,
    "lcoe_renewables_eur_per_mwh": ZERO,
    **{f"lcoe_{stem}_eur_per_mwh": ZERO for _, stem in RES_STEMS},
    "lcoe_storage_eur_per_mwh": ZERO,
    "lcoe_grid_connection_eur_per_mwh": ZERO,
    "lcoe_grid_energy_eur_per_mwh": ZERO,
    "lcoe_transmission_eur_per_mwh": ZERO,
    "lcoe_renewables_own_eur_per_mwh": UNDEFINED,
    **{f"lcoe_{stem}_own_eur_per_mwh": UNDEFINED for _, stem in RES_STEMS},
    **{f"cf_{stem}": UNDEFINED for _, stem in RES_STEMS},
    "cf_grid_connection": UNDEFINED,
    "grid_price_eur_per_mwh": UNDEFINED,
    "grid_fee_eur_per_mwh": UNDEFINED,
    "grid_connection_eur_per_mwh_imported": UNDEFINED,

    # Hydrogen: the levelised cost, its parts, and the plant behind them.
    "lcoh_eur_per_mwh_lhv": UNDEFINED,
    "lcoh_electrolyser_eur_per_mwh_lhv": UNDEFINED,
    "lcoh_electricity_eur_per_mwh_lhv": UNDEFINED,
    "lcoh_h2_storage_eur_per_mwh_lhv": UNDEFINED,
    "electrolyser_gw": ZERO,
    "electrolyser_utilization": UNDEFINED,
    "h2_buffer_gwh": ZERO,
    "h2_buffer_hours_dri": UNDEFINED,
    "dri_h2_mw_lhv": ZERO,

    # The steel chain: what each step cost per tonne, how big it is, how hard it
    # runs, and the stores that let it run at a different rate from the plant.
    **{f"plant_{link}_eur_per_t": ZERO for link in PROCESS_LINKS},
    "ore_eur_per_t_steel": ZERO,
    "consumables_eur_per_t_steel": ZERO,
    **{f"{link}_t_per_h_opt": ZERO for link in PROCESS_LINKS},
    **{f"{link}_utilization": UNDEFINED for link in PROCESS_LINKS},
    "iron_from_h2_share": UNDEFINED,
    "iron_store_kt": ZERO,
    "iron_store_hours_steel": UNDEFINED,
    "steel_store_kt": ZERO,
    "steel_store_hours_steel": UNDEFINED,

    # Built capacity.
    **{f"{tech}_gw_opt": ZERO for tech, _ in RES_STEMS},
    "grid_import_gw_opt": ZERO,
    "gas_supply_gw_opt": ZERO,
    "battery_gw_opt": ZERO,
    "battery_mwh_opt": ZERO,
    "transmission_total_annual_cost_meur": ZERO,

    # Which inputs produced the run (see common/_provenance.py).
    "inputs_hash": UNDEFINED,
}

ZERO_FILLED = tuple(field for field, fill in REPORT_FIELDS.items() if fill == ZERO)

# The fields that say what the run was rather than what it cost. They lead the
# file, so everything below them is numbers and a reader can skip straight to it.
IDENTITY_FIELDS = ("scenario", "area", "country", "route", "start_date", "end_date",
                   *(f"{tech}_variant" for tech in INPUT_TECHS),
                   "best_in_country", "lco_output_unit", "inputs_hash")

# Fields only the diagnostic carries: the report has already acted on the flag,
# so a frame without it is not missing anything.
DIAGNOSTIC_FIELDS = ("best_in_country",)

FIELD_ORDER = tuple(IDENTITY_FIELDS) + tuple(
    field for field in REPORT_FIELDS if field not in IDENTITY_FIELDS
)


def apply_schema(frame: pd.DataFrame) -> pd.DataFrame:
    """Put a frame of runs on the declared fields, in the declared order.

    Every run then writes the same fields whatever route it took, and a field the
    run had no value for says which of the two it means: `0` where it adds into a
    total, blank where it is undefined. Fields the schema does not declare — a
    multi-site run names a generator per candidate site — keep their place after
    the declared ones.
    """
    extra = [field for field in frame.columns if field not in REPORT_FIELDS]
    if extra:
        log.info(f"report carries {len(extra)} run-specific field(s): {extra}")
    declared = [field for field in FIELD_ORDER
                if field in frame.columns or field not in DIAGNOSTIC_FIELDS]
    on_schema = frame.reindex(columns=declared + extra)
    on_schema[list(ZERO_FILLED)] = on_schema[list(ZERO_FILLED)].fillna(0.0)
    return on_schema


def write_report_file(frame: pd.DataFrame, path: Path) -> None:
    """Write one report file: a column per run, a row per field.

    Runs are named `{scenario}_{n}` in the order they were compiled. The number
    counts within the file, so a run keeps the same name however many other zones
    were solved alongside it — which of two files a column came from is answered
    by the identity rows, not by matching numbers across them.
    """
    on_schema = apply_schema(frame)
    numeric = on_schema.select_dtypes("number").columns
    on_schema[numeric] = on_schema[numeric].round(2)
    counter = on_schema.groupby("scenario").cumcount() + 1
    on_schema.index = on_schema["scenario"].astype(str) + "_" + counter.astype(str)
    path.parent.mkdir(parents=True, exist_ok=True)
    on_schema.T.to_csv(path, index_label="field")


def read_report(path) -> pd.DataFrame:
    """A report file back as one row per run, numbers as numbers.

    The one place that undoes the on-disk layout: everything that reads a report
    — dashboards, taxonomies, plots — goes through here and works in runs-as-rows
    from then on, so how the file is laid out is this module's business alone.
    """
    stored = pd.read_csv(path, index_col="field")
    runs = stored.T
    runs.index.name = "run"
    measured = [field for field in runs.columns if field not in IDENTITY_FIELDS]
    runs[measured] = runs[measured].apply(pd.to_numeric)
    # The flag is a flag: read the diagnostic back and it is boolean again.
    if "best_in_country" in runs.columns:
        runs["best_in_country"] = runs["best_in_country"] == "True"
    return runs
