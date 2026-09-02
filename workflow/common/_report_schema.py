"""The report's columns: which ones there are, in what order, and what a blank means.

Every run writes the same columns whatever route it took, so the report can be
read as a table rather than probed column by column. Without this the header was
the union of whatever the runs in one scenario happened to produce: solving only
`h2-only` left out the steel chain entirely, and adding a `moe-eaf` run grew the
header of every other row. Two runs of the same model gave CSVs a reader could
not treat alike, and a consumer asking for a column that was never written got
silence rather than an error.

A blank cell therefore means one thing: the quantity is undefined for this run.
That is a ratio with nothing in its denominator — the capacity factor of a
turbine that was not built, the levelised cost of hydrogen on a route that makes
none. Anything that adds into a total is written as `0`, because zero is what it
contributed. So blank is never "we did not look", and `0` is never "missing".

The steel chain's per-tonne columns follow the totals rather than the ratios:
they are shares of the levelised cost of steel and have to stack up to it, so a
route without a MOE cell reads `0` there. A run that makes no steel at all
(`h2-only`) reads `0` down the whole chain, matching its `steel_produced_mt`.

`REPORT_COLUMNS` maps each column to how it is filled when the run did not
produce it, in the order the CSV writes them. Columns a run produces that are
not declared here (a multi-site run names a generator per candidate site) follow
the declared ones, and `inputs_hash` closes the row.
"""

# The fill for a column a run did not produce.
ZERO = 0.0  # an amount, and this run's is nothing
UNDEFINED = None  # a ratio or a label with no meaning for this run

# Techs that contribute an input series to a run, named as `config/scenarios.csv`
# names them. The variant each was solved with rides along in the report.
INPUT_TECHS = ("solar", "wind-onshore", "wind-offshore", "grid")

# Renewable tech -> the stem the report writes its per-tech columns under. The
# tech is hyphenated like the wildcard; the stem is not, because the columns
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

REPORT_COLUMNS = {
    # What this row is a result for.
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

    # Which inputs produced the row (see common/_provenance.py). Last, because it
    # identifies the row rather than describing it.
    "inputs_hash": UNDEFINED,
}

ZERO_FILLED = tuple(col for col, fill in REPORT_COLUMNS.items() if fill == ZERO)
