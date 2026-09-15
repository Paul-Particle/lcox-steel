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

Field names separate words with `_` throughout, including the parts that name a
tech or a link. The network and `config/scenarios.csv` hyphenate those ids —
`wind-onshore`, `dri-h2` — so every field built from one goes through
`field_stem`, and the report never mixes the two spellings the way
`wind-onshore_gw_opt` next to `lcoe_wind_onshore_eur_per_mwh` once did.

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


def field_stem(element_id: str) -> str:
    """The report's spelling of a network or wildcard id: `wind-onshore` → `wind_onshore`.

    Every field name a run's own parts go into is built through here, so the
    report separates words one way and a reader never has to know which spelling
    a given field inherited.
    """
    return element_id.replace("-", "_")


# Techs that contribute an input series to a run, named as `config/scenarios.csv`
# names them. The variant each was solved with rides along in the report.
INPUT_TECHS = ("solar", "wind-onshore", "wind-offshore", "grid")

# The renewables a run can build, named as the network names them.
RES_TECHS = ("solar", "wind-onshore", "wind-offshore")

# The steel chain's links, by the id `build_network` gives them.
PROCESS_LINKS = ("dri-h2", "dri-ng", "dri-mix", "eaf", "moe", "ew", "briquetting")

# Annualised cost groups, in the order `compile_report._cost_breakdown` builds
# them. Together they are the total annual cost.
COST_GROUPS = ("res", "battery", "grid", "gas", "electrolyser", "h2_buffer",
               "process", "ore_consumables", "iron_store", "steel_store",
               "transmission", "transport", "destination_power")

# The links that draw electricity, by the id `build_network` gives them. Each
# takes it either on bus0 (electricity in, product out) or on bus2 (a by-draw
# alongside its main conversion); which one is read off the network rather than
# listed here, so a link that moves between the two needs no change up here.
# This is the list the report has columns for; `compile_report` reads the
# drawing links off the network and refuses to report a run that has one this
# list does not name, rather than leaving its draw out of the total.
ELECTRICITY_USERS = (*PROCESS_LINKS, "electrolyser", "reductant-h2")

# The legs a run can ship over, by the id `build_network` gives them. Only one
# ever exists: an export route moves its iron, every other route its steel.
# Delivery is outside the accounting boundary below, so these are not steps.
FREIGHT_LEGS = ("iron_transport", "steel_transport")

# Everything a run can emit through, in the order `compile_report`'s breakdown
# builds them: the electricity users, the one link that burns gas and no power,
# and the two losses that belong to no step in particular. Together they are the
# run's total, so the shares stack to 100 %.
EMISSION_STEPS = (*ELECTRICITY_USERS, "reductant-ng",
                  "battery_losses", "transmission_losses")

# The links that buy ore, by the id `build_network` gives them: each carries an
# `ore_eur_per_t` quote on its marginal cost. The furnace buys consumables
# instead and the briquetting press buys nothing, so neither is here.
ORE_LINKS = ("dri-h2", "dri-ng", "dri-mix", "moe", "ew")

# The links whose electricity turns ore into iron — the three shafts, the two
# cells, and the preheat on the blended shaft's hydrogen feed. What separates
# them from the rest of the drawing links is what their power was for, which is
# how ELECTRICITY_JOBS below divides the electricity bill.
REDUCTION_LINKS = ("dri-h2", "dri-ng", "dri-mix", "moe", "ew", "reductant-h2")

# The finest split of the levelised cost of steel the report carries: one leaf
# per priced thing, each €/t steel, together the whole of it. The `cost_*_meur`
# groups above are what these roll up into — every leaf belongs to exactly one
# group, and `compile_report._leaf_breakdown` checks they still stack to the
# total rather than trusting that they do.
#
# Two leaves per plant, per renewable and per electrolyser, because a component
# whose annual cost is `annuity x capex + fixed opex` is two decisions and a
# reader comparing routes wants them apart. The network carries only their sum,
# so the share is taken from the quotes that priced it.
#
# Each leaf with the parent group it rolls up into, in the stack order the
# cost-breakdown page draws them. One list, so a leaf cannot be in the report
# without a group to belong to. What a group is called and what colour it is
# drawn in stays viz/cost_taxonomy.GROUPS' business; which group a leaf is in is
# structural and lives here.
LEAF_COSTS_BY_GROUP = (
    ("ore", "feedstock"),
    ("consumables", "feedstock"),
    *((f"{field_stem(link)}_{half}", "process")
      for link in PROCESS_LINKS for half in ("capex", "fom")),
    ("gas_fuel", "gas"),
    ("gas_carbon", "gas"),
    ("electrolyser_capex", "hydrogen"),
    ("electrolyser_fom", "hydrogen"),
    ("electrolyser_water", "hydrogen"),
    ("h2_buffer", "hydrogen"),
    *((f"res_{field_stem(tech)}_{half}", "electricity")
      for tech in RES_TECHS for half in ("capex", "fom")),
    # A geography may build a generator none of the three techs names; its cost
    # lands here undivided rather than going missing from the stack.
    ("res_other", "electricity"),
    ("battery_power", "electricity"),
    ("battery_energy", "electricity"),
    ("grid_connection_capex", "electricity"),
    ("grid_capacity_fee", "electricity"),
    ("grid_market", "electricity"),
    ("grid_fee", "electricity"),
    ("transmission", "electricity"),
    ("destination_power", "electricity"),
    ("iron_store", "storage"),
    ("steel_store", "storage"),
    ("transport", "storage"),
)
LEAF_COSTS = tuple(leaf for leaf, _ in LEAF_COSTS_BY_GROUP)
LEAF_GROUP = dict(LEAF_COSTS_BY_GROUP)
# The parent groups, in the order their leaves first appear.
LEAF_PARENTS = tuple(dict.fromkeys(group for _, group in LEAF_COSTS_BY_GROUP))

# The jobs the electricity did, which is the one cut of the levelised cost the
# leaves above cannot make. They price what was bought; these divide what one of
# those purchases was *for*: making the hydrogen, making the iron, melting it,
# and the handling and losses around them.
#
# Every megawatt-hour is priced at what the system that supplied it cost, and the
# system's whole cost is the megawatt-hours it delivered, so the four close on
# `cost_electricity_eur_per_t` exactly. They are the same money as the
# `electricity` leaves, cut a second way, which is why they are not `cost_*` and
# must never be added to one: doing so counts the electricity bill twice.
ELECTRICITY_JOBS = ("reduction", "melt", "hydrogen", "handling_losses")

REPORT_FIELDS = {
    # What this run is a result for.
    "scenario": UNDEFINED,
    "area": UNDEFINED,
    "country": UNDEFINED,
    "route": UNDEFINED,
    "start_date": UNDEFINED,
    "end_date": UNDEFINED,
    **{f"{field_stem(tech)}_variant": UNDEFINED for tech in INPUT_TECHS},
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
    **{f"lcoe_{field_stem(tech)}_eur_per_mwh": ZERO for tech in RES_TECHS},
    "lcoe_storage_eur_per_mwh": ZERO,
    "lcoe_grid_connection_eur_per_mwh": ZERO,
    "lcoe_grid_energy_eur_per_mwh": ZERO,
    "lcoe_transmission_eur_per_mwh": ZERO,
    "lcoe_destination_power_eur_per_mwh": ZERO,
    "lcoe_renewables_own_eur_per_mwh": UNDEFINED,
    **{f"lcoe_{field_stem(tech)}_own_eur_per_mwh": UNDEFINED for tech in RES_TECHS},
    **{f"cf_{field_stem(tech)}": UNDEFINED for tech in RES_TECHS},
    "cf_grid_connection": UNDEFINED,
    "grid_price_eur_per_mwh": UNDEFINED,
    "grid_fee_eur_per_mwh": UNDEFINED,
    "grid_connection_eur_per_mwh_imported": UNDEFINED,
    # Blank on every route that melts its iron where it made it: there is no
    # furnace in another market, so there is no price it paid there.
    "destination_price_eur_per_mwh": UNDEFINED,

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
    **{f"plant_{field_stem(link)}_eur_per_t": ZERO for link in PROCESS_LINKS},
    "ore_eur_per_t_steel": ZERO,
    "consumables_eur_per_t_steel": ZERO,
    **{f"{field_stem(link)}_t_per_h_opt": ZERO for link in PROCESS_LINKS},
    **{f"{field_stem(link)}_utilization": UNDEFINED for link in PROCESS_LINKS},
    "iron_from_h2_share": UNDEFINED,
    "iron_store_kt": ZERO,
    "iron_store_hours_steel": UNDEFINED,
    # Zero on a route that melts its iron where it made it and stops at the
    # plant gate. Only one of the two commodities is ever shipped.
    "transport_km": ZERO,
    "iron_shipped_kt": ZERO,
    "steel_shipped_kt": ZERO,
    "steel_store_kt": ZERO,
    "steel_store_hours_steel": UNDEFINED,

    # The levelised cost of steel by what was bought, €/t steel and stacking to
    # `lcos_eur_per_t`: the leaves (LEAF_COSTS) and the parent groups they roll
    # up into. A leaf a route has none of reads 0 — it contributed nothing, the
    # same way its cost group does.
    #
    # The two `cost_*_eur_per_t` blocks share a prefix but sit at different
    # levels of the same tree, so they are not to be added together: LEAF_COSTS
    # and LEAF_PARENTS are the lists that say which is which, and
    # `cost_{parent}_eur_per_t` is the sum of its own leaves.
    **{f"cost_{leaf}_eur_per_t": ZERO for leaf in LEAF_COSTS},
    **{f"cost_{parent}_eur_per_t": ZERO for parent in LEAF_PARENTS},
    # The electricity leaf group cut a second way, by the job each euro of it
    # paid for (ELECTRICITY_JOBS). These four stack to `cost_electricity_eur_per_t`
    # and to nothing else, so they belong in no `cost_*` sum.
    **{f"el_{job}_eur_per_t": ZERO for job in ELECTRICITY_JOBS},
    # The same capital/upkeep split on the two carriers, in their own units, so
    # each pair stacks to the reported part above it that it divides.
    "lcoe_res_capex_eur_per_mwh": ZERO,
    "lcoe_res_fom_eur_per_mwh": ZERO,
    "lcoh_electrolyser_capex_eur_per_mwh_lhv": UNDEFINED,
    "lcoh_electrolyser_fom_eur_per_mwh_lhv": UNDEFINED,
    "lcoh_electrolyser_water_eur_per_mwh_lhv": UNDEFINED,

    # What a tonne of steel took, for reading a cost against the thing behind
    # it. Each is measured off the run rather than taken from the coefficient
    # that priced it, so a furnace that ran hot in cheap hours says so.
    "eaf_el_mwh_per_t_steel": ZERO,
    "reduction_el_mwh_per_t_steel": ZERO,
    "electrolyser_el_mwh_per_t_steel": ZERO,
    "h2_kg_per_t_steel": ZERO,
    "gas_mwh_per_t_steel": ZERO,
    # Blank rather than zero on a run that built no battery: it is the store's
    # hours over its inverter rating, and neither exists to divide.
    "battery_duration_hours": UNDEFINED,

    # The capital recovery factor each plant's capex was annuitised at, and the
    # ore quote that applied to it — both fixed by the scenario rather than
    # chosen by the solve, but carried per run so the leaf split above can be
    # checked against its inputs without opening the assumptions file.
    **{f"annuity_factor_{field_stem(link)}": UNDEFINED for link in PROCESS_LINKS},
    **{f"ore_quote_{field_stem(link)}_eur_per_t": UNDEFINED for link in ORE_LINKS},
    # How much iron the furnace charged per tonne of steel, which is the ore
    # quote's other half: it depends on what the route feeds it.
    "eaf_iron_t_per_t_steel": UNDEFINED,

    # Built capacity.
    **{f"{field_stem(tech)}_gw_opt": ZERO for tech in RES_TECHS},
    "grid_import_gw_opt": ZERO,
    "gas_supply_gw_opt": ZERO,
    "battery_gw_opt": ZERO,
    "battery_mwh_opt": ZERO,
    "transmission_total_annual_cost_meur": ZERO,

    # Emissions from the energy a run used. Accounting only — none of it reaches
    # the objective, so these never move a cost. And only the energy: the process
    # steps' own direct emissions (electrodes, carbon injection, pellet carbon,
    # carbonate fluxes) sit outside the model boundary, as does delivering the
    # steel, which is why a figure here is not a CBAM or an ETS number. The
    # freight is measured all the same, at the foot of this block.
    #
    # `emissions_basis` says which of the three factor bases produced them, so no
    # number here can be read on the wrong footing; a grid run's intensity always
    # comes from its own `variant: emissions` generation series, never from a
    # stand-in.
    #
    # In kg rather than t throughout: the report rounds to two decimals, and a
    # clean route's tonne of steel lands near a thousandth of a tonne of CO2e.
    "emissions_basis": UNDEFINED,
    "emissions_kt_co2e_per_year": ZERO,
    # Blank rather than zero on a run that makes no steel, the same way
    # `lcos_eur_per_t` is: the per-step fields below are shares of a tonne and
    # stack, but this is the ratio they are shares of, and a route with no
    # denominator has no such tonne. `h2-only` is in every `all-routes`
    # scenario, and a 0.00 here would read as the cleanest steel in the table.
    "emissions_kg_co2e_per_t_steel": UNDEFINED,
    "emissions_kg_co2e_per_kg_h2": UNDEFINED,
    # The two things inside the boundary that emit, as annual totals. Together
    # they are `emissions_kt_co2e_per_year`.
    "emissions_electricity_kt_co2e_per_year": ZERO,
    "emissions_gas_kt_co2e_per_year": ZERO,
    # Per step: what it added to a tonne of steel, and what share of the tonne
    # that was. The shares stack to 100 %.
    **{f"emissions_{field_stem(step)}_kg_co2e_per_t_steel": ZERO
       for step in EMISSION_STEPS},
    **{f"emissions_{field_stem(step)}_pct": UNDEFINED for step in EMISSION_STEPS},
    # What delivering the steel added, which none of the totals above counts:
    # the boundary is the plant, and the distance to a customer is a question
    # asked of a route rather than a property of it. Add these to
    # `emissions_kg_co2e_per_t_steel` for a delivered figure. Only the
    # diagnostic carries them (see DIAGNOSTIC_FIELDS), and unlike everything
    # else here they survive an unknown grid mix, being a distance and a factor.
    "emissions_freight_kt_co2e_per_year": ZERO,
    **{f"emissions_{field_stem(leg)}_kg_co2e_per_t_steel": ZERO
       for leg in FREIGHT_LEGS},
    # Electricity by who drew it, and how dirty their own hours were. The system
    # average is what a MWh cost the run on average; a user above it bought the
    # dirty hours, one below it chased the clean ones.
    "emissions_kg_co2e_per_mwh_el": UNDEFINED,
    **{f"el_{field_stem(user)}_gwh": ZERO for user in ELECTRICITY_USERS},
    # The electricity nobody drew — the battery's round trip and the lines'
    # losses. With it the draws above account for the whole of what was
    # generated, which is what lets the `el_*_eur_per_t` jobs close on
    # the electricity bill.
    "el_losses_gwh": ZERO,
    **{f"emissions_{field_stem(user)}_kg_co2e_per_mwh_el": UNDEFINED
       for user in ELECTRICITY_USERS},

    # Why the emission fields above are blank, when they are. A market that
    # publishes prices but no per-carrier generation leaves a run that imports
    # from it with no intensity to report, which is neither zero nor undefined —
    # it is unknown, and the schema had no other way to say so. Blank here means
    # the emission fields mean what they usually mean.
    "emissions_unavailable_reason": UNDEFINED,

    # Which inputs produced the run (see common/_provenance.py).
    "inputs_hash": UNDEFINED,
}

ZERO_FILLED = tuple(field for field, fill in REPORT_FIELDS.items() if fill == ZERO)

# The fields that say what the run was rather than what it cost. They lead the
# file, so everything below them is numbers and a reader can skip straight to it.
IDENTITY_FIELDS = ("scenario", "area", "country", "route", "start_date", "end_date",
                   *(f"{field_stem(tech)}_variant" for tech in INPUT_TECHS),
                   "best_in_country", "lco_output_unit",
                   # Prose, not a measurement: everything outside this tuple is
                   # coerced to a number when a report is read back.
                   "emissions_basis", "emissions_unavailable_reason", "inputs_hash")

# Fields only the diagnostic carries, for two separate reasons. The flag,
# because the report has already acted on it, so a frame without it is not
# missing anything. The freight, because it is outside the boundary the rest of
# the emission fields are inside, and a figure that no total counts reads as one
# that some total does.
DIAGNOSTIC_FIELDS = ("best_in_country",
                     "emissions_freight_kt_co2e_per_year",
                     *(f"emissions_{field_stem(leg)}_kg_co2e_per_t_steel"
                       for leg in FREIGHT_LEGS))

FIELD_ORDER = tuple(IDENTITY_FIELDS) + tuple(
    field for field in REPORT_FIELDS if field not in IDENTITY_FIELDS
)


def apply_schema(frame: pd.DataFrame, hold_back: tuple = ()) -> pd.DataFrame:
    """Put a frame of runs on the declared fields, in the declared order.

    Every run then writes the same fields whatever route it took, and a field the
    run had no value for says which of the two it means: `0` where it adds into a
    total, blank where it is undefined. Fields the schema does not declare — a
    multi-site run names a generator per candidate site — keep their place after
    the declared ones.

    `hold_back` names the fields to leave out, which is how the report is written
    without what only the diagnostic carries. A held-back field is gone whether
    or not the frame had it, so which file a field appears in is this tuple's
    answer alone and not a matter of what the runs happened to produce.
    """
    extra = [field for field in frame.columns if field not in REPORT_FIELDS]
    if extra:
        log.info(f"report carries {len(extra)} run-specific field(s): {extra}")
    declared = [field for field in FIELD_ORDER if field not in hold_back]
    on_schema = frame.reindex(columns=declared + extra)
    # A run whose emission intensity is unknown keeps its emission fields blank.
    # Zero-filling them would say the run emitted nothing, which is the one
    # reading that is certainly wrong. Every other field fills as it always did,
    # including the MWh each user drew — that is known whatever the mix was.
    unknown = on_schema["emissions_unavailable_reason"].notna()
    # Against the columns in hand rather than the whole declaration: a held-back
    # field is not there to fill.
    zero_filled = [field for field in ZERO_FILLED if field in on_schema.columns]
    emission_fields = [field for field in zero_filled if field.startswith("emissions_")]
    other_fields = [field for field in zero_filled if not field.startswith("emissions_")]
    on_schema[other_fields] = on_schema[other_fields].fillna(0.0)
    on_schema.loc[~unknown, emission_fields] = (
        on_schema.loc[~unknown, emission_fields].fillna(0.0)
    )
    return on_schema


def write_report_file(frame: pd.DataFrame, path: Path, hold_back: tuple = ()) -> None:
    """Write one report file: a column per run, a row per field.

    Runs are named `{scenario}_{n}` in the order they were compiled. The number
    counts within the file, so a run keeps the same name however many other zones
    were solved alongside it — which of two files a column came from is answered
    by the identity rows, not by matching numbers across them.
    """
    on_schema = apply_schema(frame, hold_back)
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
