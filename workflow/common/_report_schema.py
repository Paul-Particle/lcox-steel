"""The report's fields: which ones there are, in what order, and what a blank means.

A report file is one column per run and one row per field, because there are far
more fields than runs and a reader wants to compare two runs down a page rather
than across one. The identity fields lead — what the run was, in text — so the
rest of the file is a numeric block that can be read past them.

Every run writes the same fields whatever route it took, so the report can be
read as a table rather than probed field by field. The shape is this schema's
alone and not the union of what the runs in a scenario happened to produce: a
file of `ng-dri-eaf` runs carries the electrolyser's rows too, adding a
`moe-eaf` run leaves every other row the same width, and a consumer asking for a
declared field always gets a cell rather than silence.

A blank cell therefore means one thing: the quantity is undefined for this run.
That is a ratio with nothing in its denominator — the capacity factor of a
turbine that was not built, the levelised cost of hydrogen on a route that makes
none. Anything that adds into a total is written as `0`, because zero is what it
contributed. So blank is never "we did not look", and `0` is never "missing".

The steel chain's per-tonne fields follow the totals rather than the ratios:
they are shares of the levelised cost of steel and have to stack up to it, so a
route without a MOE cell reads `0` there.

Field names separate words with `_` throughout, including the parts that name a
tech or a link. The network and the scenario tables hyphenate those ids —
`wind-onshore`, `dri-h2` — so every field built from one goes through
`field_stem`, so every field name separates words with `_` whatever spelling the
id behind it uses.

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

    A few ids say what the network does rather than what a reader wants to
    know, and are reported under a name of their own (`REPORT_NAMES`).
    """
    return REPORT_NAMES.get(element_id, element_id.replace("-", "_"))


# Ids reported under a name other than their own. `reductant-h2` is the link
# that feeds hydrogen into the blended shaft, and the only thing it does to a
# report is draw the power that preheats that hydrogen.
REPORT_NAMES = {"reductant-h2": "h2_preheat"}


# Techs that contribute an input series to a run, named as the scenario tables
# name them. The variant each was solved with rides along in the report.
INPUT_TECHS = ("solar", "wind-onshore", "wind-offshore", "grid")


def variant_field(tech: str) -> str:
    """The report field a tech's variant lands in: `wind-onshore` → `wind_onshore_variant`.

    The grid's says whose connection it is, because an export run has two and
    the report carries only the origin's (the destination's is always
    `emissions`, and the `destination` field already names its market).
    """
    return "grid_variant_origin" if tech == "grid" else f"{field_stem(tech)}_variant"

# The renewables a run can build, named as the network names them.
RES_TECHS = ("solar", "wind-onshore", "wind-offshore")

# The steel chain's links, by the id `build_network` gives them.
PROCESS_LINKS = ("dri-h2", "dri-ng", "dri-mix", "eaf", "moe", "ew", "briquetting")

# Annualised cost groups, in the order `compile_report._cost_breakdown` builds
# them. Together they are the total annual cost.
COST_GROUPS = ("res", "battery", "grid", "gas", "electrolyser", "h2_buffer",
               "process", "ore", "eaf_consumables", "iron_or_hbi_store", "steel_store",
               "transmission", "transport", "destination_electricity")

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
# drawn in stays viz/dashboard/cost_taxonomy.GROUPS' business; which group a leaf is in is
# structural and lives here.
LEAF_COSTS_BY_GROUP = (
    ("ore", "feedstock"),
    ("eaf_consumables", "feedstock"),
    # The ironmaking plant is one pair whatever the route — the shaft or cells it
    # built, the others being none — so the row reads the same on every route.
    # The furnace and the press after it keep a pair each.
    ("ironmaking_plant_capex", "process"),
    ("ironmaking_plant_fom", "process"),
    *((f"{field_stem(link)}_{half}", "process")
      for link in PROCESS_LINKS if link not in ORE_LINKS for half in ("capex", "fom")),
    ("gas_fuel", "gas"),
    ("gas_carbon", "gas"),
    ("electrolyser_capex", "hydrogen"),
    ("electrolyser_fom", "hydrogen"),
    ("electrolyser_water", "hydrogen"),
    ("h2_buffer", "hydrogen"),
    *((f"res_{field_stem(tech)}_{half}", "electricity")
      for tech in RES_TECHS for half in ("capex", "fom")),
    # Capital only: the battery carries no fixed O&M yet.
    ("battery_power_capex", "electricity"),
    ("battery_energy_capex", "electricity"),
    ("grid_connection_capex", "electricity"),
    ("grid_capacity_fee", "electricity"),
    ("grid_market", "electricity"),
    ("grid_fee", "electricity"),
    ("transmission", "electricity"),
    # The destination furnace's grid connection on an export route, split as the
    # origin's is. Their total is the `destination_electricity` subtotal.
    ("destination_grid_connection_capex", "electricity"),
    ("destination_grid_capacity_fee", "electricity"),
    ("destination_grid_market", "electricity"),
    ("destination_grid_fee", "electricity"),
    # The stockpile between ironmaking and what follows it, whatever it holds:
    # the iron on a route that stores it, the briquettes on one that presses its
    # DRI for export. The whole of the `iron_or_hbi_store` cost group.
    ("iron_or_hbi_store", "storage"),
    ("steel_store", "storage"),
    ("transport", "storage"),
)
LEAF_COSTS = tuple(leaf for leaf, _ in LEAF_COSTS_BY_GROUP)
LEAF_GROUP = dict(LEAF_COSTS_BY_GROUP)
# The parent groups, in the order their leaves first appear.
LEAF_PARENTS = tuple(dict.fromkeys(group for _, group in LEAF_COSTS_BY_GROUP))

# Totals of a few leaves inside one parent, for a reader who wants a component's
# cost without adding its halves up. Each is reported as `cost_{name}_eur_per_t_steel`
# beside the leaves, never among them: a stack that drew a subtotal and its own
# leaves would count them twice, so nothing that stacks leaves reads these.
LEAF_SUBTOTALS = {
    **{f"res_{field_stem(tech)}": (f"res_{field_stem(tech)}_capex", f"res_{field_stem(tech)}_fom")
       for tech in RES_TECHS},
    "res": tuple(f"res_{field_stem(tech)}_{half}" for tech in RES_TECHS for half in ("capex", "fom")),
    "battery": ("battery_power_capex", "battery_energy_capex"),
    # The route's ironmaking plant, capital and upkeep: whichever shaft or cell it
    # built, the others being 0, so one row reads the same on every route.
    "ironmaking_plant": ("ironmaking_plant_capex", "ironmaking_plant_fom"),
    # The electrolyser alone; with the buffer it is the `hydrogen` parent.
    "electrolyser": ("electrolyser_capex", "electrolyser_fom", "electrolyser_water"),
    # The origin's grid connection and the energy bought through it.
    "grid": ("grid_connection_capex", "grid_capacity_fee", "grid_market", "grid_fee"),
    # The same four for the destination furnace's connection on an export route.
    "destination_electricity": ("destination_grid_connection_capex", "destination_grid_capacity_fee",
                                "destination_grid_market", "destination_grid_fee"),
}
assert all(leaf in LEAF_COSTS for leaves in LEAF_SUBTOTALS.values() for leaf in leaves)
assert not set(LEAF_SUBTOTALS) & (set(LEAF_COSTS) | set(LEAF_PARENTS))

# The electricity leaves of the origin's own system: everything in the
# `electricity` group but the destination furnace's connection. Ironmaking
# always happens at the origin, so these are what the iron's electricity is cut
# from (`cost_electricity_*_eur_per_t_iron`).
ORIGIN_ELECTRICITY_LEAVES = tuple(
    leaf for leaf in LEAF_COSTS
    if LEAF_GROUP[leaf] == "electricity" and not leaf.startswith("destination_")
)
# The rows of that cut, in the order the report lists them: each subtotal ahead
# of its own leaves, as a reader goes from the whole to its parts.
# `battery_fom` is a placeholder, blank until the battery's fixed O&M reaches the
# solve (see `cost_battery_fom_eur_per_t_steel`).
IRON_ELECTRICITY_PARTS = (
    "res",
    *(part for tech in RES_TECHS
      for part in (f"res_{field_stem(tech)}", f"res_{field_stem(tech)}_capex",
                   f"res_{field_stem(tech)}_fom")),
    "battery", "battery_power_capex", "battery_energy_capex", "battery_fom",
    "grid", "grid_connection_capex", "grid_capacity_fee", "grid_market", "grid_fee",
    "transmission",
)
assert set(ORIGIN_ELECTRICITY_LEAVES) <= set(IRON_ELECTRICITY_PARTS)
assert all(part in ORIGIN_ELECTRICITY_LEAVES or part in LEAF_SUBTOTALS or part == "battery_fom"
           for part in IRON_ELECTRICITY_PARTS)

# The jobs the electricity did, which is the one cut of the levelised cost the
# leaves above cannot make. They price what was bought; these divide what one of
# those purchases was *for*: making the hydrogen, making the iron, melting it,
# moving it about, and the losses around all of that.
#
# `briquetting` is the press that cold-forms hot DRI before it ships, and only
# that: a new drawing link is not swept into it but needs a job of its own,
# and the report stops until it has one. `losses` is power nobody drew at all: the battery's
# round trip and the lines'. They are separate jobs because only one of them is
# electricity anyone used, and `el_losses_gwh` above measures the second in its
# own units.
#
# Every megawatt-hour is priced at what the system that supplied it cost, and the
# system's whole cost is the megawatt-hours it delivered, so the five close on
# `cost_electricity_eur_per_t_steel` exactly. They are the same money as the
# `electricity` leaves, cut a second way, and are reported as
# `cost_electricity_{job}_eur_per_t_steel` to say so: a breakdown of that one
# total, which must never be added to the leaves it re-cuts, or the electricity
# bill is counted twice.
# Named after the equipment that drew the power, as the `el_*` energy rows are,
# so a job's cost and its MWh read as a pair: the furnace's melt is `eaf`, the
# hydrogen's making `electrolyser`.
ELECTRICITY_JOBS = ("reduction", "eaf", "electrolyser", "briquetting", "losses")
# The jobs that make the iron, for the same cut per tonne of iron: every job but
# the furnace, the losses being ironmaking's share of them.
# Ironmaking stops at the DRI: the press after it is a step of its own, like the
# furnace, so neither is one of these.
IRON_ELECTRICITY_JOBS = tuple(job for job in ELECTRICITY_JOBS if job not in ("eaf", "briquetting"))

REPORT_FIELDS = {
    # What this run is a result for.
    "scenario": UNDEFINED,
    "run": UNDEFINED,
    # Where the iron is made.
    "origin": UNDEFINED,
    # Where the steel is made: the origin again, unless the route exports.
    "destination": UNDEFINED,
    "route": UNDEFINED,
    "start_date": UNDEFINED,
    "end_date": UNDEFINED,
    **{variant_field(tech): UNDEFINED for tech in INPUT_TECHS},

    # Annual system cost and the groups it is made of.
    "total_annual_cost_meur": ZERO,
    **{f"cost_{group}_meur": ZERO for group in COST_GROUPS},

    # What the run produced, and what a unit of it cost.
    "steel_produced_mt": ZERO,
    # What the shaft or cells made in the year: the denominator of every
    # `*_eur_per_t_iron` row, and `ironmaking_capacity_mt_opt` times
    # `ironmaking_utilization`.
    "iron_produced_mt": ZERO,
    # What the briquetting press made of it, on a route that briquettes; 0 elsewhere.
    "hbi_produced_mt": ZERO,
    "h2_produced_kt": ZERO,
    "ng_gwh_lhv": ZERO,
    "lcos_eur_per_t": UNDEFINED,
    # The levelised cost of the iron at the plant gate, per tonne of
    # `iron_produced_mt`: everything inside a line drawn round the ironmaking
    # asset — ore, the shaft or cells, the hydrogen, its electricity, the gas, and
    # the iron's own store where the iron can be stored. Not the press, freight,
    # furnace or steel side.
    "lcoi_eur_per_t_iron": UNDEFINED,
    # And of the briquettes, on a route that presses its DRI: the iron above plus
    # the press and the HBI stockpile, per tonne of `hbi_produced_mt`.
    "lcohbi_eur_per_t_hbi": UNDEFINED,
    "lcoh_eur_per_kg": UNDEFINED,

    # Electricity: the system LCOE, the contributions that sum to it, then each
    # technology's cost over its own generation (the unit after `per_mwh` names
    # whose MWh).
    "lcoe_eur_per_mwh": UNDEFINED,
    "lcoe_renewables_eur_per_mwh": ZERO,
    **{f"lcoe_{field_stem(tech)}_eur_per_mwh": ZERO for tech in RES_TECHS},
    "lcoe_storage_eur_per_mwh": ZERO,
    "lcoe_grid_connection_eur_per_mwh": ZERO,
    "lcoe_grid_energy_eur_per_mwh": ZERO,
    "lcoe_transmission_eur_per_mwh": ZERO,
    "lcoe_destination_electricity_eur_per_mwh": ZERO,
    "lcoe_renewables_eur_per_mwh_renewables": UNDEFINED,
    **{f"lcoe_{field_stem(tech)}_eur_per_mwh_{field_stem(tech)}": UNDEFINED for tech in RES_TECHS},
    # The battery's yearly cost over what it gave back. Not what a stored MWh
    # cost in full: the power that charged it is paid for in the rows above.
    "lcoe_storage_eur_per_mwh_discharged": UNDEFINED,
    **{f"cf_{field_stem(tech)}": UNDEFINED for tech in RES_TECHS},
    "cf_grid_connection": UNDEFINED,
    "grid_price_eur_per_mwh": UNDEFINED,
    "grid_fee_eur_per_mwh": UNDEFINED,
    "grid_connection_eur_per_mwh_imported": UNDEFINED,
    # Blank on every route that melts its iron where it made it: there is no
    # furnace in another market, so there is no price it paid there.
    "destination_price_eur_per_mwh": UNDEFINED,

    # Hydrogen: the levelised cost, its parts, and the plant behind them. The
    # parts are per kg, the unit hydrogen is quoted in; they stack to
    # `lcoh_eur_per_kg`, and the total is also given per MWh LHV.
    "lcoh_eur_per_mwh_lhv": UNDEFINED,
    "lcoh_electrolyser_eur_per_kg": UNDEFINED,
    "lcoh_electricity_eur_per_kg": UNDEFINED,
    "lcoh_h2_storage_eur_per_kg": UNDEFINED,
    # The hydrogen's yearly cost in M€, to sit beside `cost_electrolyser_meur` and
    # `cost_h2_buffer_meur` as one set: the electrolyser's electricity, and the
    # three together (the LCOH times the hydrogen made). Neither is a cost group.
    # The electricity is a slice of `cost_res_meur`, `cost_battery_meur`,
    # `cost_grid_meur` and `cost_transmission_meur`, already paid there, so
    # neither row may be added to the cost groups or to `total_annual_cost_meur`:
    # that counts the electrolyser's power twice. Add them only to each other.
    "cost_electricity_electrolyser_meur": ZERO,
    "cost_h2_production_meur": ZERO,
    # The same total per tonne of steel. Its electricity part per tonne of steel
    # is `cost_electricity_electrolyser_eur_per_t_steel`, and the same warning
    # holds: it is already inside `cost_electricity_eur_per_t_steel`, so this
    # total is never added to the six group totals or to `lcos_eur_per_t`.
    "cost_h2_production_eur_per_t_steel": ZERO,
    # `cost_electrolyser_meur` in its three parts: the annualised capital, the
    # fixed O&M, and the water and consumables. They add up to it, so they are
    # never added to it or to the cost groups beside it.
    "cost_electrolyser_capex_meur": ZERO,
    "cost_electrolyser_fom_meur": ZERO,
    "cost_electrolyser_water_meur": ZERO,
    "electrolyser_gw": ZERO,
    "electrolyser_utilization": UNDEFINED,
    "h2_buffer_gwh": ZERO,
    # The same size as hydrogen mass, the unit the buffer's flows are in.
    "h2_buffer_kt": ZERO,
    # How long a full buffer would feed the shaft on its own: its size over the
    # shaft's mean hydrogen draw, on either the H2 or the blend shaft.
    "h2_buffer_duration_hours": UNDEFINED,
    # The hydrogen taken out of the buffer in a year. The buffer is lossless, so
    # it is also what went in. Cycles are equivalent full cycles over its size;
    # the share is of the hydrogen the electrolyser made, so how much of the
    # supply went through the buffer rather than straight to its users.
    "h2_buffer_discharged_kt": ZERO,
    "h2_buffer_cycles_per_year": UNDEFINED,
    "h2_buffer_share_of_h2_pct": UNDEFINED,

    # The steel chain: how big each step is, how hard it runs, and the stores
    # that let it run at a different rate from the plant. What each step cost per
    # tonne is in its `cost_{step}_capex` / `_fom` leaves above.
    # The ironmaking step's nameplate: the iron it could put out in a year at full
    # utilisation, in Mt — the size of the asset built, not what it produced. One
    # row for every route, whichever shaft or cell makes its iron.
    "ironmaking_capacity_mt_opt": ZERO,
    # How hard that step ran: its mean output over its nameplate, so what it made
    # is the capacity above times this.
    "ironmaking_utilization": UNDEFINED,
    # The other two steps' sizes, in t/h of their own output.
    **{f"{field_stem(link)}_t_per_h_opt": ZERO for link in PROCESS_LINKS if link not in ORE_LINKS},
    **{f"{field_stem(link)}_utilization": UNDEFINED for link in PROCESS_LINKS if link not in ORE_LINKS},
    # Hydrogen's share of the reducing agent, two ways. By reduction: the share of
    # the iron reduced with hydrogen, each fuel counted as the hydrogen it stands
    # in for. By energy: hydrogen's share of the reducing agent's energy content
    # (LHV) fed into the shaft, where the gas weighs more, needing more energy per
    # tonne. 1 on the H2 shaft, 0 on the gas shaft, the model's choice on the blend.
    "reductant_h2_share": UNDEFINED,
    "reductant_h2_share_energy": UNDEFINED,
    # The stockpile between ironmaking and what follows it, iron or briquettes:
    # what it holds, and how many hours of steel production that much would feed
    # the furnace for, at the furnace's own iron or briquettes per tonne of steel.
    "iron_or_hbi_store_kt": ZERO,
    "iron_or_hbi_store_hours_steel": UNDEFINED,
    # What it took out in a year (lossless, so also what went in), its equivalent
    # full cycles, and the share of the iron or briquettes made that went through
    # it rather than straight on, as the hydrogen buffer's rows are.
    "iron_or_hbi_store_discharged_kt": ZERO,
    "iron_or_hbi_store_cycles_per_year": UNDEFINED,
    "iron_or_hbi_store_share_pct": UNDEFINED,
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
    # The two `cost_*_eur_per_t_steel` blocks share a prefix but sit at different
    # levels of the same tree, so they are not to be added together: LEAF_COSTS
    # and LEAF_PARENTS are the lists that say which is which, and
    # `cost_{parent}_eur_per_t_steel` is the sum of its own leaves.
    **{f"cost_{leaf}_eur_per_t_steel": ZERO for leaf in LEAF_COSTS},
    **{f"cost_{parent}_eur_per_t_steel": ZERO for parent in LEAF_PARENTS},
    # The subtotals (LEAF_SUBTOTALS): a component's leaves added up, beside them.
    **{f"cost_{name}_eur_per_t_steel": ZERO for name in LEAF_SUBTOTALS},
    # The battery's fixed O&M, held open until the model carries it: the
    # assumption exists in the input workbook but does not reach the solve, so
    # the battery is priced on capital alone. Blank, not 0, because 0 would say
    # the battery costs nothing to run. Once it is wired it becomes a leaf
    # (`battery_fom`, parent `electricity`) and joins the `battery` subtotal.
    "cost_battery_fom_eur_per_t_steel": UNDEFINED,
    # The electricity leaf group cut a second way, by the job each euro of it
    # paid for (ELECTRICITY_JOBS). These stack to `cost_electricity_eur_per_t_steel`
    # and to nothing else, so no stack of leaves or subtotals may add them in.
    **{f"cost_electricity_{job}_eur_per_t_steel": ZERO for job in ELECTRICITY_JOBS},
    # Ironmaking's share of the origin's electricity cost, per tonne of iron:
    # the total, then cut as the steel's electricity is. Ironmaking is the shaft
    # or cells, the H2 preheat and the electrolyser, up to the DRI; the share is
    # their MWh over everything drawn at the origin, losses shared in proportion;
    # every MWh is priced at the system's average cost.
    "cost_electricity_eur_per_t_iron": ZERO,
    **{f"cost_electricity_{part}_eur_per_t_iron": UNDEFINED if part == "battery_fom" else ZERO
       for part in IRON_ELECTRICITY_PARTS},
    # The same total cut a second way, by the job each MWh did, as the steel's
    # electricity is by ELECTRICITY_JOBS. Never to be added to the cut above.
    **{f"cost_electricity_{job}_eur_per_t_iron": ZERO for job in IRON_ELECTRICITY_JOBS},
    # The briquetting press's electricity per tonne of DRI it pressed. Beside the
    # iron's electricity, not in it: the press is a step after ironmaking, as
    # the furnace is, so this is never added to `cost_electricity_eur_per_t_iron`.
    "cost_electricity_briquetting_eur_per_t_iron": ZERO,
    # The hydrogen per tonne of iron, the whole of it, since every kg goes into
    # the shaft: the total, the electrolyser and its three parts, the buffer.
    # The electricity part is `cost_electricity_electrolyser_eur_per_t_iron`
    # above, already inside `cost_electricity_eur_per_t_iron`, so the total is
    # never added to that or to any other iron or steel total.
    "cost_h2_production_eur_per_t_iron": ZERO,
    "cost_electrolyser_eur_per_t_iron": ZERO,
    "cost_electrolyser_capex_eur_per_t_iron": ZERO,
    "cost_electrolyser_fom_eur_per_t_iron": ZERO,
    "cost_electrolyser_water_eur_per_t_iron": ZERO,
    "cost_h2_buffer_eur_per_t_iron": ZERO,
    # `lcoi_eur_per_t_iron` cut into what it bought. With the electrolyser and
    # buffer rows above and `cost_electricity_eur_per_t_iron` they stack to it
    # exactly; `cost_h2_production_eur_per_t_iron` is never one of the parts, as
    # it holds the electrolyser's power a second time. The plant is the route's
    # one shaft or cell, as in `ironmaking_capacity_mt_opt`; the carbon is the
    # price on the gas's CO2, kept apart from the gas itself.
    "cost_ore_eur_per_t_iron": ZERO,
    # The plant's two halves together: a subtotal beside them, never a third part.
    "cost_ironmaking_plant_eur_per_t_iron": ZERO,
    "cost_ironmaking_plant_capex_eur_per_t_iron": ZERO,
    "cost_ironmaking_plant_fom_eur_per_t_iron": ZERO,
    "cost_gas_fuel_eur_per_t_iron": ZERO,
    "cost_gas_carbon_eur_per_t_iron": ZERO,
    "cost_iron_store_eur_per_t_iron": ZERO,
    # The cost of the iron per tonne of steel: the same yearly euros as
    # `lcoi_eur_per_t_iron` over the steel, with the two parts the steel's own
    # breakdown has no row for. Its other parts per tonne of steel are the ore,
    # the ironmaking plant, gas, electrolyser and buffer rows above. The
    # electricity is ironmaking's share — its jobs and its share of the losses —
    # and is never added to the `cost_electricity_{job}` rows, which hold it
    # already; the store is the iron's own, never the briquettes'.
    "cost_ironmaking_eur_per_t_steel": ZERO,
    "cost_electricity_ironmaking_eur_per_t_steel": ZERO,
    "cost_ironmaking_store_eur_per_t_steel": ZERO,
    # The same capital/upkeep split on the two carriers, in their own units, so
    # each pair stacks to the reported part above it that it divides.
    "lcoe_res_capex_eur_per_mwh": ZERO,
    "lcoe_res_fom_eur_per_mwh": ZERO,
    "lcoh_electrolyser_capex_eur_per_kg": UNDEFINED,
    "lcoh_electrolyser_fom_eur_per_kg": UNDEFINED,
    "lcoh_electrolyser_water_eur_per_kg": UNDEFINED,

    # What a tonne of steel took, for reading a cost against the thing behind
    # it. Each is measured off the run rather than taken from the coefficient
    # that priced it, so a furnace that ran hot in cheap hours says so.
    "el_eaf_mwh_per_t_steel": ZERO,
    "el_reduction_mwh_per_t_steel": ZERO,
    "el_electrolyser_mwh_per_t_steel": ZERO,
    "h2_produced_kg_per_t_steel": ZERO,
    "gas_mwh_per_t_steel": ZERO,
    # Blank rather than zero on a run that built no battery: it is the store's
    # hours over its inverter rating, and neither exists to divide.
    "battery_duration_hours": UNDEFINED,
    # The same ratio the other way up: the fastest the built battery can charge
    # or discharge, as a share of its energy per hour.
    "max_battery_c_rate": UNDEFINED,

    # The capital recovery factor each plant's capex was annuitised at — fixed by
    # the scenario rather than chosen by the solve, but carried per run so the
    # leaf split above can be checked against it. The feedstock quote is not: it
    # is per tonne of iron out, so `cost_ore_eur_per_t_iron` always equals it,
    # and the run's own assumptions file beside the report holds it as given.
    **{f"annuity_factor_{field_stem(link)}": UNDEFINED for link in PROCESS_LINKS},
    # The furnace's charge per tonne of steel is not a row either: it is an input
    # the route picks, and `iron_produced_mt` (or `hbi_produced_mt` where the
    # route briquettes) over `steel_produced_mt` gives it from the run itself.

    # Built capacity.
    **{f"{field_stem(tech)}_gw_opt": ZERO for tech in RES_TECHS},
    # What it delivered in a year, and what it could have but was switched off
    # for. Generated is the denominator of `lcoe_{tech}_eur_per_mwh_{tech}`.
    **{f"{field_stem(tech)}_generated_gwh": ZERO for tech in RES_TECHS},
    **{f"{field_stem(tech)}_curtailed_gwh": ZERO for tech in RES_TECHS},
    "grid_import_gw_opt": ZERO,
    # What the origin's connection carried in a year; `cf_grid_connection` is
    # this over the connection running flat out.
    "grid_imported_gwh": ZERO,
    "gas_supply_gw_opt": ZERO,
    # The destination furnace's own connection, which only an export route
    # builds. Declared like the other two so it sits with them rather than after
    # every field the schema names.
    "destination_supply_gw_opt": ZERO,
    "battery_gw_opt": ZERO,
    "battery_mwh_opt": ZERO,
    # What it gave back in a year, after the discharge loss.
    "battery_discharged_gwh": ZERO,
    # Equivalent full cycles: what it gave back over its energy capacity, so two
    # half cycles count as one. Blank without a battery, as there is nothing to
    # divide by. The model sets a lifetime in years but no cycle limit, so this
    # is what says whether that lifetime is plausible.
    "battery_cycles_per_year": UNDEFINED,

    # Emissions from the energy a run used. Accounting only — none of it reaches
    # the objective, so these never move a cost. And only the energy: the process
    # steps' own direct emissions (electrodes, carbon injection, pellet carbon,
    # carbonate fluxes) sit outside the model boundary, as does delivering the
    # steel, which is why a figure here is not a CBAM or an ETS number. The
    # freight is measured all the same, at the foot of this block.
    #
    # `emissions_basis` says which factor basis produced them, so no
    # number here can be read on the wrong footing; a grid run's intensity always
    # comes from its own `variant: emissions` generation series, never from a
    # stand-in.
    #
    # In kg rather than t throughout: the report rounds to two decimals, and a
    # clean route's tonne of steel lands near a thousandth of a tonne of CO2e.
    "emissions_basis": UNDEFINED,
    "emissions_kt_co2e_per_year": ZERO,
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
    # All of them together: the electricity the users drew, losses excluded.
    "el_users_gwh": ZERO,
    **{f"el_{field_stem(user)}_gwh": ZERO for user in ELECTRICITY_USERS},
    # The electricity nobody drew — the battery's round trip and the lines'
    # losses. With it the draws above account for the whole of what was
    # generated, which is what lets the `cost_electricity_*_eur_per_t_steel` jobs close on
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

    # Which inputs produced the run, and which code (see common/_provenance.py).
    # Two runs that agree on both agree on everything the model could control.
    "inputs_hash": UNDEFINED,
    "code_hash": UNDEFINED,
}

ZERO_FILLED = tuple(field for field, fill in REPORT_FIELDS.items() if fill == ZERO)

# The fields that say what the run was rather than what it cost. They lead the
# file, so everything below them is numbers and a reader can skip straight to it.
IDENTITY_FIELDS = ("scenario", "run", "origin", "destination", "route", "start_date", "end_date",
                   *(variant_field(tech) for tech in INPUT_TECHS),
                   # Prose, not a measurement: everything outside this tuple is
                   # coerced to a number when a report is read back.
                   "emissions_basis", "emissions_unavailable_reason", "inputs_hash",
                   "code_hash")

# Fields only the diagnostic carries, for three separate reasons. The run,
# because the report's column is named after it. The freight, because it is
# outside the boundary the rest of the emission fields are inside, and a figure
# that no total counts reads as one that some total does. The battery C-rate,
# because the report already carries it as `battery_duration_hours`; the
# diagnostic states it in the terms the solve's plausibility warning uses.
DIAGNOSTIC_FIELDS = ("run",
                     "emissions_freight_kt_co2e_per_year",
                     *(f"emissions_{field_stem(leg)}_kg_co2e_per_t_steel"
                       for leg in FREIGHT_LEGS),
                     "max_battery_c_rate")

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

    A column is named `{scenario}_{run}`, the name of the run's own folder.
    """
    column_names = frame["scenario"].astype(str) + "_" + frame["run"].astype(str)
    on_schema = apply_schema(frame, hold_back)
    numeric = on_schema.select_dtypes("number").columns
    # Adding 0.0 turns a negative zero into a plain one: a flow the solver left at
    # -1e-12, or an output whose sign was flipped when it was nothing, would
    # otherwise round to a "-0.0" that reads like a value.
    on_schema[numeric] = on_schema[numeric].round(2) + 0.0
    on_schema.index = column_names
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
    runs.index.name = "column"
    measured = [field for field in runs.columns if field not in IDENTITY_FIELDS]
    runs[measured] = runs[measured].apply(pd.to_numeric)
    return runs
