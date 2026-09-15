#!/usr/bin/env python3
"""Read the finely cut levelised steel cost the report carries, for the charts.

The report cuts each of its thirteen cost groups (`cost_*_meur`) as finely as
the model allows and writes the result out as `cost_*_eur_per_t`, one column per
priced thing, stacking to `lcos_eur_per_t`. Beside them `el_*_eur_per_t` divides
one of those groups — the electricity — by the job each euro of it paid for.
This module reads them, keys them the way the chart bands are keyed, and puts
the hover lines together.

The by-purpose stack is assembled here rather than reported: nine of its ten
bands are leaves of the cost tree added up, and only the electricity ones say
anything the tree cannot. Reporting all ten meant a second family of fields that
restated the first under different names, and two taxonomies that had to be kept
agreeing with each other by hand.

It used to compute the split instead, out of the coarse groups plus the quotes
in `config/assumptions.yaml`: a plant's fixed O&M as a config-fixed fraction of
its annual cost, the battery's two halves as a ratio of capex quotes, the
renewables and the grid apportioned by their levelised contributions. Every one
of those numbers is now the cost the solve actually incurred, taken off the
solved network in `compile_report._leaf_breakdown` — which is also the only
place that has to know how `build_network` composed a `capital_cost`, and the
only place that can check the leaves still stack to the total.

What stays here is display: what each leaf and each parent group is called and
what colour it is drawn in (`GROUPS`, `PROCESS_PLANTS`, `RES_TECH_LABELS`), and
`spec()` — the config quotes behind each leaf, the same for every run of a
scenario and so published once with the payload rather than repeated in every
record. Which group a leaf belongs to is structural rather than cosmetic and
comes from the schema's `LEAF_GROUP`, the same map the report groups by.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "workflow"))

from common._report_schema import (  # noqa: E402
    ELECTRICITY_JOBS,
    LEAF_COSTS,
    LEAF_GROUP,
    LEAF_PARENTS,
    PROCESS_LINKS,
    field_stem,
)

# Parent groups, in stack order (bottom -> top), with the colour the coarse
# taxonomies already use for that role.
GROUPS = [
    ("feedstock",   "Ore & consumables",  "#E2B681"),
    ("process",     "Process plant",      "#33434D"),
    ("gas",         "Natural gas",        "#525F6A"),
    ("hydrogen",    "Hydrogen",           "#91C096"),
    ("electricity", "Electricity system", "#0A5680"),
    ("storage",     "Solid stores",       "#D75674"),
]

# The steel chain's plants, in stack order: the link id `build_network` gives
# each one — which is also the assumptions block it is quoted in — and its
# label. `field_stem` turns the id into the stem the report and the leaves use.
PROCESS_PLANTS = [
    ("dri-h2", "H2-DRI shaft"),
    ("dri-ng", "NG-DRI shaft"),
    ("dri-mix", "Blended-reductant shaft"),
    ("moe",    "MOE cell"),
    ("ew",     "Electrowinning"),
    ("briquetting", "Briquetting press"),
    ("eaf",    "EAF"),
]

# Renewable technology -> label. The tech names the assumptions block under
# `res` and the network's carrier; `field_stem` gives the report's stem.
RES_TECH_LABELS = [
    ("solar",         "Solar"),
    ("wind-onshore",  "Wind onshore"),
    ("wind-offshore", "Wind offshore"),
]

# ---- reading the reported split ------------------------------------------

def _value(row: pd.Series, field: str):
    """A report field as a float, or None where the run left it blank (undefined)."""
    if pd.isna(row[field]):
        return None
    return float(row[field])


def _present(row: pd.Series, fields: dict) -> dict:
    """The fields that are worth a band, keyed as the page keys them.

    A leaf a route has none of reads 0 in the report — it contributed nothing —
    and drawing it would put a legend entry in front of a reader for a thing the
    route does not have. Blank means the same here: undefined, nothing to draw.
    """
    out = {}
    for key, field in fields.items():
        value = _value(row, field)
        if value:
            out[key] = value
    return out


def leaf_costs(row: pd.Series) -> dict:
    """Every leaf cost for one run, €/t steel, keyed as the chart bands are.

    Straight off the report: `compile_report._leaf_breakdown` cut the cost groups
    from the solved network and checked the leaves still stack to LCOS, so there
    is nothing to derive here and no assumptions to derive it from.
    """
    return _present(row, {leaf: f"cost_{leaf}_eur_per_t" for leaf in LEAF_COSTS})


def parent_costs(row: pd.Series) -> dict:
    """What each parent group of leaves comes to for one run, €/t steel.

    Reported rather than summed from the rounded leaves in the payload, so a
    hover quoting a leaf's share of its group divides by the whole group.
    """
    return _present(row, {parent: f"cost_{parent}_eur_per_t"
                          for parent in LEAF_PARENTS})


def leaf_inputs(row: pd.Series) -> dict:
    """The quantities behind each leaf, as (leading text, formatted value) pairs.

    What the run itself did — what it built, how hard it ran it, what it drew and
    what it paid — for the hover under each band. The config quotes that these
    multiply are the same for every run of a scenario, so they travel once in
    `spec()` rather than being repeated per record.
    """
    inputs: dict[str, list] = {}

    def add(key: str, *lines) -> None:
        kept = [line for line in lines if line is not None]
        if kept:
            inputs[key] = kept

    def quantity(field: str, label: str, template: str):
        """One line, or None where the run has no such quantity."""
        value = _value(row, field)
        return None if value is None else (label, template.format(value))

    # -- feedstock. The quote that applies is the one on whichever reduction step
    # was built, and a route melting DRI in a furnace also pays for the iron it
    # loses in the melt.
    ore_lines = []
    for link, label in PROCESS_PLANTS:
        stem = field_stem(link)
        quote_field = f"ore_quote_{stem}_eur_per_t"
        if quote_field not in row.index or not _value(row, f"plant_{stem}_eur_per_t"):
            continue
        ore_lines.append((f"ore quote, {label.lower()}",
                          f"{_value(row, quote_field):,.0f} €/t output"))
    if _value(row, "plant_eaf_eur_per_t"):
        ore_lines.append(quantity("eaf_iron_t_per_t_steel", "iron per t steel",
                                  "{:.2f} t (gangue and melting loss)"))
    add("ore", *ore_lines)

    # -- process plants. The same two lines under either half of a plant's annual
    # cost, plus the factor its capital was annuitised at under the capital half.
    for link, _ in PROCESS_PLANTS:
        stem = field_stem(link)
        built = quantity(f"{stem}_t_per_h_opt", "built", "{:,.0f} t/h output")
        utilisation = quantity(f"{stem}_utilization", "utilisation", "{:.0%}")
        add(f"{stem}_capex", built, utilisation,
            quantity(f"annuity_factor_{stem}", "annuity factor", "{:.4f}"))
        add(f"{stem}_fom", built, utilisation)

    # -- gas. Both halves of the bill ride on the same burnt megawatt-hours.
    gas_burnt = quantity("gas_mwh_per_t_steel", "gas burnt",
                         "{:,.2f} MWh LHV / t steel")
    add("gas_fuel", gas_burnt)
    add("gas_carbon", gas_burnt)

    # -- hydrogen
    electrolyser_built = quantity("electrolyser_gw", "built", "{:,.3f} GW input")
    electrolyser_use = quantity("electrolyser_utilization", "utilisation", "{:.0%}")
    h2_made = quantity("h2_kg_per_t_steel", "H₂ produced", "{:,.0f} kg H₂ / t steel")
    add("electrolyser_capex", electrolyser_built, electrolyser_use, h2_made)
    add("electrolyser_fom", electrolyser_built, electrolyser_use)
    add("electrolyser_water", h2_made,
        quantity("electrolyser_el_mwh_per_t_steel", "electricity drawn",
                 "{:,.2f} MWh / t steel"))
    add("h2_buffer",
        quantity("h2_buffer_gwh", "size", "{:,.2f} GWh LHV"),
        quantity("h2_buffer_hours_dri", "cover", "{:,.0f} h of DRI demand"))

    # -- the electricity system
    for tech, _ in RES_TECH_LABELS:
        stem = field_stem(tech)
        own = quantity(f"lcoe_{stem}_own_eur_per_mwh", "own LCOE", "{:,.1f} €/MWh")
        capacity_factor = quantity(f"cf_{stem}", "capacity factor", "{:.0%}")
        add(f"res_{stem}_capex", own, capacity_factor)
        add(f"res_{stem}_fom", own, capacity_factor)

    battery_lines = (quantity("battery_mwh_opt", "energy built", "{:,.0f} MWh"),
                     quantity("battery_duration_hours", "duration chosen", "{:,.1f} h"))
    add("battery_power", *battery_lines)
    add("battery_energy", *battery_lines)

    connection_lines = (
        quantity("lcoe_grid_connection_eur_per_mwh", "levelised",
                 "{:,.1f} €/MWh delivered"),
        quantity("cf_grid_connection", "utilisation", "{:.0%}"),
    )
    add("grid_connection_capex", *connection_lines)
    add("grid_capacity_fee", *connection_lines)
    add("grid_market", quantity("grid_price_eur_per_mwh", "average price paid",
                                "{:,.1f} €/MWh imported"))
    add("transmission", quantity("lcoe_transmission_eur_per_mwh", "levelised",
                                 "{:,.1f} €/MWh delivered"))

    # -- getting the iron to a furnace that is somewhere else
    add("transport",
        quantity("transport_km", "distance", "{:,.0f} km"),
        quantity("iron_shipped_kt", "iron shipped", "{:,.0f} kt/yr"),
        quantity("steel_shipped_kt", "steel shipped", "{:,.0f} kt/yr"))

    # -- the solid stores that make turndown possible
    for group in ("iron_store", "steel_store"):
        add(group,
            quantity(f"{group}_kt", "size", "{:,.1f} kt"),
            quantity(f"{group}_hours_steel", "cover", "{:,.0f} h of demand"))
    return inputs


def purpose_bands(row: pd.Series) -> dict:
    """The cost of steel by what each euro was spent for (€/t steel), summing to LCOS.

    Built from the cost leaves and the electricity jobs, which between them hold
    every euro once: the plant's capital and its upkeep are its leaves gathered
    two ways, and the electricity bands are the jobs. The underscored keys are
    hover detail rather than bands, which is what keeps them out of the stack.
    """
    plants = [field_stem(link) for link in PROCESS_LINKS]
    electrolyser = ["electrolyser_capex", "electrolyser_fom", "electrolyser_water"]
    composed = {
        "ore": ["cost_feedstock_eur_per_t"],
        "capex": [f"cost_{plant}_capex_eur_per_t" for plant in plants],
        "fixed_om": [f"cost_{plant}_fom_eur_per_t" for plant in plants],
        "hydrogen": ([f"cost_{leaf}_eur_per_t" for leaf in electrolyser]
                     + ["cost_h2_buffer_eur_per_t", "el_hydrogen_eur_per_t"]),
        "gas": ["cost_gas_eur_per_t"],
        "transport": ["cost_transport_eur_per_t"],
        "store": ["cost_iron_store_eur_per_t", "cost_steel_store_eur_per_t"],
        "_hydrogen_electrolyser": [f"cost_{leaf}_eur_per_t" for leaf in electrolyser],
        "_hydrogen_buffer": ["cost_h2_buffer_eur_per_t"],
        "_hydrogen_electricity": ["el_hydrogen_eur_per_t"],
        "_electricity_total": ["cost_electricity_eur_per_t"],
    }
    bands = {}
    for key, fields in composed.items():
        total = sum(_value(row, field) or 0.0 for field in fields if field in row)
        if total:
            bands[key] = total
    bands.update(_present(row, {job: f"el_{job}_eur_per_t"
                                for job in ELECTRICITY_JOBS if job != "hydrogen"}))
    # What the two electricity bands that have a natural per-tonne figure took,
    # measured over the year each ran rather than from the coefficient that
    # priced it.
    for key, field in (("_melt_mwh_per_t", "eaf_el_mwh_per_t_steel"),
                       ("_reduction_mwh_per_t", "reduction_el_mwh_per_t_steel")):
        drawn = _value(row, field)
        if drawn:
            bands[key] = drawn
    return bands


def carrier_splits(row: pd.Series) -> dict:
    """The LCOE and LCOH bands cut the finer way, in the reports' own units.

    Same totals as the reported parts either way; this cut just separates
    capital from fixed O&M on the renewables, and capital from fixed O&M from
    water on the electrolyser — which the report carries as its own columns.
    """
    lcoe = _present(row, {
        "res_capex": "lcoe_res_capex_eur_per_mwh",
        "res_fom": "lcoe_res_fom_eur_per_mwh",
        "storage": "lcoe_storage_eur_per_mwh",
        "grid_connection": "lcoe_grid_connection_eur_per_mwh",
        "grid_energy": "lcoe_grid_energy_eur_per_mwh",
        "transmission": "lcoe_transmission_eur_per_mwh",
        "destination_power": "lcoe_destination_power_eur_per_mwh",
    })
    lcoh = _present(row, {
        "electrolyser_capex": "lcoh_electrolyser_capex_eur_per_mwh_lhv",
        "electrolyser_fom": "lcoh_electrolyser_fom_eur_per_mwh_lhv",
        "electrolyser_water": "lcoh_electrolyser_water_eur_per_mwh_lhv",
        "storage": "lcoh_h2_storage_eur_per_mwh_lhv",
        "electricity": "lcoh_electricity_eur_per_mwh_lhv",
    })
    return {"lcoe": lcoe, "lcoh": lcoh}

# ---- the published leaf specification ------------------------------------

def spec(assumptions: dict) -> list:
    """[(key, label, parent group, colour, [(constant, value)])] in stack order.

    The constants are the config quotes behind each leaf — the same for every
    scenario the dashboard covers, so they travel once with the payload instead of
    being repeated in each record. Which parent group a leaf is in comes from the
    schema's `LEAF_GROUP`, which is also what the report groups its columns by,
    so the two cannot put the same leaf in different places.
    """
    wacc = assumptions["finance"]["default_wacc"]
    wacc_line = ("WACC", f"{wacc * 100:.1f}%")
    entries = []

    def add(key, label, colour, constants=()):
        entries.append([key, label, LEAF_GROUP[key], colour,
                        [list(pair) for pair in constants if pair is not None]])

    # Ore is priced per tonne of the reduction step's own output, and the ore grade
    # each route tolerates differs — so the quote that applies is per scenario, and
    # travels with the record rather than here.
    add("ore", "Iron ore", "#E2B681")
    add("consumables", "EAF consumables", "#C99A5E",
        [("electrodes, fluxes, alloys, carbon",
          f"{assumptions['eaf']['consumables_eur_per_t']:,.0f} €/t steel")])

    # Process plants: two leaves each, capital then fixed O&M.
    # One entry per PROCESS_PLANTS link: capital first, then its fixed O&M a shade
    # lighter. The blended shaft sits between the two single-fuel ones, as it does
    # everywhere else.
    shades = {"dri-h2": ("#33434D", "#5A6B77"), "dri-ng": ("#3D4E59", "#687985"),
              "dri-mix": ("#374852", "#63747F"),
              "moe": ("#2B3A44", "#54656F"),
              "ew": ("#25333B", "#4C5D66"),
              "briquetting": ("#44555F", "#77888F"),
              "eaf": ("#3A4A54", "#6E7F89")}
    for link, label in PROCESS_PLANTS:
        plant = assumptions[link]
        stem = field_stem(link)
        capex_colour, om_colour = shades[link]
        add(f"{stem}_capex", f"{label} — capital", capex_colour,
            [("capex quote", f"{plant['capex_per_t_per_year_eur']:,.0f} €/(t·yr)"),
             ("lifetime", f"{plant['lifetime_years']:.0f} y"), wacc_line])
        add(f"{stem}_fom", f"{label} — fixed O&M", om_colour,
            [("fixed opex", f"{plant['opex_per_t_per_year_eur']:,.1f} €/(t·yr)"),
             ("as a share of capex",
              f"{plant['opex_per_t_per_year_eur'] / plant['capex_per_t_per_year_eur'] * 100:.1f}% / yr")])

    gas = assumptions["natural_gas"]
    add("gas_fuel", "Natural gas — fuel", "#525F6A",
        [("price", f"{gas['price_eur_per_mwh']:,.1f} €/MWh LHV")])
    # A run priced at zero carbon has no carbon band: it would be worth 0 €/t on
    # every route, so it draws nothing and only puts the words "carbon price" in
    # front of a reader of a run that does not have one.
    if gas.get("co2_price_eur_per_t", 0.0) > 0:
        add("gas_carbon", "Natural gas — CO₂ price", "#7A8792",
            [("carbon price", f"{gas['co2_price_eur_per_t']:,.0f} €/t CO₂"),
             ("emission factor", f"{gas['co2_t_per_mwh']:.2f} t CO₂ / MWh LHV")])

    el_cfg = assumptions["electrolyser"]
    add("electrolyser_capex", "Electrolyser — capital", "#6FA875",
        [("capex quote", f"{el_cfg['capex_per_mw_eur'] / 1e6:,.2f} M€/MW"),
         ("lifetime", f"{el_cfg['lifetime_years']:.0f} y"), wacc_line,
         ("efficiency", f"{el_cfg['efficiency_kwh_per_kg']:,.0f} kWh/kg H₂")])
    add("electrolyser_fom", "Electrolyser — fixed O&M", "#91C096",
        [("fixed opex", f"{el_cfg['opex_per_mw_per_year_eur'] / 1e3:,.0f} k€/MW·yr")])
    add("electrolyser_water", "Electrolyser — water & variable opex", "#B4D4B8",
        [("variable opex", f"{el_cfg['varopex_eur_per_mwh_el']:,.2f} €/MWh electricity")])
    # Quoted in €/MWh rather than k€/MWh: the salt-cavern sensitivity drops this from
    # 10,000 to 350, which rounds away entirely on the larger unit.
    buffer = assumptions["h2_buffer"]
    add("h2_buffer", "H₂ buffer store", "#70D2F0",
        [("capex quote", f"{buffer['capex_per_mwh_eur']:,.0f} €/MWh LHV"),
         ("lifetime", f"{buffer['lifetime_years']:.0f} y"), wacc_line])

    res_colours = {"solar": ("#0A5680", "#3E7CA3"), "wind_onshore": ("#0E6FA4", "#4B93BF"),
                   "wind_offshore": ("#0293D2", "#57B6E2")}
    for tech, label in RES_TECH_LABELS:
        cfg = assumptions["res"].get(tech)
        if cfg is None:
            continue
        stem = field_stem(tech)
        capex_colour, om_colour = res_colours[stem]
        add(f"res_{stem}_capex", f"{label} — capital", capex_colour,
            [("capex quote", f"{cfg['capex_per_mw_eur'] / 1e6:,.2f} M€/MW"),
             ("lifetime", f"{cfg['lifetime_years']:.0f} y"), wacc_line])
        add(f"res_{stem}_fom", f"{label} — fixed O&M", om_colour,
            [("fixed opex", f"{cfg['opex_per_mw_per_year_eur'] / 1e3:,.0f} k€/MW·yr")])
    add("res_other", "Renewables (not itemised)", "#7FA8C0")

    battery = assumptions["battery"]
    add("battery_power", "Battery — power capacity", "#D75674",
        [("capex quote", f"{battery['capex_per_mw_eur'] / 1e6:,.2f} M€/MW"),
         ("lifetime", f"{battery['lifetime_years']:.0f} y"), wacc_line])
    add("battery_energy", "Battery — energy capacity", "#E58AA0",
        [("capex quote", f"{battery['capex_per_mwh_eur'] / 1e6:,.2f} M€/MWh"),
         ("duration", "sized by the optimiser"),
         ("round-trip efficiency", f"{battery['efficiency_roundtrip'] * 100:.0f}%")])

    grid = assumptions["grid"]
    add("grid_connection_capex", "Grid connection — capital", "#71828F",
        [("capex quote", f"{grid['connection_capex_eur_per_mw'] / 1e3:,.0f} k€/MW"),
         ("lifetime", f"{grid['connection_lifetime_years']:.0f} y"), wacc_line])
    add("grid_capacity_fee", "Grid connection — capacity charge", "#98A5AE",
        [("capacity charge (Leistungspreis)",
          f"{grid['fee_eur_per_mw_per_year'] / 1e3:,.0f} k€/MW·yr")])
    add("grid_market", "Grid energy — market price", "#B7C1C8")
    add("grid_fee", "Grid energy — volumetric fee", "#D3DAE0",
        [("volumetric charge (Arbeitspreis)", f"{grid['fee_eur_per_mwh']:,.1f} €/MWh imported")])
    add("transmission", "Transmission (HVDC)", "#83D1DD")
    destination = assumptions["destination"]
    add("destination_power", "Destination power", "#0293D2",
        [("market", f"{destination['area']}, hourly day-ahead")])
    freight = assumptions["transport"]
    add("transport", "Freight", "#BDCCD9",
        [(f"{mode}, {commodity}",
          f"{freight[mode][commodity]['eur_per_t']:,.0f} €/t "
          f"+ {freight[mode][commodity]['eur_per_t_km'] * 1000:,.2f} €/t per 1000 km")
         for mode in ("sea", "rail") for commodity in ("iron", "steel")])

    add("iron_store", "Iron stockpile", "#D75674",
        [("capex quote", f"{assumptions['iron_store']['capex_per_t_eur']:,.0f} €/t"),
         ("lifetime", f"{assumptions['iron_store']['lifetime_years']:.0f} y")])
    add("steel_store", "Steel inventory", "#EE8DA3",
        [("capex quote", f"{assumptions['steel_store']['capex_per_t_eur']:,.0f} €/t"),
         ("lifetime", f"{assumptions['steel_store']['lifetime_years']:.0f} y")])
    return entries
