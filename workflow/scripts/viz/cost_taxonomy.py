#!/usr/bin/env python3
"""Cut a scenario's levelised steel cost finer than the report's own cost groups.

The pipeline reports eleven cost groups (`cost_*_meur`). Several of them bundle
things a reader wants to see apart — a plant's annualised capital and its fixed
O&M, the gas bill and the carbon price, the electrolyser's capital and the water
it consumes, renewables by technology. This module takes those apart.

Nothing here needs a new reported quantity. Every split is exact, because of how
`build_network` assembles the costs it hands to PyPSA:

  * a plant's `capital_cost` is `annuity x capex_per_t_per_year +
    opex_per_t_per_year`, so fixed O&M is a config-fixed fraction of the annual
    cost, applied to the reported per-plant `plant_*_eur_per_t`;
  * likewise for renewables and the electrolyser, per MW;
  * the battery folds energy capex into its per-MW cost at a fixed duration, so
    the power/energy split is `capex_per_mw : max_hours x capex_per_mwh`;
  * the gas bill and any carbon price both live on the gas_supply generator's
    marginal cost, and the reported gas energy separates them;
  * the electrolyser's variable opex is `varopex_eur_per_mwh_el` over the
    electricity it drew, which the reported H2 production fixes exactly;
  * renewables and grid costs are apportioned by the per-technology and
    per-component LCOE contributions the report already publishes, which are
    shares of the same annual cost.

Two levels come out of `leaf_costs`:

  * `leaves` — every leaf cost in €/t steel, summing to LCOS;
  * `inputs` — the scenario-varying quantities behind each leaf (built capacity,
    utilisation, prices, energy), for the hover. The config constants that go
    with them are the same for every scenario, so `spec()` publishes them once
    rather than repeating them per record.

`GROUPS` names each leaf's parent so a hover can quote a leaf's share of the
group it belongs to as well as its share of LCOS.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "workflow"))

from common._report_schema import field_stem  # noqa: E402
from common._runs import EAF_CHARGE  # noqa: E402
from scripts.solve._helpers_solve import annuity_factor  # noqa: E402

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


def _value(row: pd.Series, field: str):
    """A report field as a float, or None where the run left it blank (undefined)."""
    if pd.isna(row[field]):
        return None
    return float(row[field])


def _group_eur_per_t(row: pd.Series, group: str, steel_t: float) -> float:
    """One reported cost group in €/t steel; groups a scenario lacks read as zero."""
    value = _value(row, f"cost_{group}_meur")
    return 0.0 if value is None else value * 1e6 / steel_t


def _split_by_shares(total: float, shares: dict) -> dict:
    """Apportion `total` by `shares` (any positive scale), or evenly if all are zero.

    Used where the report publishes each part's contribution to a levelised metric
    rather than its annual cost: the contributions are shares of the same annual
    total, so their ratio carries over exactly.
    """
    positive = {key: value for key, value in shares.items() if value and value > 0}
    if not positive:
        return {}
    scale = sum(positive.values())
    return {key: total * value / scale for key, value in positive.items()}


# ---- fixed-O&M fractions (config constants) ------------------------------

def _fixed_om_fraction(annual_capex: float, fixed_om: float) -> float:
    return fixed_om / (annual_capex + fixed_om) if (annual_capex + fixed_om) else 0.0


def process_om_fractions(assumptions: dict) -> dict:
    """Fixed O&M's share of each process plant's annual cost."""
    wacc = assumptions["finance"]["default_wacc"]
    fractions = {}
    for link, _ in PROCESS_PLANTS:
        plant = assumptions[link]
        annual_capex = (annuity_factor(wacc, plant["lifetime_years"])
                        * plant["capex_per_t_per_year_eur"])
        fractions[link] = _fixed_om_fraction(annual_capex, plant["opex_per_t_per_year_eur"])
    return fractions


def res_om_fractions(assumptions: dict) -> dict:
    """Fixed O&M's share of each renewable technology's annual cost, per MW."""
    wacc = assumptions["finance"]["default_wacc"]
    fractions = {}
    for tech, _ in RES_TECH_LABELS:
        cfg = assumptions["res"].get(tech)
        if cfg is None:
            continue
        annual_capex = annuity_factor(wacc, cfg["lifetime_years"]) * cfg["capex_per_mw_eur"]
        fractions[tech] = _fixed_om_fraction(annual_capex, cfg["opex_per_mw_per_year_eur"])
    return fractions


def electrolyser_om_fraction(assumptions: dict) -> float:
    """Fixed O&M's share of the electrolyser's annual *capital* cost (water excluded)."""
    cfg = assumptions["electrolyser"]
    wacc = assumptions["finance"]["default_wacc"]
    annual_capex = annuity_factor(wacc, cfg["lifetime_years"]) * cfg["capex_per_mw_eur"]
    return _fixed_om_fraction(annual_capex, cfg["opex_per_mw_per_year_eur"])


def battery_power_fraction(assumptions: dict) -> float:
    """The power half of the battery's per-MW cost, which folds in `max_hours` of energy."""
    cfg = assumptions["battery"]
    energy = cfg["capex_per_mwh_eur"] * cfg["max_hours"]
    total = cfg["capex_per_mw_eur"] + energy
    return cfg["capex_per_mw_eur"] / total if total else 0.0


# ---- the leaf split ------------------------------------------------------

def leaf_costs(row: pd.Series, assumptions: dict) -> tuple[dict, dict]:
    """(leaves, inputs) for one scenario row — leaves in €/t steel, summing to LCOS.

    `inputs` maps a leaf key to the scenario-varying quantities behind it, each a
    (label, formatted value) pair. Config constants are not repeated here; they
    travel once in `spec()`.
    """
    steel_t = (_value(row, "steel_produced_mt") or 0.0) * 1e6
    if steel_t <= 0:
        return {}, {}

    leaves: dict[str, float] = {}
    inputs: dict[str, list] = {}
    wacc = assumptions["finance"]["default_wacc"]

    def add(key: str, value: float, lines: list | None = None) -> None:
        if value and value > 1e-9:
            leaves[key] = value
            if lines:
                inputs[key] = [line for line in lines if line is not None]

    # -- feedstock: ore on the reduction/electrolysis links, consumables on the EAF
    # The ore quote that applies is the one on whichever reduction step was built,
    # and a route melting DRI in an EAF also pays for the yield loss in melting.
    ore_quotes = [
        (f"ore quote, {label.lower()}", f"{assumptions[link]['ore_eur_per_t']:,.0f} €/t output")
        for link, label in PROCESS_PLANTS
        if "ore_eur_per_t" in assumptions[link] and _value(row, f"plant_{field_stem(link)}_eur_per_t")
    ]
    if _value(row, "plant_eaf_eur_per_t"):
        _, iron_source = EAF_CHARGE[row["route"]]
        ore_quotes.append(("iron per t steel",
                           f"{assumptions[iron_source]['iron_t_per_t_steel']:.2f} t "
                           "(gangue and melting loss)"))
    add("ore", _value(row, "ore_eur_per_t_steel") or 0.0, ore_quotes)
    add("consumables", _value(row, "consumables_eur_per_t_steel") or 0.0)

    # -- process plants, each cut into annualised capital and fixed O&M
    om_fractions = process_om_fractions(assumptions)
    for link, _ in PROCESS_PLANTS:
        stem = field_stem(link)
        plant_cost = _value(row, f"plant_{stem}_eur_per_t")
        if not plant_cost:
            continue
        fraction = om_fractions[link]
        capacity = _value(row, f"{stem}_t_per_h_opt")
        utilisation = _value(row, f"{stem}_utilization")
        lines = [
            ("built", f"{capacity:,.0f} t/h output" if capacity is not None else None),
            ("utilisation", f"{utilisation * 100:.0f}%" if utilisation is not None else None),
            ("annuity factor",
             f"{annuity_factor(wacc, assumptions[link]['lifetime_years']):.4f}"),
        ]
        lines = [line for line in lines if line[1] is not None]
        add(f"{stem}_capex", plant_cost * (1.0 - fraction), lines)
        add(f"{stem}_fom", plant_cost * fraction, lines[:2])

    # -- natural gas: the fuel bill and, separately, any carbon price on it
    gas_total = _group_eur_per_t(row, "gas", steel_t)
    if gas_total > 0:
        gas_cfg = assumptions["natural_gas"]
        gas_mwh = (_value(row, "ng_gwh_lhv") or 0.0) * 1e3
        price = gas_cfg["price_eur_per_mwh"]
        carbon = gas_cfg.get("co2_price_eur_per_t", 0.0) * gas_cfg["co2_t_per_mwh"]
        # Both ride on the same marginal cost, so their config ratio splits the bill.
        parts = _split_by_shares(gas_total, {"gas_fuel": price, "gas_carbon": carbon})
        # The price and emission factor are config quotes and travel with the spec,
        # which is built per scenario — the gas price is overlaid per geography.
        energy_line = ("gas burnt", f"{gas_mwh / max(steel_t, 1):,.2f} MWh LHV / t steel")
        add("gas_fuel", parts.get("gas_fuel", 0.0), [energy_line])
        add("gas_carbon", parts.get("gas_carbon", 0.0), [energy_line])

    # -- hydrogen: electrolyser capital, its fixed O&M, its water, and the buffer
    electrolyser_total = _group_eur_per_t(row, "electrolyser", steel_t)
    if electrolyser_total > 0:
        el_cfg = assumptions["electrolyser"]
        h2_kg = (_value(row, "h2_produced_kt") or 0.0) * 1e6
        # The link's variable opex is per MWh of electricity drawn, and the reported
        # H2 output fixes that exactly through the nominal conversion efficiency.
        el_mwh = h2_kg * el_cfg["efficiency_kwh_per_kg"] / 1000.0
        water = el_cfg["varopex_eur_per_mwh_el"] * el_mwh / steel_t
        capital = max(electrolyser_total - water, 0.0)
        fraction = electrolyser_om_fraction(assumptions)
        capacity = _value(row, "electrolyser_gw")
        utilisation = _value(row, "electrolyser_utilization")
        lines = [
            ("built", f"{capacity * 1000:,.0f} MW input" if capacity is not None else None),
            ("utilisation", f"{utilisation * 100:.0f}%" if utilisation is not None else None),
            ("H₂ produced", f"{h2_kg / max(steel_t, 1):,.0f} kg H₂ / t steel"),
        ]
        lines = [line for line in lines if line[1] is not None]
        add("electrolyser_capex", capital * (1.0 - fraction), lines)
        add("electrolyser_fom", capital * fraction, lines[:2])
        # The variable-opex quote itself is a config constant, so it travels in spec().
        add("electrolyser_water", water,
            [lines[-1], ("electricity drawn", f"{el_mwh / max(steel_t, 1):,.2f} MWh / t steel")])

    buffer_hours = _value(row, "h2_buffer_hours_dri")
    buffer_gwh = _value(row, "h2_buffer_gwh")
    add("h2_buffer", _group_eur_per_t(row, "h2_buffer", steel_t),
        [("size", f"{buffer_gwh:,.2f} GWh LHV" if buffer_gwh is not None else None),
         ("cover", f"{buffer_hours:,.0f} h of DRI demand" if buffer_hours is not None else None)])

    # -- electricity system: renewables by technology, battery, grid, HVDC
    res_total = _group_eur_per_t(row, "res", steel_t)
    if res_total > 0:
        fractions = res_om_fractions(assumptions)
        shares = {tech: _value(row, f"lcoe_{field_stem(tech)}_eur_per_mwh") or 0.0
                  for tech, _ in RES_TECH_LABELS}
        per_tech = _split_by_shares(res_total, shares)
        # A geography may build a technology the LCOE columns do not itemise; fall
        # back to one undivided renewables leaf rather than dropping the cost.
        if not per_tech:
            add("res_other", res_total)
        for tech, _ in RES_TECH_LABELS:
            tech_total = per_tech.get(tech)
            if not tech_total:
                continue
            stem = field_stem(tech)
            own = _value(row, f"lcoe_{stem}_own_eur_per_mwh")
            capacity_factor = _value(row, f"cf_{stem}")
            lines = [
                ("own LCOE", f"{own:,.1f} €/MWh" if own is not None else None),
                ("capacity factor",
                 f"{capacity_factor * 100:.0f}%" if capacity_factor is not None else None),
            ]
            lines = [line for line in lines if line[1] is not None]
            fraction = fractions.get(tech, 0.0)
            add(f"res_{stem}_capex", tech_total * (1.0 - fraction), lines)
            add(f"res_{stem}_fom", tech_total * fraction, lines)

    battery_total = _group_eur_per_t(row, "battery", steel_t)
    if battery_total > 0:
        power_fraction = battery_power_fraction(assumptions)
        battery_mwh = _value(row, "battery_mwh_opt")
        lines = [("energy built",
                  f"{battery_mwh:,.0f} MWh" if battery_mwh is not None else None)]
        lines = [line for line in lines if line[1] is not None]
        add("battery_power", battery_total * power_fraction, lines)
        add("battery_energy", battery_total * (1.0 - power_fraction), lines)

    grid_total = _group_eur_per_t(row, "grid", steel_t)
    if grid_total > 0:
        connection_share = _value(row, "lcoe_grid_connection_eur_per_mwh") or 0.0
        energy_share = _value(row, "lcoe_grid_energy_eur_per_mwh") or 0.0
        parts = _split_by_shares(grid_total, {"connection": connection_share,
                                              "energy": energy_share})
        capacity_factor = _value(row, "cf_grid_connection")
        connection_lines = [
            ("levelised", f"{connection_share:,.1f} €/MWh delivered"),
            ("utilisation",
             f"{capacity_factor * 100:.0f}%" if capacity_factor is not None else None),
        ]
        # The connection pays annuitised capex plus a yearly capacity charge per MW,
        # both on the same built MW — so their config ratio splits the reported cost.
        grid_cfg = assumptions["grid"]
        connection_capex = (annuity_factor(wacc, grid_cfg["connection_lifetime_years"])
                            * grid_cfg["connection_capex_eur_per_mw"])
        connection_parts = _split_by_shares(
            parts.get("connection", 0.0),
            {"capex": connection_capex, "fee": grid_cfg["fee_eur_per_mw_per_year"]})
        add("grid_connection_capex", connection_parts.get("capex", 0.0), connection_lines)
        add("grid_capacity_fee", connection_parts.get("fee", 0.0), connection_lines)
        # The imported energy is priced at the day-ahead market plus a flat
        # volumetric fee, which the report separates per imported MWh.
        price = _value(row, "grid_price_eur_per_mwh") or 0.0
        fee = _value(row, "grid_fee_eur_per_mwh") or 0.0
        energy_parts = _split_by_shares(parts.get("energy", 0.0),
                                        {"market": price, "fee": fee})
        # The market price is genuinely per scenario (it is an average over the hours
        # this plant chose to import); the volumetric fee is a flat config quote and
        # travels with the spec instead of being repeated here.
        add("grid_market", energy_parts.get("market", 0.0),
            [("average price paid", f"{price:,.1f} €/MWh imported")])
        add("grid_fee", energy_parts.get("fee", 0.0))

    transmission = _value(row, "lcoe_transmission_eur_per_mwh")
    add("transmission", _group_eur_per_t(row, "transmission", steel_t),
        [("levelised", f"{transmission:,.1f} €/MWh delivered")
         if transmission is not None else None])

    # -- getting the iron to a furnace that is somewhere else
    iron_kt = _value(row, "iron_shipped_kt")
    steel_kt = _value(row, "steel_shipped_kt")
    distance = _value(row, "transport_km")
    add("transport", _group_eur_per_t(row, "transport", steel_t),
        [("distance", f"{distance:,.0f} km" if distance else None),
         ("iron shipped", f"{iron_kt:,.0f} kt/yr" if iron_kt else None),
         ("steel shipped", f"{steel_kt:,.0f} kt/yr" if steel_kt else None)])
    add("destination_power", _group_eur_per_t(row, "destination_power", steel_t))

    # -- the solid stores that make turndown possible
    for group, hours_column, label in (("iron_store", "iron_store_hours_steel", "iron"),
                                       ("steel_store", "steel_store_hours_steel", "steel")):
        hours = _value(row, hours_column)
        kilotonnes = _value(row, f"{group}_kt")
        add(group, _group_eur_per_t(row, group, steel_t),
            [("size", f"{kilotonnes:,.1f} kt" if kilotonnes is not None else None),
             ("cover", f"{hours:,.0f} h of demand" if hours is not None else None)])

    return leaves, inputs


# ---- the alternative taxonomy's coarse bands -----------------------------

def alternative_lcos_bands(row: pd.Series, assumptions: dict,
                           h2_lhv_kwh_per_kg: float) -> dict:
    """The circulated taxonomy's LCOS bands (€/t steel), which sum to LCOS.

    It cuts the plant differently from the reports in two places. Process plant
    separates capital from fixed O&M. Electricity is divided three ways — the share
    that made the hydrogen, the share the EAF melts with, and the rest — with the
    rest taken as a residual so the stack cannot drift from LCOS.
    """
    steel_t = (_value(row, "steel_produced_mt") or 0.0) * 1e6
    if steel_t <= 0:
        return {}

    process_total = _group_eur_per_t(row, "process", steel_t)
    om_fractions = process_om_fractions(assumptions)
    fixed_om = 0.0
    for link, _ in PROCESS_PLANTS:
        plant_cost = _value(row, f"plant_{field_stem(link)}_eur_per_t") or 0.0
        fixed_om += plant_cost * om_fractions[link]

    electricity_total = sum(
        _group_eur_per_t(row, group, steel_t)
        for group in ("res", "grid", "battery", "transmission", "destination_power")
    )
    lcoe = _value(row, "lcoe_eur_per_mwh") or 0.0
    h2_mwh_lhv = (_value(row, "h2_produced_kt") or 0.0) * 1e6 * h2_lhv_kwh_per_kg / 1000.0
    hydrogen_electricity = ((_value(row, "lcoh_electricity_eur_per_mwh_lhv") or 0.0)
                            * h2_mwh_lhv / steel_t)
    eaf_built = (_value(row, "eaf_t_per_h_opt") or 0.0) > 0
    charge_state, _ = EAF_CHARGE[row["route"]]
    eaf_mwh_per_t = (assumptions["eaf"]["charge"][charge_state]["el_mwh_per_t"]
                     if eaf_built else 0.0)
    # An export route's furnace melts on bought power at the destination, so
    # its band is priced there rather than at the origin's levelised cost.
    if str(row["route"]).endswith("-export"):
        eaf_lcoe = (assumptions["destination"]["price_eur_per_mwh"]
                    + assumptions["grid"]["fee_eur_per_mwh"])
    else:
        eaf_lcoe = lcoe
    eaf_electricity = eaf_mwh_per_t * eaf_lcoe

    return {
        "ore": _group_eur_per_t(row, "ore_consumables", steel_t),
        "capex": process_total - fixed_om,
        "fixed_om": fixed_om,
        "hydrogen": (_group_eur_per_t(row, "electrolyser", steel_t)
                     + _group_eur_per_t(row, "h2_buffer", steel_t)
                     + hydrogen_electricity),
        "eaf": eaf_electricity,
        "rest": max(electricity_total - hydrogen_electricity - eaf_electricity, 0.0),
        "gas": _group_eur_per_t(row, "gas", steel_t),
        "transport": _group_eur_per_t(row, "transport", steel_t),
        "store": (_group_eur_per_t(row, "iron_store", steel_t)
                  + _group_eur_per_t(row, "steel_store", steel_t)),
        # Carried for the hover: what the single hydrogen block is actually made of,
        # and the electricity total the residual is taken from.
        "_hydrogen_electricity": hydrogen_electricity,
        "_electrolyser": _group_eur_per_t(row, "electrolyser", steel_t),
        "_h2_buffer": _group_eur_per_t(row, "h2_buffer", steel_t),
        "_electricity_total": electricity_total,
        "_eaf_mwh_per_t": eaf_mwh_per_t,
    }


def water_per_mwh_lhv(assumptions: dict, h2_lhv_kwh_per_kg: float) -> float:
    """The electrolyser's variable opex per MWh of H2 LHV — a pure config constant.

    Variable opex is charged per MWh of *electricity* drawn, and the nominal
    conversion efficiency fixes how much that is per MWh of hydrogen out.
    """
    cfg = assumptions["electrolyser"]
    return (cfg["varopex_eur_per_mwh_el"]
            * cfg["efficiency_kwh_per_kg"] / h2_lhv_kwh_per_kg)


def alternative_carrier_bands(row: pd.Series, assumptions: dict,
                              h2_lhv_kwh_per_kg: float) -> dict:
    """The circulated taxonomy's LCOE and LCOH bands, in the reports' own units.

    Same totals as the reported parts; it just separates capital from fixed O&M on
    the renewables, and capital / fixed O&M / water on the electrolyser.
    """
    bands: dict[str, dict] = {"lcoe": {}, "lcoh": {}}

    renewables = _value(row, "lcoe_renewables_eur_per_mwh")
    if renewables:
        # Blend the per-technology O&M fractions by what each contributes to LCOE,
        # so the pair of rows carries the mix this scenario actually built.
        fractions = res_om_fractions(assumptions)
        weights = {tech: _value(row, f"lcoe_{field_stem(tech)}_eur_per_mwh") or 0.0
                   for tech, _ in RES_TECH_LABELS}
        weighted = sum(weights[tech] * fractions.get(tech, 0.0) for tech in weights)
        total_weight = sum(weights.values())
        fraction = weighted / total_weight if total_weight else 0.0
        bands["lcoe"]["res_fom"] = renewables * fraction
        bands["lcoe"]["res_capex"] = renewables * (1.0 - fraction)
    for column, key in (("lcoe_storage_eur_per_mwh", "storage"),
                        ("lcoe_grid_connection_eur_per_mwh", "grid_connection"),
                        ("lcoe_grid_energy_eur_per_mwh", "grid_energy"),
                        ("lcoe_transmission_eur_per_mwh", "transmission")):
        value = _value(row, column)
        if value:
            bands["lcoe"][key] = value

    electrolyser = _value(row, "lcoh_electrolyser_eur_per_mwh_lhv")
    if electrolyser:
        water = min(water_per_mwh_lhv(assumptions, h2_lhv_kwh_per_kg), electrolyser)
        capital = electrolyser - water
        fraction = electrolyser_om_fraction(assumptions)
        bands["lcoh"]["electrolyser_water"] = water
        bands["lcoh"]["electrolyser_fom"] = capital * fraction
        bands["lcoh"]["electrolyser_capex"] = capital * (1.0 - fraction)
    for column, key in (("lcoh_h2_storage_eur_per_mwh_lhv", "storage"),
                        ("lcoh_electricity_eur_per_mwh_lhv", "electricity")):
        value = _value(row, column)
        if value:
            bands["lcoh"][key] = value
    return bands


# ---- the published leaf specification ------------------------------------

def spec(assumptions: dict) -> list:
    """[(key, label, parent group, colour, [(constant, value)])] in stack order.

    The constants are the config quotes behind each leaf — the same for every
    scenario the dashboard covers, so they travel once with the payload instead of
    being repeated in each record.
    """
    wacc = assumptions["finance"]["default_wacc"]
    wacc_line = ("WACC", f"{wacc * 100:.1f}%")
    entries = []

    def add(key, label, group, colour, constants=()):
        entries.append([key, label, group, colour,
                        [list(pair) for pair in constants if pair is not None]])

    # Ore is priced per tonne of the reduction step's own output, and the ore grade
    # each route tolerates differs — so the quote that applies is per scenario, and
    # travels with the record rather than here.
    add("ore", "Iron ore", "feedstock", "#E2B681")
    add("consumables", "EAF consumables", "feedstock", "#C99A5E",
        [("electrodes, fluxes, alloys, carbon",
          f"{assumptions['eaf']['consumables_eur_per_t']:,.0f} €/t steel")])

    # Process plants: two leaves each, capital then fixed O&M.
    shades = {"dri-h2": ("#33434D", "#5A6B77"), "dri-ng": ("#3D4E59", "#687985"),
              "moe": ("#2B3A44", "#54656F"),
              "ew": ("#25333B", "#4C5D66"), "eaf": ("#3A4A54", "#6E7F89")}
    for link, label in PROCESS_PLANTS:
        plant = assumptions[link]
        stem = field_stem(link)
        capex_colour, om_colour = shades[link]
        add(f"{stem}_capex", f"{label} — capital", "process", capex_colour,
            [("capex quote", f"{plant['capex_per_t_per_year_eur']:,.0f} €/(t·yr)"),
             ("lifetime", f"{plant['lifetime_years']:.0f} y"), wacc_line])
        add(f"{stem}_fom", f"{label} — fixed O&M", "process", om_colour,
            [("fixed opex", f"{plant['opex_per_t_per_year_eur']:,.1f} €/(t·yr)"),
             ("as a share of capex",
              f"{plant['opex_per_t_per_year_eur'] / plant['capex_per_t_per_year_eur'] * 100:.1f}% / yr")])

    gas = assumptions["natural_gas"]
    add("gas_fuel", "Natural gas — fuel", "gas", "#525F6A",
        [("price", f"{gas['price_eur_per_mwh']:,.1f} €/MWh LHV")])
    add("gas_carbon", "Natural gas — CO₂ price", "gas", "#7A8792",
        [("carbon price", f"{gas.get('co2_price_eur_per_t', 0.0):,.0f} €/t CO₂"),
         ("emission factor", f"{gas['co2_t_per_mwh']:.2f} t CO₂ / MWh LHV")])

    el_cfg = assumptions["electrolyser"]
    add("electrolyser_capex", "Electrolyser — capital", "hydrogen", "#6FA875",
        [("capex quote", f"{el_cfg['capex_per_mw_eur'] / 1e6:,.2f} M€/MW"),
         ("lifetime", f"{el_cfg['lifetime_years']:.0f} y"), wacc_line,
         ("efficiency", f"{el_cfg['efficiency_kwh_per_kg']:,.0f} kWh/kg H₂")])
    add("electrolyser_fom", "Electrolyser — fixed O&M", "hydrogen", "#91C096",
        [("fixed opex", f"{el_cfg['opex_per_mw_per_year_eur'] / 1e3:,.0f} k€/MW·yr")])
    add("electrolyser_water", "Electrolyser — water & variable opex", "hydrogen", "#B4D4B8",
        [("variable opex", f"{el_cfg['varopex_eur_per_mwh_el']:,.2f} €/MWh electricity")])
    # Quoted in €/MWh rather than k€/MWh: the salt-cavern sensitivity drops this from
    # 10,000 to 350, which rounds away entirely on the larger unit.
    buffer = assumptions["h2_buffer"]
    add("h2_buffer", "H₂ buffer store", "hydrogen", "#70D2F0",
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
        add(f"res_{stem}_capex", f"{label} — capital", "electricity", capex_colour,
            [("capex quote", f"{cfg['capex_per_mw_eur'] / 1e6:,.2f} M€/MW"),
             ("lifetime", f"{cfg['lifetime_years']:.0f} y"), wacc_line])
        add(f"res_{stem}_fom", f"{label} — fixed O&M", "electricity", om_colour,
            [("fixed opex", f"{cfg['opex_per_mw_per_year_eur'] / 1e3:,.0f} k€/MW·yr")])
    add("res_other", "Renewables (not itemised)", "electricity", "#7FA8C0")

    battery = assumptions["battery"]
    add("battery_power", "Battery — power capacity", "electricity", "#D75674",
        [("capex quote", f"{battery['capex_per_mw_eur'] / 1e6:,.2f} M€/MW"),
         ("lifetime", f"{battery['lifetime_years']:.0f} y"), wacc_line])
    add("battery_energy", "Battery — energy capacity", "electricity", "#E58AA0",
        [("capex quote", f"{battery['capex_per_mwh_eur'] / 1e6:,.2f} M€/MWh"),
         ("duration", f"{battery['max_hours']:.0f} h at rated power"),
         ("round-trip efficiency", f"{battery['efficiency_roundtrip'] * 100:.0f}%")])

    grid = assumptions["grid"]
    add("grid_connection_capex", "Grid connection — capital", "electricity", "#71828F",
        [("capex quote", f"{grid['connection_capex_eur_per_mw'] / 1e3:,.0f} k€/MW"),
         ("lifetime", f"{grid['connection_lifetime_years']:.0f} y"), wacc_line])
    add("grid_capacity_fee", "Grid connection — capacity charge", "electricity", "#98A5AE",
        [("capacity charge (Leistungspreis)",
          f"{grid['fee_eur_per_mw_per_year'] / 1e3:,.0f} k€/MW·yr")])
    add("grid_market", "Grid energy — market price", "electricity", "#B7C1C8")
    add("grid_fee", "Grid energy — volumetric fee", "electricity", "#D3DAE0",
        [("volumetric charge (Arbeitspreis)", f"{grid['fee_eur_per_mwh']:,.1f} €/MWh imported")])
    add("transmission", "Transmission (HVDC)", "electricity", "#83D1DD")
    destination = assumptions["destination"]
    add("destination_power", "Destination power", "electricity", "#0293D2",
        [("country", f"{destination['country']}"),
         ("flat price", f"{destination['price_eur_per_mwh']:,.0f} €/MWh")])
    freight = assumptions["transport"]
    add("transport", "Freight", "storage", "#BDCCD9",
        [(f"{mode}, {commodity}",
          f"{freight[mode][commodity]['eur_per_t']:,.0f} €/t "
          f"+ {freight[mode][commodity]['eur_per_t_km'] * 1000:,.2f} €/t per 1000 km")
         for mode in ("sea", "rail") for commodity in ("iron", "steel")])

    add("iron_store", "Iron stockpile", "storage", "#D75674",
        [("capex quote", f"{assumptions['iron_store']['capex_per_t_eur']:,.0f} €/t"),
         ("lifetime", f"{assumptions['iron_store']['lifetime_years']:.0f} y")])
    add("steel_store", "Steel inventory", "storage", "#EE8DA3",
        [("capex quote", f"{assumptions['steel_store']['capex_per_t_eur']:,.0f} €/t"),
         ("lifetime", f"{assumptions['steel_store']['lifetime_years']:.0f} y")])
    return entries
