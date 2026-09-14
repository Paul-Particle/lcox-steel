#!/usr/bin/env python3
"""Build the assumptions page: every input the run was priced on, and where a
scenario overlay moves one.

Values are read out of `config/assumptions.yaml` rather than written here, so a
number on the page cannot drift from the number that solved. What this file
holds is the reading of them — what each one is, what it is worth knowing about,
and which of them do not reach this run at all.

Two things keep that honest. `SKIPPED` names the assumptions the page leaves out
and why, and every leaf of the file must be either rendered or skipped, or the
build stops — so an assumption added later cannot go quietly unlisted. The
overlays are read the same way: `config/assumptions_{scenario}.yaml` for every
scenario in `config/scenarios.csv`, whose leaves are exactly the deviations.

Output: results/html/assumptions.html — body-only, the hub's second tab.
"""
import fnmatch
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))    # sibling build_dashboard

from build_dashboard import DASHBOARD_SCENARIOS, HTML_DIR, font_css   # noqa: E402

REPO = Path(__file__).resolve().parents[3]
CONFIG = REPO / "config"
RESULTS = REPO / "results"
TEMPLATE_HTML = Path(__file__).with_name("assumptions_template.html")
OUT = HTML_DIR / "assumptions.html"

# What the page leaves out, and the reason it gives for leaving it out. Patterns
# are fnmatch over the dotted path. Every one of these was checked against the
# run rather than assumed: the scenario table names no offshore tech and no
# multi-site overlay exists, so neither the offshore costs nor the transmission
# block can reach a network here.
SKIPPED = {
    "res.wind-offshore.*":
        "No offshore row in <code>config/scenarios.csv</code>, so no offshore "
        "generator is built anywhere in this run.",
    "anchor.min_capacity_mw":
        "A floor on what the anchor cell must build, and it binds only on a "
        "multi-site run. Every run here is single-site best-site P95.",
    "transmission.*":
        "Prices the HVDC link between a remote RES site and the plant. Only a "
        "multi-site run adds one, and none of these are multi-site.",
    "plant.dri_mt_per_year":
        "Sizes the electrolyser for the hydrogen-only route, which makes no "
        "steel and is not plotted.",
    "plant.availability_target":
        "The same: the hydrogen-only route's availability floor. A steel route "
        "leaves electrolyser sizing to the optimiser.",
    "natural_gas.co2_price_eur_per_t":
        "Zero in every scenario here, so no route pays anything for its carbon "
        "and no cost on any chart comes from it.",
    "emissions.electricity_t_co2e_per_mwh.*.delegated_act":
        "A second emission basis the run did not use.",
    "emissions.electricity_t_co2e_per_mwh.*.lifecycle":
        "A third emission basis the run did not use.",
    "emissions.gas_upstream_t_co2e_per_mwh.delegated_act": "As above.",
    "emissions.gas_upstream_t_co2e_per_mwh.lifecycle": "As above.",
    "emissions.freight_kg_co2e_per_t_km.*.delegated_act": "As above.",
    "emissions.freight_kg_co2e_per_t_km.*.lifecycle": "As above.",
}

# Plain rows: (path, label, unit, note). The path is both where the value comes
# from and what an overlay is matched against.
SECTIONS = [
    ("Money, and how big the plant is",
     "One discount rate annuitises every capital cost in the model, and one flat "
     "steel load sets the size of everything downstream of it.",
     [("finance.default_wacc", "Discount rate (WACC)", "",
       "Applied to every capex through the standard annuity over that item's own "
       "lifetime. One rate for every technology and every country."),
      ("plant.steel_mt_per_year", "Steel output", "Mt/yr",
       "A flat hourly load on the steel bus. Demand is not flexible; the steel "
       "store is what lets production move around it."),
      ("plant.h2_intensity_kg_per_t_dri", "Hydrogen per tonne of DRI", "kg/t",
       "Converts hydrogen into iron in both the H2 shaft and the blend shaft.")]),

    ("Renewables, electrolyser, battery, hydrogen store",
     "The supply side. Capacities are all extendable — the optimiser sizes them, "
     "so what is set here is only what a MW or a MWh costs.",
     [("res.wind-onshore.capex_per_mw_eur", "Onshore wind capex", "€/MW",
       "IRENA 2025, one global figure for every geography, as with every other "
       "cross-country input in the file."),
      ("res.wind-onshore.opex_per_mw_per_year_eur", "Onshore wind fixed opex", "€/MW/yr", ""),
      ("res.wind-onshore.lifetime_years", "Onshore wind lifetime", "years", ""),
      ("res.solar.capex_per_mw_eur", "Solar capex", "€/MW", "IRENA 2025, as above."),
      ("res.solar.opex_per_mw_per_year_eur", "Solar fixed opex", "€/MW/yr", ""),
      ("res.solar.lifetime_years", "Solar lifetime", "years", ""),
      ("electrolyser.capex_per_mw_eur", "Electrolyser capex", "€/MW",
       "Alkaline. The register carries alkaline and PEM side by side; this run "
       "is priced on the alkaline column."),
      ("electrolyser.opex_per_mw_per_year_eur", "Electrolyser fixed opex", "€/MW/yr",
       "2 % of capex."),
      ("electrolyser.lifetime_years", "Electrolyser lifetime", "years", ""),
      ("electrolyser.efficiency_kwh_per_kg", "Electrolyser efficiency", "kWh/kg H2",
       "Electricity in per kg of hydrogen out."),
      ("electrolyser.varopex_eur_per_mwh_el", "Electrolyser variable opex", "€/MWh el",
       "Water and consumables, charged on the electricity drawn."),
      ("battery.capex_per_mw_eur", "Battery power capex", "€/MW",
       "One bidirectional inverter, charging and discharging at that rating."),
      ("battery.capex_per_mwh_eur", "Battery energy capex", "€/MWh",
       "Priced and sized separately from power, so the optimiser picks the "
       "duration rather than being handed one."),
      ("battery.lifetime_years", "Battery lifetime", "years", ""),
      ("battery.efficiency_roundtrip", "Battery round-trip efficiency", "", ""),
      ("h2_buffer.capex_per_mwh_eur", "Hydrogen store capex", "€/MWh LHV",
       "Salt cavern, and deliberately a lower bound — the register's "
       "tank-with-compressor figure is 60 050 €/MWh. Size is optimised."),
      ("h2_buffer.lifetime_years", "Hydrogen store lifetime", "years", "")]),

    ("Stores on the iron and steel side",
     "Both are cheap on purpose: what they are for is letting production move "
     "away from delivery, not earning their own keep.",
     [("iron_store.capex_per_t_eur", "Iron stockpile capex", "€/t",
       "Electrowon plates, or briquettes waiting on a ship. A route that hands "
       "the furnace hot or liquid iron has nowhere to put it and gets none."),
      ("iron_store.lifetime_years", "Iron stockpile lifetime", "years", ""),
      ("steel_store.capex_per_t_eur", "Steel inventory capex", "€/t",
       "Covered yard and handling. This is the supply-side way of representing "
       "periodic rather than hour-by-hour delivery."),
      ("steel_store.lifetime_years", "Steel inventory lifetime", "years", ""),
      ("steel_store.max_weeks", "Steel inventory cap", "weeks of output",
       "The binding limit, not the cost — it stops the flexibility degenerating "
       "into whole-year arbitrage.")]),

    ("Natural gas",
     "One flat price for every geography, which is also what keeps the fossil "
     "benchmark comparable between countries.",
     [("natural_gas.price_eur_per_mwh", "Gas price", "€/MWh LHV",
       "Midpoint of the register's low/high benchmark."),
      ("natural_gas.co2_t_per_mwh", "Gas combustion CO2", "t/MWh LHV",
       "Accounting only. Nothing in this run puts a price on it.")]),

    ("Grid connection",
     "What a grid-connected run pays to be connected, on top of the hourly "
     "wholesale price it pays for the energy itself. The connection is sized by "
     "the optimiser. One set of charges for every area.",
     [("grid.connection_capex_eur_per_mw", "Connection capex", "€/MW", ""),
      ("grid.connection_lifetime_years", "Connection lifetime", "years", ""),
      ("grid.fee_eur_per_mw_per_year", "Capacity fee", "€/MW/yr",
       "At this plant's roughly 8 000 full-load hours it works out near "
       "19 €/MWh, between the German HV band-load position and Hydrogen "
       "Europe's 29.3 €/MWh EU average."),
      ("grid.fee_eur_per_mwh", "Volumetric fee", "€/MWh",
       "Zero, so every geography is compared on the single capacity charge.")]),

    ("Where the export routes melt their iron, and what the freight costs",
     "An <i>export</i> route ships iron and melts it at the destination, whose "
     "power is bought at that market's own hourly price. Distances are an "
     "assumption rather than geography — they depend on the port and the routing.",
     [("destination.area", "Destination", "",
       "An area from the registry, priced and counted against its own hourly "
       "series — the same download a grid run uses, so the bill and the carbon "
       "come from one file."),
      ("transport.deliver_finished_steel", "Deliver finished steel too", "",
       "Off: the model stops at the plant gate and only an export route pays "
       "freight. On, every route would also deliver its steel over the same legs."),
      ("transport.sea.iron.eur_per_t_km", "Sea freight, iron", "€/t·km",
       "Iron travels in fitted holds under inert gas, so it carries a premium "
       "over ore."),
      ("transport.sea.steel.eur_per_t_km", "Sea freight, steel", "€/t·km", ""),
      ("transport.sea.iron.eur_per_t", "Sea freight, iron, per tonne", "€/t", ""),
      ("transport.sea.steel.eur_per_t", "Sea freight, steel, per tonne", "€/t", ""),
      ("transport.rail.iron.eur_per_t_km", "Rail freight, iron", "€/t·km",
       "An order of magnitude dearer per t·km than sea, which is why iron moves "
       "between continents by ship and barely moves overland."),
      ("transport.rail.steel.eur_per_t_km", "Rail freight, steel", "€/t·km", ""),
      ("transport.rail.iron.eur_per_t", "Rail freight, iron, per tonne", "€/t",
       "Zero because the European rates this is anchored on are quoted all-in."),
      ("transport.rail.steel.eur_per_t", "Rail freight, steel, per tonne", "€/t", "")]),
]

# The process steps, as a matrix: a column per step, a row per field. Blank where
# a step has no such field — a press has no ore bill, a furnace burns no gas.
PROCESS_STEPS = [("dri-h2", "H2 shaft"), ("dri-mix", "Blend shaft"),
                 ("dri-ng", "NG shaft"), ("eaf", "EAF"), ("moe", "MOE cell"),
                 ("ew", "Electrowinning"), ("briquetting", "Briquetting")]
PROCESS_FIELDS = [
    ("capex_per_t_per_year_eur", "Capex", "€/(t/yr)"),
    ("opex_per_t_per_year_eur", "Fixed opex", "€/(t/yr)"),
    ("lifetime_years", "Lifetime", "years"),
    ("el_mwh_per_t", "Electricity", "MWh/t"),
    ("h2_preheat_el_mwh_per_t", "H2 preheat electricity", "MWh/t"),
    ("gas_mwh_per_t", "Gas", "MWh LHV/t"),
    ("ore_eur_per_t", "Ore", "€/t"),
    ("iron_t_per_t_steel", "Iron per t steel", "t/t"),
    ("consumables_eur_per_t", "Consumables", "€/t"),
    ("yield_t_per_t", "Yield", "t/t"),
    ("p_min_pu", "Minimum load while built", "of capacity"),
    ("reductant_store_hours", "Reductant buffer", "hours"),
]

# How hot the iron arrives, and what melting it therefore costs. Nearly all of an
# EAF's electricity is the melt, so iron that shows up hot has had part of that
# bill paid upstream.
CHARGE_STATES = [
    ("briquettes", "Briquettes (HBI), 25 °C",
     "The dearest charge: dense, and carrying the gangue the shaft could not remove."),
    ("plates", "Electrowon plates or cold MOE iron, 25 °C",
     "Cold like the briquettes, but nearly gangue-free."),
    ("hot", "Sponge iron, ~650 °C",
     "Walked straight from the shaft to the furnace."),
    ("liquid", "Liquid iron, ~1550 °C",
     "From the MOE cell, so the melt is already done — what is left is transfer, "
     "holding and superheat."),
]


def leaves(node, prefix=""):
    """Every scalar in a nested dict, as (dotted path, value)."""
    out = []
    for key, value in node.items():
        path = f"{prefix}{key}"
        if isinstance(value, dict):
            out.extend(leaves(value, f"{path}."))
        else:
            out.append((path, value))
    return out


def read(assumptions: dict, path: str):
    """One value out of the nested assumptions by its dotted path."""
    node = assumptions
    for key in path.split("."):
        node = node[key]
    return node


def fmt(value) -> str:
    """A number as the page prints it: grouped thousands, no trailing noise."""
    if isinstance(value, bool):
        return "no" if value is False else "yes"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        whole = value == int(value) and abs(value) >= 1000
        return f"{value:,.0f}" if whole else f"{value:g}"
    return str(value)


def main() -> None:
    """Read the assumptions, the overlays and the scenario table; write the page."""
    assumptions = yaml.safe_load((CONFIG / "assumptions.yaml").read_text())
    config = yaml.safe_load((CONFIG / "config.yaml").read_text())
    scenarios = pd.read_csv(CONFIG / "scenarios.csv", comment="#")
    run_scenarios = list(dict.fromkeys(scenarios["scenario"]))

    # The overlays, which are exactly the deviations: an overlay file holds only
    # what it moves, so its leaves are the answer without a diff.
    moved = {}
    for overlay_path in sorted(CONFIG.glob("assumptions_*.yaml")):
        scenario = overlay_path.stem.removeprefix("assumptions_")
        if scenario not in run_scenarios:
            continue
        for path, value in leaves(yaml.safe_load(overlay_path.read_text())):
            moved.setdefault(path, {})[scenario] = value

    shown = set()
    html = []

    # ---- header ---------------------------------------------------------
    html.append(
        '<div class="as-brand"><div class="as-dot"></div>'
        '<span>FC Architects · green steel model</span></div>'
        '<h1 class="as-title">What the run was priced on</h1>'
        '<div class="as-intro">'
        '<p class="as-sub">Every input behind the numbers on the other tabs, read '
        'straight out of <b>config/assumptions.yaml</b> at build time. A scenario '
        'can move any of them through an overlay file of its own, deep-merged over '
        'this one; a value that some scenario moves is marked '
        '<span class="as-moved">moved</span> and collected in <b>Where the '
        'scenarios differ</b> below. Assumptions that do not reach this run are '
        'left out, and listed at the end with the reason.</p>'
        '<p class="as-warn"><b>These are placeholders, not research.</b> Even the '
        'sourced figures are single points on scales that want a sensitivity run '
        'of their own, and several — the pre-commercial cell capexes above all — '
        'are order-of-magnitude eyeballs. Read a level as a starting position, and '
        'a comparison between routes as the thing the model is actually for.</p>'
        '</div>')

    # ---- what ran -------------------------------------------------------
    rows = []
    for scenario in run_scenarios:
        block = scenarios[scenarios["scenario"] == scenario]
        routes = " · ".join(dict.fromkeys(block["route"])).replace("|", " · ")
        areas = " · ".join(dict.fromkeys(block["area"]))
        period = f"{block['start_date'].iloc[0]}–{block['end_date'].iloc[0]}"
        overlay = ("<br>".join(f'<code>{path}</code> → {fmt(by[scenario])}'
                               for path, by in sorted(moved.items())
                               if scenario in by)
                   or '<span style="opacity:.5">—</span>')
        on_dash = "yes" if scenario in DASHBOARD_SCENARIOS else "no"
        rows.append(f'<tr><td class="k">{scenario}</td><td class="n">{routes}</td>'
                    f'<td class="n">{areas}</td><td class="u">{period[:4]}</td>'
                    f'<td class="n">{overlay}</td><td class="u">{on_dash}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>The runs these assumptions were used for</h2>'
        '<p class="lead">Fifteen scenarios. <code>all-routes</code> is every route '
        'in the model — five steel routes, their five export twins, and a '
        'hydrogen-only run that makes no steel and is therefore absent from the '
        'charts. <code>all-areas</code> resolves through the area registry, which '
        'is what sends a grid run to Brazil\'s four submarkets and Australia\'s '
        'five NEM regions instead of the country. Every run is hourly over the '
        'whole year — 8 760 snapshots — and solved with HiGHS. Only the five '
        'scenarios marked below are browsable on the other tabs; the ten sweep '
        'scenarios were solved and reported but are not plotted there.</p>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Scenario</th>'
        '<th>Routes</th><th>Areas</th><th>Year</th><th>Overlay</th>'
        '<th>On the dashboard</th></tr>' + "".join(rows) + '</table></div></div>')

    # ---- plain sections -------------------------------------------------
    for title, lead, rows_spec in SECTIONS:
        rows = []
        for path, label, unit, note in rows_spec:
            shown.add(path)
            mark = ('<span class="as-moved">moved</span>' if path in moved else "")
            rows.append(f'<tr><td>{label}{mark}</td>'
                        f'<td class="v num">{fmt(read(assumptions, path))}</td>'
                        f'<td class="u">{unit}</td><td class="n">{note}</td>'
                        f'<td class="k">{path}</td></tr>')
        html.append(f'<div class="as-sec"><h2>{title}</h2><p class="lead">{lead}</p>'
                    '<div class="as-scroll"><table class="as-t">'
                    '<tr><th>Assumption</th><th class="num">Value</th><th>Unit</th>'
                    '<th>Note</th><th>Key</th></tr>'
                    + "".join(rows) + '</table></div></div>')

    # ---- the process steps, as a matrix ---------------------------------
    head = "".join(f'<th class="num">{label}</th>' for _, label in PROCESS_STEPS)
    rows = []
    for field, label, unit in PROCESS_FIELDS:
        cells = []
        for step, _ in PROCESS_STEPS:
            path = f"{step}.{field}"
            if field not in assumptions[step]:
                cells.append('<td class="num" style="opacity:.3">·</td>')
                continue
            shown.add(path)
            mark = ('<span class="as-moved">*</span>' if path in moved else "")
            cells.append(f'<td class="v num">{fmt(read(assumptions, path))}{mark}</td>')
        rows.append(f'<tr><td>{label}</td><td class="u">{unit}</td>'
                    + "".join(cells) + '</tr>')
    html.append(
        '<div class="as-sec"><h2>The process steps</h2>'
        '<p class="lead">Capex is quoted per tonne of annual output capacity, the '
        'basis the industry quotes on; fixed opex likewise, and it is all-in — '
        'labour, maintenance and overhead, which is why it is an eighth of capex '
        'rather than the few per cent a maintenance-only figure would be. Every '
        'route melts in the same EAF, so its column is shared. A dot means the '
        'field does not apply to that step, and a '
        '<span class="as-moved">*</span> that some scenario moves the value.</p>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Assumption</th>'
        f'<th>Unit</th>{head}</tr>' + "".join(rows) + '</table></div></div>')

    # ---- the EAF's charge states ----------------------------------------
    rows = []
    for state, label, note in CHARGE_STATES:
        path = f"eaf.charge.{state}.el_mwh_per_t"
        shown.add(path)
        rows.append(f'<tr><td>{label}</td>'
                    f'<td class="v num">{fmt(read(assumptions, path))}</td>'
                    f'<td class="u">MWh/t</td><td class="n">{note}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>How the iron reaches the furnace</h2>'
        '<p class="lead">The EAF\'s electricity is set by how hot its charge '
        'arrives, and the route decides which state that is. An export route\'s '
        'iron crossed an ocean, so it always arrives cold — as briquettes if a '
        'shaft made it, as plates if a cell did.</p>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Charge</th>'
        '<th class="num">Electricity</th><th>Unit</th><th>Note</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- freight legs ---------------------------------------------------
    rows = []
    for origin, legs in read(assumptions, "transport.distance_km").items():
        for mode, km in legs.items():
            shown.add(f"transport.distance_km.{origin}.{mode}")
        drawn = " + ".join(f"{fmt(km)} km {mode}" for mode, km in legs.items())
        rows.append(f'<tr><td class="k">{origin}</td><td class="v">{drawn}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>Freight legs to the destination</h2>'
        '<p class="lead">Keyed by the producing country alone, so they are drawn '
        'to whichever area <code>destination.area</code> names — here Spain. '
        'Alberta is landlocked and reaches the sea at Vancouver; Ontario goes '
        'down the St. Lawrence. Moving the destination means redrawing all of '
        'them.</p><div class="as-scroll"><table class="as-t">'
        '<tr><th>From</th><th>Legs</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- emissions accounting -------------------------------------------
    basis = read(assumptions, "emissions.basis")
    shown.add("emissions.basis")
    shown.add("emissions.gas_upstream_t_co2e_per_mwh." + basis)
    # Twenty-two carriers down one column would run the section a screen deep and
    # leave the width empty, so they are dealt into three tables side by side.
    carriers = sorted(read(assumptions, "emissions.electricity_t_co2e_per_mwh").items())
    for carrier, _ in carriers:
        shown.add(f"emissions.electricity_t_co2e_per_mwh.{carrier}.{basis}")
    per_column = -(-len(carriers) // 3)
    columns = []
    for start in range(0, len(carriers), per_column):
        rows_here = "".join(f'<tr><td class="k">{carrier}</td>'
                            f'<td class="v num">{fmt(by_basis[basis])}</td></tr>'
                            for carrier, by_basis in carriers[start:start + per_column])
        columns.append('<table class="as-t" style="width:auto">'
                       '<tr><th>Generation carrier</th>'
                       f'<th class="num">t CO2e/MWh</th></tr>{rows_here}</table>')
    freight = []
    for mode in read(assumptions, "emissions.freight_kg_co2e_per_t_km"):
        shown.add(f"emissions.freight_kg_co2e_per_t_km.{mode}.{basis}")
        value = read(assumptions, f"emissions.freight_kg_co2e_per_t_km.{mode}.{basis}")
        freight.append(f'<tr><td class="k">{mode}</td>'
                       f'<td class="v num">{fmt(value)}</td>'
                       f'<td class="u">kg CO2e/t·km</td></tr>')
    gas_upstream = fmt(read(assumptions, f"emissions.gas_upstream_t_co2e_per_mwh.{basis}"))
    html.append(
        '<div class="as-sec"><h2>Emissions accounting</h2>'
        f'<p class="lead">The basis is <b>{basis}</b> — what comes out of the '
        'stack and nothing else, so renewables and nuclear are zero and an '
        'islanded run reads zero too. None of this reaches the objective: a '
        'number here never changed a solve. It is also not a CBAM or an ETS '
        'figure, because the model\'s boundary stops at the melt and the process '
        'steps\' own direct emissions — electrodes, carbon injection, the carbon '
        'in DR-grade pellets, fluxes — are outside it. What is counted is the '
        'energy and the freight. The <i>Estimated emissions</i> readout on the '
        'cost-breakdown tab is built from exactly this. Burned gas is counted at '
        f'its own combustion figure; its supply chain adds {gas_upstream} t '
        'CO2e/MWh on this basis, which is what a stack-only accounting means.</p>'
        '<div class="as-scroll" style="display:flex;gap:36px;align-items:flex-start;'
        'flex-wrap:wrap;">'
        + "".join(columns) +
        '<table class="as-t" style="width:auto"><tr><th>Freight</th>'
        f'<th class="num">Factor</th><th>Unit</th></tr>{"".join(freight)}</table>'
        '</div></div>')

    # ---- geographies ----------------------------------------------------
    solved_areas = set()
    for report in sorted(RESULTS.glob(".report_*_diag.csv")):
        solved_areas |= set(pd.read_csv(report, index_col="field").T["area"])
    rows = []
    for area, spec in config["areas"].items():
        if area not in solved_areas:
            continue
        market = spec.get("market", "—")
        until = spec.get("market_until")
        note = f"price series ends {until}" if until else ""
        rows.append(f'<tr><td class="k">{area}</td><td class="u">{spec["iso3"]}</td>'
                    f'<td class="v">{market}</td><td class="n">{note}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>Geographies, and where their prices come from</h2>'
        '<p class="lead">Weather is an ERA5 cutout per area, turned into hourly '
        'capacity factors by the best-site P95 method — the profile of the best '
        'cells in the area rather than its average. Hourly prices come from the '
        'market operator named here, and a grid run also pulls that market\'s '
        'generation mix, which is what the emission factors above are weighted '
        'by. An islanded run buys no grid power and needs neither.</p>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Area</th>'
        '<th>Country</th><th>Market</th><th>Note</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- where the scenarios differ -------------------------------------
    rows = []
    for path, by_scenario in sorted(moved.items()):
        by_value = {}
        for scenario, value in by_scenario.items():
            by_value.setdefault(value, []).append(scenario)
        values = " · ".join(
            f"{fmt(value)} <span style='opacity:.6'>({', '.join(sorted(names))})</span>"
            for value, names in sorted(by_value.items()))
        rows.append(f'<tr><td class="k">{path}</td>'
                    f'<td class="v num">{fmt(read(assumptions, path))}</td>'
                    f'<td class="n">{values}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>Where the scenarios differ</h2>'
        '<p class="lead">Three values, and nothing else: an overlay file holds '
        'only what it moves, so this is the whole of the deviation between '
        'scenarios. Everything not listed here is identical in all fifteen runs. '
        'Of these, only the MOE turndown is browsable on the cost-breakdown tab, '
        'as the <i>Sensitivity</i> control; the two sweeps were solved and '
        'reported but are not plotted.</p>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Assumption</th>'
        '<th class="num">Base</th><th>Moved to, and by which scenario</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- what is not on the page ----------------------------------------
    # Anything neither rendered nor skipped stops the build: an assumption added
    # later has to be given a reading or a reason, and cannot quietly vanish.
    all_paths = {path for path, _ in leaves(assumptions)}
    accounted = set(shown)
    for pattern in SKIPPED:
        accounted |= {p for p in all_paths if fnmatch.fnmatch(p, pattern)}
    missing = sorted(all_paths - accounted)
    if missing:
        raise SystemExit(
            f"assumptions not on the page and not in SKIPPED: {missing}\n"
            f"Add each to a section in {Path(__file__).name}, or to SKIPPED with "
            f"the reason it does not reach the run."
        )
    items = "".join(f'<li><code>{pattern}</code> — {reason}</li>'
                    for pattern, reason in SKIPPED.items())
    html.append(
        '<div class="as-sec"><h2>Left out, and why</h2>'
        '<p class="lead">These are in the assumptions file but do not reach this '
        'run, so putting them on the page would suggest they priced something '
        'here. The build fails if an assumption is neither shown above nor listed '
        'below, so nothing can go missing without someone deciding it should.</p>'
        f'<ul class="as-skip">{items}</ul>'
        '<p class="lead">Two things are also out of scope by intent rather than '
        'by irrelevance. The workflow settings in <code>config/config.yaml</code> '
        '— cutout caching, download polling, solver threads, file paths — are '
        'plumbing and change no result. And the weather itself is not an '
        'assumption but an input: ERA5 for the year named, through the turbine '
        'and panel models the config selects.</p></div>')

    html.append('<p class="as-foot">Generated from config/assumptions.yaml, its '
                'per-scenario overlays and config/scenarios.csv.</p>')

    page = (TEMPLATE_HTML.read_text()
            .replace("/*FONT_CSS*/", font_css())
            .replace("<!--CONTENT-->", "\n".join(html)))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page)
    print(f"wrote {OUT} ({OUT.stat().st_size/1e6:.2f} MB) — "
          f"{len(shown)} assumptions shown, {len(all_paths) - len(shown)} left out, "
          f"{len(moved)} moved by an overlay")


if __name__ == "__main__":
    main()
