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
# Where the assumptions are actually chosen and checked. The share link only —
# the Teams client appends the opening user's address and a client fingerprint to
# it, and neither belongs on a page other people open.
WORKING_DOC = ("https://futurecleantecharchitects.sharepoint.com/:x:/s/"
               "FutureCleantechArchitects/"
               "IQD-pQ-HYO2dS6UmBlA4bIbxAfsv6JAC66F21PJRRPw5v2w?e=b0390D")
OUT = HTML_DIR / "assumptions.html"

# What the page leaves out, and why — a ledger for whoever edits this file, not
# page copy. Patterns are fnmatch over the dotted path. Every entry was checked
# against the run rather than assumed.
SKIPPED = {
    "anchor.min_capacity_mw":
        "A floor on what the anchor cell must build, and it binds only on a "
        "multi-site run. Every run here is single-site best-site P95.",
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

# Plain rows: (path, label, unit) under a heading. The path is both where the
# value comes from and what an overlay is matched against.
SECTIONS = [
    ("Finance and plant",
     [("finance.default_wacc", "Discount rate (WACC)", ""),
      ("plant.steel_mt_per_year", "Steel output", "Mt/yr"),
      ("plant.h2_intensity_kg_per_t_dri", "Hydrogen per tonne of DRI", "kg/t")]),

    ("Renewables, electrolyser, battery, hydrogen storage",
     [("res.wind-onshore.capex_per_mw_eur", "Onshore wind capex", "€/MW"),
      ("res.wind-onshore.opex_per_mw_per_year_eur", "Onshore wind fixed opex", "€/MW/yr"),
      ("res.wind-onshore.lifetime_years", "Onshore wind lifetime", "years"),
      ("res.solar.capex_per_mw_eur", "Solar capex", "€/MW"),
      ("res.solar.opex_per_mw_per_year_eur", "Solar fixed opex", "€/MW/yr"),
      ("res.solar.lifetime_years", "Solar lifetime", "years"),
      ("res.wind-offshore.capex_per_mw_eur", "Offshore wind capex", "€/MW"),
      ("res.wind-offshore.opex_per_mw_per_year_eur", "Offshore wind fixed opex", "€/MW/yr"),
      ("res.wind-offshore.lifetime_years", "Offshore wind lifetime", "years"),
      ("electrolyser.capex_per_mw_eur", "Electrolyser capex", "€/MW"),
      ("electrolyser.opex_per_mw_per_year_eur", "Electrolyser fixed opex", "€/MW/yr"),
      ("electrolyser.lifetime_years", "Electrolyser lifetime", "years"),
      ("electrolyser.efficiency_kwh_per_kg", "Electrolyser efficiency", "kWh/kg H2"),
      ("electrolyser.varopex_eur_per_mwh_el", "Electrolyser variable opex", "€/MWh el"),
      ("battery.capex_per_mw_eur", "Battery power capex", "€/MW"),
      ("battery.capex_per_mwh_eur", "Battery energy capex", "€/MWh"),
      ("battery.lifetime_years", "Battery lifetime", "years"),
      ("battery.efficiency_roundtrip", "Battery round-trip efficiency", ""),
      ("h2_buffer.capex_per_mwh_eur", "Hydrogen store capex", "€/MWh LHV"),
      ("h2_buffer.lifetime_years", "Hydrogen store lifetime", "years")]),

    ("Iron and steel stores",
     [("iron_store.capex_per_t_eur", "Iron stockpile capex", "€/t"),
      ("iron_store.lifetime_years", "Iron stockpile lifetime", "years"),
      ("steel_store.capex_per_t_eur", "Steel inventory capex", "€/t"),
      ("steel_store.lifetime_years", "Steel inventory lifetime", "years"),
      ("steel_store.max_weeks", "Steel inventory cap", "weeks of output")]),

    ("Transmission",
     [("transmission.cost_per_mw_per_km_eur", "HVDC capex", "€/MW/km"),
      ("transmission.lifetime_years", "HVDC lifetime", "years"),
      ("transmission.losses_pct_per_1000km", "HVDC losses", "% per 1000 km"),
      ("transmission.indirect_route_factor", "Route factor", "")]),

    ("Natural gas",
     [("natural_gas.price_eur_per_mwh", "Gas price", "€/MWh LHV"),
      ("natural_gas.co2_t_per_mwh", "Gas combustion CO2", "t/MWh LHV")]),

    ("Grid connection",
     [("grid.connection_capex_eur_per_mw", "Connection capex", "€/MW"),
      ("grid.connection_lifetime_years", "Connection lifetime", "years"),
      ("grid.fee_eur_per_mw_per_year", "Capacity fee", "€/MW/yr"),
      ("grid.fee_eur_per_mwh", "Volumetric fee", "€/MWh")]),

    ("Export destination and freight rates",
     [("destination.area", "Destination", ""),
      ("transport.deliver_finished_steel", "Deliver finished steel too", ""),
      ("transport.sea.iron.eur_per_t_km", "Sea freight, iron", "€/t·km"),
      ("transport.sea.steel.eur_per_t_km", "Sea freight, steel", "€/t·km"),
      ("transport.sea.iron.eur_per_t", "Sea freight, iron, per tonne", "€/t"),
      ("transport.sea.steel.eur_per_t", "Sea freight, steel, per tonne", "€/t"),
      ("transport.rail.iron.eur_per_t_km", "Rail freight, iron", "€/t·km"),
      ("transport.rail.steel.eur_per_t_km", "Rail freight, steel", "€/t·km"),
      ("transport.rail.iron.eur_per_t", "Rail freight, iron, per tonne", "€/t"),
      ("transport.rail.steel.eur_per_t", "Rail freight, steel, per tonne", "€/t")]),
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
    ("briquettes", "Briquettes (HBI), 25 °C"),
    ("plates", "Electrowon plates or cold MOE iron, 25 °C"),
    ("hot", "Sponge iron, ~650 °C"),
    ("liquid", "Liquid iron, ~1550 °C"),
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
    # The scenarios the dashboard browses, and only those. The scenario table
    # also holds the capex and gas-price sweeps, which were solved and reported
    # but are plotted nowhere — an assumption of theirs on this page would read
    # as an input behind a chart the reader can reach, and none of them is.
    scenarios = scenarios[scenarios["scenario"].isin(DASHBOARD_SCENARIOS)]
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
        f'<a class="as-doc" href="{WORKING_DOC}" target="_blank" rel="noopener">'
        '<span class="as-doc-mark">Working document</span>'
        'Assumption selection and validation</a>'
        '<div class="as-intro">'
        '<p class="as-sub">Every input behind the numbers on the other tabs, read '
        'straight out of <b>config/assumptions.yaml</b> at build time. A scenario '
        'can move any of them through an overlay file of its own, deep-merged over '
        'this one; a value that some scenario moves is marked '
        '<span class="as-moved">moved</span> and collected under <b>Sensitivity</b> '
        'below. Not all assumptions are listed here.</p>'
        '<p class="as-warn"><b>This is a mix of validated and roughly estimated '
        'placeholders.</b> Central values, ranges, sensitivity analysis, country '
        'specific values are not fixed and tbd.</p>'
        '</div>')

    # ---- plain sections -------------------------------------------------
    for title, rows_spec in SECTIONS:
        rows = []
        for path, label, unit in rows_spec:
            shown.add(path)
            mark = ('<span class="as-moved">moved</span>' if path in moved else "")
            rows.append(f'<tr><td>{label}{mark}</td>'
                        f'<td class="v num">{fmt(read(assumptions, path))}</td>'
                        f'<td class="u">{unit}</td>'
                        f'<td class="k">{path}</td></tr>')
        html.append(f'<div class="as-sec"><h2>{title}</h2>'
                    '<div class="as-scroll"><table class="as-t">'
                    '<tr><th>Assumption</th><th class="num">Value</th><th>Unit</th>'
                    '<th>Key</th></tr>'
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
            mark = ('<span class="as-moved" title="moved by the sensitivity">*</span>'
                    if path in moved else "")
            cells.append(f'<td class="v num">{fmt(read(assumptions, path))}{mark}</td>')
        rows.append(f'<tr><td>{label}</td><td class="u">{unit}</td>'
                    + "".join(cells) + '</tr>')
    html.append(
        '<div class="as-sec"><h2>Process steps</h2>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Assumption</th>'
        f'<th>Unit</th>{head}</tr>' + "".join(rows) + '</table></div></div>')

    # ---- the EAF's charge states ----------------------------------------
    rows = []
    for state, label in CHARGE_STATES:
        path = f"eaf.charge.{state}.el_mwh_per_t"
        shown.add(path)
        rows.append(f'<tr><td>{label}</td>'
                    f'<td class="v num">{fmt(read(assumptions, path))}</td>'
                    f'<td class="u">MWh/t</td></tr>')
    html.append(
        '<div class="as-sec"><h2>EAF electricity by charge state</h2>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Charge</th>'
        '<th class="num">Electricity</th><th>Unit</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- freight legs ---------------------------------------------------
    rows = []
    for origin, legs in read(assumptions, "transport.distance_km").items():
        for mode, km in legs.items():
            shown.add(f"transport.distance_km.{origin}.{mode}")
        drawn = " + ".join(f"{fmt(km)} km {mode}" for mode, km in legs.items())
        rows.append(f'<tr><td class="k">{origin}</td><td class="v">{drawn}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>Freight distances to the destination</h2>'
        '<div class="as-scroll"><table class="as-t">'
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
        '<div class="as-scroll" style="display:flex;gap:36px;align-items:flex-start;'
        'flex-wrap:wrap;">'
        '<table class="as-t" style="width:auto"><tr><th>Assumption</th>'
        '<th class="num">Value</th><th>Unit</th></tr>'
        f'<tr><td>Basis</td><td class="v num">{basis}</td><td class="u"></td></tr>'
        f'<tr><td>Gas supply chain</td><td class="v num">{gas_upstream}</td>'
        '<td class="u">t CO2e/MWh LHV</td></tr></table>'
        + "".join(columns) +
        '<table class="as-t" style="width:auto"><tr><th>Freight</th>'
        f'<th class="num">Factor</th><th>Unit</th></tr>{"".join(freight)}</table>'
        '</div></div>')

    # ---- geographies ----------------------------------------------------
    solved_areas = set()
    for report in sorted(RESULTS.glob("*/.report_*_diag.csv")):
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
        '<div class="as-sec"><h2>Geographies and price sources</h2>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Area</th>'
        '<th>Country</th><th>Market</th><th>Price series</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- where the scenarios differ -------------------------------------
    rows = []
    for path, by_scenario in sorted(moved.items()):
        by_value = {}
        for scenario, value in by_scenario.items():
            by_value.setdefault(value, []).append(scenario)
        values = " · ".join(fmt(value) for value in sorted(by_value))
        rows.append(f'<tr><td class="k">{path}</td>'
                    f'<td class="v num">{fmt(read(assumptions, path))}</td>'
                    f'<td class="n">{values}</td></tr>')
    html.append(
        '<div class="as-sec"><h2>Sensitivity</h2>'
        '<div class="as-scroll"><table class="as-t"><tr><th>Assumption</th>'
        '<th class="num">Base</th><th>Moved to</th></tr>'
        + "".join(rows) + '</table></div></div>')

    # ---- coverage --------------------------------------------------------
    # The page does not list what it leaves out, but the build still checks that
    # something decided to leave it out: anything neither rendered nor named in
    # SKIPPED stops the build, so an assumption added later cannot go unread.
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
