#!/usr/bin/env python3
"""Build the cost-taxonomy hub tab: two ways of cutting the same levelised cost.

Writes `results/html/cost_taxonomies.html` — a body-only, theme-aware, self-contained
page that puts the circulated chart mockup's *category split* next to the split
the pipeline already reports, both populated from the same results, for any
scenario in the matrix.

Five panels, all with the model's five routes as the bars:

  * LCOS, LCOE and LCOH, each shown both ways side by side;
  * built capacity, for reference;
  * the same LCOS cut as finely as the model allows — up to 40 leaves, each hover
    naming the inputs behind that number.

The page renders client-side from the shared dashboard payload (the same one
`build_dashboard_v2.py` uses), so the scenario matrix can move geography, weather
year, grid mode, capacity-factor basis and sensitivity variant without a rebuild.
Every number on it is a report column: `compile_report` cuts the cost groups off
the solved networks, `cost_taxonomy.py` reads them back and keys them the way the
bands are keyed, and this module only assembles them into the payload and fills
the template. So a report compiled before those columns existed has nothing to
plot, and `attach` says so rather than failing on a missing column.

Run it directly — it is not a pipeline rule.
"""
import json
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[3]        # workflow/scripts/viz/ -> repo root
sys.path.insert(0, str(Path(__file__).parent))    # sibling build_dashboard, cost_taxonomy
sys.path.insert(0, str(REPO / "workflow"))        # common.*, scripts.*

import cost_taxonomy                                                      # noqa: E402
from build_dashboard import HTML_DIR, _axes, build_html, scenario_files      # noqa: E402
from common._constants import H2_LHV_KWH_PER_KG                           # noqa: E402
from common._report_schema import field_stem, read_report                 # noqa: E402

OUT_PATH = HTML_DIR / "cost_taxonomies.html"
TEMPLATE_HTML = Path(__file__).with_name("cost_taxonomies_template.html")
CONFIG_DIR = REPO / "config"

# ---- band specs, published to the client -------------------------------------
# (key, label, colour) bottom -> top. The two taxonomies share a colour per role,
# so a band means the same thing whichever side of the pair you read.
C_ORE, C_CAPEX, C_OM = "#E2B681", "#33434D", "#BDCCD9"
C_HYDROGEN, C_H2_STORE, C_ELEC = "#91C096", "#70D2F0", "#0A5680"
C_ELEC_EAF, C_GAS, C_STORE = "#0293D2", "#525F6A", "#D75674"
C_GRID_CONN, C_GRID_NRG, C_TRANSM = "#71828F", "#B7C1C8", "#83D1DD"

# The electricity is four bands rather than two and a remainder: making the
# hydrogen, making the iron, melting it, and the handling and losses around
# them. `rest of plant` used to carry the third of those, which on moe-eaf and
# ew-eaf is the cell that makes the iron and most of the route's electricity.
PURPOSE_BANDS = [
    ["ore",             "Ore & consumables",                C_ORE],
    ["capex",           "CAPEX (process plant)",            C_CAPEX],
    ["fixed_om",        "O&M (process plant)",              C_OM],
    ["hydrogen",        "Hydrogen (all-in)",                C_HYDROGEN],
    ["reduction",       "Electricity — making the iron",    C_ELEC],
    ["melt",            "Electricity — EAF melt",           C_ELEC_EAF],
    ["handling_losses", "Electricity — handling & losses",  C_GRID_NRG],
    ["gas",             "Gas + CO₂ (fossil routes)",        C_GAS],
    ["store",           "Storage (iron/steel)",             C_STORE],
    ["transport",       "Freight",                          "#BDCCD9"],
]
LCOE_SPLIT_BANDS = [
    ["res_capex",         "Renewables — capital",     C_ELEC],
    ["res_fom",           "Renewables — fixed O&M",   "#4B93BF"],
    ["storage",           "Battery / storage",        C_STORE],
    ["grid_connection",   "Grid connection",          C_GRID_CONN],
    ["grid_energy",       "Grid energy",              C_GRID_NRG],
    ["transmission",      "Transmission (HVDC)",      C_TRANSM],
    # Only an export route has one; on every other it is zero and draws nothing.
    ["destination_power", "Destination power",        "#0293D2"],
]
LCOE_GROUP_BANDS = [
    ["renewables",        "Renewables (capex+opex)",  C_ELEC],
    ["storage",           "Battery / storage",        C_STORE],
    ["grid_connection",   "Grid connection",          C_GRID_CONN],
    ["grid_energy",       "Grid energy",              C_GRID_NRG],
    ["transmission",      "Transmission (HVDC)",      C_TRANSM],
    ["destination_power", "Destination power",        "#0293D2"],
]
LCOH_SPLIT_BANDS = [
    ["electrolyser_capex", "Electrolyser — capital",   C_CAPEX],
    ["electrolyser_fom",   "Electrolyser — fixed O&M", C_OM],
    ["electrolyser_water", "Water & variable opex",    "#B4D4B8"],
    ["storage",            "H₂ storage (buffer)",      C_H2_STORE],
    ["electricity",        "Electricity (@ LCOE)",     C_ELEC],
]
LCOH_GROUP_BANDS = [
    ["electrolyser", "Electrolyser (capex+opex)", C_CAPEX],
    ["storage",      "H₂ storage (buffer)",       C_H2_STORE],
    ["electricity",  "Electricity (@ LCOE)",      C_ELEC],
]


def _assumptions(scenario: str) -> dict:
    """Base assumptions with the scenario's overlay applied on top.

    An overlay is thin — it carries only the keys it bumps — so the physical
    coefficients come from the base file. Reading the overlay alone silently
    yields nothing.
    """
    merged = yaml.safe_load((CONFIG_DIR / "assumptions.yaml").read_text())
    overlay_path = CONFIG_DIR / f"assumptions_{scenario}.yaml"
    if overlay_path.exists():
        for key, value in (yaml.safe_load(overlay_path.read_text()) or {}).items():
            if isinstance(value, dict) and isinstance(merged.get(key), dict):
                merged[key].update(value)
            else:
                merged[key] = value
    return merged


def _round_map(values: dict, places: int) -> dict:
    """Drop the negligible and round the rest, so the payload carries no noise."""
    floor = 0.5 * 10 ** -places
    return {key: round(value, places) for key, value in values.items()
            if value is not None and abs(value) >= floor}


def attach(payload: dict, cases: dict) -> None:
    """Add the taxonomy split to every record in `cases`, and the specs alongside.

    Walks the reports a second time — `build_payload` does not keep the raw rows —
    and matches each row onto its record with the same scenario parser, so a record
    can never be given another scenario's split.
    """
    payload["leaf_groups"] = [list(group) for group in cost_taxonomy.GROUPS]
    # The plants, by the stem their leaves are keyed on, so a hover over the
    # capital or the upkeep band can name each one and reach its own half of
    # that plant's annual cost.
    payload["process_plants"] = [[field_stem(link), label]
                                 for link, label in cost_taxonomy.PROCESS_PLANTS]
    payload["purpose_bands"] = PURPOSE_BANDS
    payload["lcoe_split_bands"] = LCOE_SPLIT_BANDS
    payload["lcoe_group_bands"] = LCOE_GROUP_BANDS
    payload["lcoh_split_bands"] = LCOH_SPLIT_BANDS
    payload["lcoh_group_bands"] = LCOH_GROUP_BANDS
    payload["h2_lhv_kwh_per_kg"] = H2_LHV_KWH_PER_KG

    # The config quotes a hover shows are per scenario, not global: the generated
    # overlays move the gas price for 120 scenarios and the H2-buffer capex for 109
    # (that is the salt-cavern sensitivity), among others. Publishing one spec built
    # from a single scenario would print the base quote next to salt-cavern numbers.
    # So each scenario gets its own spec, deduplicated — there are only a handful of
    # distinct ones, and a record just points at the one it uses.
    specs: dict[str, int] = {}
    payload["leaf_specs"] = []

    def spec_index(assumptions: dict) -> int:
        entry = cost_taxonomy.spec(assumptions)
        fingerprint = json.dumps(entry, sort_keys=True)
        if fingerprint not in specs:
            specs[fingerprint] = len(payload["leaf_specs"])
            payload["leaf_specs"].append(entry)
        return specs[fingerprint]

    attached = missing = 0
    # The same reports build_dashboard read, walked again to attach the leaf
    # split to the records it already built. A row finds its record through the
    # axes it declares, not through a name that has to be taken apart.
    for scenario, report in scenario_files().items():
        # The overlay is named by the scenario the file is for, which is not
        # what the record is filed under: `standard-grid` and `standard-islanded`
        # share the `base` pill but keep their own names on disk. Only `spec()`
        # needs it now — the split itself is read off the report, which
        # compile_report cut against this same merged overlay.
        assumptions = _assumptions(scenario)
        runs = read_report(report)
        # The split is the report's now, so a report compiled before it existed
        # has nothing for this page to plot. Say that, rather than failing on a
        # missing column three frames further in.
        probe = f"cost_{cost_taxonomy.LEAF_COSTS[0]}_eur_per_t"
        if probe not in runs.columns:
            raise SystemExit(
                f"{report.name} carries no cost leaves ({probe} is not in it), so "
                f"there is nothing to break down. Recompile the reports — "
                f"`snakemake results/report_{scenario}.csv --force` — which needs "
                f"the scenario's solved networks under results/{scenario}/."
            )
        for _, row in runs.iterrows():
            axes = _axes(row)
            project = f"{axes['geo']}-{axes['year']}-{axes['grid']}"
            record = (cases.get(project, {})
                           .get(axes["route"], {})
                           .get(axes["variant"], {})
                           .get(axes["cf"]))
            if record is None:
                continue
            bands = cost_taxonomy.purpose_bands(row)
            if not bands:
                missing += 1
                continue
            carriers = cost_taxonomy.carrier_splits(row)

            record["spec"] = spec_index(assumptions)
            record["purpose"] = _round_map(bands, 2)
            record["lcoe_split"] = _round_map(carriers["lcoe"], 2)
            record["lcoh_split"] = _round_map(carriers["lcoh"], 2)
            record["leaves"] = _round_map(cost_taxonomy.leaf_costs(row), 2)
            record["groups"] = _round_map(cost_taxonomy.parent_costs(row), 2)
            record["leaf_inputs"] = {key: [list(pair) for pair in lines]
                                     for key, lines in cost_taxonomy.leaf_inputs(row).items()
                                     if key in record["leaves"]}
            attached += 1

    print(f"  attached the taxonomy split to {attached} records"
          + (f" ({missing} had no steel load)" if missing else "")
          + f", over {len(payload['leaf_specs'])} distinct assumption sets")


def main() -> None:
    html, cases, geos = build_html(TEMPLATE_HTML, augment=attach)
    # The single-taxonomy build fills this in; here it leaves the class list clean.
    html = html.replace(" /*PURPOSE_ONLY_CLASS*/", "")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(html, encoding="utf-8")
    print(f"wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1e6:.2f} MB, "
          f"{html.count(chr(0))} NULs) — {len(cases)} projects, {len(geos)} geos, "
          "body-only hub fragment")


if __name__ == "__main__":
    main()
