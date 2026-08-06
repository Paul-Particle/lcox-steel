#!/usr/bin/env python3
"""Two cost-breakdown taxonomies, both drawn from the model's own results.

Builds `results/cost_breakdown_mockups.html` — a body-only, theme-aware,
self-contained Plotly page that puts the circulated chart mockup's *category
split* next to the split the pipeline already reports, with **both** populated
from the same DE 2023 results the scenario-comparison tab reads.

  * ALTERNATIVE — the circulated taxonomy: one "Hydrogen (all-in)" block, process
    plant cut into annualised capex and fixed O&M, and electricity divided into
    an EAF-melt share and everything else.
  * DASHBOARD — the pipeline's own cost groups, unchanged; what the rest of the
    dashboard already shows.

Both stacks close exactly on LCOS, so the two differ only in how the same total
is cut. That is the point: the alternative's hydrogen block turns out to swallow
most of the H2 route's electricity, which is visible only once real numbers are
in it.

Reconstructing the alternative needs two things the report does not name directly.

Electricity is partitioned three ways — two parts attributed and the third a
residual, so the stack cannot drift from LCOS:

  * H2 share       = LCOH's electricity component x the H2 actually produced
  * EAF melt share = the EAF's own MWh/t at that scenario's LCOE
  * rest of plant  = electricity system total - the two above

Process plant is cut into annualised capex and fixed O&M from the per-plant
`plant_*_eur_per_t` columns: a plant's PyPSA capital_cost is
`annuity x capex_per_t_per_year + opex_per_t_per_year`, so the ratio is a config
constant and the split is exact (see `_fixed_om_shares`).

Segment hovers carry the same information as the scenario-comparison tab's cost
chart — a headline with the segment's share of the levelised total, then the
per-plant, per-technology and per-carrier detail lines for that bar.

LCOE and LCOH get a single chart each: there the two taxonomies coincide apart
from sub-splits the reports fold into their parent line, so a second identical
panel would say nothing.

Output is a hub fragment: body-only, its own <style>/<script>, plotly.js and the
brand font inlined. `build_dashboard_hub.py` embeds it as the "Cost breakdown
options" tab. Nothing here is a pipeline rule — run it directly.
"""
import copy
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import yaml

REPO = Path(__file__).resolve().parents[3]        # workflow/scripts/viz/ -> repo root
sys.path.insert(0, str(Path(__file__).parent))    # sibling build_dashboard
sys.path.insert(0, str(REPO / "workflow"))        # common.*, scripts.*

from build_dashboard import COST_GROUPS, RESULTS, ROUTE_LABEL, font_css  # noqa: E402
from common._constants import H2_LHV_KWH_PER_KG                          # noqa: E402
from scripts.h2_dri._helpers import annuity_factor                       # noqa: E402
from scripts.viz.style import PLOTLY_CONFIG, fca_template                # noqa: E402

OUT_PATH = RESULTS / "cost_breakdown_mockups.html"
CONFIG_DIR = REPO / "config"

# The case both taxonomies are drawn from — the same reports the scenario tab reads.
PROJECT = "DE-2023-grid"
PROJECT_OFFGRID = "DE-2023-nogrid"
LCOS_SCENARIOS = ["h2-dri-eaf-avg", "moe-avg", "ew-avg", "ng-dri-eaf", "mix-dri-eaf-avg"]

# Charts sit on a constant white "print" surface in both themes — the same
# convention the v2 dashboard uses, so a screenshot is publication-ready either
# way and no re-render is needed when the viewer toggles the theme.
PRINT = dict(paper="#FFFFFF", ink="#1B242B", muted="#525F6A",
             grid="#E4E9ED", border="#CDD5DB")
CONFIG = {**PLOTLY_CONFIG, "responsive": True}
PLOT_HEIGHT = 440

# fca_template pins width=960/height=540 for static PNG export, and layout
# properties fall back to the template — so `fig.layout.width = None` cannot clear
# it and `to_html` would emit a fixed 960px container. These figures are sized by
# the page grid instead, so clear the pinned size on a copy of the template.
responsive_template = copy.deepcopy(fca_template)
responsive_template.layout.autosize = True
responsive_template.layout.width = None
responsive_template.layout.height = None

# ---- One colour per cost role, shared by both taxonomies ------------------
# The dashboard taxonomy's LCOS colours come from build_dashboard.COST_GROUPS;
# these mirror the same palette so equivalent roles read alike across the pair.
C_ORE = "#E2B681"        # sand yellow    — ore & consumables
C_CAPEX = "#33434D"      # blue black     — plant capex
C_OM = "#BDCCD9"         # light blue gray — plant fixed O&M: the same family as the
#                          capex band above it, and the one blue-gray the dashboard
#                          taxonomy does not use, so the paired legends never show
#                          one swatch twice.
C_HYDROGEN = "#91C096"   # green          — hydrogen / electrolyser
C_H2_STORE = "#70D2F0"   # light blue     — H2 buffer
C_ELEC = "#0A5680"       # fca blue       — electricity / renewables
C_ELEC_EAF = "#0293D2"   # highlight blue — the EAF-melt slice of the electricity
C_GAS = "#525F6A"        # very dark gray — natural gas + CO2
C_STORE = "#D75674"      # magenta red    — storage (iron/steel, battery)
C_GRID_CONN = "#71828F"  # dark gray      — grid connection
C_GRID_NRG = "#B7C1C8"   # gray           — grid energy
C_TRANSM = "#83D1DD"     # turquois       — HVDC transmission

# Report plant key -> (assumptions block, block whose lifetime annuitises it,
# label). Ladle metallurgy is folded into the MOE link by `_add_moe_link` and
# annuitised over the MOE lifetime, so it borrows MOE's lifetime here too.
PROCESS_PLANTS = [
    ("dri",            "dri",            "dri",            "H2-DRI shaft"),
    ("dri_ng",         "dri_ng",         "dri_ng",         "NG-DRI shaft"),
    ("moe",            "moe",            "moe",            "MOE"),
    ("ladle",          "ladle",          "moe",            "Ladle metallurgy"),
    ("electrowinning", "electrowinning", "electrowinning", "Electrowinning"),
    ("eaf",            "eaf",            "eaf",            "EAF"),
]


# ---- Reading the pipeline reports ----------------------------------------

def _report(project: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"report_{project}.csv").set_index("scenario")


def _value(row: pd.Series, column: str):
    """A report column as a float, or None when this scenario does not have it."""
    if column not in row.index or pd.isna(row[column]):
        return None
    return float(row[column])


def _assumptions(project: str, scenario: str) -> dict:
    """Base assumptions with the scenario's generated overlay applied on top.

    The generated `assumptions_{project}_{scenario}.yaml` files are thin overlays —
    typically just `route:` — so the physical coefficients come from the base file.
    Reading the overlay alone silently yields nothing.
    """
    merged = yaml.safe_load((CONFIG_DIR / "assumptions.yaml").read_text())
    overlay_path = CONFIG_DIR / f"assumptions_{project}_{scenario}.yaml"
    if overlay_path.exists():
        for key, value in (yaml.safe_load(overlay_path.read_text()) or {}).items():
            if isinstance(value, dict) and isinstance(merged.get(key), dict):
                merged[key].update(value)
            else:
                merged[key] = value
    return merged


def _eur_per_t(row: pd.Series, group: str, steel_t: float) -> float:
    """One cost group in €/t steel; groups a scenario does not have read as zero."""
    col = f"cost_{group}_meur"
    if col not in row.index or pd.isna(row[col]):
        return 0.0
    return float(row[col]) * 1e6 / steel_t


# ---- Splitting what the report bundles -----------------------------------

def _fixed_om_shares(assumptions: dict) -> dict:
    """Fixed O&M's share of each process plant's annual cost.

    `build_network._process_capital_cost` sets a plant's PyPSA capital_cost to
    `annuity x capex_per_t_per_year + opex_per_t_per_year`, so the reported
    `plant_*_eur_per_t` figures bundle annualised capital with fixed O&M in a ratio
    fixed entirely by the config. The split is therefore exact and needs no new
    reported quantity — at ~3 % of capex per year it comes out near a quarter of
    each plant's annual cost.
    """
    wacc = assumptions["finance"]["default_wacc"]
    shares = {}
    for key, block, lifetime_block, _ in PROCESS_PLANTS:
        plant = assumptions[block]
        annual_capex = (annuity_factor(wacc, assumptions[lifetime_block]["lifetime_years"])
                        * plant["capex_per_t_per_year_eur"])
        fixed_om = plant["opex_per_t_per_year_eur"]
        shares[key] = fixed_om / (annual_capex + fixed_om)
    return shares


def _split_process(row: pd.Series, scenario: str, steel_t: float) -> dict:
    """Process-plant cost (€/t steel) cut into annualised capex and fixed O&M.

    Per-plant figures come from the report's `plant_*_eur_per_t` columns, which sum
    to the `process` cost group; each is cut by that plant's config O&M share.
    """
    shares = _fixed_om_shares(_assumptions(PROJECT, scenario))
    plants = {}
    for key, _, _, label in PROCESS_PLANTS:
        plant_cost = _value(row, f"plant_{key}_eur_per_t")
        if plant_cost:
            plants[key] = dict(label=label, total=plant_cost,
                               fixed_om=plant_cost * shares[key],
                               capex=plant_cost * (1.0 - shares[key]))

    total = _eur_per_t(row, "process", steel_t)
    drift = sum(plant["total"] for plant in plants.values()) - total
    if abs(drift) > 0.01:
        raise ValueError(
            f"{scenario}: the per-plant costs miss the reported process group by "
            f"{drift:+.3f} €/t, so the capex/O&M split cannot be trusted."
        )
    return dict(total=total, plants=plants,
                capex=sum(plant["capex"] for plant in plants.values()),
                fixed_om=sum(plant["fixed_om"] for plant in plants.values()))


def _split_electricity(row: pd.Series, scenario: str, steel_t: float) -> dict:
    """Partition the electricity-system cost (€/t steel) the way the alternative wants it.

    Returns the H2 share, the EAF melt share and the remainder. The remainder is
    a residual, so the three always sum back to the electricity-system total and
    the alternative's stack closes on LCOS.
    """
    total = sum(_eur_per_t(row, group, steel_t)
                for group in ("res", "grid", "battery", "transmission"))
    lcoe = _value(row, "lcoe_eur_per_mwh") or 0.0

    h2_kt = _value(row, "h2_produced_kt") or 0.0
    lcoh_electricity = _value(row, "lcoh_electricity_eur_per_mwh_lhv") or 0.0
    h2_mwh_lhv = h2_kt * 1e6 * H2_LHV_KWH_PER_KG / 1000.0
    hydrogen = lcoh_electricity * h2_mwh_lhv / steel_t

    # The EAF melts only on routes that build one; MOE pours liquid steel directly
    # and its assumptions file carries no `eaf` block at all.
    eaf_built = (_value(row, "eaf_t_per_h_opt") or 0.0) > 0
    eaf_mwh_per_t = _assumptions(PROJECT, scenario)["eaf"]["el_mwh_per_t"] if eaf_built else 0.0
    eaf = eaf_mwh_per_t * lcoe

    rest = total - hydrogen - eaf
    if rest < -0.01:
        raise ValueError(
            f"{scenario}: attributing {hydrogen:.1f} €/t to hydrogen and {eaf:.1f} €/t to the "
            f"EAF over-claims an electricity system of {total:.1f} €/t. The EAF melt share is "
            "reconstructed from el_mwh_per_t x LCOE, so it is an estimate — revisit it before "
            "trusting this split."
        )
    return dict(hydrogen=hydrogen, eaf=eaf, rest=max(rest, 0.0), total=total,
                lcoe=lcoe, eaf_mwh_per_t=eaf_mwh_per_t)


# ---- Hover text ----------------------------------------------------------
# The same shape the scenario-comparison tab uses: a headline naming the segment,
# its value and its share of the levelised total, then indented detail lines built
# from that bar's own components — a bar missing a component omits the line rather
# than printing "n/a".

def _amount(value, unit: str, decimals: int = 0) -> str:
    """"41 €/t" — or "n/a" where a scenario has no such component."""
    return "n/a" if value is None else f"{value:,.{decimals}f} {unit}"


def _amount_of(value, total, basis: str, unit: str, decimals: int = 0, shown=None) -> str:
    """"41 €/t (6% of LCOS)".

    `shown` overrides the figure printed while `value` still drives the share, so a
    technology's own unit cost can carry its contribution's share of the system
    total — the convention the scenario tab uses for renewables.
    """
    text = _amount(shown if shown is not None else value, unit, decimals)
    if total and value is not None:
        text += f" ({round(value / total * 100)}% of {basis})"
    return text


def _lcoe_tech_lines(row: pd.Series, *, headline: bool) -> list:
    """Renewables detail: each technology's own LCOE, with its share of system LCOE."""
    lcoe = _value(row, "lcoe_eur_per_mwh")
    lines = []
    for own_col, part_col, label in (
        ("lcoe_renewables_own_eur_per_mwh", "lcoe_renewables_eur_per_mwh", "renewables:"),
        ("lcoe_solar_own_eur_per_mwh", "lcoe_solar_eur_per_mwh", "· solar"),
        ("lcoe_wind_onshore_own_eur_per_mwh", "lcoe_wind_onshore_eur_per_mwh", "· wind onshore"),
        ("lcoe_wind_offshore_own_eur_per_mwh", "lcoe_wind_offshore_eur_per_mwh",
         "· wind offshore"),
    ):
        if label == "renewables:" and not headline:
            continue
        part = _value(row, part_col)
        if part is None:
            continue
        lines.append((label, _amount_of(part, lcoe, "LCOE", "€/MWh", 1,
                                        shown=_value(row, own_col))))
    return lines


def _grid_lines(row: pd.Series) -> list:
    """Delivered grid price split into market price, volumetric fee and connection."""
    lcoe = _value(row, "lcoe_eur_per_mwh")
    price = _value(row, "grid_price_eur_per_mwh")
    fee = _value(row, "grid_fee_eur_per_mwh")
    connection = _value(row, "grid_connection_eur_per_mwh_imported")
    if all(part is None for part in (price, fee, connection)):
        return []
    delivered = sum(part for part in (price, fee, connection) if part is not None)
    contribution = sum(_value(row, col) or 0.0 for col in
                       ("lcoe_grid_connection_eur_per_mwh", "lcoe_grid_energy_eur_per_mwh"))
    lines = [("grid electricity:", _amount_of(contribution or None, lcoe, "LCOE",
                                              "€/MWh", 1, shown=delivered))]
    for part, label in ((price, "· market price"), (fee, "· volumetric fee"),
                        (connection, "· connection")):
        if part is not None:
            lines.append((label, _amount(part, "€/MWh", 1)))
    return lines


def _lcoh_lines(row: pd.Series, *, headline: bool) -> list:
    """Hydrogen detail: LCOH and the parts that sum to it, all €/MWh LHV."""
    lcoh = _value(row, "lcoh_eur_per_mwh_lhv")
    lines = []
    if headline and lcoh is not None:
        lines.append(("hydrogen LCOH:", _amount(lcoh, "€/MWh LHV", 1)))
    for col, label in (("lcoh_electrolyser_eur_per_mwh_lhv", "· electrolyser"),
                       ("lcoh_electricity_eur_per_mwh_lhv", "· electricity"),
                       ("lcoh_h2_storage_eur_per_mwh_lhv", "· H₂ storage")):
        part = _value(row, col)
        if part is not None:
            lines.append((label, _amount_of(part, lcoh, "LCOH", "€/MWh LHV", 1)))
    return lines


def _ore_lines(row: pd.Series, lcos: float) -> list:
    lines = []
    for col, label in (("ore_eur_per_t_steel", "· iron ore"),
                       ("consumables_eur_per_t_steel",
                        "· EAF consumables (electrodes, fluxes, alloys)")):
        part = _value(row, col)
        if part is not None:
            lines.append((label, _amount_of(part, lcos, "LCOS", "€/t")))
    return lines


def _plant_lines(process: dict, key: str, lcos: float) -> list:
    """Per-plant contribution to one of the process bands (`capex` or `fixed_om`)."""
    return [(f"· {plant['label']}", _amount_of(plant[key], lcos, "LCOS", "€/t"))
            for plant in process["plants"].values() if plant[key] > 0]


def _gas_lines(row: pd.Series, scenario: str) -> list:
    if not (_value(row, "ng_gwh_lhv") or 0.0) > 0:
        return []
    gas = _assumptions(PROJECT, scenario)["natural_gas"]
    return [("gas price:", _amount(gas["price_eur_per_mwh"], "€/MWh LHV"))]


STORES = [("iron_store", "iron_store_hours_steel", "· iron stockpile"),
          ("steel_store", "steel_store_hours_steel", "· steel inventory")]


def _store_lines(row: pd.Series, steel_t: float, groups=("iron_store", "steel_store")) -> list:
    """The solid stores, with their size — they cost too little to round to a €/t band."""
    lines = []
    for group, hours_col, label in STORES:
        if group not in groups:
            continue
        cost = _eur_per_t(row, group, steel_t)
        if cost <= 0:
            continue
        hours = _value(row, hours_col)
        text = _amount(cost, "€/t", 2)
        if hours is not None:
            text += f" for {hours:,.0f} h of demand"
        lines.append((label, text))
    return lines


def _levelised_totals_lines(row: pd.Series) -> list:
    """The cross-links every bar carries: its LCOE and LCOH, as on the scenario tab."""
    return [("electricity (LCOE):", _amount(_value(row, "lcoe_eur_per_mwh"), "€/MWh", 1)),
            ("hydrogen (LCOH):", _amount(_value(row, "lcoh_eur_per_mwh_lhv"),
                                         "€/MWh LHV", 1))]


# ---- The two LCOS taxonomies ---------------------------------------------

def alternative_lcos():
    """The circulated taxonomy, reconstructed from the model's own cost groups."""
    dataframe = _report(PROJECT)
    rows = dataframe.loc[LCOS_SCENARIOS]
    bars = [ROUTE_LABEL[scenario.removesuffix("-avg")] for scenario in LCOS_SCENARIOS]

    keys = ("ore", "capex", "fixed_om", "hydrogen", "eaf", "rest", "gas", "store")
    values = {key: [] for key in keys}
    lines = {key: [] for key in keys}
    totals_lines = []

    for scenario in LCOS_SCENARIOS:
        row = dataframe.loc[scenario]
        steel_t = float(row["steel_produced_mt"]) * 1e6
        lcos = float(row["lcos_eur_per_t"])
        electricity = _split_electricity(row, scenario, steel_t)
        process = _split_process(row, scenario, steel_t)
        electrolyser = _eur_per_t(row, "electrolyser", steel_t)
        h2_buffer = _eur_per_t(row, "h2_buffer", steel_t)

        values["ore"].append(_eur_per_t(row, "ore_consumables", steel_t))
        lines["ore"].append(_ore_lines(row, lcos))

        values["capex"].append(process["capex"])
        lines["capex"].append(_plant_lines(process, "capex", lcos))
        values["fixed_om"].append(process["fixed_om"])
        lines["fixed_om"].append(_plant_lines(process, "fixed_om", lcos))

        values["hydrogen"].append(electrolyser + h2_buffer + electricity["hydrogen"])
        hydrogen_detail = _lcoh_lines(row, headline=True)[:1]
        for part, label in ((electrolyser, "· electrolyser plant"),
                            (h2_buffer, "· H₂ buffer"),
                            (electricity["hydrogen"], "· electricity to make the H₂")):
            if part > 0:
                hydrogen_detail.append((label, _amount_of(part, lcos, "LCOS", "€/t")))
        lines["hydrogen"].append(hydrogen_detail)

        values["eaf"].append(electricity["eaf"])
        lines["eaf"].append(
            [("electricity LCOE:", _amount(electricity["lcoe"], "€/MWh", 1)),
             ("· furnace demand", f"{electricity['eaf_mwh_per_t']:.2f} MWh/t steel")]
            if electricity["eaf"] > 0 else [])

        values["rest"].append(electricity["rest"])
        lines["rest"].append([
            ("electricity LCOE:", _amount(electricity["lcoe"], "€/MWh", 1)),
            ("· electricity system, all in",
             _amount_of(electricity["total"], lcos, "LCOS", "€/t")),
            ("· less the H₂ and EAF melt shares",
             _amount(electricity["hydrogen"] + electricity["eaf"], "€/t")),
        ])

        values["gas"].append(_eur_per_t(row, "gas", steel_t))
        lines["gas"].append(_gas_lines(row, scenario))

        values["store"].append(_eur_per_t(row, "iron_store", steel_t)
                               + _eur_per_t(row, "steel_store", steel_t))
        lines["store"].append(_store_lines(row, steel_t))

        totals_lines.append(_levelised_totals_lines(row))

    series = [
        ("Ore & consumables",           C_ORE,      "ore"),
        ("CAPEX (process plant)",       C_CAPEX,    "capex"),
        ("O&M (process plant)",         C_OM,       "fixed_om"),
        ("Hydrogen (all-in)",           C_HYDROGEN, "hydrogen"),
        ("Electricity — EAF melt",      C_ELEC_EAF, "eaf"),
        ("Electricity — rest of plant", C_ELEC,     "rest"),
        ("Gas + CO₂ (fossil routes)",   C_GAS,      "gas"),
        ("Storage (iron/steel)",        C_STORE,    "store"),
    ]
    return (bars,
            [(label, colour, values[key], lines[key]) for label, colour, key in series],
            rows["lcos_eur_per_t"].tolist(), totals_lines)


def dashboard_lcos():
    """The pipeline's own LCOS cost groups, unchanged, for the same scenarios."""
    dataframe = _report(PROJECT)
    rows = dataframe.loc[LCOS_SCENARIOS]
    bars = [ROUTE_LABEL[scenario.removesuffix("-avg")] for scenario in LCOS_SCENARIOS]

    detail = {group: [] for group, _, _ in COST_GROUPS}
    totals_lines = []
    for scenario in LCOS_SCENARIOS:
        row = dataframe.loc[scenario]
        steel_t = float(row["steel_produced_mt"]) * 1e6
        lcos = float(row["lcos_eur_per_t"])
        process = _split_process(row, scenario, steel_t)
        # The process group bundles annualised capital with fixed O&M; this taxonomy
        # keeps the band whole, so the split rides in the hover instead.
        process_lines = [
            ("· annualised capex", _amount_of(process["capex"], lcos, "LCOS", "€/t")),
            ("· fixed O&M", _amount_of(process["fixed_om"], lcos, "LCOS", "€/t")),
        ] + [(f"· {plant['label']}", _amount_of(plant["total"], lcos, "LCOS", "€/t"))
             for plant in process["plants"].values()]

        for group in detail:
            if group == "process":
                detail[group].append(process_lines)
            elif group == "ore_consumables":
                detail[group].append(_ore_lines(row, lcos))
            elif group == "gas":
                detail[group].append(_gas_lines(row, scenario))
            elif group in ("electrolyser", "h2_buffer"):
                detail[group].append(_lcoh_lines(row, headline=True))
            elif group == "res":
                detail[group].append(_lcoe_tech_lines(row, headline=True))
            elif group == "battery":
                detail[group].append(
                    [("storage:", _amount_of(_value(row, "lcoe_storage_eur_per_mwh"),
                                             _value(row, "lcoe_eur_per_mwh"),
                                             "LCOE", "€/MWh", 1))])
            elif group == "transmission":
                detail[group].append(
                    [("transmission:", _amount_of(_value(row, "lcoe_transmission_eur_per_mwh"),
                                                  _value(row, "lcoe_eur_per_mwh"),
                                                  "LCOE", "€/MWh", 1))])
            elif group == "grid":
                detail[group].append(_grid_lines(row))
            elif group in ("iron_store", "steel_store"):
                detail[group].append(_store_lines(row, steel_t, (group,)))
            else:
                detail[group].append([])
        totals_lines.append(_levelised_totals_lines(row))

    steel_t = rows["steel_produced_mt"] * 1e6
    series = []
    for group, label, colour in COST_GROUPS:
        column = f"cost_{group}_meur"
        if column not in rows.columns:
            continue
        values = (rows[column].fillna(0.0) * 1e6 / steel_t).tolist()
        if all(value <= 0.005 for value in values):
            continue
        series.append((label, colour, values, detail[group]))
    return bars, series, rows["lcos_eur_per_t"].tolist(), totals_lines


# ---- The carriers underneath ---------------------------------------------

# The three supply cases the LCOE and LCOH charts compare, all DE 2023,
# H2-DRI-EAF. `bestsite-p95` is the 95th-percentile area-weighted capacity-factor
# cell, so it is a better site than the mean and comes out cheaper.
SUPPLY_CASES = [
    (PROJECT,         "h2-dri-eaf-avg", "Grid-connected"),
    (PROJECT_OFFGRID, "h2-dri-eaf-p95", "Off-grid · p95 site"),
    (PROJECT_OFFGRID, "h2-dri-eaf-avg", "Off-grid · mean site"),
]
LCOE_PARTS = [
    ("lcoe_renewables_eur_per_mwh",      "Renewables (capex+opex)", C_ELEC),
    ("lcoe_storage_eur_per_mwh",         "Battery / storage",       C_STORE),
    ("lcoe_grid_connection_eur_per_mwh", "Grid connection",         C_GRID_CONN),
    ("lcoe_grid_energy_eur_per_mwh",     "Grid energy",             C_GRID_NRG),
    ("lcoe_transmission_eur_per_mwh",    "Transmission (HVDC)",     C_TRANSM),
]
LCOH_PARTS = [
    ("lcoh_electrolyser_eur_per_mwh_lhv", "Electrolyser (capex+opex)", C_CAPEX),
    ("lcoh_h2_storage_eur_per_mwh_lhv",   "H₂ storage (buffer)",       C_H2_STORE),
    ("lcoh_electricity_eur_per_mwh_lhv",  "Electricity (@ LCOE)",      C_ELEC),
]


def _lcoe_part_lines(column: str, row: pd.Series) -> list:
    """Detail under one LCOE band — the bands are the parts, so this goes one level down."""
    if column == "lcoe_renewables_eur_per_mwh":
        own = _value(row, "lcoe_renewables_own_eur_per_mwh")
        lines = [("· own LCOE, over its own generation", _amount(own, "€/MWh", 1))] if own else []
        return lines + _lcoe_tech_lines(row, headline=False)
    # Connection and energy are their own bands here, so each takes only the part of
    # the delivered grid price it is actually paying for.
    if column == "lcoe_grid_energy_eur_per_mwh":
        return [line for line in _grid_lines(row)
                if line[0] in ("· market price", "· volumetric fee")]
    if column == "lcoe_grid_connection_eur_per_mwh":
        return [line for line in _grid_lines(row) if line[0] == "· connection"]
    return []


def _lcoh_part_lines(column: str, row: pd.Series) -> list:
    """Detail under one LCOH band, whose values are €/kg — restate them per MWh LHV."""
    lines = [("· on an LHV basis", _amount(_value(row, column), "€/MWh LHV", 1))]
    if column == "lcoh_electricity_eur_per_mwh_lhv":
        lines.append(("· priced at LCOE", _amount(_value(row, "lcoe_eur_per_mwh"), "€/MWh", 1)))
    return lines


def levelised_parts(parts, lines_for, scale: float = 1.0):
    """Stack one levelised metric across SUPPLY_CASES.

    `parts` is [(report column, legend label, colour)] bottom→top; absent columns
    and all-zero rows drop out, matching how the reports omit parts a scenario
    does not have. `lines_for(column, row)` supplies each band's hover detail.
    """
    frames = {project: _report(project) for project, _, _ in SUPPLY_CASES}
    case_rows = [frames[project].loc[scenario] for project, scenario, _ in SUPPLY_CASES]
    bars = [label for _, _, label in SUPPLY_CASES]
    series = []
    for column, label, colour in parts:
        values = [(_value(row, column) or 0.0) * scale for row in case_rows]
        if all(value <= 0 for value in values):
            continue
        series.append((label, colour, values, [lines_for(column, row) for row in case_rows]))
    return bars, series, case_rows


# ---- Figure assembly -----------------------------------------------------

def stacked_figure(bars, series, *, unit: str, value_fmt: str, total_fmt: str,
                   basis: str, bar_totals, total_label: str, totals_lines) -> go.Figure:
    """A stacked FCA bar chart; `series` is [(label, colour, values, lines)] bottom→top.

    `lines` holds one list of (leading text, formatted value) per bar — the hover
    detail. `bar_totals` are the reported levelised totals, which give each segment
    its "% of BASIS" share. Totals ride above each bar as a text-only trace (the
    house pattern from plot_lcos_bars). The legend is horizontal below the plot so
    the figure keeps its width in a half-page panel.
    """
    fig = go.Figure()
    for label, colour, values, lines in series:
        templates, customdata = [], []
        for index, value in enumerate(values):
            detail = lines[index]
            share = round(value / bar_totals[index] * 100) if bar_totals[index] else 0
            template = f"{label}: %{{y:{value_fmt}}} {unit} (%{{customdata[0]}}% of {basis})"
            for position, (text, _) in enumerate(detail):
                template += f"<br>{text} %{{customdata[{position + 1}]}}"
            templates.append(template + "<extra></extra>")
            customdata.append([share] + [shown for _, shown in detail])
        # Bars differ in how many detail lines they have; pad to a rectangle so
        # customdata stays a plain 2-D array, and each bar's own template only ever
        # reaches the indices it filled.
        width = max(len(row) for row in customdata)
        fig.add_trace(go.Bar(
            x=bars, y=values, name=label, marker_color=colour,
            customdata=[row + [""] * (width - len(row)) for row in customdata],
            hovertemplate=templates,
        ))

    totals = [sum(values[i] for _, _, values, _ in series) for i in range(len(bars))]
    templates, customdata = [], []
    for index in range(len(bars)):
        detail = totals_lines[index]
        template = f"{total_label}: %{{y:{total_fmt}}} {unit}"
        for position, (text, _) in enumerate(detail):
            template += f"<br>{text} %{{customdata[{position}]}}"
        templates.append(template + "<extra></extra>")
        customdata.append([shown for _, shown in detail])

    fig.add_trace(go.Scatter(
        x=bars, y=totals, mode="text",
        text=[format(total, total_fmt) for total in totals],
        textposition="top center",
        textfont=dict(size=13, color=PRINT["ink"], family="Titillium Web"),
        showlegend=False, cliponaxis=False,
        customdata=customdata, hovertemplate=templates,
    ))

    fig.update_layout(
        template=responsive_template, barmode="stack", bargap=0.4,
        autosize=True, height=PLOT_HEIGHT,
        paper_bgcolor=PRINT["paper"], plot_bgcolor=PRINT["paper"],
        font=dict(family="Titillium Web", size=12, color=PRINT["ink"]),
        margin=dict(l=52, r=16, t=28, b=118),
        legend=dict(orientation="h", x=0, y=-0.16, xanchor="left", yanchor="top",
                    traceorder="reversed", bgcolor="rgba(0,0,0,0)",
                    font=dict(family="Titillium Web", size=11.5, color=PRINT["muted"])),
    )
    fig.update_xaxes(type="category", automargin=True, linecolor=PRINT["border"],
                     tickfont=dict(size=12, color=PRINT["muted"]))
    fig.update_yaxes(rangemode="tozero", gridcolor=PRINT["grid"], zeroline=False,
                     tickfont=dict(size=11.5, color=PRINT["muted"]))
    return fig


def plot_div(fig: go.Figure, div_id: str) -> str:
    return fig.to_html(include_plotlyjs=False, full_html=False, config=CONFIG,
                       div_id=div_id, default_width="100%",
                       default_height=f"{PLOT_HEIGHT}px")


# ---- Page copy -----------------------------------------------------------

CARRIER_COLUMNS = [
    dict(
        key="lcoe",
        title="LCOE — levelised cost of electricity",
        unit="€ / MWh delivered",
        src="report columns",
        note="Transmission is absent — no HVDC in the DE cases. p95 is the 95th-percentile "
             "area-weighted capacity-factor cell, so it beats the mean site.",
    ),
    dict(
        key="lcoh",
        title="LCOH — levelised cost of hydrogen",
        unit="€ / kg H₂",
        src="report columns",
        note=f"Converted at {H2_LHV_KWH_PER_KG} kWh/kg LHV. Electricity is ~72% of the total — "
             "the same fact the alternative's hydrogen block folds away above.",
    ),
]

CSS = """
  .cbm-root{
    --paper:#EAEEF1; --panel:#FFFFFF; --panel-2:#F4F6F8;
    --ink:#1B242B; --muted:#525F6A; --faint:#8CA5B7;
    --border:#D6DBDF; --hair:#E4E9ED;
    --accent:#0A5680; --accent-2:#0293D2;
    --shadow:0 1px 2px rgba(27,36,43,.05),0 10px 28px rgba(27,36,43,.07);
    --tw:"Titillium Web",system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
    --mono:"SF Mono",ui-monospace,"Cascadia Mono","Roboto Mono",Menlo,Consolas,monospace;
  }
  @media (prefers-color-scheme:dark){
    .cbm-root{--paper:#10171C; --panel:#1A232A; --panel-2:#212C34;
      --ink:#E9EEF1; --muted:#A7B6C2; --faint:#6C7E8C; --border:#2B3841; --hair:#26323A;
      --accent:#0293D2; --accent-2:#70D2F0;
      --shadow:0 1px 2px rgba(0,0,0,.35),0 12px 34px rgba(0,0,0,.42);}
  }
  :root[data-theme="dark"] .cbm-root{--paper:#10171C; --panel:#1A232A; --panel-2:#212C34;
    --ink:#E9EEF1; --muted:#A7B6C2; --faint:#6C7E8C; --border:#2B3841; --hair:#26323A;
    --accent:#0293D2; --accent-2:#70D2F0;
    --shadow:0 1px 2px rgba(0,0,0,.35),0 12px 34px rgba(0,0,0,.42);}
  :root[data-theme="light"] .cbm-root{--paper:#EAEEF1; --panel:#FFFFFF; --panel-2:#F4F6F8;
    --ink:#1B242B; --muted:#525F6A; --faint:#8CA5B7; --border:#D6DBDF; --hair:#E4E9ED;
    --accent:#0A5680; --accent-2:#0293D2;
    --shadow:0 1px 2px rgba(27,36,43,.05),0 10px 28px rgba(27,36,43,.07);}

  .cbm-root{background:var(--paper);color:var(--ink);font-family:var(--tw);
    line-height:1.5;min-height:100vh;-webkit-font-smoothing:antialiased;}
  .cbm-root *{box-sizing:border-box;}
  .cbm-wrap{max-width:1320px;margin:0 auto;padding:30px 24px 68px;}

  .cbm-head{display:flex;justify-content:space-between;align-items:flex-end;
    gap:24px;flex-wrap:wrap;padding-bottom:4px;}
  .cbm-brand{display:flex;align-items:center;gap:11px;margin-bottom:10px;}
  .cbm-dot{width:12px;height:12px;border-radius:50%;background:var(--accent);flex:none;}
  .cbm-brand span{font-family:var(--mono);font-size:11px;letter-spacing:.2em;
    text-transform:uppercase;color:var(--accent-2);font-weight:600;}
  .cbm-title{font-size:clamp(24px,3.2vw,34px);font-weight:700;letter-spacing:-.02em;
    margin:0;text-wrap:balance;line-height:1.06;}
  .cbm-sub{color:var(--muted);font-size:14px;margin:9px 0 0;max-width:70ch;}
  .cbm-meta{font-family:var(--mono);font-size:11px;color:var(--faint);
    text-align:right;line-height:1.7;white-space:nowrap;}

  .cbm-legend{display:flex;gap:26px;flex-wrap:wrap;margin:22px 0 4px;}
  .cbm-key{display:flex;gap:10px;align-items:flex-start;max-width:46ch;}
  .cbm-key .chip{font-family:var(--mono);font-size:9.5px;font-weight:700;
    letter-spacing:.12em;text-transform:uppercase;border-radius:5px;
    padding:3px 8px;flex:none;line-height:1.6;}
  .cbm-key.alt .chip{background:var(--accent);color:#fff;}
  .cbm-key.dash .chip{background:var(--accent-2);color:#062430;}
  .cbm-key p{margin:0;font-size:12.5px;color:var(--muted);}

  .cbm-panel{background:var(--panel);border:1px solid var(--border);
    border-radius:12px;box-shadow:var(--shadow);margin-top:20px;overflow:hidden;}
  .cbm-panel-h{padding:16px 18px 2px;}
  .cbm-panel-h h2{font-size:15px;font-weight:700;margin:0;letter-spacing:-.01em;}
  .cbm-panel-h .unit{font-size:11.5px;font-weight:400;color:var(--muted);letter-spacing:0;}
  .cbm-panel-h p{font-size:12.5px;color:var(--muted);margin:4px 0 0;max-width:98ch;}

  .cbm-pair{display:grid;grid-template-columns:1fr 1fr;gap:4px 20px;padding:8px 18px 18px;}
  @media (max-width:980px){.cbm-pair{grid-template-columns:1fr;}}
  .cbm-col{min-width:0;}
  .cbm-col-h{display:flex;align-items:baseline;gap:9px;flex-wrap:wrap;
    padding:8px 0 2px;border-bottom:1px solid var(--hair);margin-bottom:6px;}
  .cbm-col-h .chip{font-family:var(--mono);font-size:9.5px;font-weight:700;
    letter-spacing:.12em;text-transform:uppercase;border-radius:5px;
    padding:3px 8px;flex:none;line-height:1.6;}
  .cbm-col.alt .chip{background:var(--accent);color:#fff;}
  .cbm-col.dash .chip{background:var(--accent-2);color:#062430;}
  .cbm-col-h .nm{font-size:13px;font-weight:700;letter-spacing:-.01em;}
  .cbm-col-h .unit{font-size:11.5px;font-weight:400;color:var(--muted);letter-spacing:0;}
  .cbm-col-h .src{font-family:var(--mono);font-size:10px;color:var(--faint);
    margin-left:auto;white-space:nowrap;}
  .cbm-note{font-size:12px;color:var(--muted);margin:8px 2px 0;line-height:1.45;}

  .cbm-root code{font-family:var(--mono);font-size:.88em;background:var(--panel-2);
    border:1px solid var(--hair);border-radius:4px;padding:0 4px;}
  .cbm-foot{margin-top:26px;padding-top:14px;border-top:1px solid var(--border);
    font-size:11.5px;color:var(--faint);line-height:1.6;}
"""


def _column(side: str, chip: str, name: str, src: str, div: str, note: str) -> str:
    """One half of the alternative | dashboard pair."""
    return f"""      <div class="cbm-col {side}">
        <div class="cbm-col-h"><span class="chip">{chip}</span>
          <span class="nm">{name}</span><span class="src">{src}</span></div>
        {div}
        <p class="cbm-note">{note}</p>
      </div>"""


def _metric_column(spec: dict, div: str) -> str:
    """One half of the LCOE | LCOH pair — a metric, not a taxonomy, so no chip."""
    return f"""      <div class="cbm-col">
        <div class="cbm-col-h"><span class="nm">{spec['title']}</span>
          <span class="unit">· {spec['unit']}</span><span class="src">{spec['src']}</span></div>
        {div}
        <p class="cbm-note">{spec['note']}</p>
      </div>"""


def build_core() -> str:
    """The body-only page content, ready to embed as a hub tab."""
    alt_bars, alt_series, alt_lcos, alt_totals = alternative_lcos()
    dash_bars, dash_series, dash_lcos, dash_totals = dashboard_lcos()

    # Both cuts of the same total must land on the reported LCOS, or the page is
    # quietly lying about one of them.
    for name, series, lcos in (("alternative", alt_series, alt_lcos),
                               ("dashboard", dash_series, dash_lcos)):
        totals = [sum(values[i] for _, _, values, _ in series) for i in range(len(lcos))]
        drift = max(abs(total - reported) for total, reported in zip(totals, lcos))
        if drift > 0.05:
            raise ValueError(
                f"the {name} stack drifts {drift:.3f} €/t from the reported LCOS")

    lcoe_bars, lcoe_series, lcoe_rows = levelised_parts(LCOE_PARTS, _lcoe_part_lines)
    lcoh_bars, lcoh_series, lcoh_rows = levelised_parts(
        LCOH_PARTS, _lcoh_part_lines, scale=H2_LHV_KWH_PER_KG / 1000.0)

    lcos_kwargs = dict(unit="€/t", value_fmt=",.0f", total_fmt=",.0f", basis="LCOS",
                       total_label="levelised cost of steel (LCOS)")
    divs = {
        "lcos-alt": plot_div(stacked_figure(alt_bars, alt_series, bar_totals=alt_lcos,
                                           totals_lines=alt_totals, **lcos_kwargs),
                             "cbm-lcos-alt"),
        "lcos-dash": plot_div(stacked_figure(dash_bars, dash_series, bar_totals=dash_lcos,
                                            totals_lines=dash_totals, **lcos_kwargs),
                              "cbm-lcos-dash"),
        "lcoe": plot_div(stacked_figure(
            lcoe_bars, lcoe_series, unit="€/MWh", value_fmt=",.1f", total_fmt=",.1f",
            basis="LCOE", total_label="levelised cost of electricity (LCOE)",
            bar_totals=[_value(row, "lcoe_eur_per_mwh") for row in lcoe_rows],
            totals_lines=[[("steel (LCOS):", _amount(_value(row, "lcos_eur_per_t"), "€/t")),
                           ("hydrogen (LCOH):", _amount(
                               _value(row, "lcoh_eur_per_mwh_lhv"), "€/MWh LHV", 1))]
                          for row in lcoe_rows]), "cbm-lcoe"),
        "lcoh": plot_div(stacked_figure(
            lcoh_bars, lcoh_series, unit="€/kg", value_fmt=",.2f", total_fmt=",.2f",
            basis="LCOH", total_label="levelised cost of hydrogen (LCOH)",
            bar_totals=[(_value(row, "lcoh_eur_per_mwh_lhv") or 0.0)
                        * H2_LHV_KWH_PER_KG / 1000.0 for row in lcoh_rows],
            totals_lines=[[("on an LHV basis:", _amount(
                               _value(row, "lcoh_eur_per_mwh_lhv"), "€/MWh LHV", 1)),
                           ("electricity (LCOE):", _amount(
                               _value(row, "lcoe_eur_per_mwh"), "€/MWh", 1))]
                          for row in lcoh_rows]), "cbm-lcoh"),
    }

    carrier_columns = "\n".join(_metric_column(spec, divs[spec["key"]])
                                for spec in CARRIER_COLUMNS)

    return f"""<style>
{font_css()}
{CSS}</style>

<div class="cbm-root">
  <div class="cbm-wrap">
    <header class="cbm-head">
      <div>
        <div class="cbm-brand"><span class="cbm-dot"></span><span>Cost breakdown · two taxonomies</span></div>
        <h1 class="cbm-title">Two ways to cut the same number</h1>
        <p class="cbm-sub">The circulated chart mockup's category split next to the split the
          pipeline already reports — both drawn from the same DE 2023 results the scenario tab
          reads, and both closing exactly on LCOS. All that differs is where the cuts fall, which
          is what makes the choice between them legible. Hover any segment for the same per-plant
          and per-carrier detail the scenario-comparison tab shows.</p>
      </div>
      <div class="cbm-meta">built {date.today().strftime('%-d %b %Y')}<br>
        DE 2023 · grid-connected<br>mean-CF site</div>
    </header>

    <div class="cbm-legend">
      <div class="cbm-key alt"><span class="chip">Alternative</span>
        <p>The circulated taxonomy: one hydrogen block, process plant split into capex and O&amp;M,
          electricity split into an EAF-melt share and the rest. Rebuilt from the model's cost
          groups — no placeholder numbers left.</p></div>
      <div class="cbm-key dash"><span class="chip">Dashboard</span>
        <p>The pipeline's own cost groups, as the rest of the dashboard already shows them. Same
          bars, same totals, different cuts — renewables, electrolyser and H₂ buffer each stand on
          their own.</p></div>
    </div>

    <section class="cbm-panel">
      <div class="cbm-panel-h">
        <h2>LCOS — levelised cost of steel <span class="unit">· € / t liquid steel</span></h2>
        <p>Where the choice between the two actually changes the picture. Both stacks sum to the
          same LCOS per route; the alternative folds the electrolyser, the H₂ buffer and the
          electricity that made the hydrogen into a single block, while the dashboard keeps them
          apart.</p>
      </div>
      <div class="cbm-pair">
{_column('alt', 'Alternative', 'circulated split', 'reconstructed', divs['lcos-alt'],
         "Hydrogen (all-in) = electrolyser + H₂ buffer + the electricity that made the H₂. "
         "The EAF melt share is the furnace's own MWh/t valued at that scenario's LCOE; "
         "&ldquo;rest of plant&rdquo; is the residual, so the stack closes on LCOS. Capex and "
         "O&amp;M split each plant's annual cost by its config ratio — fixed O&amp;M runs at "
         "≈3&nbsp;% of capex per year, so it lands near a quarter of the band.")}
{_column('dash', 'Dashboard', "model's own split", 'report columns', divs['lcos-dash'],
         "Straight from the <code>cost_*_meur</code> columns, converted to €/t. Its process band "
         "keeps capex and O&amp;M together as the pipeline groups them — the split is in the "
         "hover. The iron stockpile and steel inventory keep their legend entries but draw no "
         "visible band; both sit near 0.05 €/t.")}
      </div>
    </section>

    <section class="cbm-panel">
      <div class="cbm-panel-h">
        <h2>The carriers underneath <span class="unit">· LCOE and LCOH</span></h2>
        <p>Here the two taxonomies very nearly coincide, so there is one chart of each rather than
          a pair. The alternative adds a renewables capex/O&amp;M row to LCOE and an
          electrolyser-O&amp;M row to LCOH; both are derivable the same way the process split is,
          but the reports quote each carrier's parent line and these charts follow the reports.
          Everything else maps one-to-one. Both charts show the same three supply cases, and the
          LCOH on the right is built on the LCOE on the left.</p>
      </div>
      <div class="cbm-pair">
{carrier_columns}
      </div>
    </section>

    <p class="cbm-foot">Cross-chart logic, shared by both: the electricity segment in
      LCOS/LCOH is priced at LCOE; the hydrogen segment in LCOS is LCOH × kg&nbsp;H₂/t; LCOI (not
      shown) is LCOS minus the EAF stage. Two rows of the alternative stay empty because the
      model does not carry them — ore/iron freight and scrap or HBI purchase; freight interacts
      with the ore-grade work, since lower-grade ore moves more tonnes per tonne of iron. Sources:
      <code>results/report_DE-2023-grid.csv</code>,
      <code>results/report_DE-2023-nogrid.csv</code> and
      <code>config/assumptions.yaml</code> with the per-scenario
      <code>config/assumptions_DE-2023-grid_*.yaml</code> overlays. Regenerate with
      <code>python workflow/scripts/viz/build_cost_breakdown_mockups.py</code>.</p>
  </div>
</div>
"""


def main() -> None:
    core = build_core()
    plotly_js = (Path(go.__file__).parents[1] / "package_data" / "plotly.min.js").read_text()
    page = f"<script>{plotly_js}</script>\n{core}"

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(page, encoding="utf-8")
    print(f"wrote {OUT_PATH} ({OUT_PATH.stat().st_size/1e6:.2f} MB, "
          f"{page.count(chr(0))} NULs) — body-only hub fragment")


if __name__ == "__main__":
    main()
