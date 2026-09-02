"""Build a self-contained, interactive FCA dashboard from the pipeline reports.

Rather than embedding static pipeline PNGs, this emits a compact per-scenario
data payload (LCOS cost breakdown in €/t and optimised capacities, computed with
the exact pipeline reshaping in plot_lcos_bars / plot_capacity_bars) and renders
the cost-breakdown and capacity charts client-side. That lets the page toggle:

  * sensitivity — which scenario the run belongs to (a scenario is a name plus
    its optional assumptions overlay), and
  * capacity-factor method — the `variant` each RES tech was solved with.

Both axes are report columns, so the page is assembled by grouping rows rather
than by parsing names.

The controls drive the charts, the LCOS table and the summary cards; the
cross-geography LCOS overview stays on the baseline, cheapest-CF numbers.
Offline, no external hosts — safe to publish as an Artifact.
"""
import base64
import json
import sys
import types
from pathlib import Path

import pandas as pd
import yaml
from plotly.utils import PlotlyJSONEncoder

REPO = Path(__file__).resolve().parents[3]  # workflow/scripts/viz/ -> repo root
RESULTS = REPO / "results"
OUT = RESULTS / "dashboard.html"
TEMPLATE_HTML = Path(__file__).with_name("dashboard_template.html")
CONFIG_DIR = REPO / "config"

# Display names for the areas the registry can name. An area with no entry shows
# its own code, which is what a newly added zone should do until it is named.
GEO_NAMES = {
    "DEU": "Germany", "ESP": "Spain", "FRA": "France",
    "AUS": "Australia", "BRA": "Brazil",
    "NSW1": "New South Wales", "QLD1": "Queensland", "SA1": "South Australia",
    "TAS1": "Tasmania", "VIC1": "Victoria",
}

# Display names for the capacity-factor variants. A variant with no entry shows
# its own id, so a new CF method reaches the page without a code change.
CF_NAMES = {
    "area-average": "Area average",
    "bestsite-p95": "Best-site P95",
}
# Zero-fossil routes only. MIX-DRI-EAF burns natural gas (optimiser picks the H2/gas
# blend), so it is a transitional route, not clean; NG-DRI is the fossil benchmark.
# Route ids as common/_runs.py spells them. h2-only is absent on purpose: it
# makes hydrogen, so it has no cost of steel to compare.
CLEAN_ROUTES = ["h2-dri-eaf", "moe-eaf", "ew-eaf"]
ROUTE_ORDER = ["h2-dri-eaf", "moe-eaf", "ew-eaf", "mix-dri-eaf", "ng-dri-eaf"]
ROUTE_LABEL = {"h2-dri-eaf": "H2-DRI-EAF", "moe-eaf": "MOE", "ew-eaf": "Electrowinning",
               "mix-dri-eaf": "NG-H2-DRI-EAF", "ng-dri-eaf": "NG-DRI-EAF"}
# Blues = clean; sand = transitional (partial gas); red = fossil.
ROUTE_COLOR = {"h2-dri-eaf": "#0A5680", "moe-eaf": "#0293D2", "ew-eaf": "#83D1DD",
               "mix-dri-eaf": "#E2B681", "ng-dri-eaf": "#D75674"}
CO2_T_PER_MWH = 0.20
# Below this iron-from-H2 share the "MIX" route has effectively rejected hydrogen
# and is economically the fossil NG-DRI route — flagged as such in the UI.
H2_MIN = 0.02

# Cost-group stack order (bottom → top), label and colour — mirrors COST_GROUPS in
# plot_lcos_bars.py. Keys match the cost_{group}_meur columns in the report.
COST_GROUPS = [
    ["process",         "Process plant (capex+opex)", "#33434D"],
    ["ore_consumables", "Ore & consumables",          "#E2B681"],
    ["gas",             "Natural gas (fuel)",         "#525F6A"],
    ["iron_store",      "Iron stockpile",             "#B7C1C8"],
    ["steel_store",     "Steel inventory",            "#8CA5B7"],
    ["electrolyser",    "Electrolyser",               "#91C096"],
    ["h2_buffer",       "H₂ buffer",                  "#70D2F0"],
    ["res",             "Renewables (capex+opex)",    "#0A5680"],
    ["battery",         "Battery",                    "#D75674"],
    ["transmission",    "Transmission (HVDC)",        "#83D1DD"],
    ["grid",            "Grid (connection + energy)", "#71828F"],
]

# Capacity panels: (title, unit, [(key, label, colour), ...]). Keys are columns
# emitted by plot_capacity_bars.build_plot_data (solar/wind summed across any
# orientation/tech sub-columns during extraction).
CAP_PANELS = [
    ["Power capacity (MW)", "MW", [
        ["solar_mw",         "Solar",                "#E2B681"],
        ["wind-onshore_mw",  "Wind onshore",         "#0A5680"],
        ["wind-offshore_mw", "Wind offshore",        "#3E7CB1"],
        ["battery_mw",       "Battery",              "#D75674"],
        ["electrolyser_mw",  "Electrolyser (input)", "#91C096"],
        ["grid_mw",          "Grid connection",      "#71828F"],
    ]],
    ["Process capacity (t/h output)", "t/h", [
        ["dri-h2_t_per_h", "H2-DRI shaft",   "#8CA5B7"],
        ["dri-ng_t_per_h", "NG-DRI shaft",   "#BDCCD9"],
        ["eaf_t_per_h",    "EAF",            "#525F6A"],
        ["moe_t_per_h",    "MOE",            "#0293D2"],
        ["ew_t_per_h",     "Electrowinning", "#83D1DD"],
    ]],
    ["Storage (hours of demand)", "h", [
        ["h2_buffer_hours_dri",     "H₂ buffer",       "#70D2F0"],
        ["iron_store_hours_steel",  "Iron store",      "#B7C1C8"],
        ["steel_store_hours_steel", "Steel inventory", "#8CA5B7"],
    ]],
]
_CAP_KEYS = [b[0] for _, _, bars in CAP_PANELS for b in bars]


def _seed_stub():
    """Provide a working `snakemake` so the plot modules import cleanly standalone."""
    sm = types.SimpleNamespace(
        input=types.SimpleNamespace(report="_"), output=["_"],
        log=["/dev/null"], wildcards=types.SimpleNamespace(scenario="_"),
        config={}, params=types.SimpleNamespace(), threads=1,
    )
    mod = types.ModuleType("common._stubs")
    mod.snakemake = sm
    sys.modules["common._stubs"] = mod


def _num(v, default=0.0):
    try:
        f = float(v)
        return f if pd.notna(f) else default
    except (TypeError, ValueError):
        return default


def _opt(v):
    """Round to 1 dp, or None if missing/NaN (so absent LCOH stays absent, not 0)."""
    try:
        f = float(v)
        return round(f, 1) if pd.notna(f) else None
    except (TypeError, ValueError):
        return None


def _axes(row):
    """The six axes the page browses, read off one report row.

    Each is a column now: the area, the weather year, whether a grid price was
    among the run's inputs, the route, the scenario it belongs to, and the CF
    method its renewables were solved with. The old build parsed all but the
    first two out of a project and scenario name.
    """
    variants = {
        col[: -len("_variant")]: row[col]
        for col in row.index
        if col.endswith("_variant") and pd.notna(row[col])
    }
    cf_methods = sorted({v for tech, v in variants.items() if tech != "grid"})
    return {
        "geo": row["area"],
        "year": str(row["start_date"])[:4],
        "grid": "grid" if "grid" in variants else "nogrid",
        "route": row["route"],
        "variant": row["scenario"],
        # A run whose techs were solved with different methods is named by all of
        # them, rather than being filed under whichever came first.
        "cf": "+".join(cf_methods) or "na",
    }


def _baseline_scenario(scenarios):
    """The scenario the overview rests on and the others fall back to.

    A scenario is a name plus an optional assumptions overlay, so the one with no
    overlay is the base case. If every scenario overrides something the first by
    name stands in — the comparison has to be against something.
    """
    plain = [s for s in scenarios if not (CONFIG_DIR / f"assumptions_{s}.yaml").exists()]
    baseline = sorted(plain or scenarios)[0]
    return baseline


def _record(row, lcos_row, cap_row):
    """Assemble one scenario's payload record (scalars + cost breakdown + caps)."""
    costs = {}
    for group, _, _ in COST_GROUPS:
        if group in lcos_row:
            v = _num(lcos_row[group])
            if v > 0.05:
                costs[group] = round(v, 1)
    caps = {}
    solar = sum(_num(cap_row[c]) for c in cap_row.index if str(c).startswith("solar"))
    wind = sum(_num(cap_row[c]) for c in cap_row.index if str(c).startswith("wind-onshore"))
    wind_off = sum(_num(cap_row[c]) for c in cap_row.index if str(c).startswith("wind-offshore"))
    for key in _CAP_KEYS:
        if key == "solar_mw":
            v = solar
        elif key == "wind-onshore_mw":
            v = wind
        elif key == "wind-offshore_mw":
            v = wind_off
        else:
            v = _num(cap_row[key]) if key in cap_row.index else 0.0
        if v > 0.001:
            caps[key] = round(v, 1)
    # LCOE breakdown in €/MWh (parts sum to LCOE); absent parts stay absent.
    lcoe_parts = {}
    for col, key in (
        ("lcoe_renewables_eur_per_mwh", "renewables"),
        ("lcoe_solar_eur_per_mwh", "solar"),
        ("lcoe_wind_onshore_eur_per_mwh", "wind_onshore"),
        ("lcoe_wind_offshore_eur_per_mwh", "wind_offshore"),
        ("lcoe_storage_eur_per_mwh", "storage"),
        ("lcoe_grid_connection_eur_per_mwh", "grid_connection"),
        ("lcoe_grid_energy_eur_per_mwh", "grid_energy"),
        ("lcoe_transmission_eur_per_mwh", "transmission"),
    ):
        v = _opt(row.get(col))
        if v is not None:
            lcoe_parts[key] = v
    # Per-tech LCOE over each tech's own generation (the unit cost of the
    # renewable itself), distinct from the lcoe_parts contribution above.
    lcoe_own = {}
    for col, key in (
        ("lcoe_renewables_own_eur_per_mwh", "renewables"),
        ("lcoe_solar_own_eur_per_mwh", "solar"),
        ("lcoe_wind_onshore_own_eur_per_mwh", "wind_onshore"),
        ("lcoe_wind_offshore_own_eur_per_mwh", "wind_offshore"),
    ):
        v = _opt(row.get(col))
        if v is not None:
            lcoe_own[key] = v
    # LCOH breakdown in €/MWh LHV (parts sum to LCOH); absent parts stay absent.
    lcoh_parts = {}
    for col, key in (
        ("lcoh_electrolyser_eur_per_mwh_lhv", "electrolyser"),
        ("lcoh_electricity_eur_per_mwh_lhv", "electricity"),
        ("lcoh_h2_storage_eur_per_mwh_lhv", "storage"),
    ):
        v = _opt(row.get(col))
        if v is not None:
            lcoh_parts[key] = v
    # €/t-steel cost detail: per-plant capex (sums to the "process" group) plus the
    # ore/consumables split (sums to "ore_consumables"). For the segment hovers.
    cost_detail = {}
    for col, key in (
        ("plant_dri-h2_eur_per_t", "dri"),
        ("plant_dri-ng_eur_per_t", "dri_ng"),
        ("plant_eaf_eur_per_t", "eaf"),
        ("plant_moe_eur_per_t", "moe"),
        ("plant_ew_eur_per_t", "electrowinning"),
        ("ore_eur_per_t_steel", "ore"),
        ("consumables_eur_per_t_steel", "consumables"),
    ):
        v = _opt(row.get(col))
        if v is not None:
            cost_detail[key] = v
    # Capacity factors as integer percent (0–1 fraction × 100 — 1 dp would collapse
    # a 13% CF to 10%), keyed by the capacity-panel bar key they annotate.
    cf = {}
    for capkey, col in (
        ("solar_mw", "cf_solar"),
        ("wind-onshore_mw", "cf_wind_onshore"),
        ("wind-offshore_mw", "cf_wind_offshore"),
        ("grid_mw", "cf_grid_connection"),
        ("electrolyser_mw", "electrolyser_utilization"),
        ("dri-h2_t_per_h", "dri-h2_utilization"),
        ("dri-ng_t_per_h", "dri-ng_utilization"),
        ("eaf_t_per_h", "eaf_utilization"),
        ("moe_t_per_h", "moe_utilization"),
        ("ew_t_per_h", "ew_utilization"),
    ):
        raw = row.get(col)
        if raw is not None and pd.notna(raw):
            cf[capkey] = round(float(raw) * 100)
    return {
        "lcos": round(_num(row.get("lcos_eur_per_t")), 0),
        "lcoe": _opt(row.get("lcoe_eur_per_mwh")),
        "lcoh": _opt(row.get("lcoh_eur_per_mwh_lhv")),
        "steel_mt": round(_num(row.get("steel_produced_mt")), 4),
        "ng_gwh": round(_num(row.get("ng_gwh_lhv")), 1),
        "h2_share": round(_num(row.get("iron_from_h2_share")), 3),
        "costs": costs,
        "caps": caps,
        "lcoe_parts": lcoe_parts,
        "lcoe_own": lcoe_own,
        "lcoh_parts": lcoh_parts,
        "grid_price": _opt(row.get("grid_price_eur_per_mwh")),
        "grid_fee": _opt(row.get("grid_fee_eur_per_mwh")),
        "grid_conn_imported": _opt(row.get("grid_connection_eur_per_mwh_imported")),
        "cost_detail": cost_detail,
        "cf": cf,
    }


def _gas_price(scenario):
    """The €/MWh natural-gas price a scenario assumed (its overlay, else the base).

    One overlay covers a whole scenario, so unlike the old per-project lookup
    there is only one file that can carry the price.
    """
    base = yaml.safe_load((CONFIG_DIR / "assumptions.yaml").read_text()) or {}
    price = base.get("natural_gas", {}).get("price_eur_per_mwh")
    overlay_path = CONFIG_DIR / f"assumptions_{scenario}.yaml"
    if overlay_path.exists():
        overlay = yaml.safe_load(overlay_path.read_text()) or {}
        price = overlay.get("natural_gas", {}).get("price_eur_per_mwh", price)
    return price


def _co2_t_per_mwh():
    """Gas combustion CO2 intensity (t/MWh LHV) from assumptions; constant as fallback."""
    base = yaml.safe_load((REPO / "config" / "assumptions.yaml").read_text()) or {}
    return base.get("natural_gas", {}).get("co2_t_per_mwh", CO2_T_PER_MWH)


def _default_view(cases, baseline, cf_options):
    """The pair of scenarios the page opens on, chosen from what was solved.

    The old build named a geography, year and sensitivity outright, so the page
    opened on an "unavailable" notice whenever that one run had not been solved.
    B differs from A by route alone; a null axis means B tracks A.
    """
    project = sorted(cases)[0]
    geo, year, grid = project.rsplit("-", 2)
    routes = [route for route in ROUTE_ORDER if route in cases[project]]
    primary = "h2-dri-eaf" if "h2-dri-eaf" in routes else routes[0]
    secondary = next((route for route in routes if route != primary), primary)
    scenarios = cases[project][primary]
    view = {
        "a": {
            "geo": geo, "year": year, "grid": grid, "route": primary,
            "cf": cf_options[0][0],
            "variant": baseline if baseline in scenarios else sorted(scenarios)[0],
        },
        "b": {"geo": None, "year": None, "grid": None, "route": secondary,
              "cf": None, "variant": None},
    }
    return view


def build_payload(report_paths):
    """Return (cases, synth, gas, axis_options) from every scenario report.

    cases[project][route][scenario][cf] = record — the full interactive data,
    where project is the `{geo}-{year}-{grid}` key the page browses by.
    synth[geo][year][grid][route] = baseline cheapest-CF record — for the overview.

    One report holds one scenario's runs, spanning areas, years and routes, so a
    single file can contribute to several projects. That is the inversion from
    the old build, where a file *was* a project.
    """
    sys.path.insert(0, str(REPO / "workflow"))
    _seed_stub()
    import scripts.viz.plot_lcos_bars as L
    import scripts.viz.plot_capacity_bars as C
    from _run_display import run_label

    scenarios = sorted(p.stem[len("report_"):] for p in report_paths)
    baseline = _baseline_scenario(scenarios)

    cases, synth, gas = {}, {}, {}
    geos, years, cf_methods = set(), set(), set()
    for report_path in sorted(report_paths):
        df = pd.read_csv(report_path)
        lcos_df = L.build_plot_data(df)   # €/t cost groups, indexed by run label
        cap_df = C.build_plot_data(df)    # capacities, indexed by run label

        for _, row in df.iterrows():
            label = run_label(row)
            if label not in lcos_df.index:   # no LCOS (h2-only) — nothing to compare
                continue
            axes = _axes(row)
            geo, year, grid = axes["geo"], axes["year"], axes["grid"]
            geos.add(geo)
            years.add(year)
            cf_methods.add(axes["cf"])
            gas.setdefault(geo, {}).setdefault(year, {})[grid] = _gas_price(axes["variant"])

            rec = _record(row, lcos_df.loc[label], cap_df.loc[label])
            project = f"{geo}-{year}-{grid}"
            (cases.setdefault(project, {})
                  .setdefault(axes["route"], {})
                  .setdefault(axes["variant"], {}))[axes["cf"]] = rec

            if axes["variant"] == baseline and rec["lcos"] > 0:
                bucket = synth.setdefault(geo, {}).setdefault(year, {}).setdefault(grid, {})
                previous = bucket.get(axes["route"])
                if previous is None or rec["lcos"] < previous["lcos"]:
                    bucket[axes["route"]] = {
                        "lcos": rec["lcos"], "steel_mt": rec["steel_mt"],
                        "ng_gwh": rec["ng_gwh"], "h2_share": rec["h2_share"],
                        "lcoe": rec["lcoe"], "lcoh": rec["lcoh"],
                    }

    # Which routes a scenario actually solved — the page greys out the rest
    # rather than offering a combination that has no record behind it.
    variant_routes = {}
    for case in cases.values():
        for route, by_scenario in case.items():
            for scenario in by_scenario:
                variant_routes.setdefault(scenario, set()).add(route)

    axis_options = {
        "geos": sorted(geos),
        "years": sorted(years),
        "base_variant": baseline,
        "variant_label": {s: s for s in scenarios},
        "variant_routes": {s: sorted(r) for s, r in variant_routes.items()},
        "cf_options": [[c, CF_NAMES.get(c, c)] for c in sorted(cf_methods)],
    }
    if not cases:
        raise ValueError(
            "no run with a cost of steel in any report — the page compares steel "
            "routes, and h2-only produces hydrogen"
        )
    axis_options["default_view"] = _default_view(
        cases, baseline, axis_options["cf_options"]
    )
    return cases, synth, gas, axis_options


def font_css():
    """Inline the bundled Titillium Web woff2 subsets as @font-face data-URIs."""
    assets = REPO / "workflow" / "scripts" / "viz" / "assets"
    faces = []
    for weight, name in ((400, "TitilliumWeb-Regular.woff2"), (600, "TitilliumWeb-SemiBold.woff2")):
        p = assets / name
        if p.exists():
            b64 = base64.b64encode(p.read_bytes()).decode("ascii")
            faces.append(f"@font-face{{font-family:'Titillium Web';font-style:normal;font-weight:{weight};"
                         f"src:url('data:font/woff2;base64,{b64}') format('woff2');}}")
    return "".join(faces)


def build_html(template_path: Path, augment=None):
    """Assemble a self-contained dashboard HTML from a template, returning (html, cases, geos).

    The payload (cases/synth/gas + FCA template + label/colour maps) is identical
    for every template — v1 and the scenario-comparison v2 differ only in markup
    and client JS. Templates carry the /*FONT_CSS*/, /*PLOTLY_JS*/ and
    /*PAYLOAD_JSON*/ placeholders.

    `augment(payload, cases)` lets a template that needs more than the shared
    payload add to it before serialisation — the cost-taxonomy page attaches its
    leaf-level split that way, rather than growing every page's payload.
    """
    # Whatever has been compiled: one report per scenario, discovered on disk
    # rather than listed anywhere. The hidden `_diag.csv` siblings are not
    # reports — they hold the zones the report already decided against.
    report_paths = sorted(RESULTS.glob("report_*.csv"))
    if not report_paths:
        raise FileNotFoundError(
            f"no scenario reports in {RESULTS} — run `snakemake` to compile some"
        )
    cases, synth, gas, axis_options = build_payload(report_paths)

    # FCA plotly template (shared by the overview and the client-rendered charts).
    sys.path.insert(0, str(REPO / "workflow"))
    _seed_stub()
    import scripts.viz.plot_lcos_bars as L
    tpl = L.fca_template
    template = tpl.to_plotly_json() if hasattr(tpl, "to_plotly_json") else tpl

    plotly_js = (Path(__import__("plotly").__file__).parent / "package_data" / "plotly.min.js").read_text()

    payload = {
        "cases": cases, "synth": synth, "gas": gas, "template": template,
        "geo_names": GEO_NAMES,
        "clean_routes": CLEAN_ROUTES, "route_order": ROUTE_ORDER,
        "route_label": ROUTE_LABEL, "route_color": ROUTE_COLOR,
        "co2_t_per_mwh": _co2_t_per_mwh(), "h2_min": H2_MIN,
        "cost_groups": COST_GROUPS, "cap_panels": CAP_PANELS,
        # geos, years, base_variant, variant_label, variant_routes, cf_options —
        # every one of them read off the reports rather than declared here, so a
        # new area, year, scenario or CF method reaches the page on its own.
        **axis_options,
    }
    if augment is not None:
        augment(payload, cases)

    html = (template_path.read_text()
            .replace("/*FONT_CSS*/", font_css())
            .replace("/*PLOTLY_JS*/", plotly_js)
            .replace("/*PAYLOAD_JSON*/", json.dumps(payload, cls=PlotlyJSONEncoder)))
    return html, cases, axis_options["geos"]


def main():
    html, cases, geos = build_html(TEMPLATE_HTML)
    OUT.write_text(html)
    print(f"wrote {OUT} ({OUT.stat().st_size/1e6:.2f} MB) — {len(cases)} projects, {len(geos)} geos")


if __name__ == "__main__":
    main()
