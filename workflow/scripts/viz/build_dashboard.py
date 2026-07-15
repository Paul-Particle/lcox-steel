"""Build a self-contained, interactive FCA dashboard from the pipeline reports.

Rather than embedding static pipeline PNGs, this emits a compact per-scenario
data payload (LCOS cost breakdown in €/t and optimised capacities, computed with
the exact pipeline reshaping in plot_lcos_bars / plot_capacity_bars) and renders
the cost-breakdown and capacity charts client-side. That lets the page toggle:

  * sensitivity variant — Baseline / MOE turndown 70% / H₂ salt cavern
    (swap-in-place: only the affected route's bars change), and
  * capacity-factor mode — average CF / best-site P95 / compare both.

The controls drive the charts, the LCOS table and the summary cards; the
cross-geography LCOS overview stays on the baseline, cheapest-CF numbers.
Offline, no external hosts — safe to publish as an Artifact.
"""
import base64
import json
import re
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
PROJ_RE = re.compile(r"^(DE|ES|FR|AUS|BRA)-(2023|2024|2025)-(grid|nogrid)$")

GEO_NAMES = {"DE": "Germany", "ES": "Spain", "FR": "France", "AUS": "Australia", "BRA": "Brazil"}
# Zero-fossil routes only. MIX-DRI-EAF burns natural gas (optimiser picks the H2/gas
# blend), so it is a transitional route, not clean; NG-DRI is the fossil benchmark.
CLEAN_ROUTES = ["h2-dri-eaf", "moe", "ew"]
ROUTE_ORDER = ["h2-dri-eaf", "moe", "ew", "mix-dri-eaf", "ng-dri-eaf"]
ROUTE_LABEL = {"h2-dri-eaf": "H2-DRI-EAF", "moe": "MOE", "ew": "Electrowinning",
               "mix-dri-eaf": "MIX-DRI-EAF", "ng-dri-eaf": "NG-DRI-EAF"}
# Blues = clean; sand = transitional (partial gas); red = fossil.
ROUTE_COLOR = {"h2-dri-eaf": "#0A5680", "moe": "#0293D2", "ew": "#83D1DD",
               "mix-dri-eaf": "#E2B681", "ng-dri-eaf": "#D75674"}
CO2_T_PER_MWH = 0.20
# Below this iron-from-H2 share the "MIX" route has effectively rejected hydrogen
# and is economically the fossil NG-DRI route — flagged as such in the UI.
H2_MIN = 0.02

# Which routes each sensitivity variant swaps out. "base" leaves every route alone.
VARIANT_ROUTES = {"base": [], "moe70": ["moe"], "salt": ["h2-dri-eaf", "mix-dri-eaf"]}
VARIANT_LABEL = {"base": "Baseline", "moe70": "MOE turndown 70%", "salt": "H₂ salt cavern"}

# Cost-group stack order (bottom → top), label and colour — mirrors COST_GROUPS in
# plot_lcos_bars.py. Keys match the cost_{group}_meur columns in the report.
COST_GROUPS = [
    ["process",         "Process plant (capex+opex)", "#33434D"],
    ["ore_consumables", "Ore & consumables",          "#E2B681"],
    ["gas",             "Natural gas (fuel + CO₂)",    "#525F6A"],
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
        ["battery_mw",       "Battery",              "#D75674"],
        ["electrolyser_mw",  "Electrolyser (input)", "#91C096"],
        ["grid_mw",          "Grid connection",      "#71828F"],
    ]],
    ["Process capacity (t/h output)", "t/h", [
        ["dri_t_per_h",            "H2-DRI shaft",   "#8CA5B7"],
        ["dri_ng_t_per_h",         "NG-DRI shaft",   "#BDCCD9"],
        ["eaf_t_per_h",            "EAF",            "#525F6A"],
        ["moe_t_per_h",            "MOE",            "#0293D2"],
        ["electrowinning_t_per_h", "Electrowinning", "#83D1DD"],
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
        log=["/dev/null"], wildcards=types.SimpleNamespace(project="_"),
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


def _parse_scenario(scen):
    """(route, variant, cf) for a report scenario name.

    ng-dri-eaf has no CF variants (no RES) → cf 'na'. A trailing -avg/-p95 is the
    CF; the remaining stem carries the sensitivity variant: '…-t70' → moe70 (on the
    moe route), '…-salt' → salt, otherwise baseline.
    """
    if scen == "ng-dri-eaf":
        return "ng-dri-eaf", "base", "na"
    cf, stem = "na", scen
    for suf in ("-avg", "-p95"):
        if scen.endswith(suf):
            cf, stem = suf[1:], scen[: -len(suf)]
            break
    if stem.endswith("-t70"):
        return stem[: -len("-t70")], "moe70", cf
    if "-salt" in stem:
        return stem.replace("-salt", ""), "salt", cf
    return stem, "base", cf


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
    for key in _CAP_KEYS:
        if key == "solar_mw":
            v = solar
        elif key == "wind-onshore_mw":
            v = wind
        else:
            v = _num(cap_row[key]) if key in cap_row.index else 0.0
        if v > 0.001:
            caps[key] = round(v, 1)
    return {
        "lcos": round(_num(row.get("lcos_eur_per_t")), 0),
        "lcoe": _opt(row.get("lcoe_eur_per_mwh")),
        "lcoh": _opt(row.get("lcoh_eur_per_mwh_lhv")),
        "steel_mt": round(_num(row.get("steel_produced_mt")), 4),
        "ng_gwh": round(_num(row.get("ng_gwh_lhv")), 1),
        "h2_share": round(_num(row.get("iron_from_h2_share")), 3),
        "costs": costs,
        "caps": caps,
    }


def _gas_price(project):
    """The €/MWh natural-gas price assumed for a project (per-geo overlay, else base)."""
    base = yaml.safe_load((REPO / "config" / "assumptions.yaml").read_text()) or {}
    fallback = base.get("natural_gas", {}).get("price_eur_per_mwh")
    for scen in ("ng-dri-eaf", "mix-dri-eaf-avg", "mix-dri-eaf-p95"):
        p = REPO / "config" / f"assumptions_{project}_{scen}.yaml"
        if p.exists():
            ov = yaml.safe_load(p.read_text()) or {}
            price = ov.get("natural_gas", {}).get("price_eur_per_mwh")
            if price is not None:
                return price
    return fallback


def build_payload(projects):
    """Return (cases, synth, gas, geos, years).

    cases[project][route][variant][cf] = record — the full interactive data.
    synth[geo][year][grid][route] = baseline cheapest-CF record — for the overview.
    """
    sys.path.insert(0, str(REPO / "workflow"))
    _seed_stub()
    import scripts.viz.plot_lcos_bars as L
    import scripts.viz.plot_capacity_bars as C

    cases, synth, gas = {}, {}, {}
    geos, years = set(), set()
    for project in projects:
        rp = RESULTS / f"report_{project}.csv"
        m = PROJ_RE.match(project)
        if not m or not rp.exists():
            continue
        geo, year, grid = m.groups()
        geos.add(geo)
        years.add(year)
        df = pd.read_csv(rp)
        lcos_df = L.build_plot_data(df)     # €/t cost groups, indexed by scenario
        cap_df = C.build_plot_data(df)      # capacities, indexed by scenario
        gas.setdefault(geo, {}).setdefault(year, {})[grid] = _gas_price(project)

        case = cases.setdefault(project, {})
        ov_bucket = synth.setdefault(geo, {}).setdefault(year, {}).setdefault(grid, {})
        for _, row in df.iterrows():
            scen = str(row["scenario"])
            if scen not in lcos_df.index:      # no LCOS (e.g. h2_only) — skip
                continue
            route, variant, cf = _parse_scenario(scen)
            rec = _record(row, lcos_df.loc[scen], cap_df.loc[scen])
            case.setdefault(route, {}).setdefault(variant, {})[cf] = rec

            if variant == "base" and rec["lcos"] > 0:   # overview: baseline cheapest CF
                prev = ov_bucket.get(route)
                if prev is None or rec["lcos"] < prev["lcos"]:
                    ov_bucket[route] = {
                        "lcos": rec["lcos"], "steel_mt": rec["steel_mt"],
                        "ng_gwh": rec["ng_gwh"], "h2_share": rec["h2_share"],
                        "lcoe": rec["lcoe"], "lcoh": rec["lcoh"],
                    }
    return cases, synth, gas, sorted(geos), sorted(years)


def main():
    projects = sorted(
        p for p in {r["project"] for r in _read_projects()}
        if PROJ_RE.match(p) and (RESULTS / f"report_{p}.csv").exists()
    )
    cases, synth, gas, geos, years = build_payload(projects)

    # FCA plotly template (shared by the overview and the client-rendered charts).
    sys.path.insert(0, str(REPO / "workflow"))
    _seed_stub()
    import scripts.viz.plot_lcos_bars as L
    tpl = L.fca_template
    template = tpl.to_plotly_json() if hasattr(tpl, "to_plotly_json") else tpl

    faces = []
    assets = REPO / "workflow" / "scripts" / "viz" / "assets"
    b64 = {}
    for name in ("TitilliumWeb-Regular.woff2", "TitilliumWeb-SemiBold.woff2"):
        p = assets / name
        if p.exists():
            b64[name] = base64.b64encode(p.read_bytes()).decode("ascii")

    def face(family, weight, name):
        return (f"@font-face{{font-family:'{family}';font-style:normal;font-weight:{weight};"
                f"src:url('data:font/woff2;base64,{b64[name]}') format('woff2');}}")

    if "TitilliumWeb-Regular.woff2" in b64:
        faces.append(face("Titillium Web", 400, "TitilliumWeb-Regular.woff2"))
    if "TitilliumWeb-SemiBold.woff2" in b64:
        faces.append(face("Titillium Web", 600, "TitilliumWeb-SemiBold.woff2"))
    font_css = "".join(faces)

    plotly_js = (Path(__import__("plotly").__file__).parent / "package_data" / "plotly.min.js").read_text()

    payload = {
        "cases": cases, "synth": synth, "gas": gas, "template": template,
        "geos": geos, "years": years, "geo_names": GEO_NAMES,
        "clean_routes": CLEAN_ROUTES, "route_order": ROUTE_ORDER,
        "route_label": ROUTE_LABEL, "route_color": ROUTE_COLOR,
        "co2_t_per_mwh": CO2_T_PER_MWH, "h2_min": H2_MIN,
        "cost_groups": COST_GROUPS, "cap_panels": CAP_PANELS,
        "variant_routes": VARIANT_ROUTES, "variant_label": VARIANT_LABEL,
    }

    template_html = TEMPLATE_HTML.read_text()
    html = (template_html
            .replace("/*FONT_CSS*/", font_css)
            .replace("/*PLOTLY_JS*/", plotly_js)
            .replace("/*PAYLOAD_JSON*/", json.dumps(payload, cls=PlotlyJSONEncoder)))
    OUT.write_text(html)
    print(f"wrote {OUT} ({OUT.stat().st_size/1e6:.2f} MB) — {len(cases)} projects, {len(geos)} geos")


def _read_projects():
    import csv
    with open(REPO / "config" / "projects.csv") as f:
        return list(csv.DictReader(f))


if __name__ == "__main__":
    main()
