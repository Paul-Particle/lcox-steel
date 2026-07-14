"""Build a self-contained FCA dashboard wrapping the EXACT viz-pipeline figures.

Captures each project's real plot_lcos_bars / plot_capacity_bars figure in-process
(monkeypatching save_figure so nothing is written and no kaleido is needed), then
emits one HTML that shares a single inlined plotly.js and the FCA template, adds a
cross-geography synthesis (summary cards + LCOS overview) on top, and drives the
per-project pipeline figures from geography/year/grid selectors. Offline, no
external hosts — safe to publish as an Artifact.
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
ROUTE_LABEL = {"h2-dri-eaf": "H2-DRI-EAF", "moe": "MOE", "ew": "Electrowinning",
               "mix-dri-eaf": "MIX-DRI-EAF", "ng-dri-eaf": "NG-DRI-EAF"}
# Blues = clean; sand = transitional (partial gas); red = fossil.
ROUTE_COLOR = {"h2-dri-eaf": "#0A5680", "moe": "#0293D2", "ew": "#83D1DD",
               "mix-dri-eaf": "#E2B681", "ng-dri-eaf": "#D75674"}
CO2_T_PER_MWH = 0.20
# Below this iron-from-H2 share the "MIX" route has effectively rejected hydrogen
# and is economically the fossil NG-DRI route — flagged as such in the UI.
H2_MIN = 0.02


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


def _project_label(df):
    return ", ".join(dict.fromkeys(df["project"].astype(str)))


def capture_figures(projects):
    """Return {project: {'lcos': figdict, 'capacity': figdict}} and the FCA template dict."""
    sys.path.insert(0, str(REPO / "workflow"))
    _seed_stub()
    import scripts.viz.plot_lcos_bars as L
    import scripts.viz.plot_capacity_bars as C

    grabbed = {}
    def _grab(fig, out_dir, stem, **kw):
        grabbed["fig"] = fig
        return []
    L.save_figure = _grab
    C.save_figure = _grab

    tpl = L.fca_template
    template = tpl.to_plotly_json() if hasattr(tpl, "to_plotly_json") else tpl

    figs = {}
    for project in projects:
        rp = RESULTS / f"report_{project}.csv"
        if not rp.exists():
            continue
        df = pd.read_csv(rp)
        label = _project_label(df)
        entry = {}
        for key, mod in (("lcos", L), ("capacity", C)):
            grabbed.clear()
            try:
                mod.plot(mod.build_plot_data(df), Path(project), label)
            except Exception as e:  # a project may lack a plottable scenario set
                print(f"  ! {project} {key}: {type(e).__name__}: {e}")
                continue
            fig = grabbed.get("fig")
            if fig is None:
                continue
            d = fig.to_dict()
            d.get("layout", {}).pop("template", None)  # template shared globally, not per-fig
            # Keep the figure's designed pixel width/height. Rendering it at native
            # size (and letting the panel scroll horizontally) preserves the exact
            # pipeline geometry — the header dot/logo, margins and 45° labels are all
            # pixel-placed for that width, so shrinking to the panel collides them.
            _tune_layout(d, key)
            entry[key] = d
        if entry:
            figs[project] = entry
    return figs, template


def _tune_layout(d, key):
    """Dashboard-display tweaks on the captured pipeline figure (no pipeline change)."""
    lay = d["layout"]
    if key == "lcos":
        # 45° scenario labels get clipped for long names — let plotly grow the margin.
        lay.setdefault("xaxis", {})["automargin"] = True
    elif key == "capacity":
        # Move the shared 45° scenario labels from the bottom to the top of panel 1,
        # so they read right where the panels start. Hide them on the bottom axis.
        xkeys = [k for k in lay if re.match(r"xaxis\d*$", k)]
        xnum = lambda k: int(k[5:]) if k[5:] else 1
        bottom = max(xkeys, key=xnum)
        lay.setdefault("xaxis", {}).update(
            side="top", showticklabels=True, tickangle=-45, automargin=True)
        if bottom != "xaxis":
            lay.setdefault(bottom, {})["showticklabels"] = False
        m = lay.setdefault("margin", {})
        m["t"] = max(m.get("t", 0), 200)


def _num(v):
    try:
        f = float(v)
        return f if pd.notna(f) else 0.0
    except (TypeError, ValueError):
        return 0.0


def _opt(v):
    """Round to 1 dp, or None if missing/NaN (so absent LCOH stays absent, not 0)."""
    try:
        f = float(v)
        return round(f, 1) if pd.notna(f) else None
    except (TypeError, ValueError):
        return None


def _gas_price(project):
    """The €/MWh natural-gas price assumed for a project (per-geo overlay, else base).

    The fossil routes (ng-dri-eaf on grid, mix on islanded) carry the gas-price
    override; read it from whichever overlay exists, falling back to the base
    assumptions value.
    """
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


def load_synthesis(projects):
    """Per-(geo,year,grid,route) min-LCOS + fossil data for the cards & overview."""
    data, gas, geos, years = {}, {}, set(), set()
    for project in projects:
        m = PROJ_RE.match(project)
        rp = RESULTS / f"report_{project}.csv"
        if not m or not rp.exists():
            continue
        geo, year, grid = m.groups()
        geos.add(geo); years.add(year)
        df = pd.read_csv(rp)
        gas.setdefault(geo, {}).setdefault(year, {})[grid] = _gas_price(project)
        bucket = data.setdefault(geo, {}).setdefault(year, {}).setdefault(grid, {})
        for _, row in df.iterrows():
            scen = str(row["scenario"])
            route = scen[:-4] if scen.endswith(("-avg", "-p95")) else scen
            lcos = _num(row.get("lcos_eur_per_t"))
            if lcos <= 0:
                continue
            rec = bucket.get(route)
            cand = {"lcos": round(lcos, 1),
                    "steel_mt": round(_num(row.get("steel_produced_mt")), 4),
                    "ng_gwh": round(_num(row.get("ng_gwh_lhv")), 1),
                    "h2_share": round(_num(row.get("iron_from_h2_share")), 3),
                    "lcoe": _opt(row.get("lcoe_eur_per_mwh")),
                    "lcoh": _opt(row.get("lcoh_eur_per_mwh_lhv"))}
            if rec is None or cand["lcos"] < rec["lcos"]:  # keep cheapest variant per route
                bucket[route] = cand
    return data, gas, sorted(geos), sorted(years)


def main():
    projects = sorted(
        p for p in {r["project"] for r in _read_projects()}
        if PROJ_RE.match(p) and (RESULTS / f"report_{p}.csv").exists()
    )
    figs, template = capture_figures(projects)
    synth, gas, geos, years = load_synthesis(projects)

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
        # The pipeline titles set family="Titillium Web SemiBold" literally, so register
        # that as its own family (mapping every weight to the SemiBold face) or the
        # headline falls back to a placeholder.
        faces.append(face("Titillium Web SemiBold", 400, "TitilliumWeb-SemiBold.woff2"))
        faces.append(face("Titillium Web SemiBold", 600, "TitilliumWeb-SemiBold.woff2"))
    font_css = "".join(faces)

    plotly_js = (Path(__import__("plotly").__file__).parent / "package_data" / "plotly.min.js").read_text()

    payload = {
        "figs": figs, "template": template, "synth": synth, "gas": gas,
        "geos": geos, "years": years, "geo_names": GEO_NAMES,
        "clean_routes": CLEAN_ROUTES, "route_label": ROUTE_LABEL,
        "route_color": ROUTE_COLOR, "co2_t_per_mwh": CO2_T_PER_MWH,
        "h2_min": H2_MIN,
    }

    template_html = TEMPLATE_HTML.read_text()
    html = (template_html
            .replace("/*FONT_CSS*/", font_css)
            .replace("/*PLOTLY_JS*/", plotly_js)
            .replace("/*PAYLOAD_JSON*/", json.dumps(payload, cls=PlotlyJSONEncoder)))
    OUT.write_text(html)
    print(f"wrote {OUT} ({OUT.stat().st_size/1e6:.2f} MB) — {len(figs)} projects, {len(geos)} geos")


def _read_projects():
    import csv
    with open(REPO / "config" / "projects.csv") as f:
        return list(csv.DictReader(f))


if __name__ == "__main__":
    main()
