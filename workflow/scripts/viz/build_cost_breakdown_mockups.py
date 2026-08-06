#!/usr/bin/env python3
"""Cost-breakdown taxonomy drafts, both drawn from the model's own results.

Builds `results/cost_breakdown_mockups.html` — a body-only, theme-aware,
self-contained Plotly page that puts the circulated chart mockup's *category
split* next to the split the pipeline already reports, with **both** populated
from the same DE 2023 results the scenario-comparison tab reads.

  * Draft A — the circulated taxonomy: one "Hydrogen (all-in)" block, and
    electricity divided into an EAF-melt share and everything else.
  * Draft B — the pipeline's own cost groups, unchanged.

Both stacks close exactly on LCOS, so the drafts differ only in how the same
total is cut. That is the point: Draft A's hydrogen block turns out to swallow
most of the H2 route's electricity, which is visible only once real numbers are
in it.

Reconstructing Draft A needs the electricity-system cost partitioned three ways.
Two parts are attributed and the third is a residual, so the stack cannot drift
from LCOS:

  * H2 share      = LCOH's electricity component x the H2 actually produced
  * EAF melt share = the EAF's own MWh/t at that scenario's LCOE
  * rest of plant  = electricity system total - the two above

LCOE and LCOH get a single chart each: there the two taxonomies coincide apart
from sub-splits (renewables capex vs O&M, electrolyser capex vs O&M vs water)
that the reports do not separate, so a second identical panel would say nothing.

Output is a hub fragment: body-only, its own <style>/<script>, plotly.js and the
brand font inlined. `build_dashboard_hub.py` embeds it as the "Cost breakdown
drafts" tab. Nothing here is a pipeline rule — run it directly.
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
sys.path.insert(0, str(REPO / "workflow"))        # common.*, scripts.viz.*

from build_dashboard import COST_GROUPS, RESULTS, ROUTE_LABEL, font_css  # noqa: E402
from common._constants import H2_LHV_KWH_PER_KG                          # noqa: E402
from scripts.viz.style import PLOTLY_CONFIG, fca_template, lighten       # noqa: E402

OUT_PATH = RESULTS / "cost_breakdown_mockups.html"
CONFIG_DIR = REPO / "config"

# The case both drafts are drawn from — the same reports the scenario tab reads.
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

# ---- One colour per cost role, shared by both drafts ----------------------
# Draft B's LCOS colours come from build_dashboard.COST_GROUPS; these mirror the
# same palette so equivalent roles read alike across the A|B pair.
C_ORE = "#E2B681"        # sand yellow    — ore & consumables
C_CAPEX = "#33434D"      # blue black     — plant capex
C_HYDROGEN = "#91C096"   # green          — hydrogen / electrolyser
C_H2_STORE = "#70D2F0"   # light blue     — H2 buffer
C_ELEC = "#0A5680"       # fca blue       — electricity / renewables
# The EAF melt share is a *slice of* the electricity bundle, so it reads as a tint
# of the same hue rather than a different colour. Plain highlight-blue (#0293D2)
# sits too close to fca-blue to separate the two bands in a stacked bar.
C_ELEC_EAF = lighten(C_ELEC, 0.55)
C_GAS = "#525F6A"        # very dark gray — natural gas + CO2
C_STORE = "#D75674"      # magenta red    — storage (iron/steel, battery)
C_GRID_CONN = "#71828F"  # dark gray      — grid connection
C_GRID_NRG = "#B7C1C8"   # gray           — grid energy
C_TRANSM = "#83D1DD"     # turquois       — HVDC transmission


# ---- Reading the pipeline reports ----------------------------------------

def _report(project: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"report_{project}.csv").set_index("scenario")


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


def _split_electricity(row: pd.Series, scenario: str, steel_t: float) -> dict:
    """Partition the electricity-system cost (€/t steel) the way Draft A wants it.

    Returns the H2 share, the EAF melt share and the remainder. The remainder is
    a residual, so the three always sum back to the electricity-system total and
    Draft A's stack closes on LCOS.
    """
    total = sum(_eur_per_t(row, g, steel_t)
                for g in ("res", "grid", "battery", "transmission"))
    lcoe = float(row["lcoe_eur_per_mwh"]) if pd.notna(row.get("lcoe_eur_per_mwh")) else 0.0

    h2_kt = float(row["h2_produced_kt"]) if pd.notna(row.get("h2_produced_kt")) else 0.0
    lcoh_electricity = (float(row["lcoh_electricity_eur_per_mwh_lhv"])
                        if pd.notna(row.get("lcoh_electricity_eur_per_mwh_lhv")) else 0.0)
    h2_mwh_lhv = h2_kt * 1e6 * H2_LHV_KWH_PER_KG / 1000.0
    hydrogen = lcoh_electricity * h2_mwh_lhv / steel_t

    # The EAF melts only on routes that build one; MOE pours liquid steel directly
    # and its assumptions file carries no `eaf` block at all.
    eaf_built = pd.notna(row.get("eaf_t_per_h_opt")) and float(row.get("eaf_t_per_h_opt", 0)) > 0
    eaf = (_assumptions(PROJECT, scenario)["eaf"]["el_mwh_per_t"] * lcoe) if eaf_built else 0.0

    rest = total - hydrogen - eaf
    if rest < -0.01:
        raise ValueError(
            f"{scenario}: attributing {hydrogen:.1f} €/t to hydrogen and {eaf:.1f} €/t to the "
            f"EAF over-claims an electricity system of {total:.1f} €/t. The EAF melt share is "
            "reconstructed from el_mwh_per_t x LCOE, so it is an estimate — revisit it before "
            "trusting this split."
        )
    return dict(hydrogen=hydrogen, eaf=eaf, rest=max(rest, 0.0), total=total, lcoe=lcoe)


def draft_a_lcos():
    """The circulated taxonomy, reconstructed from the model's own cost groups."""
    df = _report(PROJECT)
    rows = df.loc[LCOS_SCENARIOS]
    bars = [ROUTE_LABEL[s.removesuffix("-avg")] for s in LCOS_SCENARIOS]

    columns = {k: [] for k in ("ore", "capex", "hydrogen", "eaf", "rest", "gas", "store")}
    for scenario in LCOS_SCENARIOS:
        row = df.loc[scenario]
        steel_t = float(row["steel_produced_mt"]) * 1e6
        electricity = _split_electricity(row, scenario, steel_t)
        columns["ore"].append(_eur_per_t(row, "ore_consumables", steel_t))
        columns["capex"].append(_eur_per_t(row, "process", steel_t))
        columns["hydrogen"].append(
            _eur_per_t(row, "electrolyser", steel_t)
            + _eur_per_t(row, "h2_buffer", steel_t)
            + electricity["hydrogen"])
        columns["eaf"].append(electricity["eaf"])
        columns["rest"].append(electricity["rest"])
        columns["gas"].append(_eur_per_t(row, "gas", steel_t))
        columns["store"].append(_eur_per_t(row, "iron_store", steel_t)
                                + _eur_per_t(row, "steel_store", steel_t))

    series = [
        ("Ore & consumables",           C_ORE,      columns["ore"]),
        ("CAPEX (process plant)",       C_CAPEX,    columns["capex"]),
        ("Hydrogen (all-in)",           C_HYDROGEN, columns["hydrogen"]),
        ("Electricity — EAF melt",      C_ELEC_EAF, columns["eaf"]),
        ("Electricity — rest of plant", C_ELEC,     columns["rest"]),
        ("Gas + CO₂ (fossil routes)",   C_GAS,      columns["gas"]),
        ("Storage (iron/steel)",        C_STORE,    columns["store"]),
    ]
    return bars, series, rows["lcos_eur_per_t"].tolist()


def draft_b_lcos():
    """The pipeline's own LCOS cost groups, unchanged, for the same scenarios."""
    df = _report(PROJECT)
    rows = df.loc[LCOS_SCENARIOS]
    steel_t = rows["steel_produced_mt"] * 1e6
    bars = [ROUTE_LABEL[s.removesuffix("-avg")] for s in LCOS_SCENARIOS]
    series = []
    for group, label, colour in COST_GROUPS:
        col = f"cost_{group}_meur"
        if col not in rows.columns:
            continue
        values = (rows[col].fillna(0.0) * 1e6 / steel_t).tolist()
        if all(v <= 0.005 for v in values):
            continue
        series.append((label, colour, values))
    return bars, series, rows["lcos_eur_per_t"].tolist()


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


def levelised_parts(parts, scale: float = 1.0):
    """Stack one levelised metric across SUPPLY_CASES.

    `parts` is [(report column, legend label, colour)] bottom→top; absent columns
    and all-zero rows drop out, matching how the reports omit parts a scenario
    does not have.
    """
    frames = {project: _report(project) for project, _, _ in SUPPLY_CASES}
    bars = [label for _, _, label in SUPPLY_CASES]
    series = []
    for col, label, colour in parts:
        values = []
        for project, scenario, _ in SUPPLY_CASES:
            row = frames[project].loc[scenario]
            values.append(float(row[col]) * scale
                          if col in row.index and pd.notna(row.get(col)) else 0.0)
        if all(v <= 0 for v in values):
            continue
        series.append((label, colour, values))
    return bars, series


# ---- Figure assembly -----------------------------------------------------

def stacked_figure(bars, series, *, unit: str, value_fmt: str, total_fmt: str) -> go.Figure:
    """A stacked FCA bar chart; `series` is [(label, colour, values)] bottom→top.

    Totals ride above each bar as a text-only trace (the house pattern from
    plot_lcos_bars). The legend is horizontal below the plot so the figure keeps
    its width in a half-page panel.
    """
    fig = go.Figure()
    for label, colour, values in series:
        fig.add_trace(go.Bar(
            x=bars, y=values, name=label, marker_color=colour,
            hovertemplate=f"{label}: %{{y:{value_fmt}}} {unit}<extra></extra>",
        ))

    totals = [sum(values[i] for _, _, values in series) for i in range(len(bars))]
    fig.add_trace(go.Scatter(
        x=bars, y=totals, mode="text",
        text=[format(t, total_fmt) for t in totals],
        textposition="top center",
        textfont=dict(size=13, color=PRINT["ink"], family="Titillium Web"),
        showlegend=False, hoverinfo="skip", cliponaxis=False,
    ))

    fig.update_layout(
        template=responsive_template, barmode="stack", bargap=0.4,
        autosize=True, height=PLOT_HEIGHT,
        paper_bgcolor=PRINT["paper"], plot_bgcolor=PRINT["paper"],
        font=dict(family="Titillium Web", size=12, color=PRINT["ink"]),
        margin=dict(l=52, r=16, t=28, b=104),
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

DRAFT_DELTAS = [
    ("Ore is the biggest single line, and bigger than the draft assumed",
     "The draft put feedstock at 220 €/t for H2-DRI and 205 for NG-DRI; the model reports "
     "<b>260 €/t</b> of ore &amp; consumables for both. It also had the MOE and electrowinning "
     "figures crossed — <code>assumptions.yaml</code> has <code>ore_eur_per_t</code> 150 for "
     "<code>moe</code> and 160 for <code>electrowinning</code>. Real ore &amp; consumables: "
     "<b>160</b> for MOE, <b>205</b> for electrowinning."),
    ("Storage is ~70x smaller than the draft showed",
     "The draft gave storage 6–8 €/t. The model's iron and steel stores come to "
     "<b>0.01–0.12 €/t</b> — invisible on a 700 €/t bar, which is the finding rather than a "
     "rendering problem. They are what makes turndown feasible, at essentially no cost."),
    ("The off-grid p95 and mean-site bars were inverted",
     "The draft priced p95 above the mean (86 vs 61 €/MWh). <code>bestsite-p95</code> is the "
     "95th-percentile area-weighted capacity-factor cell — a <em>better</em> site — so it must "
     "come out cheaper. The model agrees: <b>59.3 €/MWh at p95 vs 69.7 at the mean site</b>."),
    ("Process capex and O&M cannot be separated",
     "The draft splits the process plant into CAPEX and O&amp;M rows. The reports give "
     "<code>cost_process_meur</code> as annualised plant capital only, with the variable side "
     "living in <code>cost_ore_consumables_meur</code>. The two rows are merged here; splitting "
     "them would need a new reported quantity."),
    ("Two rows have no data at all",
     "Ore/iron transport and scrap/HBI purchase are not in the model, so they are absent rather "
     "than shown at a placeholder value. Freight interacts with the ore-grade work — lower-grade "
     "ore moves more tonnes per tonne of iron."),
]

OPEN_QUESTIONS = [
    ("The hydrogen block swallows most of the H₂ route's electricity",
     "This is what the real numbers make undeniable. H2-DRI-EAF's “Hydrogen (all-in)” segment is "
     "<b>290 €/t, of which 207 €/t is electricity</b>. The two visible electricity segments add to "
     "just 71 €/t, so Draft A appears to show a route using a quarter of MOE's power — when the "
     "true electricity totals are <b>278 €/t for H2-DRI against 400 for MOE</b>. Draft B shows the "
     "same total with renewables, electrolyser and H₂ buffer separate. If Draft A's form factor is "
     "kept, the hydrogen block needs to be visibly nested or hatched."),
    ("Pricing the EAF at LCOE makes an identical furnace look route-dependent",
     "Every route here melts at the same <b>0.65 MWh/t</b>. Valued at each scenario's own LCOE, "
     "that same furnace reads <b>44.6 €/t on H2-DRI and 70.5 €/t on NG-DRI</b>, purely because "
     "NG-DRI's small electricity system carries a higher LCOE (108 vs 69 €/MWh). The convention is "
     "internally consistent, but the EAF segment is not comparable across bars."),
    ("Average LCOE hides the flexibility the optimiser is actually buying",
     "Both drafts price electricity at system LCOE. That is the right average-cost basis for a "
     "levelised total, but it says nothing about <em>which hours</em> a route chooses to run — the "
     "mechanism behind the turndown work. Worth deciding whether this chart answers “what does it "
     "cost” or “why does the optimiser prefer this route”; those want different charts."),
    ("Should the tab follow the scenario picker?",
     "Both drafts are pinned to DE 2023, grid-connected, mean-CF site — one case out of the "
     "scenario tab's full matrix. Wiring the same geography/year/CF controls onto this tab is a "
     "straightforward next step; it was left out so the taxonomy comparison stays the subject."),
]

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
             "the same fact that Draft A's hydrogen block hides above.",
    ),
]

CSS = """
  .cbm-root{
    --paper:#EAEEF1; --panel:#FFFFFF; --panel-2:#F4F6F8;
    --ink:#1B242B; --muted:#525F6A; --faint:#8CA5B7;
    --border:#D6DBDF; --hair:#E4E9ED;
    --accent:#0A5680; --accent-2:#0293D2;
    --fix:#4E8C6A; --fix-bg:rgba(78,140,106,.09); --fix-line:rgba(78,140,106,.30);
    --ask:#C08A3A; --ask-bg:rgba(192,138,58,.09); --ask-line:rgba(192,138,58,.30);
    --shadow:0 1px 2px rgba(27,36,43,.05),0 10px 28px rgba(27,36,43,.07);
    --tw:"Titillium Web",system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
    --mono:"SF Mono",ui-monospace,"Cascadia Mono","Roboto Mono",Menlo,Consolas,monospace;
  }
  @media (prefers-color-scheme:dark){
    .cbm-root{--paper:#10171C; --panel:#1A232A; --panel-2:#212C34;
      --ink:#E9EEF1; --muted:#A7B6C2; --faint:#6C7E8C; --border:#2B3841; --hair:#26323A;
      --accent:#0293D2; --accent-2:#70D2F0;
      --fix:#91C096; --fix-bg:rgba(145,192,150,.13); --fix-line:rgba(145,192,150,.34);
      --ask:#E2B681; --ask-bg:rgba(226,182,129,.13); --ask-line:rgba(226,182,129,.34);
      --shadow:0 1px 2px rgba(0,0,0,.35),0 12px 34px rgba(0,0,0,.42);}
  }
  :root[data-theme="dark"] .cbm-root{--paper:#10171C; --panel:#1A232A; --panel-2:#212C34;
    --ink:#E9EEF1; --muted:#A7B6C2; --faint:#6C7E8C; --border:#2B3841; --hair:#26323A;
    --accent:#0293D2; --accent-2:#70D2F0;
    --fix:#91C096; --fix-bg:rgba(145,192,150,.13); --fix-line:rgba(145,192,150,.34);
    --ask:#E2B681; --ask-bg:rgba(226,182,129,.13); --ask-line:rgba(226,182,129,.34);
    --shadow:0 1px 2px rgba(0,0,0,.35),0 12px 34px rgba(0,0,0,.42);}
  :root[data-theme="light"] .cbm-root{--paper:#EAEEF1; --panel:#FFFFFF; --panel-2:#F4F6F8;
    --ink:#1B242B; --muted:#525F6A; --faint:#8CA5B7; --border:#D6DBDF; --hair:#E4E9ED;
    --accent:#0A5680; --accent-2:#0293D2;
    --fix:#4E8C6A; --fix-bg:rgba(78,140,106,.09); --fix-line:rgba(78,140,106,.30);
    --ask:#C08A3A; --ask-bg:rgba(192,138,58,.09); --ask-line:rgba(192,138,58,.30);
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
  .cbm-key.a .chip{background:var(--accent);color:#fff;}
  .cbm-key.b .chip{background:var(--accent-2);color:#062430;}
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
  .cbm-col.a .chip{background:var(--accent);color:#fff;}
  .cbm-col.b .chip{background:var(--accent-2);color:#062430;}
  .cbm-col-h .nm{font-size:13px;font-weight:700;letter-spacing:-.01em;}
  .cbm-col-h .unit{font-size:11.5px;font-weight:400;color:var(--muted);letter-spacing:0;}
  .cbm-col-h .src{font-family:var(--mono);font-size:10px;color:var(--faint);
    margin-left:auto;white-space:nowrap;}
  .cbm-note{font-size:12px;color:var(--muted);margin:8px 2px 0;line-height:1.45;}

  .cbm-call{border:1px solid var(--border);border-radius:12px;margin-top:22px;
    padding:16px 18px 18px 20px;position:relative;background:var(--panel);
    box-shadow:var(--shadow);overflow:hidden;}
  .cbm-call::before{content:"";position:absolute;left:0;top:0;bottom:0;width:4px;}
  .cbm-call.fix{background:var(--fix-bg);border-color:var(--fix-line);}
  .cbm-call.fix::before{background:var(--fix);}
  .cbm-call.ask{background:var(--ask-bg);border-color:var(--ask-line);}
  .cbm-call.ask::before{background:var(--ask);}
  .cbm-call h2{font-size:13px;font-weight:800;letter-spacing:.09em;
    text-transform:uppercase;margin:0 0 2px;}
  .cbm-call.fix h2{color:var(--fix);}
  .cbm-call.ask h2{color:var(--ask);}
  .cbm-call > p{font-size:12.5px;color:var(--muted);margin:0 0 12px;max-width:94ch;}
  .cbm-call ol{margin:0;padding-left:22px;display:flex;flex-direction:column;gap:11px;}
  .cbm-call li{font-size:13px;line-height:1.5;}
  .cbm-call li b{font-weight:700;}
  .cbm-call li span{display:block;color:var(--muted);font-size:12.5px;margin-top:2px;}
  .cbm-root code{font-family:var(--mono);font-size:.88em;background:var(--panel-2);
    border:1px solid var(--hair);border-radius:4px;padding:0 4px;}
  .cbm-foot{margin-top:26px;padding-top:14px;border-top:1px solid var(--border);
    font-size:11.5px;color:var(--faint);line-height:1.6;}
"""


def _column(side: str, name: str, src: str, div: str, note: str) -> str:
    """One half of the Draft A | Draft B pair."""
    return f"""      <div class="cbm-col {side}">
        <div class="cbm-col-h"><span class="chip">Draft {side.upper()}</span>
          <span class="nm">{name}</span><span class="src">{src}</span></div>
        {div}
        <p class="cbm-note">{note}</p>
      </div>"""


def _metric_column(spec: dict, div: str) -> str:
    """One half of the LCOE | LCOH pair — a metric, not a draft, so no chip."""
    return f"""      <div class="cbm-col">
        <div class="cbm-col-h"><span class="nm">{spec['title']}</span>
          <span class="unit">· {spec['unit']}</span><span class="src">{spec['src']}</span></div>
        {div}
        <p class="cbm-note">{spec['note']}</p>
      </div>"""


def _ordered_list(items) -> str:
    return "\n".join(
        f"      <li><b>{head}</b><span>{body}</span></li>" for head, body in items
    )


def build_core() -> str:
    """The body-only page content, ready to embed as a hub tab."""
    a_bars, a_series, a_lcos = draft_a_lcos()
    b_bars, b_series, b_lcos = draft_b_lcos()

    # Both cuts of the same total must land on the reported LCOS, or the page is
    # quietly lying about one of them.
    for label, series, lcos in (("A", a_series, a_lcos), ("B", b_series, b_lcos)):
        totals = [sum(v[i] for _, _, v in series) for i in range(len(lcos))]
        drift = max(abs(t - l) for t, l in zip(totals, lcos))
        if drift > 0.05:
            raise ValueError(f"Draft {label} stack drifts {drift:.3f} €/t from reported LCOS")

    lcoe_bars, lcoe_series = levelised_parts(LCOE_PARTS)
    lcoh_bars, lcoh_series = levelised_parts(
        LCOH_PARTS, scale=H2_LHV_KWH_PER_KG / 1000.0)

    divs = {
        "lcos-a": plot_div(stacked_figure(a_bars, a_series, unit="€/t",
                                          value_fmt=",.1f", total_fmt=",.0f"), "cbm-lcos-a"),
        "lcos-b": plot_div(stacked_figure(b_bars, b_series, unit="€/t",
                                          value_fmt=",.1f", total_fmt=",.0f"), "cbm-lcos-b"),
        "lcoe": plot_div(stacked_figure(lcoe_bars, lcoe_series, unit="€/MWh",
                                        value_fmt=",.1f", total_fmt=",.1f"), "cbm-lcoe"),
        "lcoh": plot_div(stacked_figure(lcoh_bars, lcoh_series, unit="€/kg",
                                        value_fmt=",.2f", total_fmt=",.2f"), "cbm-lcoh"),
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
        <div class="cbm-brand"><span class="cbm-dot"></span><span>Cost breakdown · taxonomy drafts</span></div>
        <h1 class="cbm-title">Two ways to cut the same number</h1>
        <p class="cbm-sub">The circulated chart mockup's category split next to the split the
          pipeline already reports — both drawn from the same DE 2023 results the scenario tab
          reads, and both closing exactly on LCOS. The drafts differ only in where they put the
          cuts, which is what makes the disagreement legible.</p>
      </div>
      <div class="cbm-meta">built {date.today().strftime('%-d %b %Y')}<br>
        DE 2023 · grid-connected<br>mean-CF site</div>
    </header>

    <div class="cbm-legend">
      <div class="cbm-key a"><span class="chip">Draft A</span>
        <p>The circulated taxonomy: one hydrogen block, electricity split into an EAF-melt share
          and the rest. Rebuilt from the model's cost groups — no placeholder numbers left.</p></div>
      <div class="cbm-key b"><span class="chip">Draft B</span>
        <p>The pipeline's own cost groups, unchanged. Same bars, same totals, different cuts —
          renewables, electrolyser and H₂ buffer each stand on their own.</p></div>
    </div>

    <section class="cbm-panel">
      <div class="cbm-panel-h">
        <h2>LCOS — levelised cost of steel <span class="unit">· € / t liquid steel</span></h2>
        <p>Where the two taxonomies genuinely disagree. Both stacks sum to the same LCOS per route;
          Draft A folds the electrolyser, the H₂ buffer and the electricity that made the hydrogen
          into a single block, while Draft B keeps them apart.</p>
      </div>
      <div class="cbm-pair">
{_column('a', 'circulated split', 'reconstructed', divs['lcos-a'],
         "Hydrogen (all-in) = electrolyser + H₂ buffer + the electricity that made the H₂. "
         "The EAF melt share is the furnace's own MWh/t valued at that scenario's LCOE; "
         "&ldquo;rest of plant&rdquo; is the residual, so the stack closes on LCOS.")}
{_column('b', "model's own split", 'report columns', divs['lcos-b'],
         "Straight from the <code>cost_*_meur</code> columns, converted to €/t. The iron "
         "stockpile and steel inventory keep their legend entries but draw no visible band — "
         "both sit near 0.05 €/t.")}
      </div>
    </section>

    <section class="cbm-panel">
      <div class="cbm-panel-h">
        <h2>The carriers underneath <span class="unit">· LCOE and LCOH</span></h2>
        <p>Here the two taxonomies nearly agree, so there is one chart of each rather than a pair.
          The circulated split adds a renewables capex/O&amp;M row to LCOE, and electrolyser-O&amp;M
          plus water/variable-opex rows to LCOH; the reports fold those into the parent line, so
          they cannot be filled without a new reported quantity. Everything else maps one-to-one.
          Both charts show the same three supply cases, and the LCOH on the right is built on the
          LCOE on the left.</p>
      </div>
      <div class="cbm-pair">
{carrier_columns}
      </div>
    </section>

    <section class="cbm-call fix">
      <h2>What changed when the draft met real numbers</h2>
      <p>The circulated version carried illustrative figures. Putting the model's own results into
        the same categories moved several of them a long way.</p>
      <ol>
{_ordered_list(DRAFT_DELTAS)}
      </ol>
    </section>

    <section class="cbm-call ask">
      <h2>Open presentation questions</h2>
      <p>Deliberately not resolved here — each one changes what the chart argues, not just how it
        looks.</p>
      <ol>
{_ordered_list(OPEN_QUESTIONS)}
      </ol>
    </section>

    <p class="cbm-foot">Cross-chart logic, unchanged from the draft: the electricity segment in
      LCOS/LCOH is priced at LCOE; the hydrogen segment in LCOS is LCOH × kg&nbsp;H₂/t; LCOI (not
      shown) is LCOS minus the EAF stage. Sources:
      <code>results/report_DE-2023-grid.csv</code>,
      <code>results/report_DE-2023-nogrid.csv</code> and the per-scenario
      <code>config/assumptions_DE-2023-grid_*.yaml</code>. Regenerate with
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
