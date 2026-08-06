#!/usr/bin/env python3
"""Cost-breakdown taxonomy drafts: the circulated chart mockup beside the model's own split.

Builds `results/cost_breakdown_mockups.html` — a body-only, theme-aware,
self-contained Plotly page that shows the same three levelised metrics (LCOS,
LCOE, LCOH) under two competing category taxonomies:

  * Draft A — the split from the circulated mockup, on its own illustrative
    numbers. Reproduced here in the FCA house style so the two drafts can be
    compared on taxonomy rather than on styling. Three arithmetic/labelling
    errors in the original are corrected (see CORRECTIONS below); the remaining
    open presentation questions are listed on the page rather than silently
    resolved.
  * Draft B — the taxonomy the model actually emits, on real DE 2023 results
    read from the pipeline reports.

Both drafts share one colour per cost role, so a segment keeps its meaning
across the two panels.

Output is a hub fragment: body-only, its own <style>/<script>, plotly.js and the
brand font inlined. `build_dashboard_hub.py` embeds it as the "Cost breakdown"
tab. Nothing here is a pipeline rule — run it directly.
"""
import copy
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

REPO = Path(__file__).resolve().parents[3]        # workflow/scripts/viz/ -> repo root
sys.path.insert(0, str(Path(__file__).parent))    # sibling build_dashboard
sys.path.insert(0, str(REPO / "workflow"))        # common.*, scripts.viz.*

from build_dashboard import COST_GROUPS, RESULTS, ROUTE_LABEL, font_css  # noqa: E402
from common._constants import H2_LHV_KWH_PER_KG                          # noqa: E402
from scripts.viz.style import PLOTLY_CONFIG, fca_template                # noqa: E402

OUT_PATH = RESULTS / "cost_breakdown_mockups.html"

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
C_TRANSPORT = "#BDCCD9"  # light blue-gray — inbound/outbound freight
C_OM = "#8CA5B7"         # blue gray      — operations & maintenance
C_CAPEX = "#33434D"      # blue black     — plant capex
C_HYDROGEN = "#91C096"   # green          — hydrogen / electrolyser
C_H2_STORE = "#70D2F0"   # light blue     — H2 buffer
C_ELEC = "#0A5680"       # fca blue       — electricity / renewables
C_ELEC_EAF = "#0293D2"   # highlight blue — the EAF melt share of electricity
C_GAS = "#525F6A"        # very dark gray — natural gas + CO2
C_STORE = "#D75674"      # magenta red    — storage (iron/steel, battery)
C_GRID_CONN = "#71828F"  # dark gray      — grid connection
C_GRID_NRG = "#B7C1C8"   # gray           — grid energy
C_TRANSM = "#83D1DD"     # turquois       — HVDC transmission

# ---- Draft A: the circulated mockup, corrected ---------------------------
# Values are the mockup's own illustrative figures. Three corrections are applied
# relative to the circulated version; each is listed in CORRECTIONS and rendered
# on the page so the change is visible rather than silent.
ROUTES_A = ["H2-DRI-EAF", "MOE", "EW", "NG-DRI-EAF"]
DRAFT_A_LCOS = [
    ("Ore & consumables",           C_ORE,       [220, 150, 160, 205]),
    ("Transport (ore/iron)*",       C_TRANSPORT, [15, 12, 12, 15]),
    ("O&M (process plant)",         C_OM,        [25, 30, 40, 22]),
    ("CAPEX (process plant)",       C_CAPEX,     [60, 95, 130, 60]),
    ("Hydrogen (all-in)",           C_HYDROGEN,  [260, 0, 0, 0]),
    ("Electricity — EAF melt",      C_ELEC_EAF,  [55, 0, 55, 55]),
    ("Electricity — rest of plant", C_ELEC,      [35, 380, 245, 40]),
    ("Gas + CO₂ (fossil routes)",   C_GAS,       [0, 0, 0, 150]),
    ("Storage (iron/steel)",        C_STORE,     [8, 6, 8, 7]),
]

BARS_LCOE = ["Grid-connected", "Off-grid · p95 site", "Off-grid · mean site"]
DRAFT_A_LCOE = [
    ("RES capex",           C_ELEC,      [30, 40, 55]),
    ("RES O&M",             C_OM,        [5, 6, 8]),
    ("Battery / storage",   C_STORE,     [8, 12, 18]),
    ("Grid connection",     C_GRID_CONN, [6, 0, 0]),
    ("Grid energy",         C_GRID_NRG,  [20, 0, 0]),
    ("Transmission (HVDC)", C_TRANSM,    [0, 3, 5]),
]

BARS_LCOH_A = ["Off-grid green H₂", "Grid-mix green H₂"]
DRAFT_A_LCOH = [
    ("Electrolyser CAPEX",   C_CAPEX,    [1.2, 1.0]),
    ("Electrolyser O&M",     C_OM,       [0.3, 0.25]),
    ("Water / var-opex",     C_ORE,      [0.1, 0.1]),
    ("H₂ storage (buffer)",  C_H2_STORE, [0.4, 0.2]),
    ("Electricity (@ LCOE)", C_ELEC,     [3.0, 2.2]),
]

CORRECTIONS = [
    ("Ore cost swapped between MOE and electrowinning",
     "The draft put 160 €/t on MOE and 150 €/t on EW. "
     "<code>config/assumptions.yaml</code> has <code>ore_eur_per_t</code> = 150 for "
     "<code>moe</code> and 160 for <code>electrowinning</code> — the two were crossed. "
     "Draft A now reads 150 / 160."),
    ("The starred feedstock row over-claimed what is missing",
     "“Feedstock (ore/scrap/HBI)*” marked the whole row as absent from the model, "
     "but ore <em>is</em> modelled — it is the <code>ore_consumables</code> cost group. "
     "The row is now “Ore &amp; consumables”, unstarred; the genuinely missing "
     "scrap/HBI purchase is called out under the open questions instead."),
    ("Off-grid p95 and mean sites were the wrong way round",
     "The draft priced the p95 bar above the mean bar (86 vs 61 €/MWh). In this model "
     "<code>bestsite-p95</code> is the 95th-percentile area-weighted capacity-factor "
     "cell — a <em>good</em> site — so it must come out cheaper. Real DE 2023 off-grid: "
     "59.3 €/MWh at p95 vs 69.7 at the mean site. The two bars are swapped."),
]

OPEN_QUESTIONS = [
    ("The hydrogen segment hides that H₂-DRI is also an electricity route",
     "Draft A defines “Hydrogen (all-in)” as electrolyser capex+opex plus H₂ storage "
     "<em>plus the electricity to make the H₂</em>. H2-DRI-EAF then shows 90 €/t of "
     "electricity against MOE's 380, even though most of its 260 €/t hydrogen segment "
     "is electricity too. Read off the bar heights, that says MOE is four times as "
     "power-hungry — the opposite of the point. Options: nest the electricity inside "
     "the hydrogen segment with a hatch, or add a second view stacked by physical "
     "driver (electricity / ore / plant / fuel) alongside the contractual one."),
    ("Storage is invisible at this scale, and that is itself the finding",
     "Draft A gives storage 6–8 €/t; the model's real iron and steel stores are "
     "0.01–0.08 M€/yr, i.e. well under 0.1 €/t. Both round to nothing on a 700 €/t "
     "bar. But the stores are what make turndown feasible — they are load-bearing at "
     "near-zero cost. A stacked bar cannot say that; a footnote or a separate "
     "“enabling capacity” panel can."),
    ("What price should the electricity segment carry?",
     "Both drafts price it at system LCOE. That is the right average-cost basis for a "
     "levelised total, but it hides the flexibility value the optimiser is actually "
     "exploiting — the marginal hours a route chooses to run. Worth deciding whether "
     "the breakdown is meant to answer “what does it cost” or “why does the optimiser "
     "prefer this route”; those want different charts."),
    ("Two Draft A rows still need inputs the model does not have",
     "Ore/iron transport is shown with placeholder values and scrap/HBI purchase is "
     "absent altogether. Both need new data before they can appear on a real chart. "
     "Freight in particular interacts with the ore-grade work — a low-grade ore moves "
     "more tonnes per tonne of iron."),
    ("Draft A omits the NG-H2-DRI-EAF blend route",
     "The model solves five routes; the draft shows four. The blend route is the one "
     "that shows the optimiser rejecting hydrogen at current prices, so leaving it out "
     "removes the most argumentative bar on the chart."),
]


# ---- Real numbers from the pipeline reports ------------------------------

def _report(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name).set_index("scenario")


def draft_b_lcos() -> tuple[list[str], list[tuple[str, str, list[float]]]]:
    """The model's own LCOS split for DE 2023, grid-connected, mean-CF site.

    steel_produced_mt is 1.0 for these scenarios, so the report's M€/yr cost
    groups are already €/t; divide by it anyway rather than relying on that.
    """
    df = _report("report_DE-2023-grid.csv")
    scenarios = ["h2-dri-eaf-avg", "moe-avg", "ew-avg", "ng-dri-eaf", "mix-dri-eaf-avg"]
    rows = df.loc[scenarios]
    steel_t = rows["steel_produced_mt"] * 1e6
    bars = [ROUTE_LABEL[s.removesuffix("-avg")] for s in scenarios]
    series = []
    for group, label, colour in COST_GROUPS:
        col = f"cost_{group}_meur"
        if col not in rows.columns:
            continue
        values = (rows[col].fillna(0.0) * 1e6 / steel_t).tolist()
        if all(v <= 0.05 for v in values):
            continue
        series.append((label, colour, values))
    return bars, series


def _levelised_parts(cases, parts, scale: float = 1.0):
    """Stack one levelised metric across `cases` = [(report, scenario, bar label)].

    `parts` is [(report column, legend label, colour)] bottom→top; absent columns
    and all-zero rows drop out, matching how the pipeline reports omit parts a
    scenario does not have.
    """
    frames = {name: _report(name) for name in {c[0] for c in cases}}
    bars = [label for _, _, label in cases]
    series = []
    for col, label, colour in parts:
        values = []
        for name, scenario, _ in cases:
            row = frames[name].loc[scenario]
            values.append(float(row[col]) * scale if col in row.index and pd.notna(row.get(col)) else 0.0)
        if all(v <= 0 for v in values):
            continue
        series.append((label, colour, values))
    return bars, series


LCOE_CASES_B = [
    ("report_DE-2023-grid.csv",   "h2-dri-eaf-avg", "Grid-connected"),
    ("report_DE-2023-nogrid.csv", "h2-dri-eaf-p95", "Off-grid · p95 site"),
    ("report_DE-2023-nogrid.csv", "h2-dri-eaf-avg", "Off-grid · mean site"),
]
LCOE_PARTS_B = [
    ("lcoe_renewables_eur_per_mwh",      "Renewables (capex+opex)", C_ELEC),
    ("lcoe_storage_eur_per_mwh",         "Battery / storage",       C_STORE),
    ("lcoe_grid_connection_eur_per_mwh", "Grid connection",         C_GRID_CONN),
    ("lcoe_grid_energy_eur_per_mwh",     "Grid energy",             C_GRID_NRG),
    ("lcoe_transmission_eur_per_mwh",    "Transmission (HVDC)",     C_TRANSM),
]
LCOH_PARTS_B = [
    ("lcoh_electrolyser_eur_per_mwh_lhv", "Electrolyser (capex+opex)", C_CAPEX),
    ("lcoh_h2_storage_eur_per_mwh_lhv",   "H₂ storage (buffer)",       C_H2_STORE),
    ("lcoh_electricity_eur_per_mwh_lhv",  "Electricity (@ LCOE)",      C_ELEC),
]


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
        text=[format(t, total_fmt.lstrip(":")) for t in totals],
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


# ---- Page assembly -------------------------------------------------------

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
  .cbm-sub{color:var(--muted);font-size:14px;margin:9px 0 0;max-width:68ch;}
  .cbm-meta{font-family:var(--mono);font-size:11px;color:var(--faint);
    text-align:right;line-height:1.7;white-space:nowrap;}

  .cbm-legend{display:flex;gap:26px;flex-wrap:wrap;margin:22px 0 4px;}
  .cbm-key{display:flex;gap:10px;align-items:flex-start;max-width:44ch;}
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
  .cbm-panel-h p{font-size:12.5px;color:var(--muted);margin:4px 0 0;max-width:96ch;}

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
  .cbm-call > p{font-size:12.5px;color:var(--muted);margin:0 0 12px;max-width:92ch;}
  .cbm-call ol{margin:0;padding-left:22px;display:flex;flex-direction:column;gap:11px;}
  .cbm-call li{font-size:13px;line-height:1.5;}
  .cbm-call li b{font-weight:700;}
  .cbm-call li span{display:block;color:var(--muted);font-size:12.5px;margin-top:2px;}
  .cbm-root code{font-family:var(--mono);font-size:.88em;background:var(--panel-2);
    border:1px solid var(--hair);border-radius:4px;padding:0 4px;}
  .cbm-foot{margin-top:26px;padding-top:14px;border-top:1px solid var(--border);
    font-size:11.5px;color:var(--faint);line-height:1.6;}
"""

SECTIONS = [
    dict(
        key="lcos",
        title="LCOS — levelised cost of steel",
        unit="€ / t liquid steel",
        blurb="Draft A splits by commercial category and bundles all hydrogen cost into one "
              "segment. Draft B splits by the cost groups the optimiser actually reports, so "
              "the electrolyser, the H₂ buffer and the renewables that feed them stay separate.",
        note_a="Illustrative numbers from the circulated draft, three corrections applied.",
        note_b="Real result: DE 2023, grid-connected, mean-CF site.",
        src_a="draft figures",
        src_b="report_DE-2023-grid.csv",
    ),
    dict(
        key="lcoe",
        title="LCOE — levelised cost of electricity",
        unit="€ / MWh delivered",
        blurb="The two taxonomies nearly agree here. Draft A splits renewables into capex and "
              "O&M, which the reports do not currently separate; everything else maps one-to-one.",
        note_a="Illustrative; the p95 and mean-site bars have been swapped (see corrections).",
        note_b="Real result: DE 2023, H2-DRI-EAF. Transmission is absent — no HVDC in the DE cases.",
        src_a="draft figures",
        src_b="report_DE-2023-{grid,nogrid}.csv",
    ),
    dict(
        key="lcoh",
        title="LCOH — levelised cost of hydrogen",
        unit="€ / kg H₂",
        blurb="Draft A adds a water / variable-opex line the reports do not break out, and "
              "compares an off-grid case against a grid-mix one. Draft B carries the model's "
              "three parts across the same three supply cases as the LCOE panel.",
        note_a="Illustrative numbers from the circulated draft.",
        note_b=f"Real result: DE 2023, H2-DRI-EAF; converted at {H2_LHV_KWH_PER_KG} kWh/kg LHV.",
        src_a="draft figures",
        src_b="report_DE-2023-{grid,nogrid}.csv",
    ),
]


def _column(side: str, name: str, src: str, div: str, note: str) -> str:
    return f"""      <div class="cbm-col {side}">
        <div class="cbm-col-h"><span class="chip">Draft {side.upper()}</span>
          <span class="nm">{name}</span><span class="src">{src}</span></div>
        {div}
        <p class="cbm-note">{note}</p>
      </div>"""


def _ordered_list(items) -> str:
    return "\n".join(
        f"      <li><b>{head}</b><span>{body}</span></li>" for head, body in items
    )


def build_core() -> str:
    """The body-only page content, ready to embed as a hub tab."""
    b_lcos_bars, b_lcos_series = draft_b_lcos()
    b_lcoe_bars, b_lcoe_series = _levelised_parts(LCOE_CASES_B, LCOE_PARTS_B)
    b_lcoh_bars, b_lcoh_series = _levelised_parts(
        LCOE_CASES_B, LCOH_PARTS_B, scale=H2_LHV_KWH_PER_KG / 1000.0)

    figures = {
        ("lcos", "a"): stacked_figure(ROUTES_A, DRAFT_A_LCOS,
                                      unit="€/t", value_fmt=",.0f", total_fmt=",.0f"),
        ("lcos", "b"): stacked_figure(b_lcos_bars, b_lcos_series,
                                      unit="€/t", value_fmt=",.0f", total_fmt=",.0f"),
        ("lcoe", "a"): stacked_figure(BARS_LCOE, DRAFT_A_LCOE,
                                      unit="€/MWh", value_fmt=",.1f", total_fmt=",.0f"),
        ("lcoe", "b"): stacked_figure(b_lcoe_bars, b_lcoe_series,
                                      unit="€/MWh", value_fmt=",.1f", total_fmt=",.1f"),
        ("lcoh", "a"): stacked_figure(BARS_LCOH_A, DRAFT_A_LCOH,
                                      unit="€/kg", value_fmt=",.2f", total_fmt=",.2f"),
        ("lcoh", "b"): stacked_figure(b_lcoh_bars, b_lcoh_series,
                                      unit="€/kg", value_fmt=",.2f", total_fmt=",.2f"),
    }
    divs = {k: plot_div(f, f"cbm-{k[0]}-{k[1]}") for k, f in figures.items()}

    panels = []
    for section in SECTIONS:
        key = section["key"]
        panels.append(f"""    <section class="cbm-panel">
      <div class="cbm-panel-h">
        <h2>{section['title']} <span class="unit">· {section['unit']}</span></h2>
        <p>{section['blurb']}</p>
      </div>
      <div class="cbm-pair">
{_column('a', 'circulated split', section['src_a'], divs[(key, 'a')], section['note_a'])}
{_column('b', "model's own split", section['src_b'], divs[(key, 'b')], section['note_b'])}
      </div>
    </section>""")

    return f"""<style>
{font_css()}
{CSS}</style>

<div class="cbm-root">
  <div class="cbm-wrap">
    <header class="cbm-head">
      <div>
        <div class="cbm-brand"><span class="cbm-dot"></span><span>Cost breakdown · taxonomy drafts</span></div>
        <h1 class="cbm-title">Two ways to split a levelised cost</h1>
        <p class="cbm-sub">The circulated chart mockup redrawn in the house style, next to the
          category split the model already emits — same three metrics, same colour per cost role,
          so the comparison is about taxonomy rather than styling.</p>
      </div>
      <div class="cbm-meta">built {date.today().strftime('%-d %b %Y')}<br>
        draft A: illustrative<br>draft B: DE 2023 results</div>
    </header>

    <div class="cbm-legend">
      <div class="cbm-key a"><span class="chip">Draft A</span>
        <p>The circulated split, on its own illustrative numbers. Corrected where it was
          arithmetically or definitionally wrong; left alone where the question is one of
          presentation.</p></div>
      <div class="cbm-key b"><span class="chip">Draft B</span>
        <p>The taxonomy the pipeline reports today, on real DE 2023 output. Totals land within
          ~10% of Draft A's placeholders, so the bars are directly comparable.</p></div>
    </div>

{chr(10).join(panels)}

    <section class="cbm-call fix">
      <h2>Corrections applied to the circulated draft</h2>
      <p>Each of these was wrong against the model or its config rather than a matter of taste,
        so Draft A above shows the corrected version.</p>
      <ol>
{_ordered_list(CORRECTIONS)}
      </ol>
    </section>

    <section class="cbm-call ask">
      <h2>Open presentation questions</h2>
      <p>Deliberately not resolved here — each one changes what the chart argues, not just how
        it looks.</p>
      <ol>
{_ordered_list(OPEN_QUESTIONS)}
      </ol>
    </section>

    <p class="cbm-foot">Cross-chart logic, unchanged from the draft: the electricity segment in
      LCOS/LCOH is priced at LCOE; the hydrogen segment in LCOS is LCOH × kg&nbsp;H₂/t; LCOI (not
      shown) is LCOS minus the EAF stage. Segments marked <code>*</code> need inputs the model does
      not have yet. Regenerate with
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
