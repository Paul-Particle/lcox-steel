"""Stacked LCOS cost-breakdown bars, one bar per steel-route scenario.

Snakemake rule: plot_lcos_bars (in viz.smk).

The core route-comparison figure: each scenario's levelised cost of steel
(€/t) stacked by cost group (process capex, ore & consumables, renewables,
grid, electrolyser, storage, transmission), with the LCOS total on top of
each bar. Scenarios without a steel load (h2_only) are skipped — they have
no LCOS; requesting this plot for a project with no steel scenarios at all
is an error.
"""

import logging
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.viz.style import (
    apply_header,
    blue_black,
    dark_gray,
    fca_blue,
    fca_template,
    gray,
    green,
    light_blue,
    magenta_red,
    sand_yellow,
    save_figure,
    turquois,
    very_dark_gray,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

_REPORT_PATH = Path(snakemake.input.report)
_OUT = Path(snakemake.output[0])

# Stack order (bottom → top) with legend label and color. Keys match the
# cost_{group}_meur columns written by compile_report.
COST_GROUPS = [
    ("process",         "Process plant (capex+opex)",   blue_black),
    ("ore_consumables", "Ore & consumables",            sand_yellow),
    ("gas",             "Natural gas (fuel + CO₂)",     very_dark_gray),
    ("iron_store",      "Iron stockpile",               gray),
    ("electrolyser",    "Electrolyser",                 green),
    ("h2_buffer",       "H₂ buffer",                    light_blue),
    ("res",             "Renewables (capex+opex)",      fca_blue),
    ("battery",         "Battery",                      magenta_red),
    ("transmission",    "Transmission (HVDC)",          turquois),
    ("grid",            "Grid (connection + energy)",   dark_gray),
]


def build_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    """Per-scenario cost-group columns in €/t steel, indexed by scenario label.

    Keeps only scenarios with an LCOS (steel routes); converts each
    cost_{group}_meur column from M€/yr to €/t via the scenario's steel output.
    """
    steel = df[df["lcos_eur_per_t"].notna()] if "lcos_eur_per_t" in df.columns else df.iloc[:0]
    skipped = len(df) - len(steel)
    if skipped:
        log.info(f"skipping {skipped} scenario(s) without a steel load (no LCOS)")
    if steel.empty:
        raise ValueError(
            f"{_REPORT_PATH.name} has no steel-route scenarios — "
            "plot_lcos_bars only applies to projects with route != h2_only"
        )

    steel_t_per_year = steel["steel_produced_mt"] * 1e6
    plot_df = pd.DataFrame({"label": steel["scenario"].astype(str)})
    for group, _, _ in COST_GROUPS:
        col = f"cost_{group}_meur"
        if col in steel.columns:
            plot_df[group] = steel[col].fillna(0.0) * 1e6 / steel_t_per_year.values
    plot_df["lcos_total"] = steel["lcos_eur_per_t"].values
    return plot_df.set_index("label")


def plot(plot_df: pd.DataFrame, out: Path, project_label: str) -> None:
    """Assemble the stacked LCOS bar chart and write it to PNG + HTML."""
    fig = go.Figure()
    for group, label, color in COST_GROUPS:
        if group not in plot_df.columns or (plot_df[group] == 0).all():
            continue
        fig.add_trace(go.Bar(
            x=plot_df.index,
            y=plot_df[group],
            name=label,
            marker_color=color,
            hovertemplate=f"{label}: %{{y:.0f}} €/t<extra></extra>",
        ))

    fig.add_trace(go.Scatter(
        x=plot_df.index,
        y=plot_df["lcos_total"],
        mode="text",
        text=[f"{v:,.0f}" for v in plot_df["lcos_total"]],
        textposition="top center",
        textfont=dict(size=13, color=blue_black),
        showlegend=False,
        hoverinfo="skip",
    ))

    n_sc = len(plot_df)
    fig.update_layout(
        template=fca_template,
        barmode="stack",
        bargap=0.45,
        # No rotated y-axis title — the unit lives in the subtitle.
        yaxis_title=None,
        xaxis_title=None,
        legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top", traceorder="reversed"),
    )
    fig.update_xaxes(type="category")
    fig.update_yaxes(rangemode="tozero")

    apply_header(
        fig,
        title=f"{project_label} levelised cost of steel",
        subtitle="€/t liquid steel, by cost group; label = LCOS total",
        fig_width=max(720, 200 * n_sc + 320), fig_height=600,
        margin_l=80, margin_r=280, margin_t=110, margin_b=80,
    )
    saved = save_figure(fig, out.parent, out.stem)
    log.info(f"saved {' + '.join(saved)}")


def main() -> None:
    """Load the project report and render its LCOS breakdown chart."""
    df = pd.read_csv(_REPORT_PATH)
    project_label = ", ".join(dict.fromkeys(df["project"].astype(str)))
    plot(build_plot_data(df), _OUT, project_label)


if __name__ == "__main__":
    main()
