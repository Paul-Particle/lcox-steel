"""
Bar chart comparing scenario capacities across projects.

Snakemake rule: plot_capacity_bars  (in viz.smk)

Shows per scenario:
  DRI H2 demand (MW LHV) · Electrolyser (MW) · H2 buffer (MW·days LHV)
  Solar stacked by orientation (MW) · Wind stacked by tech (MW)
  Grid annual-average import (MW)
"""

import logging
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.viz._helpers import (
    PLOTLY_CONFIG,
    blue_black,
    dark_gray,
    fca_blue,
    fca_colormap,
    fca_template,
    green,
    magenta_red,
    sand_yellow,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

_REPORT_PATH = Path(snakemake.input.report)
_OUT = Path(snakemake.output[0])

# Stacked sub-bar palettes, sampled from the FCA colormap.
# Solar: high end of the colormap (golden/sand tones — close to sand_yellow).
# Wind:  low/mid end of the colormap (blue/teal tones — close to fca_blue).
SOLAR_PALETTE_RANGE = (0.78, 1.00)
WIND_PALETTE_RANGE  = (0.02, 0.40)

# Solo-trace colors when only one sub-bar is present.
SOLAR_SOLO_COLOR = sand_yellow
WIND_SOLO_COLOR  = fca_blue

# Slot positions (in "category widths") for the six grouped bars.
_SLOTS = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
_BAR_WIDTH = 0.13


def load_report(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _solar_cols(df):
    return [c for c in df.columns if c.startswith("solar_az") and c.endswith("_gw_opt")]


def _wind_cols(df):
    return [c for c in df.columns if
            (c.startswith("wind-onshore") or c.startswith("wind-offshore"))
            and c.endswith("_gw_opt")]


def build_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    solar_cols = _solar_cols(df)
    wind_cols  = _wind_cols(df)
    rows = []
    for _, r in df.iterrows():
        row = {"label": str(r["scenario"])}
        row["dri_h2_mw_lhv"]    = r.get("dri_h2_mw_lhv", float("nan"))
        row["electrolyser_mw"]  = r.get("electrolyser_gw", 0) * 1e3
        row["battery_mw"]       = r.get("battery_gw_opt", 0) * 1e3
        # Buffer reported as hours of DRI H₂ demand by compile_report.
        row["h2_buffer_hours_dri"] = r.get("h2_buffer_hours_dri", 0)
        # If no orientation-resolved solar columns are present, fall back to the
        # plain solar_gw_opt column so single-orientation runs (e.g. DE baseline)
        # still show a solar bar.
        if solar_cols:
            for col in solar_cols:
                az = col.replace("solar_", "").replace("_gw_opt", "")
                row[f"solar_{az}_mw"] = r.get(col, 0) * 1e3
        elif "solar_gw_opt" in df.columns:
            row["solar_mw"] = r.get("solar_gw_opt", 0) * 1e3
        for col in wind_cols:
            row[f"{col.replace('_gw_opt','')}_mw"] = r.get(col, 0) * 1e3
        row["grid_avg_mw"] = r.get("grid_avg_mw", 0)
        rows.append(row)
    return pd.DataFrame(rows).set_index("label")


def _sample_colormap(cmap, frac: float) -> str:
    """Sample the FCA colormap at a fraction in [0,1]. cmap is the
    [[frac, 'rgb(...)']] structure defined in _helpers."""
    frac = max(0.0, min(1.0, frac))
    nearest = min(cmap, key=lambda fc: abs(fc[0] - frac))
    return nearest[1]


def _pretty(col: str) -> str:
    """Humanize a plot_df column name for the legend, e.g.
    'solar_mw' → 'Solar (MW nominal)',
    'solar_az180_mw' → 'Solar az180 (MW nominal)',
    'wind-onshore_mw' → 'Wind onshore (MW nominal)'."""
    base = col.removesuffix("_mw").replace("-", " ").replace("_", " ")
    return f"{base.capitalize()} (MW nominal)"


def _stacked_traces(
    plot_df: pd.DataFrame,
    cols: list[str],
    cmap_range: tuple[float, float],
    solo_color: str,
    slot: int,
    unit: str = "MW",
) -> list[go.Bar]:
    """Return a list of Bar traces that stack on top of each other at the same
    x-offset slot. Plotly's stackgroup is for scatter; for bars we set base
    manually."""
    offset = _SLOTS[slot] * _BAR_WIDTH - _BAR_WIDTH / 2
    base = pd.Series(0.0, index=plot_df.index)
    traces: list[go.Bar] = []
    if len(cols) == 1:
        colors = [solo_color]
    else:
        lo, hi = cmap_range
        fracs = [lo + (hi - lo) * i / (len(cols) - 1) for i in range(len(cols))]
        colors = [_sample_colormap(fca_colormap, f) for f in fracs]
    for col, color in zip(cols, colors):
        vals = plot_df[col].fillna(0)
        label = _pretty(col)
        traces.append(go.Bar(
            x=plot_df.index,
            y=vals,
            base=base.copy(),
            width=_BAR_WIDTH,
            offset=offset,
            marker_color=color,
            name=label,
            hovertemplate=f"{label}: %{{y:.1f}} {unit}<extra></extra>",
        ))
        base = base + vals.values
    return traces


def _single_bar(
    plot_df: pd.DataFrame,
    col: str,
    slot: int,
    color: str,
    name: str,
    unit: str = "MW",
) -> go.Bar:
    offset = _SLOTS[slot] * _BAR_WIDTH - _BAR_WIDTH / 2
    return go.Bar(
        x=plot_df.index,
        y=plot_df[col].fillna(0),
        width=_BAR_WIDTH,
        offset=offset,
        marker_color=color,
        name=name,
        hovertemplate=f"{name}: %{{y:.1f}} {unit}<extra></extra>",
    )


def plot(plot_df: pd.DataFrame, out: Path, project_label: str) -> None:
    solar_mw_cols = [c for c in plot_df.columns if c.startswith("solar_")]
    wind_mw_cols  = [c for c in plot_df.columns if c.startswith("wind_")]

    # Slot order (left → right): solar, wind, battery, electrolyser, buffer, DRI demand.
    fig = go.Figure()

    if solar_mw_cols:
        for tr in _stacked_traces(plot_df, solar_mw_cols, SOLAR_PALETTE_RANGE,
                                  SOLAR_SOLO_COLOR, 0):
            fig.add_trace(tr)

    if wind_mw_cols:
        for tr in _stacked_traces(plot_df, wind_mw_cols, WIND_PALETTE_RANGE,
                                  WIND_SOLO_COLOR, 1):
            fig.add_trace(tr)

    fig.add_trace(_single_bar(plot_df, "battery_mw", 2, magenta_red, "Battery (MW)"))
    fig.add_trace(_single_bar(plot_df, "electrolyser_mw", 3, green,
                              "Electrolyser (MW input)", "MW input"))
    fig.add_trace(_single_bar(plot_df, "h2_buffer_hours_dri", 4, dark_gray,
                              "H₂ buffer (hours of DRI demand)", "h of DRI demand"))
    fig.add_trace(_single_bar(plot_df, "dri_h2_mw_lhv", 5, blue_black,
                              "DRI H₂ demand (MW H₂ LHV)", "MW H₂ LHV"))

    n_sc = len(plot_df)
    fig.update_layout(
        template=fca_template,
        title=dict(
            text=f"{project_label} capacity breakdown",
            subtitle=dict(
                text="MW; H₂ buffer in hours of DRI demand",
                font=dict(family="Titillium Web", size=18, color=blue_black),
            ),
            font=dict(family="Titillium Web", size=26, color=blue_black),
            x=0, xref="paper", xanchor="left",
            y=1.0, yref="container", yanchor="top",
            pad=dict(t=14, l=0),
        ),
        barmode="overlay",
        bargap=0.1,
        width=max(720, 220 * n_sc + 280),
        height=600,
        yaxis_title=None,
        xaxis_title=None,
        legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top"),
        margin=dict(l=80, r=260, t=110, b=80),
    )
    fig.update_xaxes(type="category")
    fig.update_yaxes(rangemode="tozero")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(out, scale=2)
    fig.write_html(out.with_suffix(".html"), config=PLOTLY_CONFIG, include_plotlyjs="cdn")
    log.info("saved %s (+ .html)", out)


def main() -> None:
    df = load_report(_REPORT_PATH)
    project_label = ", ".join(dict.fromkeys(df["project"].astype(str)))
    plot(build_plot_data(df), _OUT, project_label)


if __name__ == "__main__":
    main()
