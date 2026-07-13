"""Bar chart comparing scenario capacities across projects.

Snakemake rule: plot_capacity_bars (in viz.smk).

Shows per scenario (slots with no data in the whole project are dropped):
  Solar stacked by orientation (MW) · Wind stacked by tech (MW) · Battery (MW)
  Electrolyser (MW) · H2 buffer (hours of DRI demand) · DRI H2 demand (MW LHV)
  Grid connection (MW)
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
    fca_colormap,
    fca_template,
    gray,
    green,
    magenta_red,
    sand_yellow,
    save_figure,
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

# Slot positions (in "category widths") for the seven grouped bars.
_SLOTS = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
_BAR_WIDTH = 0.12


def load_report(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _solar_cols(df):
    return [c for c in df.columns if c.startswith("solar_az") and c.endswith("_gw_opt")]


def _wind_cols(df):
    return [c for c in df.columns if
            (c.startswith("wind-onshore") or c.startswith("wind-offshore"))
            and c.endswith("_gw_opt")]


def build_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape a report DataFrame into per-scenario rows of plottable MW capacities.

    Converts GW columns to MW, expands orientation-resolved solar columns (falling
    back to a single `solar_mw` when none are present), and keeps the H2 buffer in
    hours of DRI demand. Indexed by scenario label.
    """
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
        row["grid_mw"] = r.get("grid_import_gw_opt", 0) * 1e3
        rows.append(row)
    return pd.DataFrame(rows).set_index("label")


def _sample_colormap(cmap, frac: float) -> str:
    """Sample the FCA colormap at a fraction in [0,1]. cmap is the
    [[frac, 'rgb(...)']] structure defined in style.py."""
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
    """Return a single Bar trace positioned at the given grouped-bar slot."""
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
    """Assemble the grouped capacity bar chart and write it to PNG + HTML.

    Lays out six slots left→right: solar (stacked by orientation), wind (stacked
    by tech), battery, electrolyser, H2 buffer, and DRI H2 demand.
    """
    def _has_data(cols) -> bool:
        """True if any scenario has a nonzero value in any of the columns.

        Slots that are empty across the whole project (e.g. electrolyser on
        MOE/EW routes, DRI H₂ demand on steel routes) are dropped entirely so
        they don't leave gaps and phantom legend entries.
        """
        cols = [c for c in cols if c in plot_df.columns]
        return bool(cols) and bool(plot_df[cols].fillna(0).to_numpy().any())

    solar_mw_cols = [c for c in plot_df.columns if c.startswith("solar")]
    wind_mw_cols  = [c for c in plot_df.columns
                     if c.startswith(("wind-onshore", "wind-offshore"))]

    # Slot order (left → right): solar, wind, battery, electrolyser, buffer,
    # DRI demand, grid connection.
    fig = go.Figure()

    if _has_data(solar_mw_cols):
        for tr in _stacked_traces(plot_df, solar_mw_cols, SOLAR_PALETTE_RANGE,
                                  SOLAR_SOLO_COLOR, 0):
            fig.add_trace(tr)

    if _has_data(wind_mw_cols):
        for tr in _stacked_traces(plot_df, wind_mw_cols, WIND_PALETTE_RANGE,
                                  WIND_SOLO_COLOR, 1):
            fig.add_trace(tr)

    if _has_data(["battery_mw"]):
        fig.add_trace(_single_bar(plot_df, "battery_mw", 2, magenta_red, "Battery (MW)"))
    if _has_data(["electrolyser_mw"]):
        fig.add_trace(_single_bar(plot_df, "electrolyser_mw", 3, green,
                                  "Electrolyser (MW input)", "MW input"))
    if _has_data(["h2_buffer_hours_dri"]):
        fig.add_trace(_single_bar(plot_df, "h2_buffer_hours_dri", 4, dark_gray,
                                  "H₂ buffer (hours of DRI demand)", "h of DRI demand"))
    if _has_data(["dri_h2_mw_lhv"]):
        fig.add_trace(_single_bar(plot_df, "dri_h2_mw_lhv", 5, blue_black,
                                  "DRI H₂ demand (MW H₂ LHV)", "MW H₂ LHV"))
    if _has_data(["grid_mw"]):
        fig.add_trace(_single_bar(plot_df, "grid_mw", 6, gray, "Grid connection (MW)"))

    n_sc = len(plot_df)
    fig.update_layout(
        template=fca_template,
        barmode="overlay",
        bargap=0.1,
        # No rotated y-axis title — the unit lives in the subtitle.
        yaxis_title=None,
        xaxis_title=None,
        legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top"),
    )
    fig.update_xaxes(type="category")
    fig.update_yaxes(rangemode="tozero")

    apply_header(
        fig,
        title=f"{project_label} capacity breakdown",
        subtitle=(
            "MW; H₂ buffer in hours of DRI demand"
            if _has_data(["h2_buffer_hours_dri"])
            else "MW"
        ),
        fig_width=max(720, 220 * n_sc + 280), fig_height=600,
        margin_l=80, margin_r=260, margin_t=110, margin_b=80,
    )
    saved = save_figure(fig, out.parent, out.stem)
    log.info(f"saved {' + '.join(saved)}")


def main() -> None:
    """Load the project report and render its capacity bar chart."""
    df = load_report(_REPORT_PATH)
    project_label = ", ".join(dict.fromkeys(df["project"].astype(str)))
    plot(build_plot_data(df), _OUT, project_label)


if __name__ == "__main__":
    main()
