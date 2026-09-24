"""Bar chart comparing route capacities within one scenario.

Snakemake rule: plot_capacity_bars (in viz.smk).

One panel per unit class, stacked vertically over a shared route axis.
Panels and bars with no data anywhere in the scenario are dropped:

  Power (MW)      solar stacked by orientation · wind stacked by tech ·
                  battery · electrolyser · grid connection ·
                  DRI H₂ demand (MW LHV, pure-H2 route only)
  Process (t/h)   DRI · EAF · MOE · electrowinning, in output units
  Storage (h)     H₂ buffer (hours of average DRI H₂ draw) ·
                  iron store (hours of steel output)
"""

import logging
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from common._report_schema import PROCESS_LINKS, field_stem, read_report
from _run_display import run_label
from scripts.viz.style import (
    apply_header,
    blue_black,
    blue_gray,
    dark_gray,
    fca_blue,
    fca_colormap,
    fca_template,
    gray,
    green,
    highlight_blue,
    light_blue,
    light_blue_gray,
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

# Stacked sub-bar palettes, sampled from the FCA colormap.
# Solar: high end of the colormap (golden/sand tones — close to sand_yellow).
# Wind:  low/mid end of the colormap (blue/teal tones — close to fca_blue).
SOLAR_PALETTE_RANGE = (0.78, 1.00)
WIND_PALETTE_RANGE  = (0.02, 0.40)

# Solo-trace colors when only one sub-bar is present.
SOLAR_SOLO_COLOR = sand_yellow
WIND_SOLO_COLOR  = fca_blue

_BAR_WIDTH = 0.12

PROCESS_BARS = [
    ("dri_h2_t_per_h", "H2-DRI shaft (t/h iron)",   blue_gray),
    ("dri_ng_t_per_h", "NG-DRI shaft (t/h iron)",   light_blue_gray),
    ("eaf_t_per_h",    "EAF (t/h steel)",           very_dark_gray),
    ("moe_t_per_h",    "MOE (t/h steel)",           highlight_blue),
    ("ew_t_per_h",     "Electrowinning (t/h iron)", turquois),
]

STORAGE_BARS = [
    ("h2_buffer_hours_dri",   "H₂ buffer (h of avg DRI H₂ draw)", light_blue),
    ("iron_store_hours_steel", "Iron store (h of steel output)",  gray),
    ("steel_store_hours_steel", "Steel inventory (h of steel output)", blue_gray),
]


def _solar_cols(df):
    # Orientation-resolved (solar_az*) or per-site (multisite: solar-c00 …) solar
    # capacities. Excludes the bare single-site `solar_gw_opt` (handled by the
    # fallback in build_plot_data) and the multisite `_total` aggregate, which
    # would otherwise stack on top of the per-site columns it already sums.
    return [c for c in df.columns
            if c.startswith("solar") and c.endswith("_gw_opt")
            and c != "solar_gw_opt" and not c.endswith("_total_gw_opt")]


def _wind_cols(df):
    # Per-tech (single-site) or per-site (multisite) wind capacities. Excludes the
    # multisite `{tech}_total_gw_opt` aggregate so it isn't double-counted against
    # the per-site columns it sums.
    return [c for c in df.columns
            if c.startswith(("wind_onshore", "wind_offshore"))
            and c.endswith("_gw_opt") and not c.endswith("_total_gw_opt")]





def build_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape a report DataFrame into per-route rows of plottable capacities.

    Converts GW columns to MW, expands orientation-resolved solar columns (falling
    back to a single `solar_mw` when none are present), and carries process
    capacities (t/h output) and stores (hours of demand) through in their native
    units. Indexed by run label.
    """
    solar_cols = _solar_cols(df)
    wind_cols  = _wind_cols(df)
    rows = []
    for _, r in df.iterrows():
        row = {"label": run_label(r)}
        row["dri_h2_mw_lhv"]    = r["dri_h2_mw_lhv"]
        row["electrolyser_mw"]  = r["electrolyser_gw"] * 1e3
        row["battery_mw"]       = r["battery_gw_opt"] * 1e3
        row["grid_mw"]          = r["grid_import_gw_opt"] * 1e3
        # If no orientation-resolved solar columns are present, fall back to the
        # plain solar_gw_opt column so single-orientation runs (e.g. DE baseline)
        # still show a solar bar.
        if solar_cols:
            for col in solar_cols:
                az = col.replace("solar_", "").replace("_gw_opt", "")
                row[f"solar_{az}_mw"] = r[col] * 1e3
        elif "solar_gw_opt" in df.columns:
            row["solar_mw"] = r["solar_gw_opt"] * 1e3
        for col in wind_cols:
            row[f"{col.replace('_gw_opt','')}_mw"] = r[col] * 1e3
        # Every report row carries `{link}_t_per_h_opt` for each of
        # PROCESS_LINKS (common/_report_schema.py), reading 0 on a route
        # without that step.
        for link in PROCESS_LINKS:
            stem = field_stem(link)
            row[f"{stem}_t_per_h"] = r[f"{stem}_t_per_h_opt"]
        row["h2_buffer_hours_dri"]    = r["h2_buffer_hours_dri"]
        row["iron_store_hours_steel"] = r["iron_store_hours_steel"]
        row["steel_store_hours_steel"] = r["steel_store_hours_steel"]
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
    'wind_onshore_mw' → 'Wind onshore (MW nominal)'."""
    base = col.removesuffix("_mw").replace("_", " ")
    return f"{base.capitalize()} (MW nominal)"


def _stacked_traces(
    plot_df: pd.DataFrame,
    cols: list[str],
    cmap_range: tuple[float, float],
    solo_color: str,
    offset: float,
    unit: str = "MW",
) -> list[go.Bar]:
    """Return a list of Bar traces that stack on top of each other at the same
    x-offset. Plotly's stackgroup is for scatter; for bars we set base
    manually."""
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
    offset: float,
    color: str,
    name: str,
    unit: str = "MW",
) -> go.Bar:
    """Return a single Bar trace at the given x-offset."""
    return go.Bar(
        x=plot_df.index,
        y=plot_df[col].fillna(0),
        width=_BAR_WIDTH,
        offset=offset,
        marker_color=color,
        name=name,
        hovertemplate=f"{name}: %{{y:.1f}} {unit}<extra></extra>",
    )


def _build_panels(plot_df: pd.DataFrame) -> list[tuple[str, list]]:
    """Assemble the non-empty panels as (title, [slot, ...]) pairs.

    A slot is one x-offset position: ("stack", cols, cmap_range, solo_color,
    unit) or ("single", col, color, name, unit). Slots empty across the whole
    scenario are dropped so they don't leave gaps and phantom legend entries;
    panels with no surviving slot are dropped entirely.
    """
    def _has_data(cols) -> bool:
        cols = [c for c in cols if c in plot_df.columns]
        return bool(cols) and bool(plot_df[cols].fillna(0).to_numpy().any())

    solar_mw_cols = [c for c in plot_df.columns if c.startswith("solar")]
    wind_mw_cols  = [c for c in plot_df.columns
                     if c.startswith(("wind_onshore", "wind_offshore"))]

    power_slots: list = []
    if _has_data(solar_mw_cols):
        power_slots.append(("stack", solar_mw_cols, SOLAR_PALETTE_RANGE,
                            SOLAR_SOLO_COLOR, "MW"))
    if _has_data(wind_mw_cols):
        power_slots.append(("stack", wind_mw_cols, WIND_PALETTE_RANGE,
                            WIND_SOLO_COLOR, "MW"))
    if _has_data(["battery_mw"]):
        power_slots.append(("single", "battery_mw", magenta_red,
                            "Battery (MW)", "MW"))
    if _has_data(["electrolyser_mw"]):
        power_slots.append(("single", "electrolyser_mw", green,
                            "Electrolyser (MW input)", "MW input"))
    if _has_data(["grid_mw"]):
        power_slots.append(("single", "grid_mw", dark_gray,
                            "Grid connection (MW)", "MW"))
    if _has_data(["dri_h2_mw_lhv"]):
        power_slots.append(("single", "dri_h2_mw_lhv", blue_black,
                            "DRI H₂ demand (MW H₂ LHV)", "MW H₂ LHV"))

    process_slots = [
        ("single", col, color, name, "t/h")
        for col, name, color in PROCESS_BARS if _has_data([col])
    ]
    storage_slots = [
        ("single", col, color, name, "h")
        for col, name, color in STORAGE_BARS if _has_data([col])
    ]

    panels = []
    if power_slots:
        panels.append(("Power capacity (MW)", power_slots))
    if process_slots:
        panels.append(("Process capacity (t/h output)", process_slots))
    if storage_slots:
        panels.append(("Storage (hours of average demand)", storage_slots))
    return panels


def plot(plot_df: pd.DataFrame, out: Path, scenario_label: str) -> None:
    """Assemble the panelled capacity bar chart and write it to PNG + HTML."""
    panels = _build_panels(plot_df)

    fig = make_subplots(
        rows=len(panels), cols=1, shared_xaxes=True,
        vertical_spacing=0.34 / len(panels),
        subplot_titles=[title for title, _ in panels],
    )
    for row_i, (_, slots) in enumerate(panels, start=1):
        n_slots = len(slots)
        for slot_i, spec in enumerate(slots):
            offset = (slot_i - n_slots / 2) * _BAR_WIDTH
            if spec[0] == "stack":
                _, cols, cmap_range, solo_color, unit = spec
                traces = _stacked_traces(plot_df, cols, cmap_range, solo_color,
                                         offset, unit)
            else:
                _, col, color, name, unit = spec
                traces = [_single_bar(plot_df, col, offset, color, name, unit)]
            for tr in traces:
                fig.add_trace(tr, row=row_i, col=1)

    n_sc = len(plot_df)
    fig.update_layout(
        template=fca_template,
        barmode="overlay",
        bargap=0.1,
        # No rotated y-axis titles — units live in the panel titles.
        yaxis_title=None,
        xaxis_title=None,
        legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top"),
    )
    fig.update_xaxes(type="category")
    fig.update_yaxes(rangemode="tozero")
    # Panel titles come out at plotly's default size; align them with the
    # template's body text and left edge.
    fig.update_annotations(font_size=15, x=0.0, xanchor="left")

    apply_header(
        fig,
        title=f"{scenario_label} capacity breakdown",
        subtitle="Optimised builds by route; one panel per unit class",
        fig_width=max(720, 220 * n_sc + 300),
        fig_height=180 + 270 * len(panels),
        margin_l=80, margin_r=280, margin_t=110, margin_b=80,
    )
    saved = save_figure(fig, out.parent, out.stem)
    log.info(f"saved {' + '.join(saved)}")


def main() -> None:
    """Load the scenario report and render its capacity bar chart."""
    df = read_report(_REPORT_PATH)
    scenario_label = ", ".join(dict.fromkeys(df["scenario"].astype(str)))
    plot(build_plot_data(df), _OUT, scenario_label)


if __name__ == "__main__":
    main()
