"""
Bar chart comparing scenario capacities across projects.

Snakemake rule: plot_capacity_bars  (in viz.smk)
Standalone:     python workflow/scripts/viz/plot_capacity_bars.py results/report_*.csv

Shows per scenario:
  DRI H2 demand (MW LHV) · Electrolyser (MW) · H2 buffer (MW·days LHV)
  Solar stacked by orientation (MW) · Wind stacked by tech (MW)
  Grid annual-average import (MW)
"""

import logging
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if "snakemake" not in globals():
    sys.path.insert(0, str(Path(__file__).parents[2]))
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

_REPORT_PATHS: list[Path] = []
_OUT = Path("results/plots/capacity_bars.png")

if "snakemake" in globals() and hasattr(snakemake, "input"):
    _REPORT_PATHS = [Path(p) for p in snakemake.input.reports]
    _OUT = Path(snakemake.output[0])

SOLAR_CMAP = plt.cm.YlOrRd
WIND_CMAP  = plt.cm.Blues


def load_reports(paths: list[Path]) -> pd.DataFrame:
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)


def _solar_cols(df):
    return [c for c in df.columns if c.startswith("solar_az") and c.endswith("_gw_opt")]

def _wind_cols(df):
    return [c for c in df.columns if
            (c.startswith("wind_onshore") or c.startswith("wind_offshore"))
            and c.endswith("_gw_opt")]


def build_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    solar_cols = _solar_cols(df)
    wind_cols  = _wind_cols(df)
    rows = []
    for _, r in df.iterrows():
        row = {"label": f"{r['project']}\n({r['scenario']})"}
        row["dri_h2_mw_lhv"]   = r.get("dri_h2_mw_lhv", float("nan"))
        row["electrolyser_mw"]  = r.get("electrolyser_gw", 0) * 1e3
        row["h2_buffer_mwdays"] = r.get("h2_buffer_mwh_lhv_opt", 0) / 24.0
        for col in solar_cols:
            az = col.replace("solar_", "").replace("_gw_opt", "")
            row[f"solar_{az}_mw"] = r.get(col, 0) * 1e3
        for col in wind_cols:
            row[f"{col.replace('_gw_opt','')}_mw"] = r.get(col, 0) * 1e3
        row["grid_avg_mw"] = r.get("grid_avg_mw", 0)
        rows.append(row)
    return pd.DataFrame(rows).set_index("label")


def plot(plot_df: pd.DataFrame, out: Path) -> None:
    solar_mw_cols = [c for c in plot_df.columns if c.startswith("solar_az")]
    wind_mw_cols  = [c for c in plot_df.columns if c.startswith("wind_")]

    n_sc    = len(plot_df)
    x       = np.arange(n_sc)
    width   = 0.13
    slots   = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
    offsets = np.array(slots) * width

    fig, ax = plt.subplots(figsize=(max(10, n_sc * 3), 6))
    handles = []

    def bar(slot, vals, color, label):
        ax.bar(x + offsets[slot], vals, width, color=color)
        handles.append(mpatches.Patch(color=color, label=label))

    bar(0, plot_df["dri_h2_mw_lhv"].values,   "#444",    "DRI H₂ demand (MW LHV)")
    bar(1, plot_df["electrolyser_mw"].values,  "#e07b39", "Electrolyser (MW)")
    bar(2, plot_df["h2_buffer_mwdays"].values, "#5b9bd5", "H₂ buffer (MW·days LHV)")

    if solar_mw_cols:
        colors  = SOLAR_CMAP(np.linspace(0.35, 0.9, len(solar_mw_cols)))
        bottoms = np.zeros(n_sc)
        for col, color in zip(solar_mw_cols, colors):
            vals = plot_df[col].fillna(0).values
            ax.bar(x + offsets[3], vals, width, bottom=bottoms, color=color)
            bottoms += vals
        handles.append(mpatches.Patch(color=SOLAR_CMAP(0.65),
                                      label="Solar (MW, stacked by orientation)"))

    if wind_mw_cols:
        colors  = WIND_CMAP(np.linspace(0.4, 0.85, len(wind_mw_cols)))
        bottoms = np.zeros(n_sc)
        for col, color in zip(wind_mw_cols, colors):
            vals = plot_df[col].fillna(0).values
            ax.bar(x + offsets[4], vals, width, bottom=bottoms, color=color)
            bottoms += vals
        handles.append(mpatches.Patch(color=WIND_CMAP(0.6),
                                      label="Wind (MW, stacked by tech)"))

    grid = plot_df["grid_avg_mw"].fillna(0).values
    if grid.any():
        bar(5, grid, "#7cb87c", "Grid avg import (MW)")

    ax.set_xticks(x)
    ax.set_xticklabels(plot_df.index.tolist(), fontsize=9)
    ax.set_ylabel("Capacity  (MW  or  MW·days for H₂ buffer)", fontsize=10)
    ax.set_title("Scenario capacity breakdown", fontsize=12)
    ax.legend(handles=handles, fontsize=8, loc="upper right")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.3)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("saved %s", out)


def main() -> None:
    paths = _REPORT_PATHS or [Path(p) for p in sys.argv[1:]]
    if not paths:
        raise SystemExit("usage: plot_capacity_bars.py report1.csv [report2.csv ...]")
    df = load_reports(paths)
    plot(build_plot_data(df), _OUT)


if __name__ == "__main__":
    main()
