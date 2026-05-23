"""
diag_plot_complementarity.py

Plots a 2-week sample of the best complementarity triplet's time series for a
given country, as a visual sanity check on the complementarity.py output.

Usage:
    python scripts/res_cf/diag_plot_complementarity.py --country DE --year 2023
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
CF_DIR       = PROJECT_ROOT / "resources" / "res_cf"
PLOT_DIR     = PROJECT_ROOT / "results" / "plots"

# ── Import reusable functions from make_bestsite_cf ──────────────────────────
_spec = importlib.util.spec_from_file_location(
    "bestsite",
    Path(__file__).parent / "make_bestsite_cf.py"
)
_bestsite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bestsite)

build_cf_year           = _bestsite.build_cf_year
extract_cell_timeseries = _bestsite.extract_cell_timeseries


def plot_best_triplet(country: str, year: int) -> None:
    cc = country.lower()
    comp_path = CF_DIR / f"{cc}_complementarity_top10_{year}.csv"
    if not comp_path.exists():
        raise FileNotFoundError(
            f"Complementarity results not found: {comp_path}\n"
            "Run complementarity.py first."
        )

    df   = pd.read_csv(comp_path)
    best = df.iloc[0]

    print(f"Loading CF grids for {country} {year}...")
    cf_on  = build_cf_year(country, "wind_onshore")
    cf_off = build_cf_year(country, "wind_offshore")
    cf_sol = build_cf_year(country, "solar")

    ts_on  = extract_cell_timeseries(cf_on,  int(best["onshore_y_idx"]),  int(best["onshore_x_idx"]),  "wind_onshore")
    ts_off = extract_cell_timeseries(cf_off, int(best["offshore_y_idx"]), int(best["offshore_x_idx"]), "wind_offshore")
    ts_sol = extract_cell_timeseries(cf_sol, int(best["solar_y_idx"]),    int(best["solar_x_idx"]),    "solar")

    hours = slice(0, 24 * 14)  # 2 weeks
    t = range(24 * 14)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(t, ts_on.values[hours],  label="wind onshore",  alpha=0.8)
    ax.plot(t, ts_off.values[hours], label="wind offshore", alpha=0.8)
    ax.plot(t, ts_sol.values[hours], label="solar",         alpha=0.8)
    ax.set_title(
        f"{country} best triplet — 2-week sample  "
        f"(score={best['score']:.3f}, coincidence={best['coincidence']:.3f}, "
        f"dist={best['dist_km']:.0f} km)"
    )
    ax.set_xlabel("Hour")
    ax.set_ylabel("Capacity factor")
    ax.legend()
    plt.tight_layout()

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOT_DIR / f"{cc}_best_triplet_{year}.png"
    plt.savefig(out, dpi=150)
    print(f"Saved → {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--country", default="DE")
    parser.add_argument("--year",    default=2023, type=int)
    args = parser.parse_args()
    plot_best_triplet(args.country, args.year)
