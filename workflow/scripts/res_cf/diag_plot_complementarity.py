"""
diag_plot_complementarity.py

Plots a 2-week sample of the best complementarity triplet's time series for a
given country, as a visual sanity check on the determine_complementarity.py output.

Usage:
    python scripts/res_cf/diag_plot_complementarity.py --country DE --year 2023
"""

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common._logging import configure_logging
from common._paths import RES_CF, RESULTS

CF_DIR   = RES_CF
PLOT_DIR = RESULTS / "plots"

from determine_bestsite_p95 import build_cf_year, extract_cell_timeseries

configure_logging(None)
log = logging.getLogger(__name__)


def plot_best_triplet(country: str, year: int) -> None:
    cc = country.lower()
    comp_path = CF_DIR / f"{cc}_complementarity_top10_{year}.parquet"
    if not comp_path.exists():
        raise FileNotFoundError(
            f"Complementarity results not found: {comp_path}\n"
            "Run determine_complementarity.py first."
        )

    df   = pd.read_parquet(comp_path)
    best = df.iloc[0]

    log.info(f"loading CF grids for {country} {year}")
    cf_on  = build_cf_year(country, "wind_onshore")
    cf_off = build_cf_year(country, "wind_offshore")
    cf_sol = build_cf_year(country, "solar")

    ts_on  = extract_cell_timeseries(cf_on,  int(best["onshore_y_idx"]),  int(best["onshore_x_idx"]))
    ts_off = extract_cell_timeseries(cf_off, int(best["offshore_y_idx"]), int(best["offshore_x_idx"]))
    ts_sol = extract_cell_timeseries(cf_sol, int(best["solar_y_idx"]),    int(best["solar_x_idx"]))

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
    log.info(f"saved → {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--country", default="DE")
    parser.add_argument("--year",    default=2023, type=int)
    args = parser.parse_args()
    plot_best_triplet(args.country, args.year)
