import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


DATA_DIR = Path("data/res_cf")
PLOTS_DIR = Path("results/plots/profiles_plots")


def load_df(country, tech, variant, year):
    cc = country.lower()

    if variant == "avg":
        filename = f"{cc}_cf_{year}.csv"
    else:
        filename = f"{cc}_cf_{year}_bestsite_p95.csv"

    df = pd.read_csv(DATA_DIR / filename, parse_dates=["time"])
    col = f"{tech}_cf"

    return df[["time", col]].rename(columns={col: "cf"})


def plot_full_year_with_rolling(avg, best, country, tech):
    plt.figure(figsize=(14, 5))

    # raw
    plt.plot(avg["time"], avg["cf"], alpha=0.3, label="avg (raw)")
    plt.plot(best["time"], best["cf"], alpha=0.3, label="bestsite (raw)")

    # rolling mean (24h)
    avg_roll = avg["cf"].rolling(24).mean()
    best_roll = best["cf"].rolling(24).mean()

    plt.plot(avg["time"], avg_roll, linewidth=2, label="avg (24h mean)")
    plt.plot(best["time"], best_roll, linewidth=2, label="bestsite (24h mean)")

    plt.title(f"{country} {tech} – full year CF (with smoothing)")
    plt.ylabel("Capacity Factor")
    plt.legend()
    plt.tight_layout()

    save_path = PLOTS_DIR / f"{country}_{tech}_full_year.png"
    plt.savefig(save_path, dpi=150)
    plt.close()

    print(f"Saved: {save_path}")


def plot_duration_curve(avg, best, country, tech):
    plt.figure(figsize=(8, 5))

    avg_sorted = np.sort(avg["cf"])[::-1]
    best_sorted = np.sort(best["cf"])[::-1]

    x = np.arange(len(avg_sorted)) / len(avg_sorted)

    plt.plot(x, avg_sorted, label="avg")
    plt.plot(x, best_sorted, label="bestsite")

    plt.title(f"{country} {tech} – duration curve")
    plt.xlabel("Fraction of time")
    plt.ylabel("Capacity Factor")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    save_path = PLOTS_DIR / f"{country}_{tech}_duration_curve.png"
    plt.savefig(save_path, dpi=150)
    plt.close()

    print(f"Saved: {save_path}")


def main():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    country = "AUS"
    tech = "wind_onshore"
    year = 2023

    avg = load_df(country, tech, "avg", year)
    best = load_df(country, tech, "bestsite_p95", year)

    plot_full_year_with_rolling(avg, best, country, tech)
    plot_duration_curve(avg, best, country, tech)


if __name__ == "__main__":
    main()
    