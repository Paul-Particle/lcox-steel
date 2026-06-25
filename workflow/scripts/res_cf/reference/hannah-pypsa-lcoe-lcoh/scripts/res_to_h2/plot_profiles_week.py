import pandas as pd
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


def main():

    # ensure folder exists
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    country = "AUS"
    tech = "wind_onshore"
    year = 2023

    avg = load_df(country, tech, "avg", year)
    best = load_df(country, tech, "bestsite_p95", year)

    # pick 1 week
    start = "2023-07-01"
    end = "2023-07-08"

    avg = avg[(avg["time"] >= start) & (avg["time"] <= end)]
    best = best[(best["time"] >= start) & (best["time"] <= end)]

    plt.figure(figsize=(12, 5))
    plt.plot(avg["time"], avg["cf"], label="avg")
    plt.plot(best["time"], best["cf"], label="bestsite")

    plt.title(f"{country} {tech} – 1 week CF comparison")
    plt.ylabel("Capacity Factor")
    plt.xlabel("Time")
    plt.legend()
    plt.tight_layout()

    # save instead of show
    filename = f"{country}_{tech}_week_comparison.png"
    save_path = PLOTS_DIR / filename

    plt.savefig(save_path, dpi=150)
    plt.close()

    print(f"Saved plot to: {save_path}")


if __name__ == "__main__":
    main()