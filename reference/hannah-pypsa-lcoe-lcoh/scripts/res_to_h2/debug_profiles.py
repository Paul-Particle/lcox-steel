import numpy as np
import pandas as pd
from pathlib import Path


DATA_DIR = Path("data/res_cf")


def load_profile(country, tech, variant, year):
    cc = country.lower()

    if variant == "avg":
        filename = f"{cc}_cf_{year}.csv"
    elif variant == "bestsite_p95":
        filename = f"{cc}_cf_{year}_bestsite_p95.csv"
    else:
        raise ValueError("unknown variant")

    df = pd.read_csv(DATA_DIR / filename)

    col = f"{tech}_cf"
    return df[col].astype(float).to_numpy()


def print_stats(name, arr):
    print(f"\n{name}")
    print("-" * 30)
    print(f"mean : {arr.mean():.4f}")
    print(f"p10  : {np.percentile(arr, 10):.4f}")
    print(f"p5   : {np.percentile(arr, 5):.4f}")
    print(f"p1   : {np.percentile(arr, 1):.4f}")
    print(f"min  : {arr.min():.4f}")


def main():

    country = "AUS"
    tech = "wind_onshore"
    year = 2023

    avg = load_profile(country, tech, "avg", year)
    best = load_profile(country, tech, "bestsite_p95", year)

    print("\n=== PROFILE COMPARISON ===")

    print_stats("AVG", avg)
    print_stats("BESTSITE", best)

    print("\n=== INTERPRETATION GUIDE ===")
    print("If BESTSITE has lower p1/p5 than AVG → deeper low-wind periods → more storage needed")


if __name__ == "__main__":
    main()