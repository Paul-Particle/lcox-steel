"""
Combines per-technology CF files from a single-period run (e.g. jan2w)
into the merged format expected by the plot script.

Use this instead of running scripts 04 + 05 when you only have one period
and just want to get to a plot quickly.

Output: data/res_cf/{COUNTRY}_cf_{YEAR}.csv
        columns: time, wind_onshore_cf, [wind_offshore_cf,] solar_cf
"""

from pathlib import Path
import pandas as pd

COUNTRY = "de"
TAG     = "oct2w"
YEAR    = 2025

DATA_DIR = Path("data/res_cf")


def read_cf(path: Path, col_name: str) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["time"])
    df = df.sort_values("time").rename(columns={"cf": col_name})
    return df[["time", col_name]]


def main():
    wind_path    = DATA_DIR / f"{COUNTRY}_wind_onshore_cf_{YEAR}_{TAG}.csv"
    solar_path   = DATA_DIR / f"{COUNTRY}_solar_cf_{YEAR}_{TAG}.csv"
    offshore_path = DATA_DIR / f"{COUNTRY}_wind_offshore_cf_{YEAR}_{TAG}.csv"

    df = read_cf(wind_path, "wind_onshore_cf").merge(
         read_cf(solar_path, "solar_cf"), on="time", how="inner")

    if offshore_path.exists():
        df = df.merge(read_cf(offshore_path, "wind_offshore_cf"), on="time", how="inner")
        df = df[["time", "wind_onshore_cf", "wind_offshore_cf", "solar_cf"]]
    else:
        df = df[["time", "wind_onshore_cf", "solar_cf"]]

    out = DATA_DIR / f"{COUNTRY}_cf_{YEAR}.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {len(df)} rows → {out}")


if __name__ == "__main__":
    main()
