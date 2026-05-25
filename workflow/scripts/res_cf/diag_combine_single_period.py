"""
Combines per-technology CF files from a single-period run (e.g. jan2w)
into the merged format expected by the plot script.

Use this instead of running scripts 04 + 05 when you only have one period
and just want to get to a plot quickly.

Output: resources/res_cf/{COUNTRY}_cf_{YEAR}.parquet
        columns: time, wind_onshore_cf, [wind_offshore_cf,] solar_cf
"""

from pathlib import Path
import pandas as pd
from common._paths import RES_CF

COUNTRY = "de"
TAG     = "oct2w"
YEAR    = 2025

DATA_DIR = RES_CF


def read_cf(path: Path, col_name: str) -> pd.DataFrame:
    df = pd.read_parquet(path).reset_index()
    return df.sort_values("time").rename(columns={"cf": col_name})[["time", col_name]]


def main():
    wind_path    = DATA_DIR / f"{COUNTRY}_wind_onshore_cf_{YEAR}_{TAG}.parquet"
    solar_path   = DATA_DIR / f"{COUNTRY}_solar_cf_{YEAR}_{TAG}.parquet"
    offshore_path = DATA_DIR / f"{COUNTRY}_wind_offshore_cf_{YEAR}_{TAG}.parquet"

    df = read_cf(wind_path, "wind_onshore_cf").merge(
         read_cf(solar_path, "solar_cf"), on="time", how="inner")

    if offshore_path.exists():
        df = df.merge(read_cf(offshore_path, "wind_offshore_cf"), on="time", how="inner")
        df = df[["time", "wind_onshore_cf", "wind_offshore_cf", "solar_cf"]]
    else:
        df = df[["time", "wind_onshore_cf", "solar_cf"]]

    out = DATA_DIR / f"{COUNTRY}_cf_{YEAR}.parquet"
    df.to_parquet(out, index=False)
    print(f"Wrote {len(df)} rows → {out}")


if __name__ == "__main__":
    main()
