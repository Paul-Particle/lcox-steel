"""
Combine yearly per-technology CF time series into a single country-year file.

Inputs (yearly, hourly; produced by Script 04):
- data/res_cf/de_wind_onshore_cf_2023.csv   (columns: time, cf)
- data/res_cf/de_wind_offshore_cf_2023.csv  (columns: time, cf)
- data/res_cf/de_solar_cf_2023.csv          (columns: time, cf)

Output:
- data/res_cf/de_cf_2023.csv
  columns: time, wind_onshore_cf, wind_offshore_cf, solar_cf

Notes:
- Uses an inner join on 'time' (should yield exactly 8760 rows for 2023).
- No aggregation; purely column-wise merge.
"""

from pathlib import Path
import pandas as pd

OUTDIR    = Path("data/res_cf")
COUNTRIES = ["de", "fr", "es", "aus", "bra"]  # standalone default
YEAR      = 2023

if "snakemake" in dir():
    COUNTRIES = [snakemake.wildcards.country.lower()]
    YEAR      = int(snakemake.wildcards.year)

def read_cf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["time"] = pd.to_datetime(df["time"])
    return df.sort_values("time")

def build_country(country: str):
    wind_path = OUTDIR / f"{country}_wind_onshore_cf_{YEAR}.csv"
    offshore_path = OUTDIR / f"{country}_wind_offshore_cf_{YEAR}.csv"
    solar_path = OUTDIR / f"{country}_solar_cf_{YEAR}.csv"

    wind = read_cf(wind_path).rename(columns={"cf": "wind_onshore_cf"})
    offshore = read_cf(offshore_path).rename(columns={"cf": "wind_offshore_cf"})
    solar = read_cf(solar_path).rename(columns={"cf": "solar_cf"})

    df = wind.merge(offshore, on="time", how="inner").merge(solar, on="time", how="inner")

    # integrity checks
    if len(df) != 8760:
        raise ValueError(f"{country}: Expected 8760 rows for {YEAR}, got {len(df)}")
    if df["time"].duplicated().any():
        raise ValueError(f"{country}: Duplicate timestamps found.")
    diffs = df["time"].diff().dropna()
    if not (diffs == pd.Timedelta(hours=1)).all():
        raise ValueError(f"{country}: Non-hourly continuity detected.")

    out = OUTDIR / f"{country}_cf_{YEAR}.csv"
    df.to_csv(out, index=False)
    print("Wrote:", out)


if __name__ == "__main__":
    for country in COUNTRIES:
        build_country(country)

