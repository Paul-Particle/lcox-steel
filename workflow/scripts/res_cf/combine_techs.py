"""
Combine yearly per-technology CF time series into a single country-year file.

Inputs (yearly, hourly; produced by concat_quarters):
- resources/res_cf/annual/de_wind_onshore_2023.parquet   (DatetimeIndex 'time', column 'cf')
- resources/res_cf/annual/de_wind_offshore_2023.parquet  (same)
- resources/res_cf/annual/de_solar_2023.parquet          (same)

Output:
- resources/res_cf/de_cf_2023.parquet
  columns: time, wind_onshore_cf, wind_offshore_cf, solar_cf

Notes:
- Uses an inner join on 'time' (8760 rows for non-leap years, 8784 for leap years).
- No aggregation; purely column-wise merge.
"""

import calendar
from pathlib import Path
import pandas as pd

from common._paths import RES_CF

if "snakemake" not in globals():
    from common._stubs import snakemake

INDIR     = RES_CF / "annual"
OUTDIR    = RES_CF
COUNTRIES = ["de"]  # standalone default
YEAR      = 2023

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRIES = [snakemake.wildcards.country.lower()]
    YEAR      = int(snakemake.wildcards.year)


def expected_hours(year: int) -> int:
    return 8784 if calendar.isleap(year) else 8760


def read_cf(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path).reset_index()
    return df.sort_values("time")


def build_country(country: str) -> None:
    wind_path     = INDIR / f"{country}_wind_onshore_{YEAR}.parquet"
    offshore_path = INDIR / f"{country}_wind_offshore_{YEAR}.parquet"
    solar_path    = INDIR / f"{country}_solar_{YEAR}.parquet"

    wind     = read_cf(wind_path).rename(columns={"cf": "wind_onshore_cf"})
    offshore = read_cf(offshore_path).rename(columns={"cf": "wind_offshore_cf"})
    solar    = read_cf(solar_path).rename(columns={"cf": "solar_cf"})

    df = wind.merge(offshore, on="time", how="inner").merge(solar, on="time", how="inner")

    # integrity checks
    expected = expected_hours(YEAR)
    if len(df) != expected:
        raise ValueError(f"{country}: Expected {expected} rows for {YEAR}, got {len(df)}")
    if df["time"].duplicated().any():
        raise ValueError(f"{country}: Duplicate timestamps found.")
    diffs = df["time"].diff().dropna()
    if not (diffs == pd.Timedelta(hours=1)).all():
        raise ValueError(f"{country}: Non-hourly continuity detected.")

    out = OUTDIR / f"{country}_cf_{YEAR}.parquet"
    df.to_parquet(out, index=False)
    print("Wrote:", out)


if __name__ == "__main__":
    for country in COUNTRIES:
        build_country(country)

