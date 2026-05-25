"""
Concatenate quarterly (Q1–Q4) RES capacity factor time series into a full-year series.

Inputs (parquet, hourly; produced by build_cf_timeseries):
- resources/res_cf/quarterly/<cc>_wind_onshore_<year>_q{1..4}.parquet
- resources/res_cf/quarterly/<cc>_wind_offshore_<year>_q{1..4}.parquet
- resources/res_cf/quarterly/<cc>_solar_<year>_q{1..4}.parquet

Expected input format:
- DatetimeIndex 'time', column 'cf'

Outputs (parquet, hourly):
- resources/res_cf/annual/<cc>_wind_onshore_<year>.parquet
- resources/res_cf/annual/<cc>_wind_offshore_<year>.parquet
- resources/res_cf/annual/<cc>_solar_<year>.parquet

Checks:
- Hourly continuity (1h steps)
- No duplicate timestamps
- Total hours = 8760 (non-leap years)
"""

from pathlib import Path
import pandas as pd

from common._paths import RES_CF

if "snakemake" not in globals():
    from common._stubs import snakemake

INDIR  = RES_CF / "quarterly"
OUTDIR = RES_CF / "annual"
YEAR = 2023
COUNTRIES = ["de", "fr", "es", "aus", "bra"]  # standalone default

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRIES = [snakemake.wildcards.country.lower()]
    YEAR      = int(snakemake.wildcards.year)

def read_series(path: Path) -> pd.Series:
    df = pd.read_parquet(path)
    return df["cf"].sort_index()

def stitch_quarters(paths, expected_hours=8760) -> pd.Series:
    s = pd.concat([read_series(p) for p in paths]).sort_index()

    if s.index.duplicated().any():
        s = s[~s.index.duplicated(keep="last")]

    diffs = s.index.to_series().diff().dropna()
    bad = diffs[diffs != pd.Timedelta(hours=1)]
    if len(bad) > 0:
        raise ValueError(f"Non-hourly steps detected:\n{bad.head(10)}")

    if len(s) != expected_hours:
        raise ValueError(f"Expected {expected_hours} hours, got {len(s)}")

    return s

def run_for_country(cc: str):
    wind_q     = [INDIR / f"{cc}_wind_onshore_{YEAR}_q{i}.parquet"  for i in [1,2,3,4]]
    offshore_q = [INDIR / f"{cc}_wind_offshore_{YEAR}_q{i}.parquet" for i in [1,2,3,4]]
    solar_q    = [INDIR / f"{cc}_solar_{YEAR}_q{i}.parquet"         for i in [1,2,3,4]]

    wind = stitch_quarters(wind_q)
    offshore = stitch_quarters(offshore_q)
    solar = stitch_quarters(solar_q)

    wind_out     = OUTDIR / f"{cc}_wind_onshore_{YEAR}.parquet"
    offshore_out = OUTDIR / f"{cc}_wind_offshore_{YEAR}.parquet"
    solar_out    = OUTDIR / f"{cc}_solar_{YEAR}.parquet"

    wind.to_frame().to_parquet(wind_out, index=True)
    offshore.to_frame().to_parquet(offshore_out, index=True)
    solar.to_frame().to_parquet(solar_out, index=True)

    print("Wrote:")
    print(" -", wind_out)
    print(" -", offshore_out)
    print(" -", solar_out)

def main():
    for cc in COUNTRIES:
        run_for_country(cc)

if __name__ == "__main__":
    main()
