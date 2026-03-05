"""
Concatenate quarterly (Q1–Q4) RES capacity factor time series into a full-year series.

Inputs (CSV, hourly; produced by Script 03):
- data/res_cf/<cc>_wind_onshore_cf_<year>_q1.csv
- data/res_cf/<cc>_wind_onshore_cf_<year>_q2.csv
- data/res_cf/<cc>_wind_onshore_cf_<year>_q3.csv
- data/res_cf/<cc>_wind_onshore_cf_<year>_q4.csv
(and same for solar)

Expected input format:
- Columns: time, cf

Outputs (CSV, hourly):
- data/res_cf/<cc>_wind_onshore_cf_<year>.csv
- data/res_cf/<cc>_solar_cf_<year>.csv

Checks:
- Hourly continuity (1h steps)
- No duplicate timestamps
- Total hours = 8760 (non-leap years)
"""

from pathlib import Path
import pandas as pd

OUTDIR = Path("data/res_cf")
YEAR = 2023

def read_series(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    df["time"] = pd.to_datetime(df["time"])
    df = df.sort_values("time")
    return pd.Series(df["cf"].values, index=df["time"], name="cf")

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
    wind_q = [OUTDIR / f"{cc}_wind_onshore_cf_{YEAR}_q{i}.csv" for i in [1,2,3,4]]
    solar_q = [OUTDIR / f"{cc}_solar_cf_{YEAR}_q{i}.csv" for i in [1,2,3,4]]

    wind = stitch_quarters(wind_q)
    solar = stitch_quarters(solar_q)

    wind_out = OUTDIR / f"{cc}_wind_onshore_cf_{YEAR}.csv"
    solar_out = OUTDIR / f"{cc}_solar_cf_{YEAR}.csv"

    wind.reset_index().rename(columns={"index":"time"}).to_csv(wind_out, index=False)
    solar.reset_index().rename(columns={"index":"time"}).to_csv(solar_out, index=False)

    print("Wrote:")
    print(" -", wind_out)
    print(" -", solar_out)

if __name__ == "__main__":
    # add/remove countries here
    for cc in ["fr", "es", "aus", "bra"]:
        run_for_country(cc)
