import pandas as pd
from pathlib import Path
from common._paths import RES_CF

OUTDIR = RES_CF
YEAR = 2023
COUNTRIES = ["fr", "es", "aus", "bra"]

def structural_qc(cc):
    path_w = OUTDIR / f"{cc}_wind_onshore_cf_{YEAR}.csv"
    path_s = OUTDIR / f"{cc}_solar_cf_{YEAR}.csv"

    for path in [path_w, path_s]:
        df = pd.read_csv(path)
        df["time"] = pd.to_datetime(df["time"])
        df = df.sort_values("time")

        print(f"\n=== {path.name} (structural) ===")
        print("rows:", len(df))
        print("start:", df["time"].iloc[0], "end:", df["time"].iloc[-1])
        print("duplicates:", df["time"].duplicated().sum())

        diffs = df["time"].diff().dropna()
        print("non-hourly steps:", (diffs != pd.Timedelta(hours=1)).sum())

if __name__ == "__main__":
    for cc in COUNTRIES:
        structural_qc(cc)
