"""
07_make_bestsite_cf_timeseries.py

Purpose
-------
Creates "best-site" hourly CF time series by applying uplift factors (P90/P95)
derived in 06_resource_spread.py to the national hourly CF series.

Method
------
- Reads baseline national hourly CFs: data/res_cf/<cc>_cf_2023.csv
- Reads spatial uplift stats: data/res_cf/resource_spread_2023.csv
- For each country + tech:
    factor_p90 = spatial_p90_mean / national_mean
    factor_p95 = spatial_p95_mean / national_mean
- Applies factors to hourly series and clips to [0, 1]

Notes
-----
- Purely climate-resource based uplift (no land-use/grid/permitting constraints)
- Keeps baseline national series untouched; writes separate scenario files

Outputs
-------
data/res_cf/<cc>_cf_2023_bestsite_p90.csv
data/res_cf/<cc>_cf_2023_bestsite_p95.csv
(columns: time, wind_onshore_cf, solar_cf)
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

YEAR = 2023
OUTDIR = Path("data/res_cf")

COUNTRIES = ["de", "fr", "es", "aus", "bra"]  # lowercase to match filenames

SPREAD_PATH = OUTDIR / f"resource_spread_{YEAR}.csv"


def load_factors(spread_path: Path) -> dict[tuple[str, str], dict[str, float]]:
    """
    Returns mapping:
      (country_iso2_upper, tech) -> {"p90": factor, "p95": factor}
    where tech in {"wind_onshore", "solar"}.
    """
    df = pd.read_csv(spread_path)
    required = {
        "country", "tech", "national_mean",
        "spatial_p90_mean", "spatial_p95_mean"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{spread_path} missing columns: {sorted(missing)}")

    factors = {}
    for _, r in df.iterrows():
        country = str(r["country"]).upper()
        tech = str(r["tech"]).strip()

        nat = float(r["national_mean"])
        p90 = float(r["spatial_p90_mean"])
        p95 = float(r["spatial_p95_mean"])

        if nat <= 0 or not np.isfinite(nat):
            f90 = np.nan
            f95 = np.nan
        else:
            f90 = p90 / nat
            f95 = p95 / nat

        factors[(country, tech)] = {"p90": f90, "p95": f95}

    return factors


def apply_factor(series: pd.Series, factor: float) -> pd.Series:
    if not np.isfinite(factor):
        # if factor missing, return unchanged
        return series
    return np.clip(series.values * factor, 0.0, 1.0)


def main() -> None:
    if not SPREAD_PATH.exists():
        raise FileNotFoundError(
            f"Missing {SPREAD_PATH}. Run 06_resource_spread.py first."
        )

    factors = load_factors(SPREAD_PATH)

    for cc in COUNTRIES:
        base_path = OUTDIR / f"{cc}_cf_{YEAR}.csv"
        if not base_path.exists():
            raise FileNotFoundError(f"Missing baseline CF file: {base_path}")

        df = pd.read_csv(base_path)

        # basic schema checks
        for col in ["time", "wind_onshore_cf", "solar_cf"]:
            if col not in df.columns:
                raise ValueError(f"{base_path} missing column: {col}")

        country_upper = cc.upper()

        # Wind factors
        f_w_p90 = factors.get((country_upper, "wind_onshore"), {}).get("p90", np.nan)
        f_w_p95 = factors.get((country_upper, "wind_onshore"), {}).get("p95", np.nan)

        # Solar factors
        f_s_p90 = factors.get((country_upper, "solar"), {}).get("p90", np.nan)
        f_s_p95 = factors.get((country_upper, "solar"), {}).get("p95", np.nan)

        # Build bestsite P90
        df_p90 = df.copy()
        df_p90["wind_onshore_cf"] = apply_factor(df_p90["wind_onshore_cf"], f_w_p90)
        df_p90["solar_cf"] = apply_factor(df_p90["solar_cf"], f_s_p90)

        out_p90 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p90.csv"
        df_p90.to_csv(out_p90, index=False)

        # Build bestsite P95
        df_p95 = df.copy()
        df_p95["wind_onshore_cf"] = apply_factor(df_p95["wind_onshore_cf"], f_w_p95)
        df_p95["solar_cf"] = apply_factor(df_p95["solar_cf"], f_s_p95)

        out_p95 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95.csv"
        df_p95.to_csv(out_p95, index=False)

        print(
            f"{country_upper}: wrote "
            f"{out_p90.name} (wind×{f_w_p90:.3f}, solar×{f_s_p90:.3f}) and "
            f"{out_p95.name} (wind×{f_w_p95:.3f}, solar×{f_s_p95:.3f})"
        )


if __name__ == "__main__":
    main()
