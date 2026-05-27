"""Build one bidding zone's full ENTSO-E frame (grid variant=full).

Wide per-zone frame: price, load, RES forecasts, generation mix, crossborder
flows and derived residual-load columns, at native resolution. Kept for ad-hoc
analysis — not consumed by the optimisation. The model price is produced
separately by process_entsoe.py (variant=dayahead).
"""

from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake


def iso(yyyymmdd: str) -> str:
    """YYYYMMDD -> YYYY-MM-DD."""
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def to_utc_naive(obj):
    """Convert a tz-aware (Europe/Brussels) index to tz-naive UTC and sort."""
    if obj.index.tz is not None:
        obj.index = obj.index.tz_convert("UTC").tz_localize(None)
    return obj.sort_index()


def read_inputs(input_paths: list[str]) -> dict[str, pd.DataFrame]:
    """Map data_type -> its single-zone DataFrame.

    Filenames are ``{data_type}_{start_date}_{end_date}.parquet``; strip the two
    trailing date segments to recover the data_type key.
    """
    return {
        Path(p).stem.rsplit("_", 2)[0]: pd.read_parquet(p) for p in input_paths
    }


def build_full_frame(by_dt: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Wide per-zone frame: base metrics + derived residual-load columns."""
    price   = by_dt["prices"].iloc[:, 0].rename("price")
    load_fc = by_dt["load_forecast"].iloc[:, 0].rename("load_forecast")
    load    = by_dt["load_actual"].iloc[:, 0].rename("load")

    # res / generation / crossborder come MultiIndexed as (area, metric); the
    # inputs are single-zone, so drop the area level down to metric columns.
    res = by_dt["res"].copy();         res.columns = res.columns.droplevel(0)
    gen = by_dt["generation"].copy();  gen.columns = gen.columns.droplevel(0)
    xb  = by_dt["crossborder"].copy(); xb.columns  = xb.columns.droplevel(0)

    df = pd.concat([price, load_fc, load, res, gen, xb], axis=1)
    df = to_utc_naive(df).ffill(limit=3).fillna(0.0)

    # Derived forecasts (fix the original 'wind_onshore_foreacast' typo that
    # silently zeroed the onshore contribution).
    df["wind_forecast"]     = df.get("wind_onshore_forecast", 0) + df.get("wind_offshore_forecast", 0)
    df["res_forecast"]      = df["wind_forecast"] + df.get("solar_forecast", 0)
    df["residual_forecast"] = df["load_forecast"] - df["res_forecast"]
    # Derived actuals from generation.
    df["wind"]     = df.get("wind_onshore", 0) + df.get("wind_offshore", 0)
    df["res"]      = df["wind"] + df.get("solar", 0)
    df["residual"] = df["load"] - df["res"]

    return df.sort_index(axis=1)


def process_data(snakemake) -> None:
    window = slice(
        iso(snakemake.wildcards.start_date),
        f"{iso(snakemake.wildcards.end_date)} 23:00",
    )
    full = build_full_frame(read_inputs(list(snakemake.input))).loc[window]
    full.index.name = "time"

    out = Path(snakemake.output[0])
    out.parent.mkdir(parents=True, exist_ok=True)
    full.to_parquet(out, index=True)
    print(f"Wrote {out} ({full.shape[0]} rows × {full.shape[1]} cols)")


if __name__ == "__main__":
    process_data(snakemake)
