"""Extract one bidding zone's hourly day-ahead price (grid variant=dayahead).

Reads only the prices parquet and emits a single "price" column (€/MWh,
tz-naive UTC, hourly) for use as a model input. ENTSO-E indexes day-ahead
prices in Europe/Brussels local time while atlite CF series are UTC-naive, so
we convert here and run.py aligns prices to the CF index directly. The rich
multi-metric frame is built separately by process_entsoe_full.py (variant=full).
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


def process_data(snakemake) -> None:
    window = slice(
        iso(snakemake.wildcards.start_date),
        f"{iso(snakemake.wildcards.end_date)} 23:00",
    )
    prices_df = pd.read_parquet(snakemake.input[0])
    price = to_utc_naive(prices_df.iloc[:, 0])
    # Bridge small gaps but never zero-fill (a fake €0 would distort the
    # optimisation); larger gaps stay NaN so run.py fails loudly.
    # variant=dayahead ⇒ hourly: no-op for hourly zones, averages sub-hourly
    # ones (e.g. post-2025 15-min MTU markets).
    price = price.ffill(limit=3).bfill(limit=3).resample("1h").mean().rename("price")
    price = price.loc[window]
    price.index.name = "time"

    out = Path(snakemake.output.prices)
    out.parent.mkdir(parents=True, exist_ok=True)
    price.to_frame().to_parquet(out, index=True)
    print(f"Wrote {out} ({len(price)} hourly prices, "
          f"{price.index.min()} .. {price.index.max()})")


if __name__ == "__main__":
    process_data(snakemake)
