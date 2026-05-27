"""Extract a single bidding zone's day-ahead price series for one date range.

Mirrors the res_cf CF outputs: one value per hourly timestamp, "time" index.

Input:  resources/entsoe/{bidding_zone}/{data_type}.parquet  (only prices is used;
        the other data_types are still downloaded and available for ad-hoc work).
Output: resources/entsoe/{bidding_zone}_{start_date}_{end_date}.parquet
        single column "price" (€/MWh), tz-naive UTC index named "time".

ENTSO-E day-ahead prices are indexed in Europe/Brussels local time. atlite CF
series are UTC-naive, so we convert to UTC here; run.py then aligns prices to the
CF index directly.
"""

from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake


def iso(yyyymmdd: str) -> str:
    """YYYYMMDD -> YYYY-MM-DD."""
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def prices_input(input_paths: list[str]) -> Path:
    for p in input_paths:
        if Path(p).stem == "prices":
            return Path(p)
    raise FileNotFoundError(f"No 'prices' parquet among inputs: {list(input_paths)}")


def process_data(snakemake) -> None:
    start_date = snakemake.params.start_date
    end_date = snakemake.params.end_date

    df = pd.read_parquet(prices_input(list(snakemake.input)))
    prices = df.iloc[:, 0]  # single area column from download_entsoe

    if prices.index.tz is not None:
        prices.index = prices.index.tz_convert("UTC").tz_localize(None)
    prices = prices.sort_index()

    # End at 23:00 so the full final day of hourly data is included.
    prices = prices.loc[iso(start_date):f"{iso(end_date)} 23:00"]

    # Bridge small gaps (DST hours, occasional missing values); leave larger gaps
    # as NaN so run.py's coverage check fails loudly rather than silently.
    prices = prices.ffill(limit=3).bfill(limit=3)

    prices = prices.rename("price")
    prices.index.name = "time"

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    prices.to_frame().to_parquet(out_path, index=True)
    print(f"Wrote {out_path} ({len(prices)} hourly prices, "
          f"{prices.index.min()} .. {prices.index.max()})")


if __name__ == "__main__":
    process_data(snakemake)
