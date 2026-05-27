"""Extract one NEM region's hourly day-ahead price (grid variant=dayahead).

Reads the raw DISPATCHPRICE table (resources/nem/raw/price_{start}_{end}.parquet),
takes the region's RRP, resamples the native 5-min series to hourly, and converts
AUD→EUR (fx param) so it matches the EUR-denominated ENTSO-E prices the
optimisation consumes.
"""

from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake


def process_data(snakemake) -> None:
    area = snakemake.wildcards.area
    eur_per_aud = snakemake.params.eur_per_aud

    # NEMOSIS writes tz-naive timestamps in NEM local time (AEST/AEDT, UTC+10/+11).
    prices = pd.read_parquet(snakemake.input[0])
    price = (prices[(area, "price")].resample("1h").mean() * eur_per_aud).rename("price")
    price.index.name = "time"

    out = Path(snakemake.output.prices)
    out.parent.mkdir(parents=True, exist_ok=True)
    price.to_frame().to_parquet(out, index=True)
    print(f"Wrote {out} ({len(price)} hourly prices, "
          f"{price.index.min()} .. {price.index.max()})")


if __name__ == "__main__":
    process_data(snakemake)
