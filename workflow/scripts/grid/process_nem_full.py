"""Build one NEM region's full frame (grid variant=full).

Reads the raw NEM tables (price, load, generation, crossborder), keeps this
region's columns, and adds derived wind/residual-load columns. Native 5-min
resolution, raw AUD prices — for ad-hoc analysis, not consumed by the model.
"""

from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake


def read_raw(input_paths: list[str]) -> dict[str, pd.DataFrame]:
    """Map table name -> its (region, metric) DataFrame.

    Filenames are ``{nem_table}_{start_date}_{end_date}.parquet``; strip the two
    trailing date segments to recover the table key.
    """
    return {
        Path(p).stem.rsplit("_", 2)[0]: pd.read_parquet(p) for p in input_paths
    }


def process_data(snakemake) -> None:
    area = snakemake.wildcards.area
    raw = read_raw(list(snakemake.input))

    df = pd.concat(
        [raw["price"], raw["load"], raw["generation"], raw["crossborder"]], axis=1
    )
    if df.index.tz is not None:
        df.index = df.index.tz_localize(None)

    load_s = df[(area, "load")]
    wind = df.get((area, "wind_onshore"), pd.Series(0, index=df.index))
    solar = df.get((area, "solar"), pd.Series(0, index=df.index))
    derived = pd.DataFrame(index=df.index)
    derived["wind"] = wind  # NEM has no offshore, so wind == onshore
    derived["residual"] = load_s - (wind + solar)

    full = pd.concat([df[area], derived], axis=1)
    full.index.name = "time"

    out = Path(snakemake.output[0])
    out.parent.mkdir(parents=True, exist_ok=True)
    full.to_parquet(out, index=True)
    print(f"Wrote {out} ({full.shape[0]} rows × {full.shape[1]} cols)")


if __name__ == "__main__":
    process_data(snakemake)
