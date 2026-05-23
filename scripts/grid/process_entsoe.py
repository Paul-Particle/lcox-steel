"""Process per-(area, data_type) ENTSO-E feathers into a single hourly EU dataset.

Inputs: a flat list of resources/entsoe/{area}/{data_type}.feather paths
        (each file already covers the full date range for its area+data_type).
Output: resources/entsoe_processed.feather — wide DataFrame, columns = (area, metric).
"""

from collections import defaultdict
from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake


def load_per_data_type(input_paths: list[str]) -> dict[str, pd.DataFrame]:
    """Group inputs by data_type and concat areas horizontally per data_type."""
    by_dt = defaultdict(list)
    for p in input_paths:
        p = Path(p)
        by_dt[p.stem].append(p)

    result = {}
    for dt, paths in by_dt.items():
        dfs = [pd.read_feather(p) for p in sorted(paths)]
        result[dt] = pd.concat(dfs, axis=1).sort_index()
    return result


def process_data(snakemake) -> None:
    print("Processing data...")
    AREAS = list(snakemake.params.areas)

    by_dt = load_per_data_type(list(snakemake.input.entsoe))

    df_price          = by_dt["prices"].ffill(limit=3).fillna(0.0)
    df_load_forecast  = by_dt["load_forecast"].ffill(limit=3).fillna(0.0)
    df_load_actual    = by_dt["load_actual"].ffill(limit=3).fillna(0.0)
    df_res_forecast   = by_dt["res"].ffill(limit=3).fillna(0.0)
    df_generation     = by_dt["generation"].ffill(limit=3).fillna(0.0)
    df_crossborder    = by_dt["crossborder"].ffill(limit=3).fillna(0.0)

    # Wrap single-level frames in a (area, metric) MultiIndex so they can be merged
    # with the already-MultiIndex generation / VRE / crossborder frames.
    df_p = df_price.copy()
    df_p.columns = pd.MultiIndex.from_tuples([(col, 'price') for col in df_p.columns])

    df_d_fc = df_load_forecast.copy()
    df_d_fc.columns = pd.MultiIndex.from_tuples([(col, 'load_forecast') for col in df_d_fc.columns])

    df_d_ac = df_load_actual.copy()
    df_d_ac.columns = pd.MultiIndex.from_tuples([(col, 'load') for col in df_d_ac.columns])

    df_EU = pd.concat([df_p, df_d_fc, df_d_ac, df_res_forecast, df_generation, df_crossborder], axis=1)

    dfs_EU = {}
    for area in [col for col in df_EU.columns.get_level_values(0).unique() if col in AREAS]:
        df_area_data = df_EU[area].copy()

        df_area_processed = (df_area_data
            .sort_index(axis=1)
            .assign(wind_forecast = lambda x: x.get('wind_onshore_foreacast', 0) + x.get('wind_offshore_forecast', 0))
            .assign(res_forecast  = lambda x: x.get('wind_forecast', 0) + x.get('solar_forecast', 0))
            .assign(residual_forecast = lambda x: x.get('load_forecast', 0) - x.get('res_forecast', 0))
            .assign(wind = lambda x: x.get('wind_onshore', 0) + x.get('wind_offshore', 0))
            .assign(res  = lambda x: x.get('wind', 0) + x.get('solar', 0))
            .assign(residual = lambda x: x.get('load', 0) - x.get('res', 0))
            .ffill(limit=3)
        )

        df_area_processed.columns = pd.MultiIndex.from_product([[area], df_area_processed.columns])
        dfs_EU[area] = df_area_processed

    df_EU_all = pd.concat(dfs_EU.values(), axis=1)
    df_EU_all.index = df_EU_all.index.tz_localize(None)  # pyright: ignore[reportAttributeAccessIssue]
    df_EU_all.to_feather(snakemake.output[0])
    print(f"Successfully saved processed data to {snakemake.output[0]}")


if __name__ == "__main__":
    process_data(snakemake)
