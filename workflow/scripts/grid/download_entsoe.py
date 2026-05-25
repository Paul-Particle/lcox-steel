"""Download ENTSO-E data for a single (area, data_type) pair.

Manages a per-month cache in data/entsoe_cache/{area}/{YYYY-MM}/{data_type}.feather.
Cached months are skipped; missing months are fetched (with 3 retries + exponential
backoff per month). Per-month failures are logged and the run continues; the rule
fails only if ZERO months succeed. To force a refresh, delete the relevant cache
files then re-run snakemake.
"""

import logging
import os
import time
from pathlib import Path
from typing import Callable

import entsoe
import pandas as pd
from dotenv import load_dotenv

if "snakemake" not in globals():
    from common._stubs import snakemake

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("download_entsoe")


def get_entsoe_client() -> entsoe.EntsoePandasClient:
    load_dotenv()
    api_key = os.environ.get("ENTSOE_API_KEY")
    if not api_key:
        raise ValueError(
            "ENTSOE_API_KEY environment variable not set. Please set it in your .env file."
        )
    return entsoe.EntsoePandasClient(api_key=api_key)  # pyright: ignore[reportPrivateImportUsage]


# ── Per-data_type fetchers (return DataFrame or raise) ────────────────────────

def fetch_prices(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_day_ahead_prices(area, start=start, end=end)
    data.name = area
    return data.to_frame()


def fetch_load_forecast(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_load_forecast(area, start=start, end=end)
    data.columns = [area]
    return data


def fetch_load_actual(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_load(area, start=start, end=end)
    data.columns = [area]
    return data


def fetch_res(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_wind_and_solar_forecast(area, start=start, end=end)
    res_names = {
        "Solar": "solar_forecast",
        "Wind Onshore": "wind_onshore_forecast",
        "Wind Offshore": "wind_offshore_forecast",
    }
    data.columns = pd.MultiIndex.from_tuples([(area, res_names[c]) for c in data.columns])
    return data


def fetch_generation(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_generation(area, start=start, end=end)
    if data.columns.nlevels == 2:
        data.columns = ["_".join(col) for col in data.columns.values]
    else:
        data.columns = ["_".join([col, "Actual Aggregated"]) for col in data.columns]
    cons_cols = data.filter(regex="Consumption$", axis=1).columns
    data.loc[:, cons_cols] = data.loc[:, cons_cols] * -1

    gen_names = {
        "Biomass": "biomass",
        "Energy storage": "energy_storage",
        "Fossil Brown coal/Lignite": "brown_coal",
        "Fossil Coal-derived gas": "coal_gas",
        "Fossil Gas": "gas",
        "Fossil Hard coal": "hard_coal",
        "Fossil Oil": "oil",
        "Fossil Oil shale": "oil_shale",
        "Fossil Peat": "peat",
        "Geothermal": "geothermal",
        "Hydro Pumped Storage": "pumped_storage",
        "Hydro Run-of-river and poundage": "hydro_river",
        "Hydro Water Reservoir": "hydro_reservoir",
        "Marine": "marine",
        "Nuclear": "nuclear",
        "Other": "other",
        "Other renewable": "other_re",
        "Solar": "solar",
        "Waste": "waste",
        "Wind Offshore": "wind_offshore",
        "Wind Onshore": "wind_onshore",
    }
    rename_map = {k + "_Actual Aggregated": v for k, v in gen_names.items()}
    rename_map |= {k + "_Actual Consumption": v + "_cons" for k, v in gen_names.items()}

    data.columns = pd.MultiIndex.from_tuples([(area, rename_map[c]) for c in data.columns])
    return data


def fetch_crossborder(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    parts = []
    data_in = client.query_physical_crossborder_allborders(
        area, start=start, end=end, export=False, per_hour=False
    )
    data_in.columns = ["from_" + c for c in data_in.columns]
    parts.append(data_in)
    data_out = client.query_physical_crossborder_allborders(
        area, start=start, end=end, export=True, per_hour=False
    )
    data_out.columns = ["to_" + c for c in data_out.columns]
    parts.append(data_out * -1)
    df = pd.concat(parts, axis=1)
    df.columns = pd.MultiIndex.from_tuples([(area, c) for c in df.columns])
    return df


FETCHERS = {
    "prices": fetch_prices,
    "load_forecast": fetch_load_forecast,
    "load_actual": fetch_load_actual,
    "res": fetch_res,
    "generation": fetch_generation,
    "crossborder": fetch_crossborder,
}


# ── Month iteration & retry helpers ──────────────────────────────────────────

def iter_months(start_date_str: str, end_date_str: str):
    """Yield (YYYY-MM, month_start_ts, month_end_ts) for each month in [start, end]."""
    start = pd.to_datetime(start_date_str, format="%Y%m%d")
    end = pd.to_datetime(end_date_str, format="%Y%m%d")
    for ts in pd.date_range(start, end, freq="MS"):
        ym = ts.strftime("%Y-%m")
        month_start = pd.Timestamp(ts, tz="Europe/Brussels")
        month_end = month_start + pd.offsets.MonthEnd(0)
        yield ym, month_start, month_end


def fetch_with_retry(
    fetcher: Callable[..., pd.DataFrame],
    client: entsoe.EntsoePandasClient,
    area: str,
    start: pd.Timestamp,
    end: pd.Timestamp,
    max_attempts: int = 3,
) -> pd.DataFrame | None:
    for attempt in range(1, max_attempts + 1):
        try:
            return fetcher(client, area, start=start, end=end)
        except Exception as e:
            if attempt == max_attempts:
                raise
            backoff = 2 ** attempt
            log.warning(
                f"  attempt {attempt}/{max_attempts} failed ({e!r}); retry in {backoff}s"
            )
            time.sleep(backoff)


# ── Main ─────────────────────────────────────────────────────────────────────

def download_data(snakemake) -> None:
    area = snakemake.wildcards.area
    data_type = snakemake.wildcards.data_type
    output_path = Path(snakemake.output[0])
    cache_dir = Path(snakemake.params.cache_dir)
    start_date = snakemake.params.start_date
    end_date = snakemake.params.end_date

    if data_type not in FETCHERS:
        raise ValueError(f"unknown data_type {data_type!r}")
    fetcher = FETCHERS[data_type]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    area_cache_dir = cache_dir / area

    client = None  # lazily instantiated; warm cache = no API call
    successful_dfs = []
    n_total = n_cached = n_fetched = n_failed = 0

    for ym, month_start, month_end in iter_months(start_date, end_date):
        n_total += 1
        month_cache_dir = area_cache_dir / ym
        cache_path = month_cache_dir / f"{data_type}.parquet"

        if cache_path.exists():
            successful_dfs.append(pd.read_parquet(cache_path))
            n_cached += 1
            continue

        if client is None:
            client = get_entsoe_client()

        try:
            log.info(f"{area}/{ym}/{data_type}: fetching")
            df = fetch_with_retry(fetcher, client, area, month_start, month_end)
        except Exception as e:
            log.error(f"{area}/{ym}/{data_type}: FAILED after retries — {e!r}")
            n_failed += 1
            continue

        month_cache_dir.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache_path, index=True)
        successful_dfs.append(df)
        n_fetched += 1

    if not successful_dfs:
        raise RuntimeError(
            f"{area}/{data_type}: zero months succeeded — refusing to write empty output. "
            f"(total={n_total}, failed={n_failed})"
        )

    combined = pd.concat(successful_dfs, axis=0)
    combined = combined[~combined.index.duplicated(keep="last")]
    combined = combined.sort_index()
    combined.to_parquet(output_path, index=True)
    log.info(
        f"{area}/{data_type}: wrote {output_path} "
        f"(cached={n_cached}, fetched={n_fetched}, failed={n_failed}, total={n_total})"
    )


if __name__ == "__main__":
    download_data(snakemake)
