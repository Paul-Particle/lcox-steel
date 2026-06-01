"""ENTSO-E download primitives — imported by retrieve_entsoe.py.

Provides the per-month raw cache layer:
  data/entsoe_cache/{area}/{YYYY-MM}/{data_type}.parquet

To force a full re-fetch, delete the relevant month directories and re-run.
"""

import logging
import os
import time
from typing import Callable

import entsoe
import pandas as pd
from dotenv import load_dotenv

# Module-level logger only — the rule script (retrieve_entsoe.py) installs handlers.
log = logging.getLogger(__name__)


def get_entsoe_client() -> entsoe.EntsoePandasClient:
    load_dotenv()
    api_key = os.environ.get("ENTSOE_API_KEY")
    if not api_key:
        raise ValueError(
            "ENTSOE_API_KEY environment variable not set. Please set it in your .env file."
        )
    return entsoe.EntsoePandasClient(api_key=api_key)  # pyright: ignore[reportPrivateImportUsage]


# ── Per-data_type fetchers (return DataFrame or raise) ────────────────────────

def download_prices(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_day_ahead_prices(area, start=start, end=end)
    data.name = area
    return data.to_frame()


def download_load_forecast(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_load_forecast(area, start=start, end=end)
    data.columns = [area]
    return data


def download_load_actual(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_load(area, start=start, end=end)
    data.columns = [area]
    return data


def download_res(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_wind_and_solar_forecast(area, start=start, end=end)
    res_names = {
        "Solar": "solar_forecast",
        "Wind Onshore": "wind_onshore_forecast",
        "Wind Offshore": "wind_offshore_forecast",
    }
    data.columns = pd.MultiIndex.from_tuples([(area, res_names[c]) for c in data.columns])
    return data


def download_generation(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    data = client.query_generation(area, start=start, end=end)
    if data.columns.nlevels == 2:
        data.columns = ["_".join(col) for col in data.columns]
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


def download_crossborder(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
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


DOWNLOADERS = {
    "prices": download_prices,
    "load_forecast": download_load_forecast,
    "load_actual": download_load_actual,
    "res": download_res,
    "generation": download_generation,
    "crossborder": download_crossborder,
}


# ── Month iteration & retry helpers ──────────────────────────────────────────

def iter_months(start_date_str: str, end_date_str: str):
    """Yield (YYYY-MM, month_start_ts, next_month_start_ts) for each month that
    overlaps [start, end].

    The query end is the *next* month's start: ENTSO-E treats the end as
    exclusive, so this fetches the full final day of each month (using MonthEnd
    landed on the last day at 00:00, dropping the last day). Starting the range
    at the first of start's month keeps a mid-month start_date from skipping it.
    """
    start = pd.to_datetime(start_date_str, format="%Y%m%d")
    end = pd.to_datetime(end_date_str, format="%Y%m%d")
    first_month = start.to_period("M").to_timestamp()
    for ts in pd.date_range(first_month, end, freq="MS"):
        ym = ts.strftime("%Y-%m")
        month_start = pd.Timestamp(ts, tz="Europe/Brussels")
        next_month_start = month_start + pd.offsets.MonthBegin(1)
        yield ym, month_start, next_month_start


def download_with_retry(
    fetcher: Callable[..., pd.DataFrame],
    client: entsoe.EntsoePandasClient,
    area: str,
    start: pd.Timestamp,
    end: pd.Timestamp,
    max_attempts: int = 3,
) -> pd.DataFrame:
    for attempt in range(1, max_attempts + 1):
        try:
            return fetcher(client, area, start=start, end=end)
        except Exception as e:
            if attempt == max_attempts:
                raise
            backoff = 2 ** attempt
            log.warning(
                "attempt %d/%d failed (%r); retry in %ds",
                attempt, max_attempts, e, backoff,
            )
            time.sleep(backoff)


