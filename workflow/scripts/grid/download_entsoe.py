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
    """Build an authenticated ENTSO-E client from ENTSOE_API_KEY (loaded via .env)."""
    load_dotenv()
    api_key = os.environ.get("ENTSOE_API_KEY")
    if not api_key:
        raise ValueError(
            "ENTSOE_API_KEY environment variable not set. Please set it in your .env file."
        )
    return entsoe.EntsoePandasClient(api_key=api_key)  # pyright: ignore[reportPrivateImportUsage]


# ── Per-data_type fetchers (return DataFrame or raise) ────────────────────────

def download_prices(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    """Fetch day-ahead prices for one area as a single-column (€/MWh) frame."""
    data = client.query_day_ahead_prices(area, start=start, end=end)
    data.name = area
    return data.to_frame()


def download_load_forecast(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    """Fetch the day-ahead load forecast for one area (single column, MW)."""
    data = client.query_load_forecast(area, start=start, end=end)
    data.columns = [area]
    return data


def download_load_actual(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    """Fetch actual (realised) load for one area (single column, MW)."""
    data = client.query_load(area, start=start, end=end)
    data.columns = [area]
    return data


def download_res(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    """Fetch the wind/solar forecast, renamed to solar/wind_onshore/wind_offshore_forecast."""
    data = client.query_wind_and_solar_forecast(area, start=start, end=end)
    res_names = {
        "Solar": "solar_forecast",
        "Wind Onshore": "wind_onshore_forecast",
        "Wind Offshore": "wind_offshore_forecast",
    }
    data.columns = pd.MultiIndex.from_tuples([(area, res_names[c]) for c in data.columns])
    return data


def download_generation(client: entsoe.EntsoePandasClient, area: str, start: pd.Timestamp, end: pd.Timestamp) -> pd.DataFrame:
    """Fetch per-carrier generation, renamed to short keys with consumption sign-flipped.

    ENTSO-E's verbose carrier labels are mapped to short keys (e.g. 'Fossil Gas'
    → 'gas'); 'Actual Consumption' columns are negated and suffixed `_cons` so
    storage/pumped-hydro consumption reads as negative generation.
    """
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
    """Fetch physical cross-border flows as signed columns (imports +, exports −)."""
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
    """Yield (YYYY-MM, month_start, next_month_start) for every month overlapping [start, end].

    ENTSO-E treats the query end as exclusive, so each month is bounded by the
    *next* month's start to capture its full final day. The range starts at the
    first of start_date's month so a mid-month start_date isn't skipped.
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
    """Call `fetcher`, retrying up to max_attempts with exponential backoff."""
    for attempt in range(1, max_attempts + 1):
        try:
            return fetcher(client, area, start=start, end=end)
        except Exception as e:
            if attempt == max_attempts:
                raise
            backoff = 2 ** attempt
            log.warning(
                f"attempt {attempt}/{max_attempts} failed ({e!r}); retry in {backoff}s"
            )
            time.sleep(backoff)


