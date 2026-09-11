"""AESO (Alberta) and IESO (Ontario) source implementations, called by retrieve_grid_data.py.

Canada runs no national wholesale market. Alberta and Ontario each operate one and
publish hourly prices; every other province is a regulated cost-of-service utility
with no hourly price series to retrieve at all. Both operators serve a whole
calendar year in a single request, cached raw under data/canada_cache/.

Both label hour endings in their own market time, and the two conventions differ in
the way that matters here: IESO holds EST all year, so every market day is 24 hours
and a fixed offset is the entire conversion, while AESO follows Mountain daylight
saving, so its days run 23, 24 or 25 hours and the repeated autumn hour arrives
labelled `02*`.

Variants
--------
dayahead  single "price" column, hourly UTC-naive, EUR/MWh
"""

import logging
from pathlib import Path

import pandas as pd
import requests

from _helpers_grid import assert_window_complete, iso

# Module-level logger only — retrieve_grid_data.py installs the handlers.
log = logging.getLogger(__name__)

AESO_POOL_PRICE_URL = (
    "http://ets.aeso.ca/ets_web/ip/Market/Reports/HistoricalPoolPriceReportServlet"
    "?contentType=csv&beginDate={begin}&endDate={end}"
)
IESO_HOEP_URL = (
    "https://reports-public.ieso.ca/public/PriceHOEPPredispOR/"
    "PUB_PriceHOEPPredispOR_{year}.csv"
)

AESO_MARKET_TZ = "America/Edmonton"
# IESO settles in EST and holds it all year, so no daylight saving applies and a
# fixed offset is the whole of the conversion to UTC.
IESO_MARKET_OFFSET = pd.Timedelta(hours=5)


def fetch_year_csv(url: str, cache_path: Path) -> Path:
    """Download `url` to `cache_path` unless it is already cached, and return the path."""
    if cache_path.exists():
        return cache_path
    log.info(f"downloading {url}")
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(response.content)
    return cache_path


def read_aeso_year(year: int, cache_dir: Path) -> pd.DataFrame:
    """One year of Alberta pool price and internal load, hourly and UTC-naive."""
    # endDate is exclusive: asking for 1 January of the next year is what returns
    # 31 December. Requesting 31 December stops the file a day short instead.
    url = AESO_POOL_PRICE_URL.format(begin=f"0101{year}", end=f"0101{year + 1}")
    raw = (
        pd.read_csv(fetch_year_csv(url, cache_dir / f"aeso_pool_price_{year}.csv"),
                    skiprows=4)
        .dropna(how="all")
    )
    raw[["market_date", "hour_ending"]] = raw["Date (HE)"].str.split(" ", expand=True)

    # The hour-ending label is wall-clock, so it skips an hour in spring (01 → 03)
    # and repeats one in autumn (02, 02*). Arithmetic on the label therefore lands
    # on hours that do not exist. Each day's rows are contiguous, so anchoring the
    # day at its own local midnight and stepping an hour per row is correct in both
    # directions. Labels are zero-padded, which puts `02*` straight after `02`.
    ordered = raw.sort_values(["market_date", "hour_ending"])
    day_start_utc = (
        pd.to_datetime(ordered["market_date"], format="%m/%d/%Y")
        .dt.tz_localize(AESO_MARKET_TZ)
        .dt.tz_convert("UTC")
        .dt.tz_localize(None)
    )
    hours_into_day = ordered.groupby("market_date").cumcount()

    prices_and_load = (
        ordered
        .set_axis(day_start_utc + pd.to_timedelta(hours_into_day, unit="h"))
        .rename(columns={"Price ($)": "price", "AIL Demand (MW)": "load"})
        .loc[:, ["price", "load"]]
        .apply(pd.to_numeric, errors="coerce")
        .sort_index()
    )
    return prices_and_load


def read_ieso_year(year: int, cache_dir: Path) -> pd.DataFrame:
    """One year of Ontario hourly energy price (HOEP), hourly and UTC-naive."""
    # Three comment lines carrying the report title sit above the header row.
    raw = pd.read_csv(
        fetch_year_csv(IESO_HOEP_URL.format(year=year), cache_dir / f"ieso_hoep_{year}.csv"),
        skiprows=3,
    )
    market_hour_start = (
        pd.to_datetime(raw["Date"]) + pd.to_timedelta(raw["Hour"] - 1, unit="h")
    )
    prices = (
        raw
        .set_axis(market_hour_start + IESO_MARKET_OFFSET)
        .rename(columns={"HOEP": "price"})
        .loc[:, ["price"]]
        .apply(pd.to_numeric, errors="coerce")
        .sort_index()
    )
    return prices


def retrieve_window(snakemake, read_year, eur_per_cad: float) -> None:
    """Assemble the requested window from annual files, convert CAD to EUR and write it."""
    variant = snakemake.wildcards.variant
    start_date = snakemake.wildcards.start_date
    end_date = snakemake.wildcards.end_date

    if variant != "dayahead":
        raise ValueError(
            f"variant {variant!r} is not available for the Canadian markets — only "
            f"'dayahead'. Per-carrier generation needs AESO's keyed API or its bulk "
            f"metered-volume files, and for Ontario a separate XML report; neither "
            f"is wired up."
        )

    # Both markets sit west of UTC, so a UTC window opens during the *previous*
    # local year: Alberta at UTC-7/-6 and Ontario at UTC-5 each leave the window's
    # first hours in the year before. (NEM pads the far end instead, being east.)
    years = range(int(start_date[:4]) - 1, int(end_date[:4]) + 1)
    cache_dir = Path("data/canada_cache")
    assembled = pd.concat([read_year(year, cache_dir) for year in years]).sort_index()

    window = slice(iso(start_date), f"{iso(end_date)} 23:59")
    out_df = assembled.loc[window, ["price"]].copy()
    # Both markets settle in CAD; every variant's `price` column is EUR/MWh.
    out_df["price"] = out_df["price"] * eur_per_cad
    out_df.index.name = "time"

    assert_window_complete(out_df, start_date, end_date, variant)

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")


def retrieve_aeso(snakemake, area: str) -> None:
    """Alberta pool price for the requested window. AESO serves one province, so `area` is AB."""
    retrieve_window(snakemake, read_aeso_year, snakemake.params.eur_per_cad)


def retrieve_ieso(snakemake, area: str) -> None:
    """Ontario HOEP for the requested window. Pre-MRP province-wide price, so `area` is ON."""
    retrieve_window(snakemake, read_ieso_year, snakemake.params.eur_per_cad)
