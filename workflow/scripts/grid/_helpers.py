"""Shared helpers for retrieve_entsoe.py and retrieve_nem.py."""

import pandas as pd


def iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def iter_months_str(start_date: str, end_date: str) -> list[str]:
    """Return a list of 'YYYY-MM' strings for every month in [start_date, end_date].

    start_date / end_date are 'YYYYMMDD' strings (Snakemake wildcard format).
    """
    start = pd.Timestamp(iso(start_date))
    end = pd.Timestamp(iso(end_date))
    return [ts.strftime("%Y-%m") for ts in pd.date_range(start=start, end=end, freq="MS")]


def to_utc_naive(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a tz-aware or fixed-offset-naive (AEST = UTC+10) index to UTC-naive."""
    if df.index.tz is not None:
        df.index = df.index.tz_convert("UTC").tz_localize(None)
    else:
        df.index = df.index.tz_localize("Australia/Brisbane").tz_convert("UTC").tz_localize(None)
    return df.sort_index()


def area_month_in_cache(cached: pd.DataFrame | None, area: str, ym: str) -> bool:
    """Return True if `cached` already holds any data for (area, month `ym`)."""
    if cached is None or area not in cached.columns.get_level_values(0):
        return False
    area_data = cached[area].dropna(how="all")
    return (area_data.index.to_period("M") == pd.Period(ym, freq="M")).any()
