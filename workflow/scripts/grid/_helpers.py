"""Shared helpers for retrieve_entsoe.py and retrieve_nem.py."""

import pandas as pd


def iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def to_utc_naive(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a tz-aware or fixed-offset-naive (AEST = UTC+10) index to UTC-naive."""
    if df.index.tz is not None:
        df.index = df.index.tz_convert("UTC").tz_localize(None)
    else:
        df.index = df.index.tz_localize("Australia/Brisbane").tz_convert("UTC").tz_localize(None)
    return df.sort_index()


def area_month_in_cache(cached: pd.DataFrame | None, area: str, ym: str) -> bool:
    if cached is None or area not in cached.columns.get_level_values(0):
        return False
    area_data = cached[area].dropna(how="all")
    return (area_data.index.to_period("M") == pd.Period(ym, freq="M")).any()
