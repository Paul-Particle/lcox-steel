"""Unit tests for NEM market-month cache membership.

NEM prices are downloaded by whole *market* months (AEST = UTC+10) but stored in
UTC. The ~10 h shift spills each market month across two UTC calendar months, so
`area_month_in_cache` must be told to match in market time (tz=NEM_MARKET_TZ) —
otherwise a neighbouring month's spillover (notably the trailing hours the pad
month contributes to the prior UTC month) makes it report a month as cached and
skip its download, dropping ~a month of data. Synthetic frames, no raw cache.
"""

import pandas as pd

import _helpers  # sys.path set by conftest

area_month_in_cache = _helpers.area_month_in_cache
NEM_MARKET_TZ = _helpers.NEM_MARKET_TZ


def _cache_for_market_months(area: str, months: list[str]) -> pd.DataFrame:
    """Build a UTC-stored cache holding exactly the given market (AEST) months."""
    frames = []
    for ym in months:
        aest = pd.date_range(
            f"{ym}-01 00:00", pd.Timestamp(f"{ym}-01") + pd.offsets.MonthEnd(0) + pd.Timedelta(hours=23),
            freq="h", tz=NEM_MARKET_TZ,
        )
        frames.append(pd.Series(1.0, index=aest.tz_convert("UTC").tz_localize(None)))
    cache = pd.concat(frames).to_frame("price")
    cache.columns = pd.MultiIndex.from_tuples([(area, "price")])
    return cache


def test_pad_month_spillover_does_not_mask_real_month():
    # Cache holds only market-January 2025 (as the pad month for a 2024 window).
    # Its UTC span reaches back into 2024-12, which must NOT make Dec look cached.
    cache = _cache_for_market_months("VIC1", ["2025-01"])
    assert area_month_in_cache(cache, "VIC1", "2024-12", tz=NEM_MARKET_TZ) is False
    assert area_month_in_cache(cache, "VIC1", "2025-01", tz=NEM_MARKET_TZ) is True


def test_real_market_month_is_detected():
    cache = _cache_for_market_months("VIC1", ["2024-12"])
    assert area_month_in_cache(cache, "VIC1", "2024-12", tz=NEM_MARKET_TZ) is True


def test_default_utc_behaviour_unchanged():
    # Without tz, matching is by UTC calendar month — the pre-existing behaviour
    # that ENTSO-E relies on. The Jan pad's UTC spillover DOES fall in 2024-12.
    cache = _cache_for_market_months("VIC1", ["2025-01"])
    assert area_month_in_cache(cache, "VIC1", "2024-12") is True


def test_missing_area_returns_false():
    cache = _cache_for_market_months("VIC1", ["2025-01"])
    assert area_month_in_cache(cache, "NSW1", "2025-01", tz=NEM_MARKET_TZ) is False
    assert area_month_in_cache(None, "VIC1", "2025-01", tz=NEM_MARKET_TZ) is False
