"""Unit tests for the ONS Brasilia-time window mapping.

ONS stamps din_instante in Brasilia time (UTC-3) while the rule output is
UTC-naive, so a UTC window's first hours come from the *previous* local month —
the mirror image of the forward pad NEM needs for AEST (UTC+10). Getting this
backwards silently drops the first three hours of every window, and at a year
boundary it reads from the wrong raw file entirely. Synthetic only, no download.
"""

import pandas as pd

import _helpers  # sys.path set by conftest
import retrieve_ons

_local_months = retrieve_ons._local_months
ONS_MARKET_TZ = _helpers.ONS_MARKET_TZ


def test_window_reaches_back_into_the_previous_local_month():
    # 2025-02-01 00:00 UTC is 2025-01-31 21:00 in Brasilia, so January is needed.
    months = _local_months("20250201", "20250228")
    assert months == ["2025-01", "2025-02"]


def test_year_boundary_pulls_the_previous_year():
    # The pad month lives in the previous raw year file, not just the previous month.
    months = _local_months("20250101", "20251231")
    assert months[0] == "2024-12"
    assert months[-1] == "2025-12"
    assert len(months) == 13


def test_mid_month_window_spans_only_its_own_month():
    assert _local_months("20250610", "20250620") == ["2025-06"]


def test_local_months_cover_the_whole_requested_utc_window():
    """Every hour of the UTC window must fall inside the local months selected."""
    start, end = "20250301", "20250331"
    months = _local_months(start, end)
    utc_hours = pd.date_range(f"2025-03-01", f"2025-03-31 23:00", freq="h")
    local = utc_hours.tz_localize("UTC").tz_convert(ONS_MARKET_TZ).tz_localize(None)
    assert set(local.to_period("M").strftime("%Y-%m")) <= set(months)


def test_pre_dst_abolition_window_is_refused():
    """Brazil observed DST until 2019; localising those timestamps is ambiguous."""
    months = _local_months("20180601", "20180630")
    assert int(months[0][:4]) < retrieve_ons.FIRST_DST_FREE_YEAR
