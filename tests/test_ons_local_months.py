"""Unit tests for the ONS Brasilia-time window mapping.

ONS stamps din_instante in Brasilia time (UTC-3) while the rule output is
UTC-naive, so a UTC window's first hours come from the *previous* local month —
the mirror image of the forward pad NEM needs for AEST (UTC+10). Getting this
backwards silently drops the first three hours of every window, and at a year
boundary it reads from the wrong raw file entirely. Synthetic only, no download.
"""

from types import SimpleNamespace

import pandas as pd
import pytest

import _helpers_grid  # sys.path set by conftest
import _ons

_local_months = _ons._local_months
ONS_MARKET_TZ = _helpers_grid.ONS_MARKET_TZ


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


def _request(start: str, end: str) -> SimpleNamespace:
    """The rule inputs `retrieve` reads before it decides it cannot serve them."""
    return SimpleNamespace(
        wildcards=SimpleNamespace(variant="dayahead", start_date=start, end_date=end),
        params=SimpleNamespace(eur_per_brl=0.16, pld_limits={}),
    )


@pytest.mark.parametrize("start, end", [
    ("20180601", "20180630"),
    # Brazil's last DST period ran to 16 February 2019, so 2019 carries a
    # transition too and a February window lands on the repeated hour.
    ("20190201", "20190228"),
])
def test_a_window_reaching_into_a_dst_year_is_refused(start, end):
    """Localising a repeated wall-clock hour is ambiguous, so say so and stop."""
    with pytest.raises(ValueError, match="daylight saving"):
        _ons.retrieve(_request(start, end), "SE")


def test_the_first_clean_year_is_the_first_one_that_localises():
    """What FIRST_DST_FREE_YEAR has to be: the year before it still repeats an hour."""
    last_dst_year = _ons.FIRST_DST_FREE_YEAR - 1
    repeated = pd.date_range(f"{last_dst_year}-02-16 22:00", periods=6, freq="h")
    with pytest.raises(ValueError, match="ambiguous"):
        repeated.tz_localize(ONS_MARKET_TZ)

    clean = pd.date_range(f"{_ons.FIRST_DST_FREE_YEAR}-02-16 22:00", periods=6, freq="h")
    assert len(clean.tz_localize(ONS_MARKET_TZ)) == 6
