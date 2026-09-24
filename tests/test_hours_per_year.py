"""A full calendar year levelises over its own hours, leap year or not."""

import pandas as pd
import pytest

from scripts.viz import compile_report


@pytest.mark.parametrize("year", [2024, 2025])
def test_a_full_year_is_not_rescaled(year):
    snapshots = pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")
    assert compile_report._hours_per_year(snapshots) == len(snapshots)


def test_a_part_year_scales_to_its_calendar_year():
    january_2024 = pd.date_range("2024-01-01", "2024-01-31 23:00", freq="h")
    assert compile_report._hours_per_year(january_2024) == 8784.0
