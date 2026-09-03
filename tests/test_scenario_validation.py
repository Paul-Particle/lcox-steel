"""Unit tests for the checks the scenario table has to pass before a DAG is built.

Both used to sit as raw Python in the Snakefile, where nothing could reach them.
Synthetic frames, no CSV on disk.
"""

import pandas as pd
import pytest

from common._runs import check_one_series_per_tech, check_run_coverage

COLUMNS = ["scenario", "route", "tech", "variant", "area", "start_date", "end_date"]


def _rows(*rows) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=COLUMNS)


def test_one_series_per_tech_accepts_distinct_techs():
    df = _rows(
        ("s", "all-routes", "wind-onshore", "bestsite-p95", "DEU", "20250101", "20251231"),
        ("s", "all-routes", "solar", "bestsite-p95", "DEU", "20250101", "20251231"),
    )
    check_one_series_per_tech(df)


def test_one_series_per_tech_rejects_a_tech_twice_in_one_run():
    """Two variants of one tech in one run would hand the solve two series for it."""
    df = _rows(
        ("s", "all-routes", "solar", "bestsite-p95", "DEU", "20250101", "20251231"),
        ("s", "all-routes", "solar", "area-average", "DEU", "20250101", "20251231"),
    )
    with pytest.raises(ValueError, match="duplicate rows for one run and tech"):
        check_one_series_per_tech(df)


def test_one_series_per_tech_allows_the_same_tech_in_another_area_or_year():
    df = _rows(
        ("s", "all-routes", "solar", "bestsite-p95", "DEU", "20250101", "20251231"),
        ("s", "all-routes", "solar", "bestsite-p95", "FRA", "20250101", "20251231"),
        ("s", "all-routes", "solar", "bestsite-p95", "DEU", "20240101", "20241231"),
    )
    check_one_series_per_tech(df)


def test_run_coverage_rejects_renewables_and_prices_that_never_meet():
    """The quiet mistake: a grid row filed under an area the CF rows never use."""
    df = _rows(
        ("s", "all-routes", "wind-onshore", "bestsite-p95", "AUS", "20250101", "20251231"),
        ("s", "all-routes", "grid", "dayahead", "VIC1", "20250101", "20251231"),
    )
    with pytest.raises(ValueError, match="never meet"):
        check_run_coverage(df)


def test_run_coverage_allows_one_area_islanded_and_another_priced():
    df = _rows(
        ("s", "all-routes", "wind-onshore", "bestsite-p95", "BRA", "20250101", "20251231"),
        ("s", "all-routes", "wind-onshore", "bestsite-p95", "DEU", "20250101", "20251231"),
        ("s", "all-routes", "grid", "dayahead", "DEU", "20250101", "20251231"),
    )
    check_run_coverage(df)
