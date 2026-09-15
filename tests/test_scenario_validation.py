"""Unit tests for the checks the scenario table has to pass before a DAG is built.

Both used to sit as raw Python in the Snakefile, where nothing could reach them.
Synthetic frames, no CSV on disk.
"""

import pandas as pd
import pytest

from common._runs import (
    build_destination_frame,
    check_one_series_per_tech,
    check_run_coverage,
)

COLUMNS = ["scenario", "route", "tech", "variant", "area", "start_date", "end_date"]
RUN_COLUMNS = ["scenario", "area", "start_date", "end_date", "route"]
AREAS = {
    "DEU": {"market": "entsoe", "market_area": "DE_LU"},
    "BRA": {},
    "VIC1": {"market": "nem", "market_area": "VIC1"},
}


def _rows(*rows) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=COLUMNS)


def _runs(*rows) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=RUN_COLUMNS)


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


def test_only_an_export_run_asks_for_the_destination_market():
    """A domestic route melts where it made its iron, so a scenario that builds
    none of the export twins needs no market download for the destination."""
    runs = _runs(
        ("s", "BRA", "20250101", "20251231", "moe-eaf"),
        ("s", "BRA", "20250101", "20251231", "moe-eaf-export"),
        ("s", "VIC1", "20250101", "20251231", "moe-eaf-export"),
    )
    destinations = build_destination_frame(runs, "DEU", AREAS)

    # Both export runs melt in the same place, so one row and one download
    # covers them — where the iron came from is not part of the key.
    assert list(destinations["route"]) == ["moe-eaf-export"]
    assert list(destinations["destination"]) == ["DEU"]
    assert "area" not in destinations.columns


def test_a_destination_that_trades_in_no_market_is_an_error():
    """Its furnace has to buy power hour by hour somewhere, and Brazil has no
    series in the model — so say so at DAG time, not after a solve."""
    runs = _runs(("s", "VIC1", "20250101", "20251231", "moe-eaf-export"))
    with pytest.raises(ValueError, match="destination.area"):
        build_destination_frame(runs, "BRA", AREAS)
