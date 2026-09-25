"""Unit tests for the checks the scenario table has to pass before a DAG is built.

Both live in `common/_runs.py` rather than in the Snakefile, so a test can reach
them. Synthetic frames, no CSV on disk; the destination tests write their
assumptions into a temporary config directory.
"""

import pandas as pd
import pytest

from common._runs import (
    build_destination_frame,
    check_one_series_per_tech,
    check_run_coverage,
    top_level_areas,
)

COLUMNS = ["scenario", "route", "tech", "variant", "area", "start_date", "end_date"]
RUN_COLUMNS = ["scenario", "area", "start_date", "end_date", "route"]
AREAS = {
    "DEU": {"market": "entsoe", "market_area": "DE_LU"},
    "BRA": {},
    "VIC1": {"market": "nem", "market_area": "VIC1"},
    "ES_SYN": {"market": "synthetic"},
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


def _config_dir(tmp_path, base_area: str, overlays: dict | None = None):
    """A config directory holding the base assumptions and any scenario overlays."""
    (tmp_path / "assumptions.yaml").write_text(f"destination:\n  area: {base_area}\n")
    for scenario, area in (overlays or {}).items():
        (tmp_path / f"assumptions_{scenario}.yaml").write_text(
            f"destination:\n  area: {area}\n"
        )
    return tmp_path


def test_only_an_export_run_asks_for_the_destination_market(tmp_path):
    """A domestic route melts where it made its iron, so a scenario that builds
    none of the export twins needs no market download for the destination."""
    runs = _runs(
        ("s", "BRA", "20250101", "20251231", "moe-eaf"),
        ("s", "BRA", "20250101", "20251231", "moe-eaf-export"),
        ("s", "VIC1", "20250101", "20251231", "moe-eaf-export"),
    )
    destinations = build_destination_frame(runs, _config_dir(tmp_path, "DEU"), AREAS)

    # Both export runs melt in the same place, so one row and one download
    # covers them — where the iron came from is not part of the key.
    assert list(destinations["route"]) == ["moe-eaf-export"]
    assert list(destinations["destination"]) == ["DEU"]
    assert "area" not in destinations.columns


def test_an_overlay_moves_only_its_own_scenario_destination(tmp_path):
    """The demo melts in the synthetic market while every other scenario keeps
    the base destination, so the DAG fetches a different series for each."""
    runs = _runs(
        ("demo", "VIC1", "20250101", "20251231", "moe-eaf-export"),
        ("study", "VIC1", "20250101", "20251231", "moe-eaf-export"),
    )
    config_dir = _config_dir(tmp_path, "DEU", {"demo": "ES_SYN"})
    destinations = build_destination_frame(runs, config_dir, AREAS).set_index("scenario")

    assert destinations.loc["demo", "destination"] == "ES_SYN"
    assert destinations.loc["study", "destination"] == "DEU"


def test_a_destination_that_trades_in_no_market_is_an_error(tmp_path):
    """Its furnace has to buy power hour by hour somewhere, and Brazil has no
    series in the model — so say so at DAG time, not after a solve."""
    runs = _runs(("s", "VIC1", "20250101", "20251231", "moe-eaf-export"))
    with pytest.raises(ValueError, match="destination.area"):
        build_destination_frame(runs, _config_dir(tmp_path, "BRA"), AREAS)


def test_all_areas_never_reaches_the_synthetic_market():
    """Its series is made up, so a real scenario naming every area must not solve it."""
    assert "ES_SYN" not in top_level_areas(AREAS)
