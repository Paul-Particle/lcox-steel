"""Unit tests for zone ranking in the report and the run labels the plots use.

A country that reaches its market through zones (Australia through its NEM
regions) is solved once per zone, so several report rows describe the same
place. `mark_best_in_country` flags the cheapest rather than dropping the rest,
and the plots hide the unflagged ones. These pin that behaviour on synthetic
frames — no cutouts, no solves, so they run anywhere. Only VIC1 is configured
today, so the multi-zone case has to be synthetic to exist at all.
"""

import pandas as pd
import pytest

import compile_report  # sys.path set by conftest
import _run_display
from common._runs import zone_parents

AREAS = {
    "DEU": {"market": "entsoe", "market_area": "DE_LU"},
    "BRA": {},
    "AUS": {"zones": ["VIC1", "NSW1", "QLD1"]},
    "VIC1": {"market": "nem", "market_area": "VIC1"},
    "NSW1": {"market": "nem", "market_area": "NSW1"},
    "QLD1": {"market": "nem", "market_area": "QLD1"},
}
PARENTS = zone_parents(AREAS)


def _report(rows):
    """A report-shaped frame: one row per run, in the column order compile_report uses."""
    return pd.DataFrame(
        rows, columns=["scenario", "area", "route", "start_date", "end_date", "lco_output"]
    )


def test_zone_parents_maps_zones_to_their_country():
    assert PARENTS["VIC1"] == "AUS"
    assert PARENTS["NSW1"] == "AUS"
    assert PARENTS["DEU"] == "DEU", "an area that is nobody's zone is its own country"


def test_cheapest_zone_of_each_country_is_flagged():
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
        ("s", "QLD1", "moe-eaf", "20250101", "20251231", 875.0),
        ("s", "DEU", "moe-eaf", "20250101", "20251231", 1010.0),
    ])
    out = compile_report.mark_best_in_country(df, PARENTS, "lco_output")

    flagged = set(out.loc[out["best_in_country"], "area"])
    assert flagged == {"NSW1", "DEU"}, "the cheapest Australian zone, and Germany on its own"
    assert list(out["country"]) == ["AUS", "AUS", "AUS", "DEU"]
    assert len(out) == len(df), "losers are flagged, never dropped"


def test_ranking_is_per_route_and_per_date_range():
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
        ("s", "VIC1", "ew-eaf", "20250101", "20251231", 700.0),
        ("s", "NSW1", "ew-eaf", "20250101", "20251231", 760.0),
        ("s", "VIC1", "moe-eaf", "20240101", "20241231", 880.0),
        ("s", "NSW1", "moe-eaf", "20240101", "20241231", 910.0),
    ])
    out = compile_report.mark_best_in_country(df, PARENTS, "lco_output")
    winners = {
        (r.route, r.start_date): r.area for r in out.itertuples() if r.best_in_country
    }
    # Each route picks its own zone, and each year is ranked separately.
    assert winners == {
        ("moe-eaf", "20250101"): "NSW1",
        ("ew-eaf", "20250101"): "VIC1",
        ("moe-eaf", "20240101"): "VIC1",
    }


def test_runs_without_the_metric_are_kept():
    """h2-only reports a cost of hydrogen, not steel; it must not be ranked away."""
    df = _report([
        ("s", "VIC1", "h2-only", "20250101", "20251231", float("nan")),
        ("s", "NSW1", "h2-only", "20250101", "20251231", float("nan")),
    ])
    out = compile_report.mark_best_in_country(df, PARENTS, "lco_output")
    assert out["best_in_country"].all()


def test_no_metric_configured_flags_everything():
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
    ])
    out = compile_report.mark_best_in_country(df, PARENTS, None)
    assert out["best_in_country"].all()
    assert "country" in out.columns, "the country is reported either way"


@pytest.mark.parametrize(
    "row, expected",
    [
        ({"area": "NSW1", "route": "moe-eaf", "start_date": "20250101", "end_date": "20251231"},
         "NSW1 moe-eaf"),
        ({"area": "DEU", "route": "h2-dri-eaf", "start_date": "20230101", "end_date": "20241231"},
         "DEU h2-dri-eaf 2023-2024"),
    ],
)
def test_run_label_identifies_the_run(row, expected):
    assert _run_display.run_label(pd.Series(row)) == expected


def test_run_labels_are_unique_across_a_multi_area_scenario():
    """The bug this guards: a scenario spanning areas rendered every one as 'moe-eaf'."""
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
        ("s", "DEU", "moe-eaf", "20250101", "20251231", 1010.0),
    ])
    labels = df.apply(_run_display.run_label, axis=1)
    assert labels.is_unique


def test_plots_hide_the_beaten_zones():
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
        ("s", "DEU", "moe-eaf", "20250101", "20251231", 1010.0),
    ])
    out = compile_report.mark_best_in_country(df, PARENTS, "lco_output")
    kept = _run_display.best_zones_only(out)
    assert set(kept["area"]) == {"NSW1", "DEU"}
