"""Unit tests for what a report row says and which rows survive into it.

A country that reaches its market through zones (Australia through its NEM
regions) is solved once per zone, so several report rows describe the same
place. `mark_best_in_country` flags the cheapest rather than dropping the rest,
and `write_report` splits that into the report viz reads and the hidden
diagnostic. `input_variants` puts back the one part of a run's identity the
solved network does not carry. These pin that behaviour on synthetic frames —
no cutouts, no solves, so they run anywhere.
"""

import pandas as pd
import pytest

import compile_report  # sys.path set by conftest
import _run_display
from common._report_schema import FIELD_ORDER, IDENTITY_FIELDS, read_report
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
    frame = pd.DataFrame(
        rows, columns=["scenario", "area", "route", "start_date", "end_date", "lco_output"]
    )
    # Every row carries the fingerprint of the inputs behind it, so a frame
    # without one is not a report frame.
    frame["inputs_hash"] = "772c9cac81fe"
    return frame


def test_a_run_is_a_column_and_a_field_is_a_row(tmp_path):
    """The file's shape: runs across the top, the identity fields leading."""
    df = _report([
        ("s", "VIC1", "h2-dri-eaf", "20250101", "20251231", 606.0),
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 914.0),
    ])
    flagged = compile_report.mark_best_in_country(df, PARENTS, "lco_output")

    report_path = tmp_path / "report_s.csv"
    compile_report.write_report(flagged, report_path, tmp_path / ".report_s_diag.csv")
    stored = pd.read_csv(report_path)

    assert list(stored.columns) == ["field", "s_1", "s_2"]
    leading = list(stored["field"])[:len(IDENTITY_FIELDS) - 1]
    # Every identity field but the diagnostic's flag, in declared order, on top.
    assert leading == [f for f in IDENTITY_FIELDS if f != "best_in_country"]
    assert stored.set_index("field").at["route", "s_2"] == "moe-eaf"


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


def test_report_holds_the_selected_areas_and_the_diagnostic_holds_all(tmp_path):
    """The seam: the report stands alone, the diagnostic answers "and the others?"."""
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
        ("s", "DEU", "moe-eaf", "20250101", "20251231", 1010.0),
    ])
    flagged = compile_report.mark_best_in_country(df, PARENTS, "lco_output")
    report_path = tmp_path / "report_s.csv"
    diagnostic_path = tmp_path / ".report_s_diag.csv"
    compile_report.write_report(flagged, report_path, diagnostic_path)

    report = read_report(report_path)
    assert set(report["area"]) == {"NSW1", "DEU"}
    assert "best_in_country" not in report.columns, "the flag is not viz's problem"

    diagnostic = read_report(diagnostic_path)
    assert set(diagnostic["area"]) == {"VIC1", "NSW1", "DEU"}
    assert set(diagnostic.loc[~diagnostic["best_in_country"], "area"]) == {"VIC1"}


def _scenario_rows(rows):
    """A scenario table as `load_scenarios` returns it: one row per input series."""
    return pd.DataFrame(
        rows, columns=["scenario", "tech", "variant", "area", "start_date", "end_date"]
    )


SCENARIOS = _scenario_rows([
    ("s", "wind-onshore", "bestsite-p95", "VIC1", "20250101", "20251231"),
    ("s", "solar",        "tilt-mix-n5",  "VIC1", "20250101", "20251231"),
    ("s", "grid",         "dayahead",     "VIC1", "20250101", "20251231"),
    ("s", "solar",        "area-average", "VIC1", "20240101", "20241231"),
    ("s", "solar",        "area-average", "NSW1", "20250101", "20251231"),
])


def test_input_variants_name_the_series_behind_a_run():
    run = {"area": "VIC1", "route": "moe-eaf",
           "start_date": "20250101", "end_date": "20251231"}
    assert compile_report.input_variants(SCENARIOS, "s", run) == {
        "wind-onshore_variant": "bestsite-p95",
        "solar_variant": "tilt-mix-n5",
        "grid_variant": "dayahead",
    }


def test_input_variants_ignore_the_route():
    """Every route of a group is solved from the same series, so they all report it."""
    group = {"area": "VIC1", "start_date": "20250101", "end_date": "20251231"}
    moe = compile_report.input_variants(SCENARIOS, "s", {**group, "route": "moe-eaf"})
    h2 = compile_report.input_variants(SCENARIOS, "s", {**group, "route": "h2-dri-eaf"})
    assert moe == h2


@pytest.mark.parametrize(
    "run, expected",
    [
        ({"area": "VIC1", "start_date": "20240101", "end_date": "20241231"},
         {"solar_variant": "area-average"}),
        ({"area": "NSW1", "start_date": "20250101", "end_date": "20251231"},
         {"solar_variant": "area-average"}),
    ],
)
def test_input_variants_do_not_leak_across_areas_or_years(run, expected):
    assert compile_report.input_variants(SCENARIOS, "s", run) == expected


def test_every_route_writes_the_same_fields(tmp_path):
    """The point of the schema: what a run reports must not depend on its route."""
    df = _report([
        ("s", "VIC1", "h2-only", "20250101", "20251231", 4.1),
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 914.0),
    ])
    # Only the MOE run has a MOE cell to report; only the H2 run has an LCOH.
    df["plant_moe_eur_per_t"] = [float("nan"), 120.0]
    df["lcoh_eur_per_kg"] = [4.1, float("nan")]
    flagged = compile_report.mark_best_in_country(df, PARENTS, "lco_output")

    report_path = tmp_path / "report_s.csv"
    compile_report.write_report(flagged, report_path, tmp_path / ".report_s_diag.csv")
    report = read_report(report_path)

    declared = [f for f in FIELD_ORDER if f != "best_in_country"]
    assert list(report.columns) == declared


def test_a_blank_means_undefined_and_a_zero_means_zero(tmp_path):
    """An amount the run has none of reads 0; a ratio with no denominator is blank."""
    df = _report([("s", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    flagged = compile_report.mark_best_in_country(df, PARENTS, "lco_output")

    report_path = tmp_path / "report_s.csv"
    compile_report.write_report(flagged, report_path, tmp_path / ".report_s_diag.csv")
    row = read_report(report_path).iloc[0]

    # No H2 chain on this route: the electrolyser it did not build is 0 GW, but
    # the cost of the hydrogen it never made is not a number.
    assert row["electrolyser_gw"] == 0.0
    assert row["cost_electrolyser_meur"] == 0.0
    assert row["plant_dri-h2_eur_per_t"] == 0.0
    assert pd.isna(row["lcoh_eur_per_mwh_lhv"])
    assert pd.isna(row["electrolyser_utilization"])
    assert pd.isna(row["cf_wind_offshore"])


def test_run_specific_fields_follow_the_declared_ones(tmp_path):
    """A multi-site run names a generator per candidate site; those are not declared."""
    df = _report([("s", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    df["solar-c00_gw_opt"] = 1.4
    flagged = compile_report.mark_best_in_country(df, PARENTS, "lco_output")

    report_path = tmp_path / "report_s.csv"
    compile_report.write_report(flagged, report_path, tmp_path / ".report_s_diag.csv")
    report = read_report(report_path)

    assert list(report.columns)[-1] == "solar-c00_gw_opt"
    assert report.at["s_1", "solar-c00_gw_opt"] == 1.4


def test_the_diagnostic_is_the_report_plus_the_flag(tmp_path):
    df = _report([
        ("s", "VIC1", "moe-eaf", "20250101", "20251231", 900.0),
        ("s", "NSW1", "moe-eaf", "20250101", "20251231", 820.0),
    ])
    flagged = compile_report.mark_best_in_country(df, PARENTS, "lco_output")
    report_path = tmp_path / "report_s.csv"
    diagnostic_path = tmp_path / ".report_s_diag.csv"
    compile_report.write_report(flagged, report_path, diagnostic_path)

    report = set(read_report(report_path).columns)
    diagnostic = set(read_report(diagnostic_path).columns)
    assert diagnostic - report == {"best_in_country"}


def test_a_report_reads_back_as_it_was_written(tmp_path):
    """read_report is the only thing that undoes the layout, so it has to be exact."""
    df = _report([
        ("s", "VIC1", "h2-dri-eaf", "20250101", "20251231", 606.78),
        ("s", "NSW1", "moe-eaf", "20240101", "20241231", 914.11),
    ])
    df["cf_solar"] = [0.19, 0.24]
    flagged = compile_report.mark_best_in_country(df, PARENTS, None)

    report_path = tmp_path / "report_s.csv"
    compile_report.write_report(flagged, report_path, tmp_path / ".report_s_diag.csv")
    report = read_report(report_path)

    assert list(report.index) == ["s_1", "s_2"]
    assert list(report["cf_solar"]) == [0.19, 0.24], "numbers come back as numbers"
    assert list(report["start_date"]) == ["20250101", "20240101"], "and dates as text"
    assert report["lco_output"].dtype == float
