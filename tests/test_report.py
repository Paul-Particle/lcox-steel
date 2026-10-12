"""Unit tests for what a report says and which run stands for its scenario.

A report is one run: compile_report writes it, and the hidden diagnostic beside
it, through `write_report_file`. A scenario usually has one run; an area on the
grid that trades through zones (Australia through its NEM regions) is run once
per zone, and `pick_best_run` copies the cheapest zone's files up to stand for
the scenario. These pin that behaviour on synthetic frames — no cutouts, no
solves, so they run anywhere.
"""

from types import SimpleNamespace

import pandas as pd

from scripts.viz import _run_display, pick_best_run
from common._report_schema import (
    DIAGNOSTIC_FIELDS,
    FIELD_ORDER,
    IDENTITY_FIELDS,
    read_report,
    write_report_file,
)

# The files pick_best_run copies up, each an input and an output of the same name.
RUN_FILES = ("report", "diagnostic", "network", "overlay", "base_assumptions", "merged_assumptions")


def _report(rows):
    """A report-shaped frame: one row per run, in the column order compile_report uses."""
    frame = pd.DataFrame(
        rows, columns=["scenario", "run", "area", "route", "start_date", "end_date", "lco_output"]
    )
    # Every row carries the fingerprint of the inputs behind it, so a frame
    # without one is not a report frame.
    frame["inputs_hash"] = "772c9cac81fe"
    return frame


def _write(df, tmp_path):
    """The report and the diagnostic, written as compile_report writes them."""
    report_path = tmp_path / "report_s.csv"
    diagnostic_path = tmp_path / ".report_s_diag.csv"
    write_report_file(df, diagnostic_path)
    write_report_file(df, report_path, hold_back=DIAGNOSTIC_FIELDS)
    return report_path, diagnostic_path


def test_a_run_is_a_column_and_a_field_is_a_row(tmp_path):
    """The file's shape: the run across the top, named after its folder, the identity fields leading."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    report_path, diagnostic_path = _write(df, tmp_path)
    stored = pd.read_csv(report_path)
    report_identity = [field for field in IDENTITY_FIELDS if field not in DIAGNOSTIC_FIELDS]

    assert list(stored.columns) == ["field", "s_VIC1"]
    assert list(stored["field"])[:len(report_identity)] == report_identity
    assert stored.set_index("field").at["route", "s_VIC1"] == "moe-eaf"
    assert "run" not in set(stored["field"]), "the column's name already says which run it is"
    assert read_report(diagnostic_path).at["s_VIC1", "run"] == "VIC1"


def _pick(tmp_path, monkeypatch, costs):
    """Run pick_best_run over one scenario whose runs cost `costs` ({run: lco_output})."""
    inputs = {name: [] for name in RUN_FILES}
    for run, cost in costs.items():
        run_dir = tmp_path / ".runs" / f"s_{run}"
        run_dir.mkdir(parents=True)
        report_path, diagnostic_path = _write(
            _report([("s", run, run, "moe-eaf", "20250101", "20251231", cost)]), run_dir
        )
        inputs["report"].append(report_path)
        inputs["diagnostic"].append(diagnostic_path)
        for name in RUN_FILES[2:]:
            (run_dir / name).write_text(f"{name} of {run}")
            inputs[name].append(run_dir / name)
    outputs = {name: tmp_path / name for name in RUN_FILES}
    monkeypatch.setattr(pick_best_run, "snakemake", SimpleNamespace(
        input=SimpleNamespace(**inputs),
        output=SimpleNamespace(**outputs, runs=tmp_path / ".runs" / "runs.csv"),
        params=SimpleNamespace(best_run_by="lco_output"),
        wildcards=SimpleNamespace(scenario="s"),
    ))
    pick_best_run.main()
    picked = pd.read_csv(tmp_path / ".runs" / "runs.csv")
    return picked


def test_the_cheapest_zone_stands_for_its_scenario(tmp_path, monkeypatch):
    picked = _pick(tmp_path, monkeypatch, {"VIC1": 900.0, "NSW1": 820.0, "QLD1": 875.0})

    assert list(picked.loc[picked["picked"], "run"]) == ["NSW1"]
    assert len(picked) == 3, "every run is listed, the losers included"
    assert (tmp_path / "network").read_text() == "network of NSW1"
    assert read_report(tmp_path / "report").iloc[0]["area"] == "NSW1"


def test_a_single_run_stands_for_its_scenario(tmp_path, monkeypatch):
    picked = _pick(tmp_path, monkeypatch, {"DEU": 1010.0})

    assert list(picked["picked"]) == [True]
    assert (tmp_path / "merged_assumptions").read_text() == "merged_assumptions of DEU"


def test_a_run_without_the_metric_is_not_picked_over_one_with_it(tmp_path, monkeypatch):
    picked = _pick(tmp_path, monkeypatch, {"VIC1": float("nan"), "NSW1": 820.0})

    assert list(picked.loc[picked["picked"], "run"]) == ["NSW1"]


def test_one_window_is_labelled_by_area_and_route_alone():
    """A chart of a single window needs no date to tell its bars apart."""
    df = _report([
        ("s", "NSW1", "NSW1", "moe-eaf", "20250101", "20251231", 900.0),
        ("t", "DEU", "DEU", "h2-dri-eaf", "20250101", "20251231", 1010.0),
    ])
    assert list(_run_display.run_labels(df)) == ["NSW1 moe-eaf", "DEU h2-dri-eaf"]


def test_two_windows_of_one_area_are_told_apart_by_year():
    """The frame decides, not the run: two single-year windows are each one year,
    so a label read off the run alone would give both bars the same name."""
    df = _report([
        ("s", "DEU", "DEU", "h2-dri-eaf", "20240101", "20241231", 900.0),
        ("t", "DEU", "DEU", "h2-dri-eaf", "20250101", "20251231", 820.0),
    ])
    labels = _run_display.run_labels(df)
    assert list(labels) == ["DEU h2-dri-eaf 2024", "DEU h2-dri-eaf 2025"]
    assert labels.is_unique


def test_a_window_that_crosses_new_year_names_both_years():
    df = _report([
        ("s", "DEU", "DEU", "h2-dri-eaf", "20230101", "20241231", 900.0),
        ("t", "DEU", "DEU", "h2-dri-eaf", "20250101", "20251231", 820.0),
    ])
    assert list(_run_display.run_labels(df)) == [
        "DEU h2-dri-eaf 2023-2024", "DEU h2-dri-eaf 2025",
    ]


def test_every_route_writes_the_same_fields(tmp_path):
    """The point of the schema: what a run reports must not depend on its route."""
    df = _report([
        ("s", "VIC1", "VIC1", "h2-only", "20250101", "20251231", 4.1),
        ("t", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0),
    ])
    # Only the MOE run has a MOE cell to report; only the H2 run has an LCOH.
    df["plant_moe_eur_per_t"] = [float("nan"), 120.0]
    df["lcoh_eur_per_kg"] = [4.1, float("nan")]
    report_path, _ = _write(df, tmp_path)
    report = read_report(report_path)

    declared = [f for f in FIELD_ORDER if f not in DIAGNOSTIC_FIELDS]
    assert list(report.columns) == declared


def test_a_blank_means_undefined_and_a_zero_means_zero(tmp_path):
    """An amount the run has none of reads 0; a ratio with no denominator is blank."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    report_path, _ = _write(df, tmp_path)
    row = read_report(report_path).iloc[0]

    # No H2 chain on this route: the electrolyser it did not build is 0 GW, but
    # the cost of the hydrogen it never made is not a number.
    assert row["electrolyser_gw"] == 0.0
    assert row["cost_electrolyser_meur"] == 0.0
    assert row["plant_dri_h2_eur_per_t"] == 0.0
    assert pd.isna(row["lcoh_eur_per_mwh_lhv"])
    assert pd.isna(row["electrolyser_utilization"])
    assert pd.isna(row["cf_wind_offshore"])


def test_run_specific_fields_follow_the_declared_ones(tmp_path):
    """A multi-site run names a generator per candidate site; those are not declared."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    df["solar_c00_gw_opt"] = 1.4
    report_path, _ = _write(df, tmp_path)
    report = read_report(report_path)

    assert list(report.columns)[-1] == "solar_c00_gw_opt"
    assert report.at["s_VIC1", "solar_c00_gw_opt"] == 1.4


def test_the_diagnostic_is_the_report_plus_the_held_back_fields(tmp_path):
    df = _report([("s", "NSW1", "NSW1", "moe-eaf", "20250101", "20251231", 820.0)])
    report_path, diagnostic_path = _write(df, tmp_path)

    report = set(read_report(report_path).columns)
    diagnostic = set(read_report(diagnostic_path).columns)
    assert diagnostic - report == set(DIAGNOSTIC_FIELDS)


def test_the_freight_a_run_reports_is_held_back_from_the_report(tmp_path):
    """A route that ships writes its freight, and only the diagnostic carries it.

    No total on either side counts it, so the two files agree on every field
    they share and differ only in that the report has no freight at all.
    """
    df = _report([("s", "VIC1", "VIC1", "moe-eaf-export", "20250101", "20251231", 900.0)])
    df["emissions_kg_co2e_per_t_steel"] = 157.4
    df["emissions_freight_kt_co2e_per_year"] = 215.0
    df["emissions_steel_transport_kg_co2e_per_t_steel"] = 215.0
    df["emissions_moe_kg_co2e_per_t_steel"] = 149.3
    report_path, diagnostic_path = _write(df, tmp_path)
    report = read_report(report_path).iloc[0]
    diagnostic = read_report(diagnostic_path).iloc[0]

    assert "emissions_steel_transport_kg_co2e_per_t_steel" not in report.index
    assert "emissions_freight_kt_co2e_per_year" not in report.index
    assert diagnostic["emissions_steel_transport_kg_co2e_per_t_steel"] == 215.0
    # The total counts neither file's freight, so the two agree on it.
    assert report["emissions_kg_co2e_per_t_steel"] == 157.4
    assert diagnostic["emissions_kg_co2e_per_t_steel"] == 157.4


def test_a_freight_field_no_run_produced_still_reaches_the_diagnostic(tmp_path):
    """The diagnostic declares every field, so which file one lands in is fixed."""
    df = _report([("s", "VIC1", "VIC1", "h2-only", "20250101", "20251231", 4.1)])
    _, diagnostic_path = _write(df, tmp_path)

    diagnostic = read_report(diagnostic_path)
    assert diagnostic.at["s_VIC1", "emissions_steel_transport_kg_co2e_per_t_steel"] == 0.0
    assert diagnostic.at["s_VIC1", "emissions_freight_kt_co2e_per_year"] == 0.0


def test_a_report_reads_back_as_it_was_written(tmp_path):
    """read_report is the only thing that undoes the layout, so it has to be exact."""
    df = _report([("s", "VIC1", "VIC1", "h2-dri-eaf", "20250101", "20251231", 606.78)])
    df["cf_solar"] = [0.19]
    report_path, _ = _write(df, tmp_path)
    report = read_report(report_path)

    assert list(report.index) == ["s_VIC1"]
    assert list(report["cf_solar"]) == [0.19], "numbers come back as numbers"
    assert list(report["start_date"]) == ["20250101"], "and dates as text"
    assert report["lco_output"].dtype == float
