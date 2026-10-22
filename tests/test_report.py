"""Unit tests for what a report says and which run stands for its scenario.

A report is one run: compile_report writes it, and the hidden diagnostic beside
it, through `write_report_file`. A scenario usually has one run; an area on the
grid that trades through zones (Australia through its NEM regions) is run once
per zone, and `pick_best_run` copies the cheapest zone's files up to stand for
the scenario. These pin that behaviour on synthetic frames — no cutouts, no
solves, so they run anywhere.
"""

import copy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytest
import yaml

from conftest import REPO_ROOT
from scripts.solve.build_network import build_network
from scripts.viz import _run_display, pick_best_run
from scripts.viz.compile_report import (
    _battery_discharged_mwh,
    _h2_buffer_discharged_mwh,
    _iron_produced_t,
    _res_energy,
    extract_summary,
    input_variants,
)
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
        rows, columns=["scenario", "run", "origin", "route", "start_date", "end_date", "lcos_eur_per_t"]
    )
    # Iron and steel made in one place, as compile_report writes a domestic run.
    frame.insert(3, "destination", frame["origin"])
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
    assert stored.set_index("field").at["destination", "s_VIC1"] == "VIC1"
    assert "run" not in set(stored["field"]), "the column's name already says which run it is"
    assert read_report(diagnostic_path).at["s_VIC1", "run"] == "VIC1"


def _pick(tmp_path, monkeypatch, costs):
    """Run pick_best_run over one scenario whose runs cost `costs` ({run: lcos_eur_per_t})."""
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
        params=SimpleNamespace(best_run_by="lcos_eur_per_t"),
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
    assert read_report(tmp_path / "report").iloc[0]["origin"] == "NSW1"


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
        ("s", "VIC1", "VIC1", "h2-dri-eaf", "20250101", "20251231", 606.78),
        ("t", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0),
    ])
    # Only the MOE run has a MOE cell to report; only the H2 run has an LCOH.
    df["plant_moe_eur_per_t_steel"] = [float("nan"), 120.0]
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
    assert row["plant_dri_h2_eur_per_t_steel"] == 0.0
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
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
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
    assert report["lcos_eur_per_t"].dtype == float


def _rows(*series):
    """A run's rows in its scenario table, as (tech, variant, area) triples."""
    return pd.DataFrame(series, columns=["tech", "variant", "area"])


def test_an_export_reports_its_own_grid_not_its_furnaces():
    """The destination's grid row prices the furnace abroad; the origin's is the plant's."""
    rows = _rows(("solar", "bestsite-p95", "AUS"), ("grid", "dayahead", "AUS"),
                 ("grid", "emissions", "DEU"))
    assert input_variants(rows, "AUS") == {"solar_variant": "bestsite-p95",
                                           "grid_variant_origin": "dayahead"}


def test_an_islanded_export_reports_no_grid():
    """A blank grid_variant_origin is what tells the dashboard the plant was islanded."""
    rows = _rows(("solar", "bestsite-p95", "AUS"), ("grid", "emissions", "DEU"))
    assert "grid_variant_origin" not in input_variants(rows, "AUS")


def test_a_domestic_run_reports_its_grid():
    rows = _rows(("wind-onshore", "area-average", "DEU"), ("grid", "emissions", "DEU"))
    assert input_variants(rows, "DEU") == {"wind_onshore_variant": "area-average",
                                           "grid_variant_origin": "emissions"}


def _res_network():
    """Two solar orientations, one wind farm and a grid connection, solved by hand.

    A day of hours at a flat half-capacity availability, so what each tech could
    have made is easy to state: 24 h x 0.5 x its capacity.
    """
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2025-01-01", periods=24, freq="h"))
    n.add("Bus", "electricity")
    for name in ("solar_az0", "solar_az30", "wind-onshore", "grid_import"):
        n.add("Generator", name, bus="electricity", p_nom_extendable=True,
              p_max_pu=0.5 if name != "grid_import" else 1.0)
    n.generators["p_nom_opt"] = [100.0, 60.0, 200.0, 80.0]
    # Solar runs flat out; wind is held back to 40 of its 100 MW available.
    n.generators_t.p = pd.DataFrame(
        {"solar_az0": 50.0, "solar_az30": 30.0, "wind-onshore": 40.0, "grid_import": 80.0},
        index=n.snapshots,
    )
    return n


def test_a_tech_s_generation_sums_its_orientations_and_leaves_the_grid_out():
    """Each orientation is a generator of its own; the tech's row is all of them."""
    fields = _res_energy(_res_network(), hours_per_year=24)
    assert fields["solar_generated_gwh"] == pytest.approx((50 + 30) * 24 / 1e3)
    assert fields["wind_onshore_generated_gwh"] == pytest.approx(40 * 24 / 1e3)
    assert fields["wind_offshore_generated_gwh"] == 0.0


def test_curtailment_is_what_was_available_but_not_delivered():
    fields = _res_energy(_res_network(), hours_per_year=24)
    assert fields["solar_curtailed_gwh"] == pytest.approx(0.0)
    assert fields["wind_onshore_curtailed_gwh"] == pytest.approx((100 - 40) * 24 / 1e3)


def test_a_short_run_reports_a_year():
    """A day's run is scaled up to a year, like every other energy in the report."""
    fields = _res_energy(_res_network(), hours_per_year=8760)
    assert fields["wind_onshore_generated_gwh"] == pytest.approx(40 * 8760 / 1e3)


def test_a_battery_is_credited_with_what_it_gave_back_not_what_it_drew():
    """The discharge loss stays out: 10 MW drawn off the store at 90 % is 9 MW delivered."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2025-01-01", periods=24, freq="h"))
    n.add("Bus", "electricity")
    n.add("Bus", "electricity_battery")
    n.add("Link", "battery_discharger", bus0="electricity_battery", bus1="electricity",
          efficiency=0.9)
    n.links_t.p0 = pd.DataFrame({"battery_discharger": 10.0}, index=n.snapshots)

    assert _battery_discharged_mwh(n, hours_per_year=24) == pytest.approx(9.0 * 24)
    assert _battery_discharged_mwh(n, hours_per_year=8760) == pytest.approx(9.0 * 8760)


def test_the_h2_shaft_reports_its_preheat_where_the_blend_shaft_does():
    """The H2 feed's preheat lands in the one preheat row on either shaft, the shaft
    row keeps only the shaft's own power, and splitting the draw moves no number
    the solve decides."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    hours = np.arange(48)
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)},
                      index=pd.date_range("2025-01-01", periods=48, freq="h"))
    run = {"origin": "X", "destination": "X", "route": "h2-dri-eaf",
           "start_date": "20250101", "end_date": "20250102"}

    def solved(a):
        n = build_network("h2-dri-eaf", a, cf)
        n.optimize(solver_name="highs", solver_options={"output_flag": False})
        return extract_summary(n, "s", run, assumptions)

    # The shaft as it was before the split: one draw carrying both.
    lumped = copy.deepcopy(assumptions)
    shaft = lumped["dri-h2"]
    shaft["el_mwh_per_t"] += shaft["h2_preheat_el_mwh_per_t"]
    shaft["h2_preheat_el_mwh_per_t"] = 0.0
    before, after = solved(lumped), solved(assumptions)

    assert after["lcos_eur_per_t"] == pytest.approx(before["lcos_eur_per_t"])
    assert after["el_users_gwh"] == pytest.approx(before["el_users_gwh"])
    assert after["el_dri_h2_gwh"] + after["el_h2_preheat_gwh"] == pytest.approx(before["el_dri_h2_gwh"])
    preheat_share = (assumptions["dri-h2"]["h2_preheat_el_mwh_per_t"]
                     / (assumptions["dri-h2"]["el_mwh_per_t"]
                        + assumptions["dri-h2"]["h2_preheat_el_mwh_per_t"]))
    assert after["el_h2_preheat_gwh"] == pytest.approx(preheat_share * before["el_dri_h2_gwh"])



def test_the_h2_buffer_counts_what_it_gave_back_not_what_it_took_in():
    """Only the hours the store fed its bus are discharge: filling is not."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2025-01-01", periods=24, freq="h"))
    n.add("Bus", "hydrogen")
    n.add("Store", "h2_buffer", bus="hydrogen", e_cyclic=True)
    # Twelve hours filling at 30 MW, twelve emptying at 30 MW.
    n.stores_t.p = pd.DataFrame({"h2_buffer": [-30.0] * 12 + [30.0] * 12}, index=n.snapshots)

    assert _h2_buffer_discharged_mwh(n, hours_per_year=24) == pytest.approx(30.0 * 12)
    assert _h2_buffer_discharged_mwh(n, hours_per_year=8760) == pytest.approx(30.0 * 12 * 365)


def test_a_zero_is_written_without_a_sign(tmp_path):
    """Solver dust and a flipped zero both round to nothing, and read as 0, not -0."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    df["h2_produced_kg_per_t_steel"] = -0.0
    df["h2_buffer_gwh"] = -1e-12
    report_path, _ = _write(df, tmp_path)
    stored = pd.read_csv(report_path).set_index("field")["s_VIC1"]

    assert stored["h2_produced_kg_per_t_steel"] == "0.0"
    assert stored["h2_buffer_gwh"] == "0.0"


def test_the_hydrogen_s_yearly_cost_is_its_lcoh_times_what_was_made():
    """The M€ set closes on itself: parts to total, and each part to its per-kg row."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    hours = np.arange(48)
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)},
                      index=pd.date_range("2025-01-01", periods=48, freq="h"))
    n = build_network("h2-dri-eaf", assumptions, cf)
    n.optimize(solver_name="highs", solver_options={"output_flag": False})
    s = extract_summary(n, "s", {"origin": "X", "destination": "X", "route": "h2-dri-eaf",
                                 "start_date": "20250101", "end_date": "20250102"}, assumptions)
    made = s["h2_produced_kt"]

    assert s["cost_h2_production_meur"] == pytest.approx(
        s["cost_electrolyser_meur"] + s["cost_h2_buffer_meur"]
        + s["cost_electricity_electrolyser_meur"])
    assert s["cost_h2_production_meur"] == pytest.approx(s["lcoh_eur_per_kg"] * made)
    assert s["cost_electricity_electrolyser_meur"] == pytest.approx(
        s["lcoh_electricity_eur_per_kg"] * made)
    assert s["cost_electrolyser_meur"] == pytest.approx(s["lcoh_electrolyser_eur_per_kg"] * made)
    assert s["cost_electrolyser_meur"] == pytest.approx(
        sum(s[f"cost_electrolyser_{part}_meur"] for part in ("capex", "fom", "water")))
    for part in ("capex", "fom", "water"):
        assert s[f"cost_electrolyser_{part}_meur"] == pytest.approx(
            s[f"lcoh_electrolyser_{part}_eur_per_kg"] * made)
    # Per tonne of steel: the same total, over the steel, and on a route that
    # melts at home its three parts are the rows the steel's breakdown carries.
    assert s["cost_h2_production_eur_per_t_steel"] == pytest.approx(
        s["cost_h2_production_meur"] / s["steel_produced_mt"])
    assert s["cost_h2_production_eur_per_t_steel"] == pytest.approx(
        s["cost_electrolyser_eur_per_t_steel"] + s["cost_h2_buffer_eur_per_t_steel"]
        + s["cost_electricity_electrolyser_eur_per_t_steel"])
    # And per tonne of iron: all of the hydrogen over the iron it reduced.
    iron_t = _iron_produced_t(n) * 8760.0 / len(n.snapshots)
    assert s["cost_h2_production_eur_per_t_iron"] == pytest.approx(
        s["cost_h2_production_meur"] * 1e6 / iron_t)
    assert s["cost_h2_production_eur_per_t_iron"] == pytest.approx(
        s["cost_electrolyser_eur_per_t_iron"] + s["cost_h2_buffer_eur_per_t_iron"]
        + s["cost_electricity_electrolyser_eur_per_t_iron"])
    assert s["cost_electrolyser_eur_per_t_iron"] == pytest.approx(
        sum(s[f"cost_electrolyser_{part}_eur_per_t_iron"] for part in ("capex", "fom", "water")))
