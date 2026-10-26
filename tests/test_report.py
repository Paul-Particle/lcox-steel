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
    EMISSION_STEPS,
    FIELD_ORDER,
    IDENTITY_FIELDS,
    field_stem,
    read_report,
    write_report_file,
)

# The files pick_best_run copies up, each an input and an output of the same name.
RUN_FILES = ("report", "diagnostic", "network", "overlay", "base_assumptions", "merged_assumptions")


def _report(rows):
    """A report-shaped frame: one row per run, in the column order compile_report uses."""
    frame = pd.DataFrame(
        rows, columns=["scenario", "run", "origin", "route", "start_date", "end_date", "lcos_eur_per_t_steel"]
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
    """Run pick_best_run over one scenario whose runs cost `costs` ({run: lcos_eur_per_t_steel})."""
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
        params=SimpleNamespace(best_run_by="lcos_eur_per_t_steel"),
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
    # Only one run has a press to report, only the other an LCOH.
    df["cost_briquetting_plant_capex_eur_per_t_steel"] = [float("nan"), 120.0]
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
    assert row["cost_briquetting_plant_capex_eur_per_t_steel"] == 0.0
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
    df["emissions_energy_t_co2e_per_t_steel"] = 0.1574
    df["emissions_freight_kt_co2e_per_year"] = 215.0
    df["emissions_steel_transport_kg_co2e_per_t_steel"] = 215.0
    df["emissions_scope2_moe_t_co2e_per_t_steel"] = 0.1493
    report_path, diagnostic_path = _write(df, tmp_path)
    report = read_report(report_path).iloc[0]
    diagnostic = read_report(diagnostic_path).iloc[0]

    assert "emissions_steel_transport_kg_co2e_per_t_steel" not in report.index
    assert "emissions_freight_kt_co2e_per_year" not in report.index
    assert diagnostic["emissions_steel_transport_kg_co2e_per_t_steel"] == 215.0
    # The total counts neither file's freight, so the two agree on it.
    assert report["emissions_energy_t_co2e_per_t_steel"] == 0.1574
    assert diagnostic["emissions_energy_t_co2e_per_t_steel"] == 0.1574


def test_a_freight_field_no_run_produced_still_reaches_the_diagnostic(tmp_path):
    """The diagnostic declares every field, so which file one lands in is fixed."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    _, diagnostic_path = _write(df, tmp_path)

    diagnostic = read_report(diagnostic_path)
    assert diagnostic.at["s_VIC1", "emissions_steel_transport_kg_co2e_per_t_steel"] == 0.0
    assert diagnostic.at["s_VIC1", "emissions_freight_kt_co2e_per_year"] == 0.0


def test_emissions_in_tonnes_keep_four_decimals(tmp_path):
    """A clean route lands in thousandths of a tonne, which two decimals would erase."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    df["emissions_energy_mt_co2e_per_year"] = 0.004812
    df["emissions_energy_t_co2e_per_t_steel"] = 0.004812
    df["emissions_scope2_kg_co2e_per_mwh_el"] = 4.8123
    report_path, _ = _write(df, tmp_path)
    report = read_report(report_path).iloc[0]

    assert report["emissions_energy_mt_co2e_per_year"] == 0.0048
    assert report["emissions_energy_t_co2e_per_t_steel"] == 0.0048
    assert report["emissions_scope2_kg_co2e_per_mwh_el"] == 4.81


def test_the_run_timestamp_reads_back_as_written_with_its_offset(tmp_path):
    """It is text, not a number, so the offset survives the round trip."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    df["run_timestamp"] = "2026-10-26T14:03:12+01:00"
    report_path, _ = _write(df, tmp_path)

    assert read_report(report_path).iloc[0]["run_timestamp"] == "2026-10-26T14:03:12+01:00"


def test_the_model_run_id_reads_back_as_text(tmp_path):
    """Hex digits that happen to be all numbers must not come back as a number."""
    df = _report([("s", "VIC1", "VIC1", "moe-eaf", "20250101", "20251231", 914.0)])
    df["model_run_id"] = "401290128830-20261026T140312-00112233"
    report_path, _ = _write(df, tmp_path)

    assert read_report(report_path).iloc[0]["model_run_id"] == "401290128830-20261026T140312-00112233"


def test_a_report_reads_back_as_it_was_written(tmp_path):
    """read_report is the only thing that undoes the layout, so it has to be exact."""
    df = _report([("s", "VIC1", "VIC1", "h2-dri-eaf", "20250101", "20251231", 606.78)])
    df["cf_solar"] = [0.19]
    report_path, _ = _write(df, tmp_path)
    report = read_report(report_path)

    assert list(report.index) == ["s_VIC1"]
    assert list(report["cf_solar"]) == [0.19], "numbers come back as numbers"
    assert list(report["start_date"]) == ["20250101"], "and dates as text"
    assert report["lcos_eur_per_t_steel"].dtype == float


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

    assert after["lcos_eur_per_t_steel"] == pytest.approx(before["lcos_eur_per_t_steel"])
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


def test_the_ironmaking_capacity_is_the_shaft_at_full_load_for_a_year():
    """The asset's size, not its output: the built t/h over a whole year, which
    is at least what the shaft made and equals it only at full utilisation."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    hours = np.arange(48)
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)},
                      index=pd.date_range("2025-01-01", periods=48, freq="h"))
    n = build_network("h2-dri-eaf", assumptions, cf)
    n.optimize(solver_name="highs", solver_options={"output_flag": False})
    s = extract_summary(n, "s", {"origin": "X", "destination": "X", "route": "h2-dri-eaf",
                                 "start_date": "20250101", "end_date": "20250102"}, assumptions)

    shaft_t_per_h = n.links.at["dri-h2", "p_nom_opt"] * n.links.at["dri-h2", "efficiency"]
    assert s["ironmaking_capacity_mt_opt"] == pytest.approx(shaft_t_per_h * 8760 / 1e6)
    assert "dri_h2_t_per_h_opt" not in s, "one ironmaking row, not one per shaft"
    # The furnace's nameplate likewise: never below the steel it made.
    furnace_t_per_h = n.links.at["eaf", "p_nom_opt"] * n.links.at["eaf", "efficiency"]
    assert s["steelmaking_capacity_mt_opt"] == pytest.approx(furnace_t_per_h * 8760 / 1e6)
    assert "eaf_t_per_h_opt" not in s, "the furnace's size is the steelmaking row"
    assert s["steelmaking_capacity_mt_opt"] * s["steelmaking_utilization"] == pytest.approx(
        s["steel_produced_mt"])
    made_mt = _iron_produced_t(n) * 8760 / len(n.snapshots) / 1e6
    assert s["ironmaking_capacity_mt_opt"] * s["ironmaking_utilization"] == pytest.approx(made_mt)
    # The feedstock is quoted per tonne of iron out, so what was spent on it over
    # the iron made is the quote itself: the reason the report no longer repeats it.
    assert s["cost_ore_eur_per_t_iron"] == pytest.approx(assumptions["dri-h2"]["ore_eur_per_t"])


def test_iron_is_counted_at_the_shaft_and_the_briquettes_beside_it():
    """On a route that briquettes, iron is what the shaft made and HBI what the
    press made of it; every per-tonne-of-iron row is over the shaft's tonnes."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    hours = np.arange(48)
    index = pd.date_range("2025-01-01", periods=48, freq="h")
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)}, index=index)
    price = pd.Series(60 + 40 * np.sin(hours / 24 * 2 * np.pi), index=index)
    legs = {"sea": 8945}
    n = build_network("h2-dri-eaf-export", assumptions, cf, transport_legs=legs,
                      destination_price=price)
    n.optimize(solver_name="highs", solver_options={"output_flag": False})
    mix = pd.DataFrame({"price": price, "gas": 100.0, "wind_onshore": 100.0}, index=index)
    s = extract_summary(n, "s", {"origin": "X", "destination": "Y", "route": "h2-dri-eaf-export",
                                 "start_date": "20250101", "end_date": "20250102"},
                        assumptions, legs, None, mix)

    assert s["iron_produced_mt"] == pytest.approx(
        s["ironmaking_capacity_mt_opt"] * s["ironmaking_utilization"])
    assert s["hbi_produced_mt"] == pytest.approx(
        s["iron_produced_mt"] * assumptions["briquetting"]["yield_t_per_t"])
    # The press's nameplate times how hard it ran is what it made.
    assert s["briquetting_capacity_mt_opt"] * s["briquetting_utilization"] == pytest.approx(
        s["hbi_produced_mt"])
    assert s["cost_h2_production_eur_per_t_iron"] * s["iron_produced_mt"] == pytest.approx(
        s["cost_h2_production_meur"])
    # The destination furnace's power per MWh it drew, in three parts: together
    # they are its whole bill, which is the destination's electricity group.
    per_mwh = (s["destination_grid_price_eur_per_mwh"] + s["destination_grid_fee_eur_per_mwh"]
               + s["destination_grid_connection_eur_per_mwh_imported"])
    assert s["destination_grid_total_eur_per_mwh_imported"] == pytest.approx(per_mwh)
    assert s["destination_grid_import_gw_opt"] == pytest.approx(
        n.generators.at["destination_supply", "p_nom_opt"] / 1e3)
    assert "destination_supply_gw_opt" not in s
    assert per_mwh * s["el_eaf_gwh"] * 1e3 == pytest.approx(s["cost_destination_electricity_meur"] * 1e6)
    # The connection's two halves, as the destination's cost leaves have them.
    assert (s["destination_grid_connection_capex_eur_per_mwh_imported"]
            + s["destination_grid_capacity_fee_eur_per_mwh_imported"]) == pytest.approx(
        s["destination_grid_connection_eur_per_mwh_imported"])
    steel_t = s["steel_produced_mt"] * 1e6
    assert s["destination_grid_capacity_fee_eur_per_mwh_imported"] * s["el_eaf_gwh"] * 1e3 == pytest.approx(
        s["cost_destination_grid_capacity_fee_eur_per_t_steel"] * steel_t)

    # The briquette stockpile in hours of steel: what it holds over what the
    # furnace charges an hour, more than a tonne of briquettes per tonne of steel.
    steel_t_per_h = s["steel_produced_mt"] * 1e6 / 8760
    assert s["iron_or_hbi_store_duration_hours_steel"] == pytest.approx(
        s["iron_or_hbi_store_kt"] * 1e3
        / (steel_t_per_h * s["hbi_produced_mt"] / s["steel_produced_mt"]))
    # What went through the stockpile: its cycles are that over its size, its
    # share that over the briquettes made, which is where it sits on this route.
    assert s["iron_or_hbi_store_cycles_per_year"] == pytest.approx(
        s["iron_or_hbi_store_discharged_kt"] / s["iron_or_hbi_store_kt"])
    assert s["iron_or_hbi_store_share_pct"] == pytest.approx(
        s["iron_or_hbi_store_discharged_kt"] / 1e3 * 100 / s["hbi_produced_mt"])
    # The steel store's rows, built the same way over the steel made.
    if s["steel_store_kt"] > 0:
        assert s["steel_store_cycles_per_year"] == pytest.approx(
            s["steel_store_discharged_kt"] / s["steel_store_kt"])
    assert s["steel_store_share_pct"] == pytest.approx(
        s["steel_store_discharged_kt"] / 1e3 * 100 / s["steel_produced_mt"])

def test_the_origin_s_grid_per_mwh_has_the_cost_leaves_four_parts():
    """Price, volumetric fee, connection capital and capacity fee per MWh imported,
    each the matching cost leaf over the MWh, and the connection their subtotal."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    hours = np.arange(48)
    index = pd.date_range("2025-01-01", periods=48, freq="h")
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)}, index=index)
    price = pd.Series(40 + 30 * np.sin(hours / 24 * 2 * np.pi), index=index)
    n = build_network("h2-dri-eaf", assumptions, cf, price_series=price)
    n.optimize(solver_name="highs", solver_options={"output_flag": False})
    s = extract_summary(n, "s", {"origin": "X", "destination": "X", "route": "h2-dri-eaf",
                                 "start_date": "20250101", "end_date": "20250102"}, assumptions)
    imported_mwh = s["grid_imported_gwh"] * 1e3
    steel_t = s["steel_produced_mt"] * 1e6
    assert imported_mwh > 0, "the test needs a run that bought power"

    assert (s["grid_connection_capex_eur_per_mwh_imported"]
            + s["grid_capacity_fee_eur_per_mwh_imported"]) == pytest.approx(
        s["grid_connection_eur_per_mwh_imported"])
    for per_mwh, leaf in (("grid_connection_capex", "grid_connection_capex"),
                          ("grid_capacity_fee", "grid_capacity_fee")):
        assert s[f"{per_mwh}_eur_per_mwh_imported"] * imported_mwh == pytest.approx(
            s[f"cost_{leaf}_eur_per_t_steel"] * steel_t)

def test_the_gas_side_reports_its_peak_draw_not_a_free_rating():
    """The gas supply costs nothing to size, so its rating is not reported; its
    busiest hour is, and that is at least the year's average draw."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    hours = np.arange(48)
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)},
                      index=pd.date_range("2025-01-01", periods=48, freq="h"))
    n = build_network("ng-dri-eaf", assumptions, cf)
    n.optimize(solver_name="highs", solver_options={"output_flag": False})
    s = extract_summary(n, "s", {"origin": "X", "destination": "X", "route": "ng-dri-eaf",
                                 "start_date": "20250101", "end_date": "20250102"}, assumptions)

    assert "gas_supply_gw_opt" not in s
    assert s["gas_peak_draw_gw"] == pytest.approx(float(n.generators_t.p["gas_supply"].max()) / 1e3)
    assert s["gas_peak_draw_gw"] * 8760 >= s["gas_consumed_gwh_lhv"] - 1e-6
    # A gas shaft emits both ways, so its split has two real parts to add up.
    assert s["emissions_scope1_gas_t_co2e_per_t_steel"] > 0
    assert (s["emissions_scope1_gas_t_co2e_per_t_steel"]
            + s["emissions_scope2_electricity_t_co2e_per_t_steel"]
            == pytest.approx(s["emissions_energy_t_co2e_per_t_steel"]))
    # The steps are the electricity's alone, and with the gas's share beside
    # them the shares make the whole.
    steps = [f"emissions_scope2_{field_stem(step)}_t_co2e_per_t_steel" for step in EMISSION_STEPS]
    assert sum(s[field] for field in steps) == pytest.approx(
        s["emissions_scope2_electricity_t_co2e_per_t_steel"])
    shares = [f"emissions_scope2_{field_stem(step)}_pct" for step in EMISSION_STEPS]
    assert sum(s.get(field, 0.0) for field in shares) + s["emissions_scope1_gas_pct"] == pytest.approx(100.0)
    # Every shaft and cell is one reduction step, and ironmaking is it plus the
    # electrolyser and its share of the losses: no less than its own steps, no
    # more than the whole of scope 2.
    assert "emissions_scope2_dri_ng_t_co2e_per_t_steel" not in s
    iron_steps = (s["emissions_scope2_reduction_t_co2e_per_t_steel"]
                  + s["emissions_scope2_electrolyser_t_co2e_per_t_steel"])
    assert iron_steps - 1e-12 <= s["emissions_scope2_ironmaking_t_co2e_per_t_steel"]
    assert (s["emissions_scope2_ironmaking_t_co2e_per_t_steel"]
            <= s["emissions_scope2_electricity_t_co2e_per_t_steel"] + 1e-12)


def test_ironmaking_s_emissions_take_their_share_of_the_losses_by_the_mwh_drawn():
    """Reduction and electrolyser, plus the losses in proportion to the MWh they
    drew among everything drawn at the origin, as their electricity bill is cut."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    # A dirty solar and wind, made up, so every step and every loss emits.
    factors = assumptions["emissions"][assumptions["emissions"]["basis"]]["electricity_t_co2e_per_mwh"]
    factors["solar"] = factors["wind_onshore"] = 0.1
    hours = np.arange(48)
    cf = pd.DataFrame({"solar": np.clip(np.sin((hours % 24 - 6) / 12 * np.pi), 0, None),
                       "wind-onshore": 0.35 + 0.3 * np.sin(hours / 7.0)},
                      index=pd.date_range("2025-01-01", periods=48, freq="h"))
    n = build_network("h2-dri-eaf", assumptions, cf)
    n.optimize(solver_name="highs", solver_options={"output_flag": False})
    s = extract_summary(n, "s", {"origin": "X", "destination": "X", "route": "h2-dri-eaf",
                                 "start_date": "20250101", "end_date": "20250102"}, assumptions)

    losses = (s["emissions_scope2_battery_losses_t_co2e_per_t_steel"]
              + s["emissions_scope2_transmission_losses_t_co2e_per_t_steel"])
    assert losses > 0, "the test needs a run that lost electricity"
    iron_gwh = s["el_dri_h2_gwh"] + s["el_h2_preheat_gwh"] + s["el_electrolyser_gwh"]
    share = iron_gwh / s["el_users_gwh"]
    assert s["emissions_scope2_ironmaking_t_co2e_per_t_steel"] == pytest.approx(
        s["emissions_scope2_reduction_t_co2e_per_t_steel"]
        + s["emissions_scope2_electrolyser_t_co2e_per_t_steel"]
        + losses * share)
    # The preheat is part of the reduction step, not a step of its own.
    assert "emissions_scope2_h2_preheat_t_co2e_per_t_steel" not in s
    # The hydrogen takes its share of the same losses by the same rule.
    h2_t = s["emissions_scope2_kg_co2e_per_kg_h2"] * s["h2_produced_kt"] / 1e3
    steel_t = s["steel_produced_mt"] * 1e6
    assert h2_t * 1e6 / steel_t == pytest.approx(
        s["emissions_scope2_electrolyser_t_co2e_per_t_steel"]
        + losses * s["el_electrolyser_gwh"] / s["el_users_gwh"])
    # A step's intensity is its CO2e over the MWh it drew, shaft and preheat
    # together for the reduction.
    assert "emissions_dri_h2_kg_co2e_per_mwh_el" not in s
    reduction_mwh_per_t = s["el_reduction_mwh_per_t_steel"]
    assert s["emissions_scope2_reduction_kg_co2e_per_mwh_el"] == pytest.approx(
        s["emissions_scope2_reduction_t_co2e_per_t_steel"] * 1e3 / reduction_mwh_per_t)
    assert s["emissions_scope2_eaf_kg_co2e_per_mwh_el"] == pytest.approx(
        s["emissions_scope2_eaf_t_co2e_per_t_steel"] * 1e3 / s["el_eaf_mwh_per_t_steel"])
