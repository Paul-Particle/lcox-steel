"""Unit tests for turning the input workbook's values into model keys.

The workbook rows are a stand-in for the `assumptions_inputs` tab of
data/assumptions/Assumptions_yaml_inputs.xlsx, indexed by its YAML key column.
"""

import json
import re
import sys

import pandas as pd
import pytest

from common import generate_inputs
from common._workbook_translation import UNREAD_ROWS
from common.generate_inputs import (
    REGIONAL_SOURCES,
    SENSITIVITY_ROWS,
    build_scenarios,
    flatten,
    read_overrides,
    read_workbook_value,
    set_nested,
)

WORKBOOK_ROWS = pd.DataFrame(
    [("plant.capex_eur_base", 1000, "EUR/t"),
     ("plant.capex_eur_high", 2000, "EUR/t"),
     ("plant.capex_usd_base", 1130, "USD/t"),
     ("plant.empty_mwh", None, "MWh/t"),
     ("plant.lifetime", 20, "years")],
    columns=["key", "value", "unit"],
).set_index("key")


read_with = lambda path, spec: read_workbook_value(WORKBOOK_ROWS, path, spec)


def test_a_key_takes_its_row():
    assert read_with("plant.lifetime_years", "plant.lifetime") == 20
    assert read_with("plant.capex_eur_per_t", "plant.capex_eur_high") == 2000


def test_a_missing_row_or_value_stops_the_build():
    with pytest.raises(KeyError, match="plant.capex_eur_per_t"):
        read_with("plant.capex_eur_per_t", "plant.capex_eur_low")
    with pytest.raises(ValueError, match="no value"):
        read_with("plant.heat_mwh_per_t", "plant.empty_mwh")


def test_a_currency_the_key_does_not_name_stops_the_build():
    """The model is in EUR: a USD row, a EUR row under a key without `eur`, or the reverse is a wrong row."""
    with pytest.raises(ValueError, match="currency"):
        read_with("plant.capex_eur_per_t", "plant.capex_usd_base")
    with pytest.raises(ValueError, match="currency"):
        read_with("plant.lifetime_years", "plant.capex_eur_base")
    with pytest.raises(ValueError, match="currency"):
        read_with("plant.capex_eur_per_t", "plant.lifetime")


def test_overrides_survive_excel_curly_quotes():
    assert read_overrides('{“finance.origin_wacc”:“+30%”}') == {"finance.origin_wacc": "+30%"}
    assert read_overrides(float("nan")) == {}


IMPORT_ROW = {
    "scenario_group": "import", "iron_location": "AUS", "steel_location": "ESP", "route": "moe-eaf",
    "power_supply": "islanded_res", "res_configuration": "single", "res_technology": "solar", "res_anchor": None,
    "site_variant": "bestsite-p95", "start_date": 20250101, "end_date": 20251231, "sensitivity_overrides_JSON": "{}",
}
NAMED_INPUTS = ({path: {"noneu": 1.0, "aus": 1.0, "global": 1.0, "bra": 3.0, "eu": 2.0} for path in REGIONAL_SOURCES}
                | {path: {"low": 5.0, "high": 7.0} for path in SENSITIVITY_ROWS})
AREAS = {"AUS": {"zones": ["NSW1", "QLD1"]}, "ESP": {"market": "entsoe"}}


def build_from_rows(changes_by_scenario: dict) -> tuple[pd.DataFrame, dict]:
    """The scenario table and flattened overlays of AUS->ESP import rows with the given changes, on a base of zeros."""
    scenario_rows = pd.DataFrame([{**IMPORT_ROW, "scenario_id": scenario, **changes}
                                  for scenario, changes in changes_by_scenario.items()])
    assumptions = {}
    for path in NAMED_INPUTS:
        set_nested(assumptions, path, 0.0)
    scenarios, overlays = build_scenarios(scenario_rows, NAMED_INPUTS, assumptions, AREAS)
    flat_overlays = {scenario: flatten(overlay) for scenario, overlay in overlays.items()}
    return scenarios, flat_overlays


def test_an_import_takes_iron_side_values_from_its_iron_location_and_the_furnace_from_its_steel_location():
    """The overrides apply on top, each to the model key it names: changing the region's value,
    switching to another region's, or setting a number."""
    _, overlays = build_from_rows({
        "wacc-up": {"sensitivity_overrides_JSON": json.dumps({"finance.origin_wacc": "+30%",
                                                              "finance.destination_wacc": "+30%"})},
        "rows-overridden": {"sensitivity_overrides_JSON": json.dumps({
            "battery.capex_per_mw_eur": "eu", "moe.opex_per_t_per_year_eur": "+10%"})},
        "destination-wacc-set": {"sensitivity_overrides_JSON": json.dumps({"finance.destination_wacc": 0.09})},
    })

    wacc_up, rows_overridden, destination_wacc_set = (
        overlays["wacc-up"], overlays["rows-overridden"], overlays["destination-wacc-set"])
    assert wacc_up["finance.origin_wacc"] == pytest.approx(1.3)
    assert wacc_up["finance.destination_wacc"] == pytest.approx(2.6)
    assert wacc_up["moe.opex_per_t_per_year_eur"] == 1.0
    assert wacc_up["battery.capex_per_mw_eur"] == 1.0
    assert wacc_up["eaf.opex_per_t_per_year_eur"] == 2.0
    assert rows_overridden["moe.opex_per_t_per_year_eur"] == pytest.approx(1.1)
    assert rows_overridden["battery.capex_per_mw_eur"] == 2.0
    assert rows_overridden["battery.capex_per_mwh_eur"] == 1.0
    assert destination_wacc_set["finance.origin_wacc"] == 1.0
    assert destination_wacc_set["finance.destination_wacc"] == 0.09


def test_an_override_naming_a_workbook_row_or_an_unknown_region_stops_the_build():
    with pytest.raises(ValueError) as error:
        build_from_rows({
            "workbook-row": {"sensitivity_overrides_JSON": json.dumps({
                "iron.fom_eur_per_t_capacity_per_year_bra": 80})},
            "unknown-region": {"sensitivity_overrides_JSON": json.dumps({"moe.opex_per_t_per_year_eur": "brazil"})}})
    assert "override key 'iron.fom_eur_per_t_capacity_per_year_bra' is not a model key" in str(error.value)
    assert "moe.opex_per_t_per_year_eur has no sensitivity level or region 'brazil'" in str(error.value)


def test_the_terms_decide_the_route_and_the_series_each_run_draws_on():
    """Islanded runs draw on the site's renewables; on the grid, an area without a market
    of its own runs once per zone. An import also draws on the destination's grid."""
    scenarios, _ = build_from_rows({
        "islanded-import": {"route": "ewin-eaf", "site_variant": "national-average"},
        "grid-domestic": {"scenario_group": "domestic", "steel_location": "AUS", "power_supply": "res_grid"},
    })

    lines = lambda scenario: list(scenarios.query("scenario == @scenario")[["run", "route", "tech", "variant", "area"]]
                                  .itertuples(index=False, name=None))
    assert lines("islanded-import") == [("AUS-ESP", "ew-eaf-export", "solar", "area-average", "AUS"),
                                        ("AUS-ESP", "ew-eaf-export", "grid", "emissions", "ESP")]
    assert lines("grid-domestic") == [("NSW1", "moe-eaf", "solar", "bestsite-p95", "NSW1"),
                                      ("NSW1", "moe-eaf", "grid", "emissions", "NSW1"),
                                      ("QLD1", "moe-eaf", "solar", "bestsite-p95", "QLD1"),
                                      ("QLD1", "moe-eaf", "grid", "emissions", "QLD1")]


def test_a_grid_only_run_draws_on_its_grids_alone_and_names_no_renewables():
    """Run per zone as res_grid is; a renewables column on such a row is a contradiction."""
    no_renewables = {"power_supply": "grid", "res_configuration": None, "res_technology": None,
                     "site_variant": None}
    scenarios, _ = build_from_rows({"grid-import": no_renewables})

    lines = list(scenarios[["run", "route", "tech", "variant", "area"]].itertuples(index=False, name=None))
    assert lines == [("NSW1-ESP", "moe-eaf-export", "grid", "emissions", "NSW1"),
                     ("NSW1-ESP", "moe-eaf-export", "grid", "emissions", "ESP"),
                     ("QLD1-ESP", "moe-eaf-export", "grid", "emissions", "QLD1"),
                     ("QLD1-ESP", "moe-eaf-export", "grid", "emissions", "ESP")]
    with pytest.raises(ValueError, match="res_technology on power_supply 'grid', which has no renewables"):
        build_from_rows({"grid-with-solar": {**no_renewables, "res_technology": "solar"}})


def test_a_mix_is_built_around_its_anchor_with_the_partner_cells_it_names_and_skipped_without():
    mix = {"res_configuration": "mix", "res_technology": None, "res_anchor": "wind-onshore"}
    scenarios, _ = build_from_rows({"named": {**mix, "res_partner_cells": 3},
                                    "unnamed": {**mix, "res_partner_cells": None}})

    assert set(scenarios["scenario"]) == {"named"}
    renewables = scenarios.query("tech != 'grid'")
    assert renewables[["tech", "variant"]].values.tolist() == [["wind-onshore", "anchor-colo-n3"]]


def test_unknown_or_contradicting_terms_stop_the_build():
    with pytest.raises(ValueError) as error:
        build_from_rows({"typo": {"route": "moe-eef"},
                         "grid-only": {"power_supply": "grid_only"},
                         "domestic-abroad": {"scenario_group": "domestic"},
                         "national-zones": {"scenario_group": "domestic", "steel_location": "AUS",
                                            "power_supply": "res_grid", "site_variant": "national-average"}})
    message = str(error.value)
    assert "unknown route 'moe-eef'" in message
    assert "unknown power_supply 'grid_only'" in message
    assert "scenario_group 'domestic' with iron in AUS and steel in ESP" in message
    assert "national-average on the grid of an area that runs per zone" in message


@pytest.fixture
def outputs_in_tmp(tmp_path, monkeypatch):
    """Run the generator on the real workbook, writing into tmp_path instead of config/."""
    overlays_dir, generated_dir = tmp_path / "overlays", tmp_path / "generated"
    overlays_dir.mkdir()
    monkeypatch.setattr(generate_inputs, "OVERLAYS_DIR", overlays_dir)
    monkeypatch.setattr(generate_inputs, "GENERATED_DIR", generated_dir)
    monkeypatch.setattr(generate_inputs, "ASSUMPTIONS_PATH", generated_dir / "assumptions.yaml")
    monkeypatch.setattr(generate_inputs, "SCENARIOS_PATH", generated_dir / "scenarios.csv")
    monkeypatch.setattr(sys, "argv", ["generate_inputs"])
    return tmp_path


def test_an_overlay_kept_by_hand_under_a_generated_scenario_name_stops_the_build(outputs_in_tmp):
    """Whatever the generated overlay holds, the solve would otherwise read the hand-kept one."""
    (outputs_in_tmp / "overlays" / "domestic-ESP-moe-eaf-solar-p95.yaml").write_text("finance:\n  origin_wacc: 0.05\n")

    with pytest.raises(ValueError, match=r"\['domestic-ESP-moe-eaf-solar-p95'\] .* are kept by hand"):
        generate_inputs.main()
    assert not (outputs_in_tmp / "generated").exists()


def test_the_template_values_standing_in_for_the_workbook_are_named_with_the_scenarios_they_reach(
        outputs_in_tmp, caplog):
    """A key every scenario's overlay sets (a regional one) is no stand-in, so no special case is needed."""
    generate_inputs.main()

    warning = next(record.message for record in caplog.records if "as the workbook gives none" in record.message)
    assert re.search(r"briquetting\.yield_t_per_t: \d+", warning)
    assert "finance.origin_wacc" not in warning


def test_every_workbook_row_is_read_or_listed_with_why_it_is_not(outputs_in_tmp, monkeypatch):
    """A row nobody accounted for stops the build, as does a listing the workbook or the reading has outgrown."""
    without_battery_hours = {row: reason for row, reason in UNREAD_ROWS.items() if row != "battery.max_hours"}
    monkeypatch.setattr(generate_inputs, "UNREAD_ROWS", without_battery_hours)
    with pytest.raises(ValueError, match=r"\['battery.max_hours'\] are neither read nor listed"):
        generate_inputs.main()

    with_read_and_absent_rows = UNREAD_ROWS | {"moe.el_mwh_per_t": "", "moe.no_such_row": ""}
    monkeypatch.setattr(generate_inputs, "UNREAD_ROWS", with_read_and_absent_rows)
    with pytest.raises(ValueError, match=r"lists \['moe.el_mwh_per_t', 'moe.no_such_row'\]"):
        generate_inputs.main()
    assert not (outputs_in_tmp / "generated").exists()


def test_an_override_naming_a_sensitivity_row_takes_that_level():
    _, overlays = build_from_rows({"solar-low-gas-high": {"sensitivity_overrides_JSON": json.dumps({
        "res.solar.capex_per_mw_eur": "low", "natural_gas.price_eur_per_mwh": "high"})}})

    assert overlays["solar-low-gas-high"]["res.solar.capex_per_mw_eur"] == 5.0
    assert overlays["solar-low-gas-high"]["natural_gas.price_eur_per_mwh"] == 7.0
