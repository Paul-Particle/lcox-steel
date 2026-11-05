"""Unit tests for turning the input workbook's values into model keys.

The workbook rows are a stand-in for the `assumptions_inputs` tab of
data/assumptions/Assumptions_yaml_inputs.xlsx, indexed by its YAML key column.
"""

import json

import pandas as pd
import pytest

from common._workbook_translation import WORKBOOK_SOURCES, translate_override
from common.generate_inputs import (
    REGIONAL_SOURCES,
    build_scenarios,
    flatten,
    read_overrides,
    resolve_sources,
    set_nested,
)

WORKBOOK_ROWS = pd.DataFrame(
    [("plant.capex_eur_base", 1000, "EUR/t"),
     ("plant.capex_eur_high", 2000, "EUR/t"),
     ("plant.capex_usd_base", 1130, "USD/t"),
     ("plant.heat_mwh", 0.5, "MWh/t"),
     ("plant.drive_mwh", 0.25, "MWh/t"),
     ("plant.lifetime", 20, "years")],
    columns=["key", "value", "unit"],
).set_index("key")


def resolve_with(sources):
    return resolve_sources(WORKBOOK_ROWS, sources)


def test_values_are_combined_and_kept_with_their_alternatives():
    resolved = resolve_with({
        "plant.capex_eur_per_t": {"base": "plant.capex_eur_base", "high": "plant.capex_eur_high"},
        "plant.opex_eur_per_t": {"mean": ["plant.capex_eur_base", "plant.capex_eur_high"]},
        "plant.total_mwh_per_t": {"sum": ["plant.heat_mwh", "plant.drive_mwh"]},
        "plant.lifetime_years": "plant.lifetime",
    })

    assert resolved["plant.capex_eur_per_t"] == pytest.approx({"base": 1000, "high": 2000})
    assert resolved["plant.opex_eur_per_t"] == pytest.approx({"default": 1500})
    assert resolved["plant.total_mwh_per_t"] == pytest.approx({"default": 0.75})
    assert resolved["plant.lifetime_years"] == {"default": 20}


def test_a_missing_row_stops_the_build():
    with pytest.raises(KeyError, match="plant.capex_eur_per_t"):
        resolve_with({"plant.capex_eur_per_t": "plant.capex_eur_low"})


def test_rows_in_different_units_stop_the_build():
    with pytest.raises(ValueError, match="one unit"):
        resolve_with({"plant.capex_eur_per_t": {"sum": ["plant.capex_eur_base", "plant.heat_mwh"]}})


def test_a_currency_the_key_does_not_name_stops_the_build():
    """The model is in EUR: a USD row, a EUR row under a key without `eur`, or the reverse is a wrong row."""
    with pytest.raises(ValueError, match="currency"):
        resolve_with({"plant.capex_eur_per_t": "plant.capex_usd_base"})
    with pytest.raises(ValueError, match="currency"):
        resolve_with({"plant.lifetime_years": "plant.capex_eur_base"})
    with pytest.raises(ValueError, match="currency"):
        resolve_with({"plant.capex_eur_per_t": "plant.lifetime"})


def test_overrides_survive_excel_curly_quotes():
    assert read_overrides('{“finance.default_wacc”:“+30%”}') == {"finance.default_wacc": "+30%"}
    assert read_overrides(float("nan")) == {}


def test_workbook_override_keys_translate_to_model_keys_and_alternatives():
    assert translate_override("res.solar.capex_eur_per_mw", "res.solar.capex_eur_per_mw_low", WORKBOOK_SOURCES) == [
        ("res.solar.capex_per_mw_eur", "low")]
    assert translate_override("electrolyser.capex_eur_per_mw", "electrolyser.capex_eur_per_mw_low",
                              WORKBOOK_SOURCES) == [("electrolyser.capex_per_mw_eur", "electrolyser.capex_eur_per_mw_low")]


def test_one_workbook_key_reaches_every_model_key_it_feeds():
    assert translate_override("finance.default_wacc", "+30%", WORKBOOK_SOURCES) == [
        ("finance.origin_wacc", "+30%"), ("finance.destination_wacc", "+30%")]
    assert [path for path, _ in translate_override("dri.capex_eur_per_t", 400, WORKBOOK_SOURCES)] == [
        "dri-h2.capex_per_t_per_year_eur", "dri-mix.capex_per_t_per_year_eur", "dri-ng.capex_per_t_per_year_eur"]


EXPORT_ROW = {
    "iron_location": "AUS", "steel_location": "ESP", "route": "moe-eaf", "power_supply": "islanded",
    "res_configuration": "single", "res_technology": "solar", "res_anchor": None,
    "site_variant": "bestsite-p95", "start_date": 20250101, "end_date": 20251231,
}
REGIONAL_ALTERNATIVES_BY_PATH = {path: {"noneu": 1.0, "aus": 1.0, "global": 1.0, "bra": 3.0, "eu": 2.0}
                                 for path in REGIONAL_SOURCES}
EXPORT_AREAS = {"AUS": {}, "ESP": {"market": "entsoe"}}


def build_export_overlays(overrides_by_scenario: dict) -> dict:
    """The overlays of AUS->ESP exports with the given overrides, on a base of zeros."""
    scenario_rows = pd.DataFrame([{**EXPORT_ROW, "scenario_id": scenario, "sensitivity_overrides_JSON": json.dumps(overrides)}
                                  for scenario, overrides in overrides_by_scenario.items()])
    assumptions = {}
    for path in REGIONAL_SOURCES:
        set_nested(assumptions, path, 0.0)
    _, overlays = build_scenarios(scenario_rows, REGIONAL_ALTERNATIVES_BY_PATH, assumptions, EXPORT_AREAS)
    flat_overlays = {scenario: flatten(overlay) for scenario, overlay in overlays.items()}
    return flat_overlays


def test_an_export_takes_iron_side_values_from_its_iron_location_and_the_furnace_from_its_steel_location():
    """The overrides apply on top: to the region's value, or switching to another region's row.
    A regional row reaches only the keys the scenario takes from it."""
    overlays = build_export_overlays({
        "wacc-up": {"finance.default_wacc": "+30%"},
        "rows-overridden": {"battery.capex_eur_per_mw_global": "battery.capex_eur_per_mw_eu",
                            "iron.fom_eur_per_t_capacity_per_year_aus": "+10%"},
        "eu-wacc-set": {"finance.other.eu_wacc": 0.09},
    })

    wacc_up, rows_overridden, eu_wacc_set = overlays["wacc-up"], overlays["rows-overridden"], overlays["eu-wacc-set"]
    assert wacc_up["finance.origin_wacc"] == pytest.approx(1.3)
    assert wacc_up["finance.destination_wacc"] == pytest.approx(2.6)
    assert wacc_up["moe.opex_per_t_per_year_eur"] == 1.0
    assert wacc_up["battery.capex_per_mw_eur"] == 1.0
    assert wacc_up["eaf.opex_per_t_per_year_eur"] == 2.0
    assert rows_overridden["moe.opex_per_t_per_year_eur"] == pytest.approx(1.1)
    assert rows_overridden["battery.capex_per_mw_eur"] == 2.0
    assert rows_overridden["battery.capex_per_mwh_eur"] == 1.0
    assert eu_wacc_set["finance.origin_wacc"] == 1.0
    assert eu_wacc_set["finance.destination_wacc"] == 0.09


def test_an_override_of_a_row_the_scenario_takes_no_key_from_stops_the_build():
    with pytest.raises(ValueError, match="a row of a region this scenario does not take"):
        build_export_overlays({"brazil-fom": {"iron.fom_eur_per_t_capacity_per_year_bra": 80}})
