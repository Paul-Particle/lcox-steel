"""Unit tests for turning the input workbook's values into model keys.

The workbook rows are a stand-in for the `assumptions_inputs` tab of
data/assumptions/Assumptions_yaml_inputs.xlsx, indexed by its YAML key column.
"""

import pandas as pd
import pytest

from common._workbook_translation import translate_override
from common.generate_inputs import read_overrides, resolve_sources

WORKBOOK_ROWS = pd.DataFrame(
    [("plant.capex_usd_base", 1.13 * 1000, "USD/t"),
     ("plant.capex_usd_high", 1.13 * 2000, "USD/t"),
     ("plant.heat_mwh", 0.5, "MWh/t"),
     ("plant.drive_mwh", 0.25, "MWh/t"),
     ("plant.lifetime", 20, "years")],
    columns=["key", "value", "unit"],
).set_index("key")


def resolve_with(sources):
    return resolve_sources(WORKBOOK_ROWS, sources)


def test_values_are_converted_combined_and_kept_with_their_alternatives():
    resolved = resolve_with({
        "plant.capex_eur_per_t": {"base": "plant.capex_usd_base", "high": "plant.capex_usd_high"},
        "plant.opex_eur_per_t": {"mean": ["plant.capex_usd_base", "plant.capex_usd_high"]},
        "plant.total_mwh_per_t": {"sum": ["plant.heat_mwh", "plant.drive_mwh"]},
        "plant.lifetime_years": "plant.lifetime",
    })

    assert resolved["plant.capex_eur_per_t"] == pytest.approx({"base": 1000, "high": 2000})
    assert resolved["plant.opex_eur_per_t"] == pytest.approx({"default": 1500})
    assert resolved["plant.total_mwh_per_t"] == pytest.approx({"default": 0.75})
    assert resolved["plant.lifetime_years"] == {"default": 20}


def test_a_missing_row_stops_the_build():
    with pytest.raises(KeyError, match="plant.capex_eur_per_t"):
        resolve_with({"plant.capex_eur_per_t": "plant.capex_usd_low"})


def test_rows_in_different_units_stop_the_build():
    with pytest.raises(ValueError, match="one unit"):
        resolve_with({"plant.capex_eur_per_t": {"sum": ["plant.capex_usd_base", "plant.heat_mwh"]}})


def test_a_currency_the_key_does_not_name_stops_the_build():
    """A USD row under a key without `eur` (or the reverse) is a wrong row, not a value."""
    with pytest.raises(ValueError, match="currency"):
        resolve_with({"plant.lifetime_years": "plant.capex_usd_base"})
    with pytest.raises(ValueError, match="currency"):
        resolve_with({"plant.capex_eur_per_t": "plant.lifetime"})


def test_overrides_survive_excel_curly_quotes():
    assert read_overrides('{“finance.default_wacc”:“+30%”}') == {"finance.default_wacc": "+30%"}
    assert read_overrides(float("nan")) == {}


def test_workbook_override_keys_translate_to_model_keys_and_alternatives():
    assert translate_override("res.solar.capex_usd_per_mw", "res.solar.capex_usd_per_mw_low") == [
        ("res.solar.capex_per_mw_eur", "low")]
    assert translate_override("electrolyser.capex_usd_per_mw", "electrolyser.capex_usd_per_mw_low") == [
        ("electrolyser.capex_per_mw_eur", "electrolyser.capex_usd_per_mw_low")]


def test_one_workbook_key_reaches_every_model_key_it_feeds():
    assert translate_override("finance.default_wacc", "+30%") == [
        ("finance.origin_wacc", "+30%"), ("finance.destination_wacc", "+30%")]
    assert [path for path, _ in translate_override("dri.capex_usd_per_t", 400)] == [
        "dri-h2.capex_per_t_per_year_eur", "dri-mix.capex_per_t_per_year_eur", "dri-ng.capex_per_t_per_year_eur"]
