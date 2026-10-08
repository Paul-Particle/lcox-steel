"""Unit tests for writing config/assumptions.yaml from the template and the input workbook.

The workbook is a stand-in laid out like data/assumptions/Assumptions_yaml_inputs.xlsx:
an `assumptions_inputs` tab with one row per value and its register row ID.
"""

import pandas as pd
import pytest
import yaml

from common.generate_inputs import (
    USD_PER_EUR,
    format_number,
    read_template_references,
    read_workbook_inputs,
    render_assumptions,
    resolve_reference,
)

WORKBOOK_ROWS = [
    # YAML key, value, unit, confirmed row ID
    ("plant.capex_usd_base", 1.13 * 1000, "USD/t", 1),
    ("plant.capex_usd_high", 1.13 * 2000, "USD/t", 2),
    ("plant.heat_mwh", 0.5, "MWh/t", 3),
    ("plant.drive_mwh", 0.25, "MWh/t", 4),
    ("plant.lifetime", 20, "years", 5),
    ("plant.p_min_pu", 0, "0", None),
]

TEMPLATE = """\
# Explains the template, and is not copied over.

plant:
  capex_eur_per_t: !register {base: 1, high: 2}   # kept comment
  opex_eur_per_t: !register {mean: [1, 2]}
  total_mwh_per_t: !register {sum: [3, 4]}
  lifetime_years: !register 5
  freight: { eur_per_t: 0, eur_per_t_km: !register 1 }
  untouched: 7
"""


@pytest.fixture
def workbook_inputs(tmp_path):
    workbook_path = tmp_path / "inputs.xlsx"
    (
        pd.DataFrame(WORKBOOK_ROWS, columns=["YAML key", "Value", "Unit", "Confirmed row ID"])
        .assign(**{"Old YAML key": None})
        .to_excel(workbook_path, sheet_name="assumptions_inputs", index=False)
    )
    return read_workbook_inputs(workbook_path)


def resolve_template(template, workbook_inputs):
    references, _ = read_template_references(template)
    for reference in references:
        resolve_reference(reference, workbook_inputs)
    return references


def test_values_are_converted_combined_and_written(workbook_inputs):
    rendered_text = render_assumptions(TEMPLATE, resolve_template(TEMPLATE, workbook_inputs))
    plant = yaml.safe_load(rendered_text)["plant"]

    assert USD_PER_EUR == 1.13
    assert plant["capex_eur_per_t"] == pytest.approx(1000)
    assert plant["opex_eur_per_t"] == pytest.approx(1500)
    assert plant["total_mwh_per_t"] == pytest.approx(0.75)
    assert plant["lifetime_years"] == 20
    assert plant["freight"] == {"eur_per_t": 0, "eur_per_t_km": pytest.approx(1000)}
    assert plant["untouched"] == 7


def test_alternatives_are_kept_beside_the_default(workbook_inputs):
    references = {r.path: r for r in resolve_template(TEMPLATE, workbook_inputs)}

    assert references["plant.capex_eur_per_t"].values == pytest.approx({"base": 1000, "high": 2000})


def test_comments_carry_the_rows_and_the_alternatives(workbook_inputs):
    rendered_text = render_assumptions(TEMPLATE, resolve_template(TEMPLATE, workbook_inputs))

    assert "Explains the template" not in rendered_text
    assert "# register base = 1: plant.capex_usd_base 1130, plant.capex_usd_high 2260 USD/t" in rendered_text
    assert "# alternatives in EUR: high = 2000" in rendered_text
    assert "# register sum of 3, 4: plant.heat_mwh 0.5, plant.drive_mwh 0.25 MWh/t" in rendered_text
    assert "# kept comment" in rendered_text


def test_untagged_values_are_listed_as_backstop_values():
    _, backstop_paths = read_template_references(TEMPLATE)

    assert backstop_paths == ["plant.freight.eur_per_t", "plant.untouched"]


def test_an_unknown_row_stops_the_build(workbook_inputs):
    with pytest.raises(KeyError, match="plant.capex_eur_per_t"):
        resolve_template("plant:\n  capex_eur_per_t: !register 99\n", workbook_inputs)


def test_rows_in_different_units_stop_the_build(workbook_inputs):
    with pytest.raises(ValueError, match="different units"):
        resolve_template("plant:\n  capex_eur_per_t: !register {sum: [1, 3]}\n", workbook_inputs)


def test_a_currency_the_key_does_not_name_stops_the_build(workbook_inputs):
    """A USD row under a key without `eur` (or the reverse) is a wrong ID, not a value."""
    with pytest.raises(ValueError, match="currency"):
        resolve_template("plant:\n  lifetime_years: !register 1\n", workbook_inputs)
    with pytest.raises(ValueError, match="currency"):
        resolve_template("plant:\n  capex_eur_per_t: !register 5\n", workbook_inputs)


def test_numbers_are_never_written_in_exponent_notation():
    assert format_number(2593805.3) == "2593810"
    assert format_number(0.0000343037) == "0.0000343037"
    assert format_number(20.0) == "20"
