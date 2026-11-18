"""Unit tests for turning a tariff sheet into LP regimes.

The check that matters: for any plant size and annual energy, the cheapest
regime whose (h, E) box holds the plant bills exactly what the tariff says.
The rows are a synthetic copy of the German tariff, so the tests need no
workbook.
"""

import math

import pandas as pd
import pytest

from scripts.grid import _tariff

COLUMNS = ["year", "component", "model_variable", "value", "unit", "period", "h_min", "h_max",
           "e_min_mwh", "tier_from_mwh", "tier_to_mwh", "offpeak_min", "offpeak_max",
           "consumer", "applies_to", "params"]
NA = float("nan")


def _row(component, variable, value, h=(NA, NA), e_min=NA, tier=(NA, NA), consumer=NA,
         applies_to=NA, params=NA):
    return [2026, component, variable, value, _tariff.UNITS[variable], NA, h[0], h[1], e_min,
            tier[0], tier[1], NA, NA, consumer, applies_to, params]


def _sheet(*rows):
    return pd.DataFrame(list(rows), columns=COLUMNS)


GERMAN = _sheet(
    _row("industriestrompreis", "subsidy_eur_per_mwh", -18.72, consumer="steel", params="one_of=ets_or_isp"),
    _row("strompreiskompensation", "subsidy_eur_per_mwh", -30.79, consumer="steel", params="one_of=ets_or_isp"),
    _row("baukostenzuschuss", "connection_capex_eur_per_mw", 0.0),
    _row("leistungspreis", "fee_eur_per_mw_per_year", 11390, h=(0, 2500)),
    _row("leistungspreis", "fee_eur_per_mw_per_year", 53060, h=(2500, NA)),
    _row("arbeitspreis", "fee_eur_per_mwh", 23.6, h=(0, 2500)),
    _row("arbeitspreis", "fee_eur_per_mwh", 6.9, h=(2500, NA)),
    _row("s19_discount", "multiplier", 0.2, h=(7000, 7500), e_min=10000, applies_to="leistungspreis, arbeitspreis"),
    _row("s19_discount", "multiplier", 0.15, h=(7500, 8000), e_min=10000, applies_to="leistungspreis, arbeitspreis"),
    _row("s19_discount", "multiplier", 0.1, h=(8000, NA), e_min=10000, applies_to="leistungspreis, arbeitspreis"),
    _row("messstellenbetrieb", "fee_eur_per_year", 3629),
    _row("kwkg_surcharge", "fee_eur_per_mwh", 4.46, tier=(NA, 1000)),
    _row("kwkg_surcharge", "fee_eur_per_mwh", 0.669, tier=(1000, NA)),
    _row("s19_surcharge", "fee_eur_per_mwh", 15.59, tier=(NA, 1000)),
    _row("s19_surcharge", "fee_eur_per_mwh", 0.25, tier=(1000, NA), consumer="steel"),
    _row("s19_surcharge", "fee_eur_per_mwh", 0.5, tier=(1000, NA), consumer="standard"),
    _row("stromsteuer", "fee_eur_per_mwh", 0.0, consumer="steel"),
    _row("stromsteuer", "fee_eur_per_mwh", 20.5, consumer="standard"),
)


def _german_bill(p_mw, e_mwh):
    """The steel plant's bill, written out from the tariff by hand."""
    h = e_mwh / p_mw
    capacity, energy = (11390, 23.6) if h < 2500 else (53060, 6.9)
    if e_mwh >= 10000:
        for lo, factor in ((8000, 0.1), (7500, 0.15), (7000, 0.2)):
            if h >= lo:
                capacity, energy = capacity * factor, energy * factor
                break
    tiered = sum(first * min(e_mwh, 1000) + rest * max(e_mwh - 1000, 0)
                 for first, rest in ((4.46, 0.669), (15.59, 0.25)))
    return capacity * p_mw + energy * e_mwh + tiered + 3629


def _cheapest(regimes, p_mw, e_mwh):
    h = e_mwh / p_mw
    fits = [r for r in regimes
            if r.h_min <= h <= r.h_max and r.e_min_mwh <= e_mwh <= r.e_max_mwh]
    return min(r.capacity_eur_per_mw_per_year * p_mw + r.energy_eur_per_mwh * e_mwh
               + r.fixed_eur_per_year for r in fits)


def _regimes(sheet):
    return _tariff.build_regimes(_tariff._select_rows(sheet, "steel", "test"), "test")


@pytest.mark.parametrize("p_mw, e_mwh", [
    (1, 500), (1, 1000), (1, 2400), (1, 2600),
    (10, 20000), (10, 69000), (10, 72000), (10, 77000), (10, 85000),
    (1.2, 8500), (100, 600000), (100, 790000), (500, 4_380_000),
])
def test_the_cheapest_fitting_regime_bills_the_tariff(p_mw, e_mwh):
    assert _cheapest(_regimes(GERMAN), p_mw, e_mwh) == pytest.approx(_german_bill(p_mw, e_mwh))


def test_the_discount_waits_for_its_energy_threshold():
    """8500 h below 10 GWh pays the full rate: no §19 regime fits."""
    regimes = _regimes(GERMAN)
    names = [r.name for r in regimes if r.h_min <= 8500 <= r.h_max and r.e_min_mwh <= 8500 <= r.e_max_mwh]
    assert not any(name.startswith("s19_discount") for name in names)


def test_the_regimes_are_the_bands_the_tiers_and_the_discounts():
    names = [r.name for r in _regimes(GERMAN)]
    assert names == [
        "h0-2500_e0-1000", "h0-2500_e1000-", "h2500-_e0-1000", "h2500-_e1000-",
        "s19_discount_h7000-7500_e10000-", "s19_discount_h7500-8000_e10000-",
        "s19_discount_h8000-_e10000-",
    ]


def test_the_larger_of_two_exclusive_subsidies_is_kept():
    rows = _tariff._select_rows(GERMAN, "steel", "test")
    assert _tariff._subsidy(rows) == (-30.79, ["strompreiskompensation"])


def test_only_the_consumers_rows_are_billed():
    standard = _tariff.build_regimes(_tariff._select_rows(GERMAN, "standard", "test"), "test")
    steel = _regimes(GERMAN)
    # 0.5 instead of 0.25 on the §19 surcharge, and 20.5 of electricity tax.
    assert standard[1].energy_eur_per_mwh - steel[1].energy_eur_per_mwh == pytest.approx(20.75)


def test_a_multiplier_above_one_is_refused():
    sheet = _sheet(
        _row("leistungspreis", "fee_eur_per_mw_per_year", 100.0),
        _row("penalty", "multiplier", 1.5, h=(0, 1000), applies_to="leistungspreis"),
    )
    with pytest.raises(ValueError, match="only discounts"):
        _regimes(sheet)


def test_a_multiplier_on_an_unknown_component_is_refused():
    sheet = _sheet(
        _row("leistungspreis", "fee_eur_per_mw_per_year", 100.0),
        _row("s19_discount", "multiplier", 0.5, h=(7000, NA), applies_to="leistungpreis"),
    )
    with pytest.raises(ValueError, match="unknown components"):
        _regimes(sheet)


def test_a_gap_in_a_components_hour_bands_is_refused():
    sheet = _sheet(
        _row("leistungspreis", "fee_eur_per_mw_per_year", 100.0, h=(0, 2500)),
        _row("leistungspreis", "fee_eur_per_mw_per_year", 50.0, h=(3000, NA)),
    )
    with pytest.raises(ValueError, match="no rate for 2500–3000 h"):
        _regimes(sheet)


def test_a_wrong_unit_is_refused():
    sheet = _sheet(_row("arbeitspreis", "fee_eur_per_mwh", 6.9))
    sheet.loc[0, "unit"] = "ct/kWh"
    with pytest.raises(ValueError, match="expected 'EUR/MWh'"):
        _regimes(sheet)


def test_an_off_peak_window_is_not_silently_dropped():
    sheet = _sheet(_row("arbeitspreis", "fee_eur_per_mwh", 6.9))
    sheet.loc[0, "offpeak_min"] = 22
    with pytest.raises(NotImplementedError, match="offpeak_min"):
        _regimes(sheet)


def test_open_bounds_are_written_as_null():
    regime = _regimes(GERMAN)[-1].to_dict()
    assert regime["h_max"] is None and regime["e_max_mwh"] is None
    assert regime["capacity_eur_per_mw_per_year"] == pytest.approx(5306)
    assert not math.isinf(regime["e_min_mwh"])
