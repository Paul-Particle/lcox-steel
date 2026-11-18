"""The solve runs once per grid tariff regime and keeps the cheapest; the report
charges what the area charges on consumption, with a grid or without.

Built on the two-day grid-only network of test_grid_only_run: the plant buys
all its power, so its utilisation hours and annual energy are what the regimes
are cut on. The charges are dicts in the shape build_electricity_charges writes.
"""

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import REPO_ROOT
from scripts.solve import _grid_tariff
from scripts.solve.build_network import build_network
from scripts.solve.solve_network import _solve_cheapest_regime, _tie_battery_inverter
from scripts.viz import compile_report

pd.options.mode.string_storage = "python"
TEMPLATE = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
INDEX = pd.date_range("2025-03-01", periods=48, freq="h")
PRICE = pd.Series(np.where(INDEX.hour < 12, 30.0, 120.0), index=INDEX)
GRID = TEMPLATE["grid"]
MIX = pd.DataFrame({"gas": 400.0, "hard_coal": 100.0, "solar": np.where(INDEX.hour < 12, 500.0, 0.0)},
                   index=INDEX)
RUN = {"run": "DEU", "origin": "DEU", "destination": "DEU", "route": "h2-dri-eaf",
       "start_date": "20250301", "end_date": "20250302"}


def _regime(name, capacity, energy, fixed=0.0, h=(None, None), e=(None, None), components=None):
    return {"name": name, "h_min": h[0], "h_max": h[1], "e_min_mwh": e[0], "e_max_mwh": e[1],
            "capacity_eur_per_mw_per_year": capacity, "energy_eur_per_mwh": energy,
            "fixed_eur_per_year": fixed, "components": components or {}}


def _charges(area, *regimes, consumption=0.0, reference=None, reference_consumption=0.0):
    """An area's electricity charges: grid regimes, a consumption charge, and a
    reference consumer if `reference` gives its regimes."""
    return {
        "area": area, "source": "test", "consumer": "steel",
        "reference_consumer": "standard" if reference is not None else None,
        "component_groups": {"leistungspreis": "grid_fees", "arbeitspreis": "grid_fees",
                             "messstellenbetrieb": "metering", "strompreiskompensation": "subsidy"},
        "grid": {"regimes": list(regimes), "reference": reference},
        "consumption": {
            "eur_per_mwh": consumption,
            "components": {"strompreiskompensation": consumption} if consumption else {},
            "reference": ({"eur_per_mwh": reference_consumption, "components": {}}
                          if reference is not None else None),
        },
    }


def _solve(route="h2-dri-eaf", origin=None, destination=None):
    is_export = route.endswith("-export")
    n = build_network(route, TEMPLATE, pd.DataFrame(index=INDEX), price_series=PRICE,
                      transport_legs={"sea": 17_965} if is_export else None,
                      destination_price=PRICE * 1.1 if is_export else None)
    tariffs = {_grid_tariff.GRID_IMPORT: origin, _grid_tariff.DESTINATION_SUPPLY: destination}
    prices = {_grid_tariff.GRID_IMPORT: PRICE, _grid_tariff.DESTINATION_SUPPLY: PRICE * 1.1}
    waccs = {_grid_tariff.GRID_IMPORT: TEMPLATE["finance"]["origin_wacc"],
             _grid_tariff.DESTINATION_SUPPLY: TEMPLATE["finance"]["destination_wacc"]}
    result = _solve_cheapest_regime(n, "test", {}, GRID, tariffs, prices, waccs)
    n.meta["grid_tariff"] = result
    return n, result


# -- the solve: one run per regime, the cheapest kept

def test_an_area_without_a_workbook_is_billed_as_before():
    n, result = _solve()
    reference = build_network("h2-dri-eaf", TEMPLATE, pd.DataFrame(index=INDEX), price_series=PRICE)
    reference.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)

    bill = result["meters"]["grid_import"]
    assert bill["regime"] == "config" and len(result["candidates"]) == 1
    assert float(n.objective) == pytest.approx(float(reference.objective), rel=1e-6)


def test_the_cheapest_regime_is_kept_and_its_box_holds():
    # The discount needs ≥ 6000 h; the plant is flat enough to reach it.
    charges = _charges("DEU",
                       _regime("full", 50_000, 10.0),
                       _regime("discount", 5_000, 1.0, fixed=1_000.0, h=(6000, None)))
    n, result = _solve(origin=charges)
    bill = result["meters"]["grid_import"]

    assert bill["regime"] == "discount"
    assert bill["utilisation_hours"] >= 6000 * (1 - 1e-6)
    totals = {c["regimes"]: c["total_eur_per_year"] for c in result["candidates"]}
    assert totals["grid_import=discount"] < totals["grid_import=full"]
    # The fixed €/yr counts in the comparison, not in the objective.
    assert totals["grid_import=discount"] == pytest.approx(float(n.objective) + 1_000.0)


def test_a_regime_the_plant_cannot_reach_is_skipped():
    charges = _charges("DEU",
                       _regime("tiny", 0.0, 0.0, e=(None, 1.0)),
                       _regime("full", 50_000, 10.0))
    _, result = _solve(origin=charges)

    assert result["meters"]["grid_import"]["regime"] == "full"
    assert [c["total_eur_per_year"] is None for c in result["candidates"]] == [True, False]


def test_the_annual_energy_counts_a_short_window_as_a_year():
    n, result = _solve()
    bill = result["meters"]["grid_import"]
    two_days = float(n.generators_t.p["grid_import"].sum())
    assert bill["energy_mwh_per_year"] == pytest.approx(two_days * 8760 / 48)


def test_the_consumption_charge_stays_out_of_the_objective():
    plain, _ = _solve(origin=_charges("DEU", _regime("full", 50_000, 10.0)))
    subsidised, _ = _solve(origin=_charges("DEU", _regime("full", 50_000, 10.0), consumption=-30.0))
    assert float(subsidised.objective) == pytest.approx(float(plain.objective), rel=1e-9)


def test_an_export_tries_every_pair_of_regimes():
    origin = _charges("AUS", _regime("a1", 50_000, 10.0), _regime("a2", 60_000, 5.0))
    destination = _charges("DEU", _regime("d1", 50_000, 10.0), _regime("d2", 40_000, 20.0),
                           _regime("d3", 45_000, 15.0))
    _, result = _solve("moe-eaf-export", origin, destination)

    assert len(result["candidates"]) == 6
    assert set(result["meters"]) == {"grid_import", "destination_supply"}
    best = min(c["total_eur_per_year"] for c in result["candidates"])
    chosen = f"grid_import={result['meters']['grid_import']['regime']}, " \
             f"destination_supply={result['meters']['destination_supply']['regime']}"
    assert {c["regimes"]: c["total_eur_per_year"] for c in result["candidates"]}[chosen] == best


# -- the report: the grid bill and the consumption charge reach the stack and the table

def _billed(consumption):
    """Charges with a component per kind of grid charge, and a dearer reference consumer."""
    regime = _regime("full", 50_000, 10.0, fixed=1_000.0, components={
        "leistungspreis": {"capacity_eur_per_mw_per_year": 50_000},
        "arbeitspreis": {"energy_eur_per_mwh": 10.0},
        "messstellenbetrieb": {"fixed_eur_per_year": 1_000.0}})
    return _charges("DEU", regime, consumption=consumption,
                    reference=[_regime("standard", 60_000, 30.0, fixed=1_000.0)])


def _report(consumption):
    charges = _billed(consumption)
    n, _ = _solve(origin=charges)
    s = compile_report.extract_summary(n, "tariff", RUN, TEMPLATE, None, MIX, None,
                                       {"origin": charges, "destination": None})
    return n, s, charges


def test_the_grid_bill_is_the_capacity_fee_and_the_fee():
    _, s, _ = _report(consumption=-20.0)
    assert s["grid_tariff_regime"] == "full"
    assert (s["cost_grid_capacity_fee_eur_per_t_steel"] + s["cost_grid_fee_eur_per_t_steel"]
            == pytest.approx(s["grid_charges_eur_per_t_steel"]))
    # The consumption charge is a line of its own, and the leaves still close on LCOS.
    assert s["cost_consumption_charges_eur_per_t_steel"] < 0
    leaves = sum(s[f"cost_{leaf}_eur_per_t_steel"] for leaf in compile_report.LEAF_COSTS)
    assert leaves == pytest.approx(s["lcos_eur_per_t_steel"])


def test_the_subsidy_lowers_the_lcos_by_what_it_paid_on_consumption():
    n, with_subsidy, _ = _report(consumption=-20.0)
    _, without, _ = _report(consumption=0.0)
    steel_t = with_subsidy["steel_produced_mt"] * 1e6
    consumed = with_subsidy["el_users_gwh"] * 1e3
    assert without["lcos_eur_per_t_steel"] - with_subsidy["lcos_eur_per_t_steel"] == pytest.approx(
        20.0 * consumed / steel_t)
    assert with_subsidy["cost_consumption_charges_eur_per_t_steel"] == pytest.approx(
        -20.0 * consumed / steel_t)


def test_the_subsidy_is_paid_on_what_the_plant_drew_not_on_what_it_imported():
    """The processes' draw, as the report counts it: the battery's losses are not
    consumption, so a plant that cycles its battery imports more than it consumed."""
    n, s, _ = _report(consumption=-20.0)
    consumed = compile_report._consumption_mwh(n, at_destination=False)
    assert consumed <= n.meta["grid_tariff"]["meters"]["grid_import"]["energy_mwh_per_year"] * (1 + 1e-9)
    assert consumed == pytest.approx(s["el_users_gwh"] * 1e3, rel=1e-6)


def test_the_reference_consumer_is_billed_on_the_same_profile():
    n, s, _ = _report(consumption=-20.0)
    meter = n.meta["grid_tariff"]["meters"]["grid_import"]
    steel_t = s["steel_produced_mt"] * 1e6
    expected = 60_000 * meter["p_nom_mw"] + 30.0 * meter["energy_mwh_per_year"] + 1_000.0
    assert s["grid_charges_reference_eur_per_t_steel"] * steel_t == pytest.approx(expected)
    assert s["grid_charges_reference_eur_per_t_steel"] > s["grid_charges_eur_per_t_steel"]
    # The reference consumer is charged nothing on consumption: no subsidy.
    assert s["consumption_charges_reference_eur_per_t_steel"] == 0.0


def test_the_charges_table_lists_every_component_and_closes_on_each_total():
    n, _, charges = _report(consumption=-20.0)
    table = compile_report.electricity_charges(n, "tariff", "DEU", {"origin": charges})
    own = table[table["consumer"] == "steel"]
    for base, components in (("grid", ["leistungspreis", "arbeitspreis", "messstellenbetrieb"]),
                             ("consumption", ["strompreiskompensation"])):
        rows = own[own["base"] == base]
        parts = rows[rows["component"] != "total"]
        assert list(parts["component"]) == components
        assert parts["eur_per_year"].sum() == pytest.approx(
            rows.loc[rows["component"] == "total", "eur_per_year"].item())
    # One reference total per base.
    assert sorted(table.loc[table["consumer"] == "standard", "base"]) == ["consumption", "grid"]


def test_an_export_charges_each_side_on_its_own_consumption():
    """The furnace abroad is the destination's consumption, everything else the
    origin's; together they are every user's draw."""
    origin = _charges("AUS", _regime("a", 50_000, 10.0))
    destination = _charges("DEU", _regime("d", 50_000, 10.0), consumption=-20.0)
    n, _ = _solve("moe-eaf-export", origin, destination)
    run = {**RUN, "run": "AUS-DEU", "origin": "AUS", "route": "moe-eaf-export"}
    s = compile_report.extract_summary(n, "tariff", run, TEMPLATE, {"sea": 17_965}, MIX, MIX,
                                       {"origin": origin, "destination": destination})

    at_origin = compile_report._consumption_mwh(n, at_destination=False)
    at_destination = compile_report._consumption_mwh(n, at_destination=True)
    assert at_origin + at_destination == pytest.approx(s["el_users_gwh"] * 1e3, rel=1e-6)
    assert at_destination == pytest.approx(s["el_eaf_gwh"] * 1e3, rel=1e-6)
    # Only the destination's area charges consumption, so only its side pays.
    steel_t = s["steel_produced_mt"] * 1e6
    assert s["cost_consumption_charges_eur_per_t_steel"] == 0.0
    assert s["cost_destination_consumption_charges_eur_per_t_steel"] == pytest.approx(
        -20.0 * at_destination / steel_t)


def test_a_plant_on_its_own_renewables_is_charged_on_its_consumption_too():
    """No grid, so no tariff regime: the area's consumption charge still applies."""
    sun = pd.DataFrame({"solar": np.where((INDEX.hour >= 6) & (INDEX.hour < 18), 0.7, 0.0)},
                       index=INDEX)
    n = build_network("h2-dri-eaf", TEMPLATE, sun)
    status, condition = n.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)
    assert condition == "optimal"
    charges = _charges("DEU", _regime("full", 50_000, 10.0), consumption=-20.0)
    s = compile_report.extract_summary(n, "res-only", RUN, TEMPLATE, None, None, None,
                                       {"origin": charges, "destination": None})

    steel_t = s["steel_produced_mt"] * 1e6
    assert "grid_import" not in n.generators.index
    assert s["cost_consumption_charges_eur_per_t_steel"] == pytest.approx(
        -20.0 * s["el_users_gwh"] * 1e3 / steel_t)
    leaves = sum(s[f"cost_{leaf}_eur_per_t_steel"] for leaf in compile_report.LEAF_COSTS)
    assert leaves == pytest.approx(s["lcos_eur_per_t_steel"])
