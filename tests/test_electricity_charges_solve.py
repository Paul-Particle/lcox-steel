"""The solve runs once per grid tariff regime and keeps the cheapest, and prices
the grid connection from the workbook.

Built on the two-day grid-only network of test_grid_only_run: the plant buys
all its power, so its utilisation hours and annual energy are what the regimes
are cut on. The charges are dicts in the shape `_tariff.for_consumer` gives the
solve. test_electricity_charges_report builds on the helpers here.
"""

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import REPO_ROOT
from scripts.solve import _grid_tariff
from scripts.solve.build_network import build_network
from scripts.solve.solve_network import _solve_cheapest_regime, _tie_battery_inverter

pd.options.mode.string_storage = "python"
TEMPLATE = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
INDEX = pd.date_range("2025-03-01", periods=48, freq="h")
PRICE = pd.Series(np.where(INDEX.hour < 12, 30.0, 120.0), index=INDEX)
GRID = TEMPLATE["grid"]


def _regime(name, capacity, energy, fixed=0.0, h=(None, None), e=(None, None), components=None):
    return {"name": name, "h_min": h[0], "h_max": h[1], "e_min_mwh": e[0], "e_max_mwh": e[1],
            "capacity_eur_per_mw_per_year": capacity, "energy_eur_per_mwh": energy,
            "fixed_eur_per_year": fixed, "components": components or {}}


def _charges(area, *regimes, consumption=0.0, connection=None, consumer="iron_steel"):
    """An area's electricity charges as one consumer pays them, the shape
    `_tariff.for_consumer` gives the solve and the report: grid regimes, a
    consumption charge, and the connection's one-off costs if `connection`
    gives them (else the config's)."""
    return {
        "area": area, "source": "test", "consumer": consumer,
        "component_groups": {"leistungspreis": "grid_fees", "arbeitspreis": "grid_fees",
                             "messstellenbetrieb": "metering", "strompreiskompensation": "subsidy",
                             "netzanschlusskosten": "grid_connection",
                             "baukostenzuschuss": "grid_connection"},
        "connection": connection,
        "grid": {"regimes": list(regimes)},
        "consumption": {
            "eur_per_mwh": consumption,
            "components": {"strompreiskompensation": consumption} if consumption else {},
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
    n, result = _solve_cheapest_regime(n, "test", {}, GRID, tariffs, prices, waccs)
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


def test_the_winner_is_returned_when_it_was_not_solved_last():
    """The cheapest network so far is kept as a copy, so the winner comes back
    without solving it again, whatever ran after it."""
    charges = _charges("DEU",
                       _regime("discount", 5_000, 1.0, fixed=1_000.0, h=(6000, None)),
                       _regime("full", 50_000, 10.0))
    n, result = _solve(origin=charges)
    totals = {c["regimes"]: c["total_eur_per_year"] for c in result["candidates"]}
    wacc = TEMPLATE["finance"]["origin_wacc"]

    assert result["meters"]["grid_import"]["regime"] == "discount"
    assert totals["grid_import=discount"] < totals["grid_import=full"]
    # The network is the winner's: its capital cost and its objective, not those
    # of "full", which was solved after it.
    assert float(n.generators.at["grid_import", "capital_cost"]) == pytest.approx(
        _annuity(wacc) * GRID["connection_capex_eur_per_mw"] + 5_000)
    assert float(n.objective) + 1_000.0 == pytest.approx(totals["grid_import=discount"])
    assert result["meters"]["grid_import"]["utilisation_hours"] >= 6000 * (1 - 1e-6)


# -- the connection: a part per MW in the LP, a fixed sum as a constant


CONNECTION = {"capex_eur_per_mw": 30_000.0, "capex_eur": 2_000_000.0, "components": {
    "netzanschlusskosten": {"capex_eur_per_mw": 20_000.0, "capex_eur": 2_000_000.0},
    "baukostenzuschuss": {"capex_eur_per_mw": 10_000.0, "capex_eur": 0.0}}}


def _annuity(wacc):
    return _grid_tariff._connection_annuity(GRID, wacc)


def test_the_workbook_connection_replaces_the_configs():
    charges = _charges("DEU", _regime("full", 50_000, 10.0), connection=CONNECTION)
    n, result = _solve(origin=charges)
    wacc = TEMPLATE["finance"]["origin_wacc"]

    # The €/MW part is in the capital cost the LP sized the connection with.
    assert float(n.generators.at["grid_import", "capital_cost"]) == pytest.approx(
        _annuity(wacc) * 30_000.0 + 50_000)
    # The fixed sum is not: a constant per year, as a grid run always builds a
    # connection, added to the regime's cost for the comparison.
    connection = result["meters"]["grid_import"]["connection"]
    assert connection["lump_sum_eur_per_year"] == pytest.approx(_annuity(wacc) * 2_000_000.0)
    assert result["candidates"][0]["total_eur_per_year"] == pytest.approx(
        float(n.objective) + connection["lump_sum_eur_per_year"])


def test_what_the_lp_leaves_out_is_the_regimes_fixed_costs_and_the_lump_sum():
    """Metering, the tier offsets and the connection's lump sum: constants of
    the regime, as every grid run builds a connection. An area without a
    workbook has none."""
    regime = _regime("full", 50_000, 10.0, fixed=30_000.0)
    wacc = TEMPLATE["finance"]["origin_wacc"]
    assert _grid_tariff.outside_lp_eur_per_year(regime, CONNECTION, GRID, wacc) == pytest.approx(
        30_000.0 + _annuity(wacc) * 2_000_000.0)
    config = _grid_tariff.config_regime(GRID)
    assert _grid_tariff.outside_lp_eur_per_year(config, _grid_tariff.connection_of(None, GRID), GRID, wacc) == 0.0
