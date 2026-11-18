"""The solve runs once per tariff regime and keeps the cheapest.

Built on the two-day grid-only network of test_grid_only_run: the plant buys
all its power, so its utilisation hours and annual energy are what the regimes
are cut on.
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


def _regime(name, capacity, energy, fixed=0.0, h=(None, None), e=(None, None)):
    return {"name": name, "h_min": h[0], "h_max": h[1], "e_min_mwh": e[0], "e_max_mwh": e[1],
            "capacity_eur_per_mw_per_year": capacity, "energy_eur_per_mwh": energy,
            "fixed_eur_per_year": fixed, "components": {}}


def _tariff(area, *regimes, subsidy=0.0):
    return {"area": area, "source": "test", "regimes": list(regimes),
            "post_solve": {"subsidy_eur_per_mwh": subsidy, "subsidy_components": ["s"] if subsidy else []}}


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
    return n, result


def test_an_area_without_a_workbook_is_billed_as_before():
    n, result = _solve()
    reference = build_network("h2-dri-eaf", TEMPLATE, pd.DataFrame(index=INDEX), price_series=PRICE)
    reference.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)

    bill = result["meters"]["grid_import"]
    assert bill["regime"] == "config" and len(result["candidates"]) == 1
    assert float(n.objective) == pytest.approx(float(reference.objective), rel=1e-6)


def test_the_cheapest_regime_is_kept_and_its_box_holds():
    # The discount needs ≥ 6000 h; the plant is flat enough to reach it.
    tariff = _tariff("DEU",
                     _regime("full", 50_000, 10.0),
                     _regime("discount", 5_000, 1.0, fixed=1_000.0, h=(6000, None)))
    n, result = _solve(origin=tariff)
    bill = result["meters"]["grid_import"]

    assert bill["regime"] == "discount"
    assert bill["utilisation_hours"] >= 6000 * (1 - 1e-6)
    totals = {c["regimes"]: c["total_eur_per_year"] for c in result["candidates"]}
    assert totals["grid_import=discount"] < totals["grid_import=full"]
    # The fixed €/yr counts in the comparison, not in the objective.
    assert totals["grid_import=discount"] == pytest.approx(float(n.objective) + 1_000.0)


def test_a_regime_the_plant_cannot_reach_is_skipped():
    tariff = _tariff("DEU",
                     _regime("tiny", 0.0, 0.0, e=(None, 1.0)),
                     _regime("full", 50_000, 10.0))
    _, result = _solve(origin=tariff)

    assert result["meters"]["grid_import"]["regime"] == "full"
    assert [c["total_eur_per_year"] is None for c in result["candidates"]] == [True, False]


def test_the_annual_energy_counts_a_short_window_as_a_year():
    n, result = _solve()
    bill = result["meters"]["grid_import"]
    two_days = float(n.generators_t.p["grid_import"].sum())
    assert bill["energy_mwh_per_year"] == pytest.approx(two_days * 8760 / 48)


def test_the_subsidy_stays_out_of_the_objective():
    plain, _ = _solve(origin=_tariff("DEU", _regime("full", 50_000, 10.0)))
    subsidised, result = _solve(origin=_tariff("DEU", _regime("full", 50_000, 10.0), subsidy=-30.0))
    bill = result["meters"]["grid_import"]

    assert float(subsidised.objective) == pytest.approx(float(plain.objective), rel=1e-9)
    assert bill["subsidy_eur_per_year"] == pytest.approx(-30.0 * bill["energy_mwh_per_year"])


def test_an_export_tries_every_pair_of_regimes():
    origin = _tariff("AUS", _regime("a1", 50_000, 10.0), _regime("a2", 60_000, 5.0))
    destination = _tariff("DEU", _regime("d1", 50_000, 10.0), _regime("d2", 40_000, 20.0),
                          _regime("d3", 45_000, 15.0))
    _, result = _solve("moe-eaf-export", origin, destination)

    assert len(result["candidates"]) == 6
    assert set(result["meters"]) == {"grid_import", "destination_supply"}
    best = min(c["total_eur_per_year"] for c in result["candidates"])
    chosen = f"grid_import={result['meters']['grid_import']['regime']}, " \
             f"destination_supply={result['meters']['destination_supply']['regime']}"
    assert {c["regimes"]: c["total_eur_per_year"] for c in result["candidates"]}[chosen] == best
