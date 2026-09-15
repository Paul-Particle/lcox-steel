"""The battery's charger and discharger are one inverter, rated and paid once.

`build_network` prices only the charger, so without the tie the solver sizes
discharge power for nothing and the battery comes out cheaper than it is. These
build the smallest network that wants the two ratings to differ — slow charge,
fast discharge — and check the constraint is what stops it.
"""

import pypsa
import pytest

from solve_network import _tie_battery_inverter

CHARGER_CAPEX = 1.0


def _battery_network() -> pypsa.Network:
    """One MW of supply over two hours, two MW of demand in the third."""
    n = pypsa.Network()
    n.set_snapshots(range(3))
    n.add("Bus", "electricity", carrier="AC")
    n.add("Bus", "electricity_battery", carrier="battery")
    n.add("Generator", "supply", bus="electricity", carrier="AC",
          p_nom=1.0, marginal_cost=0.0, p_max_pu=[1.0, 1.0, 0.0])
    n.add("Load", "demand", bus="electricity", p_set=[0.0, 0.0, 2.0])
    n.add("Store", "battery", bus="electricity_battery", carrier="battery",
          e_nom_extendable=True, capital_cost=0.0)
    n.add("Link", "battery_charger", bus0="electricity", bus1="electricity_battery",
          carrier="battery", efficiency=1.0,
          p_nom_extendable=True, capital_cost=CHARGER_CAPEX)
    n.add("Link", "battery_discharger", bus0="electricity_battery", bus1="electricity",
          carrier="battery", efficiency=1.0,
          p_nom_extendable=True, capital_cost=0.0)
    return n


def test_the_two_ratings_would_differ_without_the_tie():
    """The premise: the demand really does pull the two ratings apart."""
    n = _battery_network()
    n.optimize(solver_name="highs")

    assert n.links.at["battery_charger", "p_nom_opt"] == pytest.approx(1.0)
    assert n.links.at["battery_discharger", "p_nom_opt"] == pytest.approx(2.0)


def test_the_tie_holds_the_two_ratings_together():
    n = _battery_network()
    n.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)

    charger = n.links.at["battery_charger", "p_nom_opt"]
    discharger = n.links.at["battery_discharger", "p_nom_opt"]
    assert charger == pytest.approx(discharger)
    # Sized for the harder of the two directions, and paid for once.
    assert charger == pytest.approx(2.0)
