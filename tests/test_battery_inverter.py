"""The battery's charger and discharger are one inverter, rated and paid once.

`build_network` prices only the charger, so without the tie the solver sizes
discharge power for nothing and the battery comes out cheaper than it is. These
build the smallest network that wants the two ratings to differ — slow charge,
fast discharge — and check the constraint is what stops it, with the rating
measured on the grid side in both directions.
"""

import pypsa
import pytest

from scripts.solve.solve_network import _tie_battery_inverter

CHARGER_CAPEX = 1.0
ETA = 0.9
PEAK_DEMAND_MW = 2.0


def _battery_network() -> pypsa.Network:
    """Two MW of supply over two hours, two MW of demand in the third."""
    n = pypsa.Network()
    n.set_snapshots(range(3))
    n.add("Bus", "electricity", carrier="AC")
    n.add("Bus", "electricity_battery", carrier="battery")
    n.add("Generator", "supply", bus="electricity", carrier="AC",
          p_nom=2.0, marginal_cost=0.0, p_max_pu=[1.0, 1.0, 0.0])
    n.add("Load", "demand", bus="electricity", p_set=[0.0, 0.0, PEAK_DEMAND_MW])
    n.add("Store", "battery", bus="electricity_battery", carrier="battery",
          e_nom_extendable=True, capital_cost=0.0)
    n.add("Link", "battery_charger", bus0="electricity", bus1="electricity_battery",
          carrier="battery", efficiency=ETA,
          p_nom_extendable=True, capital_cost=CHARGER_CAPEX)
    n.add("Link", "battery_discharger", bus0="electricity_battery", bus1="electricity",
          carrier="battery", efficiency=ETA,
          p_nom_extendable=True, capital_cost=0.0)
    return n


def test_the_two_ratings_would_differ_without_the_tie():
    """The premise: the demand really does pull the two ratings apart."""
    n = _battery_network()
    n.optimize(solver_name="highs")

    charge_over_two_hours = PEAK_DEMAND_MW / ETA**2 / 2
    assert n.links.at["battery_charger", "p_nom_opt"] == pytest.approx(charge_over_two_hours)
    assert n.links.at["battery_discharger", "p_nom_opt"] == pytest.approx(PEAK_DEMAND_MW / ETA)


def test_the_tie_rates_both_directions_at_the_grid():
    n = _battery_network()
    n.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)

    charger_grid_mw = n.links.at["battery_charger", "p_nom_opt"]
    discharger_grid_mw = ETA * n.links.at["battery_discharger", "p_nom_opt"]
    assert charger_grid_mw == pytest.approx(discharger_grid_mw)
    # Sized for the harder of the two directions, and paid for once.
    assert charger_grid_mw == pytest.approx(PEAK_DEMAND_MW)


SUPPLY_PRICE_EUR_PER_MWH = -10.0
INVERTER_MAX_MW = 1.0


def _battery_paid_to_waste_network() -> pypsa.Network:
    """Power that pays to be taken and nowhere to put it but the battery's losses.

    Charging and discharging at once burns (1 - ETA**2) of what goes round, and
    every MWh burnt is a MWh more bought at a negative price.
    """
    n = pypsa.Network()
    n.set_snapshots(range(2))
    n.add("Bus", "electricity", carrier="AC")
    n.add("Bus", "electricity_battery", carrier="battery")
    n.add("Generator", "supply", bus="electricity", carrier="AC",
          p_nom=10.0, marginal_cost=SUPPLY_PRICE_EUR_PER_MWH)
    n.add("Store", "battery", bus="electricity_battery", carrier="battery",
          e_nom_extendable=True, e_cyclic=True, capital_cost=0.0)
    n.add("Link", "battery_charger", bus0="electricity", bus1="electricity_battery",
          carrier="battery", efficiency=ETA, p_nom_extendable=True,
          p_nom_max=INVERTER_MAX_MW, capital_cost=CHARGER_CAPEX)
    n.add("Link", "battery_discharger", bus0="electricity_battery", bus1="electricity",
          carrier="battery", efficiency=ETA, p_nom_extendable=True,
          p_nom_max=INVERTER_MAX_MW / ETA, capital_cost=0.0)
    return n


def _grid_side_flows(n: pypsa.Network) -> tuple:
    charge = n.links_t.p0["battery_charger"]
    discharge = ETA * n.links_t.p0["battery_discharger"]
    return charge, discharge


def test_without_the_shared_rating_the_battery_charges_and_discharges_at_full_power_at_once():
    """The premise: paid to waste power, the two links each run at the full rating."""
    n = _battery_paid_to_waste_network()
    n.optimize(solver_name="highs")

    charge, discharge = _grid_side_flows(n)
    assert (charge + discharge).max() == pytest.approx(2 * INVERTER_MAX_MW)


def test_charging_and_discharging_together_stay_within_the_one_inverter_rating():
    n = _battery_paid_to_waste_network()
    n.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)

    charge, discharge = _grid_side_flows(n)
    rating_mw = n.links.at["battery_charger", "p_nom_opt"]
    assert rating_mw == pytest.approx(INVERTER_MAX_MW)
    assert (charge + discharge).max() == pytest.approx(rating_mw)
