"""Unit tests for which freight a run pays and where its furnace is financed.

A run with a destination elsewhere ships: an `-export` route its iron, to a
furnace financed at the destination's WACC, any other route its steel, from a
furnace financed at the origin's.
"""

import numpy as np
import pandas as pd
import pytest
import yaml

from scripts.solve.build_network import build_network
from conftest import REPO_ROOT

HOURS = 24
LEGS = {"sea": 17_965}


@pytest.fixture(scope="module")
def assumptions() -> dict:
    """The template's assumptions, with the destination dearer to borrow in."""
    template = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    template["finance"]["destination_wacc"] = template["finance"]["origin_wacc"] + 0.05
    return template


def _built(route: str, assumptions: dict, legs: dict | None):
    index = pd.date_range("2025-01-01", periods=HOURS, freq="h")
    cf = pd.DataFrame({"solar": np.full(HOURS, 0.3)}, index=index)
    price = pd.Series(np.full(HOURS, 50.0), index=index)
    network = build_network(
        route, assumptions, cf, price_series=price, transport_legs=legs,
        destination_price=price if route.endswith("-export") else None,
    )
    return network


def test_a_run_without_a_destination_elsewhere_stops_at_the_plant_gate(assumptions):
    n = _built("moe-eaf", assumptions, legs=None)

    assert "steel_transport" not in n.links.index
    assert n.loads.at["steel_load", "bus"] == "steel"


def test_a_domestic_route_with_a_destination_ships_its_steel(assumptions):
    home = _built("moe-eaf", assumptions, legs=None)
    shipped = _built("moe-eaf", assumptions, legs=LEGS)

    assert shipped.links.at["steel_transport", "length"] == 17_965
    assert "iron_transport" not in shipped.links.index
    assert shipped.loads.at["steel_load", "bus"] != "steel"
    # The furnace stays at home, financed as it is there.
    assert shipped.links.at["eaf", "capital_cost"] == home.links.at["eaf", "capital_cost"]


def test_an_export_ships_its_iron_to_a_furnace_financed_at_the_destination(assumptions):
    home = _built("moe-eaf", assumptions, legs=None)
    exported = _built("moe-eaf-export", assumptions, legs=LEGS)

    assert "iron_transport" in exported.links.index
    assert "steel_transport" not in exported.links.index
    assert exported.stores.at["steel_store", "capital_cost"] > home.stores.at["steel_store", "capital_cost"]
