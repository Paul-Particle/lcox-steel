"""A NaN or inf that would reach the optimisation stops the build.

PyPSA replaces a NaN with the attribute's default when a component is added (a NaN
capital cost makes the asset free, a NaN efficiency becomes 1), and an inf ends the
solve with status ok, condition unknown and every capacity at zero.
"""

import copy

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import REPO_ROOT
from scripts.solve.build_network import build_network

TEMPLATE = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
INDEX = pd.date_range("2025-03-01", periods=24, freq="h")
CF = pd.DataFrame({"solar": np.linspace(0.0, 1.0, 24)}, index=INDEX)
PRICE = pd.Series(80.0, index=INDEX)


def build_with(route: str, section: str, key: str, value: float):
    """The route's network on the template with one value replaced."""
    assumptions = copy.deepcopy(TEMPLATE)
    assumptions[section][key] = value
    network = build_network(route, assumptions, CF, price_series=PRICE)
    return network


def test_a_nan_stops_the_build_where_the_component_is_added():
    with pytest.raises(ValueError, match=r"Link electrolyser: NaN in \['capital_cost'\]"):
        build_with("h2-dri-eaf", "electrolyser", "capex_eur_per_mw", np.nan)
    with pytest.raises(ValueError, match=r"Link dri-h2: NaN in \['efficiency3'\]"):
        build_with("h2-dri-eaf", "dri-h2", "h2_preheat_el_mwh_per_t", np.nan)


def test_an_inf_stops_the_build_once_the_network_stands():
    with pytest.raises(ValueError, match=r"Link capital_cost.*electrolyser"):
        build_with("h2-dri-eaf", "electrolyser", "capex_eur_per_mw", np.inf)


def test_a_nan_the_route_does_not_read_is_no_problem():
    """A key no workbook scenario requires is NaN, and only the runs that read it stop."""
    build_with("moe-eaf", "dri-h2", "h2_preheat_el_mwh_per_t", np.nan)
