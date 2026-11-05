"""A grid-only run (power_supply `grid`) solves and reports: no renewables, the grid
import with the battery and, on a hydrogen route, the H2 buffer beside it."""

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import REPO_ROOT
from scripts.solve.build_network import build_network
from scripts.viz import compile_report

# As solve_network sets it: xarray, under PyPSA's optimise, cannot index pandas 3's Arrow strings.
pd.options.mode.string_storage = "python"
TEMPLATE = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
INDEX = pd.date_range("2025-03-01", periods=48, freq="h")
PRICE = pd.Series(np.where(INDEX.hour < 12, 30.0, 120.0), index=INDEX)
MIX = pd.DataFrame({"gas": 400.0, "hard_coal": 100.0, "solar": np.where(INDEX.hour < 12, 500.0, 0.0)}, index=INDEX)


@pytest.mark.parametrize("route", ["h2-dri-eaf", "moe-eaf-export"])
def test_a_grid_only_run_solves_and_reports_without_renewables(route):
    is_export = route.endswith("-export")
    legs = {"sea": 17_965} if is_export else None
    no_renewables = pd.DataFrame(index=INDEX)
    n = build_network(route, TEMPLATE, no_renewables, price_series=PRICE, transport_legs=legs,
                      destination_price=PRICE * 1.1 if is_export else None)

    status, condition = n.optimize(solver_name="highs")
    assert (status, condition) == ("ok", "optimal")
    assert not n.generators.carrier.isin(TEMPLATE["res"]).any()
    assert "battery" in n.storage_units.index or "battery" in n.stores.index
    assert ("h2_buffer" in n.stores.index) == route.startswith("h2-")

    run = {"run": "AUS-ESP" if is_export else "ESP", "origin": "AUS" if is_export else "ESP",
           "destination": "ESP", "route": route, "start_date": "20250301", "end_date": "20250302"}
    summary = compile_report.extract_summary(n, "grid-only", run, TEMPLATE, legs, MIX,
                                             MIX if is_export else None)
    assert np.isfinite(summary["lcos_eur_per_t_steel"]) and summary["lcos_eur_per_t_steel"] > 0
    assert summary["emissions_energy_t_co2_per_t_steel"] > 0
