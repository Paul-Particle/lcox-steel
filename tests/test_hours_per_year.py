"""A full calendar year levelises over its own hours, leap year or not, and a solve
sizes the plant the same whatever length of window stands for that year."""

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import REPO_ROOT
from scripts.solve._helpers_solve import mean_year_hours
from scripts.solve.build_network import build_network

# linopy cannot index pandas' Arrow strings; solve_network sets the same.
pd.options.mode.string_storage = "python"


@pytest.mark.parametrize("year", [2024, 2025])
def test_a_full_year_is_not_rescaled(year):
    snapshots = pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")
    assert mean_year_hours(snapshots) == len(snapshots)


def test_a_part_year_scales_to_its_calendar_year():
    january_2024 = pd.date_range("2024-01-01", "2024-01-31 23:00", freq="h")
    assert mean_year_hours(january_2024) == 8784.0


def _solved_solar_mw(days: int) -> float:
    """Solar built for a grid-connected MOE plant over `days` repeats of one sunny day."""
    assumptions = yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())
    index = pd.date_range("2025-03-01", periods=24 * days, freq="h")
    daylight = np.clip(np.sin(np.linspace(0, 2 * np.pi, 24, endpoint=False) - np.pi / 2), 0, None)
    cf = pd.DataFrame({"solar": np.tile(daylight, days)}, index=index)
    price = pd.Series(80.0, index=index)
    n = build_network("moe-eaf", assumptions, cf, price_series=price)
    n.optimize(solver_name="highs")
    solar_mw = float(n.generators.at["solar", "p_nom_opt"])
    return solar_mw


def test_a_week_and_a_fortnight_of_the_same_days_size_the_same_plant():
    """The objective counts a window's hours as one year against the annuities.
    Unweighted, a week's grid bill is too small to be worth any solar at all."""
    week_mw, fortnight_mw = _solved_solar_mw(7), _solved_solar_mw(14)
    assert week_mw > 0
    assert fortnight_mw == pytest.approx(week_mw, rel=1e-4)
