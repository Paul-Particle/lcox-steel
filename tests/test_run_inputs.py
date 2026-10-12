"""Unit tests for how a run's rows say what each input series is for.

A run's first area is its origin and its second its destination; the origin's
grid series prices its own imports and the destination's the furnace an export
melts in.
"""

import pandas as pd
import pytest

from scripts.solve._helpers_solve import find_transport_legs, read_run_inputs

ROWS = pd.DataFrame(
    [("s", "DEU", "h2-dri-eaf", "solar", "bestsite-p95", "DEU"),
     ("s", "DEU", "h2-dri-eaf", "grid", "emissions", "DEU"),
     ("s", "AUS", "h2-dri-eaf-export", "solar", "bestsite-p95", "AUS"),
     ("s", "AUS", "h2-dri-eaf-export", "grid", "emissions", "DEU"),
     ("s", "NSW1", "h2-dri-eaf-export", "solar", "bestsite-p95", "NSW1"),
     ("s", "NSW1", "h2-dri-eaf-export", "grid", "emissions", "NSW1"),
     ("s", "NSW1", "h2-dri-eaf-export", "grid", "emissions", "DEU"),
     ("s", "BAD", "moe-eaf-export", "solar", "bestsite-p95", "AUS"),
     ("s", "BAD", "moe-eaf-export", "grid", "emissions", "DEU"),
     ("s", "BAD", "moe-eaf-export", "grid", "emissions", "FRA")],
    columns=["scenario", "run", "route", "tech", "variant", "area"],
).assign(start_date="20250101", end_date="20251231")


def _inputs(tmp_path, run):
    table = tmp_path / "scenarios.csv"
    ROWS.to_csv(table, index=False)
    names = (ROWS["area"] + "_" + ROWS["tech"] + "_" + ROWS["variant"] + "_20250101_20251231.parquet")
    run_inputs = read_run_inputs(table, "s", run, [f"resources/timeseries/{name}" for name in names])
    return run_inputs


def test_a_domestic_run_is_its_own_destination(tmp_path):
    run_inputs = _inputs(tmp_path, "DEU")

    assert (run_inputs["origin"], run_inputs["destination"]) == ("DEU", "DEU")
    assert run_inputs["grid_path"].name == "DEU_grid_emissions_20250101_20251231.parquet"
    assert run_inputs["destination_path"] is None


def test_an_islanded_export_prices_only_its_furnace(tmp_path):
    run_inputs = _inputs(tmp_path, "AUS")

    assert (run_inputs["origin"], run_inputs["destination"]) == ("AUS", "DEU")
    assert run_inputs["route"] == "h2-dri-eaf-export"
    assert [path.name for path in run_inputs["cf_paths"]] == ["AUS_solar_bestsite-p95_20250101_20251231.parquet"]
    assert run_inputs["grid_path"] is None
    assert run_inputs["destination_path"].name == "DEU_grid_emissions_20250101_20251231.parquet"


def test_a_grid_export_keeps_the_two_markets_apart(tmp_path):
    run_inputs = _inputs(tmp_path, "NSW1")

    assert run_inputs["grid_path"].name == "NSW1_grid_emissions_20250101_20251231.parquet"
    assert run_inputs["destination_path"].name == "DEU_grid_emissions_20250101_20251231.parquet"


def test_a_third_area_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="at most one destination"):
        _inputs(tmp_path, "BAD")


DISTANCES = {"default": {"sea": 3141.593}, "AUS": {"DEU": {"sea": 21_947}}}
AREAS = {"AUS": {"zones": ["NSW1", "VIC1"]}, "NSW1": {"market": "nem"}, "DEU": {"market": "entsoe"}}


def test_a_zone_ships_as_its_country():
    assert find_transport_legs(DISTANCES, "NSW1", "DEU", AREAS) == {"sea": 21_947}


def test_a_pair_not_listed_gets_the_default():
    assert find_transport_legs(DISTANCES, "AUS", "FRA", AREAS) == {"sea": 3141.593}
    assert find_transport_legs(DISTANCES, "BRA", "DEU", AREAS) == {"sea": 3141.593}
