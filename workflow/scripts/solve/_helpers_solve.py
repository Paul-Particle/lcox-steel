"""Financial, sizing and config utilities for the steel model."""

from pathlib import Path

import numpy as np
import pandas as pd


def find_transport_legs(distances: dict, origin: str, destination: str, areas: dict) -> dict:
    """A run's freight legs, km by mode: its origin country's to its destination, else the default.

    A market zone ships as the country whose `zones` list it is in.
    """
    country = next((parent for parent, entry in areas.items() if origin in entry.get("zones", [])), origin)
    legs = distances.get(country, {}).get(destination, distances["default"])
    return legs


def read_run_inputs(scenario_table: Path, scenario: str, run: str, timeseries: list) -> dict:
    """One run's rows of its scenario table, and the part each of its input series plays.

    The first area in the rows is the run's origin and the second its destination;
    a domestic run names one area, which is both. The origin's grid series prices
    its own imports, and the destination's prices the furnace an export melts in.
    """
    rows = (
        pd.read_csv(scenario_table, comment="#", dtype=str)
        .query("scenario == @scenario and run == @run")
    )
    areas = list(dict.fromkeys(rows["area"]))
    if len(areas) > 2:
        raise ValueError(f"{scenario}/{run} names {len(areas)} areas {areas}: "
                         f"a run has an origin and at most one destination")
    path_by_name = {Path(raw).name: Path(raw) for raw in timeseries}
    rows = rows.assign(path=(
        rows["area"] + "_" + rows["tech"] + "_" + rows["variant"] + "_"
        + rows["start_date"] + "_" + rows["end_date"] + ".parquet"
    ).map(path_by_name))
    origin, destination = areas[0], areas[-1]
    grid_rows = rows[rows["tech"] == "grid"]
    run_inputs = {
        "rows": rows,
        "origin": origin,
        "destination": destination,
        "route": rows["route"].iloc[0],
        "cf_paths": list(rows.loc[rows["tech"] != "grid", "path"]),
        "grid_path": next(iter(grid_rows.loc[grid_rows["area"] == origin, "path"]), None),
        "destination_path": (next(iter(grid_rows.loc[grid_rows["area"] == destination, "path"]), None)
                             if destination != origin else None),
    }
    return run_inputs


def deep_merge(base: dict, overlay: dict) -> dict:
    """Recursively merge `overlay` into `base` (neither input mutated).

    Overlay leaves replace base leaves; dict branches are merged key-by-key.
    """
    out = dict(base)
    for k, v in overlay.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def annuity_factor(wacc: float, lifetime_years: float) -> float:
    """Capital recovery factor: annualises a lump-sum capex over a project life."""
    if lifetime_years <= 0:
        raise ValueError("lifetime_years must be > 0")
    if wacc == 0:
        return 1.0 / lifetime_years
    return (wacc * (1.0 + wacc) ** lifetime_years) / (
        (1.0 + wacc) ** lifetime_years - 1.0
    )


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """Great-circle distance (km) between two points given in degrees.

    A scalar copy of res_cf/_helpers_res_cf.haversine_distance_km (which is array-oriented
    and lives in the res_cf package); duplicated here to keep the solve package free of a
    res_cf import for the one place it needs an inter-site distance.
    """
    earth_radius_km = 6371.0
    lon1_rad, lat1_rad = np.deg2rad(lon1), np.deg2rad(lat1)
    lon2_rad, lat2_rad = np.deg2rad(lon2), np.deg2rad(lat2)
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2.0) ** 2
    )
    return float(earth_radius_km * 2.0 * np.arcsin(np.sqrt(a)))
