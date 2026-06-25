"""Thematic helpers shared across WIP res_cf analysis scripts.

Lives next to its consumers (07_make_bestsite_cf_timeseries, 08_complementarity_screen,
06_resource_spread, 100_plot_bestsite_locations) rather than in a top-level
common/ module. Not used by the active Snakemake pipeline
(01_build_regions, 01b_build_offshore_regions, 02_make_cutouts, 03_build_cf_timeseries).
"""

from pathlib import Path

import numpy as np
import yaml

from common._paths import CUTOUTS, REPO_ROOT


def load_res_cf_cfg() -> dict:
    """Read the `res_cf` config block from config/config.yaml (standalone-mode default).

    Snakemake-driven runs receive this via snakemake.config instead.
    """
    with open(REPO_ROOT / "config/config.yaml") as f:
        return yaml.safe_load(f)["res_cf"]


def annual_cutout_path(cf_area: str, year: int) -> Path:
    """Return the path to the single annual atlite cutout for (cf_area, year).

    Matches the output pattern of the `download_cutout` rule:
    cutouts/{cf_area}_{year}0101_{year}1231.nc
    """
    return CUTOUTS / f"{cf_area.lower()}_{year}0101_{year}1231.nc"


def haversine_distance_km(
    lon1: float,
    lat1: float,
    lon2: np.ndarray,
    lat2: np.ndarray,
) -> np.ndarray:
    """Great-circle distance (km) from one target point to arrays of points.

    `lon1`/`lat1` are the target point and `lon2`/`lat2` the candidate-point
    arrays, all in degrees.
    """
    # --- Previous docstring (kept for reference) below ---
    # Great-circle distance (km) between one target point and arrays of points.
    #
    # lon1, lat1: target point in degrees.
    # lon2, lat2: arrays of candidate point coordinates in degrees.
    earth_radius_km = 6371.0

    lon1_rad = np.deg2rad(lon1)
    lat1_rad = np.deg2rad(lat1)
    lon2_rad = np.deg2rad(lon2)
    lat2_rad = np.deg2rad(lat2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2.0) ** 2
    )
    c = 2.0 * np.arcsin(np.sqrt(a))

    return earth_radius_km * c
