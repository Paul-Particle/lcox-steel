"""Thematic helpers shared across WIP res_cf analysis scripts.

Lives next to its consumers (determine_bestsite_p95, determine_complementarity,
determine_resource_spread, diag_plot_bestsite_p95) rather than in a top-level
common/ module — these helpers are specific to the quarterly-cutout workflow
used by those scripts and are not used by the active Snakemake pipeline
(build_regions, build_offshore_regions, download_cutout, build_res_cf_profile).
"""

from pathlib import Path

import numpy as np
import yaml

from common._paths import CUTOUTS, REPO_ROOT

QUARTERS = ["q1", "q2", "q3", "q4"]


def load_res_cf_cfg() -> dict:
    """Read the `res_cf` config block from config/config.yaml (standalone-mode default).

    Snakemake-driven runs receive this via snakemake.config instead.
    """
    # --- Previous docstring (kept for reference) below ---
    # Read the `res_cf` block from config/config.yaml. Standalone-mode default;
    # Snakemake-driven runs go through snakemake.config instead.
    with open(REPO_ROOT / "config/config.yaml") as f:
        return yaml.safe_load(f)["res_cf"]


def cutout_path(country: str, year: int, quarter: str) -> Path:
    """Return the canonical path to an atlite cutout for (country, year, quarter)."""
    # --- Previous docstring (kept for reference) below ---
    # Canonical path to an atlite cutout.
    return CUTOUTS / f"{country.lower()}_{year}_{quarter}.nc"


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
