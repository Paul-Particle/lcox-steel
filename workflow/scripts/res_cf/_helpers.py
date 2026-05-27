"""
Thematic helpers shared across res_cf pipeline scripts.

Lives next to its consumers (make_cutout, determine_bestsite_p95, determine_complementarity,
resource_spread, diag_plot_bestsite_locations) rather than in a top-level
common/ module — the helpers are specific to the atlite/quarterly-cutout
workflow.
"""

from pathlib import Path

import numpy as np
import yaml

from common._paths import CUTOUTS, REPO_ROOT

QUARTERS = ["q1", "q2", "q3", "q4"]


def load_res_cf_cfg() -> dict:
    """Read the `res_cf` block from config/config.yaml. Standalone-mode default;
    Snakemake-driven runs go through snakemake.config instead."""
    with open(REPO_ROOT / "config/config.yaml") as f:
        return yaml.safe_load(f)["res_cf"]

QUARTER_DATES = {
    "q1": ("-01-01", "-03-31 23:00"),
    "q2": ("-04-01", "-06-30 23:00"),
    "q3": ("-07-01", "-09-30 23:00"),
    "q4": ("-10-01", "-12-31 23:00"),
}


def cutout_path(country: str, year: int, quarter: str) -> Path:
    """Canonical path to an atlite cutout."""
    return CUTOUTS / f"{country.lower()}_{year}_{quarter}.nc"


def haversine_distance_km(
    lon1: float,
    lat1: float,
    lon2: np.ndarray,
    lat2: np.ndarray,
) -> np.ndarray:
    """
    Great-circle distance (km) between one target point and arrays of points.

    lon1, lat1: target point in degrees.
    lon2, lat2: arrays of candidate point coordinates in degrees.
    """
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
