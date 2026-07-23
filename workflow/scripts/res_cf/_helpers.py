"""Thematic helpers shared across res_cf scripts.

Consumed by active Snakemake-pipeline scripts (03_build_cf_timeseries,
07b_make_anchor_colocated_cf_timeseries) as well as the WIP analysis script
07_make_bestsite_cf_timeseries. Lives next to its consumers rather than in the
top-level common/ package.
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


def cos_lat_weights(cutout) -> np.ndarray:
    """Per-cell cos(latitude) weights, flat in cutout-grid order.

    atlite's indicatormatrix returns land-coverage fractions computed in
    EPSG:4326, so every cell counts equally regardless of latitude. A lat/lon
    cell's physical area is proportional to cos(lat), so multiplying the
    fraction weights by these values turns a degree-area average into a
    physical-area average (issue #37). Unlike a region-specific equal-area CRS
    (e.g. EPSG:3035, Europe-only), cos(lat) is valid for every cf_area.

    The ordering matches cutout.grid rows, which is also the column order of
    cutout.indicatormatrix(...) and the ravel of a (y, x) weight grid — so the
    result can scale either directly.
    """
    lats = cutout.grid["y"].to_numpy()
    return np.cos(np.deg2rad(lats))


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
