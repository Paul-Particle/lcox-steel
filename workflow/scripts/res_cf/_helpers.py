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


def weighted_percentile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Weighted percentile of `values` (nonnegative `weights`), q in [0, 1].

    Entries with a non-finite value, non-finite weight, or zero weight are
    dropped. Returns NaN if nothing survives.
    """
    if not (0.0 <= q <= 1.0):
        raise ValueError("q must be in [0, 1].")

    m = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    v = values[m]
    w = weights[m]

    if v.size == 0:
        return np.nan

    order = np.argsort(v)
    v = v[order]
    w = w[order]

    cw = np.cumsum(w)
    cw /= cw[-1]

    idx = np.searchsorted(cw, q, side="left")
    idx = min(idx, v.size - 1)
    return float(v[idx])


def geom_area_weights(cutout, geom) -> np.ndarray:
    """Physical-area cell weights for `geom` as a (y, x) grid.

    atlite's indicatormatrix gives each cell's fractional overlap with `geom`
    (computed in EPSG:4326); scaling by cos_lat_weights turns those degree-area
    fractions into physical-area weights (issue #37). Cells outside `geom` get
    weight 0. Shared by scripts 03b, 07, 07b and viz/plot_cf_map so the whole
    pipeline uses one definition of area weighting.
    """
    indicator = cutout.indicatormatrix([geom]).tocsr()
    weights_1d = np.asarray(indicator[0, :].todense()).ravel() * cos_lat_weights(cutout)
    n_y = cutout.data.sizes["y"]
    n_x = cutout.data.sizes["x"]
    return weights_1d.reshape(n_y, n_x)


def pick_p95_cell(cell_mean, weights: np.ndarray, q: float = 0.95) -> tuple[int, int]:
    """(y, x) index of the in-region cell whose annual-mean CF is closest to the
    area-weighted q-percentile.

    `cell_mean` is a 2-D (y, x) grid of annual-mean CF (xr.DataArray or ndarray);
    `weights` is the matching (y, x) grid from geom_area_weights. Cells with
    weight 0 (outside the region) are excluded. This is the single definition of
    "the P95 cell" shared by scripts 03b, 07, 07b and viz/plot_cf_map.
    """
    vals2d = np.asarray(getattr(cell_mean, "values", cell_mean))
    w = np.asarray(weights)

    p = weighted_percentile(vals2d.ravel(), w.ravel(), q)

    valid = w.ravel() > 0
    vals = np.where(valid, vals2d.ravel(), np.nan)
    dist = np.abs(vals - p)
    idx_flat = np.nanargmin(dist)
    y_idx, x_idx = np.unravel_index(idx_flat, vals2d.shape)
    return int(y_idx), int(x_idx)


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
