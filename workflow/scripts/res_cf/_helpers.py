"""Thematic helpers shared across res_cf scripts.

Lives next to its consumers (03_build_cf_timeseries, 07_make_bestsite_cf_timeseries,
07b_make_anchor_colocated_cf_timeseries) rather than in a top-level common/ module.
The eligibility helpers below are used by the active Snakemake pipeline (03/07/07b);
the config/haversine helpers are also used by standalone analysis scripts.
"""

from pathlib import Path

import numpy as np
import yaml

from common._paths import CUTOUTS, REPO_ROOT


# ── Land-sea eligibility (issue #41) ─────────────────────────────────────────
# Coarse ERA5 (~27-31 km) coastal cells blend smooth-sea and rough-land wind, so
# onshore wind CF is inflated wherever a cell contains sea. atlite's
# `indicatormatrix` gives each cell's land-coverage fraction (area overlap with
# the region polygon, EPSG:4326) but a mostly-sea coastal cell still carries its
# inflated CF and can dominate P95 / best-site / max *selection*. These helpers
# drop cells whose land fraction is below a threshold (relative to the region's
# max), reproducing the validated notebook fix
# (notebooks/res_cf_coastal_cf_exploration.ipynb: full-land = w >= 0.98 * wmax).
#
# PLUGGABLE SOURCE: this is the single definition of per-cell eligibility shared
# by 03/07/07b. The land-fraction *source* is selectable via the config value
# `res_cf.eligibility_source` ("indicatormatrix" | "availabilitymatrix"):
#   - "indicatormatrix" (default): coarse land/sea outline, no extra data. This
#     is what the notebook validated and what ships enabled today.
#   - "availabilitymatrix": a finer land-sea fraction from atlite's
#     `availabilitymatrix(shapes, ExclusionContainer)` over a high-res land polygon.
#     Sharper coastline than NE-110m. Wired (see `cell_availability_fraction`) and
#     cross-checked against indicatormatrix — agrees in 5/6 areas. Needs the
#     ne_10m_land download; consumers are identical either way.
# See docs/land_eligibility_design.md.
#
# ONSHORE ONLY: callers pass min_land_fraction=0 for offshore wind (which is
# sited on sea by design), making every helper a no-op there.

ELIGIBILITY_SOURCES = ("indicatormatrix", "availabilitymatrix")


def cell_land_fraction(cutout, geom) -> np.ndarray:
    """Flat (n_cells,) land-coverage fraction for one geometry.

    From atlite's `indicatormatrix` — the area fraction of each grid cell that
    overlaps `geom` (EPSG:4326). Cells outside the region are 0.
    """
    ind = cutout.indicatormatrix([geom]).tocsr()
    return np.asarray(ind[0, :].todense()).ravel()


_LAND_GEOM_CACHE: dict = {}


def _load_land_geoms(path):
    """Cached high-res land polygons (GeoDataFrame, EPSG:4326).

    Read once per path — 07/07b call the availability source several times
    (per tech), and re-reading the shapefile each time would be wasteful.
    """
    key = str(path)
    if key not in _LAND_GEOM_CACHE:
        import geopandas as gpd

        _LAND_GEOM_CACHE[key] = gpd.read_file(path).to_crs(4326)
    return _LAND_GEOM_CACHE[key]


def cell_availability_fraction(cutout, geom) -> np.ndarray:
    """Flat (n_cells,) land fraction from a finer coastline (availabilitymatrix).

    The `eligibility_source='availabilitymatrix'` path. Excludes sea using a
    high-res land polygon (config `res_cf.availability.land_shapes`, e.g. Natural
    Earth 10 m land) rasterised at `res` m in the equal-area `crs`, then returns
    atlite's per-cell available fraction — 'fraction of cell that is land', the
    finer-coastline analogue of `cell_land_fraction` (which uses the coarse
    NE-110m outline). Same (y, x) C-order as `indicatormatrix`, so it drops into
    the same downstream slot unchanged. Onshore only (offshore passes
    min_land_fraction=0, so this is never invoked for offshore).
    """
    import geopandas as gpd
    from shapely.geometry import box
    from atlite.gis import ExclusionContainer

    cfg = load_res_cf_cfg().get("availability", {})
    land_path = cfg.get("land_shapes", "data/shapes/ne_10m_land/ne_10m_land.zip")
    crs = cfg.get("crs", 6933)
    res = cfg.get("res", 200)

    p = Path(land_path)
    land = _load_land_geoms(p if p.is_absolute() else REPO_ROOT / p)

    # Clip land to the region's padded bbox so we only rasterise near the country
    # (the raw NE-10m land is whole continents — expensive to rasterise globally).
    minx, miny, maxx, maxy = geom.bounds
    pad = 1.0
    land_clip = land.clip(box(minx - pad, miny - pad, maxx + pad, maxy + pad))

    excluder = ExclusionContainer(crs=crs, res=res)
    # invert=False here yields available = land within the region (verified
    # empirically against indicatormatrix); atlite composes the geometry with the
    # region mask so the effect is "keep land, drop sea", not the literal reading
    # of the add_geometry docstring.
    excluder.add_geometry(land_clip.geometry, invert=False)

    shapes = gpd.GeoSeries([geom], crs=4326)
    avail = cutout.availabilitymatrix(shapes, excluder, disable_progressbar=True)
    return np.asarray(avail.values[0]).ravel()


def cell_eligible_fraction(cutout, geom, source: str = "indicatormatrix") -> np.ndarray:
    """Flat (n_cells,) eligible fraction from the configured `source`.

    Dispatches to `cell_land_fraction` (indicatormatrix, default) or
    `cell_availability_fraction` (availabilitymatrix). The one place the
    eligibility ingredient is chosen; the threshold and consumers are unaffected.
    """
    if source == "indicatormatrix":
        return cell_land_fraction(cutout, geom)
    if source == "availabilitymatrix":
        return cell_availability_fraction(cutout, geom)
    raise ValueError(
        f"Unknown eligibility_source {source!r}; expected one of {ELIGIBILITY_SOURCES}"
    )


def eligibility_mask(land_fraction: np.ndarray, min_land_fraction: float) -> np.ndarray:
    """Boolean keep-mask: cells with land fraction >= min_land_fraction * max.

    Threshold is relative to the region's maximum land fraction (`wmax`), matching
    the validated notebook method and staying robust for small regions with no
    fully-interior cell (there `wmax < 1`). The max-land cell is always kept, so
    at least one cell survives. Returns all-True when `min_land_fraction <= 0`.
    """
    if min_land_fraction <= 0:
        return np.ones_like(land_fraction, dtype=bool)
    wmax = land_fraction.max()
    if wmax <= 0:
        return np.zeros_like(land_fraction, dtype=bool)
    return land_fraction >= min_land_fraction * wmax


def eligibility_fraction_1d(cutout, geom, min_land_fraction: float = 0.0,
                            source: str = "indicatormatrix") -> np.ndarray:
    """Flat (n_cells,) eligible fraction with sub-threshold cells zeroed.

    The core of the fix: take the per-cell eligible fraction from `source` and
    zero every cell below the land-sea cutoff. Shared by `eligibility_weights`
    (07/07b) and `eligibility_matrix` (03) so all consumers apply one rule.
    """
    lf = cell_eligible_fraction(cutout, geom, source)
    return np.where(eligibility_mask(lf, min_land_fraction), lf, 0.0)


def eligibility_weights(cutout, geom, min_land_fraction: float = 0.0,
                        source: str = "indicatormatrix") -> np.ndarray:
    """Per-cell eligibility weight as a (y, x) array — 07/07b's P95/max/mean
    weights. Cells outside the region are already 0, so they stay excluded.
    """
    frac = eligibility_fraction_1d(cutout, geom, min_land_fraction, source)
    n_y = cutout.data.sizes["y"]
    n_x = cutout.data.sizes["x"]
    return frac.reshape(n_y, n_x)


def eligibility_matrix(cutout, geom, min_land_fraction: float = 0.0,
                       source: str = "indicatormatrix"):
    """(1, n_cells) sparse weight matrix for atlite's per_unit aggregation (03's
    country mean). Row values are the thresholded eligible fraction per cell, so
    `cutout.wind(matrix=..., per_unit=True)` averages over eligible land only.
    """
    from scipy.sparse import csr_matrix

    frac = eligibility_fraction_1d(cutout, geom, min_land_fraction, source)
    return csr_matrix(frac.reshape(1, -1))


def eligibility_mask_2d(cutout, geom, min_land_fraction: float = 0.0,
                        source: str = "indicatormatrix") -> np.ndarray:
    """Boolean (y, x) keep-mask for `geom` — the 2-D form of `eligibility_mask`,
    for AND-ing into cell-validity checks (e.g. 07b candidate selection).
    """
    lf = cell_eligible_fraction(cutout, geom, source)
    n_y = cutout.data.sizes["y"]
    n_x = cutout.data.sizes["x"]
    return eligibility_mask(lf, min_land_fraction).reshape(n_y, n_x)


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
