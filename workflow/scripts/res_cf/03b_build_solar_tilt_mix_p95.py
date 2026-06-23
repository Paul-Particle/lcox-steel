"""Build N solar CF time series spanning the east–west orientation trade-off.

For a southern-hemisphere site, panels can be oriented anywhere from due west
(azimuth 270°) through north (0°, maximum annual yield) to due east (90°,
maximum morning shift). This script:

  1. Identifies the bestsite cell — P95 of annual-mean latitude_optimal CF —
     so all orientations compare at the same high-quality location. The
     site-selection helpers are copied from `07_make_bestsite_cf_timeseries.py`.
  2. Sweeps N equally spaced azimuths from 270° (west) → 0° (north) → 90° (east).
  3. For each azimuth, finds the slope that maximises annual plane-of-array (POA)
     irradiance using the Hay-Davies model and the cell's pre-computed ERA5 solar
     position data — pure numpy, no full-grid atlite call in the optimization loop.
  4. Computes the full hourly CF at each (optimal_slope, azimuth) via atlite.pv().
  5. Writes a single multi-column parquet: columns named solar_az{azimuth}.

The downstream network builder treats each column as an independent extendable
solar generator, letting the capacity optimiser choose the orientation mix.

Orientation convention (atlite / clockwise from North):
  0°   = North  (equator-facing for SH; returned by latitude_optimal for VIC)
  90°  = East
  180° = South  (equator-facing for NH)
  270° = West
"""

import logging
from pathlib import Path

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd

from common._paths import CUTOUTS, RES_CF, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_CF_AREA = "aus"
_N_STEPS = 7
_START_DATE = "20250101"
_END_DATE = "20251231"
_CUTOUT_PATH = CUTOUTS / "aus_20250101_20251231.nc"
_REGIONS_PATH = SHAPES_RES / "aus_geo.parquet"
_REGION = "AUS"
_PV_PANEL = "CSi"
_OUT = RES_CF / "aus_solar_tilt-mix-n7_20250101_20251231.parquet"

# Slope search resolution for the optimization loop (degrees).
# 1° is pure numpy so it runs in milliseconds regardless of resolution.
_SLOPE_STEP = 1

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _CF_AREA = snakemake.wildcards.cf_area
    _N_STEPS = int(snakemake.wildcards.n_steps)
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _CUTOUT_PATH = Path(snakemake.input.cutout)
    _REGIONS_PATH = Path(snakemake.input.regions)
    _REGION = snakemake.params.region
    _PV_PANEL = snakemake.params.pv_panel
    _OUT = Path(snakemake.output[0])


# ---------------------------------------------------------------------------
# Region and P95 cell helpers
# ---------------------------------------------------------------------------

def _get_region_geometry(path: Path, region: str):
    """Return the geometry of `region` from the GeoParquet at `path` (raises if missing)."""
    gdf = gpd.read_parquet(path).to_crs(4326)
    row = gdf.loc[gdf["region"] == region]
    if row.empty:
        raise ValueError(f"Region '{region}' not found in {path}")
    return row.geometry.iloc[0]


def _mask_cells_inside(cell_mean, geom) -> np.ndarray:
    """Return a boolean grid marking cutout cells whose centre lies within `geom`."""
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)


def _find_p95_cell(cf_grid, geom) -> tuple[int, int]:
    """Return the (y, x) index of the in-region cell closest to the P95 annual-mean CF."""
    # --- Previous docstring (kept for reference) below ---
    # Return (y_idx, x_idx) of the cell closest to the P95 annual-mean CF.
    cell_mean = cf_grid.mean("time")
    inside = _mask_cells_inside(cell_mean, geom)
    vals = np.where(inside, cell_mean.values, np.nan)
    valid = np.isfinite(vals)
    p95 = np.nanpercentile(vals[valid], 95)
    dist = np.abs(np.where(valid, vals, np.nan) - p95)
    idx_flat = np.nanargmin(dist)
    return tuple(int(i) for i in np.unravel_index(idx_flat, vals.shape))


# ---------------------------------------------------------------------------
# Fast analytical slope optimisation (Hay-Davies, replicates atlite internals)
# ---------------------------------------------------------------------------

def _annual_poa(
    direct: np.ndarray,
    diffuse: np.ndarray,
    influx_toa: np.ndarray,
    albedo: np.ndarray,
    sun_alt_rad: np.ndarray,
    sun_az_rad: np.ndarray,
    slope_deg: float,
    azimuth_deg: float,
) -> float:
    """Annual sum of plane-of-array irradiance via the Hay-Davies model.

    Replicates atlite's TiltedDirectIrrad + TiltedDiffuseIrrad + TiltedGroundIrrad
    (atlite.pv.irradiation) on pre-extracted single-cell numpy arrays. Solar
    position angles are in radians (as stored in the ERA5 cutout); slope and
    azimuth are in degrees, converted here.
    """
    # --- Previous docstring (kept for reference) below ---
    # Annual sum of plane-of-array irradiance using the Hay-Davies model.
    #
    # Replicates atlite's TiltedDirectIrrad + TiltedDiffuseIrrad + TiltedGroundIrrad
    # (see atlite.pv.irradiation) on pre-extracted single-cell numpy arrays.
    #
    # All solar position angles are in radians (as stored in the ERA5 cutout).
    # Slope and azimuth are in degrees (converted here, matching atlite's
    # make_constant convention).
    slope = np.radians(slope_deg)
    azim  = np.radians(azimuth_deg)

    # Cosine of incidence angle (SurfaceOrientation formula)
    cosinc = np.clip(
        np.sin(slope) * np.cos(sun_alt_rad) * np.cos(azim - sun_az_rad)
        + np.cos(slope) * np.sin(sun_alt_rad),
        0.0, None,
    )

    sinalt = np.clip(np.sin(sun_alt_rad), 1e-3, None)
    R_b = cosinc / sinalt  # geometric factor

    # Direct component
    direct_t = R_b * direct

    # Hay-Davies diffuse (matches atlite's TiltedDiffuseIrrad exactly)
    influx = direct + diffuse
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.sqrt(np.where(influx > 0, direct / influx, 0.0))   # brightening
        A = np.where(influx_toa > 0, direct / influx_toa, 0.0)    # anisotropy
    diffuse_t = np.clip(
        ((1.0 - A) * (1.0 + np.cos(slope)) / 2.0
         * (1.0 + f * np.sin(slope / 2.0) ** 3)
         + A * R_b) * diffuse,
        0.0, None,
    )

    # Ground reflection
    ground_t = albedo * influx * (1.0 - np.cos(slope)) / 2.0

    return float((direct_t + diffuse_t + ground_t).sum())


def _optimal_slope(
    direct, diffuse, influx_toa, albedo, sun_alt_rad, sun_az_rad, azimuth_deg
) -> int:
    """Return the integer slope (0–90°) maximising annual POA at the given azimuth."""
    # --- Previous docstring (kept for reference) below ---
    # Integer slope (0–90°) maximising annual POA at the given azimuth.
    slopes = np.arange(0, 91, _SLOPE_STEP)
    poa = np.array([
        _annual_poa(direct, diffuse, influx_toa, albedo,
                    sun_alt_rad, sun_az_rad, float(s), azimuth_deg)
        for s in slopes
    ])
    return int(slopes[np.argmax(poa)])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Sweep east–west solar orientations at the P95 cell and write per-azimuth CF columns.

    Finds the P95 latitude-optimal cell, then for N azimuths (west→north→east)
    picks the POA-maximising slope and computes the hourly CF, writing one
    `solar_az{az}` column per orientation.
    """
    _OUT.parent.mkdir(parents=True, exist_ok=True)

    cutout = atlite.Cutout(str(_CUTOUT_PATH))
    geom = _get_region_geometry(_REGIONS_PATH, _REGION)

    # Step 1: identify the P95 cell using latitude_optimal as reference
    log.info("computing latitude_optimal CF grid to identify P95 cell")
    cf_ref = cutout.pv(
        panel=_PV_PANEL,
        orientation="latitude_optimal",
        capacity_factor_timeseries=True,
    )
    y_idx, x_idx = _find_p95_cell(cf_ref, geom)
    cell_lat = float(cutout.data.coords["y"].isel(y=y_idx))
    cell_lon = float(cutout.data.coords["x"].isel(x=x_idx))
    annual_mean_ref = float(cf_ref.isel(y=y_idx, x=x_idx).mean())
    log.info(
        f"P95 cell: y={y_idx} x={x_idx}  lat={cell_lat:.2f} lon={cell_lon:.2f}  "
        f"annual_mean_ref={annual_mean_ref:.3f}"
    )

    # Step 2: extract pre-computed irradiance and solar position at P95 cell
    # (atlite pre-computes solar_altitude / solar_azimuth in radians during
    # cutout preparation, so no solar-position call is needed here)
    cell = cutout.data.isel(y=y_idx, x=x_idx).compute()
    direct     = np.clip(cell["influx_direct"].values,  0.0, cell["influx_toa"].values)
    diffuse    = np.clip(cell["influx_diffuse"].values, 0.0, cell["influx_toa"].values)
    influx_toa = cell["influx_toa"].values
    albedo     = cell["albedo"].values
    sun_alt    = cell["solar_altitude"].values   # radians
    sun_az     = cell["solar_azimuth"].values    # radians, clockwise from North

    # Step 3: N equally spaced azimuths, west (270°) → north (0°) → east (90°)
    raw = np.linspace(-90.0, 90.0, _N_STEPS)
    azimuths = np.where(raw < 0, raw + 360.0, raw).round().astype(int)

    # Step 4: optimize slope for each azimuth (pure numpy, milliseconds)
    opt_slopes: dict[int, int] = {}
    for az in azimuths:
        slope = _optimal_slope(direct, diffuse, influx_toa, albedo, sun_alt, sun_az, float(az))
        opt_slopes[az] = slope
        log.info(f"azimuth={az:3d}°  optimal_slope={slope:2d}°")

    # Step 5: compute final CF at optimal (slope, azimuth) via atlite.pv()
    results: dict[str, pd.Series] = {}
    for az, slope in opt_slopes.items():
        log.info(f"computing CF: azimuth={az:3d}° slope={slope:2d}°")
        cf_grid = cutout.pv(
            panel=_PV_PANEL,
            orientation={"slope": float(slope), "azimuth": float(az)},
            capacity_factor_timeseries=True,
        )
        ts = cf_grid.isel(y=y_idx, x=x_idx).to_pandas()
        ts.index = pd.to_datetime(ts.index)
        ts = ts.clip(0.0, 1.0)
        col = f"solar_az{az}"
        results[col] = ts
        log.info(f"  annual_mean_cf={float(ts.mean()):.3f}")

    df = pd.DataFrame(results)
    df.index.name = "time"
    df.to_parquet(_OUT, index=True)
    log.info(f"wrote {_OUT} ({len(df)} rows, {len(results)} orientations)")


if __name__ == "__main__":
    main()
