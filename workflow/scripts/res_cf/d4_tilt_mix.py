"""Build N solar CF time series spanning the east–west orientation trade-off.

Panels can be oriented anywhere from due west (azimuth 270°) through the
equator-facing direction (maximum annual yield) to due east (90°, maximum
morning shift). The equator-facing direction depends on hemisphere: north (0°)
in the southern hemisphere, south (180°) in the northern hemisphere. This
script:

  1. Identifies the bestsite cell — P95 of annual-mean latitude_optimal CF —
     so all orientations compare at the same high-quality location. The
     site-selection helpers are copied from `d2_bestsite_p95.py`.
  2. Sweeps N equally spaced azimuths spanning ±90° around the cell's
     equator-facing direction — i.e. the full east–west range, passing
     through equator-facing at the midpoint. Use an odd N so that midpoint
     (the maximum annual-yield orientation) is actually sampled.
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

from common._paths import CUTOUTS, SHAPES_RES, TIMESERIES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.res_cf._helpers import geom_area_weights, pick_p95_cell

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_AREA = "aus"
_N_STEPS = 7
_START_DATE = "20250101"
_END_DATE = "20251231"
_CUTOUT_PATH = CUTOUTS / "aus_20250101_20251231.nc"
_REGIONS_PATH = SHAPES_RES / "aus_geo.parquet"
_REGION = "AUS"
_PV_PANEL = "CSi"
_OUT = TIMESERIES / "AUS_solar_tilt-mix-n7_20250101_20251231.parquet"

# Slope search resolution for the optimization loop (degrees).
# 1° is pure numpy so it runs in milliseconds regardless of resolution.
_SLOPE_STEP = 1

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _AREA = snakemake.wildcards.area
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

    Finds the P95 latitude-optimal cell, then for N azimuths spanning ±90°
    around the cell's equator-facing direction (the full east–west range)
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
    weights = geom_area_weights(cutout, geom)
    y_idx, x_idx = pick_p95_cell(cf_ref.mean("time"), weights)
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

    # Step 3: N equally spaced azimuths spanning ±90° around the cell's
    # equator-facing direction — north (0°) in the southern hemisphere, south
    # (180°) in the northern hemisphere. Covers the full east–west range, with
    # equator-facing at the midpoint (sampled only when N is odd).
    equator_facing_az = 0.0 if cell_lat < 0 else 180.0
    log.info(
        f"cell_lat={cell_lat:.2f} → equator-facing azimuth {equator_facing_az:.0f}° "
        f"({'north, SH' if cell_lat < 0 else 'south, NH'})"
    )
    if _N_STEPS % 2 == 0:
        log.warning(
            f"n_steps={_N_STEPS} is even → the equator-facing azimuth "
            f"({equator_facing_az:.0f}°, maximum annual yield) is not sampled; "
            f"use an odd n_steps to include it"
        )
    raw = np.linspace(equator_facing_az - 90.0, equator_facing_az + 90.0, _N_STEPS)
    azimuths = np.mod(raw, 360.0).round().astype(int)

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
