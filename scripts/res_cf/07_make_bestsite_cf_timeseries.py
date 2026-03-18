"""
07_make_bestsite_cf_timeseries.py

Purpose
-------
Creates "best-site" hourly CF time series by extracting physically consistent,
site-level CF profiles from Atlite CF grids for wind onshore, wind offshore, and solar.

This replaces the previous uplift-factor-based approach with direct extraction
from the underlying climate data.

Method
------
- Reads Atlite cutouts for each country segment (Q1–Q4)
- Computes hourly CF per grid cell for each technology using:
    - wind: smoothed turbine power curve (smooth=True)
    - solar: standard PV conversion
- Computes annual mean CF per grid cell
- Identifies the representative P95 grid cell (mean-based ranking)
- Extracts that cell's hourly CF time series
- Applies 3×3 spatial averaging for wind technologies to approximate
  wind farm-scale variability and avoid artificial saturation at CF ≈ 1
- Writes the resulting hourly series as the best-site P95 scenario

Key updates vs previous implementation
-------------------------------------
- ❌ Removed uplift-factor scaling (P90/P95 × national CF)
- ❌ Removed dependence on resource_spread outputs for time series generation
- ✅ Direct extraction from full CF grid (cf_year)
- ✅ Introduced wind power-curve smoothing (smooth=True)
- ✅ Added 3×3 spatial averaging for wind to approximate site-level behaviour
- ✅ Uses mean-based P95 selection (median-based selection tested and rejected)

Notes
-----
- Purely climate-resource based (no land-use, grid, or permitting constraints)
- Best-site CFs represent high-quality project locations within national boundaries
- National CFs remain unchanged and are generated separately (Script 03–05)
- Internal consistency check (not enforced in code):
    best-site mean CF ≥ national mean CF

Outputs
-------
data/res_cf/<cc>_cf_2023_bestsite_p95.csv

(columns: time, wind_onshore_cf, wind_offshore_cf, solar_cf)
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import atlite
import geopandas as gpd
from shapely.geometry import box

YEAR = 2023
OUTDIR = Path("data/res_cf")

COUNTRIES = ["de", "fr", "es", "aus", "bra"]  # lowercase to match filenames

CUTOUT_DIR = Path("data/cutouts")
REGIONS_PATH = Path("data/shapes/regions.geojson")
OFFSHORE_REGIONS_PATH = Path("data/shapes/offshore_regions.geojson")

COUNTRY_SEGMENTS = {
    "DE": ["q1", "q2", "q3", "q4"],
    "FR": ["q1", "q2", "q3", "q4"],
    "ES": ["q1", "q2", "q3", "q4"],
    "AUS": ["q1", "q2", "q3", "q4"],
    "BRA": ["q1", "q2", "q3", "q4"],
}

TECHS = ["wind_onshore", "wind_offshore", "solar"]

WIND_TURBINE = "Vestas_V112_3MW"
WIND_OFFSHORE_TURBINE = "NREL_ReferenceTurbine_5MW_offshore"
PV_PANEL = "CSi"
PV_ORIENTATION = "latitude_optimal"


def extract_cell_timeseries(cf_year, y_idx, x_idx, tech):

    if tech.startswith("wind"):
        ts = cf_year.isel(
            y=slice(max(0, y_idx-1), y_idx+2),
            x=slice(max(0, x_idx-1), x_idx+2)
        ).mean(dim=("y", "x"))
    else:
        ts = cf_year.isel(y=y_idx, x=x_idx)

    s = ts.to_pandas()
    s.index = pd.to_datetime(s.index)
    s.name = "cf"

    return s.clip(0, 1)

def build_cf_year(country_upper: str, tech: str) -> xr.DataArray:
    parts = []

    for seg in COUNTRY_SEGMENTS[country_upper]:
        cutout_path = CUTOUT_DIR / f"{country_upper.lower()}_{YEAR}_{seg}.nc"
        co = atlite.Cutout(path=str(cutout_path))

        if tech == "wind_onshore":
            cf = co.wind(
                turbine=WIND_TURBINE,
                capacity_factor_timeseries=True,
                smooth=True,
                add_cutout_windspeed=True,
            )
        elif tech == "wind_offshore":
            cf = co.wind(
                turbine=WIND_OFFSHORE_TURBINE,
                capacity_factor_timeseries=True,
                smooth=True,
                add_cutout_windspeed=True,
            )
        elif tech == "solar":
            cf = co.pv(
                panel=PV_PANEL,
                orientation=PV_ORIENTATION,
                capacity_factor_timeseries=True,
            )
        else:
            raise ValueError(f"Unknown tech: {tech}")

        if isinstance(cf, xr.Dataset):
            cf = cf[list(cf.data_vars)[0]]

        parts.append(cf)

    return xr.concat(parts, dim="time")

def load_land_geometry(iso2: str):
    gdf = gpd.read_file(REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"{iso2} not found in {REGIONS_PATH}")
    return row.geometry.iloc[0]


def load_offshore_geometry(iso2: str):
    gdf = gpd.read_file(OFFSHORE_REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"{iso2} not found in {OFFSHORE_REGIONS_PATH}")
    return row.geometry.iloc[0]

def geometry_for_tech(iso2: str, tech: str):
    return load_offshore_geometry(iso2) if tech == "wind_offshore" else load_land_geometry(iso2)


def mask_cells_inside(cell_mean: xr.DataArray, geom) -> np.ndarray:
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)

def find_p95_cell(cf_year: xr.DataArray, geom) -> tuple[int, int]:
    """
    Identify the representative P95 grid cell based on annual mean CF.

    Returns
    -------
    (y_idx, x_idx)
        Index of the selected P95 grid cell.
    """
    cell_mean = cf_year.mean("time")

    inside = mask_cells_inside(cell_mean, geom)

    vals = np.where(inside, cell_mean.values, np.nan)

    valid = np.isfinite(vals)
    v = vals[valid]

    p95 = np.nanpercentile(v, 95)

    # representative P95 cell = valid cell whose value is closest to the P95 threshold
    dist = np.abs(np.where(valid, vals, np.nan) - p95)
    idx_flat = np.nanargmin(dist)
    y_idx, x_idx = np.unravel_index(idx_flat, vals.shape)

    return int(y_idx), int(x_idx)

def to_cf_series(x, name="cf"):
    obj = x.to_pandas()

    if isinstance(obj, pd.DataFrame):
        if obj.shape[1] == 1:
            s = obj.iloc[:, 0]
        else:
            s = obj.iloc[:, 0]
    else:
        s = obj

    s.index = pd.to_datetime(s.index)
    s.index.name = "time"
    return s.rename(name).clip(0, 1)



def main() -> None:


    for cc in COUNTRIES:
        country_upper = cc.upper()

        results = {}

        for tech in TECHS:
            cf_year = build_cf_year(country_upper, tech)

            geom = geometry_for_tech(country_upper, tech)

            cell_mean = cf_year.mean("time")
            inside = mask_cells_inside(cell_mean, geom)

            mask = xr.DataArray(inside, coords={"y": cf_year.y, "x": cf_year.x}, dims=("y", "x"))

            cf_national = cf_year.where(mask).mean(dim=("y", "x"))
            nat_mean = float(cf_national.mean().item())

            y_idx, x_idx = find_p95_cell(cf_year, geom)

            ts = extract_cell_timeseries(cf_year, y_idx, x_idx, tech)

            best_mean = float(ts.mean())

            print(f"{country_upper} | {tech}")
            print(f"  national_mean = {nat_mean:.3f}")
            print(f"  best_mean     = {best_mean:.3f}")

            results[tech] = ts


        df_p95 = pd.DataFrame({
            "time": results["wind_onshore"].index,
            "wind_onshore_cf": results["wind_onshore"].values,
            "wind_offshore_cf": results["wind_offshore"].values,
            "solar_cf": results["solar"].values,
        })

        out_p95 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95.csv"
        df_p95.to_csv(out_p95, index=False)

        print(f"{country_upper}: wrote {out_p95.name}")




if __name__ == "__main__":
    main()
