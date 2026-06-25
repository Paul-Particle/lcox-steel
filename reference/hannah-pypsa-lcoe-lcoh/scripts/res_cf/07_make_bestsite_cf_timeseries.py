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
import yaml
from shapely.geometry import box

YEAR = 2023
OUTDIR = Path("data/res_cf")
CONFIG_PATH = Path("config_hannah.yaml")
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

def haversine_distance_km(lon1, lat1, lon2, lat2):
    """
    Compute great-circle distance in km between one target point and arrays of points.

    lon1, lat1: target point in degrees
    lon2, lat2: arrays of candidate point coordinates in degrees
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

def find_nearest_valid_cell(
    target_x,
    target_y,
    cf_year: xr.DataArray,
    geom,
    max_radius_km: float = 100.0,
    quality_floor_fraction: float = 0.90,
) -> tuple[int, int]:
    """
    Select a deterministic 'best nearby counterpart' cell inside `geom`.

    Rule:
    1. Build all valid cells inside `geom` with finite annual mean CF > 0.
    2. Restrict to cells within `max_radius_km` of the target point.
    3. If nearby candidates exist:
       - find the best nearby annual mean CF
       - keep cells with CF >= quality_floor_fraction * best_nearby_cf
       - choose the nearest among them
    4. If no nearby candidates exist:
       - fall back to the nearest valid cell overall
    """
    cell_mean = cf_year.mean("time")
    inside = mask_cells_inside(cell_mean, geom)

    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)

    mean_vals = cell_mean.values
    valid_mask = inside & np.isfinite(mean_vals) & (mean_vals > 0)

    if not np.any(valid_mask):
        raise ValueError("No valid cells found for counterpart matching.")

    dist_km = haversine_distance_km(target_x, target_y, xx, yy)

    nearby_mask = valid_mask & (dist_km <= max_radius_km)

    if np.any(nearby_mask):
        nearby_cf = np.where(nearby_mask, mean_vals, np.nan)
        best_nearby_cf = np.nanmax(nearby_cf)

        good_nearby_mask = nearby_mask & (
            mean_vals >= quality_floor_fraction * best_nearby_cf
        )

        candidate_dist_km = np.where(good_nearby_mask, dist_km, np.nan)
    else:
        candidate_dist_km = np.where(valid_mask, dist_km, np.nan)

    idx_flat = np.nanargmin(candidate_dist_km)
    y_idx, x_idx = np.unravel_index(idx_flat, dist_km.shape)

    return int(y_idx), int(x_idx)

def get_cell_coords(cf_year: xr.DataArray, y_idx: int, x_idx: int) -> tuple[float, float]:
    x = float(cf_year.x.values[x_idx])
    y = float(cf_year.y.values[y_idx])
    return x, y

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

def load_spatial_matching_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    spatial_cfg = config.get("spatial_matching_res_mix", {})

    return {
        "max_radius_km": float(spatial_cfg.get("max_radius_km", 100.0)),
        "quality_floor_fraction": float(spatial_cfg.get("quality_floor_fraction", 0.90)),
    }

def load_bestsite_res_mix_scenarios() -> list[dict]:
    """
    Load RES-mix scenarios from config_hannah.yaml that require
    scenario-specific bestsite co-located CF files.
    """
    with CONFIG_PATH.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    scenarios = config.get("scenarios", [])

    return [
        s for s in scenarios
        if s.get("tech") == "res_mix"
        and s.get("variant") == "bestsite_p95"
        and isinstance(s.get("res_mix"), dict)
        and s.get("mix_anchor_tech") in TECHS
    ]

def format_res_mix_label(res_mix: dict) -> str:
    """
    Create a deterministic filename-safe label for RES-mix weights.
    Example:
    {"wind_onshore": 0.8, "solar": 0.2}
    -> solar-0p2_wind_onshore-0p8
    """
    parts = []

    for tech in sorted(res_mix.keys()):
        weight = res_mix[tech]
        weight_str = str(float(weight)).replace(".", "p")
        parts.append(f"{tech}-{weight_str}")

    return "_".join(parts)

def build_anchor_colocated_profiles(
    country_upper: str,
    anchor_tech: str,
    res_mix: dict[str, float] | None = None,
) -> dict[str, pd.Series]:
    """
    Build co-located hourly CF profiles for all techs using one anchor-defined site basis.

    Logic:
    - land anchor (wind_onshore, solar):
        - anchor tech uses its own P95 land cell
        - other land tech uses the same land cell
        - offshore uses a matched offshore counterpart cell
    - offshore anchor (wind_offshore):
        - anchor tech uses its own P95 offshore cell
        - one shared land counterpart cell is selected near the offshore anchor
          using a weighted land-tech score based on res_mix
        - both wind_onshore and solar use that same selected land cell for now
          (this is a simplification and could later be extended so each land tech
          gets its own selected counterpart cell)
    """
    if res_mix is None:
        raise ValueError("res_mix must be provided for anchor-based co-located profiles.")

    spatial_cfg = load_spatial_matching_config()

    # Build CF grids once per tech
    cf_year_by_tech = {tech: build_cf_year(country_upper, tech) for tech in TECHS}

    # 1) Find anchor cell in its own valid geometry
    anchor_geom = geometry_for_tech(country_upper, anchor_tech)
    anchor_cf_year = cf_year_by_tech[anchor_tech]
    anchor_y_idx, anchor_x_idx = find_p95_cell(anchor_cf_year, anchor_geom)
    anchor_x, anchor_y = get_cell_coords(anchor_cf_year, anchor_y_idx, anchor_x_idx)

    selected_land_y_idx = None
    selected_land_x_idx = None

    if anchor_tech == "wind_offshore":
        land_geom = geometry_for_tech(country_upper, "wind_onshore")

        wind_onshore_mean = cf_year_by_tech["wind_onshore"].mean("time")
        solar_mean = cf_year_by_tech["solar"].mean("time")

        inside_land = mask_cells_inside(wind_onshore_mean, land_geom)

        xs = wind_onshore_mean.coords["x"].values
        ys = wind_onshore_mean.coords["y"].values
        xx, yy = np.meshgrid(xs, ys)

        dist_km = haversine_distance_km(anchor_x, anchor_y, xx, yy)

        w_onshore = float(res_mix.get("wind_onshore", 0.0))
        w_solar = float(res_mix.get("solar", 0.0))

        if (w_onshore + w_solar) <= 0:
            raise ValueError(
                "Offshore anchor requires at least one land tech in res_mix "
                "(wind_onshore and/or solar)."
            )

        land_score = np.zeros_like(wind_onshore_mean.values, dtype=float)
        valid_land = inside_land.copy()

        if w_onshore > 0:
            onshore_vals = wind_onshore_mean.values
            valid_land = valid_land & np.isfinite(onshore_vals) & (onshore_vals > 0)
            land_score = land_score + w_onshore * onshore_vals

        if w_solar > 0:
            solar_vals = solar_mean.values
            valid_land = valid_land & np.isfinite(solar_vals) & (solar_vals > 0)
            land_score = land_score + w_solar * solar_vals

        nearby_land = valid_land & (dist_km <= spatial_cfg["max_radius_km"])

        if np.any(nearby_land):
            nearby_scores = np.where(nearby_land, land_score, np.nan)
            best_nearby_score = np.nanmax(nearby_scores)

            good_nearby_land = nearby_land & (
                land_score >= spatial_cfg["quality_floor_fraction"] * best_nearby_score
            )

            candidate_dist_km = np.where(good_nearby_land, dist_km, np.nan)
        else:
            candidate_dist_km = np.where(valid_land, dist_km, np.nan)

        idx_flat = np.nanargmin(candidate_dist_km)
        selected_land_y_idx, selected_land_x_idx = np.unravel_index(idx_flat, dist_km.shape)

        selected_land_y_idx = int(selected_land_y_idx)
        selected_land_x_idx = int(selected_land_x_idx)

        selected_land_x, selected_land_y = get_cell_coords(
            cf_year_by_tech["wind_onshore"],
            selected_land_y_idx,
            selected_land_x_idx,
        )

        print(
            f"{country_upper} | offshore anchor -> selected shared land cell "
            f"(x={selected_land_x:.2f}, y={selected_land_y:.2f})"
        )

    results = {}
    selected_cells = {}

    for tech in TECHS:
        tech_cf_year = cf_year_by_tech[tech]

        if tech == anchor_tech:
            y_idx, x_idx = anchor_y_idx, anchor_x_idx

        elif anchor_tech in {"wind_onshore", "solar"} and tech in {"wind_onshore", "solar"}:
            # same land cell basis for land-anchor cases
            y_idx, x_idx = anchor_y_idx, anchor_x_idx

        elif anchor_tech in {"wind_onshore", "solar"} and tech == "wind_offshore":
            # match offshore counterpart to the anchor land cell
            offshore_geom = geometry_for_tech(country_upper, "wind_offshore")
            y_idx, x_idx = find_nearest_valid_cell(
                anchor_x,
                anchor_y,
                tech_cf_year,
                offshore_geom,
                max_radius_km=spatial_cfg["max_radius_km"],
                quality_floor_fraction=spatial_cfg["quality_floor_fraction"],
            )

        elif anchor_tech == "wind_offshore" and tech in {"wind_onshore", "solar"}:
            # For now, both land techs use the same selected land counterpart cell.
            # This could later be extended so solar and wind_onshore each get their
            # own selected land counterpart cell.
            y_idx, x_idx = selected_land_y_idx, selected_land_x_idx

        else:
            raise ValueError(f"Unsupported anchor/counterpart combination: {anchor_tech} -> {tech}")

        # store coordinates of selected cell
        x, y = get_cell_coords(tech_cf_year, y_idx, x_idx)

        selected_cells[tech] = {
            "x": x,
            "y": y,
            "x_idx": int(x_idx),
            "y_idx": int(y_idx),
        }

        # extract time series (unchanged)
        ts = extract_cell_timeseries(tech_cf_year, y_idx, x_idx, tech)
        results[tech] = ts

    return results, selected_cells

def add_location_metadata(
    df: pd.DataFrame,
    anchor_tech: str | None,
    mix_label: str | None,
    selected_cells: dict[str, dict],
) -> pd.DataFrame:
    """
    Add selected-cell metadata to every row of a CF output file.

    Metadata is repeated because CSV is a rectangular table.
    Script 08 can read the first row to recover plotting locations.
    """
    df = df.copy()

    df["anchor_tech"] = anchor_tech if anchor_tech is not None else ""
    df["mix_label"] = mix_label if mix_label is not None else ""

    for tech in TECHS:
        cell = selected_cells[tech]
        df[f"{tech}_x"] = cell["x"]
        df[f"{tech}_y"] = cell["y"]
        df[f"{tech}_x_idx"] = cell["x_idx"]
        df[f"{tech}_y_idx"] = cell["y_idx"]

    return df


def main() -> None:


    for cc in COUNTRIES:
        country_upper = cc.upper()

        results = {}
        selected_cells = {}

        for tech in TECHS:
            cf_year = build_cf_year(country_upper, tech)

            geom = geometry_for_tech(country_upper, tech)

            cell_mean = cf_year.mean("time")
            inside = mask_cells_inside(cell_mean, geom)

            mask = xr.DataArray(inside, coords={"y": cf_year.y, "x": cf_year.x}, dims=("y", "x"))

            cf_national = cf_year.where(mask).mean(dim=("y", "x"))
            nat_mean = float(cf_national.mean().item())

            y_idx, x_idx = find_p95_cell(cf_year, geom)

            x, y = get_cell_coords(cf_year, y_idx, x_idx)

            selected_cells[tech] = {
                "x": x,
                "y": y,
                "x_idx": int(x_idx),
                "y_idx": int(y_idx),
            }

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

        df_p95 = add_location_metadata(
            df_p95,
            anchor_tech=None,
            mix_label=None,
            selected_cells=selected_cells,
        )

        out_p95 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95.csv"
        df_p95.to_csv(out_p95, index=False)

        print(f"{country_upper}: wrote {out_p95.name}")

        # Additional anchor-specific co-located bestsite files
        # For offshore anchor, one shared land counterpart cell is selected for both land techs.
        for anchor_tech in ["wind_onshore", "solar"]:
            res_mix_for_anchor = {anchor_tech: 1.0}

            anchor_results, anchor_selected_cells = build_anchor_colocated_profiles(
                country_upper,
                anchor_tech,
                res_mix=res_mix_for_anchor,
            )

            df_anchor = pd.DataFrame({
                "time": anchor_results["wind_onshore"].index,
                "wind_onshore_cf": anchor_results["wind_onshore"].values,
                "wind_offshore_cf": anchor_results["wind_offshore"].values,
                "solar_cf": anchor_results["solar"].values,
            })

            df_anchor = add_location_metadata(
                df_anchor,
                anchor_tech=anchor_tech,
                mix_label=None,
                selected_cells=anchor_selected_cells,
            )

            out_anchor = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95_anchor-{anchor_tech}.csv"
            df_anchor.to_csv(out_anchor, index=False)

            print(f"{country_upper}: wrote {out_anchor.name}")

        # Scenario-specific co-located bestsite files for res_mix + bestsite_p95
        # These use the actual res_mix weights from config_hannah.yaml.
        res_mix_scenarios = load_bestsite_res_mix_scenarios()

        for scenario in res_mix_scenarios:
            if scenario["country"].upper() != country_upper:
                continue

            anchor_tech = scenario["mix_anchor_tech"]
            res_mix = scenario["res_mix"]
            mix_label = format_res_mix_label(res_mix)

            if anchor_tech in {"wind_onshore", "solar"}:
                print(
                    f"{country_upper}: skipping scenario-specific file for land anchor "
                    f"{anchor_tech} because land-anchor mix files are identical to generic anchor files."
                )
                continue

            scenario_results, scenario_selected_cells = build_anchor_colocated_profiles(
                country_upper,
                anchor_tech,
                res_mix=res_mix,
            )

            df_scenario_anchor = pd.DataFrame({
                "time": scenario_results["wind_onshore"].index,
                "wind_onshore_cf": scenario_results["wind_onshore"].values,
                "wind_offshore_cf": scenario_results["wind_offshore"].values,
                "solar_cf": scenario_results["solar"].values,
            })

            df_scenario_anchor = add_location_metadata(
                df_scenario_anchor,
                anchor_tech=anchor_tech,
                mix_label=mix_label,
                selected_cells=scenario_selected_cells,
            )

            out_scenario_anchor = OUTDIR / (
                f"{cc}_cf_{YEAR}_bestsite_p95_"
                f"anchor-{anchor_tech}_mix-{mix_label}.csv"
            )

            df_scenario_anchor.to_csv(out_scenario_anchor, index=False)

            print(f"{country_upper}: wrote {out_scenario_anchor.name}")




if __name__ == "__main__":
    main()
