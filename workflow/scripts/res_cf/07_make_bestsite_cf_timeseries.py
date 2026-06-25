"""
07_make_bestsite_cf_timeseries.py

Purpose
-------
Creates "best-site" hourly CF time series by extracting physically consistent,
site-level CF profiles from Atlite CF grids for wind onshore, wind offshore, and solar.

This replaces the previous uplift-factor-based approach with direct extraction
from the underlying climate data.

Two families of output are produced:

1. Per-tech best-site P95 — each technology is taken from its own P95 grid cell
   (not co-located across technologies).
2. Anchor-co-located profiles — one "anchor" technology fixes the site, and the
   other technologies are evaluated at that same location basis (land techs share
   the anchor's land cell; offshore is matched to the nearest high-quality
   offshore counterpart, and vice-versa). These feed RES-mix scenarios.

Method
------
- Reads Atlite cutouts for each country segment (Q1-Q4)
- Computes hourly CF per grid cell for each technology using:
    - wind: smoothed turbine power curve (smooth=True)
    - solar: standard PV conversion
- Computes annual mean CF per grid cell
- Identifies the representative P95 grid cell (mean-based ranking)
- Extracts that cell's hourly CF time series
- For anchor co-location, matches counterpart cells within `max_radius_km`
  subject to a `quality_floor_fraction` of the best nearby CF

Single-cell vs spatial averaging
--------------------------------
The earlier iteration of this script applied a 3x3 spatial average to wind cells
to mimic farm-scale variability. The active pipeline dropped that ("single-cell
sampling is sufficient"), so this version samples single cells. If farm-scale
smoothing is wanted back, reintroduce it in `extract_cell_timeseries`.

Notes
-----
- Purely climate-resource based (no land-use, grid, or permitting constraints)
- Best-site CFs represent high-quality project locations within national boundaries
- National CFs remain unchanged and are generated separately (Script 03-05)

Outputs (parquet, columns: time, wind_onshore_cf, wind_offshore_cf, solar_cf,
plus selected-cell metadata columns)
-------
resources/res_cf/<cc>_cf_2023_bestsite_p95.parquet
resources/res_cf/<cc>_cf_2023_bestsite_p95_anchor-<anchor>.parquet
resources/res_cf/<cc>_cf_2023_bestsite_p95_anchor-<anchor>_mix-<label>.parquet
"""

from __future__ import annotations

import logging
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import atlite
import geopandas as gpd

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from common._paths import RES_CF, SHAPES_RES
from scripts.res_cf._helpers import (
    annual_cutout_path,
    haversine_distance_km,
    load_res_cf_cfg,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)


RES_CF_CFG = load_res_cf_cfg()
YEAR       = 2023
OUTDIR     = RES_CF
COUNTRIES  = ["de"]  # lowercase to match filenames; standalone default

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRIES  = [snakemake.wildcards.country.lower()]
    RES_CF_CFG = snakemake.config["res_cf"]

REGIONS_PATH          = SHAPES_RES / "regions.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"

TECHS                 = ["wind_onshore", "wind_offshore", "solar"]
WIND_ONSHORE_TURBINE  = RES_CF_CFG["wind_onshore_turbine"]
WIND_OFFSHORE_TURBINE = RES_CF_CFG["wind_offshore_turbine"]
PV_PANEL              = RES_CF_CFG["pv_panel"]
PV_ORIENTATION        = RES_CF_CFG["pv_orientation"]
WIND_CF_CFG           = RES_CF_CFG.get("wind_cf", {})
WIND_SMOOTH           = WIND_CF_CFG.get("smooth", True)
WIND_ADD_CUTOUT_WS    = WIND_CF_CFG.get("add_cutout_windspeed", True)

# Spatial matching for anchor co-location (max search radius + how close to the
# best nearby cell a candidate must be to qualify).
_SPATIAL_CFG          = RES_CF_CFG.get("spatial_matching_res_mix", {})
MAX_RADIUS_KM         = float(_SPATIAL_CFG.get("max_radius_km", 100.0))
QUALITY_FLOOR_FRAC    = float(_SPATIAL_CFG.get("quality_floor_fraction", 0.90))

# RES-mix scenarios needing scenario-specific co-located files. Each entry:
#   {country, mix_anchor_tech, res_mix: {tech: weight, ...}}  (weights sum to 1).
# Empty by default; populate the `res_cf.bestsite_res_mix_scenarios` config list.
RES_MIX_SCENARIOS     = RES_CF_CFG.get("bestsite_res_mix_scenarios", [])


def extract_cell_timeseries(
    cf_year: xr.DataArray,
    y_idx: int,
    x_idx: int,
) -> pd.Series:
    """Extract one grid cell's hourly CF series (clipped to [0, 1], 'time'-indexed)."""
    s = cf_year.isel(y=y_idx, x=x_idx).to_pandas()
    s.index = pd.to_datetime(s.index)
    s.name = "cf"
    return s.clip(0, 1)


def build_cf_year(cutout_path: Path, tech: str) -> xr.DataArray:
    """Build the full-year per-cell CF grid for `tech` from a single annual cutout."""
    co = atlite.Cutout(path=str(cutout_path))

    if tech == "wind_onshore":
        cf = co.wind(
            turbine=WIND_ONSHORE_TURBINE,
            capacity_factor_timeseries=True,
            smooth=WIND_SMOOTH,
            add_cutout_windspeed=WIND_ADD_CUTOUT_WS,
        )
    elif tech == "wind_offshore":
        cf = co.wind(
            turbine=WIND_OFFSHORE_TURBINE,
            capacity_factor_timeseries=True,
            smooth=WIND_SMOOTH,
            add_cutout_windspeed=WIND_ADD_CUTOUT_WS,
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

    return cf


def _load_land_geometry(iso2: str):
    gdf = gpd.read_parquet(REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"{iso2} not found in {REGIONS_PATH}")
    return row.geometry.iloc[0]


def _load_offshore_geometry(iso2: str):
    gdf = gpd.read_parquet(OFFSHORE_REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"{iso2} not found in {OFFSHORE_REGIONS_PATH}")
    return row.geometry.iloc[0]


def geometry_for_tech(iso2: str, tech: str):  # returns shapely geometry
    return _load_offshore_geometry(iso2) if tech == "wind_offshore" else _load_land_geometry(iso2)


def mask_cells_inside(cell_mean: xr.DataArray, geom) -> np.ndarray:
    """Return a boolean grid marking cells whose centre lies within `geom`."""
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)


def _find_p95_cell(cf_year: xr.DataArray, geom) -> tuple[int, int]:
    """Return the (y, x) index of the in-region cell closest to the P95 annual-mean CF."""
    cell_mean = cf_year.mean("time")
    inside = mask_cells_inside(cell_mean, geom)
    vals = np.where(inside, cell_mean.values, np.nan)
    valid = np.isfinite(vals)
    p95 = np.nanpercentile(vals[valid], 95)
    # representative P95 cell = valid cell whose value is closest to the P95 threshold
    dist = np.abs(np.where(valid, vals, np.nan) - p95)
    idx_flat = np.nanargmin(dist)
    y_idx, x_idx = np.unravel_index(idx_flat, vals.shape)
    return int(y_idx), int(x_idx)


def _cell_coords(cf_year: xr.DataArray, y_idx: int, x_idx: int) -> tuple[float, float]:
    """Return the (lon, lat) centre of grid cell (y_idx, x_idx)."""
    return float(cf_year.x.values[x_idx]), float(cf_year.y.values[y_idx])


def find_nearest_valid_cell(
    target_x: float,
    target_y: float,
    cf_year: xr.DataArray,
    geom,
    max_radius_km: float = MAX_RADIUS_KM,
    quality_floor_fraction: float = QUALITY_FLOOR_FRAC,
) -> tuple[int, int]:
    """Select a deterministic 'best nearby counterpart' cell inside `geom`.

    1. Build all valid cells inside `geom` with finite annual-mean CF > 0.
    2. Restrict to cells within `max_radius_km` of the target point.
    3. If nearby candidates exist: keep those with CF >= quality_floor_fraction
       * best-nearby-CF, then pick the nearest of them.
    4. Otherwise: fall back to the nearest valid cell overall.
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
        best_nearby_cf = np.nanmax(np.where(nearby_mask, mean_vals, np.nan))
        good_nearby = nearby_mask & (mean_vals >= quality_floor_fraction * best_nearby_cf)
        candidate_dist_km = np.where(good_nearby, dist_km, np.nan)
    else:
        candidate_dist_km = np.where(valid_mask, dist_km, np.nan)

    idx_flat = np.nanargmin(candidate_dist_km)
    y_idx, x_idx = np.unravel_index(idx_flat, dist_km.shape)
    return int(y_idx), int(x_idx)


def build_anchor_colocated_profiles(
    cutout_path: Path,
    country_upper: str,
    anchor_tech: str,
    res_mix: dict[str, float],
) -> tuple[dict[str, pd.Series], dict[str, dict]]:
    """Build co-located hourly CF profiles for all techs from one anchor-defined site.

    - Land anchor (wind_onshore / solar): the anchor tech uses its own P95 land
      cell; the other land tech reuses that cell; offshore is matched to the
      nearest high-quality offshore counterpart.
    - Offshore anchor (wind_offshore): the anchor uses its own P95 offshore cell;
      a single shared land counterpart cell is chosen near the offshore anchor by
      a res_mix-weighted land-CF score, and both land techs use it (a deliberate
      simplification — could later give each land tech its own counterpart).

    Returns (profiles_by_tech, selected_cells_by_tech).
    """
    cf_year_by_tech = {tech: build_cf_year(cutout_path, tech) for tech in TECHS}

    # 1) Anchor cell in its own geometry.
    anchor_cf_year = cf_year_by_tech[anchor_tech]
    anchor_geom = geometry_for_tech(country_upper, anchor_tech)
    anchor_y_idx, anchor_x_idx = _find_p95_cell(anchor_cf_year, anchor_geom)
    anchor_x, anchor_y = _cell_coords(anchor_cf_year, anchor_y_idx, anchor_x_idx)

    # 2) For an offshore anchor, pick one shared land counterpart cell weighted by res_mix.
    selected_land_y_idx = selected_land_x_idx = None
    if anchor_tech == "wind_offshore":
        land_geom = geometry_for_tech(country_upper, "wind_onshore")
        wind_onshore_mean = cf_year_by_tech["wind_onshore"].mean("time")
        solar_mean = cf_year_by_tech["solar"].mean("time")

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

        valid_land = mask_cells_inside(wind_onshore_mean, land_geom)
        land_score = np.zeros_like(wind_onshore_mean.values, dtype=float)
        if w_onshore > 0:
            onshore_vals = wind_onshore_mean.values
            valid_land = valid_land & np.isfinite(onshore_vals) & (onshore_vals > 0)
            land_score = land_score + w_onshore * onshore_vals
        if w_solar > 0:
            solar_vals = solar_mean.values
            valid_land = valid_land & np.isfinite(solar_vals) & (solar_vals > 0)
            land_score = land_score + w_solar * solar_vals

        nearby_land = valid_land & (dist_km <= MAX_RADIUS_KM)
        if np.any(nearby_land):
            best_nearby_score = np.nanmax(np.where(nearby_land, land_score, np.nan))
            good_nearby = nearby_land & (land_score >= QUALITY_FLOOR_FRAC * best_nearby_score)
            candidate_dist_km = np.where(good_nearby, dist_km, np.nan)
        else:
            candidate_dist_km = np.where(valid_land, dist_km, np.nan)

        idx_flat = np.nanargmin(candidate_dist_km)
        selected_land_y_idx, selected_land_x_idx = (
            int(i) for i in np.unravel_index(idx_flat, dist_km.shape)
        )
        sl_x, sl_y = _cell_coords(
            cf_year_by_tech["wind_onshore"], selected_land_y_idx, selected_land_x_idx
        )
        log.info(
            f"{country_upper} | offshore anchor -> shared land cell (x={sl_x:.2f}, y={sl_y:.2f})"
        )

    results: dict[str, pd.Series] = {}
    selected_cells: dict[str, dict] = {}
    for tech in TECHS:
        tech_cf_year = cf_year_by_tech[tech]

        if tech == anchor_tech:
            y_idx, x_idx = anchor_y_idx, anchor_x_idx
        elif anchor_tech in {"wind_onshore", "solar"} and tech in {"wind_onshore", "solar"}:
            y_idx, x_idx = anchor_y_idx, anchor_x_idx  # shared land cell
        elif anchor_tech in {"wind_onshore", "solar"} and tech == "wind_offshore":
            y_idx, x_idx = find_nearest_valid_cell(
                anchor_x, anchor_y, tech_cf_year, geometry_for_tech(country_upper, "wind_offshore")
            )
        elif anchor_tech == "wind_offshore" and tech in {"wind_onshore", "solar"}:
            y_idx, x_idx = selected_land_y_idx, selected_land_x_idx
        else:
            raise ValueError(f"Unsupported anchor/counterpart combination: {anchor_tech} -> {tech}")

        x, y = _cell_coords(tech_cf_year, y_idx, x_idx)
        selected_cells[tech] = {"x": x, "y": y, "x_idx": int(x_idx), "y_idx": int(y_idx)}
        results[tech] = extract_cell_timeseries(tech_cf_year, y_idx, x_idx)

    return results, selected_cells


def _format_res_mix_label(res_mix: dict[str, float]) -> str:
    """Deterministic filename-safe label, e.g. {wind_onshore:0.8, solar:0.2} -> solar-0p2_wind_onshore-0p8."""
    parts = [f"{tech}-{str(float(res_mix[tech])).replace('.', 'p')}" for tech in sorted(res_mix)]
    return "_".join(parts)


def _to_dataframe(profiles: dict[str, pd.Series]) -> pd.DataFrame:
    """Assemble the standard (time + per-tech CF) dataframe from co-located profiles."""
    return pd.DataFrame({
        "time": profiles["wind_onshore"].index,
        "wind_onshore_cf": profiles["wind_onshore"].values,
        "wind_offshore_cf": profiles["wind_offshore"].values,
        "solar_cf": profiles["solar"].values,
    })


def _add_location_metadata(
    df: pd.DataFrame,
    anchor_tech: str | None,
    mix_label: str | None,
    selected_cells: dict[str, dict],
) -> pd.DataFrame:
    """Repeat selected-cell metadata on every row so 08/100 can recover plotting locations."""
    df = df.copy()
    df["anchor_tech"] = anchor_tech or ""
    df["mix_label"] = mix_label or ""
    for tech in TECHS:
        cell = selected_cells[tech]
        df[f"{tech}_x"] = cell["x"]
        df[f"{tech}_y"] = cell["y"]
        df[f"{tech}_x_idx"] = cell["x_idx"]
        df[f"{tech}_y_idx"] = cell["y_idx"]
    return df


def main() -> None:
    """Write per-tech P95, generic anchor, and scenario-specific RES-mix bestsite parquets."""
    for cc in COUNTRIES:
        country_upper = cc.upper()
        cutout = annual_cutout_path(cc, YEAR)

        # 1) Per-tech best-site P95 (each tech from its own P95 cell).
        results: dict[str, pd.Series] = {}
        selected_cells: dict[str, dict] = {}
        for tech in TECHS:
            cf_year = build_cf_year(cutout, tech)
            geom = geometry_for_tech(country_upper, tech)

            cell_mean = cf_year.mean("time")
            inside = mask_cells_inside(cell_mean, geom)
            mask = xr.DataArray(inside, coords={"y": cf_year.y, "x": cf_year.x}, dims=("y", "x"))
            nat_mean = float(cf_year.where(mask).mean(dim=("y", "x")).mean().item())

            y_idx, x_idx = _find_p95_cell(cf_year, geom)
            x, y = _cell_coords(cf_year, y_idx, x_idx)
            selected_cells[tech] = {"x": x, "y": y, "x_idx": int(x_idx), "y_idx": int(y_idx)}
            ts = extract_cell_timeseries(cf_year, y_idx, x_idx)
            results[tech] = ts
            log.info(f"{country_upper} | {tech}: national_mean={nat_mean:.3f} best_mean={float(ts.mean()):.3f}")

        df_p95 = _add_location_metadata(_to_dataframe(results), None, None, selected_cells)
        out_p95 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95.parquet"
        df_p95.to_parquet(out_p95, index=False)
        log.info(f"{country_upper}: wrote {out_p95.name}")

        # 2) Generic anchor-co-located files (one land anchor per file).
        for anchor_tech in ["wind_onshore", "solar"]:
            profiles, cells = build_anchor_colocated_profiles(
                cutout, country_upper, anchor_tech, res_mix={anchor_tech: 1.0}
            )
            df_anchor = _add_location_metadata(_to_dataframe(profiles), anchor_tech, None, cells)
            out_anchor = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95_anchor-{anchor_tech}.parquet"
            df_anchor.to_parquet(out_anchor, index=False)
            log.info(f"{country_upper}: wrote {out_anchor.name}")

        # 3) Scenario-specific RES-mix files (config-driven; offshore anchors etc.).
        for scenario in RES_MIX_SCENARIOS:
            if scenario["country"].upper() != country_upper:
                continue
            anchor_tech = scenario["mix_anchor_tech"]
            res_mix = scenario["res_mix"]
            if anchor_tech in {"wind_onshore", "solar"}:
                log.info(
                    f"{country_upper}: skipping land-anchor mix file for {anchor_tech} "
                    "(identical to the generic anchor file)."
                )
                continue
            mix_label = _format_res_mix_label(res_mix)
            profiles, cells = build_anchor_colocated_profiles(cutout, country_upper, anchor_tech, res_mix)
            df_mix = _add_location_metadata(_to_dataframe(profiles), anchor_tech, mix_label, cells)
            out_mix = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95_anchor-{anchor_tech}_mix-{mix_label}.parquet"
            df_mix.to_parquet(out_mix, index=False)
            log.info(f"{country_upper}: wrote {out_mix.name}")


if __name__ == "__main__":
    main()
