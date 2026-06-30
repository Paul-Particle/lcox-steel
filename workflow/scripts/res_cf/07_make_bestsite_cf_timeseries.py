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
- Reads an annual Atlite cutout
- Computes hourly CF per grid cell for each technology using:
    - wind: smoothed turbine power curve (smooth=True)
    - solar: standard PV conversion
- Computes annual mean CF per grid cell
- Identifies the representative P95 grid cell (mean-based ranking)
- Extracts that cell's hourly CF time series
- For anchor co-location, matches counterpart cells within `max_radius_km`
  subject to a `quality_floor_fraction` of the best nearby CF

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

# NOTE: `from __future__ import annotations` (present in the reference script)
# cannot be used here — Snakemake prepends its preamble before running the
# script, so this would no longer be the file's first statement (SyntaxError).
# The type hints below are valid on the pinned Python (3.12) without it.

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import atlite
import geopandas as gpd
import yaml
from shapely.geometry import box
import logging
import re
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

YEAR = 2023
OUTDIR = RES_CF
CONFIG_PATH = load_res_cf_cfg()
COUNTRIES = ["de"]  # lowercase to match filenames

CUTOUT_DIR = None
REGIONS_PATH = SHAPES_RES / "regions.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"

SM_VARIANT:            str | None = None   # "bestsite-p95" or "anchored-w{W}-s{S}"
SM_ANCHOR:             str | None = None   # anchor tech in snake_case (from tech wildcard)
SM_RES_MIX:            dict[str, float] | None = None  # parsed from variant; offshore anchor only
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRIES             = [snakemake.wildcards.cf_area.lower()]
    CONFIG_PATH            = snakemake.config["res_cf"]
    REGIONS_PATH          = Path(snakemake.input.regions)
    OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    CUTOUT_DIR           = Path(snakemake.input.cutout)
    SM_VARIANT            = snakemake.wildcards.variant
    SM_ANCHOR             = snakemake.wildcards.tech.replace("-", "_")
    if SM_VARIANT.startswith("anchored-"):
        _m = re.match(r"anchored-w(\d+)-s(\d+)$", SM_VARIANT)
        assert _m, f"Cannot parse variant '{SM_VARIANT}'"
        _w = int(_m.group(1)) / 10
        _s = int(_m.group(2)) / 10
        SM_RES_MIX = {"wind_onshore": _s, "solar": max(0.0, round(1.0 - _w - _s, 10))}

TECHS = ["wind_onshore", "wind_offshore", "solar"]

WIND_TURBINE          = CONFIG_PATH["wind_onshore_turbine"]
WIND_OFFSHORE_TURBINE = CONFIG_PATH["wind_offshore_turbine"]
PV_PANEL              = CONFIG_PATH["pv_panel"]
PV_ORIENTATION        = CONFIG_PATH["pv_orientation"]

WIND_ONSHORE_TURBINE  = WIND_TURBINE
WIND_CF_CFG           = CONFIG_PATH.get("wind_cf", {})
WIND_SMOOTH           = WIND_CF_CFG.get("smooth", True)
WIND_ADD_CUTOUT_WS    = WIND_CF_CFG.get("add_cutout_windspeed", True)

# Spatial matching for anchor co-location (max search radius + how close to the
# best nearby cell a candidate must be to qualify).
_SPATIAL_CFG          = CONFIG_PATH.get("spatial_matching_res_mix", {})
MAX_RADIUS_KM         = float(_SPATIAL_CFG.get("max_radius_km", 100.0))
QUALITY_FLOOR_FRAC    = float(_SPATIAL_CFG.get("quality_floor_fraction", 0.90))

# RES-mix scenarios needing scenario-specific co-located files. Each entry:
#   {country, mix_anchor_tech, res_mix: {tech: weight, ...}}  (weights sum to 1).
# Empty by default; populate the `res_cf.bestsite_res_mix_scenarios` config list.
RES_MIX_SCENARIOS     = CONFIG_PATH.get("bestsite_res_mix_scenarios", [])


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
    co = atlite.Cutout(path=str(cutout_path))

    if tech == "wind_onshore":
        cf = co.wind(
            turbine=WIND_TURBINE,
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


def load_land_geometry(iso2: str):
    gdf = gpd.read_parquet(REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"{iso2} not found in {REGIONS_PATH}")
    return row.geometry.iloc[0]


def load_offshore_geometry(iso2: str):
    gdf = gpd.read_parquet(OFFSHORE_REGIONS_PATH)
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
    """Return the (y, x) index of the in-region cell closest to the P95 annual-mean CF."""
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


def get_cell_coords(cf_year: xr.DataArray, y_idx: int, x_idx: int) -> tuple[float, float]:
    x = float(cf_year.x.values[x_idx])
    y = float(cf_year.y.values[y_idx])
    return x, y


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
    # Build CF grids once per tech
    cf_year_by_tech = {tech: build_cf_year(cutout_path, tech) for tech in TECHS}

    # 1) Find anchor cell in its own valid geometry.
    anchor_cf_year = cf_year_by_tech[anchor_tech]
    anchor_geom = geometry_for_tech(country_upper, anchor_tech)
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

        nearby_land = valid_land & (dist_km <= MAX_RADIUS_KM)

        if np.any(nearby_land):
            nearby_scores = np.where(nearby_land, land_score, np.nan)
            best_nearby_score = np.nanmax(nearby_scores)

            good_nearby_land = nearby_land & (
                land_score >= QUALITY_FLOOR_FRAC * best_nearby_score
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

        log.info(
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
                max_radius_km=MAX_RADIUS_KM,
                quality_floor_fraction=QUALITY_FLOOR_FRAC,
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

        # extract time series
        ts = extract_cell_timeseries(tech_cf_year, y_idx, x_idx)
        results[tech] = ts

    return results, selected_cells


def format_res_mix_label(res_mix: dict[str, float]) -> str:
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


def add_location_metadata(
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


def _write_sm_p95_output(results: dict[str, pd.Series]) -> None:
    if "snakemake" not in globals() or SM_VARIANT != "bestsite-p95":
        return
    results[SM_ANCHOR].rename(SM_ANCHOR.replace("_", "-")).to_frame().to_parquet(
        snakemake.output[0], index=True
    )


def _write_sm_anchored_output(profiles: dict[str, pd.Series], anchor_tech: str) -> None:
    if "snakemake" not in globals() or anchor_tech != SM_ANCHOR:
        return
    if not (SM_VARIANT and SM_VARIANT.startswith("anchored-")):
        return
    pd.DataFrame(
        {tech.replace("_", "-"): profiles[tech] for tech in TECHS},
    ).to_parquet(snakemake.output[0], index=True)


def main() -> None:
    for cc in COUNTRIES:
        country_upper = cc.upper()
        cutout = CUTOUT_DIR or annual_cutout_path(cc, YEAR)

        # 1) Per-tech best-site P95 (each tech from its own P95 cell).
        results: dict[str, pd.Series] = {}
        selected_cells: dict[str, dict] = {}
        for tech in TECHS:
            cf_year = build_cf_year(cutout, tech)
            geom = geometry_for_tech(country_upper, tech)

            cell_mean = cf_year.mean("time")
            inside = mask_cells_inside(cell_mean, geom)
            mask = xr.DataArray(inside, coords={"y": cf_year.y, "x": cf_year.x}, dims=("y", "x"))
            cf_national = cf_year.where(mask).mean(dim=("y", "x"))
            nat_mean = float(cf_national.mean().item())

            y_idx, x_idx = find_p95_cell(cf_year, geom)
            x, y = get_cell_coords(cf_year, y_idx, x_idx)
            selected_cells[tech] = {"x": x, "y": y, "x_idx": int(x_idx), "y_idx": int(y_idx)}
            ts = extract_cell_timeseries(cf_year, y_idx, x_idx)
            best_mean = float(ts.mean())
            results[tech] = ts
            log.info(f"{country_upper} | {tech}: national_mean={nat_mean:.3f} best_mean={best_mean:.3f}")

        df_p95 = add_location_metadata(_to_dataframe(results), None, None, selected_cells)
        out_p95 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95.parquet"
        df_p95.to_parquet(out_p95, index=False)
        log.info(f"{country_upper}: wrote {out_p95.name}")
        _write_sm_p95_output(results)

        # 2) Anchor-co-located files (one per anchor tech: onshore, solar, offshore).
        for anchor_tech in ["wind_onshore", "solar", "wind_offshore"]:
            if SM_ANCHOR == anchor_tech and SM_RES_MIX is not None:
                anchor_res_mix = SM_RES_MIX
            elif anchor_tech == "wind_offshore":
                anchor_res_mix = {"wind_onshore": 0.5, "solar": 0.5}
            else:
                anchor_res_mix = {anchor_tech: 1.0}
            profiles, cells = build_anchor_colocated_profiles(
                cutout, country_upper, anchor_tech, res_mix=anchor_res_mix
            )
            df_anchor = add_location_metadata(_to_dataframe(profiles), anchor_tech, None, cells)
            out_anchor = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95_anchor-{anchor_tech}.parquet"
            df_anchor.to_parquet(out_anchor, index=False)
            log.info(f"{country_upper}: wrote {out_anchor.name}")
            _write_sm_anchored_output(profiles, anchor_tech)

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
            mix_label = format_res_mix_label(res_mix)
            profiles, cells = build_anchor_colocated_profiles(cutout, country_upper, anchor_tech, res_mix)
            df_mix = add_location_metadata(_to_dataframe(profiles), anchor_tech, mix_label, cells)
            out_mix = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95_anchor-{anchor_tech}_mix-{mix_label}.parquet"
            df_mix.to_parquet(out_mix, index=False)
            log.info(f"{country_upper}: wrote {out_mix.name}")


if __name__ == "__main__":
    main()
