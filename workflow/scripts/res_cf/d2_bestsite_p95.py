"""
d2_bestsite_p95.py

Purpose
-------
Creates "best-site" hourly CF time series by extracting physically consistent,
site-level CF profiles from Atlite CF grids for wind onshore, wind offshore, and
solar — one independent best cell per technology (no co-location across techs).

This replaces the previous uplift-factor-based approach with direct extraction
from the underlying climate data.

Anchor-based co-location (one technology fixing the site for the others) has
moved to d3_anchor_colo.py, which also replaces the
old proximity+quality-floor counterpart matching with complementarity-based
multi-candidate selection.

Method
------
- Reads an annual Atlite cutout
- Computes hourly CF per grid cell for each technology using:
    - wind: smoothed turbine power curve (smooth=True)
    - solar: standard PV conversion
- Computes annual mean CF per grid cell
- Builds area-based cell weights via atlite's indicatormatrix (area-fraction
  overlap with the region), matching the method used in scripts 03 and 06
- Identifies the representative P95 grid cell using area-weighted percentile
  ranking (not a plain unweighted percentile — matters for large/distorted
  countries where grid-cell area varies with latitude)
- Extracts that cell's hourly CF time series
- Also computes and saves a per-country/tech resource summary: area-weighted
  national mean, P90, P95, and max CF (no uplift ratios — that was an obsolete
  metric from the old methodology, dropped when 06_resource_spread.py was
  folded into this script and then deleted)

Notes
-----
- Purely climate-resource based (no land-use, grid, or permitting constraints)
- Best-site CFs represent high-quality project locations within national boundaries
- National (area-average) CFs are generated separately by script 03
- Internal consistency check (not enforced in code):
    best-site mean CF >= national mean CF

Outputs
-------
resources/timeseries/<area>_cf_2023_bestsite_p95.parquet
    Per-tech independent P95 CF time series (time, wind_onshore_cf,
    wind_offshore_cf, solar_cf, plus selected-cell metadata columns).
resources/timeseries/<area>_cf_2023_resource_summary.parquet
    One row per country/tech: national_mean, p90, p95, max (all area-weighted),
    plus run_timestamp_utc.
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
if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from common._paths import SHAPES_RES, TIMESERIES
from scripts.res_cf._helpers_res_cf import (
    annual_cutout_path,
    eligibility_weights,
    load_res_cf_cfg,
    pick_p95_cell,
    weighted_percentile,
)
configure_logging(snakemake)
log = logging.getLogger(__name__)

YEAR = 2023
OUTDIR = TIMESERIES
RES_CF_CFG = load_res_cf_cfg()
AREAS = ["DE_LU"]          # area codes, as in config/scenarios.csv
REGION = "DE"              # region tag stored in the geo parquet

CUTOUT_DIR = None # test: Path("cutouts/vic_20250101_20251231.nc")
REGIONS_PATH = SHAPES_RES / "regions.parquet" # test: "vic_geo.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"  # test: "vic_offshore_geo.parquet"

SM_VARIANT:            str | None = None   # "bestsite-p95" or "anchored-w{W}-s{S}"
SM_ANCHOR:             str | None = None   # anchor tech in snake_case (from tech wildcard)
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    AREAS                 = [snakemake.wildcards.area]
    REGION                = snakemake.params.region
    RES_CF_CFG            = snakemake.config["res_cf"]
    REGIONS_PATH          = Path(snakemake.input.regions)
    OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    CUTOUT_DIR           = Path(snakemake.input.cutout)
    SM_VARIANT            = snakemake.wildcards.variant
    SM_ANCHOR             = snakemake.wildcards.tech.replace("-", "_")

TECHS = ["wind_onshore", "wind_offshore", "solar"]

WIND_TURBINE          = RES_CF_CFG["wind_onshore_turbine"]
WIND_OFFSHORE_TURBINE = RES_CF_CFG["wind_offshore_turbine"]
PV_PANEL              = RES_CF_CFG["pv_panel"]
PV_ORIENTATION        = RES_CF_CFG["pv_orientation"]

WIND_ONSHORE_TURBINE  = WIND_TURBINE
WIND_CF_CFG           = RES_CF_CFG.get("wind_cf", {})
WIND_SMOOTH           = WIND_CF_CFG.get("smooth", True)
WIND_ADD_CUTOUT_WS    = WIND_CF_CFG.get("add_cutout_windspeed", True)
# Land-sea eligibility cutoff (#41). Onshore only — offshore passes 0 below.
MIN_LAND_FRACTION     = float(RES_CF_CFG.get("min_land_fraction", 0.0))
ELIGIBILITY_SOURCE    = RES_CF_CFG.get("eligibility_source", "indicatormatrix")


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

def get_cell_coords(cf_year: xr.DataArray, y_idx: int, x_idx: int) -> tuple[float, float]:
    x = float(cf_year.x.values[x_idx])
    y = float(cf_year.y.values[y_idx])
    return x, y

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

def main() -> None:
    for area in AREAS:
        cutout = CUTOUT_DIR or annual_cutout_path(area, YEAR)

        # 1) Per-tech best-site P95 (each tech from its own P95 cell).
        results: dict[str, pd.Series] = {}
        selected_cells: dict[str, dict] = {}
        summary_rows: list[dict] = []
        for tech in TECHS:
            cf_year = build_cf_year(cutout, tech)
            geom = geometry_for_tech(REGION, tech)

            co = atlite.Cutout(path=str(cutout))
            # Land-sea eligibility (#41): onshore drops sea-contaminated coastal
            # border cells; offshore keeps all sea cells (min_lf = 0). Weights are
            # also cos(lat) physical-area weighted (#37), folded into eligibility_weights.
            min_lf = MIN_LAND_FRACTION if tech != "wind_offshore" else 0.0
            weights = eligibility_weights(co, geom, min_lf, ELIGIBILITY_SOURCE)

            cell_mean = cf_year.mean("time").values
            w_flat = weights.ravel()
            v_flat = cell_mean.ravel()
            nat_mean = float(np.nansum(v_flat * w_flat) / np.nansum(w_flat)) if np.nansum(w_flat) > 0 else np.nan

            p95 = weighted_percentile(v_flat, w_flat, 0.95)
            p90 = weighted_percentile(v_flat, w_flat, 0.90)

            valid_mask = w_flat > 0
            cell_max = float(np.nanmax(np.where(valid_mask, v_flat, np.nan))) if np.any(valid_mask) else np.nan

            y_idx, x_idx = pick_p95_cell(cell_mean, weights)
            x, y = get_cell_coords(cf_year, y_idx, x_idx)
            selected_cells[tech] = {"x": x, "y": y, "x_idx": int(x_idx), "y_idx": int(y_idx)}
            ts = extract_cell_timeseries(cf_year, y_idx, x_idx)
            best_mean = float(ts.mean())
            results[tech] = ts
            log.info(f"{area} | {tech}: national_mean={nat_mean:.3f} best_mean={best_mean:.3f}")
            summary_rows.append({
                "area": area,
                "tech": tech,
                "national_mean": nat_mean,
                "p90": p90,
                "p95": p95,
                "max": cell_max,
            })

        df_p95 = add_location_metadata(_to_dataframe(results), None, None, selected_cells)
        out_p95 = OUTDIR / f"{area}_cf_{YEAR}_bestsite_p95.parquet"
        df_p95.to_parquet(out_p95, index=False)
        log.info(f"{area}: wrote {out_p95.name}")
        _write_sm_p95_output(results)

        df_summary = pd.DataFrame(summary_rows)
        df_summary["run_timestamp_utc"] = pd.Timestamp.now("UTC").isoformat()
        out_summary = OUTDIR / f"{area}_cf_{YEAR}_resource_summary.parquet"
        df_summary.to_parquet(out_summary, index=False)
        log.info(f"{area}: wrote {out_summary.name}")


if __name__ == "__main__":
    main()
