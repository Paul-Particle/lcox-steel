"""
07_make_bestsite_cf_timeseries.py

Purpose
-------
Creates "best-site" hourly CF time series by extracting physically consistent,
site-level CF profiles from Atlite CF grids for wind onshore, wind offshore, and
solar — one independent best cell per technology (no co-location across techs).

This replaces the previous uplift-factor-based approach with direct extraction
from the underlying climate data.

Anchor-based co-location (one technology fixing the site for the others) has
moved to 07b_make_anchor_colocated_cf_timeseries.py, which also replaces the
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
- National (country-average) CFs are generated separately by script 03
- Internal consistency check (not enforced in code):
    best-site mean CF >= national mean CF

Outputs
-------
resources/res_cf/<cc>_cf_2023_bestsite_p95.parquet
    Per-tech independent P95 CF time series (time, wind_onshore_cf,
    wind_offshore_cf, solar_cf, plus selected-cell metadata columns).
resources/res_cf/<cc>_cf_2023_resource_summary.parquet
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
from common._paths import RES_CF, SHAPES_RES
from scripts.res_cf._helpers import (
    annual_cutout_path,
    load_res_cf_cfg,
)
configure_logging(snakemake)
log = logging.getLogger(__name__)

YEAR = 2023
OUTDIR = RES_CF
RES_CF_CFG = load_res_cf_cfg()
COUNTRIES = ["de"]  # lowercase to match filenames

CUTOUT_DIR = None # test: Path("cutouts/vic_20250101_20251231.nc")
REGIONS_PATH = SHAPES_RES / "regions.parquet" # test: "vic_geo.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"  # test: "vic_offshore_geo.parquet"

SM_VARIANT:            str | None = None   # "bestsite-p95" or "anchored-w{W}-s{S}"
SM_ANCHOR:             str | None = None   # anchor tech in snake_case (from tech wildcard)
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRIES             = [snakemake.wildcards.cf_area.lower()]
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

def weighted_percentile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Area-weighted percentile of `values` (nonnegative `weights`), q in [0, 1]."""
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

def find_p95_cell(cf_year: xr.DataArray, weights: np.ndarray) -> tuple[int, int]:
    """Return the (y, x) index of the cell closest to the area-weighted P95
    annual-mean CF.

    `weights` is a (y, x) array of area-based weights (e.g. from
    cutout.indicatormatrix). Cells outside the region already have weight 0,
    so they're naturally excluded — no separate geometry check needed here.
    """
    cell_mean = cf_year.mean("time").values

    p95 = weighted_percentile(cell_mean.ravel(), weights.ravel(), 0.95)

    valid = weights.ravel() > 0
    vals = np.where(valid, cell_mean.ravel(), np.nan)
    dist = np.abs(vals - p95)
    idx_flat = np.nanargmin(dist)
    y_idx, x_idx = np.unravel_index(idx_flat, cell_mean.shape)

    return int(y_idx), int(x_idx)

def geom_weights(cutout: atlite.Cutout, geom) -> np.ndarray:
    """Area-based weights for `geom` as a flat (y, x)-shaped array, built from
    atlite's indicatormatrix (area-fraction overlap) — same method used in
    scripts 03 and 06.
    """
    indicator = cutout.indicatormatrix([geom]).tocsr()
    indicator_1d = np.asarray(indicator[0, :].todense()).ravel()
    n_y = cutout.data.sizes["y"]
    n_x = cutout.data.sizes["x"]
    return indicator_1d.reshape(n_y, n_x)

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
    for cc in COUNTRIES:
        country_upper = cc.upper()
        cutout = CUTOUT_DIR or annual_cutout_path(cc, YEAR)

        # 1) Per-tech best-site P95 (each tech from its own P95 cell).
        results: dict[str, pd.Series] = {}
        selected_cells: dict[str, dict] = {}
        summary_rows: list[dict] = []
        for tech in TECHS:
            cf_year = build_cf_year(cutout, tech)
            geom = geometry_for_tech(country_upper, tech)

            co = atlite.Cutout(path=str(cutout))
            weights = geom_weights(co, geom)

            cell_mean = cf_year.mean("time").values
            w_flat = weights.ravel()
            v_flat = cell_mean.ravel()
            nat_mean = float(np.nansum(v_flat * w_flat) / np.nansum(w_flat)) if np.nansum(w_flat) > 0 else np.nan

            p95 = weighted_percentile(v_flat, w_flat, 0.95)
            p90 = weighted_percentile(v_flat, w_flat, 0.90)

            valid_mask = w_flat > 0
            cell_max = float(np.nanmax(np.where(valid_mask, v_flat, np.nan))) if np.any(valid_mask) else np.nan

            y_idx, x_idx = find_p95_cell(cf_year, weights)
            x, y = get_cell_coords(cf_year, y_idx, x_idx)
            selected_cells[tech] = {"x": x, "y": y, "x_idx": int(x_idx), "y_idx": int(y_idx)}
            ts = extract_cell_timeseries(cf_year, y_idx, x_idx)
            best_mean = float(ts.mean())
            results[tech] = ts
            log.info(f"{country_upper} | {tech}: national_mean={nat_mean:.3f} best_mean={best_mean:.3f}")
            summary_rows.append({
                "country": country_upper,
                "tech": tech,
                "national_mean": nat_mean,
                "p90": p90,
                "p95": p95,
                "max": cell_max,
            })

        df_p95 = add_location_metadata(_to_dataframe(results), None, None, selected_cells)
        out_p95 = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95.parquet"
        df_p95.to_parquet(out_p95, index=False)
        log.info(f"{country_upper}: wrote {out_p95.name}")
        _write_sm_p95_output(results)

        df_summary = pd.DataFrame(summary_rows)
        df_summary["run_timestamp_utc"] = pd.Timestamp.now("UTC").isoformat()
        out_summary = OUTDIR / f"{cc}_cf_{YEAR}_resource_summary.parquet"
        df_summary.to_parquet(out_summary, index=False)
        log.info(f"{country_upper}: wrote {out_summary.name}")


if __name__ == "__main__":
    main()
