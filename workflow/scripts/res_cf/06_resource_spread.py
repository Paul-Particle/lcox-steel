"""
06_resource_spread.py

Purpose
-------
Quantifies intra-country spatial resource heterogeneity ("best-site uplift")
for wind onshore, wind offshore, and solar using existing Atlite cutouts (2023).

Method
------
- Computes hourly CF per grid cell (no national aggregation)
- Derives annual mean CF per grid cell
- Calculates area-weighted spatial statistics within country polygon:
    * spatial mean (sanity check)
    * spatial P90
    * spatial P95
    * spatial max
- Computes uplift vs national mean CF

Important
---------
- Uses same technology assumptions as national pipeline
- Uses indicator-matrix-based area weights for methodological consistency
- Purely climate-resource based (no land-use, grid, slope or permitting constraints)
- Results represent theoretical upper envelope within country

Output
------
resources/res_cf/resource_spread_2023.parquet
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
from scripts.res_cf._helpers import QUARTERS, cutout_path, load_res_cf_cfg

configure_logging(snakemake)
log = logging.getLogger(__name__)


REGIONS_PATH          = SHAPES_RES / "regions.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"
NATIONAL_CF_DIR       = RES_CF

YEAR        = 2023
OUT_PATH    = RES_CF / f"resource_spread_{YEAR}.parquet"
RES_CF_CFG  = load_res_cf_cfg()

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    OUT_PATH   = Path(snakemake.output[0])
    RES_CF_CFG = snakemake.config["res_cf"]

# Explicit list of countries to process in standalone mode. The old `enabled`
# flag was removed from config.yaml, so we hardcode the target countries here
# until this script is wired into the Snakemake pipeline.
COUNTRIES = {
    info["region"]: QUARTERS
    for key, info in RES_CF_CFG["countries"].items()
    if key in ("de", "fr", "es", "aus", "bra")
}

TECHS = ["wind_onshore", "wind_offshore", "solar"]

WIND_ONSHORE_TURBINE  = RES_CF_CFG["wind_onshore_turbine"]
WIND_OFFSHORE_TURBINE = RES_CF_CFG["wind_offshore_turbine"]
PV_PANEL              = RES_CF_CFG["pv_panel"]
PV_ORIENTATION        = RES_CF_CFG["pv_orientation"]
WIND_CF_CFG           = RES_CF_CFG.get("wind_cf", {})
WIND_SMOOTH           = WIND_CF_CFG.get("smooth", True)
WIND_ADD_CUTOUT_WS    = WIND_CF_CFG.get("add_cutout_windspeed", True)

def weighted_percentile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Area-weighted percentile of `values` (nonnegative `weights`), q in [0, 1]."""
    # --- Previous docstring (kept for reference) below ---
    # Area-weighted percentile of 'values' with nonnegative 'weights'.
    # q in [0, 1].
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

def load_country_geometry(iso2: str):
    gdf = gpd.read_parquet(REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"Country {iso2} not found in {REGIONS_PATH}")
    return row.geometry.iloc[0]

def load_offshore_geometry(iso2: str):
    gdf = gpd.read_parquet(OFFSHORE_REGIONS_PATH)
    row = gdf.loc[gdf["region"] == iso2]
    if row.empty:
        raise ValueError(f"Country {iso2} not found in {OFFSHORE_REGIONS_PATH}")
    return row.geometry.iloc[0]

def compute_cf_grid(cutout: atlite.Cutout, tech: str) -> xr.DataArray:
    """Return the per-cell CF grid (time, y, x) for `tech`, matching the national pipeline's atlite assumptions."""
    # --- Previous docstring (kept for reference) below ---
    # Returns CF grid with dims (time, y, x).
    # Keeps Atlite tech assumptions identical to national pipeline.
    if tech == "wind_onshore":
        cf = cutout.wind(turbine=WIND_ONSHORE_TURBINE, capacity_factor=True,
                         smooth=WIND_SMOOTH, add_cutout_windspeed=WIND_ADD_CUTOUT_WS)

    elif tech == "wind_offshore":
        cf = cutout.wind(turbine=WIND_OFFSHORE_TURBINE, capacity_factor=True,
                         smooth=WIND_SMOOTH, add_cutout_windspeed=WIND_ADD_CUTOUT_WS)

    elif tech == "solar":
        cf = cutout.pv(panel=PV_PANEL, orientation=PV_ORIENTATION, capacity_factor=True)

    else:
        raise ValueError(f"Unknown tech: {tech}")

    # Ensure DataArray
    if isinstance(cf, xr.Dataset):
        # If Atlite returns a Dataset, take the first variable or the known one
        # Adjust if your Atlite version names it differently
        var = list(cf.data_vars)[0]
        cf = cf[var]

    return cf

def build_weights(cutout: atlite.Cutout, country_geom) -> xr.DataArray:
    """Build per-cell area weights (cell_area × in-country fraction) as a (y, x) DataArray."""
    # --- Previous docstring (kept for reference) below ---
    # Build per-cell weights = cell_area * in-country fraction.
    # Returns DataArray with dims (y, x).
    # indicator matrix: fraction of each grid cell inside the country polygon
    # In some Atlite versions this returns a scipy sparse matrix (shape: (n_shapes, n_cells))
    indicator = cutout.indicatormatrix([country_geom]).tocsr()
    indicator_1d = np.asarray(indicator[0, :].todense()).ravel()  # (n_cells,)

    # Cell area (prefer true area if available)
    if hasattr(cutout, "grid") and isinstance(cutout.grid, xr.Dataset) and "area" in cutout.grid:
        cell_area = cutout.grid["area"]
        # ensure 2D (y, x)
        for dim in list(cell_area.dims):
            if dim not in ("y", "x"):
                cell_area = cell_area.isel({dim: 0})
    else:
        # fallback proxy if area not provided by atlite version
        lat = cutout.data.coords["y"]
        cell_area = xr.DataArray(
            np.cos(np.deg2rad(lat.values))[:, None] * np.ones((lat.size, cutout.data.coords["x"].size)),
            coords={"y": cutout.data.coords["y"], "x": cutout.data.coords["x"]},
            dims=("y", "x"),
        )

    cell_area_1d = cell_area.values.ravel()
    weights_1d = indicator_1d * cell_area_1d

    weights = xr.DataArray(
        weights_1d.reshape(cell_area.shape),
        coords=cell_area.coords,
        dims=cell_area.dims,
    )
    return weights



def national_mean_from_csv(iso2: str, tech: str) -> float:
    """Return the national-mean CF for (iso2, tech) from the per-country parquet, or NaN if absent."""
    p = NATIONAL_CF_DIR / f"{iso2.lower()}_cf_{YEAR}.parquet"
    if not p.exists():
        return np.nan
    df = pd.read_parquet(p)
    if tech == "wind_onshore":
        col = "wind_onshore_cf"
    elif tech == "wind_offshore":
        col = "wind_offshore_cf"
    else:
        col = "solar_cf"
    return float(df[col].mean())

def main():
    """Compute area-weighted spatial CF stats (P90/P95/max + uplift) per country/tech and write parquet."""
    rows = []
    regions_geom_cache = {}

    for iso2, segments in COUNTRIES.items():
        if iso2 not in regions_geom_cache:
            regions_geom_cache[iso2] = load_country_geometry(iso2)
        geom = regions_geom_cache[iso2]


        for tech in TECHS:
            # Build full-year CF grid by concatenating segments along time
            geom = load_offshore_geometry(iso2) if tech == "wind_offshore" else regions_geom_cache[iso2]
            first_cutout = atlite.Cutout(path=str(cutout_path(iso2, YEAR, segments[0])))
            weights = build_weights(first_cutout, geom)
            w = weights.values.ravel()
            cf_parts = []
            for seg in segments:
                p = cutout_path(iso2, YEAR, seg)
                if not p.exists():
                    raise FileNotFoundError(f"Missing cutout: {p}")
                co = atlite.Cutout(path=str(p))
                cf_grid = compute_cf_grid(co, tech)  # (time, y, x)
                cf_parts.append(cf_grid)

            cf_year = xr.concat(cf_parts, dim="time")
            # Annual mean per cell
            cell_mean = cf_year.mean("time")  # (y, x)
            v = cell_mean.values.ravel()

            # Spatial stats (weighted by area*country-fraction)
            spatial_mean = float(np.nansum(v * w) / np.nansum(w)) if np.nansum(w) > 0 else np.nan
            p90 = weighted_percentile(v, w, 0.90)
            p95 = weighted_percentile(v, w, 0.95)

            # Max over in-country cells only
            m = np.isfinite(v) & (w > 0)
            vmax = float(np.nanmax(v[m])) if np.any(m) else np.nan

            # National mean from your authoritative national CSV (best)
            nat_mean = national_mean_from_csv(iso2, tech)

            # If CSV missing, fall back to spatial_mean (still consistent)
            if not np.isfinite(nat_mean):
                nat_mean = spatial_mean

            uplift_p95 = (p95 / nat_mean - 1.0) if np.isfinite(p95) and nat_mean > 0 else np.nan
            uplift_p90 = (p90 / nat_mean - 1.0) if np.isfinite(p90) and nat_mean > 0 else np.nan

            rows.append({
                "country": iso2,
                "tech": tech,
                "national_mean": nat_mean,
                "spatial_p90_mean": p90,
                "spatial_p95_mean": p95,
                "spatial_max_mean": vmax,
                "uplift_p95_vs_national": uplift_p95,
                "uplift_p90_vs_national": uplift_p90,
            })

            log.info(
                f"{iso2} {tech}: national_mean={nat_mean:.4f} spatial_mean={spatial_mean:.4f} p95={p95:.4f}"
            )

    out = pd.DataFrame(rows)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT_PATH, index=False)
    log.info(f"wrote {OUT_PATH}")

if __name__ == "__main__":
    main()
