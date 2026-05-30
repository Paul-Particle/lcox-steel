"""
Plot mean annual CF as a spatial heatmap for a given cutout area and tech,
with the region boundary as a white outline and the P95 bestsite marked.

Snakemake rule: plot_cf_map  (in viz.smk)
Standalone:     python workflow/scripts/viz/plot_cf_map.py [--area aus] [--tech solar]
"""

import argparse
import logging
import sys
from pathlib import Path

import atlite
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

if "snakemake" not in globals():
    # standalone: add workflow/ to path and load stub
    sys.path.insert(0, str(Path(__file__).parents[2]))
    from common._stubs import snakemake

from common._logging import configure_logging
from common._paths import CUTOUTS, SHAPES_RES

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_CF_AREA = "aus"
_TECH = "solar"
_START_DATE = "20250101"
_END_DATE = "20251231"
_CUTOUT_PATH = CUTOUTS / "aus_20250101_20251231.nc"
_REGIONS_PATH = SHAPES_RES / "aus_geo.parquet"
_REGION = "AUS"
_PV_PANEL = "CSi"
_WIND_TURBINE = "Vestas_V112_3MW"
_OUT = Path("results/plots/cf_map/aus_solar_20250101_20251231_cf_map.png")

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _CF_AREA = snakemake.wildcards.cf_area
    _TECH = snakemake.wildcards.tech
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _CUTOUT_PATH = Path(snakemake.input.cutout)
    _REGIONS_PATH = Path(snakemake.input.regions)
    _REGION = snakemake.params.region
    _PV_PANEL = snakemake.params.pv_panel
    _WIND_TURBINE = snakemake.params.wind_onshore_turbine
    _OUT = Path(snakemake.output[0])


def _mask_cells_inside(cell_mean, geom) -> np.ndarray:
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)


def _find_p95_cell(cf_grid, geom):
    cell_mean = cf_grid.mean("time")
    inside = _mask_cells_inside(cell_mean, geom)
    vals = np.where(inside, cell_mean.values, np.nan)
    valid = np.isfinite(vals)
    p95 = np.nanpercentile(vals[valid], 95)
    dist = np.abs(np.where(valid, vals, np.nan) - p95)
    idx_flat = np.nanargmin(dist)
    y_idx, x_idx = np.unravel_index(idx_flat, vals.shape)
    return int(y_idx), int(x_idx)


def main() -> None:
    _OUT.parent.mkdir(parents=True, exist_ok=True)

    log.info("loading cutout %s", _CUTOUT_PATH)
    cutout = atlite.Cutout(str(_CUTOUT_PATH))

    log.info("computing %s CF grid", _TECH)
    if _TECH == "solar":
        cf_grid = cutout.pv(
            panel=_PV_PANEL,
            orientation="latitude_optimal",
            capacity_factor_timeseries=True,
        )
    elif _TECH == "wind_onshore":
        cf_grid = cutout.wind(
            turbine=_WIND_TURBINE,
            capacity_factor_timeseries=True,
            smooth=True,
            add_cutout_windspeed=True,
        )
    else:
        raise ValueError(f"Unsupported tech: {_TECH!r}")

    cell_mean = cf_grid.mean("time")

    gdf = gpd.read_parquet(_REGIONS_PATH).to_crs(4326)
    geom = gdf.geometry.iloc[0]

    y_idx, x_idx = _find_p95_cell(cf_grid, geom)
    p95_lat = float(cutout.data.coords["y"].isel(y=y_idx))
    p95_lon = float(cutout.data.coords["x"].isel(x=x_idx))
    p95_cf  = float(cell_mean.isel(y=y_idx, x=x_idx))
    log.info("P95 cell: lat=%.2f lon=%.2f cf=%.3f", p95_lat, p95_lon, p95_cf)

    lons = cutout.data.coords["x"].values
    lats = cutout.data.coords["y"].values

    fig, ax = plt.subplots(figsize=(10, 8))
    pcm = ax.pcolormesh(lons, lats, cell_mean.values, cmap="YlOrRd",
                        shading="nearest", vmin=0)
    cbar = fig.colorbar(pcm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean annual CF (latitude-optimal)", fontsize=10)

    gdf.boundary.plot(ax=ax, color="white", linewidth=1.5, zorder=3)
    ax.scatter(p95_lon, p95_lat, marker="*", s=250, color="white",
               edgecolors="black", linewidths=0.8, zorder=5,
               label=f"P95 site (CF={p95_cf:.3f})")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.8)

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(
        f"{_CF_AREA.upper()} — {_TECH} mean annual CF ({_START_DATE[:4]})\n"
        f"latitude-optimal orientation · P95 bestsite marked",
        fontsize=11,
    )
    ax.set_aspect("equal")

    fig.savefig(_OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("saved %s", _OUT)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--area",       default=_CF_AREA)
    parser.add_argument("--tech",       default=_TECH)
    parser.add_argument("--start-date", default=_START_DATE)
    parser.add_argument("--end-date",   default=_END_DATE)
    args = parser.parse_args()
    _CF_AREA = args.area; _TECH = args.tech
    _START_DATE = args.start_date; _END_DATE = args.end_date
    _CUTOUT_PATH = CUTOUTS / f"{_CF_AREA}_{_START_DATE}_{_END_DATE}.nc"
    _REGIONS_PATH = SHAPES_RES / f"{_CF_AREA}_geo.parquet"
    _OUT = Path(f"results/plots/cf_map/{_CF_AREA}_{_TECH}_{_START_DATE}_{_END_DATE}_cf_map.png")
    main()
