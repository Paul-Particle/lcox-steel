"""
08_plot_bestsite_locations.py

Purpose
-------
Visualise the spatial distribution of annual mean CF per grid cell and
highlight:

• the single best grid cell
• the P95 resource region

for each country and technology.

Outputs
-------
results/plots/bestsite_locations/<country>_<tech>.png
"""

from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import atlite
import regionmask
import pandas as pd

CUTOUT_DIR = Path("cutouts")
REGIONS_PATH = Path("resources/shapes/regions.geojson")
OFFSHORE_REGIONS_PATH = Path("resources/shapes/offshore_regions.geojson")

OUTDIR = Path("results/plots/bestsite_locations")
OUTDIR.mkdir(parents=True, exist_ok=True)

YEAR = 2023

COUNTRIES = {
    "DE": ["q1","q2","q3","q4"],
    "FR": ["q1","q2","q3","q4"],
    "ES": ["q1","q2","q3","q4"],
    "AUS": ["q1","q2","q3","q4"],
    "BRA": ["q1","q2","q3","q4"],
}

TECHS = ["wind_onshore", "wind_offshore", "solar"]

WIND_TURBINE = "Vestas_V112_3MW"
WIND_OFFSHORE_TURBINE = "NREL_ReferenceTurbine_5MW_offshore"
PV_PANEL = "CSi"
PV_ORIENTATION = "latitude_optimal"

def get_national_mean_from_csv(iso2: str, tech: str) -> float:
    path = Path(f"resources/res_cf/{iso2.lower()}_cf_2023.csv")
    df = pd.read_csv(path)

    col_map = {
        "wind_onshore": "wind_onshore_cf",
        "wind_offshore": "wind_offshore_cf",
        "solar": "solar_cf",
    }

    return float(df[col_map[tech.strip()]].mean())

def land_mask(cell_mean):
    """
    Returns boolean mask of land cells with the same dims/coords as cell_mean.
    """
    land = regionmask.defined_regions.natural_earth_v5_0_0.land_110

    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values

    mask = land.mask(xs, ys)

    return xr.DataArray(
        ~np.isnan(mask.values),
        coords={"y": cell_mean.coords["y"], "x": cell_mean.coords["x"]},
        dims=("y", "x"),
    )

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


def cutout_path(iso2: str, seg: str) -> Path:
    return CUTOUT_DIR / f"{iso2.lower()}_{YEAR}_{seg}.nc"


def compute_cf_grid(cutout: atlite.Cutout, tech: str) -> xr.DataArray:
    if tech == "wind_onshore":
        cf = cutout.wind(
            turbine=WIND_TURBINE,
            capacity_factor_timeseries=True,
            smooth=True,
            add_cutout_windspeed=True,
        )
    elif tech == "wind_offshore":
        cf = cutout.wind(
            turbine=WIND_OFFSHORE_TURBINE,
            capacity_factor_timeseries=True,
            smooth=True,
            add_cutout_windspeed=True,
        )
    elif tech == "solar":
        cf = cutout.pv(
            panel=PV_PANEL,
            orientation=PV_ORIENTATION,
            capacity_factor_timeseries=True,
        )
    else:
        raise ValueError(f"Unknown tech: {tech}")

    if isinstance(cf, xr.Dataset):
        var = list(cf.data_vars)[0]
        cf = cf[var]

    return cf

def build_annual_mean_grid(iso2: str, tech: str) -> xr.DataArray:
    parts = []
    for seg in COUNTRIES[iso2]:
        p = cutout_path(iso2, seg)
        if not p.exists():
            raise FileNotFoundError(f"Missing cutout: {p}")
        co = atlite.Cutout(path=str(p))
        cf = compute_cf_grid(co, tech)   # (time, y, x)
        parts.append(cf)

    cf_year = xr.concat(parts, dim="time")
    return cf_year.mean("time")          # (y, x)


def geometry_for_tech(iso2: str, tech: str):
    return load_offshore_geometry(iso2) if tech == "wind_offshore" else load_land_geometry(iso2)


def mask_cells_inside(cell_mean: xr.DataArray, geom) -> np.ndarray:
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)


def best_and_p95_masks(cell_mean: xr.DataArray, geom):
    inside = mask_cells_inside(cell_mean, geom)
    vals = cell_mean.values.copy()

    valid = np.isfinite(vals) & inside
    if not np.any(valid):
        raise ValueError("No valid in-geometry cells found.")

    v = vals[valid]
    p95_threshold = np.nanpercentile(v, 95)

    # single representative P95 cell = valid cell whose value is closest to the P95 threshold
    dist = np.abs(np.where(valid, vals, np.nan) - p95_threshold)
    p95_idx_flat = np.nanargmin(dist)
    p95_mask = np.zeros(vals.shape, dtype=bool)
    p95_mask[np.unravel_index(p95_idx_flat, vals.shape)] = True

    # single best cell
    best_idx_flat = np.nanargmax(np.where(valid, vals, np.nan))
    best_mask = np.zeros(vals.shape, dtype=bool)
    best_mask[np.unravel_index(best_idx_flat, vals.shape)] = True

    return p95_mask, best_mask, p95_threshold

def plot_bestsite_map(iso2: str, tech: str) -> None:
    geom = geometry_for_tech(iso2, tech)
    cell_mean = build_annual_mean_grid(iso2, tech)

    # use a masked copy only for onshore best/P95 selection
    cell_mean_select = cell_mean
    if tech == "wind_onshore":
        land = land_mask(cell_mean)
        cell_mean_select = cell_mean.where(land)

    p95_mask, best_mask, p95_threshold = best_and_p95_masks(cell_mean_select, geom)
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values

    fig, ax = plt.subplots(figsize=(8, 6))

    # base CF heatmap
    im = ax.pcolormesh(xs, ys, cell_mean.values, shading="auto")

    # always draw land-country outline for orientation
    land_boundary = gpd.GeoSeries([load_land_geometry(iso2)], crs=4326)
    land_boundary.boundary.plot(ax=ax, color="white", linewidth=2.5)

    # for offshore maps, also draw offshore polygon boundary
    if tech == "wind_offshore":
        offshore_boundary = gpd.GeoSeries([geom], crs=4326)
        offshore_boundary.boundary.plot(ax=ax, color="white", linewidth=1.2, linestyle="--")

    # P95 cells
    yy, xx = np.where(p95_mask)
    if len(xx) > 0:
        ax.scatter(
            xs[xx],
            ys[yy],
            marker="s",
            s=80,
            facecolors="none",
            edgecolors="white",
            linewidths=2,
            label="P95 cell",
            zorder=9,
        )

        # 3x3 footprint (wind only)
    if tech.startswith("wind") and len(xx) > 0:
        y0, x0 = yy[0], xx[0]

        y_min = max(0, y0 - 1)
        y_max = min(len(ys) - 1, y0 + 1)
        x_min = max(0, x0 - 1)
        x_max = min(len(xs) - 1, x0 + 1)

        xs_foot = xs[x_min:x_max+1]
        ys_foot = ys[y_min:y_max+1]

        xx_f, yy_f = np.meshgrid(xs_foot, ys_foot)

        ax.scatter(
            xx_f.flatten(),
            yy_f.flatten(),
            marker="s",
            s=40,
            facecolors="none",
            edgecolors="yellow",
            linewidths=1,
            label="3x3 footprint",
            zorder=8,
        )

    # best cell
    yyb, xxb = np.where(best_mask)
    ax.scatter(
        xs[xxb],
        ys[yyb],
        marker="x",
        s=120,
        linewidths=3,
        color="white",
        label="Best cell",
        zorder=10,
    )

    # summary stats box
    grid_mean = float(np.nanmean(cell_mean.values[mask_cells_inside(cell_mean, geom)]))
    national_mean = get_national_mean_from_csv(iso2, tech)
    best_cf = float(cell_mean.values[best_mask].item())

    stats_text = (
        f"National mean (area-weighted): {national_mean:.3f}\n"
        f"Grid mean (unweighted): {grid_mean:.3f}\n"
        f"P95 CF: {p95_threshold:.3f}\n"
        f"Best CF: {best_cf:.3f}"
    )

    ax.text(
        0.02, 0.98,
        stats_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        color="white",
        bbox=dict(facecolor="black", alpha=0.6, edgecolor="white", boxstyle="round,pad=0.3"),
    )

    ax.set_title(f"{iso2} {tech} annual mean CF\nP95 threshold = {p95_threshold:.3f}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.legend()
    fig.colorbar(im, ax=ax, label="Annual mean CF")

    out = OUTDIR / f"{iso2.lower()}_{tech}_bestsite_new.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")

def main():
    for iso2 in COUNTRIES:
        for tech in TECHS:
            plot_bestsite_map(iso2, tech)


if __name__ == "__main__":
    main()