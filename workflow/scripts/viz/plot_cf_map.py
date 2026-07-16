"""Plot mean annual CF as a spatial heatmap for one cutout area and tech.

Draws the region boundary as a white outline and marks the P95 best-site cell.
Snakemake rule: plot_cf_map (in viz.smk).
"""

import logging
from pathlib import Path

import atlite
import pandas as pd
import geopandas as gpd
import numpy as np
import plotly.graph_objects as go
from shapely.geometry import MultiPolygon, Polygon

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.viz.style import (
    apply_header,
    blue_black,
    fca_colormap,
    fca_template,
    save_figure,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

_CF_AREA = snakemake.wildcards.cf_area
_TECH = snakemake.wildcards.tech
_START_DATE = snakemake.wildcards.start_date
_END_DATE = snakemake.wildcards.end_date
_CUTOUT_PATH = Path(snakemake.input.cutout)
_REGIONS_PATH = Path(snakemake.input.regions)
_OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
_NATIONAL_MEAN_CF_PATH = Path(snakemake.input.national_mean_cf)
_REGION = snakemake.params.region
_PV_PANEL = snakemake.params.pv_panel
_WIND_TURBINE = snakemake.params.wind_onshore_turbine
_WIND_OFFSHORE_TURBINE = snakemake.params.wind_offshore_turbine
_OUT = Path(snakemake.output[0])


def _mask_cells_inside(cell_mean, geom) -> np.ndarray:
    """Return a boolean grid marking cutout cells whose centre lies within `geom`."""
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)


def _find_p95_cell(cf_grid, geom):
    """Return the (y, x) index of the in-region cell closest to the P95 mean CF."""
    cell_mean = cf_grid.mean("time")
    inside = _mask_cells_inside(cell_mean, geom)
    vals = np.where(inside, cell_mean.values, np.nan)
    valid = np.isfinite(vals)
    p95 = np.nanpercentile(vals[valid], 95)
    dist = np.abs(np.where(valid, vals, np.nan) - p95)
    idx_flat = np.nanargmin(dist)
    y_idx, x_idx = np.unravel_index(idx_flat, vals.shape)
    return int(y_idx), int(x_idx)

def _find_best_cell(cell_mean, geom):
    """Return the (y, x) index of the single highest-CF in-region cell."""
    inside = _mask_cells_inside(cell_mean, geom)
    vals = np.where(inside, cell_mean.values, np.nan)
    idx_flat = np.nanargmax(vals)
    y_idx, x_idx = np.unravel_index(idx_flat, vals.shape)
    return int(y_idx), int(x_idx)


def _grid_mean_unweighted(cell_mean, geom) -> float:
    """Plain (point-in-polygon, unweighted) mean CF over in-region cells.
    NOT the area-weighted national mean — see _national_mean_cf."""
    inside = _mask_cells_inside(cell_mean, geom)
    return float(np.nanmean(np.where(inside, cell_mean.values, np.nan)))


def _national_mean_cf(path: Path, tech_col: str) -> float:
    """Annual mean of the country-average CF export — the number the
    LCOE/LCOH pipeline actually consumes for the country-average scenario."""
    df = pd.read_parquet(path)
    return float(df[tech_col].mean())


def _boundary_segments(geom) -> tuple[list[float], list[float]]:
    """Flatten one or more polygon exteriors+holes into x/y arrays with NaN
    separators so a single Scatter trace can render all rings."""
    xs: list[float] = []
    ys: list[float] = []

    def _push(ring):
        rx, ry = zip(*ring.coords)
        xs.extend(rx)
        ys.extend(ry)
        xs.append(np.nan)
        ys.append(np.nan)

    polys = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    for poly in polys:
        if not isinstance(poly, Polygon):
            continue
        _push(poly.exterior)
        for interior in poly.interiors:
            _push(interior)
    return xs, ys


def main() -> None:
    """Build the CF grid, locate the P95 cell, and render the heatmap to PNG + HTML."""
    _OUT.parent.mkdir(parents=True, exist_ok=True)

    log.info(f"loading cutout {_CUTOUT_PATH}")
    cutout = atlite.Cutout(str(_CUTOUT_PATH))

    log.info(f"computing {_TECH} CF grid")
    if _TECH == "solar":
        cf_grid = cutout.pv(
            panel=_PV_PANEL,
            orientation="latitude_optimal",
            capacity_factor_timeseries=True,
        )
    elif _TECH == "wind-onshore":
        cf_grid = cutout.wind(
            turbine=_WIND_TURBINE,
            capacity_factor_timeseries=True,
            smooth=True,
            add_cutout_windspeed=True,
        )
    elif _TECH == "wind-offshore":
        cf_grid = cutout.wind(
            turbine=_WIND_OFFSHORE_TURBINE,
            capacity_factor_timeseries=True,
            smooth=True,
            add_cutout_windspeed=True,
        )
    else:
        raise ValueError(f"Unsupported tech: {_TECH!r}")

    cell_mean = cf_grid.mean("time")

    land_gdf = gpd.read_parquet(_REGIONS_PATH).to_crs(4326)
    land_geom = land_gdf.geometry.iloc[0]

    if _TECH == "wind-offshore":
        offshore_gdf = gpd.read_parquet(_OFFSHORE_REGIONS_PATH).to_crs(4326)
        geom = offshore_gdf.geometry.iloc[0]
    else:
        geom = land_geom

    if _TECH == "solar":
        tech_detail = "latitude-optimal orientation"
    elif _TECH == "wind-onshore":
        tech_detail = f"{_WIND_TURBINE} turbine"
    elif _TECH == "wind-offshore":
        tech_detail = f"{_WIND_OFFSHORE_TURBINE} turbine"
    else:
        tech_detail = _TECH

    y_idx, x_idx = _find_p95_cell(cf_grid, geom)
    p95_lat = float(cutout.data.coords["y"].isel(y=y_idx))
    p95_lon = float(cutout.data.coords["x"].isel(x=x_idx))
    p95_cf  = float(cell_mean.isel(y=y_idx, x=x_idx))
    log.info(f"P95 cell: lat={p95_lat:.2f} lon={p95_lon:.2f} cf={p95_cf:.3f}")
    by_idx, bx_idx = _find_best_cell(cell_mean, geom)
    best_lat = float(cutout.data.coords["y"].isel(y=by_idx))
    best_lon = float(cutout.data.coords["x"].isel(x=bx_idx))
    best_cf  = float(cell_mean.isel(y=by_idx, x=bx_idx))
    log.info(f"Best cell: lat={best_lat:.2f} lon={best_lon:.2f} cf={best_cf:.3f}")

    national_mean_cf = _national_mean_cf(_NATIONAL_MEAN_CF_PATH, _TECH)
    grid_mean_cf = _grid_mean_unweighted(cell_mean, geom)

    lons = cutout.data.coords["x"].values
    lats = cutout.data.coords["y"].values

    land_bx, land_by = _boundary_segments(land_geom)

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        x=lons,
        y=lats,
        z=cell_mean.values,
        colorscale=fca_colormap,
        zmin=0,
        colorbar=dict(
            title=dict(text=f"Mean annual CF<br>({tech_detail})", side="right"),
            thickness=14,
            len=0.8,
        ),
        hovertemplate="lon %{x:.2f}, lat %{y:.2f}<br>CF %{z:.3f}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=land_bx, y=land_by,
        mode="lines",
        line=dict(color="white", width=1.6),
        hoverinfo="skip",
        showlegend=False,
        name="land boundary",
    ))

    if _TECH == "wind-offshore":
        offshore_bx, offshore_by = _boundary_segments(geom)
        fig.add_trace(go.Scatter(
            x=offshore_bx, y=offshore_by,
            mode="lines",
            line=dict(color="white", width=1.2, dash="dash"),
            hoverinfo="skip",
            showlegend=False,
            name="offshore region",
        ))
    fig.add_trace(go.Scatter(
        x=[p95_lon], y=[p95_lat],
        mode="markers",
        marker=dict(symbol="star", size=18, color="white",
                    line=dict(color=blue_black, width=1.2)),
        name=f"P95 site (CF={p95_cf:.3f})",
        hovertemplate=f"P95 site<br>lon %{{x:.2f}}, lat %{{y:.2f}}<br>CF {p95_cf:.3f}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=[best_lon], y=[best_lat],
        mode="markers",
        marker=dict(symbol="x", size=14, color="white",
                    line=dict(color=blue_black, width=2.5)),
        name=f"Best site (CF={best_cf:.3f})",
        hovertemplate=f"Best site<br>lon %{{x:.2f}}, lat %{{y:.2f}}<br>CF {best_cf:.3f}<extra></extra>",
    ))

    stats_text = (
        f"National mean (area-weighted): {national_mean_cf:.5f}<br>"
        f"Grid mean (unweighted): {grid_mean_cf:.5f}<br>"
        f"P95 CF: {p95_cf:.3f}<br>"
        f"Best CF: {best_cf:.3f}"
    )
    fig.add_annotation(
        text=stats_text,
        xref="paper", yref="paper",
        x=0.02, y=0.98, xanchor="left", yanchor="top",
        showarrow=False, align="left",
        font=dict(family="Titillium Web", size=13, color="white"),
        bgcolor="rgba(51,67,77,0.75)",
        bordercolor="white", borderwidth=1, borderpad=6,
    )

    fig.update_layout(
        template=fca_template,
        # Degrees are self-evident from the ticks, so no axis titles (and no
        # sideways "Latitude" title — see the no-rotated-y-title house rule).
        xaxis_title=None,
        yaxis_title=None,
        # Legend bottom-left so it clears the bottom-right brand logo.
        legend=dict(x=0.01, y=0.01, xanchor="left", yanchor="bottom",
                    bgcolor="rgba(255,255,255,0.7)"),
    )
    # Equal-aspect lat/lon. scaleratio≈1 is fine away from poles; at mid-
    # latitudes longitude visually stretches but matches the matplotlib version.
    fig.update_yaxes(scaleanchor="x", scaleratio=1)

    apply_header(
        fig,
        title=f"{_REGION} — {_TECH} mean annual CF ({_START_DATE[:4]})",
        subtitle=f"{tech_detail} · P95 and best site marked",
        fig_width=900, fig_height=720,
        margin_l=70, margin_r=120, margin_t=90, margin_b=70,
        show_logo=False,   # the colorbar already occupies the bottom-right corner
    )
    saved = save_figure(fig, _OUT.parent, _OUT.stem)
    log.info(f"saved {' + '.join(saved)}")


if __name__ == "__main__":
    main()
