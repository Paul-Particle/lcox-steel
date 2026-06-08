"""Plot mean annual CF as a spatial heatmap for one cutout area and tech.

Draws the region boundary as a white outline and marks the P95 best-site cell.
Snakemake rule: plot_cf_map (in viz.smk).
"""

import logging
from pathlib import Path

import atlite
import geopandas as gpd
import numpy as np
import plotly.graph_objects as go
from shapely.geometry import MultiPolygon, Polygon

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.viz._helpers import (
    PLOTLY_CONFIG,
    blue_black,
    fca_colormap,
    fca_template,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

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
    else:
        raise ValueError(f"Unsupported tech: {_TECH!r}")

    cell_mean = cf_grid.mean("time")

    gdf = gpd.read_parquet(_REGIONS_PATH).to_crs(4326)
    geom = gdf.geometry.iloc[0]

    y_idx, x_idx = _find_p95_cell(cf_grid, geom)
    p95_lat = float(cutout.data.coords["y"].isel(y=y_idx))
    p95_lon = float(cutout.data.coords["x"].isel(x=x_idx))
    p95_cf  = float(cell_mean.isel(y=y_idx, x=x_idx))
    log.info(f"P95 cell: lat={p95_lat:.2f} lon={p95_lon:.2f} cf={p95_cf:.3f}")

    lons = cutout.data.coords["x"].values
    lats = cutout.data.coords["y"].values

    bx, by = _boundary_segments(geom)

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        x=lons,
        y=lats,
        z=cell_mean.values,
        colorscale=fca_colormap,
        zmin=0,
        colorbar=dict(
            title=dict(text="Mean annual CF<br>(latitude-optimal)", side="right"),
            thickness=14,
            len=0.8,
        ),
        hovertemplate="lon %{x:.2f}, lat %{y:.2f}<br>CF %{z:.3f}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=bx, y=by,
        mode="lines",
        line=dict(color="white", width=1.6),
        hoverinfo="skip",
        showlegend=False,
        name="region",
    ))
    fig.add_trace(go.Scatter(
        x=[p95_lon], y=[p95_lat],
        mode="markers",
        marker=dict(symbol="star", size=18, color="white",
                    line=dict(color=blue_black, width=1.2)),
        name=f"P95 site (CF={p95_cf:.3f})",
        hovertemplate=f"P95 site<br>lon %{{x:.2f}}, lat %{{y:.2f}}<br>CF {p95_cf:.3f}<extra></extra>",
    ))

    fig.update_layout(
        template=fca_template,
        title=(f"{_REGION} — {_TECH} mean annual CF ({_START_DATE[:4]})<br>"
               "<sup>latitude-optimal orientation · P95 bestsite marked</sup>"),
        xaxis_title="Longitude",
        yaxis_title="Latitude",
        width=900,
        height=720,
        legend=dict(x=0.99, y=0.01, xanchor="right", yanchor="bottom",
                    bgcolor="rgba(255,255,255,0.7)"),
        margin=dict(l=70, r=120, t=90, b=70),
    )
    # Equal-aspect lat/lon. scaleratio≈1 is fine away from poles; at mid-
    # latitudes longitude visually stretches but matches the matplotlib version.
    fig.update_yaxes(scaleanchor="x", scaleratio=1)

    fig.write_image(_OUT, scale=2)
    fig.write_html(_OUT.with_suffix(".html"), config=PLOTLY_CONFIG, include_plotlyjs="cdn")
    log.info(f"saved {_OUT} (+ .html)")


if __name__ == "__main__":
    main()
