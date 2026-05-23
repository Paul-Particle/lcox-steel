"""
Compute and plot wind CF for a single hour across a spatial subset of the cutout.

Useful for inspecting per-cell CF values — including flat-top saturation at CF=1
that is invisible in the nationally-aggregated time series.

Run from project root:
    python scripts/res_cf/diag_single_hour_cf.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import atlite
import geopandas as gpd
import plotly.graph_objects as go

# --- configure here ---
CUTOUT_PATH = str(CUTOUTS / "de_2025_oct2w.nc")
TURBINE     = "Vestas_V112_3MW"
TARGET_HOUR = "2025-10-26 12:00"
LAT_MIN     = 46.0
LAT_MAX     = 56.0
# ----------------------

REGIONS_PATH = SHAPES_RES / "regions.geojson"


def main():
    print(f"Loading cutout: {CUTOUT_PATH}")
    cutout = atlite.Cutout(CUTOUT_PATH)

    cutout.data = cutout.data.sel(
        time=[pd.Timestamp(TARGET_HOUR)],
        y=slice(LAT_MIN, LAT_MAX),
    )
    print(f"Subset: {len(cutout.data.y)} lat cells × {len(cutout.data.x)} lon cells")

    print("Computing CF...")
    cf = cutout.wind(turbine=TURBINE, aggregate_time="mean")

    vals = cf.values.ravel()
    vals_finite = vals[np.isfinite(vals)]
    print(f"Max CF:  {vals_finite.max():.4f}")
    print(f"Mean CF: {vals_finite.mean():.4f}")
    print(f"Top 5:   {np.sort(vals_finite)[::-1][:5].round(4).tolist()}")

    idx = cf.argmax(dim=["y", "x"])
    best_lat = float(cf.y[idx["y"]])
    best_lon = float(cf.x[idx["x"]])
    print(f"Best cell: lat={best_lat:.2f}, lon={best_lon:.2f}")

    # --- plot ---
    xs = cf.coords["x"].values
    ys = cf.coords["y"].values

    # country outline for context
    shapes = []
    if REGIONS_PATH.exists():
        regions = gpd.read_file(REGIONS_PATH)
        de = regions.loc[regions["region"] == "DE", "geometry"].iloc[0]
        if de.geom_type == "Polygon":
            coords = [de.exterior.coords]
        else:
            coords = [p.exterior.coords for p in de.geoms]
        for ring in coords:
            rx, ry = zip(*ring)
            shapes.append(go.Scattergeo(
                lon=list(rx), lat=list(ry),
                mode="lines",
                line=dict(color="white", width=1.5),
                showlegend=False,
            ))

    # CF heatmap
    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=cf.values,
        x=xs,
        y=ys,
        colorscale="YlOrRd",
        zmin=0, zmax=1,
        colorbar=dict(title="CF"),
    ))

    # country outline on top
    for shape in shapes:
        fig.add_trace(go.Scatter(
            x=shape.lon, y=shape.lat,
            mode="lines",
            line=dict(color="white", width=1.5),
            showlegend=False,
        ))

    # best cell marker
    fig.add_trace(go.Scatter(
        x=[best_lon], y=[best_lat],
        mode="markers",
        marker=dict(symbol="x", size=12, color="white", line=dict(width=2)),
        name=f"Best cell (CF={vals_finite.max():.3f})",
    ))

    fig.update_layout(
        title=f"Wind CF — {TARGET_HOUR} | {TURBINE}<br>lat {LAT_MIN}–{LAT_MAX}°N",
        xaxis_title="Longitude",
        yaxis_title="Latitude",
        yaxis=dict(scaleanchor="x", scaleratio=1),
        height=500,
        hovermode="closest",
    )

    fig.show()


if __name__ == "__main__":
    main()
