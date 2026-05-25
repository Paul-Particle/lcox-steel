"""
Viz 1: wind and solar annual-mean capacity factor over the VIC bounding box,
as side-by-side plotly heatmaps. Reads results/vic_2025/cf_grids.nc.

Output: results/vic_2025/cf_maps.html (and PNG via save_fig).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common._paths import RESULTS  # noqa: E402
from utils import save_fig, fca_template  # noqa: E402


def crop_to_vic(da: xr.DataArray, in_vic: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (cropped_values, lon_crop, lat_crop) tight to the VIC bbox."""
    arr = da.values.copy()
    arr[~in_vic] = np.nan
    rows = np.any(in_vic, axis=1)
    cols = np.any(in_vic, axis=0)
    y0, y1 = np.where(rows)[0][[0, -1]]
    x0, x1 = np.where(cols)[0][[0, -1]]
    return (
        arr[y0:y1 + 1, x0:x1 + 1],
        da.x.values[x0:x1 + 1],
        da.y.values[y0:y1 + 1],
    )


def main() -> None:
    in_dir = RESULTS / "vic_2025"
    grids = xr.open_dataset(in_dir / "cf_grids.nc")
    in_vic = grids["in_vic"].values.astype(bool)

    wind_arr,  lon, lat = crop_to_vic(grids["wind_cf_annual_mean"],  in_vic)
    solar_arr, _,   _   = crop_to_vic(grids["solar_cf_annual_mean"], in_vic)

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Wind onshore — annual mean CF", "Solar — annual mean CF"))
    fig.add_trace(go.Heatmap(z=wind_arr,  x=lon, y=lat, colorscale="Viridis", colorbar=dict(title="CF", x=0.46)), row=1, col=1)
    fig.add_trace(go.Heatmap(z=solar_arr, x=lon, y=lat, colorscale="Inferno", colorbar=dict(title="CF", x=1.0)),  row=1, col=2)

    fig.update_layout(
        template=fca_template,
        title=f"VIC 2025 — per-cell annual mean capacity factor",
        height=500, width=1200,
    )
    fig.update_xaxes(title="lon", row=1, col=1)
    fig.update_xaxes(title="lon", row=1, col=2)
    fig.update_yaxes(title="lat", row=1, col=1)

    save_fig(fig, filename="cf_maps", path=in_dir, add_time=False, html=True, png=False)
    print(f"Wrote {in_dir / 'cf_maps.html'}")


if __name__ == "__main__":
    main()
