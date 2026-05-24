"""
Viz 2: per-location CF maps using the optimised wind/solar mix.

For each of (4 locations) × (captive, hybrid) = 8 panels, plot per-cell
mix-weighted CF =
    (wind_cap_opt * cf_wind_cell + solar_cap_opt * cf_solar_cell)
    / (wind_cap_opt + solar_cap_opt)
clipped to the VIC bbox. Overlay all 4 location markers on every panel; the
panel's own location is ringed.

Output: results/vic_2025/location_maps.html
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common._paths import RESULTS  # noqa: E402
from utils import save_fig, fca_template  # noqa: E402


LOCATION_ORDER = ["wind_p95", "solar_p95", "random_1", "random_2"]
MODES = ["captive", "hybrid"]


def crop_to_vic(arr: np.ndarray, in_vic: np.ndarray, lon: np.ndarray, lat: np.ndarray):
    arr = arr.copy()
    arr[~in_vic] = np.nan
    rows = np.any(in_vic, axis=1)
    cols = np.any(in_vic, axis=0)
    y0, y1 = np.where(rows)[0][[0, -1]]
    x0, x1 = np.where(cols)[0][[0, -1]]
    return arr[y0:y1 + 1, x0:x1 + 1], lon[x0:x1 + 1], lat[y0:y1 + 1]


def main() -> None:
    in_dir = RESULTS / "vic_2025"

    grids     = xr.open_dataset(in_dir / "cf_grids.nc")
    capacities = pd.read_csv(in_dir / "capacities.csv")
    locations = pd.read_csv(in_dir / "locations.csv").set_index("location")

    in_vic = grids["in_vic"].values.astype(bool)
    wind_grid  = grids["wind_cf_annual_mean"].values
    solar_grid = grids["solar_cf_annual_mean"].values
    lon = grids.x.values
    lat = grids.y.values

    fig = make_subplots(
        rows=len(LOCATION_ORDER), cols=len(MODES),
        subplot_titles=[f"{loc} — {mode}" for loc in LOCATION_ORDER for mode in MODES],
        horizontal_spacing=0.08, vertical_spacing=0.05,
    )

    for i, loc in enumerate(LOCATION_ORDER, start=1):
        for j, mode in enumerate(MODES, start=1):
            sel = capacities[(capacities["location"] == loc) & (capacities["scenario"] == mode)]
            if sel.empty:
                continue
            w = float(sel["wind_onshore_mw_opt"].iloc[0] or 0.0)
            s = float(sel["solar_mw_opt"].iloc[0] or 0.0)
            if w + s <= 0:
                mix = np.full_like(wind_grid, np.nan)
            else:
                mix = (w * wind_grid + s * solar_grid) / (w + s)

            mix_crop, lon_crop, lat_crop = crop_to_vic(mix, in_vic, lon, lat)
            fig.add_trace(
                go.Heatmap(z=mix_crop, x=lon_crop, y=lat_crop, colorscale="Viridis",
                           zmin=0, zmax=float(np.nanmax([wind_grid[in_vic].max(), solar_grid[in_vic].max()])),
                           showscale=(i == 1 and j == 2),
                           colorbar=dict(title="mix CF") if (i == 1 and j == 2) else None),
                row=i, col=j,
            )

            # Overlay all 4 location markers
            for name in LOCATION_ORDER:
                lon_l, lat_l = locations.at[name, "lon"], locations.at[name, "lat"]
                fig.add_trace(
                    go.Scatter(
                        x=[lon_l], y=[lat_l],
                        mode="markers",
                        marker=dict(
                            size=14 if name == loc else 8,
                            symbol="circle-open" if name == loc else "x",
                            line=dict(width=2 if name == loc else 1, color="red" if name == loc else "white"),
                            color="red" if name == loc else "white",
                        ),
                        name=name,
                        showlegend=(i == 1 and j == 1),
                    ),
                    row=i, col=j,
                )

    fig.update_layout(
        template=fca_template,
        title="VIC 2025 — per-cell mix CF at each location's optimised wind/solar capacities",
        height=300 * len(LOCATION_ORDER), width=1100,
    )

    save_fig(fig, filename="location_maps", path=in_dir, add_time=False, html=True, png=False)
    print(f"Wrote {in_dir / 'location_maps.html'}")


if __name__ == "__main__":
    main()
