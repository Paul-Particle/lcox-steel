"""
Viz 3: VIC1 2025 hourly price distribution (EUR/MWh), with vertical lines
showing the max hourly price the model was still willing to import grid power
at, for both grid-only and hybrid scenarios. One figure per P95 cell
(wind_p95 + solar_p95).

Output: results/vic_2025/prices.html
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common._paths import RESULTS  # noqa: E402
from utils import save_fig, fca_template  # noqa: E402


P95_LOCATIONS = ["wind_p95", "solar_p95"]
GRID_SCENARIOS = ["grid_only", "hybrid"]


def main() -> None:
    in_dir = RESULTS / "vic_2025"

    prices = pd.read_csv(in_dir / "price_series.csv", index_col=0, parse_dates=True)
    price_series = prices["price_eur_per_mwh"]
    caps = pd.read_csv(in_dir / "capacities.csv")

    fig = make_subplots(
        rows=len(P95_LOCATIONS), cols=1,
        subplot_titles=[f"VIC1 price distribution — max accepted at {loc}" for loc in P95_LOCATIONS],
        vertical_spacing=0.12,
    )

    for i, loc in enumerate(P95_LOCATIONS, start=1):
        fig.add_trace(
            go.Histogram(
                x=price_series.values,
                nbinsx=80,
                marker=dict(color="rgba(120,120,120,0.6)"),
                name=f"VIC1 hourly price",
                showlegend=(i == 1),
            ),
            row=i, col=1,
        )

        line_colors = {"grid_only": "tomato", "hybrid": "steelblue"}
        for sc in GRID_SCENARIOS:
            sel = caps[(caps["location"] == loc) & (caps["scenario"] == sc)]
            if sel.empty:
                continue
            max_p = sel["max_grid_import_price_eur_mwh"].iloc[0]
            if pd.isna(max_p):
                continue
            fig.add_vline(
                x=float(max_p),
                line=dict(color=line_colors[sc], width=2, dash="dash"),
                annotation_text=f"{sc}: {max_p:.1f} €/MWh",
                annotation_position="top",
                row=i, col=1,
            )

        fig.update_xaxes(title="price (EUR/MWh)", row=i, col=1)
        fig.update_yaxes(title="hours", row=i, col=1)

    fig.update_layout(
        template=fca_template,
        title="VIC 2025 prices — distribution + max accepted price for grid-connected scenarios",
        height=400 * len(P95_LOCATIONS), width=900,
        barmode="overlay",
    )

    save_fig(fig, filename="prices", path=in_dir, add_time=False, html=True, png=False)
    print(f"Wrote {in_dir / 'prices.html'}")


if __name__ == "__main__":
    main()
