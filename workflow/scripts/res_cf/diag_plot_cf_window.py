"""
Quick diagnostic: plot a short window of hourly wind CF for a given country.

Useful for inspecting:
- flat tops at CF=1 (rated wind speed regime)
- cut-out zeros (extreme winds, turbine shutdown)
- overall profile shape

Run from project root:
    python scripts/res_cf/diag_plot_cf_window.py
"""

import pandas as pd
import plotly.graph_objects as go
from common._paths import RES_CF

# --- configure here ---
COUNTRY   = "de"
VARIANT   = "avg"            # "avg" or "bestsite_p95"
YEAR      = 2025
START     = "2025-10-20"
END       = "2025-11-02"
# ----------------------

DATA_DIR = RES_CF


def load_cf(country: str, variant: str, year: int) -> pd.DataFrame:
    if variant == "bestsite_p95":
        path = DATA_DIR / f"{country}_cf_{year}_bestsite_p95.csv"
    elif variant == "avg":
        path = DATA_DIR / f"{country}_cf_{year}.csv"
    else:
        raise ValueError(f"Unknown variant: {variant!r}")

    if not path.exists():
        raise FileNotFoundError(
            f"No CF file found at {path}\n"
            "Run scripts 01–07 first to generate the data."
        )

    df = pd.read_csv(path, parse_dates=["time"])
    df = df.sort_values("time").set_index("time")
    return df


def main():
    df = load_cf(COUNTRY, VARIANT, YEAR)
    window = df.loc[START:END]

    if window.empty:
        raise ValueError(f"No data in window {START}–{END}. Check COUNTRY/YEAR.")

    tech_cols = [c for c in window.columns if c.endswith("_cf")]

    fig = go.Figure()

    for col in tech_cols:
        label = col.replace("_cf", "").replace("_", " ")
        fig.add_trace(go.Scatter(
            x=window.index,
            y=window[col],
            mode="lines",
            name=label,
            line=dict(width=1.2),
        ))

    fig.add_hline(
        y=1.0,
        line=dict(color="grey", width=1, dash="dot"),
        annotation_text="CF = 1",
        annotation_position="top right",
    )

    fig.update_layout(
        title=f"{COUNTRY.upper()} {VARIANT} — {START} to {END}",
        xaxis_title="Time",
        yaxis_title="Capacity factor",
        yaxis=dict(range=[-0.02, 1.05]),
        hovermode="x unified",
        height=420,
    )

    fig.show()


if __name__ == "__main__":
    main()
