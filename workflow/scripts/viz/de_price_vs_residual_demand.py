#!/usr/bin/env python3
"""Diagnostic: German day-ahead price vs. residual demand (load − wind − solar), 2025.

Residual demand = actual load − wind forecast (on+offshore) − solar forecast.
Data: ENTSO-E DE-LU bidding zone, hourly, full year 2025.

Usage (from project root):
    python workflow/scripts/viz/de_price_vs_residual_demand.py

Output: results/de_price_vs_residual_demand_2025.html

Missing cache months are auto-fetched via the ENTSO-E API (reads ENTSOE_API_KEY
from .env).
"""
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "workflow" / "scripts" / "grid"))
sys.path.insert(0, str(ROOT / "workflow"))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dotenv import load_dotenv

from _helpers import to_utc_naive
from download_entsoe import DOWNLOADERS, download_with_retry, get_entsoe_client
from scripts.viz.style import (
    apply_header,
    blue_black,
    fca_template,
    highlight_blue,
    save_figure,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

AREA = "DE_LU"
YEAR = 2025
CACHE = ROOT / "data" / "entsoe_cache" / AREA
OUT = ROOT / "results" / "de_price_vs_residual_demand_2025.html"

_DATA_TYPES = ["prices", "load_actual", "res"]
BIN_WIDTH_MW = 500
BIN_WIDTH_PRES_MW = 2000   # wider bins → smoother average line for presentation
MIN_OBS = 20  # hide bins with fewer observations (too noisy)


def _ensure_months(months: list[str]) -> None:
    """Download any missing (month, data_type) combos into the raw cache."""
    client = None
    for ym in months:
        for dt in _DATA_TYPES:
            path = CACHE / ym / f"{dt}.parquet"
            if path.exists():
                continue
            if client is None:
                load_dotenv(ROOT / ".env")
                client = get_entsoe_client()
            start = pd.Timestamp(f"{ym}-01", tz="Europe/Brussels")
            end = start + pd.offsets.MonthBegin(1)
            log.info(f"fetching {AREA}/{ym}/{dt}")
            df = download_with_retry(DOWNLOADERS[dt], client, AREA, start, end)
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(path)
            log.info(f"cached  {AREA}/{ym}/{dt}")


def load_year(year: int) -> pd.DataFrame:
    months = [ts.strftime("%Y-%m") for ts in pd.date_range(f"{year}-01-01", f"{year}-12-01", freq="MS")]
    _ensure_months(months)

    parts = []
    for ym in months:
        d = CACHE / ym
        prices_raw = pd.read_parquet(d / "prices.parquet")
        load_raw = pd.read_parquet(d / "load_actual.parquet")
        res_raw = pd.read_parquet(d / "res.parquet")

        # Prices, load and RES may be hourly or 15-min depending on the year
        # (ENTSO-E switched to 15-min granularity in 2025).  Resample all to
        # hourly mean so everything aligns cleanly.
        price = to_utc_naive(prices_raw).iloc[:, 0].resample("1h").mean().rename("price")
        load = to_utc_naive(load_raw).iloc[:, 0].resample("1h").mean().rename("load")
        res = to_utc_naive(res_raw.copy())
        res.columns = res.columns.droplevel(0)
        res_h = res.resample("1h").mean()

        df = pd.concat([price, load, res_h], axis=1, sort=False).ffill(limit=3)
        wind = df.get("wind_onshore_forecast", 0) + df.get("wind_offshore_forecast", 0)
        df["residual_demand"] = df["load"] - wind - df.get("solar_forecast", 0)
        parts.append(df[["price", "residual_demand"]].dropna())

    return pd.concat(parts).sort_index()


def plot(df: pd.DataFrame, out: Path) -> None:
    x = df["residual_demand"]
    y = df["price"]
    agg = _binned_mean(x, y, BIN_WIDTH_MW)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x, y=y,
        mode="markers",
        marker=dict(size=2.5, color="rgba(100,100,100,0.10)"),
        showlegend=False,
        hovertemplate="Residual: %{x:,.0f} MW<br>Price: %{y:.1f} €/MWh<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=agg.index, y=agg["mean"],
        mode="lines",
        line=dict(color=highlight_blue, width=2.5),
        name=f"bin avg (n ≥ {MIN_OBS}, {BIN_WIDTH_MW/1000:.0f} GW bins)",
        hovertemplate="Residual: %{x:,.0f} MW<br>Avg price: %{y:.1f} €/MWh<extra></extra>",
    ))

    fig.update_layout(
        template=fca_template,
        xaxis=dict(title="Residual demand (MW)", range=[-19_000, 69_000]),
        # No rotated y-axis title — the price quantity/unit lives in the subtitle.
        yaxis=dict(title=None, range=[-100, 250]),
        legend=dict(x=0.02, y=0.97, xanchor="left", yanchor="top",
                    bgcolor="rgba(255,255,255,0.65)"),
    )

    apply_header(
        fig,
        title="DE day-ahead price vs. residual demand — 2025",
        subtitle="Day-ahead price (€/MWh)  ·  Residual = load − wind (on/offshore) − solar  ·  DE-LU  ·  hourly",
        fig_width=960, fig_height=600,
        margin_l=80, margin_r=60, margin_t=110, margin_b=80,
    )

    saved = save_figure(fig, out.parent, out.stem, scale=4)   # 960×4 = 3840 px wide PNG
    fig.write_image(out.with_suffix(".svg"))
    log.info(f"saved {' + '.join(saved)} + {out.with_suffix('.svg')}")


def _binned_mean(x, y, bin_width: int) -> pd.Series:
    bins = np.arange(-20_000, 70_000 + bin_width, bin_width)
    mids = 0.5 * (bins[:-1] + bins[1:])
    idx = np.clip(np.digitize(x.values, bins) - 1, 0, len(mids) - 1)
    agg = (
        pd.DataFrame({"price": y.values, "bin": idx})
        .groupby("bin")["price"]
        .agg(["mean", "count"])
    )
    agg = agg[agg["count"] >= MIN_OBS]
    agg.index = mids[agg.index]
    return agg


def plot_presentation(df: pd.DataFrame, out: Path) -> None:
    x = df["residual_demand"]
    y = df["price"]
    agg = _binned_mean(x, y, BIN_WIDTH_PRES_MW)
    agg = agg[agg.index <= 48_000]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x, y=y,
        mode="markers",
        marker=dict(size=4, color="rgba(100,100,100,0.10)"),
        showlegend=False,
        hovertemplate="Residual: %{x:,.0f} MW<br>Price: %{y:.1f} €/MWh<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=agg.index, y=agg["mean"],
        mode="lines",
        line=dict(color=highlight_blue, width=4),
        name="Average",
        hovertemplate="Residual: %{x:,.0f} MW<br>Avg price: %{y:.1f} €/MWh<extra></extra>",
    ))

    fig.update_layout(
        template=fca_template,
        title=dict(
            text="Day-ahead price (€/MWh),<br>Germany, 2025",
            font=dict(family="Titillium Web", size=24, color=blue_black),
            x=0.01,
            xref="container",
            xanchor="left",
        ),
        xaxis=dict(
            title="Demand not supplied by<br>wind and solar (MW)",
            title_font=dict(size=20),
            tickfont=dict(size=20),
            dtick=10_000,
            range=[-9_000, 49_000],
        ),
        yaxis=dict(
            title=None,
            tickfont=dict(size=20),
            dtick=50,
            range=[-40, 150],
        ),
        width=393,    # 5.8" at 67.8 DPI (PowerPoint default for unit-less PNG)
        height=384,
        legend=dict(x=0.98, y=0.02, xanchor="right", yanchor="bottom", font=dict(size=18)),
        margin=dict(l=60, r=30, t=90, b=100),
    )

    # Bespoke PowerPoint-sized slide asset: keep its custom title/sizing, but
    # save it the house way (self-contained, font-embedded HTML + retina PNG).
    saved = save_figure(fig, out.parent, out.stem, scale=4)   # 393×4 = 1572 px → 5.8" at 72 DPI
    fig.write_image(out.with_suffix(".svg"))
    log.info(f"saved {' + '.join(saved)} + {out.with_suffix('.svg')}")


if __name__ == "__main__":
    df = load_year(YEAR)
    log.info(f"loaded {len(df):,} hourly observations for {YEAR}")
    plot(df, OUT)
    plot_presentation(df, OUT.with_name(OUT.stem + "_pres" + OUT.suffix))
