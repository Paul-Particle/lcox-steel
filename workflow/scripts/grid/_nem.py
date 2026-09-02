"""AEMO NEM source implementation, called by retrieve_grid_data.py.

Builds the requested window month by month out of the NEMOSIS raw cache under
data/nem_cache/, which NEMOSIS manages itself; no per-month raw files are
written here. The rule output is itself the processed cache — it carries the
whole run in its name, so Snakemake skips this work on a re-run.

Variants
--------
dayahead  price table only → single "price" column, hourly UTC-naive, EUR/MWh
full      all four tables → wide per-area frame with derived wind/residual columns
"""

import logging
from pathlib import Path

import pandas as pd

from _helpers import (
    assert_window_complete,
    iso,
    iter_months_str,
    to_utc_naive,
)
from download_nem import DOWNLOADERS

# Module-level logger only — retrieve_grid_data.py installs the handlers.
log = logging.getLogger(__name__)

FULL_TABLES = ["price", "load", "generation", "crossborder"]


def _month_range_str(ym: str) -> tuple[str, str]:
    """Return NEMOSIS-format (start, end) strings for a full calendar month."""
    ts = pd.Timestamp(ym + "-01")
    start = ts.strftime("%Y/%m/%d") + " 00:00:00"
    end   = (ts + pd.offsets.MonthEnd(0)).strftime("%Y/%m/%d") + " 23:59:59"
    return start, end


# ── Per-variant month processing ──────────────────────────────────────────────

def _process_dayahead_month(
    area: str, ym: str, cache_dir: Path, eur_per_aud: float
) -> pd.DataFrame:
    """Return one month of prices as a single hourly, UTC-naive `price` column (AUD→EUR)."""
    start_str, end_str = _month_range_str(ym)
    raw = DOWNLOADERS["price"](start_str, end_str, cache_dir, rebuild=False)
    raw = to_utc_naive(raw)
    price = (raw[(area, "price")].resample("1h").mean() * eur_per_aud).rename("price")
    return price.to_frame()


def _process_full_month(
    area: str, ym: str, cache_dir: Path, eur_per_aud: float
) -> pd.DataFrame:
    """Assemble one month's four tables into a wide per-area frame with derived columns.

    Concatenates the price, load, generation, and cross-border tables, then adds
    a `wind` aggregate and residual load (load − wind − solar) for the area.
    """
    start_str, end_str = _month_range_str(ym)
    tables = {
        t: to_utc_naive(DOWNLOADERS[t](start_str, end_str, cache_dir, rebuild=False))
        for t in FULL_TABLES
    }
    df = pd.concat(list(tables.values()), axis=1, sort=False)
    area_df = df[area].copy()

    wind  = area_df.get("wind_onshore", pd.Series(0, index=area_df.index))
    solar = area_df.get("solar",        pd.Series(0, index=area_df.index))
    area_df["wind"]     = wind
    area_df["residual"] = area_df["load"] - (wind + solar)

    return area_df


# ── Main ──────────────────────────────────────────────────────────────────────

def retrieve(snakemake, area: str) -> None:
    """Build the requested (area, variant, date range) window and write it out.

    `area` is the code this market knows the area by — the bidding zone or NEM
    region — which is not always the area code the rest of the workflow uses.

    Processes every market month the window touches (NEMOSIS manages its own raw
    download cache under data/nem_cache/) and writes the slice to the rule output.
    """
    variant     = snakemake.wildcards.variant
    start_date  = snakemake.wildcards.start_date
    end_date    = snakemake.wildcards.end_date
    eur_per_aud = snakemake.params.eur_per_aud

    cache_dir = Path("data/nem_cache")

    if variant not in ("dayahead", "full"):
        raise ValueError(f"Unknown variant {variant!r}. Expected 'dayahead' or 'full'.")

    months = iter_months_str(start_date, end_date)
    # NEM data is downloaded in market time (AEST, UTC+10) and converted to UTC,
    # which shifts each month ~10 h earlier — the requested window's final UTC hours
    # live in the *next* market-time month. Pad the download list by one month so
    # the UTC window slice below is complete.
    pad_month = (pd.Timestamp(iso(end_date)) + pd.offsets.MonthBegin(1)).strftime("%Y-%m")
    if pad_month not in months:
        months = [*months, pad_month]

    process_month = _process_dayahead_month if variant == "dayahead" else _process_full_month

    monthly_frames = []
    for ym in months:
        log.info(f"{area}/{ym}/{variant}: processing")
        monthly_frames.append(process_month(area, ym, cache_dir, eur_per_aud))

    # Market months are stored UTC, so consecutive ones meet ~10 h off the UTC month
    # boundary and can each claim the hours they straddle. Frames are chronological,
    # so keeping the last occurrence gives those hours to the month that owns them.
    assembled = pd.concat(monthly_frames)
    assembled = assembled[~assembled.index.duplicated(keep="last")].sort_index()

    window = slice(iso(start_date), f"{iso(end_date)} 23:59")
    out_df = assembled.loc[window]
    out_df.index.name = "time"

    assert_window_complete(out_df, start_date, end_date, variant)

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")
