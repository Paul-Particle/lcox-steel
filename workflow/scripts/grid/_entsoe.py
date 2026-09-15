"""ENTSO-E source implementation, called by retrieve_grid_data.py.

Builds the requested window out of the per-month raw cache under
  data/entsoe_cache/{area}/{YYYY-MM}/{data_type}.parquet

A month absent from it is downloaded once (via the download_entsoe primitives);
every month the window touches is then processed and sliced. The rule output is
itself the processed cache — it carries the whole run in its name, so Snakemake
skips this work on a re-run.

Variants
--------
dayahead  prices data_type only → single "price" column, hourly UTC-naive
full      all six data_types → wide frame with derived residual-load columns
"""

import logging
from pathlib import Path

import pandas as pd

from _helpers_grid import assert_window_complete, iso, to_utc_naive
from download_entsoe import DOWNLOADERS, download_with_retry, get_entsoe_client, iter_months

# Module-level logger only — retrieve_grid_data.py installs the handlers.
log = logging.getLogger(__name__)

FULL_DATA_TYPES = ["prices", "load_forecast", "load_actual", "res", "generation", "crossborder"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_bidding_zones(raw_cache_dir: Path) -> set[str]:
    """Return the set of recognised bidding-zone codes from the registry CSV."""
    return set(pd.read_csv(raw_cache_dir / "entsoe_bidding_zones.csv")["area"])


# ── Month selection ───────────────────────────────────────────────────────────

def _months_to_process(start_date: str, end_date: str) -> list[str]:
    """Calendar months ('YYYY-MM') covering the window, padded by one month.

    ENTSO-E raw months are fetched on Brussels-time boundaries and stored
    UTC-naive. In summer (CEST = UTC+2) a Brussels month's data ends at 22:00 UTC
    on its last day, so the requested window's final UTC hour(s) live in the *next*
    Brussels month. Padding by one month makes the window slice complete; without
    it `assert_window_complete` fails every Mar–Sep-ending window. (Mirrors the pad
    in _nem.py.)
    """
    months = [ym for ym, _, _ in iter_months(start_date, end_date)]
    pad_month = (pd.Timestamp(iso(end_date)) + pd.offsets.MonthBegin(1)).strftime("%Y-%m")
    if pad_month not in months:
        months = [*months, pad_month]
    return months


# ── Raw-cache management ──────────────────────────────────────────────────────

def _ensure_raw_months(
    area: str, months: list[str], data_types: list[str], raw_cache_dir: Path
) -> None:
    """Download any (month, data_type) pairs absent from the raw monthly cache."""
    area_cache_dir = raw_cache_dir / area
    client = None  # lazily instantiated — warm cache = no API call
    for ym in months:
        for dt in data_types:
            cache_path = area_cache_dir / ym / f"{dt}.parquet"
            if cache_path.exists():
                continue
            if client is None:
                client = get_entsoe_client()
            month_start = pd.Timestamp(ym + "-01", tz="Europe/Brussels")
            next_month = month_start + pd.offsets.MonthBegin(1)
            log.info(f"{area}/{ym}/{dt}: fetching")
            df = download_with_retry(DOWNLOADERS[dt], client, area, month_start, next_month)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(cache_path, index=True)
            log.info(f"{area}/{ym}/{dt}: cached")


# ── Per-variant month processing ──────────────────────────────────────────────

def _process_dayahead_month(area: str, ym: str, raw_cache_dir: Path) -> pd.DataFrame:
    """Return one month of prices as a single hourly, UTC-naive `price` column."""
    raw = pd.read_parquet(raw_cache_dir / area / ym / "prices.parquet")
    price = to_utc_naive(raw).iloc[:, 0]
    price = price.ffill(limit=3).bfill(limit=3).resample("1h").mean().rename("price")
    return price.to_frame()


def _process_full_month(area: str, ym: str, raw_cache_dir: Path) -> pd.DataFrame:
    """Assemble one month's six data types into a wide frame with residual-load columns.

    Joins prices, load (forecast + actual), RES forecast, generation, and
    cross-border flows on a common UTC-naive index, then derives the wind/RES
    aggregates and residual load (load − RES) for both forecast and actual.
    """
    d = raw_cache_dir / area / ym
    prices_raw  = pd.read_parquet(d / "prices.parquet")
    load_fc_raw = pd.read_parquet(d / "load_forecast.parquet")
    load_raw    = pd.read_parquet(d / "load_actual.parquet")
    res_raw     = pd.read_parquet(d / "res.parquet")
    gen_raw     = pd.read_parquet(d / "generation.parquet")
    xb_raw      = pd.read_parquet(d / "crossborder.parquet")

    price   = to_utc_naive(prices_raw).iloc[:, 0].rename("price")
    load_fc = to_utc_naive(load_fc_raw).iloc[:, 0].rename("load_forecast")
    load    = to_utc_naive(load_raw).iloc[:, 0].rename("load")
    res     = to_utc_naive(res_raw.copy());  res.columns = res.columns.droplevel(0)
    gen     = to_utc_naive(gen_raw.copy());  gen.columns = gen.columns.droplevel(0)
    xb      = to_utc_naive(xb_raw.copy());  xb.columns  = xb.columns.droplevel(0)

    df = pd.concat([price, load_fc, load, res, gen, xb], axis=1, sort=False)
    # Sort by time first: concat leaves the hourly series (price, load_forecast)
    # as a leading block, so their sub-hourly slots (:15/:30/:45) are not adjacent
    # to the :00 value. ffill only reaches them once the index is time-ordered —
    # otherwise those slots fall through to fillna(0.0) and dilute hourly means.
    df = df.sort_index()
    df = df.ffill(limit=3).fillna(0.0)

    df["wind_forecast"]     = df.get("wind_onshore_forecast", 0) + df.get("wind_offshore_forecast", 0)
    df["res_forecast"]      = df["wind_forecast"] + df.get("solar_forecast", 0)
    df["residual_forecast"] = df["load_forecast"] - df["res_forecast"]
    df["wind"]              = df.get("wind_onshore", 0) + df.get("wind_offshore", 0)
    df["res"]               = df["wind"] + df.get("solar", 0)
    df["residual"]          = df["load"] - df["res"]

    return df.sort_index(axis=1)


# ── Completeness guard ─────────────────────────────────────────────────────────

# ── Main ──────────────────────────────────────────────────────────────────────

def retrieve(snakemake, area: str) -> None:
    """Build the requested (area, variant, date range) window and write it out.

    `area` is the code this market knows the area by — the bidding zone or NEM
    region — which is not always the area code the rest of the workflow uses.

    Validates the bidding zone, downloads any month missing from the raw cache,
    processes every month the window touches, and writes the slice to the rule
    output. A warm raw cache means no API calls.
    """
    variant    = snakemake.wildcards.variant
    start_date = snakemake.wildcards.start_date
    end_date   = snakemake.wildcards.end_date

    raw_cache_dir = Path("data/entsoe_cache")

    if variant not in ("dayahead", "full"):
        raise ValueError(f"Unknown variant {variant!r}. Expected 'dayahead' or 'full'.")

    valid_zones = _load_bidding_zones(raw_cache_dir)
    if area not in valid_zones:
        raise ValueError(
            f"{area!r} is not a recognised ENTSO-E bidding zone. "
            f"See data/entsoe_cache/entsoe_bidding_zones.csv for the full list."
        )

    data_types    = ["prices"] if variant == "dayahead" else FULL_DATA_TYPES
    process_month = _process_dayahead_month if variant == "dayahead" else _process_full_month

    months = _months_to_process(start_date, end_date)
    _ensure_raw_months(area, months, data_types, raw_cache_dir)

    monthly_frames = [process_month(area, ym, raw_cache_dir) for ym in months]
    # Brussels months meet ~1-2 h off the UTC hour, so consecutive months can each
    # claim the hour they straddle. Frames are chronological, so keeping the last
    # occurrence gives that hour to the month whose calendar owns it.
    assembled = pd.concat(monthly_frames)
    assembled = assembled[~assembled.index.duplicated(keep="last")].sort_index()

    window = slice(iso(start_date), f"{iso(end_date)} 23:00")
    out_df = assembled.loc[window]
    # Cross-month boundary gaps: forward-fill only (bfill would propagate future data backward).
    # Per-month processing already handles within-month gaps including start-of-month.
    out_df = out_df.ffill(limit=3)
    if variant == "full":
        out_df = out_df.fillna(0.0)
    out_df.index.name = "time"

    assert_window_complete(out_df, start_date, end_date, variant)

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")
