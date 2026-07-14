"""Retrieve ENTSO-E grid data for a (area, variant, date-range) slice.

Maintains a persistent second-level processed cache at
  resources/entsoe/{variant}.parquet

Cache columns are a MultiIndex (area, metric) so all areas share one file per
variant. Accessing one area's data: df["DE_LU"]. On a warm-cache run the rule
just slices and writes the temp output; on a miss it downloads any absent months
to the per-month raw cache (via download_entsoe primitives), processes them, and
extends the cache before slicing.

Variants
--------
dayahead  prices data_type only → single "price" column, hourly UTC-naive
full      all six data_types → wide frame with derived residual-load columns
"""

import logging
from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from _helpers import area_month_in_cache, iso, to_utc_naive
from download_entsoe import DOWNLOADERS, download_with_retry, get_entsoe_client, iter_months

configure_logging(snakemake)
log = logging.getLogger(__name__)

FULL_DATA_TYPES = ["prices", "load_forecast", "load_actual", "res", "generation", "crossborder"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_bidding_zones(raw_cache_dir: Path) -> set[str]:
    """Return the set of recognised bidding-zone codes from the registry CSV."""
    return set(pd.read_csv(raw_cache_dir / "entsoe_bidding_zones.csv")["area"])


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

def _summarise_runs(times: pd.DatetimeIndex, step: pd.Timedelta, limit: int = 3) -> str:
    """Compress a sorted DatetimeIndex into 'start→end' run descriptions, largest first."""
    times = times.sort_values()
    runs = []
    start = prev = times[0]
    for t in times[1:]:
        if t - prev == step:
            prev = t
        else:
            runs.append((start, prev))
            start = prev = t
    runs.append((start, prev))
    runs.sort(key=lambda ab: ab[1] - ab[0], reverse=True)
    parts = [f"{a}→{b}" for a, b in runs[:limit]]
    if len(runs) > limit:
        parts.append(f"…+{len(runs) - limit} more")
    return ", ".join(parts)


def _assert_window_complete(
    out_df: pd.DataFrame, start_date: str, end_date: str, variant: str
) -> None:
    """Fail loudly if the produced window has holes over the requested date range.

    Truncated raw-cache months (a fetch that stopped mid-month) and partial
    downloads are otherwise silent: the slice simply has fewer rows, and the hole
    only surfaces far downstream as NaN after a reindex. This guard catches them
    at the source.

    dayahead is resampled to a clean hourly grid, so every hour of the window must
    be present and non-null. full spans a mixed 15-min/hourly resolution (ENTSO-E
    switched DE_LU day-ahead to 15-min in Oct 2025), so a fixed grid would
    false-positive; instead we flag any gap larger than the ffill tolerance, which
    only occurs on truncated months.
    """
    idx = out_df.index
    want_start = pd.Timestamp(iso(start_date))
    want_end = pd.Timestamp(f"{iso(end_date)} 23:00")
    problems = []

    if len(idx) < 2:
        raise ValueError(f"{variant} {start_date}–{end_date}: only {len(idx)} rows produced")
    if idx.min() > want_start:
        problems.append(f"starts at {idx.min()}, after {want_start}")
    if idx.max() < want_end:
        problems.append(f"ends at {idx.max()}, before {want_end}")

    if variant == "dayahead":
        expected = pd.date_range(want_start, want_end, freq="h")
        missing = expected.difference(idx)
        if len(missing):
            problems.append(
                f"{len(missing)} missing hours: {_summarise_runs(missing, pd.Timedelta('1h'))}"
            )
        n_nan = int(out_df.isna().any(axis=1).sum())
        if n_nan:
            problems.append(f"{n_nan} rows still NaN after fill")
    else:
        gaps = idx.to_series().diff()
        if gaps.max() > pd.Timedelta("3h"):
            problems.append(f"gap of {gaps.max()} ending {gaps.idxmax()}")

    if problems:
        raise ValueError(f"{variant} {start_date}–{end_date} incomplete: " + "; ".join(problems))


# ── Main ──────────────────────────────────────────────────────────────────────

def retrieve(snakemake) -> None:
    """Slice the requested (area, variant, date range) out of the processed cache.

    Validates the bidding zone, downloads any missing months into the raw cache,
    processes and appends them to the shared per-variant processed cache, then
    writes the requested window to the rule output. Warm-cache runs make no API
    calls.
    """
    area       = snakemake.wildcards.area
    variant    = snakemake.wildcards.variant
    start_date = snakemake.wildcards.start_date
    end_date   = snakemake.wildcards.end_date

    raw_cache_dir        = Path("data/entsoe_cache")
    processed_cache_dir  = Path("resources/entsoe")
    processed_cache_path = processed_cache_dir / f"{variant}.parquet"

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

    cached = pd.read_parquet(processed_cache_path) if processed_cache_path.exists() else None

    months = [ym for ym, _, _ in iter_months(start_date, end_date)]

    _ensure_raw_months(area, months, data_types, raw_cache_dir)

    new_frames = []
    for ym in months:
        if area_month_in_cache(cached, area, ym):
            continue
        frame = process_month(area, ym, raw_cache_dir)
        frame.columns = pd.MultiIndex.from_tuples([(area, c) for c in frame.columns])
        new_frames.append(frame)

    if new_frames:
        all_frames = ([cached] if cached is not None else []) + new_frames
        cached = pd.concat(all_frames)
        cached = cached[~cached.index.duplicated(keep="last")].sort_index()
        processed_cache_dir.mkdir(parents=True, exist_ok=True)
        cached.to_parquet(processed_cache_path, index=True)
        log.info(f"updated processed cache: {processed_cache_path} ({len(cached)} rows)")

    window = slice(iso(start_date), f"{iso(end_date)} 23:00")
    out_df = cached[area].loc[window]
    # Cross-month boundary gaps: forward-fill only (bfill would propagate future data backward).
    # Per-month processing already handles within-month gaps including start-of-month.
    out_df = out_df.ffill(limit=3)
    if variant == "full":
        out_df = out_df.fillna(0.0)
    out_df.index.name = "time"

    _assert_window_complete(out_df, start_date, end_date, variant)

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")


if __name__ == "__main__":
    retrieve(snakemake)
