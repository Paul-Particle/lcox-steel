"""Shared helpers for _entsoe.py and _nem.py."""

import pandas as pd

# NEM market time is AEST (UTC+10, no daylight saving). "Australia/Brisbane" is
# the fixed-offset zone that models it — used both to convert raw downloads to UTC
# and to test market-month cache membership (see area_month_in_cache).
NEM_MARKET_TZ = "Australia/Brisbane"


def iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def iter_months_str(start_date: str, end_date: str) -> list[str]:
    """Return a list of 'YYYY-MM' strings for every month in [start_date, end_date].

    start_date / end_date are 'YYYYMMDD' strings (Snakemake wildcard format).
    """
    start = pd.Timestamp(iso(start_date))
    end = pd.Timestamp(iso(end_date))
    return [ts.strftime("%Y-%m") for ts in pd.date_range(start=start, end=end, freq="MS")]


def to_utc_naive(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a tz-aware or fixed-offset-naive (AEST = UTC+10) index to UTC-naive."""
    if df.index.tz is not None:
        df.index = df.index.tz_convert("UTC").tz_localize(None)
    else:
        df.index = df.index.tz_localize(NEM_MARKET_TZ).tz_convert("UTC").tz_localize(None)
    return df.sort_index()


def area_month_in_cache(
    cached: pd.DataFrame | None, area: str, ym: str, tz: str | None = None
) -> bool:
    """Return True if `cached` already holds data for (area, calendar month `ym`).

    By default the cached (UTC-naive) index is matched against `ym` directly. Pass
    `tz` to interpret the index in that timezone first. NEM stores prices in UTC
    but downloads them by *market* month (AEST), and the ~10 h shift spills each
    market month across two UTC months. A plain UTC-month match then reports a
    market month as "cached" on the strength of a neighbouring month's spillover
    (e.g. the trailing hours the pad month contributes to the prior UTC month), so
    that month gets skipped and never downloaded. Matching in market time
    (tz=NEM_MARKET_TZ) attributes every hour to its own market month and avoids this.
    """
    if cached is None or area not in cached.columns.get_level_values(0):
        return False
    index = cached[area].dropna(how="all").index
    if tz is not None:
        index = index.tz_localize("UTC").tz_convert(tz).tz_localize(None)
    return bool((index.to_period("M") == pd.Period(ym, freq="M")).any())


def summarise_runs(times: pd.DatetimeIndex, step: pd.Timedelta, limit: int = 3) -> str:
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


def assert_window_complete(
    out_df: pd.DataFrame,
    start_date: str,
    end_date: str,
    variant: str,
    full_gap_tolerance: pd.Timedelta = pd.Timedelta("3h"),
) -> None:
    """Fail loudly if the produced window has holes over the requested date range.

    Truncated raw-cache months (a fetch that stopped mid-month) and partial
    downloads are otherwise silent: the slice simply has fewer rows, and the hole
    only surfaces far downstream as NaN after a reindex (e.g. a solve aligning
    prices to a full CF year). This guard catches them at the source.

    dayahead is resampled to a clean hourly grid, so every hour of the window must
    be present and non-null. full spans a mixed resolution (ENTSO-E switched DE_LU
    day-ahead to 15-min in Oct 2025; NEM tables are 5-min), so a fixed hourly grid
    would false-positive; instead we flag any gap larger than `full_gap_tolerance`,
    which only occurs on truncated months.
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
                f"{len(missing)} missing hours: {summarise_runs(missing, pd.Timedelta('1h'))}"
            )
        n_nan = int(out_df.isna().any(axis=1).sum())
        if n_nan:
            problems.append(f"{n_nan} rows are NaN")
    else:
        gaps = idx.to_series().diff()
        if gaps.max() > full_gap_tolerance:
            problems.append(f"gap of {gaps.max()} ending {gaps.idxmax()}")

    if problems:
        raise ValueError(f"{variant} {start_date}–{end_date} incomplete: " + "; ".join(problems))
