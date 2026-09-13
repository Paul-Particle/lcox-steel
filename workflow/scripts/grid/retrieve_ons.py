"""Retrieve ONS (Brazil) grid data for a (area, variant, date-range) slice.

Maintains a persistent second-level processed cache at
  resources/ons/{variant}.parquet

Cache columns are a MultiIndex (area, metric) so all areas share one file per
variant. Accessing one area's data: df["SE"]. Areas are the four SIN submarkets
(SE, S, NE, N); the balance dataset's "SIN" national total is not a market area
and is skipped.

Variants
--------
dayahead  CMO only → single "price" column, hourly UTC-naive, EUR/MWh
full      CMO + subsystem energy balance → wide per-area frame with derived
          res/residual columns, hourly UTC-naive

Gap filling
-----------
ONS drops whole days from the CMO series (the hourly balance series is complete
for 2023-2025). Missing rows are filled from the same hour one week earlier,
else one week later, and the counts are logged. There is no ffill fallback — see
_fill_from_adjacent_week.

Price caveat
------------
ONS publishes the Marginal Operating Cost (CMO) — the DESSEM model's shadow
price. The settled market price is CCEE's PLD, which is the hourly CMO clipped
to an annual floor and cap set by ANEEL. This script reconstructs PLD from CMO
by applying those limits (config `ons.pld_limits`), because the clip is
material: raw 2025 SE CMO runs from -0.08 to 2151.70 R$/MWh, well outside the
2025 band of 58.60-1542.23. The reconstruction is *not* validated against
published PLD — CCEE's portal IP-blocks this host (see README). Raw CMO is
preserved in the `full` variant as `cmo_brl` so nothing is lost to the clip.
"""

import logging
from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from _helpers import (
    ONS_MARKET_TZ,
    area_month_in_cache,
    assert_window_complete,
    iso,
    summarise_runs,
    to_utc_naive,
)
from download_ons import SUBSYSTEMS, read_area_year

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Brazil abolished daylight saving from 2019; earlier years carry DST transitions
# that make the fixed-offset localisation in to_utc_naive ambiguous.
FIRST_DST_FREE_YEAR = 2019

BALANCE_COLUMNS = {
    "val_gerhidraulica": "hydro",
    "val_gertermica": "thermal",
    "val_gereolica": "wind",
    "val_gersolar": "solar",
    "val_carga": "load",
    "val_intercambio": "intercambio",
}


# ── Window mapping ────────────────────────────────────────────────────────────

def _local_months(start_date: str, end_date: str) -> list[str]:
    """Brasilia-time months ('YYYY-MM') whose rows the requested UTC window needs.

    ONS stamps din_instante in Brasilia time (UTC-3) and the rule output is
    UTC-naive, so the window's first UTC hours come from the *previous* local
    month — the mirror image of the forward pad retrieve_nem needs for AEST
    (UTC+10). Converting both window bounds into local time and taking the months
    they span handles the year boundary too, where the previous local month lives
    in the previous year's raw file.
    """
    to_local = lambda ts: ts.tz_localize("UTC").tz_convert(ONS_MARKET_TZ).tz_localize(None)
    local_start = to_local(pd.Timestamp(iso(start_date)))
    local_end = to_local(pd.Timestamp(f"{iso(end_date)} 23:00"))
    return pd.period_range(local_start, local_end, freq="M").strftime("%Y-%m").tolist()


# ── Per-variant month processing ──────────────────────────────────────────────

def _hourly_cmo(area: str, ym: str, cache_dir: Path) -> pd.Series:
    """Return one local month of CMO as an hourly UTC-naive series, R$/MWh.

    DESSEM publishes CMO semi-hourly; PLD is formed on the hourly mean, so the
    resample happens before any clipping.
    """
    year = int(ym[:4])
    raw = read_area_year("cmo", area, year, cache_dir)
    raw = raw.loc[ym]
    return to_utc_naive(raw, naive_tz=ONS_MARKET_TZ)["val_cmo"].resample("1h").mean()


def _to_pld_proxy(cmo_brl: pd.Series, limits: dict, eur_per_brl: float) -> pd.Series:
    """Clip CMO to the year's ANEEL PLD band and convert R$/MWh → EUR/MWh."""
    floor, cap = limits["min"], limits["max_hourly"]
    n_floored = int((cmo_brl < floor).sum())
    n_capped = int((cmo_brl > cap).sum())
    if n_floored or n_capped:
        log.info(
            f"PLD clip to [{floor}, {cap}] R$/MWh: "
            f"{n_floored} hours at floor, {n_capped} at cap ({len(cmo_brl)} total)"
        )
    return cmo_brl.clip(floor, cap) * eur_per_brl


def _process_dayahead_month(
    area: str, ym: str, cache_dir: Path, eur_per_brl: float, limits: dict
) -> pd.DataFrame:
    """Return one local month of prices as a single hourly, UTC-naive `price` column."""
    cmo_brl = _hourly_cmo(area, ym, cache_dir)
    return _to_pld_proxy(cmo_brl, limits, eur_per_brl).rename("price").to_frame()


def _process_full_month(
    area: str, ym: str, cache_dir: Path, eur_per_brl: float, limits: dict
) -> pd.DataFrame:
    """Join one local month's CMO and energy balance into a wide per-area frame.

    The balance dataset is natively hourly and CMO is resampled to match, so the
    joined frame is hourly throughout — unlike the ENTSO-E and NEM `full`
    variants, which keep their sources' native sub-hourly resolution.
    """
    year = int(ym[:4])
    balance = read_area_year("balance", area, year, cache_dir).loc[ym]
    balance = to_utc_naive(balance, naive_tz=ONS_MARKET_TZ).rename(columns=BALANCE_COLUMNS)

    cmo_brl = _hourly_cmo(area, ym, cache_dir)
    balance["cmo_brl"] = cmo_brl
    balance["price"] = _to_pld_proxy(cmo_brl, limits, eur_per_brl)

    balance["res"] = balance["wind"] + balance["solar"]
    balance["residual"] = balance["load"] - balance["res"]

    return balance.sort_index(axis=1)


# ── Processed cache ───────────────────────────────────────────────────────────

def _merge_into_cache(cached: pd.DataFrame | None, new_frames: list) -> pd.DataFrame:
    """Merge freshly processed months into the shared per-variant cache.

    Merges on both axes rather than concatenating and de-duplicating. Areas share
    timestamps but occupy different columns, so concatenating them produces one
    duplicate row per area (NaN elsewhere), and dropping duplicates keeps only the
    area written most recently — silently blanking the others. With four Brazilian
    submarkets that fires on the second area. combine_first unions index and
    columns instead, letting the fresh block win wherever the two overlap.
    """
    new_block = pd.concat(new_frames)
    new_block = new_block[~new_block.index.duplicated(keep="last")]
    merged = new_block if cached is None else new_block.combine_first(cached)
    return merged.sort_index()


# ── Gap filling ───────────────────────────────────────────────────────────────

def _fill_from_adjacent_week(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Fill NaNs from the same hour one week earlier, else one week later.

    ONS drops whole days from the CMO series (the hourly balance series is
    unaffected), so the hole to cover is 24 h wide. A same-hour-of-week donor
    keeps the diurnal and weekday shape that a plain ffill across a whole day
    destroys. Returns the filled frame and an audit dict of what each step
    touched — the caller logs it.

    Deliberately no ffill fallback: if both weekly donors are missing, the outage
    spans three aligned weeks and no cheap fill is defensible. Leaving those NaN
    lets assert_window_complete fail loudly instead. (PyPSA-Brazil does fall back
    to ffill, and their own audit shows 4752 of 5895 fills went that way — almost
    all of it ffill across multi-hundred-hour outages, which is the case this
    refuses to paper over.)
    """
    audit = {"missing": int(frame.isna().any(axis=1).sum())}
    if not audit["missing"]:
        return frame, audit

    filled = frame
    for label, offset in (("prev_week", "7D"), ("next_week", "-7D")):
        before = filled.isna().any(axis=1)
        # shift(freq=+7D) relabels t -> t+7d, so at t it carries the value from t-7d.
        donor = filled.shift(freq=offset).reindex(filled.index)
        filled = filled.fillna(donor)
        audit[label] = int((before & ~filled.isna().any(axis=1)).sum())

    audit["unfilled"] = int(filled.isna().any(axis=1).sum())
    return filled, audit


# ── Main ──────────────────────────────────────────────────────────────────────

def retrieve(snakemake) -> None:
    """Slice the requested (area, variant, date range) out of the processed cache.

    Downloads any absent dataset-years into the raw cache, processes the local
    months missing from the shared per-variant processed cache, appends them,
    then writes the requested UTC window to the rule output.
    """
    area = snakemake.wildcards.area
    variant = snakemake.wildcards.variant
    start_date = snakemake.wildcards.start_date
    end_date = snakemake.wildcards.end_date
    eur_per_brl = snakemake.params.eur_per_brl
    pld_limits = snakemake.params.pld_limits

    cache_dir = Path("data/ons_cache")
    processed_cache_dir = Path("resources/ons")
    processed_cache_path = processed_cache_dir / f"{variant}.parquet"

    if variant not in ("dayahead", "full"):
        raise ValueError(f"Unknown variant {variant!r}. Expected 'dayahead' or 'full'.")
    if area not in SUBSYSTEMS:
        raise ValueError(f"{area!r} is not an ONS subsystem. Expected one of {SUBSYSTEMS}.")

    months = _local_months(start_date, end_date)

    first_year = int(months[0][:4])
    if first_year < FIRST_DST_FREE_YEAR:
        raise ValueError(
            f"Window reaches back to {first_year}. Brazil observed daylight saving "
            f"until 2019, so Brasilia-time timestamps before then are ambiguous and "
            f"this pipeline refuses them."
        )
    missing_limits = sorted({int(ym[:4]) for ym in months} - pld_limits.keys())
    if missing_limits:
        raise ValueError(
            f"No ANEEL PLD limits configured for {missing_limits}. "
            f"Add them under `ons.pld_limits` in config/config.yaml."
        )

    cached = pd.read_parquet(processed_cache_path) if processed_cache_path.exists() else None

    new_frames = []
    for ym in months:
        # Match cache membership in market time: ONS publishes whole local months
        # but the cache is stored in UTC, so a local month spills across two UTC
        # months. A plain UTC-month check would see a neighbour's spillover and
        # wrongly skip this month. Same reasoning as retrieve_nem.
        if area_month_in_cache(cached, area, ym, tz=ONS_MARKET_TZ):
            continue
        log.info(f"{area}/{ym}/{variant}: processing")
        limits = pld_limits[int(ym[:4])]
        if variant == "dayahead":
            frame = _process_dayahead_month(area, ym, cache_dir, eur_per_brl, limits)
        else:
            frame = _process_full_month(area, ym, cache_dir, eur_per_brl, limits)
        frame.columns = pd.MultiIndex.from_tuples([(area, c) for c in frame.columns])
        new_frames.append(frame)

    if new_frames:
        cached = _merge_into_cache(cached, new_frames)
        processed_cache_dir.mkdir(parents=True, exist_ok=True)
        cached.to_parquet(processed_cache_path, index=True)
        log.info(f"updated processed cache: {processed_cache_path} ({len(cached)} rows)")

    window = slice(iso(start_date), f"{iso(end_date)} 23:00")

    # Fill before slicing, on a complete hourly index spanning the area's own
    # coverage: a gap in the window's first or last week needs a donor that sits
    # outside the window, and reindexing turns a dropped day into fillable NaN
    # rows. Bounding by this area's observed span keeps other areas' months —
    # which share the cache index — from looking like gaps in this one.
    observed = cached[area].dropna(how="all").index
    area_df = cached[area].reindex(pd.date_range(observed.min(), observed.max(), freq="h"))

    area_df, audit = _fill_from_adjacent_week(area_df)
    if audit["missing"]:
        log.info(
            f"gap fill over {area} coverage: {audit['missing']} rows missing → "
            f"{audit['prev_week']} from the previous week, "
            f"{audit['next_week']} from the next, {audit['unfilled']} left NaN"
        )

    out_df = area_df.loc[window]
    out_df.index.name = "time"

    # assert_window_complete only checks `full` for index gaps, so a `full` slice
    # could otherwise carry NaN on a complete index. Anything still NaN here
    # survived both weekly donors.
    nan_times = out_df.index[out_df.isna().any(axis=1)]
    if variant == "full" and len(nan_times):
        log.warning(
            f"{len(nan_times)} rows still NaN after gap fill: "
            f"{summarise_runs(nan_times, pd.Timedelta('1h'))}"
        )

    assert_window_complete(out_df, start_date, end_date, variant)

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")


if __name__ == "__main__":
    retrieve(snakemake)
