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
import sys
from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from common._stubs import snakemake

from _helpers import area_month_in_cache, iso, to_utc_naive
from download_entsoe import FETCHERS, fetch_with_retry, get_entsoe_client, iter_months

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("retrieve_entsoe")

FULL_DATA_TYPES = ["prices", "load_forecast", "load_actual", "res", "generation", "crossborder"]


# ── Raw-cache management ──────────────────────────────────────────────────────

def ensure_raw_months(
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
            df = fetch_with_retry(FETCHERS[dt], client, area, month_start, next_month)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(cache_path, index=True)
            log.info(f"{area}/{ym}/{dt}: cached")


# ── Per-variant month processing ──────────────────────────────────────────────

def process_dayahead_month(area: str, ym: str, raw_cache_dir: Path) -> pd.DataFrame:
    raw = pd.read_parquet(raw_cache_dir / area / ym / "prices.parquet")
    price = to_utc_naive(raw).iloc[:, 0]
    price = price.ffill(limit=3).bfill(limit=3).resample("1h").mean().rename("price")
    return price.to_frame()


def process_full_month(area: str, ym: str, raw_cache_dir: Path) -> pd.DataFrame:
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
    df = df.ffill(limit=3).fillna(0.0)

    df["wind_forecast"]     = df.get("wind_onshore_forecast", 0) + df.get("wind_offshore_forecast", 0)
    df["res_forecast"]      = df["wind_forecast"] + df.get("solar_forecast", 0)
    df["residual_forecast"] = df["load_forecast"] - df["res_forecast"]
    df["wind"]              = df.get("wind_onshore", 0) + df.get("wind_offshore", 0)
    df["res"]               = df["wind"] + df.get("solar", 0)
    df["residual"]          = df["load"] - df["res"]

    return df.sort_index(axis=1)


# ── Main ──────────────────────────────────────────────────────────────────────

def retrieve(snakemake) -> None:
    area       = snakemake.wildcards.area
    variant    = snakemake.wildcards.variant
    start_date = snakemake.wildcards.start_date
    end_date   = snakemake.wildcards.end_date

    raw_cache_dir        = Path("data/entsoe_cache")
    processed_cache_dir  = Path("resources/entsoe")
    processed_cache_path = processed_cache_dir / f"{variant}.parquet"

    data_types    = ["prices"] if variant == "dayahead" else FULL_DATA_TYPES
    process_month = process_dayahead_month if variant == "dayahead" else process_full_month

    cached = pd.read_parquet(processed_cache_path) if processed_cache_path.exists() else None

    months = [ym for ym, _, _ in iter_months(start_date, end_date)]

    ensure_raw_months(area, months, data_types, raw_cache_dir)

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
        log.info(f"Updated processed cache: {processed_cache_path} ({len(cached)} rows)")

    window = slice(iso(start_date), f"{iso(end_date)} 23:00")
    out_df = cached[area].loc[window]
    # Cross-month boundary gaps: forward-fill only (bfill would propagate future data backward).
    # Per-month processing already handles within-month gaps including start-of-month.
    out_df = out_df.ffill(limit=3)
    if variant == "full":
        out_df = out_df.fillna(0.0)
    out_df.index.name = "time"

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"Wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")


if __name__ == "__main__":
    retrieve(snakemake)
