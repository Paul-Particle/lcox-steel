"""Retrieve NEM grid data for a (area, variant, date-range) slice.

Maintains a persistent second-level processed cache at
  resources/nem/{variant}.parquet

Cache columns are a MultiIndex (area, metric) so all areas share one file per
variant. Accessing one area's data: df["VIC1"]. NEMOSIS manages its own raw
cache under data/nem_cache/; no per-month raw files are written here.

Variants
--------
dayahead  price table only → single "price" column, hourly UTC-naive, EUR/MWh
full      all four tables → wide per-area frame with derived wind/residual columns
"""

import logging
from pathlib import Path

import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake

from _helpers import area_month_in_cache, iso, to_utc_naive
from download_nem import DOWNLOADERS

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("retrieve_nem")

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
    start_str, end_str = _month_range_str(ym)
    raw = DOWNLOADERS["price"](start_str, end_str, cache_dir, rebuild=False)
    raw = to_utc_naive(raw)
    price = (raw[(area, "price")].resample("1h").mean() * eur_per_aud).rename("price")
    return price.to_frame()


def _process_full_month(
    area: str, ym: str, cache_dir: Path, eur_per_aud: float
) -> pd.DataFrame:
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

def retrieve(snakemake) -> None:
    area        = snakemake.wildcards.area
    variant     = snakemake.wildcards.variant
    start_date  = snakemake.wildcards.start_date
    end_date    = snakemake.wildcards.end_date
    eur_per_aud = snakemake.params.eur_per_aud

    cache_dir            = Path("data/nem_cache")
    processed_cache_dir  = Path("resources/nem")
    processed_cache_path = processed_cache_dir / f"{variant}.parquet"

    cached = pd.read_parquet(processed_cache_path) if processed_cache_path.exists() else None

    months = [
        ts.strftime("%Y-%m")
        for ts in pd.date_range(
            start=pd.Timestamp(iso(start_date)),
            end=pd.Timestamp(iso(end_date)),
            freq="MS",
        )
    ]

    new_frames = []
    for ym in months:
        if area_month_in_cache(cached, area, ym):
            continue
        log.info(f"{area}/{ym}/{variant}: processing")
        if variant == "dayahead":
            frame = _process_dayahead_month(area, ym, cache_dir, eur_per_aud)
        else:
            frame = _process_full_month(area, ym, cache_dir, eur_per_aud)
        frame.columns = pd.MultiIndex.from_tuples([(area, c) for c in frame.columns])
        new_frames.append(frame)

    if new_frames:
        all_frames = ([cached] if cached is not None else []) + new_frames
        cached = pd.concat(all_frames)
        cached = cached[~cached.index.duplicated(keep="last")].sort_index()
        processed_cache_dir.mkdir(parents=True, exist_ok=True)
        cached.to_parquet(processed_cache_path, index=True)
        log.info(f"Updated processed cache: {processed_cache_path} ({len(cached)} rows)")

    window = slice(iso(start_date), f"{iso(end_date)} 23:59")
    out_df = cached[area].loc[window]
    out_df.index.name = "time"

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"Wrote {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")


if __name__ == "__main__":
    retrieve(snakemake)
