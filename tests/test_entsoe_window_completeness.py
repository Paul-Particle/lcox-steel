"""ENTSO-E window completeness: the month pad that keeps summer-ending windows whole.

ENTSO-E raw months are fetched on Brussels-time boundaries and stored UTC-naive. In
summer (CEST = UTC+2) a Brussels month's data ends at 22:00 UTC on its last day, so a
window's final UTC hour lives in the *next* Brussels month. _entsoe pads the
download list by one month to cover it; without the pad `assert_window_complete` fails
every Mar–Sep-ending window (verified on real DE_LU 2023 data).

The unit check runs anywhere. The end-to-end check needs the local ENTSO-E raw price
cache and skips without it.
"""
from pathlib import Path

import pandas as pd
import pytest

import _entsoe as R  # sys.path set by conftest
from _helpers_grid import assert_window_complete, iso

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_CACHE = REPO_ROOT / "data" / "entsoe_cache"
AREA = "DE_LU"


# ── unit: the pad month is appended ─────────────────────────────────────────────

def test_pad_month_appended_for_summer_end():
    assert R._months_to_process("20230101", "20230630")[-1] == "2023-07"


def test_pad_month_rolls_into_next_year():
    assert R._months_to_process("20230101", "20231231")[-1] == "2024-01"


def test_pad_month_not_duplicated():
    # A window already spanning into the pad month must not repeat it.
    months = R._months_to_process("20230101", "20230701")
    assert months.count("2023-07") == 1


# ── end-to-end: summer-ending windows reach the final UTC hour (needs raw cache) ─

def _dayahead_window(start_date: str, end_date: str) -> pd.DataFrame:
    """Rebuild the dayahead window exactly as retrieve() does, with padded months."""
    months = R._months_to_process(start_date, end_date)
    frames = [R._process_dayahead_month(AREA, ym, RAW_CACHE) for ym in months]
    cached = pd.concat(frames)
    cached = cached[~cached.index.duplicated(keep="last")].sort_index()
    return cached.loc[slice(iso(start_date), f"{iso(end_date)} 23:00")].ffill(limit=3)


# Summer-ending windows (CEST) are the ones that raised pre-pad: without the pad
# month their UTC data stops at 22:00 on the last day. Jun→pad Jul, Sep→pad Oct,
# all within a DE_LU 2023 cache, so the check is self-contained.
@pytest.mark.parametrize("end_date", ["20230630", "20230930"])
def test_summer_window_reaches_final_hour(end_date):
    months = R._months_to_process("20230101", end_date)
    missing = [ym for ym in months if not (RAW_CACHE / AREA / ym / "prices.parquet").exists()]
    if missing:
        pytest.skip(f"ENTSO-E raw price cache missing months: {missing}")
    out = _dayahead_window("20230101", end_date)
    assert str(out.index.max()) == f"{iso(end_date)} 23:00:00"  # would be 22:00 pre-pad
    assert_window_complete(out, "20230101", end_date, "dayahead")
