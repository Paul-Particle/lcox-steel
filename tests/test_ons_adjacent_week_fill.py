"""Unit tests for the ONS adjacent-week gap fill and its audit counts.

ONS drops whole days from the CMO series, so the hole to cover is 24 h wide and
a plain ffill would flatten a full diurnal cycle. `_fill_from_adjacent_week`
takes the same hour one week earlier, falling back to one week later, and
returns an audit dict the rule logs. The audit is the only record of which
values are real and which are borrowed, so it is tested as carefully as the
fill itself. Synthetic frames, no download.
"""

import numpy as np
import pandas as pd

import retrieve_ons  # sys.path set by conftest

_fill_from_adjacent_week = retrieve_ons._fill_from_adjacent_week

HOURS = pd.date_range("2025-05-01", "2025-06-15 23:00", freq="h")


def _price_frame() -> pd.DataFrame:
    """A frame whose value encodes its own timestamp, so donors are identifiable."""
    return pd.DataFrame({"price": np.arange(len(HOURS), dtype=float)}, index=HOURS)


def _blank(frame: pd.DataFrame, day: str) -> pd.DataFrame:
    frame = frame.copy()
    frame.loc[day] = np.nan
    return frame


def test_hole_is_filled_from_the_previous_week():
    frame = _blank(_price_frame(), "2025-05-16")
    filled, audit = _fill_from_adjacent_week(frame)

    assert filled.isna().sum().sum() == 0
    # Every filled hour must equal the same hour seven days earlier.
    gap = pd.date_range("2025-05-16", "2025-05-16 23:00", freq="h")
    expected = _price_frame().loc[gap - pd.Timedelta("7D"), "price"].to_numpy()
    assert (filled.loc[gap, "price"].to_numpy() == expected).all()
    assert audit == {"missing": 24, "prev_week": 24, "next_week": 0, "unfilled": 0}


def test_next_week_is_used_when_the_previous_one_is_also_missing():
    frame = _blank(_blank(_price_frame(), "2025-05-16"), "2025-05-09")
    filled, audit = _fill_from_adjacent_week(frame)

    gap = pd.date_range("2025-05-16", "2025-05-16 23:00", freq="h")
    expected = _price_frame().loc[gap + pd.Timedelta("7D"), "price"].to_numpy()
    assert (filled.loc[gap, "price"].to_numpy() == expected).all()
    # Each pass reads donors from the frame as it stood *before* that pass, so
    # 05-16 cannot borrow the 05-09 values the same pass is filling. 05-09 takes
    # 05-02 on the previous-week pass; 05-16 takes 05-23 on the next-week pass.
    assert audit == {"missing": 48, "prev_week": 24, "next_week": 24, "unfilled": 0}


def test_three_aligned_weeks_missing_is_left_nan():
    frame = _price_frame()
    for day in ("2025-05-09", "2025-05-16", "2025-05-23"):
        frame = _blank(frame, day)
    filled, audit = _fill_from_adjacent_week(frame)

    middle = pd.date_range("2025-05-16", "2025-05-16 23:00", freq="h")
    assert filled.loc[middle, "price"].isna().all()
    assert audit["unfilled"] == 24


def test_audit_counts_account_for_every_missing_row():
    frame = _price_frame()
    for day in ("2025-05-09", "2025-05-16", "2025-05-23", "2025-06-02"):
        frame = _blank(frame, day)
    _, audit = _fill_from_adjacent_week(frame)

    assert audit["missing"] == audit["prev_week"] + audit["next_week"] + audit["unfilled"]


def test_rows_outside_the_hole_are_untouched():
    original = _price_frame()
    frame = _blank(original, "2025-05-16")
    filled, _ = _fill_from_adjacent_week(frame)

    untouched = original.index.difference(pd.date_range("2025-05-16", "2025-05-16 23:00", freq="h"))
    assert filled.loc[untouched].equals(original.loc[untouched])


def test_complete_frame_is_a_no_op():
    original = _price_frame()
    filled, audit = _fill_from_adjacent_week(original)

    assert filled.equals(original)
    assert audit == {"missing": 0}


def test_full_variant_fills_every_column_independently():
    """A `full` row can have NaN price but valid load — only the NaN must move."""
    frame = _price_frame()
    frame["load"] = 1000.0
    frame.loc["2025-05-16", "price"] = np.nan
    filled, audit = _fill_from_adjacent_week(frame)

    assert (filled["load"] == 1000.0).all()
    assert filled["price"].isna().sum() == 0
    assert audit["missing"] == 24
