"""Unit tests for the ENTSO-E output completeness guard.

`_assert_window_complete` is the pipeline's defence against silently truncated
raw-cache months and partial downloads — holes that otherwise only surface far
downstream as NaN after a reindex (e.g. a solve aligning prices to a full CF
year). These tests pin the guard's behaviour on synthetic frames, so they need
no raw cache and run anywhere.
"""

import pandas as pd
import pytest

import retrieve_entsoe  # sys.path set by conftest

_assert_window_complete = retrieve_entsoe._assert_window_complete


def _dayahead_frame(start="2024-01-01", end="2024-12-31"):
    """A complete hourly single-column dayahead-shaped frame for the window."""
    idx = pd.date_range(f"{start} 00:00", f"{end} 23:00", freq="h", name="time")
    return pd.DataFrame({"price": range(len(idx))}, index=idx, dtype=float)


def _full_frame(idx):
    """A full-variant-shaped frame (several columns, no NaN) on the given index."""
    return pd.DataFrame({"price": 1.0, "load": 2.0, "res": 3.0}, index=idx)


# ── dayahead: strict hourly grid ────────────────────────────────────────────────

def test_complete_dayahead_passes():
    out = _dayahead_frame()
    _assert_window_complete(out, "20240101", "20241231", "dayahead")


def test_missing_day_raises():
    out = _dayahead_frame()
    out = out.drop(out.loc["2024-06-15"].index)  # drop 24 interior hours
    with pytest.raises(ValueError, match="missing hours"):
        _assert_window_complete(out, "20240101", "20241231", "dayahead")


def test_month_boundary_truncation_raises():
    # Reproduces the observed bug: the last ~20h of each month absent.
    out = _dayahead_frame()
    for month in range(1, 13):
        last = out.loc[f"2024-{month:02d}"].index[-1].normalize()
        out = out.drop(out.loc[f"{last:%Y-%m-%d} 03:00":f"{last:%Y-%m-%d} 23:00"].index)
    with pytest.raises(ValueError, match="missing hours"):
        _assert_window_complete(out, "20240101", "20241231", "dayahead")


def test_trailing_truncation_raises():
    # A fetch that stopped part-way through the window (cf. Nov 2025).
    out = _dayahead_frame().loc[:"2024-11-05 23:00"]
    with pytest.raises(ValueError, match="before"):
        _assert_window_complete(out, "20240101", "20241231", "dayahead")


def test_nan_after_fill_raises():
    out = _dayahead_frame()
    out.loc["2024-03-10 05:00", "price"] = float("nan")
    with pytest.raises(ValueError, match="NaN"):
        _assert_window_complete(out, "20240101", "20241231", "dayahead")


# ── full: mixed-resolution tolerant, truncation caught ──────────────────────────

def test_full_mixed_resolution_passes():
    # Hourly through Sep, 15-min from Oct — the real Oct-2025 switch. No holes,
    # so the guard must not false-positive on the resolution change.
    hourly = pd.date_range("2025-01-01 00:00", "2025-09-30 23:00", freq="h")
    quarter = pd.date_range("2025-10-01 00:00", "2025-12-31 23:45", freq="15min")
    idx = hourly.union(quarter)
    _assert_window_complete(_full_frame(idx), "20250101", "20251231", "full")


def test_full_large_gap_raises():
    idx = pd.date_range("2025-01-01 00:00", "2025-12-31 23:00", freq="h")
    idx = idx.drop(idx[(idx >= "2025-11-05") & (idx < "2025-11-30")])  # ~missing month
    with pytest.raises(ValueError, match="gap of"):
        _assert_window_complete(_full_frame(idx), "20250101", "20251231", "full")
