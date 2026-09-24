"""Regression tests for merging processed months into the shared ONS cache.

One cache file per variant holds every area, keyed on MultiIndex (area, metric)
columns. Areas share timestamps but occupy different columns, so concatenating a
new area onto the cache produces one duplicate row per area with NaN elsewhere.
De-duplicating that by timestamp keeps only the most recently written area and
silently blanks the rest — the outputs still look right, because each area's
file is written before the next area clobbers it, so the damage is invisible
until something re-reads the cache. With four Brazilian submarkets it fires on
the second area. Synthetic frames, no download.
"""

import numpy as np
import pandas as pd

import _ons  # sys.path set by conftest

_merge_into_cache = _ons._merge_into_cache

HOURS = pd.date_range("2025-01-01", "2025-01-07 23:00", freq="h")


def _area_frame(area: str, value: float, index: pd.DatetimeIndex = HOURS) -> pd.DataFrame:
    frame = pd.DataFrame({"price": np.full(len(index), value)}, index=index)
    frame.columns = pd.MultiIndex.from_tuples([(area, c) for c in frame.columns])
    return frame


def test_second_area_does_not_blank_the_first():
    cached = _merge_into_cache(None, [_area_frame("SE", 1.0)])
    cached = _merge_into_cache(cached, [_area_frame("S", 2.0)])

    assert (cached[("SE", "price")] == 1.0).all()
    assert (cached[("S", "price")] == 2.0).all()
    assert len(cached) == len(HOURS)


def test_all_four_submarkets_coexist():
    cached = None
    for i, area in enumerate(("SE", "S", "NE", "N")):
        cached = _merge_into_cache(cached, [_area_frame(area, float(i))])

    assert len(cached) == len(HOURS)
    for i, area in enumerate(("SE", "S", "NE", "N")):
        assert (cached[(area, "price")] == float(i)).all(), area


def test_reprocessing_an_area_overwrites_its_own_values():
    """A refreshed month must win over the stale one it replaces."""
    cached = _merge_into_cache(None, [_area_frame("SE", 1.0)])
    cached = _merge_into_cache(cached, [_area_frame("SE", 9.0)])

    assert (cached[("SE", "price")] == 9.0).all()
    assert len(cached) == len(HOURS)


def test_areas_with_different_spans_are_unioned():
    later = pd.date_range("2025-02-01", "2025-02-07 23:00", freq="h")
    cached = _merge_into_cache(None, [_area_frame("SE", 1.0)])
    cached = _merge_into_cache(cached, [_area_frame("S", 2.0, index=later)])

    assert len(cached) == len(HOURS) + len(later)
    # Each area keeps its own span and is NaN outside it — which is what the
    # coverage bounds in retrieve() rely on to avoid "filling" another area's months.
    assert cached.loc[HOURS, ("SE", "price")].notna().all()
    assert cached.loc[HOURS, ("S", "price")].isna().all()
    assert cached.loc[later, ("S", "price")].notna().all()


def test_index_is_sorted_and_unique():
    cached = _merge_into_cache(None, [_area_frame("SE", 1.0)])
    cached = _merge_into_cache(cached, [_area_frame("S", 2.0)])

    assert cached.index.is_monotonic_increasing
    assert cached.index.is_unique
