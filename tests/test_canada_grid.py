"""Unit tests for the Canadian grid sources' market-time to UTC conversion.

The conversion is where Alberta and Ontario differ, and where a mistake would be
silent. AESO labels hour endings in Mountain wall-clock time, so the label skips
an hour each spring (01 → 03) and repeats one each autumn (02, 02*); arithmetic on
the label lands on hours that do not exist. IESO holds EST all year, so a fixed
offset is the whole conversion and every market day is 24 hours.

`read_*_year` downloads only when its cache file is absent, so writing a fixture
at the path it would use exercises the parsing with no network.
"""

import pandas as pd
import pytest

import _canada  # sys.path set by conftest

AESO_PREAMBLE = 'Pool Price\n\n""\n\nDate (HE),Price ($),30Ravg ($),AIL Demand (MW)\n'
IESO_PREAMBLE = (
    "\\\\Yearly HOEP OR Predispatch Report,,\n\\\\Created at 2024-01-31 08:02:24,,\n"
    "\\\\For 2023,,\nDate,Hour,HOEP\n"
)


def _write_aeso_fixture(cache_dir, year, rows):
    """Place an AESO-shaped CSV where read_aeso_year would cache its download."""
    body = "".join(f'"{date} {he}","{price}","0.0","9000.0"\n' for date, he, price in rows)
    path = cache_dir / f"aeso_pool_price_{year}.csv"
    path.write_text(AESO_PREAMBLE + body)
    return path


def _write_ieso_fixture(cache_dir, year, rows):
    """Place an IESO-shaped CSV where read_ieso_year would cache its download."""
    body = "".join(f"{date},{hour},{price}\n" for date, hour, price in rows)
    path = cache_dir / f"ieso_hoep_{year}.csv"
    path.write_text(IESO_PREAMBLE + body)
    return path


# ── AESO: Mountain time, with daylight saving ─────────────────────────────────

def test_aeso_spring_forward_day_is_23_contiguous_hours(tmp_path):
    # 12 March 2023: clocks jump 02:00 MST → 03:00 MDT, and AESO omits HE 02.
    hours = [f"{h:02d}" for h in range(1, 25) if h != 2]
    _write_aeso_fixture(tmp_path, 2023, [("03/12/2023", he, "10.0") for he in hours])

    frame = _canada.read_aeso_year(2023, tmp_path)

    assert len(frame) == 23
    # Local midnight is 00:00 MST = 07:00 UTC; the 23 hours then run unbroken.
    assert frame.index[0] == pd.Timestamp("2023-03-12 07:00")
    assert frame.index[-1] == pd.Timestamp("2023-03-13 05:00")
    assert (frame.index.to_series().diff().dropna() == pd.Timedelta("1h")).all()


def test_aeso_fall_back_day_is_25_contiguous_hours(tmp_path):
    # 5 November 2023: clocks go back 02:00 MDT → 01:00 MST, so the hour ending
    # 02:00 happens twice and AESO labels the second one 02*.
    hours = ["01", "02", "02*"] + [f"{h:02d}" for h in range(3, 25)]
    _write_aeso_fixture(tmp_path, 2023, [("11/05/2023", he, "10.0") for he in hours])

    frame = _canada.read_aeso_year(2023, tmp_path)

    assert len(frame) == 25
    assert frame.index.is_unique
    # Local midnight is 00:00 MDT = 06:00 UTC; 25 unbroken hours follow it.
    assert frame.index[0] == pd.Timestamp("2023-11-05 06:00")
    assert frame.index[-1] == pd.Timestamp("2023-11-06 06:00")
    assert (frame.index.to_series().diff().dropna() == pd.Timedelta("1h")).all()


def test_aeso_starred_hour_keeps_its_own_price(tmp_path):
    # The repeated hour is a distinct settlement hour, not a duplicate to drop.
    rows = [("11/05/2023", "01", "11.0"), ("11/05/2023", "02", "22.0"),
            ("11/05/2023", "02*", "33.0"), ("11/05/2023", "03", "44.0")]
    _write_aeso_fixture(tmp_path, 2023, rows)

    prices = _canada.read_aeso_year(2023, tmp_path)["price"].tolist()

    assert prices == [11.0, 22.0, 33.0, 44.0]


# ── IESO: EST all year ────────────────────────────────────────────────────────

def test_ieso_hour_ending_one_is_local_midnight(tmp_path):
    _write_ieso_fixture(tmp_path, 2023, [("2023-01-01", 1, 14.42), ("2023-01-01", 2, 19.21)])

    frame = _canada.read_ieso_year(2023, tmp_path)

    # HE 1 covers 00:00–01:00 EST, and EST is UTC-5 the whole year.
    assert frame.index[0] == pd.Timestamp("2023-01-01 05:00")
    assert frame.index[1] == pd.Timestamp("2023-01-01 06:00")
    assert frame["price"].tolist() == [14.42, 19.21]


@pytest.mark.parametrize("date", ["2023-03-12", "2023-11-05"])
def test_ieso_keeps_24_hours_across_dst_dates(tmp_path, date):
    # IESO does not observe daylight saving, so the transition dates are ordinary
    # 24-hour days and the offset never moves.
    _write_ieso_fixture(tmp_path, 2023, [(date, h, 10.0) for h in range(1, 25)])

    frame = _canada.read_ieso_year(2023, tmp_path)

    assert len(frame) == 24
    assert frame.index[0] == pd.Timestamp(f"{date} 05:00")
    assert (frame.index.to_series().diff().dropna() == pd.Timedelta("1h")).all()


# ── Variant support ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("variant", ["emissions", "full"])
def test_unsupported_variant_names_what_is_missing(variant, tmp_path):
    class _Stub:
        wildcards = type("W", (), {"variant": variant, "start_date": "20230101",
                                   "end_date": "20231231"})()
        params = type("P", (), {"eur_per_cad": 0.64})()
        output = [str(tmp_path / "out.parquet")]

    with pytest.raises(ValueError, match="not available for the Canadian markets"):
        _canada.retrieve_aeso(_Stub(), "AB")
