"""Unit tests for the `emissions` grid variant's month processor.

The variant exists so a grid run carries the mix behind its imports, not just
the price it paid. Two things about the frame it produces are load-bearing
downstream and neither is obvious from reading it: `price` has to lead the
columns, and the index has to be a clean hourly UTC-naive grid — the solve
reindexes the price onto the CF year and `compile_report` reindexes the mix onto
the snapshots, and a 15-minute frame would silently hand both the :00 reading
instead of the hour's mean.

Synthetic raw caches, so no ENTSO-E credentials and no NEMOSIS download.
"""

import pandas as pd
import pytest

import _entsoe  # sys.path set by conftest
from _helpers_grid import assert_window_complete
from common._report_schema import field_stem

AREA = "DE_LU"
# A quarter-hourly month, which is what ENTSO-E serves DE_LU since Oct 2025.
INDEX = pd.date_range("2025-01-01", "2025-01-02", freq="15min", tz="UTC")[:-1]


@pytest.fixture
def raw_cache(tmp_path):
    """A month of prices and per-carrier generation in the raw-cache layout."""
    month_dir = tmp_path / AREA / "2025-01"
    month_dir.mkdir(parents=True)
    pd.DataFrame({"price": range(len(INDEX))}, index=INDEX).to_parquet(
        month_dir / "prices.parquet"
    )
    generation = pd.DataFrame(
        {"hard_coal": 100.0, "wind_onshore": 50.0, "solar": 0.0}, index=INDEX
    )
    generation.columns = pd.MultiIndex.from_tuples(
        [(AREA, carrier) for carrier in generation.columns]
    )
    generation.to_parquet(month_dir / "generation.parquet")
    return tmp_path


def test_price_leads_the_columns(raw_cache):
    """The solve reads the price by name, but a leading `price` keeps the frame
    readable and keeps any positional reader honest."""
    out = _entsoe._process_emissions_month(AREA, "2025-01", raw_cache)
    assert out.columns[0] == "price"


def test_the_carriers_come_through_named_as_the_factor_table_keys(raw_cache):
    """The whole point: these column names are what `emissions` looks up."""
    out = _entsoe._process_emissions_month(AREA, "2025-01", raw_cache)
    assert set(out.columns) == {"price", "hard_coal", "wind_onshore", "solar"}
    # And they are already the report's spelling — no translation step anywhere.
    assert all(field_stem(col) == col for col in out.columns)


def test_a_sub_hourly_month_is_resampled_to_hourly_means(raw_cache):
    """Both consumers are hourly, so the resolution is settled here rather than
    reindexed away twice, differently."""
    out = _entsoe._process_emissions_month(AREA, "2025-01", raw_cache)
    assert len(out) == 24
    assert (out.index.to_series().diff().dropna() == pd.Timedelta("1h")).all()
    assert out.index.tz is None
    # Means, not samples: the first hour's four quarter-hour prices are 0..3.
    assert out["price"].iloc[0] == pytest.approx(1.5)
    assert out["hard_coal"].iloc[0] == pytest.approx(100.0)


def test_the_variant_registry_names_only_what_it_fetches():
    """`emissions` is two data types, not `full`'s six — that is its reason to
    exist beside `full`, which would answer the same question at three times
    the download."""
    data_types, processor = _entsoe.VARIANTS["emissions"]
    assert data_types == ["prices", "generation"]
    assert processor is _entsoe._process_emissions_month
    assert set(data_types) < set(_entsoe.FULL_DATA_TYPES)


def test_the_variant_is_held_to_the_strict_hourly_check():
    """It is resampled to a clean hourly grid like `dayahead`, so it wants that
    guard rather than `full`'s gap tolerance — under which a month truncated by
    two hours passes."""
    hours = pd.date_range("2025-01-01", "2025-01-01 23:00", freq="h")
    complete = pd.DataFrame({"price": 50.0, "hard_coal": 100.0}, index=hours)
    assert_window_complete(complete, "20250101", "20250101", "emissions")

    with pytest.raises(ValueError, match="missing hours"):
        assert_window_complete(complete.drop(hours[5]), "20250101", "20250101", "emissions")


def test_a_hole_inside_a_carrier_column_does_not_pass_the_guard():
    """A carrier the zone never reported all window is filled with zeros before
    the guard sees it — it has no plants of that kind. A hole *inside* a column
    is a truncated fetch, and an hour where nothing was burned is not the same
    statement."""
    hours = pd.date_range("2025-01-01", "2025-01-01 23:00", freq="h")
    truncated = pd.DataFrame({"price": 50.0, "hard_coal": 100.0}, index=hours)
    truncated.loc[hours[7:12], "hard_coal"] = pd.NA

    with pytest.raises(ValueError, match="rows are NaN"):
        assert_window_complete(truncated, "20250101", "20250101", "emissions")
