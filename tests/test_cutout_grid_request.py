"""Unit tests for the cutout grid-vs-request check.

The keyed cutout cache hashes the bbox and resolution, so a cache hit is only
trustworthy if the file really is the cutout that request asked for. One route
that fills the cache — the `_backup.nc` fallback — copies whatever the backup
contained and then files it under the params of the request it stood in for, so
the check has to be enforced where the cache is written. Synthetic lattices, no
real cutouts.
"""

import numpy as np
import pytest
import xarray as xr

import common._cutout_cache as cutout_cache
from common._cutout_qc import CutoutQCError, grid_request_mismatches

# A cutout request and the file atlite actually produces for it: the domain is
# snapped to the ERA5 lattice inside the requested box, so the file's bounds sit
# a fraction of a grid step in from the request on every side.
NATIVE_REQUEST = {
    "module": "era5",
    "x0": 4.03097,
    "x1": 16.017,
    "y0": 46.30249,
    "y1": 56.51119,
    "dx": None,
    "dy": None,
    "start_date": "20230101",
    "end_date": "20231231",
}


def _lattice_file(tmp_path, x_bounds, y_bounds, step, name="cutout.nc"):
    """A minimal netCDF carrying only the x/y lattice the grid check reads."""
    x = np.round(np.arange(x_bounds[0], x_bounds[1] + step / 2, step), 5)
    y = np.round(np.arange(y_bounds[0], y_bounds[1] + step / 2, step), 5)
    path = tmp_path / name
    xr.Dataset(coords={"x": x, "y": y}).to_netcdf(path)
    return path


def test_correctly_snapped_cutout_has_no_mismatches(tmp_path):
    # Bounds land 0.78-0.88 grid steps inside the request, which is what a real
    # atlite download looks like — this must not be read as a mismatch.
    path = _lattice_file(tmp_path, (4.25, 16.0), (46.5, 56.5), 0.25)
    assert grid_request_mismatches(path, NATIVE_REQUEST) == []


def test_half_resolution_file_under_native_request_is_caught(tmp_path):
    # The VIC1 case: a 0.5 degree backup filed under a request that asked for
    # ERA5's native lattice.
    path = _lattice_file(tmp_path, (4.5, 16.0), (46.5, 56.5), 0.5)
    mismatches = grid_request_mismatches(path, NATIVE_REQUEST)
    assert any("x resolution 0.5 deg" in mismatch for mismatch in mismatches)
    assert any("y resolution 0.5 deg" in mismatch for mismatch in mismatches)


def test_bounds_short_of_the_request_are_caught(tmp_path):
    # Issue #62's case: a file built before the offshore zone was unioned into
    # the bbox, so it is ~4 grid steps short in the west and ~3 in the north.
    path = _lattice_file(tmp_path, (5.0, 16.0), (46.5, 55.75), 0.25)
    mismatches = grid_request_mismatches(path, NATIVE_REQUEST)
    assert any(mismatch.startswith("x0 bound") for mismatch in mismatches)
    assert any(mismatch.startswith("y1 bound") for mismatch in mismatches)
    assert not any("resolution" in mismatch for mismatch in mismatches)


def test_coarse_request_accepts_the_coarse_file_it_asked_for(tmp_path):
    # Same file as the resolution failure above, but now the request declares
    # dx/dy 0.5 — the pair is consistent and must pass.
    coarse_request = NATIVE_REQUEST | {"dx": 0.5, "dy": 0.5}
    path = _lattice_file(tmp_path, (4.5, 16.0), (46.5, 56.5), 0.5)
    assert grid_request_mismatches(path, coarse_request) == []


def test_store_in_cache_refuses_a_mismatched_file_and_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(cutout_cache, "CACHE_DIR", tmp_path / "cache")
    path = _lattice_file(tmp_path, (4.5, 16.0), (46.5, 56.5), 0.5)

    with pytest.raises(CutoutQCError, match="not the one the request describes"):
        cutout_cache.store_in_cache(path, "de", NATIVE_REQUEST)

    assert list((tmp_path / "cache").glob("*")) == [] or not (tmp_path / "cache").exists()


def test_store_in_cache_accepts_the_cutout_the_request_describes(tmp_path, monkeypatch):
    monkeypatch.setattr(cutout_cache, "CACHE_DIR", tmp_path / "cache")
    path = _lattice_file(tmp_path, (4.25, 16.0), (46.5, 56.5), 0.25)

    cached = cutout_cache.store_in_cache(path, "de", NATIVE_REQUEST)

    assert cached.exists()
    assert cached.with_suffix(".json").exists()
