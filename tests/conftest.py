"""Pytest fixtures for the grid-pipeline validation tests.

Puts the workflow script dirs on ``sys.path`` so the ENTSO-E processing code can
be imported, and exposes the two datasets the validation compares:

* ``pipeline_full_de_2023`` — the ``full``-variant frame for DE_LU / 2023, rebuilt
  from the local ENTSO-E raw monthly cache via the pipeline's own
  ``_process_full_month``. Skipped if the raw cache isn't present locally.
* ``energycharts_ref`` — the committed hourly energy-charts.info reference
  (``tests/data/energycharts_de_2023_hourly.csv``).
"""

from pathlib import Path
import sys

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# _entsoe imports ``common.*`` (workflow/) and its siblings
# (_helpers_grid, download_entsoe) as top-level modules
# (workflow/scripts/grid/); the viz scripts likewise import each other by bare
# name. Each scripts/ dir suffixes its helper module with its own name, so a
# second dir on this flat path cannot shadow the first's.
for _p in (
    REPO_ROOT / "workflow",
    REPO_ROOT / "workflow" / "scripts" / "grid",
    REPO_ROOT / "workflow" / "scripts" / "viz",
):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

RAW_CACHE = REPO_ROOT / "data" / "entsoe_cache"
AREA = "DE_LU"
MONTHS = [f"2023-{month:02d}" for month in range(1, 13)]


@pytest.fixture(scope="session")
def pipeline_full_de_2023() -> pd.DataFrame:
    """DE_LU 2023 ``full`` frame rebuilt from the raw cache through the pipeline."""
    probe = RAW_CACHE / AREA / MONTHS[0] / "generation.parquet"
    if not probe.exists():
        pytest.skip(
            f"ENTSO-E raw cache not found ({probe}). Build it by running the grid "
            "pipeline for DE_LU 2023 (all six data types)."
        )

    import _entsoe  # imported here so sys.path (set above) is in effect

    frames = [_entsoe._process_full_month(AREA, ym, RAW_CACHE) for ym in MONTHS]
    combined = (
        pd.concat(frames)
        .pipe(lambda df: df[~df.index.duplicated(keep="last")])
        .sort_index()
    )
    return combined


@pytest.fixture(scope="session")
def energycharts_ref() -> pd.DataFrame:
    """Committed hourly energy-charts.info reference for DE, 2023."""
    ref_path = REPO_ROOT / "tests" / "data" / "energycharts_de_2023_hourly.csv"
    return pd.read_csv(ref_path, index_col="time", parse_dates=["time"])
