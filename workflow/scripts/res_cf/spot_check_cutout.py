"""Cheap CDS spot-check: probe (country, year) download health at ~1 grid cell.

Rationale
---------
The download-stage failures we care about (partial years, ERA5/ERA5T mixing,
missing months) are properties of *time and data availability*, not of spatial
extent — ERA5 publishes each hour globally or not at all. So a single ERA5 grid
cell over a country, requested for the *full year*, exercises the identical
atlite download/GRIB/merge path across all 12 months while moving a few MB
instead of hundreds. That makes it a fast, cheap reproduction of the real bug
before committing to full country-year cutouts.

Runs standalone (no Snakemake):

    python workflow/scripts/res_cf/spot_check_cutout.py --areas de fr es aus bra --years 2024 2025

Exits non-zero if any probe fails QC.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Standalone: put workflow/ on sys.path so `common.*` imports resolve.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import atlite

from common._cutout_qc import check_cutout
from common._paths import ATLITE_CACHE, CUTOUTS

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("spot_check_cutout")

# Representative interior point per area (lon, lat). Location within a country is
# irrelevant to the availability failures we probe, so a hardcoded point avoids
# any dependency on the built geo parquets.
PROBE_POINTS = {
    "de": (10.0, 51.0),
    "fr": (2.5, 46.5),
    "es": (-3.7, 40.4),
    "aus": (145.0, -37.0),
    "bra": (-47.9, -15.8),
}

SPOT_DIR = CUTOUTS / "spotcheck"


def probe(area: str, year: int) -> bool:
    lon, lat = PROBE_POINTS[area]
    start, end = f"{year}0101", f"{year}1231"
    SPOT_DIR.mkdir(parents=True, exist_ok=True)
    ATLITE_CACHE.mkdir(parents=True, exist_ok=True)
    out = SPOT_DIR / f"{area}_{start}_{end}_spot.nc"

    # ±0.5° around the point — a handful of 0.25° ERA5 cells. Must span >1 cell
    # per axis: a single-cell cutout makes atlite divide by (n_cells - 1) == 0.
    x = slice(lon - 0.5, lon + 0.5)
    y = slice(lat - 0.5, lat + 0.5)
    log.info(f"[{area} {year}] probing single cell x={x} y={y} full year -> {out}")

    cutout = atlite.Cutout(
        path=str(out),
        module="era5",
        x=x,
        y=y,
        time=slice(f"{year}-01-01", f"{year}-12-31 23:00"),
    )
    cutout.prepare(tmpdir=str(ATLITE_CACHE))

    report = check_cutout(out, start, end)
    log.info(f"[{area} {year}] QC report:\n{report.summary()}")
    return report.ok


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--areas", nargs="+", default=list(PROBE_POINTS))
    parser.add_argument("--years", nargs="+", type=int, required=True)
    args = parser.parse_args()

    results: dict[tuple[str, int], bool] = {}
    for area in args.areas:
        for year in args.years:
            try:
                results[(area, year)] = probe(area, year)
            except Exception as exc:  # download-stage failure is itself a result
                log.error(f"[{area} {year}] FAILED during download/prepare: {exc!r}")
                results[(area, year)] = False

    log.info("=== spot-check summary ===")
    for (area, year), ok in sorted(results.items()):
        log.info(f"  {area} {year}: {'OK' if ok else 'FAILED'}")

    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
