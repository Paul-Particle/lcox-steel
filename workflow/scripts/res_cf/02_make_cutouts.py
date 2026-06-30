"""
Create Atlite ERA5 cutouts for RES capacity factors.

The cutout extent is the bounding box of the union of the onshore land
geometry (01_build_regions) and the offshore zone geometry
(01b_build_offshore_regions), padded by `bbox_pad_deg`. Unioning in the
offshore zone matters: it can reach up to `offshore_max_distance_km` (~200 km)
from the coast — far beyond a land-only bbox + 1° pad — so without it the
offshore-wind cells get clipped out of the cutout.

Backup caching: if a backup cutout file `<cutout>_backup.nc` exists it is copied in place
of hitting CDS — useful to preserve expensive downloads across re-runs.
"""

import logging
import shutil
from pathlib import Path
import atlite
import geopandas as gpd

from common._paths import ATLITE_CACHE, CUTOUTS, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_REGIONS_PATH = SHAPES_RES / "de_geo.parquet"
_OFFSHORE_REGIONS_PATH = SHAPES_RES / "de_offshore_geo.parquet"
_CF_AREA = "de"
_START_DATE = "20230101"
_END_DATE = "20231231"
_OUTPUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_COARSE = False
_BBOX_PAD_DEG = 1.0
_MONTHLY_REQUESTS = False
_ATLITE_CACHE = ATLITE_CACHE

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _REGIONS_PATH = Path(snakemake.input.regions)
    _OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    _CF_AREA = snakemake.wildcards.cf_area
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _OUTPUT_PATH = Path(snakemake.output[0])
    _COARSE = snakemake.params.coarse
    _BBOX_PAD_DEG = snakemake.params.bbox_pad_deg
    _MONTHLY_REQUESTS = snakemake.params.monthly_requests


def _iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


TMPDIR = _ATLITE_CACHE


def bounds_for(pad: float = 1.0):
    """Cutout x/y slices: bbox of the land ∪ offshore union, padded by `pad` degrees.

    Both parquets hold a single dissolved geometry for this area (the
    mainland_bbox filter and EEZ clip are already applied upstream by 01/01b),
    so we just union their geometry columns.
    """
    land = gpd.read_parquet(_REGIONS_PATH).to_crs(4326)
    offshore = gpd.read_parquet(_OFFSHORE_REGIONS_PATH).to_crs(4326)

    geom = gpd.GeoSeries(
        list(land.geometry) + list(offshore.geometry),
        crs=4326,
    ).union_all()

    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main():
    _OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)

    # check for backup/cache before we go any further
    backup = _OUTPUT_PATH.with_name(_OUTPUT_PATH.stem + "_backup" + _OUTPUT_PATH.suffix)
    if backup.exists():
        log.info(f"backup cutout found at {backup} — copying to {_OUTPUT_PATH} (skipping CDS)")
        with open(backup, 'rb') as fsrc, open(_OUTPUT_PATH, 'wb') as fdst:
            shutil.copyfileobj(fsrc, fdst)
        return

    x, y = bounds_for(pad=_BBOX_PAD_DEG)

    # End at 23:00 so the full final day of hourly data is included.
    time_range = slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00")

    cutout = atlite.Cutout(
        path=str(_OUTPUT_PATH),
        module="era5",
        x=x,
        y=y,
        **({"dx": 0.5,
            "dy": 0.5} if _COARSE else {}),
        time=time_range,
    )
    log.info(
        f"starting CDS request: x={x} y={y} time={time_range} "
        f"monthly={_MONTHLY_REQUESTS} coarse={_COARSE} → {_OUTPUT_PATH}"
    )
    cutout.prepare(tmpdir=str(TMPDIR), monthly_requests=_MONTHLY_REQUESTS)
    log.info(f"CDS request complete: {_OUTPUT_PATH}")


if __name__ == "__main__":
    main()
