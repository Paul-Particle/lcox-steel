"""
Create Atlite ERA5 cutouts for RES capacity factors.

Start small + stable (DE, Jan–Feb 2023) to avoid CDS/GRIB flakiness on Windows.
Output:
- data/cutouts/de_2023_q1.nc
"""
"""
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

from pathlib import Path
import atlite
import geopandas as gpd

REGIONS = _REGIONS_PATH
OFFSHORE_REGIONS = _OFFSHORE_REGIONS_PATH 

TMPDIR = ATLITE_CACHE


def bounds_for(region_names, pad=1.0):
    land_geom = gpd.read_parquet(REGIONS).to_crs(4326)
    offshore_geom = gpd.read_parquet(OFFSHORE_REGIONS).to_crs(4326)

    geom = gpd.GeoSeries(
        list(land_geom.geometry) + list(offshore_geom.geometry),
        crs=4326,
    ).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main():
    _OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)
    backup = _OUTPUT_PATH.with_name(_OUTPUT_PATH.stem + "_backup" + _OUTPUT_PATH.suffix)
    if backup.exists():
        log.info(f"backup cutout found at {backup} — copying to {_OUTPUT_PATH} (skipping CDS)")
        with open(backup, 'rb') as fsrc, open(_OUTPUT_PATH, 'wb') as fdst:
            shutil.copyfileobj(fsrc, fdst)
        return

    # ✅ start with ONE country/annual time slice
    region_names = REGIONS #[AUS]
    x, y = bounds_for(region_names, pad=_BBOX_PAD_DEG)
    log.info("Cutout bounds:")
    log.info(f"x ={x}")
    log.info(f"y ={y}")

    cutout = atlite.Cutout(
        path=str(_OUTPUT_PATH),
        module="era5",
        x=x,
        y=y,
        **({"dx": 0.5,
            "dy": 0.5} if _COARSE else {}),
        time=slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00"),# End at 23:00 so the full final day of hourly data is included.
    )
    time_range = slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00")
    log.info(f"starting CDS request: x={x} y={y} time={time_range} " f"monthly={_MONTHLY_REQUESTS} coarse={_COARSE} → {_OUTPUT_PATH}")

    # ✅ reduce parallelism + use fixed tmpdir (Windows-safe)
    cutout.prepare(tmpdir=str(TMPDIR), monthly_requests=_MONTHLY_REQUESTS)
    log.info(f"CDS request complete: {_OUTPUT_PATH}")

if __name__ == "__main__":
    main()

