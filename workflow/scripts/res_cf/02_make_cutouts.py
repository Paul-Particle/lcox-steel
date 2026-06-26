"""
Create Atlite ERA5 cutouts for RES capacity factors.

Start small + stable (DE, Jan–Feb 2023) to avoid CDS/GRIB flakiness on Windows.
Output:
- data/cutouts/de_2023_q1.nc

Backup caching: if a backup cutout file `<cutout>_backup.nc` exists it is copied in place
of hitting CDS — useful to preserve expensive downloads across re-runs.
"""

import logging
import shutil
from pathlib import Path
import atlite
import geopandas as gpd

from common._paths import ATLITE_CACHE, CUTOUTS, SHAPES_RAW

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_NE_ZIP = SHAPES_RAW / "ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip"
_CF_AREA = "de"
_ISO3 = "DEU"
_MAINLAND_BBOX = None
_START_DATE = "20230101"
_END_DATE = "20231231"
_OUTPUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_COARSE = False
_REGION = "DE"
_BBOX_PAD_DEG = 1.0
_MONTHLY_REQUESTS = False
_ATLITE_CACHE = ATLITE_CACHE

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _NE_ZIP = Path(snakemake.input.ne_zip)
    _CF_AREA = snakemake.wildcards.cf_area
    _ISO3 = snakemake.params.iso3
    _MAINLAND_BBOX = snakemake.params.mainland_bbox
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _OUTPUT_PATH = Path(snakemake.output[0])
    _COARSE = snakemake.params.coarse
    _REGION = snakemake.params.region
    _BBOX_PAD_DEG = snakemake.params.bbox_pad_deg
    _MONTHLY_REQUESTS = snakemake.params.monthly_requests

def _iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


REGIONS = _ISO3 #gpd.read_file("data/shapes/regions.geojson").to_crs(4326)
OFFSHORE_REGIONS = None #gpd.read_file("data/shapes/offshore_regions.geojson").to_crs()  # TODO: restore offshore geometry in bounds_for

TMPDIR = _ATLITE_CACHE


def bounds_for(region_name=REGIONS, pad=1.0):
    land_geom = gpd.read_file(str(_NE_ZIP)).to_crs(4326)
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in land_geom.columns:
            selection = land_geom.loc[land_geom[col] == region_name, "geometry"]
            if not selection.empty:
                land_geom = selection.union_all()
                break
    else:
        raise ValueError(f"ISO3 code '{region_name}' not found in Natural Earth shapefile.")
    
    # if we have configured a mainland bounding box in the config, we need to remove overseas territories
    if _MAINLAND_BBOX is not None:
        lon_min, lon_max, lat_min, lat_max = _MAINLAND_BBOX
        parts = list(land_geom.geoms) if land_geom.geom_type == "MultiPolygon" else [land_geom]
        land_geom = [p for p in parts if lon_min <= p.centroid.x <= lon_max and lat_min <= p.centroid.y <= lat_max]
        if not land_geom:
            raise ValueError(f"No polygon parts inside mainland_bbox {_MAINLAND_BBOX}.")
        geom = gpd.GeoSeries(
            land_geom,
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

    x, y = bounds_for(region_name=REGIONS,pad=_BBOX_PAD_DEG)
    

    cutout = atlite.Cutout(
        path=str(_OUTPUT_PATH),
        module="era5",
        x=x,
        y=y,
        **({"dx": 0.5, 
            "dy": 0.5} if _COARSE else {}),
        time=slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00"), # End at 23:00 so the full final day of hourly data is included.
    )
    log.info(
        f"starting CDS request: x={x} y={y} time={time_range} "
        f"monthly={_MONTHLY_REQUESTS} coarse={_COARSE} → {_OUTPUT_PATH}"
    )
    cutout.prepare(tmpdir=str(TMPDIR), monthly_requests=_MONTHLY_REQUESTS)
    log.info(f"CDS request complete: {_OUTPUT_PATH}")


if __name__ == "__main__":
    main()

