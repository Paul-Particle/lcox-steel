"""
Create one Atlite ERA5 cutout for a given area and date range.

Output: cutouts/{cf_area}_{start_date}_{end_date}.nc

Caching hack: if a sibling file named `<cutout>_backup.nc` exists, copy from
it instead of hitting CDS. This lets us preserve hand-curated or
expensively-downloaded cutouts across code-trigger reruns without giving up
the rule's reproducibility story. Drop this once res_cf inputs get the
same content-addressed caching treatment that the grid pipelines already have.

Caveat: the backup is copied without verification — partial or stale files
will silently propagate into downstream CFs. Operator's responsibility to
ensure the backup is a valid, complete cutout (e.g. the result of a
previous successful CDS run).
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


def get_bounds(pad: float) -> tuple[slice, slice]:
    world = gpd.read_file(str(_NE_ZIP)).to_crs(4326)
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in world.columns:
            sel = world.loc[world[col] == _ISO3, "geometry"]
            if not sel.empty:
                geom = sel.union_all()
                break
    else:
        raise ValueError(f"ISO3 code '{_ISO3}' not found in Natural Earth shapefile.")
    if _MAINLAND_BBOX is not None:
        lon_min, lon_max, lat_min, lat_max = _MAINLAND_BBOX
        parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
        kept = [p for p in parts if lon_min <= p.centroid.x <= lon_max and lat_min <= p.centroid.y <= lat_max]
        if not kept:
            raise ValueError(f"No polygon parts inside mainland_bbox {_MAINLAND_BBOX}.")
        geom = gpd.GeoSeries(kept, crs=4326).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main() -> None:
    _OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    _ATLITE_CACHE.mkdir(parents=True, exist_ok=True)

    backup = _OUTPUT_PATH.with_name(_OUTPUT_PATH.stem + "_backup" + _OUTPUT_PATH.suffix)
    if backup.exists():
        log.info(f"backup cutout found at {backup} — copying to {_OUTPUT_PATH} (skipping CDS)")
        shutil.copy2(backup, _OUTPUT_PATH)
        return

    x, y = get_bounds(_BBOX_PAD_DEG)
    # End at 23:00 so the full final day of hourly data is included.
    time_range = slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00")

    kwargs = dict(path=str(_OUTPUT_PATH), module="era5", x=x, y=y, time=time_range)
    if _COARSE:
        kwargs.update(dx=0.5, dy=0.5)

    cutout = atlite.Cutout(**kwargs)
    log.info(
        f"starting CDS request: x={kwargs['x']} y={kwargs['y']} time={time_range} "
        f"monthly={_MONTHLY_REQUESTS} coarse={_COARSE} → {_OUTPUT_PATH}"
    )
    # TODO: when CDS polling lands, spawn a background thread that calls
    # cdsapi `get_jobs()` every N seconds and `log.info`s status — see TODO.md.
    cutout.prepare(tmpdir=str(_ATLITE_CACHE), monthly_requests=_MONTHLY_REQUESTS)
    log.info(f"CDS request complete: {_OUTPUT_PATH}")


if __name__ == "__main__":
    main()
