"""
Create an Atlite ERA5 cutout for one (cf_area, start_date, end_date).

The cutout extent is the bounding box of the union of the onshore land geometry
(01_build_regions) and the offshore zone geometry (01b_build_offshore_regions),
padded by `bbox_pad_deg`. Unioning in the offshore zone matters: it can reach up
to `offshore_max_distance_km` (~200 km) from the coast — far beyond a land-only
bbox — so without it the offshore-wind cells get clipped out of the cutout.

Caching: the download is keyed on its actual request parameters (bbox, dx/dy,
time range) and reused from a persistent cache (common._cutout_cache) rather than
re-fetched from CDS. The legacy `<cutout>_backup.nc` sibling is still honoured as
a fallback but is superseded by the keyed cache.

QC: every cutout (freshly downloaded, cached, or from a backup) is validated
before the rule succeeds, so a partial download or ERA5/ERA5T mix fails loudly
instead of flowing into the CF pipeline (common._cutout_qc).

CDS progress: the download is wrapped in common._cds_monitor so queue/running
status is logged periodically instead of the run appearing to hang.
"""

from pathlib import Path
import logging
import shutil
from common._paths import ATLITE_CACHE, CUTOUTS, SHAPES_RES
from common._cutout_qc import validate_cutout
from common._cutout_cache import cache_params, cache_paths, link_or_copy, store_in_cache
from common._cds_monitor import cds_progress_logger
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
_CDS_POLL_INTERVAL_S = 30.0
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
    _CDS_POLL_INTERVAL_S = snakemake.params.cds_poll_interval_s
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

    # Cutout bounds = bbox of the land ∪ offshore union, padded. Computed up front
    # because the cache key is derived from the actual request parameters.
    x, y = bounds_for(REGIONS, pad=_BBOX_PAD_DEG)
    dx = dy = 0.5 if _COARSE else None
    params = cache_params("era5", x, y, dx, dy, _START_DATE, _END_DATE)
    cached_cutout, _ = cache_paths(_CF_AREA, params)
    log.info(f"cutout bounds x={x} y={y}  cache key -> {cached_cutout.name}")

    # 1) Keyed cache: reuse an identical prior download.
    if cached_cutout.exists():
        log.info(f"cache hit: {cached_cutout} — materialising {_OUTPUT_PATH} (skipping CDS)")
        link_or_copy(cached_cutout, _OUTPUT_PATH)
        report = validate_cutout(_OUTPUT_PATH, _START_DATE, _END_DATE)
        log.info(f"cutout QC passed (from cache):\n{report.summary()}")
        return

    # 2) Legacy backup sibling (superseded by the keyed cache; kept as a fallback
    #    so previously pinned `_backup.nc` files still short-circuit CDS).
    backup = _OUTPUT_PATH.with_name(_OUTPUT_PATH.stem + "_backup" + _OUTPUT_PATH.suffix)
    if backup.exists():
        log.info(f"backup cutout found at {backup} — copying to {_OUTPUT_PATH} (skipping CDS)")
        with open(backup, "rb") as fsrc, open(_OUTPUT_PATH, "wb") as fdst:
            shutil.copyfileobj(fsrc, fdst)
        report = validate_cutout(_OUTPUT_PATH, _START_DATE, _END_DATE)
        log.info(f"cutout QC passed (from backup):\n{report.summary()}")
        store_in_cache(_OUTPUT_PATH, _CF_AREA, params)
        return

    # 3) Download from CDS.
    time_range = slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00")  # 23:00 keeps the full final day
    cutout = atlite.Cutout(
        path=str(_OUTPUT_PATH),
        module="era5",
        x=x,
        y=y,
        **({"dx": dx, "dy": dy} if _COARSE else {}),
        time=time_range,
    )
    log.info(f"starting CDS request: x={x} y={y} time={time_range} "
             f"monthly={_MONTHLY_REQUESTS} coarse={_COARSE} → {_OUTPUT_PATH}")

    with cds_progress_logger(interval_s=_CDS_POLL_INTERVAL_S):
        cutout.prepare(tmpdir=str(TMPDIR), monthly_requests=_MONTHLY_REQUESTS)
    log.info(f"CDS request complete: {_OUTPUT_PATH}")

    # Gate: a partial download or ERA5/ERA5T mix must fail the rule loudly rather
    # than flow silently into the CF pipeline. Raises CutoutQCError on failure,
    # which makes Snakemake discard the broken output.
    report = validate_cutout(_OUTPUT_PATH, _START_DATE, _END_DATE)
    log.info(f"cutout QC passed:\n{report.summary()}")

    store_in_cache(_OUTPUT_PATH, _CF_AREA, params)


if __name__ == "__main__":
    main()

