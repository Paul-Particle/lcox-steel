"""
Create an Atlite ERA5 cutout for one (area, start_date, end_date).

The cutout extent is the bounding box of the union of the onshore land geometry
(make_area_geometry) and the offshore zone geometry (make_offshore_geometry),
padded by `bbox_pad_deg`. Unioning in the offshore zone matters: it can reach up
to `offshore_max_distance_km` (~200 km) from the coast — far beyond a land-only
bbox — so without it the offshore-wind cells get clipped out of the cutout. An
area's `cutout_bbox` replaces that box (the demo area uses it to match the cutout
it ships).

Caching: the download is keyed on its actual request parameters (bbox, dx/dy,
time range) and reused from a persistent cache (common._cutout_cache) rather than
re-fetched from CDS. A request no entry matches exactly is cut out of the
smallest cached cutout, of any area, whose box covers it. The legacy `<cutout>_backup.nc` sibling is still honoured as
a fallback but is superseded by the keyed cache.

QC: every cutout (freshly downloaded, cached, or from a backup) is validated
before the rule succeeds, so a partial download or ERA5/ERA5T mix fails loudly
instead of flowing into the CF pipeline (common._cutout_qc).

CDS progress: the download is wrapped in common._cds_monitor so queue/running
status is logged periodically instead of the run appearing to hang.
"""

from pathlib import Path
import json
import logging
import shutil

import pandas as pd

from common._paths import ATLITE_CACHE, CUTOUTS, SHAPES_RES
from common._cutout_qc import ERA5_NATIVE_STEP_DEG, validate_cutout
from common._cutout_cache import (
    CACHE_DIR,
    cache_params,
    cache_paths,
    link_or_copy,
    store_in_cache,
    warn_if_low_disk,
)
from common._cds_monitor import cds_progress_logger
if "snakemake" not in globals():
    from common._stubs import snakemake
from common._logging import configure_logging
configure_logging(snakemake)
log = logging.getLogger(__name__)
# Standalone defaults
_REGIONS_PATH = SHAPES_RES / "de_geo.parquet"
_OFFSHORE_REGIONS_PATH = SHAPES_RES / "de_offshore_geo.parquet"
_AREA = "DE_LU"
_START_DATE = "20230101"
_END_DATE = "20231231"
_OUTPUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_COARSE = False
_BBOX_PAD_DEG = 1.0
_CUTOUT_BBOX = None
_MONTHLY_REQUESTS = False
_CDS_POLL_INTERVAL_S = 30.0
_CACHE_WARN_SIZE_GB = 100.0
_MIN_FREE_DISK_GB = 20.0
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _REGIONS_PATH = Path(snakemake.input.regions)
    _OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    _AREA = snakemake.wildcards.area
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _OUTPUT_PATH = Path(snakemake.output[0])
    _COARSE = snakemake.params.coarse
    _BBOX_PAD_DEG = snakemake.params.bbox_pad_deg
    _CUTOUT_BBOX = snakemake.params.cutout_bbox
    _MONTHLY_REQUESTS = snakemake.params.monthly_requests
    _CDS_POLL_INTERVAL_S = snakemake.params.cds_poll_interval_s
    _CACHE_WARN_SIZE_GB = snakemake.params.cache_warn_size_gb
    _MIN_FREE_DISK_GB = snakemake.params.min_free_disk_gb
def _iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"

from pathlib import Path
import atlite
import geopandas as gpd

REGIONS = _REGIONS_PATH
OFFSHORE_REGIONS = _OFFSHORE_REGIONS_PATH 

TMPDIR = ATLITE_CACHE


def bounds_for(pad=1.0):
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

    # Cutout bounds = bbox of the land ∪ offshore union, padded, unless the area
    # names its own. Computed up front because the cache key is derived from the
    # actual request parameters.
    if _CUTOUT_BBOX:
        lon_min, lon_max, lat_min, lat_max = _CUTOUT_BBOX
        x, y = slice(lon_min, lon_max), slice(lat_min, lat_max)
    else:
        x, y = bounds_for(pad=_BBOX_PAD_DEG)
    dx = dy = 0.5 if _COARSE else None
    params = cache_params("era5", x, y, dx, dy, _START_DATE, _END_DATE)
    cached_cutout, _ = cache_paths(_AREA, params)
    log.info(f"cutout bounds x={x} y={y}  cache key -> {cached_cutout.name}")

    # 1) Keyed cache: reuse an identical prior download.
    if cached_cutout.exists():
        log.info(f"cache hit: {cached_cutout} — materialising {_OUTPUT_PATH} (skipping CDS)")
        link_or_copy(cached_cutout, _OUTPUT_PATH)
        report = validate_cutout(_OUTPUT_PATH, _START_DATE, _END_DATE)
        log.info(f"cutout QC passed (from cache):\n{report.summary()}")
        return

    # 2) A larger cached cutout of the same window and resolution, of any area,
    #    cut down to this box. One grid step of slack, because atlite snaps a
    #    request to the ERA5 lattice inside its box, so a box a fraction of a step
    #    wider asks for no cell the cached one lacks.
    step = dx or ERA5_NATIVE_STEP_DEG
    cache_entries = pd.DataFrame(
        [{**json.loads(path.read_text()), "path": path.with_suffix(".nc")}
         for path in CACHE_DIR.glob("*.json") if path.with_suffix(".nc").exists()],
        columns=[*params, "path"],
    )
    covering = (
        cache_entries[
            (cache_entries["module"] == params["module"])
            & (cache_entries["start_date"] == params["start_date"])
            & (cache_entries["end_date"] == params["end_date"])
            & (cache_entries["dx"].fillna(0) == (params["dx"] or 0))
            & (cache_entries["x0"] <= params["x0"] + step)
            & (cache_entries["x1"] >= params["x1"] - step)
            & (cache_entries["y0"] <= params["y0"] + step)
            & (cache_entries["y1"] >= params["y1"] - step)
        ]
        .assign(size=lambda df: (df["x1"] - df["x0"]) * (df["y1"] - df["y0"]))
        .sort_values("size")
    )
    if not covering.empty:
        source = covering["path"].iloc[0]
        log.info(f"covered by cached {source.name} — cutting {_OUTPUT_PATH.name} out of it (skipping CDS)")
        # The same box under another area's name (Brazil's submarkets take the
        # national one) is linked rather than rewritten: these run to gigabytes.
        same_box = all(covering[bound].iloc[0] == params[bound] for bound in ("x0", "x1", "y0", "y1"))
        if same_box:
            link_or_copy(source, _OUTPUT_PATH)
        else:
            # The cut keeps the cached cutout's compression.
            atlite.Cutout(path=str(source)).sel(x=x, y=y).to_file(_OUTPUT_PATH)
        report = validate_cutout(_OUTPUT_PATH, _START_DATE, _END_DATE)
        log.info(f"cutout QC passed (cut from cache):\n{report.summary()}")
        store_in_cache(_OUTPUT_PATH, _AREA, params, warn_size_gb=_CACHE_WARN_SIZE_GB)
        return

    # 3) Legacy backup sibling (superseded by the keyed cache; kept as a fallback
    #    so previously pinned `_backup.nc` files still short-circuit CDS).
    backup = _OUTPUT_PATH.with_name(_OUTPUT_PATH.stem + "_backup" + _OUTPUT_PATH.suffix)
    if backup.exists():
        log.info(f"backup cutout found at {backup} — copying to {_OUTPUT_PATH} (skipping CDS)")
        with open(backup, "rb") as fsrc, open(_OUTPUT_PATH, "wb") as fdst:
            shutil.copyfileobj(fsrc, fdst)
        report = validate_cutout(_OUTPUT_PATH, _START_DATE, _END_DATE)
        log.info(f"cutout QC passed (from backup):\n{report.summary()}")
        cached = store_in_cache(_OUTPUT_PATH, _AREA, params, warn_size_gb=_CACHE_WARN_SIZE_GB)
        # Backups are user-managed pins (one is git-tracked for the demo), so we
        # never delete them — but once grandfathered into the cache they are
        # redundant and can be removed by hand.
        log.info(f"backup {backup.name} now cached as {cached.name}; it is redundant and may be removed manually")
        return

    # 4) Download from CDS.
    warn_if_low_disk(_MIN_FREE_DISK_GB)
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

    store_in_cache(_OUTPUT_PATH, _AREA, params, warn_size_gb=_CACHE_WARN_SIZE_GB)


if __name__ == "__main__":
    main()

