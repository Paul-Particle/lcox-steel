"""
Create one Atlite ERA5 cutout for a given area and date range.

Output: cutouts/{cf_area}_{start_date}_{end_date}.nc
"""

from pathlib import Path

import atlite
import geopandas as gpd

from common._paths import ATLITE_CACHE, CUTOUTS, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

# Standalone defaults
_REGIONS_PATH = SHAPES_RES / "de_geo.parquet"
_CF_AREA = "de"
_START_DATE = "20230101"
_END_DATE = "20231231"
_OUTPUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_COARSE = False
_REGION = "DE"
_BBOX_PAD_DEG = 1.0
_ATLITE_CACHE = ATLITE_CACHE

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _REGIONS_PATH = Path(snakemake.input.regions)
    _CF_AREA = snakemake.wildcards.cf_area
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _OUTPUT_PATH = Path(snakemake.output[0])
    _COARSE = snakemake.params.coarse
    _REGION = snakemake.params.region
    _BBOX_PAD_DEG = snakemake.params.bbox_pad_deg


def _iso(yyyymmdd: str) -> str:
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def get_bounds(pad: float) -> tuple[slice, slice]:
    regions = gpd.read_parquet(_REGIONS_PATH).to_crs(4326)
    geom = regions.loc[regions["region"] == _REGION, "geometry"].union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main() -> None:
    _OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    _ATLITE_CACHE.mkdir(parents=True, exist_ok=True)

    x, y = get_bounds(_BBOX_PAD_DEG)
    # End at 23:00 so the full final day of hourly data is included.
    time_range = slice(_iso(_START_DATE), f"{_iso(_END_DATE)} 23:00")

    kwargs = dict(path=str(_OUTPUT_PATH), module="era5", x=x, y=y, time=time_range)
    if _COARSE:
        kwargs.update(dx=0.5, dy=0.5)

    cutout = atlite.Cutout(**kwargs)
    cutout.prepare(tmpdir=str(_ATLITE_CACHE), concurrent_requests=1)
    print(f"Prepared: {_OUTPUT_PATH}")


if __name__ == "__main__":
    main()
