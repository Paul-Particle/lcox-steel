"""
Create one Atlite ERA5 cutout for a given area and date range.

Output: cutouts/{cf_area}_{start_date}_{end_date}.nc
"""
from pathlib import Path

import atlite
import geopandas as gpd

from common._paths import SHAPES_RES, ATLITE_CACHE, CUTOUTS
from scripts.res_cf._helpers import load_res_cf_cfg

if "snakemake" not in globals():
    from common._stubs import snakemake

REGIONS_PATH          = SHAPES_RES / "regions.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"


# Defaults for standalone use
CF_AREA     = "de"
START_DATE  = "20230101"
END_DATE    = "20231231"
OUTPUT_PATH = CUTOUTS / f"{CF_AREA}_{START_DATE}_{END_DATE}.nc"
RES_CF_CFG  = load_res_cf_cfg()

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    CF_AREA     = snakemake.wildcards.cf_area.lower()
    START_DATE  = snakemake.params.start_date
    END_DATE    = snakemake.params.end_date
    OUTPUT_PATH = Path(snakemake.output[0])
    RES_CF_CFG  = snakemake.config["res_cf"]

COUNTRIES_CFG = RES_CF_CFG["countries"]
BBOX_PAD_DEG  = RES_CF_CFG.get("cutout", {}).get("bbox_pad_deg", 1.0)


def iso(yyyymmdd: str) -> str:
    """YYYYMMDD -> YYYY-MM-DD."""
    return f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}"


def get_bounds(cf_area: str, pad: float = BBOX_PAD_DEG) -> tuple[slice, slice]:
    region_tag = COUNTRIES_CFG[cf_area]["region"]
    regions = gpd.read_parquet(REGIONS_PATH).to_crs(4326)
    geoms = list(regions.loc[regions["region"] == region_tag, "geometry"])
    if OFFSHORE_REGIONS_PATH.exists():
        offshore = gpd.read_parquet(OFFSHORE_REGIONS_PATH).to_crs(4326)
        geoms += list(offshore.loc[offshore["region"] == region_tag, "geometry"])
    geom = gpd.GeoSeries(geoms, crs=4326).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    ATLITE_CACHE.mkdir(parents=True, exist_ok=True)

    x, y = get_bounds(CF_AREA)
    # End at 23:00 so the full final day of hourly data is included.
    time_range = slice(iso(START_DATE), f"{iso(END_DATE)} 23:00")

    kwargs = dict(path=str(OUTPUT_PATH), module="era5", x=x, y=y, time=time_range)
    if COUNTRIES_CFG[CF_AREA].get("coarse"):
        kwargs.update(dx=0.5, dy=0.5)

    cutout = atlite.Cutout(**kwargs)
    cutout.prepare(tmpdir=str(ATLITE_CACHE), concurrent_requests=1)
    print(f"Prepared: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
