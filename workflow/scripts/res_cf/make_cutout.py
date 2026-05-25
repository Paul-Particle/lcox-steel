"""
Create one Atlite ERA5 cutout for a given country, year, and quarter.

Output: cutouts/{country}_{year}_{quarter}.nc
"""
from pathlib import Path

import atlite
import geopandas as gpd

from common._paths import SHAPES_RES, ATLITE_CACHE
from scripts.res_cf._helpers import QUARTER_DATES, cutout_path, load_res_cf_cfg

if "snakemake" not in globals():
    from common._stubs import snakemake

REGIONS_PATH          = SHAPES_RES / "regions.geojson"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.geojson"


# Defaults for standalone use
COUNTRY     = "de"
YEAR        = 2023
QUARTER     = "q1"
OUTPUT_PATH = cutout_path(COUNTRY, YEAR, QUARTER)
RES_CF_CFG  = load_res_cf_cfg()

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRY     = snakemake.wildcards.country.lower()
    YEAR        = int(snakemake.wildcards.year)
    QUARTER     = snakemake.wildcards.quarter
    OUTPUT_PATH = Path(snakemake.output[0])
    RES_CF_CFG  = snakemake.config["res_cf"]

COUNTRIES_CFG = RES_CF_CFG["countries"]
BBOX_PAD_DEG  = RES_CF_CFG.get("cutout", {}).get("bbox_pad_deg", 1.0)


def get_bounds(country: str, pad: float = BBOX_PAD_DEG) -> tuple[slice, slice]:
    region_tag = COUNTRIES_CFG[country]["region"]
    regions = gpd.read_file(REGIONS_PATH).to_crs(4326)
    geoms = list(regions.loc[regions["region"] == region_tag, "geometry"])
    if OFFSHORE_REGIONS_PATH.exists():
        offshore = gpd.read_file(OFFSHORE_REGIONS_PATH).to_crs(4326)
        geoms += list(offshore.loc[offshore["region"] == region_tag, "geometry"])
    geom = gpd.GeoSeries(geoms, crs=4326).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    ATLITE_CACHE.mkdir(parents=True, exist_ok=True)

    x, y = get_bounds(COUNTRY)
    start_suffix, end_suffix = QUARTER_DATES[QUARTER]
    time_range = slice(f"{YEAR}{start_suffix}", f"{YEAR}{end_suffix}")

    kwargs = dict(path=str(OUTPUT_PATH), module="era5", x=x, y=y, time=time_range)
    if COUNTRIES_CFG[COUNTRY].get("coarse"):
        kwargs.update(dx=0.5, dy=0.5)

    cutout = atlite.Cutout(**kwargs)
    cutout.prepare(tmpdir=str(ATLITE_CACHE), concurrent_requests=1)
    print(f"Prepared: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
