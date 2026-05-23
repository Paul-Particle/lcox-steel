"""
Create one Atlite ERA5 cutout for a given country, year, and quarter.

Output: cutouts/{country}_{year}_{quarter}.nc
"""
from pathlib import Path
import atlite
import geopandas as gpd

if "snakemake" not in globals():
    from common._stubs import snakemake

REGIONS_PATH         = Path("resources/shapes/regions.geojson")
OFFSHORE_REGIONS_PATH = Path("resources/shapes/offshore_regions.geojson")
TMPDIR = Path(".atlite-cache")

# Countries that need coarser resolution to avoid CDS instability / memory limits
COARSE_RESOLUTION_COUNTRIES = {"aus", "bra"}

QUARTER_DATES = {
    "q1": ("-01-01", "-03-31 23:00"),
    "q2": ("-04-01", "-06-30 23:00"),
    "q3": ("-07-01", "-09-30 23:00"),
    "q4": ("-10-01", "-12-31 23:00"),
}

# Defaults for standalone use
COUNTRY = "de"
YEAR    = 2023
QUARTER = "q1"
OUTPUT_PATH = Path(f"cutouts/{COUNTRY}_{YEAR}_{QUARTER}.nc")

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRY     = snakemake.wildcards.country.lower()
    YEAR        = int(snakemake.wildcards.year)
    QUARTER     = snakemake.wildcards.quarter
    OUTPUT_PATH = Path(snakemake.output[0])


def get_bounds(country_lower: str, pad: float = 1.0):
    country_upper = country_lower.upper()
    regions = gpd.read_file(REGIONS_PATH).to_crs(4326)
    geoms = list(regions.loc[regions["region"] == country_upper, "geometry"])
    if OFFSHORE_REGIONS_PATH.exists():
        offshore = gpd.read_file(OFFSHORE_REGIONS_PATH).to_crs(4326)
        geoms += list(offshore.loc[offshore["region"] == country_upper, "geometry"])
    geom = gpd.GeoSeries(geoms, crs=4326).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main():
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)

    x, y = get_bounds(COUNTRY)
    start_suffix, end_suffix = QUARTER_DATES[QUARTER]
    time_range = slice(f"{YEAR}{start_suffix}", f"{YEAR}{end_suffix}")

    kwargs = dict(path=str(OUTPUT_PATH), module="era5", x=x, y=y, time=time_range)
    if COUNTRY in COARSE_RESOLUTION_COUNTRIES:
        kwargs.update(dx=0.5, dy=0.5)

    cutout = atlite.Cutout(**kwargs)
    cutout.prepare(tmpdir=str(TMPDIR), concurrent_requests=1)
    print(f"Prepared: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
