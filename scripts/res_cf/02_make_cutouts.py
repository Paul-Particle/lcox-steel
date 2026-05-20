"""
Create Atlite ERA5 cutouts for RES capacity factors.

Start small + stable (DE, Jan–Feb 2023) to avoid CDS/GRIB flakiness on Windows.
Output:
- data/cutouts/de_2023_q1.nc
"""

from pathlib import Path
import atlite
import geopandas as gpd

REGIONS = gpd.read_file("data/shapes/regions.geojson").to_crs(4326)
_offshore_path = Path("data/shapes/offshore_regions.geojson")
OFFSHORE_REGIONS = gpd.read_file(_offshore_path).to_crs(4326) if _offshore_path.exists() else None

TMPDIR = Path("data/tmp/atlite")


def bounds_for(region_names, pad=1.0):
    land_geom = REGIONS.loc[REGIONS["region"].isin(region_names), "geometry"]
    geoms = list(land_geom)
    if OFFSHORE_REGIONS is not None:
        offshore_geom = OFFSHORE_REGIONS.loc[OFFSHORE_REGIONS["region"].isin(region_names), "geometry"]
        geoms += list(offshore_geom)

    geom = gpd.GeoSeries(geoms, crs=4326).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main():
    Path("data/cutouts").mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)

    # diagnostic: DE, 2 weeks — minimal download for testing
    region_names = ["DE"]
    x, y = bounds_for(region_names)
    print("Cutout bounds:")
    print("x =", x)
    print("y =", y)

    cutout = atlite.Cutout(
        path="data/cutouts/de_2023_jan2w.nc",
        module="era5",
        x=x,
        y=y,
        time=slice("2023-01-01", "2023-01-14 23:00"),
    )

    # ✅ reduce parallelism + use fixed tmpdir (Windows-safe)
    cutout.prepare(tmpdir=str(TMPDIR), concurrent_requests=1)


if __name__ == "__main__":
    main()

