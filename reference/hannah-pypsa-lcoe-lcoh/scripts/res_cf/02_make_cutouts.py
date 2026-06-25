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
OFFSHORE_REGIONS = gpd.read_file("data/shapes/offshore_regions.geojson").to_crs(4326)

TMPDIR = Path("data/tmp/atlite")


def bounds_for(region_names, pad=1.0):
    land_geom = REGIONS.loc[REGIONS["region"].isin(region_names), "geometry"]
    offshore_geom = OFFSHORE_REGIONS.loc[OFFSHORE_REGIONS["region"].isin(region_names), "geometry"]

    geom = gpd.GeoSeries(
        list(land_geom) + list(offshore_geom),
        crs=4326,
    ).union_all()
    minx, miny, maxx, maxy = geom.bounds
    return slice(minx - pad, maxx + pad), slice(miny - pad, maxy + pad)


def main():
    Path("data/cutouts").mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)

    # ✅ start with ONE country/Quaterly time slice
    region_names = ["AUS"]
    x, y = bounds_for(region_names)
    print("Cutout bounds:")
    print("x =", x)
    print("y =", y)

    cutout = atlite.Cutout(
        path="data/cutouts/aus_2023_q1_dxdy.nc",
        module="era5",
        x=x,
        y=y,
        #dx=0.5, #coarser grid if needed
        #dy=0.5, #coarser grid if needed
        time=slice("2023-01-01", "2023-03-31 23:00"),  # Q1 "2023-01-01", "2023-03-31 23:00", Q2 "2023-04-01", "2023-06-30 23:00", Q3 "2023-07-01", "2023-09-30 23:00", Q4 "2023-10-01", "2023-12-31 23:00"
    )

    # ✅ reduce parallelism + use fixed tmpdir (Windows-safe)
    cutout.prepare(tmpdir=str(TMPDIR), concurrent_requests=1)


if __name__ == "__main__":
    main()

