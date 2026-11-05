"""Build a sub-national area's offshore geometry: its share of the country's waters by nearest coast.

The country's offshore area is its EEZ within offshore_max_distance_km of its
land, as for a country. Every Admin 1 unit of the country then claims the water
closest to its own coast, and the area takes what its `admin1_names` claim. So
neighbouring areas never share a cell of sea, and water nearest a unit that is
no area's (Western Australia, Nunavut) belongs to none of them.

Output: resources/shapes/{area}_offshore_geo.parquet — one row, (region, geometry).
"""

import logging
from pathlib import Path

import geopandas as gpd
import shapely

from common._paths import DATA, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
ADMIN1_SHP = DATA / "shapes/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.zip"
EEZ_SHP = DATA / "shapes/offshore_zones/eez_v12.zip"
OUT_PARQUET = SHAPES_RES / "VIC1_offshore_geo.parquet"
_ISO3 = "AUS"
_REGION = "VIC"
_ADMIN1_NAMES = ["Victoria"]
_OFFSHORE_MAX_KM = 200.0
# How finely each unit's coast is sampled for the nearest-coast split.
COAST_SPACING_KM = 2.0

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    ADMIN1_SHP = Path(snakemake.input.admin1)
    EEZ_SHP = Path(snakemake.input.offshore_zone)
    OUT_PARQUET = Path(snakemake.output[0])
    _ISO3 = snakemake.params.iso3
    _REGION = snakemake.params.region
    _ADMIN1_NAMES = snakemake.params.admin1_names
    _OFFSHORE_MAX_KM = float(snakemake.params.offshore_max_distance_km)


def main() -> None:
    units = gpd.read_file(ADMIN1_SHP).to_crs(4326).query("adm0_a3 == @_ISO3")[["name", "geometry"]]
    # Distances in metres, true from the area's centre outwards, so they hold
    # along its coast however wide the country is.
    centre = units.query("name in @_ADMIN1_NAMES").union_all().centroid
    metric_crs = f"+proj=aeqd +lat_0={centre.y} +lon_0={centre.x} +units=m"
    units = units.to_crs(metric_crs)
    country_land = units.union_all()

    eez = (
        gpd.read_file(EEZ_SHP)
        .query("ISO_TER1 == @_ISO3 and POL_TYPE == '200NM'")
        .to_crs(metric_crs)
        .union_all()
    )
    country_offshore = eez.difference(country_land).intersection(
        country_land.buffer(_OFFSHORE_MAX_KM * 1000)
    )

    # Points along every unit's outline, each labelled with its unit, on a 100 m
    # lattice: GEOS cannot build the diagram from points closer than that. A
    # point two units share lies on their land border and keeps one label, which
    # only moves the split between those two units where it already runs.
    coast_points = (
        units.set_index("name").boundary
        .segmentize(COAST_SPACING_KM * 1000)
        .get_coordinates()
        .round(-2)
        .drop_duplicates()
    )
    cells = shapely.voronoi_polygons(
        shapely.multipoints(coast_points.to_numpy()),
        extend_to=country_offshore.envelope, ordered=True,
    )
    nearest_unit = gpd.GeoDataFrame(
        {"name": coast_points.index}, geometry=list(cells.geoms), crs=metric_crs
    )
    claimed = nearest_unit.query("name in @_ADMIN1_NAMES").union_all()
    offshore = country_offshore.intersection(claimed)

    gdf = gpd.GeoDataFrame({"region": [_REGION], "geometry": [offshore]}, crs=metric_crs)
    gdf["geometry"] = gdf.buffer(0)
    gdf.to_crs(4326).to_parquet(OUT_PARQUET)
    area_km2 = gdf.to_crs(6933).area.iloc[0] / 1e6
    log.info(f"wrote {OUT_PARQUET} ({_REGION}: {area_km2:.0f} km², "
             f"{offshore.area / country_offshore.area:.0%} of {_ISO3}'s offshore area)")


if __name__ == "__main__":
    main()
