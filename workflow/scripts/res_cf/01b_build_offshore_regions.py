"""Build offshore region geometry from EEZ shapefile for one area.

Takes the country's EEZ polygon, subtracts the onshore land, and clips to
a buffer of offshore_max_distance_km km around the coast — the area where
offshore wind makes economic sense.

Output: resources/shapes/{cf_area}_offshore_geo.parquet — one row, (region, geometry).
"""

from pathlib import Path
import geopandas as gpd
from common._paths import DATA, SHAPES_RES



#CONFIG_PATH = Path("config_hannah.yaml")
EEZ_SHP = DATA / "shapes/offshore_zones/eez_v12.zip"
OUT_GEOJSON = SHAPES_RES / "de_offshore_geo.parquet"
LAND_REGIONS = SHAPES_RES / "de_geo.parquet"


import logging

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_ISO3 = "DEU"
_REGION = "DE"
_OFFSHORE_MAX_KM = 200.0

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _OFFSHORE_ZONE_ZIP = Path(snakemake.input.offshore_zone)
    LAND_REGIONS = Path(snakemake.input.regions)
    OUT_GEOJSON = Path(snakemake.output[0])
    _ISO3 = snakemake.params.iso3
    _REGION = snakemake.params.region
    _OFFSHORE_MAX_KM = float(snakemake.params.offshore_max_distance_km)



def main():

    max_distance_km = _OFFSHORE_MAX_KM
    log.info(f"Configured offshore distance limit: {max_distance_km}")
    if not EEZ_SHP.exists():
        raise FileNotFoundError(f"Offshore zone ZIP not found: {EEZ_SHP}")

    if not LAND_REGIONS.exists():
        raise FileNotFoundError(f"Land regions file not found: {LAND_REGIONS}")

    OUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"building offshore region for {_REGION} (ISO3={_ISO3}, max {max_distance_km:.0f} km)")

    eez = gpd.read_file(str(EEZ_SHP)).to_crs(4326)

    gdf = eez.loc[
        (eez["ISO_TER1"] == _ISO3) &
        (eez["POL_TYPE"] == "200NM"),
        ["ISO_TER1", "geometry"],
    ].copy()

    gdf["region"] = _REGION
    gdf = gdf[["region", "geometry"]]
    gdf = gdf.dissolve(by="region", as_index=False)


    land = gpd.read_parquet(LAND_REGIONS).to_crs(4326)

    land_m = land.to_crs(6933).copy()
    land_m["geometry"] = land_m["geometry"].buffer(max_distance_km * 1000)
    land_buffer = land_m.to_crs(4326)
    offshore_zone_geom = gdf.loc[gdf["region"] == _REGION, "geometry"].iloc[0]

    land_geom = land.loc[land["region"] == _REGION, "geometry"].iloc[0]

    buf_geom = land_buffer.loc[land_buffer["region"] == _REGION, "geometry"].iloc[0]

    offshore_geom = (
            offshore_zone_geom.difference(land_geom)
            .intersection(buf_geom)
    )

    gdf = gpd.GeoDataFrame({"region": [_REGION], "geometry": [offshore_geom]}, crs=4326)
    gdf["geometry"] = gdf["geometry"].buffer(0)
    if not gdf.is_valid.all():
        log.warning("offshore geometry is invalid after cleanup")
    # -----------------------------------

    gdf.to_parquet(OUT_GEOJSON)
    area_km2 = gdf.to_crs(6933)["geometry"].area.iloc[0] / 1e6

    log.info(f"wrote {OUT_GEOJSON} ({_REGION}: "
             f"{area_km2:.0f} km²)")



if __name__ == "__main__":
    main()
