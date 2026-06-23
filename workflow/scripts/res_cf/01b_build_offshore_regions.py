"""Build offshore region geometry from EEZ shapefile for one area.

Takes the country's EEZ polygon, subtracts the onshore land, and clips to
a buffer of offshore_max_distance_km km around the coast — the area where
offshore wind makes economic sense.

Output: resources/shapes/{cf_area}_offshore_geo.parquet — one row, (region, geometry).
"""

import logging
from pathlib import Path

import geopandas as gpd

from common._paths import DATA, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_OFFSHORE_ZONE_ZIP = DATA / "shapes/offshore_zones/eez_v12.zip"
_LAND_REGIONS = SHAPES_RES / "de_geo.parquet"
_OUT = SHAPES_RES / "de_offshore_geo.parquet"
_ISO3 = "DEU"
_REGION = "DE"
_OFFSHORE_MAX_KM = 200.0

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _OFFSHORE_ZONE_ZIP = Path(snakemake.input.offshore_zone)
    _LAND_REGIONS = Path(snakemake.input.regions)
    _OUT = Path(snakemake.output[0])
    _ISO3 = snakemake.params.iso3
    _REGION = snakemake.params.region
    _OFFSHORE_MAX_KM = float(snakemake.params.offshore_max_distance_km)


def main() -> None:
    """Clip the area's EEZ to a near-shore band minus land, and write the offshore parquet."""
    if not _OFFSHORE_ZONE_ZIP.exists():
        raise FileNotFoundError(f"Offshore zone ZIP not found: {_OFFSHORE_ZONE_ZIP}")
    if not _LAND_REGIONS.exists():
        raise FileNotFoundError(f"Land regions file not found: {_LAND_REGIONS}")

    _OUT.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"building offshore region for {_REGION} (ISO3={_ISO3}, max {_OFFSHORE_MAX_KM:.0f} km)")

    offshore_zone = gpd.read_file(str(_OFFSHORE_ZONE_ZIP)).to_crs(4326)
    gdf = offshore_zone.loc[
        (offshore_zone["ISO_TER1"] == _ISO3) & (offshore_zone["POL_TYPE"] == "200NM"),
        ["ISO_TER1", "geometry"],
    ].copy()
    gdf["region"] = _REGION
    gdf = gdf[["region", "geometry"]].dissolve(by="region", as_index=False)

    land = gpd.read_parquet(_LAND_REGIONS).to_crs(4326)

    land_m = land.to_crs(6933).copy()
    land_m["geometry"] = land_m["geometry"].buffer(_OFFSHORE_MAX_KM * 1000)
    land_buffer = land_m.to_crs(4326)

    offshore_zone_geom = gdf.loc[gdf["region"] == _REGION, "geometry"].iloc[0]
    land_geom = land.loc[land["region"] == _REGION, "geometry"].iloc[0]
    buf_geom = land_buffer.loc[land_buffer["region"] == _REGION, "geometry"].iloc[0]

    offshore_geom = offshore_zone_geom.difference(land_geom).intersection(buf_geom)

    result = gpd.GeoDataFrame({"region": [_REGION], "geometry": [offshore_geom]}, crs=4326)
    result["geometry"] = result["geometry"].buffer(0)
    if not result.is_valid.all():
        log.warning("offshore geometry is invalid after cleanup")

    result.to_parquet(_OUT)
    area_km2 = result.to_crs(6933)["geometry"].area.iloc[0] / 1e6
    log.info(f"wrote {_OUT} ({_REGION}: {area_km2:.0f} km²)")


if __name__ == "__main__":
    main()
