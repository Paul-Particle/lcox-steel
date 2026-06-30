"""Build a single country/region geometry for RES capacity factor extraction.

Source: Natural Earth Admin 0 countries (1:110m).

Output: resources/shapes/{cf_area}_geo.parquet — one row, columns (region, geometry).
If mainland_bbox is given, only polygon parts whose centroid falls inside
[lon_min, lon_max, lat_min, lat_max] are kept (drops overseas territories).

--- Previous module docstring (kept for reference) below ---

Build country/region geometries for RES capacity factor extraction with Atlite.

Source:
- Natural Earth Admin 0 countries (1:110m)

Output:
- data/shapes/regions.geojson

Regions:
- DE
- FR
- ES
- AUS
- BRA
"""

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry.base import BaseGeometry

# ---- paths ----
from common._paths import DATA, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
NE_SHP = DATA / "shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip"
OUT_GEOJSON = SHAPES_RES / "de_geo.parquet"
_CF_AREA = "de"
_ISO3 = "DEU"
_REGION = "DE"
_MAINLAND_BBOX = None

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    NE_SHP = Path(snakemake.input[0])
    OUT_GEOJSON = Path(snakemake.output[0])
    _CF_AREA = snakemake.wildcards.cf_area
    _ISO3 = snakemake.params.iso3
    _REGION = snakemake.params.region
    _MAINLAND_BBOX = snakemake.params.mainland_bbox


def _one_country(world: gpd.GeoDataFrame, iso_a3: str) -> BaseGeometry:
    # Natural Earth: ISO_A3 can be -99; ADM0_A3 / SOV_A3 are more reliable
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in world.columns:
            sel = world.loc[world[col] == iso_a3, "geometry"]
            if not sel.empty:
                # dissolve to single geometry
                return sel.union_all() if hasattr(sel, "union_all") else sel.unary_union

    raise ValueError(f"Code '{iso_a3}' not found in columns ADM0_A3/SOV_A3/ISO_A3.")
# generalzied version of France-only function
def restrict_to_bbox(region_geom: BaseGeometry, bbox: list[float]) -> BaseGeometry:
    """Return only the polygon parts whose centroid falls inside `bbox`.

    `bbox` is [lon_min, lon_max, lat_min, lat_max]; raises if nothing remains.
    """
    parts = list(region_geom.geoms) if region_geom.geom_type == "MultiPolygon" else [region_geom]
    lon_min, lon_max, lat_min, lat_max = bbox
    # Keep only polygons whose centroid is in [lon_min, lon_max, lat_min, lat_max].
    mainland_parts = [
        p for p in parts
        if lon_min <= p.centroid.x <= lon_max and lat_min <= p.centroid.y <= lat_max
    ]

    if not mainland_parts:
        raise ValueError(f"No polygon parts inside bbox {bbox}.")

    return gpd.GeoSeries(mainland_parts, crs=4326).union_all()


def main() -> None:
    if not NE_SHP.exists():
        raise FileNotFoundError(
            f"Cannot find Natural Earth shapefile at: {NE_SHP}\n"
            "Download 'Admin 0 – Countries' (1:110m) and extract into:\n"
            "  data/shapes/ne_110m_admin_0_countries/"
        )

    OUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)

    world = gpd.read_file(str(NE_SHP)).to_crs(4326)

    geom = _one_country(world, _ISO3)
    if _MAINLAND_BBOX is not None:
        geom = restrict_to_bbox(geom, _MAINLAND_BBOX)

    #de_lu = geoms["DE"].union(geoms["LU"])

    region = gpd.GeoDataFrame(
        {
            "region": [_REGION], 
            "geometry": [geom]
        }, 
        crs=4326
        )

    # Fix occasional invalid geometries
    region["geometry"] = region["geometry"].buffer(0)

    region.to_parquet(OUT_GEOJSON)
    log.info(f"wrote {OUT_GEOJSON}")
    log.info(f"({_CF_AREA} → {_REGION})")


if __name__ == "__main__":
    main()
