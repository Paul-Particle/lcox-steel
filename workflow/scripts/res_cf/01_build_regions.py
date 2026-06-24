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
from shapely.geometry.base import BaseGeometry

from common._paths import DATA, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_NE_ZIP = DATA / "shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip"
_CF_AREA = "de"
_ISO3 = "DEU"
_REGION = "DE"
_MAINLAND_BBOX = None
_OUT = SHAPES_RES / "de_geo.parquet"

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _NE_ZIP = Path(snakemake.input[0])
    _CF_AREA = snakemake.wildcards.cf_area
    _ISO3 = snakemake.params.iso3
    _REGION = snakemake.params.region
    _MAINLAND_BBOX = snakemake.params.mainland_bbox
    _OUT = Path(snakemake.output[0])


def country_geometry(world: gpd.GeoDataFrame, iso_a3: str) -> BaseGeometry:
    """Return the dissolved geometry for `iso_a3`, trying ADM0_A3/SOV_A3/ISO_A3 in turn."""
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in world.columns:
            sel = world.loc[world[col] == iso_a3, "geometry"]
            if not sel.empty:
                return sel.union_all() if hasattr(sel, "union_all") else sel.unary_union
    raise ValueError(f"Code '{iso_a3}' not found in columns ADM0_A3/SOV_A3/ISO_A3.")


def restrict_to_bbox(geom: BaseGeometry, bbox: list[float]) -> BaseGeometry:
    """Return only the polygon parts whose centroid falls inside `bbox`.

    `bbox` is [lon_min, lon_max, lat_min, lat_max]; raises if nothing remains.
    """
    # --- Previous docstring (kept for reference) below ---
    # Keep only polygon parts whose centroid falls in [lon_min, lon_max, lat_min, lat_max].
    lon_min, lon_max, lat_min, lat_max = bbox
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    kept = [
        p for p in parts
        if lon_min <= p.centroid.x <= lon_max and lat_min <= p.centroid.y <= lat_max
    ]
    if not kept:
        raise ValueError(f"No polygon parts inside bbox {bbox}.")
    return gpd.GeoSeries(kept, crs=4326).union_all()


def main() -> None:
    """Load Natural Earth, extract the area geometry (optionally bbox-clipped), write parquet."""
    if not _NE_ZIP.exists():
        raise FileNotFoundError(f"Natural Earth shapefile ZIP not found: {_NE_ZIP}")

    _OUT.parent.mkdir(parents=True, exist_ok=True)

    world = gpd.read_file(str(_NE_ZIP)).to_crs(4326)
    geom = country_geometry(world, _ISO3)
    if _MAINLAND_BBOX is not None:
        geom = restrict_to_bbox(geom, _MAINLAND_BBOX)

    out = gpd.GeoDataFrame({"region": [_REGION], "geometry": [geom]}, crs=4326)
    out["geometry"] = out["geometry"].buffer(0)
    out.to_parquet(_OUT)
    log.info(f"wrote {_OUT} ({_CF_AREA} → {_REGION})")


if __name__ == "__main__":
    main()
