"""
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

from pathlib import Path

import geopandas as gpd
import pandas as pd

if "snakemake" not in globals():
    from common._stubs import snakemake

# ---- paths ----
NE_SHP = Path("data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp")
OUT_GEOJSON = Path("resources/shapes/regions.geojson")

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    OUT_GEOJSON = Path(snakemake.output[0])


def one_country(world: gpd.GeoDataFrame, iso_a3: str) -> gpd.GeoSeries:
    # Natural Earth: ISO_A3 can be -99; ADM0_A3 / SOV_A3 are more reliable
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in world.columns:
            sel = world.loc[world[col] == iso_a3, "geometry"]
            if not sel.empty:
                # dissolve to single geometry
                return sel.union_all() if hasattr(sel, "union_all") else sel.unary_union

    raise ValueError(f"Code '{iso_a3}' not found in columns ADM0_A3/SOV_A3/ISO_A3.")

def fr_mainland_only(fr_geom):
    parts = list(fr_geom.geoms) if fr_geom.geom_type == "MultiPolygon" else [fr_geom]

    # keep only polygons whose centroid is in metropolitan Europe
    eu_parts = [
        p for p in parts
        if (-10 <= p.centroid.x <= 15) and (41 <= p.centroid.y <= 52)
    ]

    if not eu_parts:
        raise ValueError("Could not isolate metropolitan France polygon.")

    return gpd.GeoSeries(eu_parts, crs=4326).union_all()


def main() -> None:
    if not NE_SHP.exists():
        raise FileNotFoundError(
            f"Cannot find Natural Earth shapefile at: {NE_SHP}\n"
            "Download 'Admin 0 – Countries' (1:110m) and extract into:\n"
            "  data/shapes/ne_110m_admin_0_countries/"
        )

    OUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)

    world = gpd.read_file(NE_SHP).to_crs(4326)

    geoms = {
        "DE": one_country(world, "DEU"),
        "FR": fr_mainland_only(one_country(world, "FRA")), #Only European landmass
        "ES": one_country(world, "ESP"),
        "AUS": one_country(world, "AUS"),
        "BRA": one_country(world, "BRA"),
    }

    #de_lu = geoms["DE"].union(geoms["LU"])

    regions = gpd.GeoDataFrame(
        {
            "region": ["DE", "FR", "ES", "AUS", "BRA"],
            "geometry": [geoms["DE"], geoms["FR"], geoms["ES"], geoms["AUS"], geoms["BRA"]],
        },
        crs=4326,
    )

    # Fix occasional invalid geometries
    regions["geometry"] = regions["geometry"].buffer(0)

    regions.to_file(OUT_GEOJSON, driver="GeoJSON")
    print("Wrote:", OUT_GEOJSON)
    print(regions[["region"]])


if __name__ == "__main__":
    main()
