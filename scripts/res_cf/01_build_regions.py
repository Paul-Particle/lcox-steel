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

# ---- paths ----
NE_SHP = Path("data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp")
OUT_GEOJSON = Path("data/shapes/regions.geojson")


def _one_country(world: gpd.GeoDataFrame, iso_a3: str) -> gpd.GeoSeries:
    # Natural Earth: ISO_A3 can be -99; ADM0_A3 / SOV_A3 are more reliable
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in world.columns:
            sel = world.loc[world[col] == iso_a3, "geometry"]
            if not sel.empty:
                # dissolve to single geometry
                return sel.union_all() if hasattr(sel, "union_all") else sel.unary_union

    raise ValueError(f"Code '{iso_a3}' not found in columns ADM0_A3/SOV_A3/ISO_A3.")



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
        "DE": _one_country(world, "DEU"),
        "FR": _one_country(world, "FRA"),
        "ES": _one_country(world, "ESP"),
        "AUS": _one_country(world, "AUS"),
        "BRA": _one_country(world, "BRA"),
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
