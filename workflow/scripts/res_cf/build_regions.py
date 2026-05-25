"""
Build country/region geometries for RES capacity factor extraction with Atlite.

Source: Natural Earth Admin 0 countries (1:110m).

Output: resources/shapes/regions.parquet — one row per enabled country in
config.yaml, columns = (region, geometry). GeoParquet preserves CRS + geometry
metadata natively. The "region" tag is the per-country field
config["res_cf"]["countries"][cc]["region"].

If a country has a mainland_bbox, only polygon parts whose centroid falls
inside the box are kept (drops overseas territories — e.g. metropolitan France).
"""

import yaml
from pathlib import Path

import geopandas as gpd
from shapely.geometry.base import BaseGeometry

from common._paths import DATA, SHAPES_RES, REPO_ROOT

if "snakemake" not in globals():
    from common._stubs import snakemake

NATURAL_EARTH_SHAPE = DATA / "shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
OUT_PATH = SHAPES_RES / "regions.parquet"

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    OUT_PATH = Path(snakemake.output[0])


def country_geometry(world: gpd.GeoDataFrame, iso_a3: str) -> BaseGeometry:
    # Natural Earth: ISO_A3 can be -99; ADM0_A3 / SOV_A3 are more reliable
    for col in ["ADM0_A3", "SOV_A3", "ISO_A3"]:
        if col in world.columns:
            sel = world.loc[world[col] == iso_a3, "geometry"]
            if not sel.empty:
                return sel.union_all() if hasattr(sel, "union_all") else sel.unary_union
    raise ValueError(f"Code '{iso_a3}' not found in columns ADM0_A3/SOV_A3/ISO_A3.")


def restrict_to_bbox(geom: BaseGeometry, bbox: list[float]) -> BaseGeometry:
    """Keep only polygon parts whose centroid falls in [lon_min, lon_max, lat_min, lat_max]."""
    lon_min, lon_max, lat_min, lat_max = bbox
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    kept = [
        p for p in parts
        if lon_min <= p.centroid.x <= lon_max and lat_min <= p.centroid.y <= lat_max
    ]
    if not kept:
        raise ValueError(f"No polygon parts inside bbox {bbox}.")
    return gpd.GeoSeries(kept, crs=4326).union_all()


def load_enabled_countries() -> dict[str, dict]:
    """Return dict of enabled country entries keyed by lowercase ISO-2."""
    with open(REPO_ROOT / "config/config.yaml") as f:
        cfg = yaml.safe_load(f)
    return {cc: info for cc, info in cfg["res_cf"]["countries"].items() if info.get("enabled")}


def main() -> None:
    if not NATURAL_EARTH_SHAPE.exists():
        raise FileNotFoundError(
            f"Cannot find Natural Earth shapefile at: {NATURAL_EARTH_SHAPE}\n"
            "Run scripts/res_cf/check_external_data.py for setup instructions."
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    if "snakemake" in globals() and hasattr(snakemake, "config"):
        countries = {cc: info for cc, info in snakemake.config["res_cf"]["countries"].items() if info.get("enabled")}
    else:
        countries = load_enabled_countries()

    world = gpd.read_file(NATURAL_EARTH_SHAPE).to_crs(4326)

    regions, geoms = [], []
    for cc, info in countries.items():
        geom = country_geometry(world, info["iso3"])
        if "mainland_bbox" in info:
            geom = restrict_to_bbox(geom, info["mainland_bbox"])
        regions.append(info["region"])
        geoms.append(geom)

    out = gpd.GeoDataFrame({"region": regions, "geometry": geoms}, crs=4326)
    out["geometry"] = out["geometry"].buffer(0)  # fix occasional invalid geometries

    out.to_parquet(OUT_PATH)
    print("Wrote:", OUT_PATH)
    print(out[["region"]])


if __name__ == "__main__":
    main()
