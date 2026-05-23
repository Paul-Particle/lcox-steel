"""
Build offshore region geometries from EEZ shapefile + onshore regions.

For each enabled country in config.yaml, takes the country's EEZ polygon, subtracts
the onshore land, and clips to a buffer of `res_cf.offshore_max_distance_km` km
around the coast. Result: roughly the part of the EEZ within commuting distance
of shore — the area where offshore wind makes economic sense.

The iso3 → region mapping is read from config["res_cf"]["countries"][cc].
"""

from pathlib import Path

import geopandas as gpd
import yaml

from common._paths import DATA, SHAPES_RES, REPO_ROOT

if "snakemake" not in globals():
    from common._stubs import snakemake

EEZ_SHAPE      = DATA / "shapes/eez/eez_v12.shp"
LAND_REGIONS   = SHAPES_RES / "regions.geojson"
OUT_GEOJSON    = SHAPES_RES / "offshore_regions.geojson"


def load_cfg() -> dict:
    with open(REPO_ROOT / "config/config.yaml") as f:
        return yaml.safe_load(f)


if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    OUT_GEOJSON     = Path(snakemake.output[0])
    LAND_REGIONS    = Path(snakemake.input.regions)
    OFFSHORE_MAX_KM = float(snakemake.config["res_cf"]["offshore_max_distance_km"])
    COUNTRIES_CFG   = snakemake.config["res_cf"]["countries"]
else:
    _cfg            = load_cfg()
    OFFSHORE_MAX_KM = float(_cfg["res_cf"]["offshore_max_distance_km"])
    COUNTRIES_CFG   = _cfg["res_cf"]["countries"]


def main() -> None:
    iso3_to_region = {info["iso3"]: info["region"] for info in COUNTRIES_CFG.values() if info.get("enabled")}
    print(f"Configured offshore distance limit: {OFFSHORE_MAX_KM} km")
    print(f"Building offshore regions for: {sorted(iso3_to_region.values())}")

    if not EEZ_SHAPE.exists():
        raise FileNotFoundError(f"Cannot find EEZ shapefile at: {EEZ_SHAPE}")
    if not LAND_REGIONS.exists():
        raise FileNotFoundError(f"Cannot find land regions file at: {LAND_REGIONS}")

    OUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)

    eez = gpd.read_file(EEZ_SHAPE).to_crs(4326)

    gdf = eez.loc[
        (eez["ISO_TER1"].isin(iso3_to_region.keys())) &
        (eez["POL_TYPE"] == "200NM"),
        ["ISO_TER1", "geometry"]
    ].copy()

    gdf["region"] = gdf["ISO_TER1"].map(iso3_to_region)
    gdf = gdf[["region", "geometry"]].dissolve(by="region", as_index=False)

    land = gpd.read_file(LAND_REGIONS).to_crs(4326)

    land_buffer = land.to_crs(6933).copy()
    land_buffer["geometry"] = land_buffer.buffer(OFFSHORE_MAX_KM * 1000)
    land_buffer = land_buffer.to_crs(4326)
    land_buffer = land_buffer[["region", "geometry"]].rename(columns={"geometry": "geometry_buffer"})

    gdf = gdf.merge(land, on="region", suffixes=("_eez", "_land"))
    gdf = gdf.merge(land_buffer, on="region", suffixes=("", "_buffer"))

    gdf["geometry"] = gdf.apply(
        lambda r: r.geometry_eez.difference(r.geometry_land).intersection(r.geometry_buffer),
        axis=1
    )

    gdf = gdf[["region", "geometry"]]
    gdf["geometry"] = gdf["geometry"].buffer(0)
    if not gdf.is_valid.all():
        print("Warning: some offshore geometries are invalid")

    gdf = gdf.set_crs(4326)
    gdf.to_file(OUT_GEOJSON, driver="GeoJSON")

    gdf_area = gdf.to_crs(6933).copy()
    gdf_area["area_km2"] = gdf_area.geometry.area / 1e6

    print("Wrote:", OUT_GEOJSON)
    print(gdf_area[["region", "area_km2"]].sort_values("region").to_string(index=False))


if __name__ == "__main__":
    main()
