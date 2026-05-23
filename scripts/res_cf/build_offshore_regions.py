from pathlib import Path
import geopandas as gpd
import yaml


EEZ_SHP      = Path("data/shapes/eez/eez_v12.shp")
OUT_GEOJSON  = Path("resources/shapes/offshore_regions.geojson")
LAND_REGIONS = Path("resources/shapes/regions.geojson")


def _read_offshore_max_distance_km() -> float:
    with open("config/config.yaml", "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return float(cfg["res_cf"]["offshore_max_distance_km"])


if "snakemake" in dir():
    OUT_GEOJSON  = Path(snakemake.output[0])
    LAND_REGIONS = Path(snakemake.input.regions)
    _OFFSHORE_MAX_KM = float(snakemake.config["res_cf"]["offshore_max_distance_km"])
else:
    _OFFSHORE_MAX_KM = None  # resolved in main() via _read_offshore_max_distance_km()


ISO3_TO_REGION = {
    "DEU": "DE",
    "FRA": "FR",
    "ESP": "ES",
    "AUS": "AUS",
    "BRA": "BRA",
}


def main():
    max_distance_km = _OFFSHORE_MAX_KM if _OFFSHORE_MAX_KM is not None else _read_offshore_max_distance_km()
    print(f"Configured offshore distance limit: {max_distance_km} km")

    if not EEZ_SHP.exists():
        raise FileNotFoundError(f"Cannot find EEZ shapefile at: {EEZ_SHP}")
    if not LAND_REGIONS.exists():
        raise FileNotFoundError(f"Cannot find land regions file at: {LAND_REGIONS}")

    OUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)

    eez = gpd.read_file(EEZ_SHP).to_crs(4326)

    gdf = eez.loc[
        (eez["ISO_TER1"].isin(ISO3_TO_REGION.keys())) &
        (eez["POL_TYPE"] == "200NM"),
        ["ISO_TER1", "geometry"]
    ].copy()

    gdf["region"] = gdf["ISO_TER1"].map(ISO3_TO_REGION)
    gdf = gdf[["region", "geometry"]]
    gdf = gdf.dissolve(by="region", as_index=False)

    land = gpd.read_file(LAND_REGIONS).to_crs(4326)

    land_buffer = land.to_crs(6933).copy()
    land_buffer["geometry"] = land_buffer.buffer(max_distance_km * 1000)
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
