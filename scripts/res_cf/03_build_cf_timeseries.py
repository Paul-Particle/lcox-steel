from pathlib import Path
import atlite
import geopandas as gpd
import pandas as pd

REGIONS_PATH = "data/shapes/regions.geojson"
OFFSHORE_REGIONS_PATH = "data/shapes/offshore_regions.geojson"
OUTDIR = Path("data/res_cf")

# --- set these for each run ---
CUTOUT_PATH = "data/cutouts/de_2023_jan2w.nc"
TAG = "jan2w"
YEAR = 2023
COUNTRY = "de"
# ------------------------------


def to_cf_series(x, name="cf"):
    obj = x.to_pandas()

    if isinstance(obj, pd.DataFrame):
        if obj.shape[1] == 1:
            s = obj.iloc[:, 0]
        else:
            col = COUNTRY.upper() if COUNTRY.upper() in obj.columns else obj.columns[0]
            s = obj[col]
    else:
        s = obj

    s.index = pd.to_datetime(s.index)
    s.index.name = "time"
    return s.rename(name).clip(0, 1)

def get_region_gdf(path: str, region_code: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path).to_crs(4326)
    gdf = gdf.loc[gdf["region"] == region_code.upper(), ["region", "geometry"]].copy()

    if gdf.empty:
        raise ValueError(f"Region '{region_code.upper()}' not found in {path}")

    return gdf

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    gdf = get_region_gdf(REGIONS_PATH, COUNTRY)
    has_offshore = Path(OFFSHORE_REGIONS_PATH).exists()

    cutout = atlite.Cutout(CUTOUT_PATH)
    matrix = cutout.indicatormatrix(gdf)

    wind_cf = cutout.wind(
        matrix=matrix,
        turbine="Vestas_V112_3MW",
        capacity_factor=False,
        per_unit=True,
    )

    solar_cf = cutout.pv(
        matrix=matrix,
        panel="CSi",
        orientation="latitude_optimal",
        capacity_factor=False,
        per_unit=True,
    )

    wind_cf = to_cf_series(wind_cf)
    solar_cf = to_cf_series(solar_cf)

    wind_out = OUTDIR / f"{COUNTRY}_wind_onshore_cf_{YEAR}_{TAG}.csv"
    solar_out = OUTDIR / f"{COUNTRY}_solar_cf_{YEAR}_{TAG}.csv"

    wind_cf.to_csv(wind_out)
    solar_cf.to_csv(solar_out)

    print("Wrote:")
    print(" -", wind_out)
    print(" -", solar_out)

    if has_offshore:
        offshore_gdf = get_region_gdf(OFFSHORE_REGIONS_PATH, COUNTRY)
        offshore_matrix = cutout.indicatormatrix(offshore_gdf)
        offshore_wind_cf = cutout.wind(
            matrix=offshore_matrix,
            turbine="NREL_ReferenceTurbine_5MW_offshore",
            capacity_factor=False,
            per_unit=True,
        )
        offshore_wind_cf = to_cf_series(offshore_wind_cf)
        offshore_wind_out = OUTDIR / f"{COUNTRY}_wind_offshore_cf_{YEAR}_{TAG}.csv"
        offshore_wind_cf.to_csv(offshore_wind_out)
        print(" -", offshore_wind_out)
    else:
        print("No offshore regions file — skipping offshore wind.")


if __name__ == "__main__":
    main()
