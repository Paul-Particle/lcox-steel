from pathlib import Path
import atlite
import geopandas as gpd
import pandas as pd

REGIONS_PATH = "data/shapes/regions.geojson"
OUTDIR = Path("data/res_cf")

# --- set these two for each run ---
CUTOUT_PATH = "data/cutouts/bra_2023_q1.nc"   # or q2
TAG = "q1"                                   # or q2
YEAR = 2023
COUNTRY = "bra"
# ----------------------------------


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


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    regions = gpd.read_file(REGIONS_PATH).to_crs(4326)
    gdf = regions.loc[regions["region"] == COUNTRY.upper(), ["region", "geometry"]].copy()

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

    wind_cf.to_csv(wind_out)   # keeps time index in first column
    solar_cf.to_csv(solar_out)

    print("Wrote:")
    print(" -", wind_out)
    print(" -", solar_out)


if __name__ == "__main__":
    main()
