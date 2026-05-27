from pathlib import Path

import atlite
import geopandas as gpd
import pandas as pd

from common._paths import CUTOUTS, RES_CF, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

# Standalone defaults
_CF_AREA = "de"
_TECH = "wind_onshore"
_START_DATE = "20230101"
_END_DATE = "20231231"
_CUTOUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_REGIONS_PATH = SHAPES_RES / "de_geo.parquet"
_OFFSHORE_REGIONS_PATH = SHAPES_RES / "de_offshore_geo.parquet"
_OUT = RES_CF / "de_wind_onshore_country-average_20230101_20231231.parquet"
_REGION = "DE"
_WIND_ONSHORE_TURBINE = "Vestas_V112_3MW"
_WIND_OFFSHORE_TURBINE = "NREL_ReferenceTurbine_5MW_offshore"
_PV_PANEL = "CSi"
_PV_ORIENTATION = "latitude_optimal"
_WIND_CF = {"smooth": True, "add_cutout_windspeed": True}

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _CF_AREA = snakemake.wildcards.cf_area
    _TECH = snakemake.wildcards.tech
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _CUTOUT_PATH = Path(snakemake.input.cutout)
    _REGIONS_PATH = Path(snakemake.input.regions)
    _OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    _OUT = Path(snakemake.output[0])
    _REGION = snakemake.params.region
    _WIND_ONSHORE_TURBINE = snakemake.params.wind_onshore_turbine
    _WIND_OFFSHORE_TURBINE = snakemake.params.wind_offshore_turbine
    _PV_PANEL = snakemake.params.pv_panel
    _PV_ORIENTATION = snakemake.params.pv_orientation
    _WIND_CF = snakemake.params.wind_cf


_WIND_SMOOTH = _WIND_CF.get("smooth", True)
_WIND_ADD_CUTOUT_WS = _WIND_CF.get("add_cutout_windspeed", True)


def to_cf_series(x, name: str = "cf") -> pd.Series:
    obj = x.to_pandas()
    if not isinstance(obj, pd.DataFrame):
        s = obj
    elif obj.shape[1] == 1:
        s = obj.iloc[:, 0]
    elif _REGION in obj.columns:
        s = obj[_REGION]
    else:
        s = obj.iloc[:, 0]
    s.index = pd.to_datetime(s.index)
    s.index.name = "time"
    return s.rename(name).clip(0, 1)


def get_region_gdf(path: Path) -> gpd.GeoDataFrame:
    gdf = gpd.read_parquet(path).to_crs(4326)
    gdf = gdf.loc[gdf["region"] == _REGION, ["region", "geometry"]].copy()
    if gdf.empty:
        raise ValueError(f"Region '{_REGION}' not found in {path}")
    return gdf


def main() -> None:
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    cutout = atlite.Cutout(str(_CUTOUT_PATH))

    if _TECH == "solar":
        gdf = get_region_gdf(_REGIONS_PATH)
        matrix = cutout.indicatormatrix(gdf)
        cf = cutout.pv(
            matrix=matrix,
            panel=_PV_PANEL,
            orientation=_PV_ORIENTATION,
            capacity_factor=False,
            per_unit=True,
        )
        series = to_cf_series(cf)

    elif _TECH == "wind_onshore":
        gdf = get_region_gdf(_REGIONS_PATH)
        matrix = cutout.indicatormatrix(gdf)
        cf = cutout.wind(
            matrix=matrix,
            turbine=_WIND_ONSHORE_TURBINE,
            capacity_factor=False,
            per_unit=True,
            smooth=_WIND_SMOOTH,
            add_cutout_windspeed=_WIND_ADD_CUTOUT_WS,
        )
        series = to_cf_series(cf)

    elif _TECH == "wind_offshore":
        gdf = get_region_gdf(_OFFSHORE_REGIONS_PATH)
        matrix = cutout.indicatormatrix(gdf)
        cf = cutout.wind(
            matrix=matrix,
            turbine=_WIND_OFFSHORE_TURBINE,
            capacity_factor=False,
            per_unit=True,
            smooth=_WIND_SMOOTH,
            add_cutout_windspeed=_WIND_ADD_CUTOUT_WS,
        )
        series = to_cf_series(cf)

    else:
        raise ValueError(f"Unknown tech: {_TECH!r}")

    series.to_frame().to_parquet(_OUT, index=True)
    print(f"Wrote: {_OUT}  ({len(series)} rows, mean={series.mean():.3f})")


if __name__ == "__main__":
    main()
