from pathlib import Path
import atlite
import geopandas as gpd
import pandas as pd
from scipy.sparse import diags

import logging
from common._paths import CUTOUTS, SHAPES_RES, TIMESERIES
from scripts.res_cf._helpers_res_cf import cos_lat_weights, eligibility_matrix

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

REGIONS_PATH = SHAPES_RES / "de_geo.parquet" # was "data/shapes/regions.geojson"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "de_offshore_geo.parquet"
OUTDIR = TIMESERIES / "DE_LU_wind-onshore_area-average_20230101_20231231.parquet"

# --- set these for each standalone run ---
CUTOUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
COUNTRY = "DE"
# ----------------------------------

_WIND_ONSHORE_TURBINE = "Vestas_V112_3MW"
_WIND_OFFSHORE_TURBINE = "NREL_ReferenceTurbine_5MW_offshore"
_PV_PANEL = "CSi"
_PV_ORIENTATION = "latitude_optimal"
_WIND_CF = {"smooth": True, "add_cutout_windspeed": True}
_MIN_LAND_FRACTION = 0.98
_ELIGIBILITY_SOURCE = "indicatormatrix"
_TECH = "wind-onshore"
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _TECH = snakemake.wildcards.tech
    CUTOUT_PATH = Path(snakemake.input.cutout)
    REGIONS_PATH = Path(snakemake.input.regions)
    OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    OUTDIR = Path(snakemake.output[0])
    COUNTRY = snakemake.params.region
    _WIND_ONSHORE_TURBINE = snakemake.params.wind_onshore_turbine
    _WIND_OFFSHORE_TURBINE = snakemake.params.wind_offshore_turbine
    _PV_PANEL = snakemake.params.pv_panel
    _PV_ORIENTATION = snakemake.params.pv_orientation
    _WIND_CF = snakemake.params.wind_cf
    _MIN_LAND_FRACTION = snakemake.params.min_land_fraction
    _ELIGIBILITY_SOURCE = snakemake.params.eligibility_source
_WIND_SMOOTH = _WIND_CF.get("smooth", True)
_WIND_ADD_CUTOUT_WS = _WIND_CF.get("add_cutout_windspeed", True)


def to_cf_series(x, name="cf"):
    obj = x.to_pandas()

    if not isinstance(obj, pd.DataFrame):
        s = obj
    elif obj.shape[1] == 1:
        s = obj.iloc[:, 0]
    elif COUNTRY in obj.columns:
        s = obj[COUNTRY]
    else:
        s = obj.iloc[:, 0]

    s.index = pd.to_datetime(s.index)
    s.index.name = "time"
    return s.rename(name).clip(0, 1)

def get_region_gdf(path: Path) -> gpd.GeoDataFrame:
    gdf = gpd.read_parquet(path).to_crs(4326)
    gdf = gdf.loc[gdf["region"] == COUNTRY, ["region", "geometry"]].copy()

    if gdf.empty:
        raise ValueError(f"Region '{COUNTRY}' not found in {path}")

    return gdf

def main():
    OUTDIR.parent.mkdir(parents=True, exist_ok=True)

    gdf = get_region_gdf(REGIONS_PATH)
    offshore_gdf = get_region_gdf(OFFSHORE_REGIONS_PATH)

    cutout = atlite.Cutout(str(CUTOUT_PATH))
    # Onshore land weights carry the #41 land-sea eligibility cutoff (drops
    # sea-contaminated coastal border cells) via the configured source, and are
    # cos(lat) physical-area weighted (#37) — both folded into eligibility_matrix.
    # Offshore is sited on sea by design, so it keeps the raw indicatormatrix (no
    # cutoff, no source swap) but is still area-weighted by cos(lat) for #37.
    onshore_geom = gdf.geometry.iloc[0]
    matrix = eligibility_matrix(cutout, onshore_geom, _MIN_LAND_FRACTION, _ELIGIBILITY_SOURCE)
    offshore_matrix = cutout.indicatormatrix(offshore_gdf).tocsr() @ diags(cos_lat_weights(cutout))

    if _TECH == "wind-onshore":
        wind_cf = cutout.wind(
            matrix=matrix,
            turbine=_WIND_ONSHORE_TURBINE,
            capacity_factor=False,
            per_unit=True,
            smooth=_WIND_SMOOTH,
            add_cutout_windspeed=_WIND_ADD_CUTOUT_WS,
        )
        series = to_cf_series(wind_cf)

    elif _TECH == "wind-offshore":
        offshore_wind_cf = cutout.wind(
            matrix=offshore_matrix,
            turbine=_WIND_OFFSHORE_TURBINE,
            capacity_factor=False,
            per_unit=True,
            smooth=_WIND_SMOOTH,
            add_cutout_windspeed=_WIND_ADD_CUTOUT_WS,
        )
        series = to_cf_series(offshore_wind_cf)

    elif _TECH == "solar":
        solar_cf = cutout.pv(
            matrix=matrix,
            panel=_PV_PANEL,
            orientation=_PV_ORIENTATION,
            capacity_factor=False,
            per_unit=True,
        )
        series = to_cf_series(solar_cf)

    else:
        raise ValueError(f"Unknown tech: {_TECH!r}")

    # Column name = tech wildcard, so downstream scripts (solve_network)
    # can read the tech key straight off the parquet without a separate param.
    series.rename(_TECH).to_frame().to_parquet(OUTDIR, index=True)

    log.info("Wrote:")
    log.info(f"{OUTDIR}")
    log.info(f"({len(series)} rows")
    log.info(f"mean={series.mean():.3f})")


if __name__ == "__main__":
    main()
