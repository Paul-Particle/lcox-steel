from pathlib import Path

import atlite
import geopandas as gpd
import pandas as pd

from common._paths import CUTOUTS, RES_CF, SHAPES_RES
from scripts.res_cf._helpers import load_res_cf_cfg

if "snakemake" not in globals():
    from common._stubs import snakemake

REGIONS_PATH          = SHAPES_RES / "regions.parquet"
OFFSHORE_REGIONS_PATH = SHAPES_RES / "offshore_regions.parquet"


# --- defaults for standalone use ---
CF_AREA       = "de"
START_DATE    = "20230101"
END_DATE      = "20231231"
CUTOUT_PATH   = CUTOUTS / f"{CF_AREA}_{START_DATE}_{END_DATE}.nc"
WIND_ON_OUT   = RES_CF / f"{CF_AREA}_wind_onshore_country-average_{START_DATE}_{END_DATE}.parquet"
WIND_OFF_OUT  = RES_CF / f"{CF_AREA}_wind_offshore_country-average_{START_DATE}_{END_DATE}.parquet"
SOLAR_OUT     = RES_CF / f"{CF_AREA}_solar_country-average_{START_DATE}_{END_DATE}.parquet"
RES_CF_CFG    = load_res_cf_cfg()

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    CF_AREA               = snakemake.wildcards.cf_area.lower()
    START_DATE            = snakemake.wildcards.start_date
    END_DATE              = snakemake.wildcards.end_date
    CUTOUT_PATH           = snakemake.input.cutout
    OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    WIND_ON_OUT           = Path(snakemake.output.wind_onshore)
    WIND_OFF_OUT          = Path(snakemake.output.wind_offshore)
    SOLAR_OUT             = Path(snakemake.output.solar)
    RES_CF_CFG            = snakemake.config["res_cf"]

REGION_TAG            = RES_CF_CFG["countries"][CF_AREA]["region"]
WIND_ONSHORE_TURBINE  = RES_CF_CFG["wind_onshore_turbine"]
WIND_OFFSHORE_TURBINE = RES_CF_CFG["wind_offshore_turbine"]
PV_PANEL              = RES_CF_CFG["pv_panel"]
PV_ORIENTATION        = RES_CF_CFG["pv_orientation"]
WIND_CF_CFG           = RES_CF_CFG.get("wind_cf", {})
WIND_SMOOTH           = WIND_CF_CFG.get("smooth", True)
WIND_ADD_CUTOUT_WS    = WIND_CF_CFG.get("add_cutout_windspeed", True)


def to_cf_series(x, name: str = "cf") -> pd.Series:
    obj = x.to_pandas()

    # We need a series but might get a df
    if not isinstance(obj, pd.DataFrame):
        s = obj
    elif obj.shape[1] == 1:
        s = obj.iloc[:, 0]
    else:
        if REGION_TAG in obj.columns:
            s = obj[REGION_TAG]
        else:
            s = obj.iloc[:, 0]

    s.index = pd.to_datetime(s.index)
    s.index.name = "time"

    return s.rename(name).clip(0, 1)


def get_region_gdf(path: Path, region_tag: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_parquet(path).to_crs(4326)
    gdf = gdf.loc[gdf["region"] == region_tag, ["region", "geometry"]].copy()
    if gdf.empty:
        raise ValueError(f"Region '{region_tag}' not found in {path}")
    return gdf


def main() -> None:
    WIND_ON_OUT.parent.mkdir(parents=True, exist_ok=True)

    gdf = get_region_gdf(REGIONS_PATH, REGION_TAG)
    has_offshore = OFFSHORE_REGIONS_PATH.exists()

    cutout = atlite.Cutout(str(CUTOUT_PATH))
    matrix = cutout.indicatormatrix(gdf)

    wind_cf = cutout.wind(
        matrix=matrix,
        turbine=WIND_ONSHORE_TURBINE,
        capacity_factor=False,
        per_unit=True,
        smooth=WIND_SMOOTH,
        add_cutout_windspeed=WIND_ADD_CUTOUT_WS,
    )

    solar_cf = cutout.pv(
        matrix=matrix,
        panel=PV_PANEL,
        orientation=PV_ORIENTATION,
        capacity_factor=False,
        per_unit=True,
    )

    wind_cf  = to_cf_series(wind_cf)
    solar_cf = to_cf_series(solar_cf)

    wind_cf.to_frame().to_parquet(WIND_ON_OUT, index=True)
    solar_cf.to_frame().to_parquet(SOLAR_OUT, index=True)

    print("Wrote:")
    print(" -", WIND_ON_OUT)
    print(" -", SOLAR_OUT)

    if has_offshore:
        offshore_gdf = get_region_gdf(OFFSHORE_REGIONS_PATH, REGION_TAG)
        offshore_matrix = cutout.indicatormatrix(offshore_gdf)
        offshore_wind_cf = cutout.wind(
            matrix=offshore_matrix,
            turbine=WIND_OFFSHORE_TURBINE,
            capacity_factor=False,
            per_unit=True,
            smooth=WIND_SMOOTH,
            add_cutout_windspeed=WIND_ADD_CUTOUT_WS,
        )
        offshore_wind_cf = to_cf_series(offshore_wind_cf)
    else:
        # Write zeros so Snakemake output contract is satisfied; real values require EEZ file.
        offshore_wind_cf = pd.Series(0.0, index=wind_cf.index, name="cf")
        offshore_wind_cf.index.name = "time"
        print("No offshore regions file — writing zero placeholder for offshore wind.")
    offshore_wind_cf.to_frame().to_parquet(WIND_OFF_OUT, index=True)
    print(" -", WIND_OFF_OUT)


if __name__ == "__main__":
    main()
