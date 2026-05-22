import yaml
import pandas as pd

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

def get_year_months(start_date_str, end_date_str):
    start = pd.to_datetime(start_date_str, format='%Y%m%d')
    end = pd.to_datetime(end_date_str, format='%Y%m%d')
    date_range = pd.date_range(start, end, freq='MS')
    return sorted(list(set(date_range.strftime('%Y-%m'))))

DOWNLOAD_YEAR_MONTHS = get_year_months(config["entsoe_download"]["start_date"], config["entsoe_download"]["end_date"])
INTEGRATE_YEAR_MONTHS = get_year_months(config["entsoe_integration"]["start_date"], config["entsoe_integration"]["end_date"])

def get_enabled_areas():
    areas_df = pd.read_csv("areas.csv")
    return areas_df[areas_df["enabled"]]["area_code"].tolist()

enabled_areas = get_enabled_areas()

CF_COUNTRIES = config["res_cf"]["countries"]
CF_YEAR      = config["res_cf"]["year"]
CF_QUARTERS  = ["q1", "q2", "q3", "q4"]
CF_TOP_N     = config["res_cf"]["top_n"]


rule all:
    input:
        # Grid pipeline
        "data/processed_data.feather",
        "data/nem_processed.feather",
        expand("data/entsoe_cache/{area}/{year_month}/{data_type}.feather",
               area=enabled_areas,
               year_month=DOWNLOAD_YEAR_MONTHS,
               data_type=config['data_types']),
        # res_cf pipeline
        expand("data/res_cf/{country}_cf_{year}.csv",
               country=CF_COUNTRIES, year=CF_YEAR),
        expand("data/res_cf/{country}_cf_{year}_bestsite_p95.csv",
               country=CF_COUNTRIES, year=CF_YEAR),
        f"data/res_cf/resource_spread_{CF_YEAR}.csv",
        expand("data/res_cf/{country}_complementarity_top{n}_{year}.csv",
               country=CF_COUNTRIES, n=CF_TOP_N, year=CF_YEAR),


# ── Grid pipeline ──────────────────────────────────────────────────────────────

rule download_entsoe_data:
    output:
        "data/entsoe_cache/{area}/{year_month}/{data_type}.feather"
    script:
        "scripts/grid/download_entsoe.py"

rule download_nem_data:
    output:
        "data/nem_processed.feather"
    params:
        start_date=config["nem_download"]["start_date"],
        end_date=config["nem_download"]["end_date"],
        nemosis_cache_dir=config["nem_download"]["nemosis_cache_dir"],
        rebuild=config["nem_download"]["rebuild"]
    script: "scripts/grid/download_nem.py"

rule integrate_data:
    input:
        expand("data/entsoe_cache/{area}/{year_month}/{data_type}.feather",
               area=enabled_areas,
               year_month=INTEGRATE_YEAR_MONTHS,
               data_type=config['data_types'])
    output:
        expand("data/integrated/{data_type}.feather",
               data_type=config['data_types'])
    script:
        "scripts/grid/integrate_entsoe.py"

rule process_data:
    input:
        prices="data/integrated/prices.feather",
        load_forecast="data/integrated/load_forecast.feather",
        load_actual="data/integrated/load_actual.feather",
        vre="data/integrated/vre.feather",
        generation="data/integrated/generation.feather",
        crossborder="data/integrated/crossborder.feather",
        areas_config="areas.csv"
    output:
        "data/processed_data.feather",
    script:
        "scripts/grid/process_entsoe.py"


# ── res_cf pipeline ─────────────────────────────────────────────────────────────

rule build_regions:
    input: "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
    output: "data/shapes/regions.geojson"
    script: "scripts/res_cf/build_regions.py"

rule build_offshore_regions:
    input:
        regions="data/shapes/regions.geojson",
        eez="data/shapes/eez/eez_v12.shp"
    output: "data/shapes/offshore_regions.geojson"
    script: "scripts/res_cf/build_offshore_regions.py"

rule make_cutout:
    input:
        regions="data/shapes/regions.geojson",
        offshore_regions="data/shapes/offshore_regions.geojson"
    output: "data/cutouts/{country}_{year}_{quarter}.nc"
    script: "scripts/res_cf/make_cutouts.py"

rule build_cf_timeseries:
    input:
        cutout="data/cutouts/{country}_{year}_{quarter}.nc",
        regions="data/shapes/regions.geojson",
        offshore_regions="data/shapes/offshore_regions.geojson"
    output:
        wind_onshore="data/res_cf/{country}_wind_onshore_cf_{year}_{quarter}.csv",
        wind_offshore="data/res_cf/{country}_wind_offshore_cf_{year}_{quarter}.csv",
        solar="data/res_cf/{country}_solar_cf_{year}_{quarter}.csv"
    script: "scripts/res_cf/build_cf_timeseries.py"

rule concat_quarters:
    input:
        wind_onshore=expand("data/res_cf/{{country}}_wind_onshore_cf_{{year}}_{quarter}.csv",
                            quarter=CF_QUARTERS),
        wind_offshore=expand("data/res_cf/{{country}}_wind_offshore_cf_{{year}}_{quarter}.csv",
                             quarter=CF_QUARTERS),
        solar=expand("data/res_cf/{{country}}_solar_cf_{{year}}_{quarter}.csv",
                     quarter=CF_QUARTERS),
    output:
        wind_onshore="data/res_cf/{country}_wind_onshore_cf_{year}.csv",
        wind_offshore="data/res_cf/{country}_wind_offshore_cf_{year}.csv",
        solar="data/res_cf/{country}_solar_cf_{year}.csv"
    script: "scripts/res_cf/concat_quarters.py"

rule combine_techs:
    input:
        wind_onshore="data/res_cf/{country}_wind_onshore_cf_{year}.csv",
        wind_offshore="data/res_cf/{country}_wind_offshore_cf_{year}.csv",
        solar="data/res_cf/{country}_solar_cf_{year}.csv"
    output: "data/res_cf/{country}_cf_{year}.csv"
    script: "scripts/res_cf/combine_techs.py"

rule resource_spread:
    input:
        cutouts=expand("data/cutouts/{country}_{year}_{quarter}.nc",
                       country=CF_COUNTRIES, year=CF_YEAR, quarter=CF_QUARTERS),
        national_cfs=expand("data/res_cf/{country}_cf_{year}.csv",
                            country=CF_COUNTRIES, year=CF_YEAR),
        regions="data/shapes/regions.geojson",
        offshore_regions="data/shapes/offshore_regions.geojson"
    output: f"data/res_cf/resource_spread_{CF_YEAR}.csv"
    script: "scripts/res_cf/resource_spread.py"

rule make_bestsite_cf:
    input:
        cutouts=expand("data/cutouts/{{country}}_{year}_{quarter}.nc",
                       year=CF_YEAR, quarter=CF_QUARTERS),
        regions="data/shapes/regions.geojson",
        offshore_regions="data/shapes/offshore_regions.geojson"
    output: f"data/res_cf/{{country}}_cf_{CF_YEAR}_bestsite_p95.csv"
    script: "scripts/res_cf/make_bestsite_cf_timeseries.py"

rule complementarity:
    input:
        national_cf=f"data/res_cf/{{country}}_cf_{CF_YEAR}.csv",
        cutouts=expand("data/cutouts/{{country}}_{year}_{quarter}.nc",
                       year=CF_YEAR, quarter=CF_QUARTERS),
        regions="data/shapes/regions.geojson",
        offshore_regions="data/shapes/offshore_regions.geojson"
    output:
        top=f"data/res_cf/{{country}}_complementarity_top{CF_TOP_N}_{CF_YEAR}.csv",
        avg=f"data/res_cf/{{country}}_average_profiles_{CF_YEAR}.csv"
    script: "scripts/res_cf/complementarity.py"
