configfile: "config/config.yaml"

enabled_areas = [code for code, info in config["entsoe"]["areas"].items() if info.get("enabled")]
CF_COUNTRIES  = [code for code, info in config["res_cf"]["countries"].items() if info.get("enabled")]
CF_YEAR       = config["res_cf"]["year"]
CF_QUARTERS   = ["q1", "q2", "q3", "q4"]
CF_TOP_N      = config["res_cf"]["top_n"]


rule all:
    input:
        # Grid pipeline
        "resources/entsoe_processed.feather",
        "resources/nem_processed.feather",
        expand("resources/entsoe/{area}/{data_type}.feather",
               area=enabled_areas,
               data_type=config['entsoe']['data_types']),
        # res_cf pipeline
        expand("resources/res_cf/{country}_cf_{year}.csv",
               country=CF_COUNTRIES, year=CF_YEAR),
        expand("resources/res_cf/{country}_cf_{year}_bestsite_p95.csv",
               country=CF_COUNTRIES, year=CF_YEAR),
        f"resources/res_cf/resource_spread_{CF_YEAR}.csv",
        expand("resources/res_cf/{country}_complementarity_top{n}_{year}.csv",
               country=CF_COUNTRIES, n=CF_TOP_N, year=CF_YEAR),


# ── Grid pipeline ──────────────────────────────────────────────────────────────

rule download_entsoe:
    output:
        "resources/entsoe/{area}/{data_type}.feather"
    params:
        start_date=config["entsoe"]["start_date"],
        end_date=config["entsoe"]["end_date"],
        cache_dir="data/entsoe_cache",
    resources:
        entsoe_api=1
    script:
        "scripts/grid/download_entsoe.py"

rule download_nem:
    output:
        "resources/nem_processed.feather"
    params:
        start_date=config["nem_download"]["start_date"],
        end_date=config["nem_download"]["end_date"],
        cache_dir=config["nem_download"]["cache_dir"],
        rebuild=config["nem_download"]["rebuild"]
    script: "scripts/grid/download_nem.py"

rule process_entsoe:
    input:
        entsoe=expand("resources/entsoe/{area}/{data_type}.feather",
                      area=enabled_areas,
                      data_type=config['entsoe']['data_types'])
    params:
        areas=enabled_areas
    output:
        "resources/entsoe_processed.feather",
    script:
        "scripts/grid/process_entsoe.py"


# ── res_cf pipeline ─────────────────────────────────────────────────────────────

rule build_regions:
    input: "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
    output: "resources/shapes/regions.geojson"
    script: "scripts/res_cf/build_regions.py"

rule build_offshore_regions:
    input:
        regions="resources/shapes/regions.geojson",
        eez="data/shapes/eez/eez_v12.shp"
    output: "resources/shapes/offshore_regions.geojson"
    script: "scripts/res_cf/build_offshore_regions.py"

rule make_cutout:
    input:
        regions="resources/shapes/regions.geojson"
    output: "cutouts/{country}_{year}_{quarter}.nc"
    script: "scripts/res_cf/make_cutout.py"

rule build_cf_timeseries:
    input:
        cutout="cutouts/{country}_{year}_{quarter}.nc",
        regions="resources/shapes/regions.geojson"
    output:
        wind_onshore="resources/res_cf/quarterly/{country}_wind_onshore_{year}_{quarter}.csv",
        wind_offshore="resources/res_cf/quarterly/{country}_wind_offshore_{year}_{quarter}.csv",
        solar="resources/res_cf/quarterly/{country}_solar_{year}_{quarter}.csv"
    script: "scripts/res_cf/build_cf_timeseries.py"

rule concat_quarters:
    input:
        wind_onshore=expand("resources/res_cf/quarterly/{{country}}_wind_onshore_{{year}}_{quarter}.csv",
                            quarter=CF_QUARTERS),
        wind_offshore=expand("resources/res_cf/quarterly/{{country}}_wind_offshore_{{year}}_{quarter}.csv",
                             quarter=CF_QUARTERS),
        solar=expand("resources/res_cf/quarterly/{{country}}_solar_{{year}}_{quarter}.csv",
                     quarter=CF_QUARTERS),
    output:
        wind_onshore="resources/res_cf/annual/{country}_wind_onshore_{year}.csv",
        wind_offshore="resources/res_cf/annual/{country}_wind_offshore_{year}.csv",
        solar="resources/res_cf/annual/{country}_solar_{year}.csv"
    script: "scripts/res_cf/concat_quarters.py"

rule combine_techs:
    input:
        wind_onshore="resources/res_cf/annual/{country}_wind_onshore_{year}.csv",
        wind_offshore="resources/res_cf/annual/{country}_wind_offshore_{year}.csv",
        solar="resources/res_cf/annual/{country}_solar_{year}.csv"
    output: "resources/res_cf/{country}_cf_{year}.csv"
    script: "scripts/res_cf/combine_techs.py"

rule resource_spread:
    input:
        cutouts=expand("cutouts/{country}_{year}_{quarter}.nc",
                       country=CF_COUNTRIES, year=CF_YEAR, quarter=CF_QUARTERS),
        national_cfs=expand("resources/res_cf/{country}_cf_{year}.csv",
                            country=CF_COUNTRIES, year=CF_YEAR),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output: f"resources/res_cf/resource_spread_{CF_YEAR}.csv"
    script: "scripts/res_cf/resource_spread.py"

rule make_bestsite_cf:
    input:
        cutouts=expand("cutouts/{{country}}_{year}_{quarter}.nc",
                       year=CF_YEAR, quarter=CF_QUARTERS),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output: f"resources/res_cf/{{country}}_cf_{CF_YEAR}_bestsite_p95.csv"
    script: "scripts/res_cf/make_bestsite_cf.py"

rule complementarity:
    input:
        national_cf=f"resources/res_cf/{{country}}_cf_{CF_YEAR}.csv",
        cutouts=expand("cutouts/{{country}}_{year}_{quarter}.nc",
                       year=CF_YEAR, quarter=CF_QUARTERS),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output:
        top=f"resources/res_cf/{{country}}_complementarity_top{CF_TOP_N}_{CF_YEAR}.csv",
        avg=f"resources/res_cf/{{country}}_average_profiles_{CF_YEAR}.csv"
    script: "scripts/res_cf/complementarity.py"
