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

rule all:
    input:
        "data/processed_data.feather",
        "data/nem_processed.feather",
        expand("data/entsoe_cache/{area}/{year_month}/{data_type}.feather",
               area=enabled_areas,
               year_month=DOWNLOAD_YEAR_MONTHS,
               data_type=config['data_types'])

rule download_entsoe_data:
    output:
        temp("data/entsoe_cache/{area}/{year_month}/{data_type}.feather")
    input:
        # If rebuild is false, the rule depends on its own output.
        # If the output exists, Snakemake skips this job.
        lambda wildcards: [] if config["entsoe_download"]["rebuild"] else "data/entsoe_cache/{area}/{year_month}/{data_type}.feather"
    conda: "environment.yaml"
    script:
        "scripts/download_entsoe_data.py"

rule download_nem_data:
    output:
        "data/nem_processed.feather"
    params:
        start_date=config["nem_download"]["start_date"],
        end_date=config["nem_download"]["end_date"],
        nemosis_cache_dir=config["nem_download"]["nemosis_cache_dir"],
        rebuild=config["nem_download"]["rebuild"]
    conda: "environment.yaml"
    script: "scripts/download_nem_data.py"

rule integrate_data:
    input:
        expand("data/entsoe_cache/{area}/{year_month}/{data_type}.feather",
               area=enabled_areas,
               year_month=INTEGRATE_YEAR_MONTHS,
               data_type=config['data_types'])
    output:
        expand("data/integrated/{data_type}.feather",
               data_type=config['data_types'])
    conda: "environment.yaml"
    script:
        "scripts/integrate_entsoe_data.py"

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
    conda: "environment.yaml"
    script:
        "scripts/process_entsoe_data.py"