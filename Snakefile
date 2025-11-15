import yaml
import pandas as pd

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

def get_year_months(start_date_str, end_date_str):
    start = pd.to_datetime(start_date_str, format='%Y%m%d')
    end = pd.to_datetime(end_date_str, format='%Y%m%d')
    date_range = pd.date_range(start, end, freq='MS')
    return sorted(list(set(date_range.strftime('%Y-%m'))))

DOWNLOAD_YEAR_MONTHS = get_year_months(config["split_download"]["start_date"], config["split_download"]["end_date"])
INTEGRATE_YEAR_MONTHS = get_year_months(config["integration"]["start_date"], config["integration"]["end_date"])

def get_enabled_areas():
    areas_df = pd.read_csv("areas.csv")
    return areas_df[areas_df["enabled"]]["area_code"].tolist()

enabled_areas = get_enabled_areas()

rule all:
    input:
        "results/data_availability.html",

# Rule to download data for a specific area, year, and month
rule download_split:
    output:
        "data/downloads/{area}/{year}-{month}/{data_type}.feather"
    log:
        "logs/download/{area}-{year}-{month}-{data_type}.log"
    shell:
        "(python -u scripts/download_entsoe_data_split.py {wildcards.data_type} data/downloads {wildcards.area} {wildcards.year} {wildcards.month})" #&> {log}"

# Rule to integrate the downloaded data
rule integrate_data:
    input:
        expand("data/downloads/{area}/{year_month}/{data_type}.feather",
               area=enabled_areas,
               year_month=INTEGRATE_YEAR_MONTHS,
               data_type=config['data_types'])
    output:
        expand("data/integrated/{data_type}.feather", 
               data_type=config['data_types'])
    script:
        "scripts/integrate_data.py"

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
    conda: "environment.yml"
    script:
        "scripts/processing.py"


rule verify_data:
    input:
        "data/processed_data.feather"
    output:
        "results/data_availability.html"
    script:
        "scripts/verify_data.py"