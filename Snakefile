# Snakefile for the cfp_data_analysis project

import yaml

# Load configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)


# Define the final output of the workflow
rule all:
    input:
        "results/02_analysis_and_visualization.html",


# Rule to download day-ahead prices
rule download_prices:
    output:
        "data/prices.feather",
    shell:
        "python scripts/download_entsoe_data.py prices {output} --start {config[start_date]} --end {config[end_date]}"


# Rule to download load forecast data
rule download_load:
    output:
        "data/load.feather",
    shell:
        "python scripts/download_entsoe_data.py load {output} --start {config[start_date]} --end {config[end_date]}"


# Rule to download VRE data
rule download_vre:
    output:
        "data/vre.feather",
    shell:
        "python scripts/download_entsoe_data.py vre {output} --start {config[start_date]} --end {config[end_date]}"


# Rule to process the raw data
rule process_data:
    input:
        prices="data/prices.feather",
        load="data/load.feather",
        vre="data/vre.feather",
        areas_config="areas.csv"
    output:
        "data/processed_data.feather",
    conda: "environment.yml"
    shell:
        "python scripts/processing.py {input.prices} {input.load} {input.vre} {output}"


# Rule to render the notebook
rule render_notebook:
    input:
        "notebooks/02_analysis_and_visualization.ipynb",
        "data/processed_data.feather",
    output:
        "/results/02_analysis_and_visualization.html",
    conda: "environment.yml"
    shell:
        "PYTHONPATH=$(pwd) jupyter nbconvert --to html --execute {input[0]} --output {output}"
