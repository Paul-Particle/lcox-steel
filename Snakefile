# Snakefile for the cfp_data_analysis project

import yaml

# Load configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)


# Define the final output of the workflow
rule all:
    input:
        "notebooks/results/02_analysis_and_visualization.html",


# Rule to download day-ahead prices
rule download_prices:
    output:
        "data/prices.feather",
    shell:
        "python -u scripts/download_entsoe_data.py prices {output} --start {config[start_date]} --end {config[end_date]}"


# Rule to download load forecast data
rule download_load_forecast:
    output:
        "data/load_fc.feather",
    shell:
        "python -u scripts/download_entsoe_data.py load_forecast {output} --start {config[start_date]} --end {config[end_date]}"

# Rule to download actual load data
rule download_load_actual:
    output:
        "data/load_ac.feather",
    shell:
        "python -u scripts/download_entsoe_data.py load_actual {output} --start {config[start_date]} --end {config[end_date]}"

# Rule to download vre forecast data
rule download_vre:
    output:
        "data/vre.feather",
    shell:
        "python -u scripts/download_entsoe_data.py vre {output} --start {config[start_date]} --end {config[end_date]}"

# Rule to download actual generation data
rule download_gen:
    output:
        "data/gen.feather",
    shell:
        "python -u scripts/download_entsoe_data.py generation {output} --start {config[start_date]} --end {config[end_date]}"

# Rule to download crossborder export import data
rule download_xim:
    output:
        "data/xim.feather",
    shell:
        "python -u scripts/download_entsoe_data.py crossborder {output} --start {config[start_date]} --end {config[end_date]}"




# Rule to process the raw data
rule process_data:
    input:
        prices="data/prices.feather",
        load_fc="data/load_fc.feather",
        load_ac="data/load_ac.feather",
        vre="data/vre.feather",
        gen="data/gen.feather",
        xim="data/xim.feather",
        areas_config="areas.csv"
    output:
        "data/processed_data.feather",
    conda: "environment.yml"
    shell:
        "python -u scripts/processing.py {input.prices} {input.load_fc} {input.load_ac} {input.vre} {input.gen} {input.xim} {output}"


# Rule to render the notebook
rule render_notebookxim:
    input:
        notebook="notebooks/02_analysis_and_visualization.ipynb",
        data="data/processed_data.feather"
    output:
        "notebooks/results/02_analysis_and_visualization.html"
    params:
        intermediate_html=lambda wildcards, input: str(input.notebook).replace(".ipynb", ".html")
    shell:
        """
        jupyter nbconvert --to html --execute {input.notebook}
        mv {params.intermediate_html} {output}
        """