import yaml
import pandas as pd

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

if config.get('scheduler') != 'greedy':
    print("WARNING: You are not using the greedy scheduler which might cause too many API calls from parallel execution of many small downloads. Using --scheduler greedy will ensure sequential downloads. Forcing greedy.")
    config['scheduler'] = 'greedy'

def get_enabled_areas():
    areas_df = pd.read_csv("areas.csv")
    return areas_df[areas_df["enabled"]]["area_code"].tolist()

enabled_areas = get_enabled_areas()

rule all:
    input:
        "results/data_availability.html",
        # "notebooks/results/CFP_analysis_and_visualization.html",

# Rule to download data for a specific area, year, and month
rule download_split:
    output:
        "data/downloads/{area}/{year}-{month}/{data_type}.feather"
    log:
        "logs/download/{area}-{year}-{month}-{data_type}.log"
    shell:
        "(python -u scripts/download_entsoe_data_split.py {wildcards.data_type} data/downloads {wildcards.area} {wildcards.year} {wildcards.month}) &> {log}"
        "python -u scripts/download_entsoe_data_split.py {wildcards.data_type} data/downloads {wildcards.area} {wildcards.year} {wildcards.month}"

# Rule to integrate the downloaded data
rule integrate_data:
    input:
        expand("data/downloads/{area}/{year}-{month}/{data_type}.feather",
               area=enabled_areas,
               year=config['split_download']['years'],
               month=config['split_download']['months'],
               data_type=config['data_types'])
    output:
        expand("data/integrated/{data_type}.feather", data_type=config['data_types'])
    params:
        areas=' '.join(enabled_areas),
        year_months=' '.join(config['integration']['year_months'])
    shell:
        "python -u scripts/integrate_data.py --output-dir data/integrated --inputs {input}"
        "python -u scripts/integrate_data.py data/downloads data/integrated --data_types {config[data_types]} --areas {params.areas} --year_months {params.year_months}"
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
    shell:
        "python -u scripts/processing.py {input.prices} {input.load_forecast} {input.load_actual} {input.vre} {input.generation} {input.crossborder} {output}"
    script:
        "scripts/processing.py"

# rule render_notebook:
#     input:
#         notebook="notebooks/CFP_analysis_and_visualization.ipynb",
#         data="data/processed_data.feather"
#     output:
#         "notebooks/results/CFP_analysis_and_visualization.html"
#     params:
#         intermediate_html=lambda wildcards, input: str(input.notebook).replace(".ipynb", ".html")
#     shell:
#         '''
#         jupyter nbconvert --to html --execute {input.notebook}
#         mv {params.intermediate_html} {output}
#         '''

rule verify_data:
    input:ts"
    script:
        "scrip/verify_data.py
        "data/processed_data.feather"
    output:
        "results/data_availability.html"
    shell:
        "python -u scripts/verify_data.py {input} results"