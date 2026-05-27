wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    data_type=r"prices|load_forecast|load_actual|res|generation|crossborder",
    nem_table=r"price|generation|load|crossborder",


rule download_entsoe:
    output:
        "resources/entsoe/{area}/{data_type}_{start_date}_{end_date}.parquet",
    resources:
        entsoe_api=2,
    script:
        "../scripts/grid/download_entsoe.py"


rule process_entsoe:
    input:
        "resources/entsoe/{area}/prices_{start_date}_{end_date}.parquet",
    output:
        prices="resources/entsoe/{area}_grid_dayahead_{start_date}_{end_date}.parquet",
    script:
        "../scripts/grid/process_entsoe.py"


rule process_entsoe_full:
    input:
        expand(
            "resources/entsoe/{area}/{data_type}_{start_date}_{end_date}.parquet",
            data_type=config["entsoe"]["data_types"],
            allow_missing=True,
        ),
    output:
        "resources/entsoe/{area}_grid_full_{start_date}_{end_date}.parquet",
    script:
        "../scripts/grid/process_entsoe_full.py"


rule download_nem:
    output:
        "resources/nem/raw/{nem_table}_{start_date}_{end_date}.parquet",
    script:
        "../scripts/grid/download_nem.py"


rule process_nem:
    input:
        "resources/nem/raw/price_{start_date}_{end_date}.parquet",
    output:
        prices="resources/nem/{area}_grid_dayahead_{start_date}_{end_date}.parquet",
    params:
        eur_per_aud=config["fx"]["eur_per_aud"],
    script:
        "../scripts/grid/process_nem.py"


rule process_nem_full:
    input:
        expand(
            "resources/nem/raw/{nem_table}_{start_date}_{end_date}.parquet",
            nem_table=["price", "load", "generation", "crossborder"],
            allow_missing=True,
        ),
    output:
        "resources/nem/{area}_grid_full_{start_date}_{end_date}.parquet",
    script:
        "../scripts/grid/process_nem_full.py"
