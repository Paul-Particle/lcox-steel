rule download_entsoe:
    output:
        "resources/entsoe/{area}/{data_type}_{start_date}_{end_date}.parquet",
    resources:
        entsoe_api=2,
    script:
        "../scripts/grid/download_entsoe.py"


rule process_entsoe:
    input:
        expand(
            "resources/entsoe/{area}/{data_type}_{start_date}_{end_date}.parquet",
            data_type=config["entsoe"]["data_types"],
            allow_missing=True,
        ),
    output:
        prices="resources/entsoe/{area}_grid_dayahead_{start_date}_{end_date}.parquet",
        full="resources/entsoe/{area}_grid_dayahead_{start_date}_{end_date}_full.parquet",
    script:
        "../scripts/grid/process_entsoe.py"


rule process_nem:
    output:
        prices="resources/nem/{area}_grid_dayahead_{start_date}_{end_date}.parquet",
        full="resources/nem/{area}_grid_dayahead_{start_date}_{end_date}_full.parquet",
    script:
        "../scripts/grid/download_nem.py"
