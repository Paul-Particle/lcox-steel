rule download_entsoe:
    output:
        "resources/entsoe/{area}/{data_type}.parquet"
    params:
        start_date=config["entsoe"]["start_date"],
        end_date=config["entsoe"]["end_date"],
        cache_dir="data/entsoe_cache",
    resources:
        entsoe_api=2
    script:
        "../scripts/grid/download_entsoe.py"


rule download_nem:
    output:
        "resources/nem_processed.parquet"
    params:
        start_date=config["nem_download"]["start_date"],
        end_date=config["nem_download"]["end_date"],
        cache_dir=config["nem_download"]["cache_dir"],
        resample_freq=config["nem_download"].get("resample_freq"),
        rebuild=config["nem_download"]["rebuild"]
    script:
        "../scripts/grid/download_nem.py"


rule process_entsoe:
    input:
        entsoe=expand("resources/entsoe/{area}/{data_type}.parquet",
                      area=enabled_areas,
                      data_type=config["entsoe"]["data_types"])
    params:
        areas=enabled_areas
    output:
        "resources/entsoe_processed.parquet"
    script:
        "../scripts/grid/process_entsoe.py"
