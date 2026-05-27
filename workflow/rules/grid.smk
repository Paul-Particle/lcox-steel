# Download date range covers the union of all project periods.
_dl_start = min(p["start_date"] for p in config["projects"].values())
_dl_end   = max(p["end_date"]   for p in config["projects"].values())


rule download_entsoe:
    output:
        "resources/entsoe/{area}/{data_type}.parquet"
    params:
        start_date=_dl_start,
        end_date=_dl_end,
        cache_dir="data/entsoe_cache",
    resources:
        entsoe_api=2
    script:
        "../scripts/grid/download_entsoe.py"


rule process_entsoe:
    input:
        lambda wc: expand(
            "resources/entsoe/{area}/{data_type}.parquet",
            area=wc.bidding_zone,
            data_type=config["entsoe"]["data_types"],
        )
    params:
        start_date="{start_date}",
        end_date="{end_date}",
    output:
        prices="resources/entsoe/{bidding_zone}_{start_date}_{end_date}.parquet",
        full="resources/entsoe/{bidding_zone}_{start_date}_{end_date}_full.parquet",
    script:
        "../scripts/grid/process_entsoe.py"


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
