wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    variant=r"dayahead|full",


rule retrieve_grid_data:
    output:
        temp("resources/timeseries/{area}_grid_{variant}_{start_date}_{end_date}.parquet"),
    log:
        "logs/retrieve_grid_data/{area}_{variant}_{start_date}_{end_date}.log",
    params:
        eur_per_aud=config["nem"]["eur_per_aud"],
    resources:
        # ENTSO-E rate limit; harmless for NEM areas, which make no API calls.
        entsoe_api=2,
    script:
        "../scripts/grid/retrieve_grid_data.py"
