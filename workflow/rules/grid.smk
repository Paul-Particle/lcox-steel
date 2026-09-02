wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    variant=r"dayahead|full",


rule retrieve_grid_data:
    output:
        # Not temp(): this file is the processed cache. It is keyed by the whole
        # run, so a re-run reuses it and a changed date range writes a new name
        # instead of invalidating anything.
        "resources/timeseries/{area}_grid_{variant}_{start_date}_{end_date}.parquet",
    log:
        "logs/retrieve_grid_data/{area}_{variant}_{start_date}_{end_date}.log",
    params:
        # Empty for an area with no price series — the script turns that into a
        # clear error rather than picking a source.
        market=lookup(dpath="areas/{area}/market", within=config, default=""),
        market_area=lookup(dpath="areas/{area}/market_area", within=config, default=""),
        eur_per_aud=config["nem"]["eur_per_aud"],
    resources:
        # ENTSO-E rate limit; harmless for NEM areas, which make no API calls.
        entsoe_api=2,
    script:
        "../scripts/grid/retrieve_grid_data.py"
