wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    variant=r"dayahead|full",


rule retrieve_entsoe:
    output:
        temp("resources/entsoe/{area}_grid_{variant}_{start_date}_{end_date}.parquet"),
    log:
        "logs/retrieve_entsoe/{area}_{variant}_{start_date}_{end_date}.log",
    resources:
        entsoe_api=2,
    script:
        "../scripts/grid/retrieve_entsoe.py"


rule retrieve_nem:
    output:
        temp("resources/nem/{area}_grid_{variant}_{start_date}_{end_date}.parquet"),
    log:
        "logs/retrieve_nem/{area}_{variant}_{start_date}_{end_date}.log",
    params:
        eur_per_aud=config["fx"]["eur_per_aud"],
    script:
        "../scripts/grid/retrieve_nem.py"
