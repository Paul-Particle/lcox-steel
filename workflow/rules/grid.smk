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
        eur_per_aud=config["nem"]["eur_per_aud"],
    script:
        "../scripts/grid/retrieve_nem.py"


rule retrieve_ons:
    output:
        temp("resources/ons/{area}_grid_{variant}_{start_date}_{end_date}.parquet"),
    log:
        "logs/retrieve_ons/{area}_{variant}_{start_date}_{end_date}.log",
    params:
        eur_per_brl=config["ons"]["eur_per_brl"],
        pld_limits=config["ons"]["pld_limits"],
    script:
        "../scripts/grid/retrieve_ons.py"
