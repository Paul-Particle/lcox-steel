wildcard_constraints:
    scenario=r"[^/]+",
    route=ROUTE_PATTERN,


rule solve_network:
    input:
        assumptions_base="config/assumptions.yaml",
        assumptions_overlay=optional(
            "config/assumptions_{scenario}.yaml"
        ),
        sites_overlay=optional(
            "config/sites_{scenario}.yaml"
        ),
        # tech == 'grid' splits the scenario's rows into a price series and the
        # capacity-factor series; the CSV column does the classifying, so the
        # script never has to infer it from a path.
        cf_inputs=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and tech != 'grid'",
                within=scenarios_df,
            ),
        ),
        grid_input=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and tech == 'grid'",
                within=scenarios_df,
            ),
        ),
    output:
        network="results/{scenario}/{route}.nc",
    log:
        "logs/solve_network/{scenario}_{route}.log",
    script:
        "../scripts/solve/solve_network.py"
