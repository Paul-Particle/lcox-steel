wildcard_constraints:
    scenario=r"[^/_]+",
    route=ROUTE_PATTERN,


rule solve_network:
    input:
        assumptions_base="config/assumptions.yaml",
        assumptions_overlay=optional(
            "config/assumptions_{scenario}.yaml"
        ),
        # The plant sits at the area's own representative point, so multi-site
        # scenarios need no hand-placed demand site.
        regions="resources/shapes/{area}_geo.parquet",
        # tech == 'grid' splits the run's rows into a price series and the
        # capacity-factor series; the CSV column does the classifying, so the
        # script never has to infer it from a path.
        cf_inputs=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and area == '{area}' "
                      "and start_date == '{start_date}' and end_date == '{end_date}' "
                      "and tech != 'grid'",
                within=scenarios_df,
            ),
        ),
        grid_input=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and area == '{area}' "
                      "and start_date == '{start_date}' and end_date == '{end_date}' "
                      "and tech == 'grid'",
                within=scenarios_df,
            ),
        ),
    output:
        network="results/{scenario}/{area}_{route}_{start_date}_{end_date}.nc",
    log:
        "logs/solve_network/{scenario}_{area}_{route}_{start_date}_{end_date}.log",
    script:
        "../scripts/solve/solve_network.py"
