wildcard_constraints:
    scenario=r"[^/]+",
    run=r"[^/]+",


rule solve_network:
    input:
        # The run's rows are read back from here by the script: they say which
        # input series is the origin's, which the destination's, and the route.
        scenario_table=lookup(dpath="inputs/scenarios", within=config),
        assumptions_base=lookup(dpath="inputs/assumptions", within=config),
        assumptions_overlay=optional("config/overlays/{scenario}.yaml"),
        # The plant sits at the origin's own representative point, so multi-site
        # runs need no hand-placed demand site. The origin is where the run's
        # renewables are.
        regions=collect(
            "resources/shapes/{item.area}_geo.parquet",
            item=lookup(
                query="scenario == '{scenario}' and run == '{run}' and tech != 'grid'",
                within=scenarios_df,
            ),
        ),
        # Every input series of the run, in the order its rows list them.
        timeseries=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and run == '{run}'",
                within=scenarios_df,
            ),
        ),
    output:
        network="results/{scenario}/.runs/{scenario}_{run}/outputs/{scenario}_{run}.nc",
    log:
        "logs/solve_network/{scenario}/{run}.log",
    script:
        "../scripts/solve/solve_network.py"
