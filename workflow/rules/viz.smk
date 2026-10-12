rule compile_report:
    """One run's report: what it built, cost and emitted, beside what priced it."""
    input:
        network="results/{scenario}/.runs/{scenario}_{run}/outputs/{scenario}_{run}.nc",
        # The run's rows: its origin, destination and route, and the variant each
        # tech was solved with. The solve reads the parquet, not its name, so the
        # table is where the variant survives.
        scenario_table=lookup(dpath="inputs/scenarios", within=config),
        # Emission factors, freight legs, the grid fee, and the capex and
        # lifetime quotes the cost leaves are split on, merged exactly as the
        # solve merged them, so a leaf is split on the numbers that priced it.
        assumptions_base=lookup(dpath="inputs/assumptions", within=config),
        assumptions_overlay=optional("config/overlays/{scenario}.yaml"),
        # The same grid series the solve priced against, the origin's and the
        # destination's, for the generation mix behind them.
        grid_input=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and run == '{run}' and tech == 'grid'",
                within=scenarios_df,
            ),
        ),
    output:
        report="results/{scenario}/.runs/{scenario}_{run}/report_{scenario}_{run}.csv",
        # Every field, including the ones the report holds back.
        diagnostic="results/{scenario}/.runs/{scenario}_{run}/.report_{scenario}_{run}_diag.csv",
        # What priced the run, beside the run.
        overlay="results/{scenario}/.runs/{scenario}_{run}/assumptions_{scenario}_{run}.overlay.yaml",
        base_assumptions="results/{scenario}/.runs/{scenario}_{run}/config/assumptions_base.yaml",
        merged_assumptions="results/{scenario}/.runs/{scenario}_{run}/config/assumptions_{scenario}_{run}_merged.yaml",
    log:
        "logs/compile_report/{scenario}/{run}.log",
    script:
        "../scripts/viz/compile_report.py"


rule pick_best_run:
    """Give a scenario the files of the run that stands for it.

    Usually its only run. An area on the grid that trades through market zones
    is run once per zone, and the cheapest by `report.best_run_by` stands for
    the scenario; every run keeps its own files under `.runs/`.
    """
    input:
        # Each named like the output it may become, one path per run.
        report=collect(
            "results/{item.scenario}/.runs/{item.scenario}_{item.run}/report_{item.scenario}_{item.run}.csv",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
        diagnostic=collect(
            "results/{item.scenario}/.runs/{item.scenario}_{item.run}/.report_{item.scenario}_{item.run}_diag.csv",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
        network=collect(
            "results/{item.scenario}/.runs/{item.scenario}_{item.run}/outputs/{item.scenario}_{item.run}.nc",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
        overlay=collect(
            "results/{item.scenario}/.runs/{item.scenario}_{item.run}/assumptions_{item.scenario}_{item.run}.overlay.yaml",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
        base_assumptions=collect(
            "results/{item.scenario}/.runs/{item.scenario}_{item.run}/config/assumptions_base.yaml",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
        merged_assumptions=collect(
            "results/{item.scenario}/.runs/{item.scenario}_{item.run}/config/assumptions_{item.scenario}_{item.run}_merged.yaml",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
    output:
        report="results/{scenario}/report_{scenario}.csv",
        diagnostic="results/{scenario}/.report_{scenario}_diag.csv",
        network="results/{scenario}/outputs/{scenario}.nc",
        overlay="results/{scenario}/assumptions_{scenario}.overlay.yaml",
        base_assumptions="results/{scenario}/config/assumptions_base.yaml",
        merged_assumptions="results/{scenario}/config/assumptions_{scenario}_merged.yaml",
        # Every run of the scenario, its ranking value, and which one was picked.
        runs="results/{scenario}/.runs/runs.csv",
    params:
        best_run_by=lookup(dpath="report/best_run_by", within=config),
    log:
        "logs/pick_best_run/{scenario}.log",
    script:
        "../scripts/viz/pick_best_run.py"


rule plot_cf_map:
    input:
        cutout="cutouts/{area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{area}_geo.parquet",
        offshore_regions="resources/shapes/{area}_offshore_geo.parquet",
        area_average_cf="resources/timeseries/{area}_{tech}_area-average_{start_date}_{end_date}.parquet",
    output:
        png="results/diag_plots/cf_map/{area}_{tech}_{start_date}_{end_date}_cf_map.png",
        html="results/diag_plots/cf_map/{area}_{tech}_{start_date}_{end_date}_cf_map.html",
    wildcard_constraints:
        tech=r"solar|wind-onshore|wind-offshore",
    log:
        "logs/plot_cf_map/{area}_{tech}_{start_date}_{end_date}.log",
    params:
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(dpath="res_cf/wind_offshore_turbine", within=config),
        region=lookup(dpath="areas/{area}/region", within=config),
    script:
        "../scripts/viz/diagnostics/plot_cf_map.py"


rule plot_capacity_bars:
    """One PNG/HTML per scenario — its runs go on the x-axis within each plot."""
    input:
        report="results/{scenario}/report_{scenario}.csv",
    output:
        png="results/{scenario}/plots/capacity_bars.png",
        html="results/{scenario}/plots/capacity_bars.html",
    log:
        "logs/plot_capacity_bars/{scenario}.log",
    script:
        "../scripts/viz/outputs/plot_capacity_bars.py"


rule plot_lcos_bars:
    """Stacked LCOS cost-breakdown per run for one scenario, requested on demand."""
    input:
        report="results/{scenario}/report_{scenario}.csv",
    output:
        png="results/{scenario}/plots/lcos_bars.png",
        html="results/{scenario}/plots/lcos_bars.html",
    log:
        "logs/plot_lcos_bars/{scenario}.log",
    script:
        "../scripts/viz/outputs/plot_lcos_bars.py"


rule plot_siting_map:
    """Multi-site only: geographic map of chosen sites + HVDC links for one run."""
    input:
        network="results/{scenario}/.runs/{scenario}_{run}/outputs/{scenario}_{run}.nc",
    output:
        png="results/{scenario}/.runs/{scenario}_{run}/plots/siting_map_{scenario}_{run}.png",
        html="results/{scenario}/.runs/{scenario}_{run}/plots/siting_map_{scenario}_{run}.html",
    log:
        "logs/plot_siting_map/{scenario}/{run}.log",
    script:
        "../scripts/viz/diagnostics/plot_siting_map.py"


rule plot_site_capacity_bars:
    """Multi-site only: per-site built capacity + HVDC link MW for one run."""
    input:
        network="results/{scenario}/.runs/{scenario}_{run}/outputs/{scenario}_{run}.nc",
    output:
        png="results/{scenario}/.runs/{scenario}_{run}/plots/site_capacity_bars_{scenario}_{run}.png",
        html="results/{scenario}/.runs/{scenario}_{run}/plots/site_capacity_bars_{scenario}_{run}.html",
    log:
        "logs/plot_site_capacity_bars/{scenario}/{run}.log",
    script:
        "../scripts/viz/outputs/plot_site_capacity_bars.py"
