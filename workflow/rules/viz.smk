rule compile_report:
    input:
        # The scenario table, for the variant each tech was solved with: the
        # solve reads the parquet, not its name, so this is where that survives.
        scenarios="config/scenarios.csv",
        # Emission factors, freight legs, the grid fee, and the capex and
        # lifetime quotes the cost leaves are split on. The scenario's own
        # overlay is merged over the base here exactly as the solve merged it,
        # so a leaf is split on the numbers that priced it: the gas-price and
        # EW-capex sweeps and the salt-cavern buffer all move a quote the split
        # reads, and the base file alone would divide a solved cost by a
        # quote the run never saw.
        assumptions_base="config/assumptions.yaml",
        assumptions_overlay=optional(
            "config/assumptions_{scenario}.yaml"
        ),
        networks=collect(
            "results/{item.scenario}/outputs/{item.area}_{item.route}_{item.start_date}_{item.end_date}.nc",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
        # The same grid series the solve priced its imports against, for the
        # generation mix behind them. Empty for an islanded scenario, and
        # carrying prices only unless the row asks for `variant: full`.
        grid_input=collect(
            "resources/timeseries/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}' and tech == 'grid'",
                within=scenarios_df,
            ),
        ),
        # And the destination market's series, for the same reason at the other
        # end of an `-export` route: the furnace there carries its own market's
        # mix, not the one that made the iron. Empty for a scenario that builds
        # no export twin.
        destination_input=collect(
            "resources/timeseries/{item.destination}_grid_emissions_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="scenario == '{scenario}'",
                within=destinations_df,
            ),
        ),
    output:
        # The report stands on its own: one row per reported place, the zone
        # ranking already resolved. The diagnostic keeps every zone and the
        # `best_in_country` flag, and is hidden because it answers a follow-up
        # question rather than being the thing to read.
        report="results/{scenario}/report_{scenario}.csv",
        diagnostic="results/{scenario}/.report_{scenario}_diag.csv",
        # What priced the run, beside the run. The overlay is written even when
        # the scenario has none, because "this changed nothing" is a fact worth
        # recording and an absent file cannot say it.
        overlay="results/{scenario}/assumptions_{scenario}.overlay.yaml",
        base_assumptions="results/{scenario}/config/assumptions_base.yaml",
        merged_assumptions="results/{scenario}/config/assumptions_{scenario}_merged.yaml",
    params:
        best_zone_by=lookup(dpath="report/best_zone_by", within=config, default=""),
    log:
        "logs/compile_report/{scenario}.log",
    script:
        "../scripts/viz/compile_report.py"


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
        "../scripts/viz/plot_cf_map.py"


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
        "../scripts/viz/plot_capacity_bars.py"


rule plot_lcos_bars:
    """Steel routes only: stacked LCOS cost-breakdown per run for one scenario.

    Errors on scenarios whose only route is h2-only (no LCOS), so it is
    requested on demand rather than fanned out in `rule all`."""
    input:
        report="results/{scenario}/report_{scenario}.csv",
    output:
        png="results/{scenario}/plots/lcos_bars.png",
        html="results/{scenario}/plots/lcos_bars.html",
    log:
        "logs/plot_lcos_bars/{scenario}.log",
    script:
        "../scripts/viz/plot_lcos_bars.py"


rule plot_siting_map:
    """Multi-site only: geographic map of chosen sites + HVDC links for one run."""
    input:
        network="results/{scenario}/outputs/{area}_{route}_{start_date}_{end_date}.nc",
    output:
        png="results/{scenario}/plots/siting_map_{area}_{route}_{start_date}_{end_date}.png",
        html="results/{scenario}/plots/siting_map_{area}_{route}_{start_date}_{end_date}.html",
    log:
        "logs/plot_siting_map/{scenario}_{area}_{route}_{start_date}_{end_date}.log",
    script:
        "../scripts/viz/plot_siting_map.py"


rule plot_site_capacity_bars:
    """Multi-site only: per-site built capacity + HVDC link MW for one run."""
    input:
        network="results/{scenario}/outputs/{area}_{route}_{start_date}_{end_date}.nc",
    output:
        png="results/{scenario}/plots/site_capacity_bars_{area}_{route}_{start_date}_{end_date}.png",
        html="results/{scenario}/plots/site_capacity_bars_{area}_{route}_{start_date}_{end_date}.html",
    log:
        "logs/plot_site_capacity_bars/{scenario}_{area}_{route}_{start_date}_{end_date}.log",
    script:
        "../scripts/viz/plot_site_capacity_bars.py"
