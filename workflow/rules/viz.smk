rule compile_report:
    input:
        networks=collect(
            "results/{item.scenario}/{item.area}_{item.route}_{item.start_date}_{item.end_date}.nc",
            item=lookup(query="scenario == '{scenario}'", within=runs_df),
        ),
    output:
        "results/report_{scenario}.csv",
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
        png="results/plots/cf_map/{area}_{tech}_{start_date}_{end_date}_cf_map.png",
        html="results/plots/cf_map/{area}_{tech}_{start_date}_{end_date}_cf_map.html",
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
        report="results/report_{scenario}.csv",
    output:
        png="results/plots/capacity_bars/{scenario}.png",
        html="results/plots/capacity_bars/{scenario}.html",
    log:
        "logs/plot_capacity_bars/{scenario}.log",
    script:
        "../scripts/viz/plot_capacity_bars.py"


rule plot_lcos_bars:
    """Steel routes only: stacked LCOS cost-breakdown per run for one scenario.

    Errors on scenarios whose only route is h2-only (no LCOS), so it is
    requested on demand rather than fanned out in `rule all`."""
    input:
        report="results/report_{scenario}.csv",
    output:
        png="results/plots/lcos_bars/{scenario}.png",
        html="results/plots/lcos_bars/{scenario}.html",
    log:
        "logs/plot_lcos_bars/{scenario}.log",
    script:
        "../scripts/viz/plot_lcos_bars.py"


rule plot_siting_map:
    """Multi-site only: geographic map of chosen sites + HVDC links for one run."""
    input:
        network="results/{scenario}/{area}_{route}_{start_date}_{end_date}.nc",
    output:
        png="results/plots/siting_map/{scenario}_{area}_{route}_{start_date}_{end_date}.png",
        html="results/plots/siting_map/{scenario}_{area}_{route}_{start_date}_{end_date}.html",
    log:
        "logs/plot_siting_map/{scenario}_{area}_{route}_{start_date}_{end_date}.log",
    script:
        "../scripts/viz/plot_siting_map.py"


rule plot_site_capacity_bars:
    """Multi-site only: per-site built capacity + HVDC link MW for one run."""
    input:
        network="results/{scenario}/{area}_{route}_{start_date}_{end_date}.nc",
    output:
        png="results/plots/site_capacity_bars/{scenario}_{area}_{route}_{start_date}_{end_date}.png",
        html="results/plots/site_capacity_bars/{scenario}_{area}_{route}_{start_date}_{end_date}.html",
    log:
        "logs/plot_site_capacity_bars/{scenario}_{area}_{route}_{start_date}_{end_date}.log",
    script:
        "../scripts/viz/plot_site_capacity_bars.py"
