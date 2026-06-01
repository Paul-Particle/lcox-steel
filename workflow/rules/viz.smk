rule compile_report:
    input:
        networks=collect(
            "results/{item.project}/{item.scenario}.nc",
            item=lookup(query="project == '{project}'", within=projects_df),
        ),
    output:
        "results/report_{project}.csv",
    log:
        "logs/compile_report/{project}.log",
    script:
        "../scripts/viz/compile_report.py"


rule plot_cf_map:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
    output:
        png="results/plots/cf_map/{cf_area}_{tech}_{start_date}_{end_date}_cf_map.png",
        html="results/plots/cf_map/{cf_area}_{tech}_{start_date}_{end_date}_cf_map.html",
    wildcard_constraints:
        tech=r"solar|wind-onshore|wind-offshore",
    log:
        "logs/plot_cf_map/{cf_area}_{tech}_{start_date}_{end_date}.log",
    params:
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
    script:
        "../scripts/viz/plot_cf_map.py"


rule plot_capacity_bars:
    """One PNG/HTML per project — scenarios go on the x-axis within each plot."""
    input:
        reports=["results/report_{project}.csv"],
    output:
        png="results/plots/capacity_bars/{project}.png",
        html="results/plots/capacity_bars/{project}.html",
    log:
        "logs/plot_capacity_bars/{project}.log",
    script:
        "../scripts/viz/plot_capacity_bars.py"
