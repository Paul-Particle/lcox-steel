wildcard_constraints:
    project=r"[^/]+",
    scenario=r"[^/]+",


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
        "../scripts/h2_dri/compile_report.py"


rule h2_dri_optimize:
    input:
        assumptions="config/assumptions.yaml",
        tech_inputs=collect(
            "resources/{item.pipeline}/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="project == '{project}' and scenario == '{scenario}'",
                within=projects_df,
            ),
        ),
    params:
        techs=collect(
            "{item.tech}",
            item=lookup(
                query="project == '{project}' and scenario == '{scenario}'",
                within=projects_df,
            ),
        ),
    output:
        network="results/{project}/{scenario}.nc",
    log:
        "logs/h2_dri_optimize/{project}_{scenario}.log",
    script:
        "../scripts/h2_dri/build_and_solve_network.py"
