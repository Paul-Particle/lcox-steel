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
    script:
        "../scripts/h2_dri/compile_report.py"


rule h2_dri_optimize:
    input:
        tech_inputs=collect(
            "resources/{item.pipeline}/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
            item=lookup(
                query="project == '{project}' and scenario == '{scenario}'",
                within=projects_df,
            ),
        ),
        assumptions="config/assumptions.yaml",
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
    script:
        "../scripts/h2_dri/run.py"
