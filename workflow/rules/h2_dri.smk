wildcard_constraints:
    project=r"[^/]+",
    scenario=r"[^/]+",


rule h2_dri_optimize:
    input:
        assumptions_base="config/assumptions.yaml",
        assumptions_overlay=optional(
            "config/assumptions_{project}_{scenario}.yaml"
        ),
        tech_inputs=collect(
            "resources/{item.pipeline}/{item.area}_{item.tech}_{item.variant}_{item.start_date}_{item.end_date}.parquet",
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
