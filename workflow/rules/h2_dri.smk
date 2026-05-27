rule compile_report:
    input:
        networks=lambda wc: [
            f"results/{wc.project}/{s}.nc"
            for s in projects["projects"][wc.project]["scenarios"]
        ],
        assumptions="config/assumptions.yaml",
    output:
        "results/report_{project}.csv",
    script:
        "../scripts/h2_dri/compile_report.py"


rule h2_dri_optimize:
    input:
        tech_inputs=collect(
            lookup(
                dpath="projects/{project}/scenarios/{scenario}/tech_inputs/{tech}",
                within=projects,
                default=lookup(
                    dpath="projects/{project}/tech_inputs/{tech}",
                    within=projects,
                ),
            ),
            tech=lookup(
                dpath="projects/{project}/scenarios/{scenario}/techs",
                within=projects,
            ),
        ),
        assumptions="config/assumptions.yaml",
        projects="config/projects.yaml",
    output:
        network="results/{project}/{scenario}.nc",
        summary="results/{project}/{scenario}_summary.csv",
    script:
        "../scripts/h2_dri/run.py"
