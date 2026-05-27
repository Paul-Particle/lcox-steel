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
        # Map each scenario tech to its resolved registry path (order preserved so
        # run.py can zip techs ↔ files). tech_inputs templates were pre-resolved
        # per project in common.smk.
        tech_inputs=lambda wc: [
            projects["projects"][wc.project]["tech_inputs"][tech]
            for tech in projects["projects"][wc.project]["scenarios"][wc.scenario]["techs"]
        ],
        assumptions="config/assumptions.yaml",
        projects="config/projects.yaml",
    output:
        network="results/{project}/{scenario}.nc",
        summary="results/{project}/{scenario}_summary.csv",
    script:
        "../scripts/h2_dri/run.py"
