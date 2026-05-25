rule h2_dri_optimize:
    input:
        unpack(h2_dri_inputs)
    output:
        network="results/{project}/{scenario}.nc",
        summary="results/{project}/{scenario}_summary.csv",
    script:
        "../scripts/h2_dri/run.py"
