wildcard_constraints:
    cf_area=r"[a-z]{2,3}",
    tech=r"wind-onshore|wind-offshore|solar",


rule build_regions:
    input:
        "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip",
    output:
        "resources/shapes/{cf_area}_geo.parquet",
    log:
        "logs/build_regions/{cf_area}.log",
    params:
        iso3=lookup(dpath="res_cf/countries/{cf_area}/iso3", within=config),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        mainland_bbox=lookup(
            dpath="res_cf/countries/{cf_area}/mainland_bbox",
            within=config,
            default=None,
        ),
    script:
        "../scripts/res_cf/01_build_regions.py"


rule build_offshore_regions:
    input:
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_zone="data/shapes/offshore_zones/eez_v12.zip",
    output:
        "resources/shapes/{cf_area}_offshore_geo.parquet",
    log:
        "logs/build_offshore_regions/{cf_area}.log",
    params:
        iso3=lookup(dpath="res_cf/countries/{cf_area}/iso3", within=config),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        offshore_max_distance_km=lookup(
            dpath="res_cf/offshore_max_distance_km", within=config
        ),
    script:
        "../scripts/res_cf/01b_build_offshore_regions.py"


rule download_cutout:
    input:
        ne_zip=ancient("data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip"),
    output:
        # Not protected(): the 02_make_cutouts.py script preserves expensive
        # downloads via the `_backup.nc` sibling-file convention.
        "cutouts/{cf_area}_{start_date}_{end_date}.nc",
    log:
        "logs/download_cutout/{cf_area}_{start_date}_{end_date}.log",
    params:
        iso3=lookup(dpath="res_cf/countries/{cf_area}/iso3", within=config),
        mainland_bbox=lookup(
            dpath="res_cf/countries/{cf_area}/mainland_bbox",
            within=config,
            default=None,
        ),
        coarse=lookup(
            dpath="res_cf/countries/{cf_area}/coarse", within=config, default=False
        ),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        bbox_pad_deg=lookup(dpath="res_cf/cutout/bbox_pad_deg", within=config),
        monthly_requests=lookup(dpath="res_cf/cutout/monthly_requests", within=config),
    script:
        "../scripts/res_cf/02_make_cutouts.py"


rule build_solar_tilt_mix_p95:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_solar_tilt-mix-n{n_steps}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        n_steps=r"\d+",
    log:
        "logs/build_solar_tilt_mix_p95/{cf_area}_n{n_steps}_{start_date}_{end_date}.log",
    params:
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
    script:
        "../scripts/res_cf/03b_build_solar_tilt_mix_p95.py"


rule build_res_cf_profile:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_country-average_{start_date}_{end_date}.parquet",
    log:
        "logs/build_res_cf_profile/{cf_area}_{tech}_{start_date}_{end_date}.log",
    params:
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(
            dpath="res_cf/wind_offshore_turbine", within=config
        ),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
    script:
        "../scripts/res_cf/03_build_cf_timeseries.py"


rule build_bestsite_p95:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        variant=r"bestsite-p95",
    log:
        "logs/build_bestsite_p95/{cf_area}_{tech}_{start_date}_{end_date}.log",
    params:
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(dpath="res_cf/wind_offshore_turbine", within=config),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
        spatial_matching_res_mix=lookup(
            dpath="res_cf/spatial_matching_res_mix", within=config
        ),
    script:
        "../scripts/res_cf/07_make_bestsite_cf_timeseries.py"


rule build_anchored_cf:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        variant=r"anchored-w\d+-s\d+",
    log:
        "logs/build_anchored_cf/{cf_area}_{tech}_{variant}_{start_date}_{end_date}.log",
    params:
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(dpath="res_cf/wind_offshore_turbine", within=config),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
        spatial_matching_res_mix=lookup(
            dpath="res_cf/spatial_matching_res_mix", within=config
        ),
    script:
        "../scripts/res_cf/07_make_bestsite_cf_timeseries.py"


rule build_res_cf_candidates:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_grid-n{n_cells}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        n_cells=r"\d+",
    log:
        "logs/build_res_cf_candidates/{cf_area}_{tech}_grid-n{n_cells}_{start_date}_{end_date}.log",
    params:
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(
            dpath="res_cf/wind_offshore_turbine", within=config
        ),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
    script:
        "../scripts/res_cf/03c_build_res_cf_candidates.py"
