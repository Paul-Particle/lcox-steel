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
        # Cutout bounds are the bbox of the land ∪ offshore union (see
        # 02_make_cutouts.py), so the cutout reaches far enough offshore to
        # cover the offshore-wind zone. Both geometries are read from the
        # pre-built parquets — the mainland_bbox filter / EEZ clip already
        # applied by 01/01b — so no NE/EEZ zips are needed here.
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        # Not protected(): the 02_make_cutouts.py script preserves expensive
        # downloads via the `_backup.nc` sibling-file convention.
        "cutouts/{cf_area}_{start_date}_{end_date}.nc",
    log:
        "logs/download_cutout/{cf_area}_{start_date}_{end_date}.log",
    params:
        coarse=lookup(
            dpath="res_cf/countries/{cf_area}/coarse", within=config, default=False
        ),
        bbox_pad_deg=lookup(dpath="res_cf/cutout/bbox_pad_deg", within=config),
        monthly_requests=lookup(dpath="res_cf/cutout/monthly_requests", within=config),
        cds_poll_interval_s=lookup(dpath="res_cf/cutout/cds_poll_interval_s", within=config),
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


rule build_country_average_cf:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_country-average_{start_date}_{end_date}.parquet",
    log:
        "logs/build_country_average_cf/{cf_area}_{tech}_{start_date}_{end_date}.log",
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
        "logs/build_bestsite_p95/{cf_area}_{tech}_{variant}_{start_date}_{end_date}.log",
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


rule build_anchored_res_mix_cfs:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        variant=r"anchored-w\d+-s\d+",
    log:
        "logs/build_anchored_res_mix_cfs/{cf_area}_{tech}_{variant}_{start_date}_{end_date}.log",
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


rule build_complementarity_screen:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        screen="resources/res_cf/{cf_area}_{tech}_complementarity-screen_{start_date}_{end_date}.parquet",
        cf="resources/res_cf/{cf_area}_{tech}_complementarity-top1_{start_date}_{end_date}.parquet",
    log:
        "logs/build_complementarity_screen/{cf_area}_{tech}_{start_date}_{end_date}.log",
    params:
        top_n=lookup(dpath="res_cf/complementarity/top_n", within=config),
        coincidence_threshold=lookup(
            dpath="res_cf/complementarity/coincidence_threshold", within=config
        ),
        w_coincidence=lookup(dpath="res_cf/complementarity/w_coincidence", within=config),
        w_correlation=lookup(dpath="res_cf/complementarity/w_correlation", within=config),
        max_radius_km=lookup(dpath="res_cf/complementarity/max_radius_km", within=config),
        quality_floor=lookup(dpath="res_cf/complementarity/quality_floor", within=config),
        max_triplets_brute_force=lookup(
            dpath="res_cf/complementarity/max_triplets_brute_force", within=config
        ),
    script:
        "../scripts/res_cf/08_complementarity_screen.py"


rule build_multi_site_cfs:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_multi-n{n_cells}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        n_cells=r"\d+",
    log:
        "logs/build_multi_site_cfs/{cf_area}_{tech}_multi-n{n_cells}_{start_date}_{end_date}.log",
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
