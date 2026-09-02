wildcard_constraints:
    area=r"[A-Z][A-Z0-9_]{1,5}",
    tech=r"wind-onshore|wind-offshore|solar",


rule make_area_geometry:
    wildcard_constraints:
        area=COUNTRY_AREAS,
    input:
        "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip",
    output:
        "resources/shapes/{area}_geo.parquet",
    log:
        "logs/make_area_geometry/{area}.log",
    params:
        iso3=lookup(dpath="areas/{area}/iso3", within=config),
        region=lookup(dpath="areas/{area}/region", within=config),
        mainland_bbox=lookup(
            dpath="areas/{area}/mainland_bbox",
            within=config,
            default=None,
        ),
    script:
        "../scripts/res_cf/a_make_area_geometry.py"


rule make_offshore_geometry:
    input:
        regions="resources/shapes/{area}_geo.parquet",
        offshore_zone="data/shapes/offshore_zones/eez_v12.zip",
    output:
        "resources/shapes/{area}_offshore_geo.parquet",
    log:
        "logs/make_offshore_geometry/{area}.log",
    params:
        iso3=lookup(dpath="areas/{area}/iso3", within=config),
        region=lookup(dpath="areas/{area}/region", within=config),
        offshore_max_distance_km=lookup(
            dpath="res_cf/offshore_max_distance_km", within=config
        ),
    script:
        "../scripts/res_cf/b_make_offshore_geometry.py"


rule retrieve_area_cutout:
    input:
        # Cutout bounds are the bbox of the land ∪ offshore union (see
        # c_retrieve_area_cutout.py), so the cutout reaches far enough offshore
        # to cover the offshore-wind zone. Both geometries are read from the
        # pre-built parquets — the mainland_bbox filter / EEZ clip already
        # applied upstream — so no NE/EEZ zips are needed here.
        regions="resources/shapes/{area}_geo.parquet",
        offshore_regions="resources/shapes/{area}_offshore_geo.parquet",
    output:
        # Not protected(): the script preserves expensive downloads via the
        # `_backup.nc` sibling-file convention.
        "cutouts/{area}_{start_date}_{end_date}.nc",
    log:
        "logs/retrieve_area_cutout/{area}_{start_date}_{end_date}.log",
    params:
        coarse=lookup(
            dpath="areas/{area}/coarse", within=config, default=False
        ),
        bbox_pad_deg=lookup(dpath="res_cf/cutout/bbox_pad_deg", within=config),
        monthly_requests=lookup(dpath="res_cf/cutout/monthly_requests", within=config),
        cds_poll_interval_s=lookup(dpath="res_cf/cutout/cds_poll_interval_s", within=config),
        cache_warn_size_gb=lookup(dpath="res_cf/cutout/cache_warn_size_gb", within=config),
        min_free_disk_gb=lookup(dpath="res_cf/cutout/min_free_disk_gb", within=config),
    script:
        "../scripts/res_cf/c_retrieve_area_cutout.py"


rule area_average:
    input:
        cutout="cutouts/{area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{area}_geo.parquet",
        offshore_regions="resources/shapes/{area}_offshore_geo.parquet",
    output:
        "resources/timeseries/{area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        variant=r"area-average",
    log:
        "logs/area_average/{area}_{tech}_{variant}_{start_date}_{end_date}.log",
    params:
        region=lookup(dpath="areas/{area}/region", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(
            dpath="res_cf/wind_offshore_turbine", within=config
        ),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
        min_land_fraction=lookup(dpath="res_cf/min_land_fraction", within=config),
        eligibility_source=lookup(dpath="res_cf/eligibility_source", within=config),
    script:
        "../scripts/res_cf/d1_area_average.py"


rule bestsite_p95:
    input:
        cutout="cutouts/{area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{area}_geo.parquet",
        offshore_regions="resources/shapes/{area}_offshore_geo.parquet",
    output:
        "resources/timeseries/{area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        variant=r"bestsite-p95",
    log:
        "logs/bestsite_p95/{area}_{tech}_{variant}_{start_date}_{end_date}.log",
    params:
        region=lookup(dpath="areas/{area}/region", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(dpath="res_cf/wind_offshore_turbine", within=config),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
        spatial_matching_res_mix=lookup(
            dpath="res_cf/spatial_matching_res_mix", within=config
        ),
    script:
        "../scripts/res_cf/d2_bestsite_p95.py"


rule anchor_colo:
    input:
        cutout="cutouts/{area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{area}_geo.parquet",
        offshore_regions="resources/shapes/{area}_offshore_geo.parquet",
    output:
        "resources/timeseries/{area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        variant=r"anchor-colo-n\d+",
    log:
        "logs/anchor_colo/{area}_{tech}_{variant}_{start_date}_{end_date}.log",
    params:
        region=lookup(dpath="areas/{area}/region", within=config),
        anchor_colocation=lookup(dpath="res_cf/anchor_colocation", within=config),
    script:
        "../scripts/res_cf/d3_anchor_colo.py"


rule tilt_mix:
    input:
        cutout="cutouts/{area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{area}_geo.parquet",
    output:
        "resources/timeseries/{area}_solar_tilt-mix-n{n_steps}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        n_steps=r"\d+",
    log:
        "logs/tilt_mix/{area}_n{n_steps}_{start_date}_{end_date}.log",
    params:
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        region=lookup(dpath="areas/{area}/region", within=config),
    script:
        "../scripts/res_cf/d4_tilt_mix.py"


rule multi:
    input:
        cutout="cutouts/{area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{area}_geo.parquet",
        offshore_regions="resources/shapes/{area}_offshore_geo.parquet",
    output:
        "resources/timeseries/{area}_{tech}_multi-n{n_cells}_{start_date}_{end_date}.parquet",
    wildcard_constraints:
        n_cells=r"\d+",
    log:
        "logs/multi/{area}_{tech}_multi-n{n_cells}_{start_date}_{end_date}.log",
    params:
        region=lookup(dpath="areas/{area}/region", within=config),
        wind_onshore_turbine=lookup(dpath="res_cf/wind_onshore_turbine", within=config),
        wind_offshore_turbine=lookup(
            dpath="res_cf/wind_offshore_turbine", within=config
        ),
        pv_panel=lookup(dpath="res_cf/pv_panel", within=config),
        pv_orientation=lookup(dpath="res_cf/pv_orientation", within=config),
        wind_cf=lookup(dpath="res_cf/wind_cf", within=config),
    script:
        "../scripts/res_cf/d5_multi.py"
