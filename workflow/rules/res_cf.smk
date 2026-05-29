wildcard_constraints:
    cf_area=r"[a-z]{2,3}",
    tech=r"wind_onshore|wind_offshore|solar",


rule extract_offshore_zone_shapefile:
    input:
        "data/shapes/offshore_zones/eez_v12.zip",
    output:
        shp="data/shapes/offshore_zones/eez_v12.shp",
        shx="data/shapes/offshore_zones/eez_v12.shx",
        dbf="data/shapes/offshore_zones/eez_v12.dbf",
        prj="data/shapes/offshore_zones/eez_v12.prj",
    script:
        "../scripts/res_cf/extract_shapefile.py"


rule extract_ne_countries_shapefile:
    input:
        "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip",
    output:
        shp="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
        shx="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shx",
        dbf="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.dbf",
        prj="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.prj",
    script:
        "../scripts/res_cf/extract_shapefile.py"


rule build_regions:
    input:
        "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
    output:
        "resources/shapes/{cf_area}_geo.parquet",
    params:
        iso3=lookup(dpath="res_cf/countries/{cf_area}/iso3", within=config),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        mainland_bbox=lookup(
            dpath="res_cf/countries/{cf_area}/mainland_bbox",
            within=config,
            default=None,
        ),
    script:
        "../scripts/res_cf/build_regions.py"


rule build_offshore_regions:
    input:
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_zone="data/shapes/offshore_zones/eez_v12.shp",
    output:
        "resources/shapes/{cf_area}_offshore_geo.parquet",
    params:
        iso3=lookup(dpath="res_cf/countries/{cf_area}/iso3", within=config),
        region=lookup(dpath="res_cf/countries/{cf_area}/region", within=config),
        offshore_max_distance_km=lookup(
            dpath="res_cf/offshore_max_distance_km", within=config
        ),
    script:
        "../scripts/res_cf/build_offshore_regions.py"


rule download_cutout:
    input:
        ne_shp="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
    output:
        protected("cutouts/{cf_area}_{start_date}_{end_date}.nc"),
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
        "../scripts/res_cf/download_cutout.py"


rule build_cf_timeseries:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/{cf_area}_geo.parquet",
        offshore_regions="resources/shapes/{cf_area}_offshore_geo.parquet",
    output:
        "resources/res_cf/{cf_area}_{tech}_country-average_{start_date}_{end_date}.parquet",
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
        "../scripts/res_cf/build_cf_timeseries.py"
