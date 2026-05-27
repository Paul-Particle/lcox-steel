rule extract_eez_shapefile:
    input:
        "data/shapes/eez/eez_v12.zip"
    output:
        shp="data/shapes/eez/eez_v12.shp",
        shx="data/shapes/eez/eez_v12.shx",
        dbf="data/shapes/eez/eez_v12.dbf",
        prj="data/shapes/eez/eez_v12.prj",
    script:
        "../scripts/res_cf/extract_shapefile.py"


rule extract_ne_countries_shapefile:
    input:
        "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.zip"
    output:
        shp="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
        shx="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shx",
        dbf="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.dbf",
        prj="data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.prj",
    script:
        "../scripts/res_cf/extract_shapefile.py"


rule build_regions:
    input:
        "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
    output:
        "resources/shapes/regions.parquet"
    script:
        "../scripts/res_cf/build_regions.py"


rule build_offshore_regions:
    input:
        regions="resources/shapes/regions.parquet",
        eez="data/shapes/eez/eez_v12.shp"
    output:
        "resources/shapes/offshore_regions.parquet"
    script:
        "../scripts/res_cf/build_offshore_regions.py"


rule make_cutout:
    input:
        regions="resources/shapes/regions.parquet"
    output:
        "cutouts/{cf_area}_{start_date}_{end_date}.nc"
    params:
        start_date="{start_date}",
        end_date="{end_date}",
    script:
        "../scripts/res_cf/make_cutout.py"


rule build_cf_timeseries:
    input:
        cutout="cutouts/{cf_area}_{start_date}_{end_date}.nc",
        regions="resources/shapes/regions.parquet",
        offshore_regions="resources/shapes/offshore_regions.parquet",
    output:
        wind_onshore= "resources/res_cf/{cf_area}_wind_onshore_{start_date}_{end_date}.parquet",
        wind_offshore="resources/res_cf/{cf_area}_wind_offshore_{start_date}_{end_date}.parquet",
        solar=        "resources/res_cf/{cf_area}_solar_{start_date}_{end_date}.parquet",
    script:
        "../scripts/res_cf/build_cf_timeseries.py"
