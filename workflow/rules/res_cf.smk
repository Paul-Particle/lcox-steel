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
        "resources/shapes/regions.geojson"
    script:
        "../scripts/res_cf/build_regions.py"


rule build_offshore_regions:
    input:
        regions="resources/shapes/regions.geojson",
        eez="data/shapes/eez/eez_v12.shp"
    output:
        "resources/shapes/offshore_regions.geojson"
    script:
        "../scripts/res_cf/build_offshore_regions.py"


rule make_cutout:
    input:
        regions="resources/shapes/regions.geojson"
    output:
        "cutouts/{country}_{year}_{quarter}.nc"
    script:
        "../scripts/res_cf/make_cutout.py"


rule build_cf_timeseries:
    input:
        cutout="cutouts/{country}_{year}_{quarter}.nc",
        regions="resources/shapes/regions.geojson"
    output:
        wind_onshore="resources/res_cf/quarterly/{country}_wind_onshore_{year}_{quarter}.parquet",
        wind_offshore="resources/res_cf/quarterly/{country}_wind_offshore_{year}_{quarter}.parquet",
        solar="resources/res_cf/quarterly/{country}_solar_{year}_{quarter}.parquet"
    script:
        "../scripts/res_cf/build_cf_timeseries.py"


rule concat_quarters:
    input:
        wind_onshore=expand("resources/res_cf/quarterly/{{country}}_wind_onshore_{{year}}_{quarter}.parquet",
                            quarter=CF_QUARTERS),
        wind_offshore=expand("resources/res_cf/quarterly/{{country}}_wind_offshore_{{year}}_{quarter}.parquet",
                             quarter=CF_QUARTERS),
        solar=expand("resources/res_cf/quarterly/{{country}}_solar_{{year}}_{quarter}.parquet",
                     quarter=CF_QUARTERS),
    output:
        wind_onshore="resources/res_cf/annual/{country}_wind_onshore_{year}.parquet",
        wind_offshore="resources/res_cf/annual/{country}_wind_offshore_{year}.parquet",
        solar="resources/res_cf/annual/{country}_solar_{year}.parquet"
    script:
        "../scripts/res_cf/concat_quarters.py"


rule combine_techs:
    input:
        wind_onshore="resources/res_cf/annual/{country}_wind_onshore_{year}.parquet",
        wind_offshore="resources/res_cf/annual/{country}_wind_offshore_{year}.parquet",
        solar="resources/res_cf/annual/{country}_solar_{year}.parquet"
    output:
        "resources/res_cf/{country}_cf_{year}.parquet"
    script:
        "../scripts/res_cf/combine_techs.py"


rule resource_spread:
    input:
        cutouts=expand("cutouts/{country}_{{year}}_{quarter}.nc",
                       country=CF_COUNTRIES, quarter=CF_QUARTERS),
        national_cfs=expand("resources/res_cf/{country}_cf_{{year}}.parquet",
                            country=CF_COUNTRIES),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output:
        "resources/res_cf/resource_spread_{year}.parquet"
    script:
        "../scripts/res_cf/resource_spread.py"


rule make_bestsite_cf:
    input:
        cutouts=expand("cutouts/{{country}}_{{year}}_{quarter}.nc",
                       quarter=CF_QUARTERS),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output:
        "resources/res_cf/{country}_cf_{year}_bestsite_p95.parquet"
    script:
        "../scripts/res_cf/make_bestsite_cf.py"


rule complementarity:
    input:
        national_cf="resources/res_cf/{country}_cf_{year}.parquet",
        cutouts=expand("cutouts/{{country}}_{{year}}_{quarter}.nc",
                       quarter=CF_QUARTERS),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output:
        # top_n is interpolated as a literal (not a wildcard) because the
        # paired `avg` output doesn't depend on N — Snakemake requires all
        # outputs of a rule to share the same wildcards.
        top=f"resources/res_cf/{{country}}_complementarity_top{CF_TOP_N}_{{year}}.parquet",
        avg="resources/res_cf/{country}_average_profiles_{year}.parquet"
    script:
        "../scripts/res_cf/complementarity.py"
