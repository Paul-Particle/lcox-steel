import subprocess
import sys
import yaml

configfile: "config/config.yaml"

enabled_areas = [code for code, info in config["entsoe"]["areas"].items() if info.get("enabled")]
CF_COUNTRIES  = [code for code, info in config["res_cf"]["countries"].items() if info.get("enabled")]
CF_YEAR       = config["res_cf"]["year"]
CF_QUARTERS   = ["q1", "q2", "q3", "q4"]
CF_TOP_N      = config["res_cf"]["top_n"]

# h2_dri target expansion — driven by projects.yaml (one optimization per
# (project, scenario) pair). Adding a project/scenario is a config edit, not a
# Snakefile edit.
with open("config/projects.yaml") as _f:
    PROJECTS_CFG = yaml.safe_load(_f)

H2_TARGETS = [
    f"results/{p['name']}/{s['name']}_summary.csv"
    for p in PROJECTS_CFG["projects"]
    for s in p["scenarios"]
]

# Constrain `year` globally to 4 digits so that e.g. `combine_techs`
# (output: {country}_cf_{year}.csv) does not greedily match a bestsite_p95
# filename and collide with the dedicated `make_bestsite_cf` rule.
wildcard_constraints:
    year=r"\d{4}",
    quarter=r"q[1-4]",


onstart:
    # Refuse to start if required external shapefiles are missing. The check
    # script also auto-extracts any ZIPs found in the data/shapes/* dirs.
    rc = subprocess.call([sys.executable, "scripts/res_cf/check_external_data.py"])
    if rc != 0:
        raise WorkflowError(
            "External data check failed. Drop the EEZ and Natural Earth ZIPs "
            "into data/shapes/<dataset>/ and re-run."
        )


rule all:
    input:
        # Grid pipeline
        "resources/entsoe_processed.feather",
        "resources/nem_processed.feather",
        expand("resources/entsoe/{area}/{data_type}.feather",
               area=enabled_areas,
               data_type=config['entsoe']['data_types']),
        # res_cf pipeline
        expand("resources/res_cf/{country}_cf_{year}.csv",
               country=CF_COUNTRIES, year=CF_YEAR),
        expand("resources/res_cf/{country}_cf_{year}_bestsite_p95.csv",
               country=CF_COUNTRIES, year=CF_YEAR),
        f"resources/res_cf/resource_spread_{CF_YEAR}.csv",
        expand("resources/res_cf/{country}_complementarity_top{n}_{year}.csv",
               country=CF_COUNTRIES, n=CF_TOP_N, year=CF_YEAR),
        # h2_dri pipeline (one (project, scenario) per target)
        H2_TARGETS,


# ── Grid pipeline ──────────────────────────────────────────────────────────────

rule download_entsoe:
    output:
        "resources/entsoe/{area}/{data_type}.feather"
    params:
        start_date=config["entsoe"]["start_date"],
        end_date=config["entsoe"]["end_date"],
        cache_dir="data/entsoe_cache",
    resources:
        entsoe_api=1
    script:
        "scripts/grid/download_entsoe.py"

rule download_nem:
    output:
        "resources/nem_processed.feather"
    params:
        start_date=config["nem_download"]["start_date"],
        end_date=config["nem_download"]["end_date"],
        cache_dir=config["nem_download"]["cache_dir"],
        resample_freq=config["nem_download"].get("resample_freq"),
        rebuild=config["nem_download"]["rebuild"]
    script: "scripts/grid/download_nem.py"

rule process_entsoe:
    input:
        entsoe=expand("resources/entsoe/{area}/{data_type}.feather",
                      area=enabled_areas,
                      data_type=config['entsoe']['data_types'])
    params:
        areas=enabled_areas
    output:
        "resources/entsoe_processed.feather",
    script:
        "scripts/grid/process_entsoe.py"


# ── res_cf pipeline ─────────────────────────────────────────────────────────────

rule build_regions:
    input: "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp"
    output: "resources/shapes/regions.geojson"
    script: "scripts/res_cf/build_regions.py"

rule build_offshore_regions:
    input:
        regions="resources/shapes/regions.geojson",
        eez="data/shapes/eez/eez_v12.shp"
    output: "resources/shapes/offshore_regions.geojson"
    script: "scripts/res_cf/build_offshore_regions.py"

rule make_cutout:
    input:
        regions="resources/shapes/regions.geojson"
    output: "cutouts/{country}_{year}_{quarter}.nc"
    script: "scripts/res_cf/make_cutout.py"

rule build_cf_timeseries:
    input:
        cutout="cutouts/{country}_{year}_{quarter}.nc",
        regions="resources/shapes/regions.geojson"
    output:
        wind_onshore="resources/res_cf/quarterly/{country}_wind_onshore_{year}_{quarter}.csv",
        wind_offshore="resources/res_cf/quarterly/{country}_wind_offshore_{year}_{quarter}.csv",
        solar="resources/res_cf/quarterly/{country}_solar_{year}_{quarter}.csv"
    script: "scripts/res_cf/build_cf_timeseries.py"

rule concat_quarters:
    input:
        wind_onshore=expand("resources/res_cf/quarterly/{{country}}_wind_onshore_{{year}}_{quarter}.csv",
                            quarter=CF_QUARTERS),
        wind_offshore=expand("resources/res_cf/quarterly/{{country}}_wind_offshore_{{year}}_{quarter}.csv",
                             quarter=CF_QUARTERS),
        solar=expand("resources/res_cf/quarterly/{{country}}_solar_{{year}}_{quarter}.csv",
                     quarter=CF_QUARTERS),
    output:
        wind_onshore="resources/res_cf/annual/{country}_wind_onshore_{year}.csv",
        wind_offshore="resources/res_cf/annual/{country}_wind_offshore_{year}.csv",
        solar="resources/res_cf/annual/{country}_solar_{year}.csv"
    script: "scripts/res_cf/concat_quarters.py"

rule combine_techs:
    input:
        wind_onshore="resources/res_cf/annual/{country}_wind_onshore_{year}.csv",
        wind_offshore="resources/res_cf/annual/{country}_wind_offshore_{year}.csv",
        solar="resources/res_cf/annual/{country}_solar_{year}.csv"
    output: "resources/res_cf/{country}_cf_{year}.csv"
    script: "scripts/res_cf/combine_techs.py"

rule resource_spread:
    input:
        cutouts=expand("cutouts/{country}_{year}_{quarter}.nc",
                       country=CF_COUNTRIES, year=CF_YEAR, quarter=CF_QUARTERS),
        national_cfs=expand("resources/res_cf/{country}_cf_{year}.csv",
                            country=CF_COUNTRIES, year=CF_YEAR),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output: f"resources/res_cf/resource_spread_{CF_YEAR}.csv"
    script: "scripts/res_cf/resource_spread.py"

rule make_bestsite_cf:
    input:
        cutouts=expand("cutouts/{{country}}_{year}_{quarter}.nc",
                       year=CF_YEAR, quarter=CF_QUARTERS),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output: f"resources/res_cf/{{country}}_cf_{CF_YEAR}_bestsite_p95.csv"
    script: "scripts/res_cf/make_bestsite_cf.py"

rule complementarity:
    input:
        national_cf=f"resources/res_cf/{{country}}_cf_{CF_YEAR}.csv",
        cutouts=expand("cutouts/{{country}}_{year}_{quarter}.nc",
                       year=CF_YEAR, quarter=CF_QUARTERS),
        regions="resources/shapes/regions.geojson",
        offshore_regions="resources/shapes/offshore_regions.geojson"
    output:
        top=f"resources/res_cf/{{country}}_complementarity_top{CF_TOP_N}_{CF_YEAR}.csv",
        avg=f"resources/res_cf/{{country}}_average_profiles_{CF_YEAR}.csv"
    script: "scripts/res_cf/complementarity.py"


# ── h2_dri pipeline ─────────────────────────────────────────────────────────────

def _find(items, name):
    for it in items:
        if it["name"] == name:
            return it
    raise KeyError(name)


def h2_dri_inputs(wildcards):
    """Resolve CF + (optional) prices inputs from projects.yaml for a (project, scenario).

    Single-cell shape today; the (project, scenario) wildcard pair is the seam
    where per-location wildcards land later without restructuring rules.
    """
    proj = _find(PROJECTS_CFG["projects"], wildcards.project)
    scen = _find(proj["scenarios"], wildcards.scenario)

    cc = proj["country"].lower()
    variant = proj.get("cf_variant", "avg")
    cf_name = f"{cc}_cf_{proj['year']}.csv" if variant == "avg" else f"{cc}_cf_{proj['year']}_{variant}.csv"

    inputs = {
        "cf":          f"resources/res_cf/{cf_name}",
        "assumptions": "config/assumptions.yaml",
        "projects":    "config/projects.yaml",
    }
    if scen.get("grid_connected"):
        inputs["prices"] = "resources/entsoe_processed.feather"
    return inputs


rule h2_dri_optimize:
    input: unpack(h2_dri_inputs)
    output:
        network = "results/{project}/{scenario}.nc",
        summary = "results/{project}/{scenario}_summary.csv",
    script: "scripts/h2_dri/run.py"
