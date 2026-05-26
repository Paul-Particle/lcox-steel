import sys
from pathlib import Path

# Put `workflow/` on sys.path so scripts under workflow/scripts/ can
# `from common._paths import …` (the common/ package lives at workflow/common/).
sys.path.insert(0, str(Path(workflow.basedir)))


enabled_areas = [code for code, info in config["entsoe"]["areas"].items() if info.get("enabled")]
CF_COUNTRIES  = [code for code, info in config["res_cf"]["countries"].items() if info.get("enabled")]
CF_YEAR       = config["res_cf"]["year"]
CF_QUARTERS   = ["q1", "q2", "q3", "q4"]
CF_TOP_N      = config["res_cf"]["top_n"]


# Constrain wildcards globally so e.g. `combine_techs` (output:
# {country}_cf_{year}.csv) cannot greedily match a bestsite_p95 filename and
# collide with `make_bestsite_cf`. `country` covers the lowercase iso2/iso3
# keys in config["res_cf"]["countries"].
wildcard_constraints:
    year=r"\d{4}",
    quarter=r"q[1-4]",
    country=r"[a-z]{2,3}",


def h2_dri_targets():
    """One `results/{project}/{scenario}_summary.csv` per (project, scenario) pair."""
    pairs = [(p["name"], s["name"])
             for p in config["projects"]
             for s in p["scenarios"]]
    projects, scenarios = zip(*pairs)
    return expand(
        "results/{project}/{scenario}_summary.csv",
        zip, project=projects, scenario=scenarios,
    )


def _find_project(name):
    for p in config["projects"]:
        if p["name"] == name:
            return p
    raise KeyError(name)


def _find_scenario(project, name):
    for s in project["scenarios"]:
        if s["name"] == name:
            return s
    raise KeyError(name)


def h2_dri_inputs(wildcards):
    """Resolve CF + (optional) prices inputs for a (project, scenario) target.

    The (project, scenario) wildcard pair is the seam where per-location
    wildcards land later without restructuring rules.
    """
    proj = _find_project(wildcards.project)
    scen = _find_scenario(proj, wildcards.scenario)

    cc   = proj["cf_area"]
    year = proj["time_period"]["start_date"][:4]
    sel  = scen.get("site_selection", "country-average")
    suffix   = "" if sel == "country-average" else f"_{sel}"
    cf_stem  = f"{cc}_cf_{year}{suffix}"

    inputs = {
        "cf":          f"resources/res_cf/{cf_stem}.parquet",
        "assumptions": "config/assumptions.yaml",
        "projects":    "config/projects.yaml",
    }
    if "grid" in scen.get("techs", []):
        inputs["prices"] = "resources/entsoe_processed.parquet"
    return inputs
