import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir)))

# Pre-resolve tech_inputs path templates for each project so that collect()
# in h2_dri.smk needs only a tech= kwarg — no cf_area/start_date/end_date kwargs.
# After this loop, projects["projects"][name]["tech_inputs"][tech] is a fully
# resolved path string with no remaining placeholders.
for _pname, _proj in config["projects"].items():
    _proj["tech_inputs"] = {
        tech: template.format(
            cf_area=_proj.get("cf_area", ""),
            start_date=_proj["start_date"],
            end_date=_proj["end_date"],
            bidding_zone=_proj.get("bidding_zone", ""),
        )
        for tech, template in config["tech_inputs"].items()
    }

enabled_areas = [
    code for code, info in config["entsoe"]["areas"].items() if info.get("enabled")
]

wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    cf_area=r"[a-z]{2,3}",
    project=r"[^/]+",
    scenario=r"[^/]+",
