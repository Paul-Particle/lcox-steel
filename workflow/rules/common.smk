import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir)))

enabled_areas = [
    code for code, info in config["entsoe"]["areas"].items() if info.get("enabled")
]

wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    cf_area=r"[a-z]{2,3}",
    project=r"[^/]+",
    scenario=r"[^/]+",
