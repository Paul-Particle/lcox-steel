import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir)))


wildcard_constraints:
    start_date=r"\d{8}",
    end_date=r"\d{8}",
    cf_area=r"[a-z]{2,3}",
    data_type=r"prices|load_forecast|load_actual|res|generation|crossborder",
    nem_table=r"price|generation|load|crossborder",
    project=r"[^/]+",
    scenario=r"[^/]+",
