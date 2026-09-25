"""Single source of truth for repo-relative path roots.

Scripts import this with:

    from common._paths import DATA, RESOURCES, TIMESERIES, CUTOUTS

`common` and `scripts` are importable because the repo is installed as a
package (see "Setup" in README.md).
"""

import os
from pathlib import Path

# workflow/common/_paths.py → workflow/common/ → workflow/ → repo root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Where the expensive inputs live. Defaults to the repo root, so a lone checkout
# behaves as it always has. Set `LCOX_STORE` to a checkout's path and every git
# worktree then reads and writes that one copy: a cutout is gigabytes and its CDS
# download queues for hours, so a per-worktree copy costs far more than it saves.
# Only what is expensive to fetch moves — `resources/` and `results/` stay per
# worktree, since they are derived and one branch's output is not another's.
STORE_ROOT = Path(os.environ.get("LCOX_STORE", REPO_ROOT))

# Raw / external / expensive (won't be re-fetched on rebuild)
DATA = STORE_ROOT / "data"
SHAPES_RAW = DATA / "shapes"                   # ne_110m, offshore_zone (eez_v12)

# Derived (Snakemake-tracked, reproducible from raw + scripts + config)
RESOURCES = REPO_ROOT / "resources"
# Every CF and grid-price series a scenario can consume, in one flat namespace
# keyed {area}_{tech}_{variant}_{start_date}_{end_date} — see config/scenarios.csv.
TIMESERIES = RESOURCES / "timeseries"
SHAPES_RES = RESOURCES / "shapes"              # regions.geojson, offshore_regions.geojson

# Atlite weather cutouts — derived in principle but expensive enough to treat as
# raw; lives at the repo root per PyPSA-Eur convention.
CUTOUTS = STORE_ROOT / "cutouts"

# Working directories
ATLITE_CACHE = STORE_ROOT / ".atlite-cache"    # atlite scratch (gitignored)
RESULTS = REPO_ROOT / "results"
