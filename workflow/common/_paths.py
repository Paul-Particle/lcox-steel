"""Single source of truth for repo-relative path roots.

Scripts import this with:

    from common._paths import DATA, RESOURCES, TIMESERIES, CUTOUTS

`common` and `scripts` are importable because the repo is installed as a
package (see "Setup" in README.md).
"""

import os
from pathlib import Path

from dotenv import dotenv_values

# workflow/common/_paths.py → workflow/common/ → workflow/ → repo root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# The checkout every git worktree of this repo hangs off. A worktree's `.git` is
# a file reading `gitdir: <main>/.git/worktrees/<name>`; the main checkout's is
# the directory itself, and a copy without `.git` counts as its own main.
_git_entry = REPO_ROOT / ".git"
MAIN_CHECKOUT = (
    (REPO_ROOT / _git_entry.read_text().removeprefix("gitdir:").strip()).resolve().parents[2]
    if _git_entry.is_file()
    else REPO_ROOT
)

# Where the expensive inputs live, shared by every worktree: a cutout is gigabytes
# and its CDS download queues for hours, so a per-worktree copy costs far more than
# it saves. `LCOX_STORE` names it, from the environment or else from the main
# checkout's `.env`, so one line there serves all its worktrees; unset, the store
# is the main checkout itself. Only what is expensive to fetch moves —
# `resources/` and `results/` stay per worktree, since they are derived and one
# branch's output is not another's.
STORE_ROOT = Path(
    os.environ.get("LCOX_STORE")
    or dotenv_values(MAIN_CHECKOUT / ".env").get("LCOX_STORE")
    or MAIN_CHECKOUT
)

# Raw / external / expensive (won't be re-fetched on rebuild)
DATA = STORE_ROOT / "data"
SHAPES_RAW = DATA / "shapes"                   # ne_110m, offshore_zone (eez_v12)

# Derived (Snakemake-tracked, reproducible from raw + scripts + config)
RESOURCES = REPO_ROOT / "resources"
# Every CF and grid-price series a scenario can consume, in one flat namespace
# keyed {area}_{tech}_{variant}_{start_date}_{end_date}, one per row of a scenario table.
TIMESERIES = RESOURCES / "timeseries"
SHAPES_RES = RESOURCES / "shapes"              # regions.geojson, offshore_regions.geojson

# Atlite weather cutouts — derived in principle but expensive enough to treat as
# raw; lives at the repo root per PyPSA-Eur convention.
CUTOUTS = STORE_ROOT / "cutouts"

# Working directories
ATLITE_CACHE = STORE_ROOT / ".atlite-cache"    # atlite scratch (gitignored)
RESULTS = REPO_ROOT / "results"
