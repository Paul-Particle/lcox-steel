"""Single source of truth for repo-relative path roots.

Scripts import this with:

    from common._paths import DATA, RESOURCES, RES_CF, CUTOUTS

This works because `workflow/` is on sys.path (added by
`workflow/rules/common.smk` under Snakemake; scripts that run standalone add it
themselves). For standalone runs from elsewhere, add
`sys.path.insert(0, "<repo>/workflow")` first.
"""

from pathlib import Path

# workflow/common/_paths.py → workflow/common/ → workflow/ → repo root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Raw / external / expensive (won't be re-fetched on rebuild)
DATA = REPO_ROOT / "data"
SHAPES_RAW = DATA / "shapes"                   # ne_110m, offshore_zone (eez_v12)

# Derived (Snakemake-tracked, reproducible from raw + scripts + config)
RESOURCES = REPO_ROOT / "resources"
RES_CF = RESOURCES / "res_cf"
SHAPES_RES = RESOURCES / "shapes"              # regions.geojson, offshore_regions.geojson

# Atlite weather cutouts — derived in principle but expensive enough to treat as
# raw; lives at the repo root per PyPSA-Eur convention.
CUTOUTS = REPO_ROOT / "cutouts"

# Working directories
ATLITE_CACHE = REPO_ROOT / ".atlite-cache"     # atlite scratch (gitignored)
RESULTS = REPO_ROOT / "results"
