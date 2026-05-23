"""Single source of truth for repo-relative path roots.

Scripts in `scripts/<subdir>/` import this with:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from _paths import DATA, RESOURCES, RES_CF  # noqa: E402

`parents[2]` resolves to the repo root in both standalone runs
(scripts/<subdir>/foo.py → scripts/<subdir> → scripts → repo) and
under Snakemake's tmp script copy (.snakemake/scripts/tmp.py → scripts → .snakemake → repo).
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# Raw / external / expensive (won't be re-fetched on rebuild)
DATA = REPO_ROOT / "data"
SHAPES_RAW = DATA / "shapes"                   # ne_110m, eez

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
