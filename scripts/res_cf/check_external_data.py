"""
Validate that required external data files are present before running the
res_cf pipeline. Run this once after cloning / setting up the repo.

External files that Snakemake cannot download automatically:

  data/shapes/eez/eez_v12.shp
      World EEZ v12 from Marine Regions (marineregions.org).
      Free to download after brief registration:
        https://www.marineregions.org/eez.php
      → Download "World EEZ v12 (2023)" as a shapefile (.zip)
      → Extract the contents so that eez_v12.shp (and its companion
        .dbf / .shx / .prj files) live at data/shapes/eez/

  data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp
      Natural Earth 1:110m Admin-0 countries shapefile.
      Download from naturalearthdata.com or via:
        https://www.naturalearthdata.com/downloads/110m-cultural-vectors/
      → Extract to data/shapes/ne_110m_admin_0_countries/

Both datasets are stable; any version released within the last 2–3 years
will work. If you already have them from a different project, symlinks are
fine.
"""

from pathlib import Path
import sys

REQUIRED = {
    "EEZ (Marine Regions v12)": "data/shapes/eez/eez_v12.shp",
    "Natural Earth 110m countries": (
        "data/shapes/ne_110m_admin_0_countries/"
        "ne_110m_admin_0_countries.shp"
    ),
}

ROOT = Path(__file__).parent.parent.parent

ok = True
for label, rel_path in REQUIRED.items():
    p = ROOT / rel_path
    if p.exists():
        print(f"  [ok]  {label}")
    else:
        print(f"  [!!]  {label}")
        print(f"        Missing: {rel_path}")
        ok = False

if ok:
    print("\nAll external data files present.")
else:
    print("\nSome files are missing — see the docstring at the top of this")
    print("script or the README for download instructions.")
    sys.exit(1)


if __name__ == "__main__":
    pass
