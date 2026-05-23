"""
Validate (and set up) required external data files for the res_cf pipeline.

For each dataset:
  1. Auto-creates the target directory if missing.
  2. Extracts any .zip in that directory (a shapefile is .shp + .shx + .dbf
     + .prj that all need to live together; dropping the ZIP and letting this
     script extract it is the easiest path).
  3. Reports OK / missing.

Workflow after cloning:
  1. `python scripts/res_cf/check_external_data.py`   (creates the dirs)
  2. Download:
       - World EEZ v12 from https://www.marineregions.org/downloads.php
       - Natural Earth 110m countries from
         https://www.naturalearthdata.com/downloads/110m-cultural-vectors/
  3. Drop the two ZIPs into the directories created in step 1.
  4. Re-run this script — it will extract and verify.
"""

from pathlib import Path
import shutil
import sys
import zipfile

ROOT = Path(__file__).parent.parent.parent

DATASETS = [
    {
        "label": "EEZ (Marine Regions v12)",
        "dir":   ROOT / "data/shapes/eez",
        "shp":   ROOT / "data/shapes/eez/eez_v12.shp",
        "stem":  "eez_v12",
    },
    {
        "label": "Natural Earth 110m countries",
        "dir":   ROOT / "data/shapes/ne_110m_admin_0_countries",
        "shp":   ROOT / "data/shapes/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
        "stem":  "ne_110m_admin_0_countries",
    },
]


def extract_zips(target_dir: Path, stem: str) -> list[Path]:
    """
    Extract every .zip in target_dir. If the ZIP nests shapefile components
    inside a subfolder, move them up so the .shp ends up directly at
    target_dir/{stem}.shp.
    """
    extracted = []
    for zp in sorted(target_dir.glob("*.zip")):
        with zipfile.ZipFile(zp) as zf:
            zf.extractall(target_dir)
        extracted.append(zp)

    # Flatten: if {stem}.shp is nested, move it + companions up
    expected = target_dir / f"{stem}.shp"
    if not expected.exists():
        found = next(iter(target_dir.rglob(f"{stem}.shp")), None)
        if found is not None and found.parent != target_dir:
            for sibling in list(found.parent.iterdir()):
                if sibling.is_file() and sibling.stem == stem:
                    shutil.move(str(sibling), target_dir / sibling.name)
            try:
                found.parent.rmdir()
            except OSError:
                pass

    return extracted


def main() -> int:
    ok = True
    for ds in DATASETS:
        if not ds["dir"].exists():
            ds["dir"].mkdir(parents=True)
            print(f"  [mkdir] {ds['dir'].relative_to(ROOT)}/")

        if not ds["shp"].exists():
            for zp in extract_zips(ds["dir"], ds["stem"]):
                print(f"  [unzip] {zp.relative_to(ROOT)}")

        if ds["shp"].exists():
            print(f"  [ok]    {ds['label']}")
        else:
            print(f"  [!!]    {ds['label']}")
            print(f"          Missing: {ds['shp'].relative_to(ROOT)}")
            print(f"          Drop the ZIP into {ds['dir'].relative_to(ROOT)}/ and re-run.")
            ok = False

    if not ok:
        return 1
    print("\nAll external data files present.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
