"""
Extract a shapefile ZIP into its target directory.

Driven by Snakemake: input[0] is the ZIP, output.shp is the canonical .shp path
inside the same directory. If the ZIP nests the shapefile components under a
subfolder (the marineregions EEZ ZIP does this), flatten them up so the
components land directly next to the ZIP — matching what the downstream rules
declare as input.
"""

from pathlib import Path
import shutil
import zipfile

zip_path = Path(snakemake.input[0])
shp_path = Path(snakemake.output.shp)
target_dir = shp_path.parent
stem = shp_path.stem

with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(target_dir)

if not shp_path.exists():
    nested = next(iter(target_dir.rglob(f"{stem}.shp")), None)
    if nested is None:
        raise FileNotFoundError(
            f"{zip_path} did not contain {stem}.shp (or any file with that stem)."
        )
    for sibling in list(nested.parent.iterdir()):
        if sibling.is_file() and sibling.stem == stem:
            shutil.move(str(sibling), target_dir / sibling.name)
    try:
        nested.parent.rmdir()
    except OSError:
        pass
