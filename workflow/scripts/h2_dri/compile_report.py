"""
Compile per-scenario summaries into a single project-level report CSV.

Invoked by Snakemake's `script:` directive (compile_report rule).
"""

import sys
from pathlib import Path

import pandas as pd
import pypsa
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from costs import extract_summary

if "snakemake" not in globals():
    raise RuntimeError("compile_report.py must be run via Snakemake")

assumptions_path = Path(snakemake.input.assumptions)
with assumptions_path.open("r", encoding="utf-8") as f:
    assumptions = yaml.safe_load(f)

h2_lhv_kwh_per_kg = assumptions["h2"]["lhv_kwh_per_kg"]
project_name = snakemake.wildcards.project

rows = []
# networks may contain duplicates (collect fans out per tech row); dedupe
# while preserving order so each scenario appears once.
for nc_path in dict.fromkeys(snakemake.input.networks):
    nc_path = Path(nc_path)
    scenario_name = nc_path.stem
    n = pypsa.Network()
    n.import_from_netcdf(nc_path)
    rows.append(extract_summary(n, project_name, scenario_name, h2_lhv_kwh_per_kg))

out_path = Path(snakemake.output[0])
out_path.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out_path, index=False)
