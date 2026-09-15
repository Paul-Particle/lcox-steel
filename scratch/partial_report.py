"""Compile a report from whatever has solved so far, without waiting for the rest.

`compile_report` needs every run in a scenario before it writes anything, which
is right for the real report and useless while a long run is still going. This
does the same per-run extraction over the networks that exist right now and
writes results/partial/report_{scenario}.csv beside a one-page progress summary,
so there is something to read at any point.

Reuses compile_report's own extraction rather than reimplementing it, so a
partial row and a final row are the same row.

Run from the repo root with the project env; takes no arguments.
"""
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pypsa
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "workflow"))
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts"))
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts/viz"))

from common._runs import build_runs_frame, load_scenarios, zone_parents
from viz.compile_report import (
    extract_summary,
    input_variants,
    mark_best_in_country,
    write_report,
)

config = yaml.safe_load((REPO_ROOT / "config/config.yaml").read_text())
assumptions = yaml.safe_load((REPO_ROOT / "config/assumptions.yaml").read_text())
scenarios = load_scenarios(REPO_ROOT / "config/scenarios.csv", config["areas"])
runs = build_runs_frame(scenarios, config["areas"])
parents = zone_parents(config["areas"])

series_dir = REPO_ROOT / "resources/timeseries"
out_dir = REPO_ROOT / "results/partial"
out_dir.mkdir(parents=True, exist_ok=True)

grid_series = {}
for path in series_dir.glob("*_grid_*.parquet"):
    area_tech_variant, start, end = path.stem.rsplit("_", 2)
    grid_series[(area_tech_variant.rsplit("_", 2)[0], start, end)] = pd.read_parquet(path)

destination_area = assumptions["destination"]["area"]
destination_series = {}
for path in series_dir.glob(f"{destination_area}_grid_emissions_*.parquet"):
    _, start, end = path.stem.rsplit("_", 2)
    destination_series[(start, end)] = pd.read_parquet(path)

lines = [f"# Test run progress — {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}", ""]

for scenario, expected in runs.groupby("scenario"):
    solved = sorted((REPO_ROOT / "results" / scenario).glob("*.nc"))
    lines.append(f"## {scenario}: {len(solved)} of {len(expected)} runs solved")
    if not solved:
        lines.append("")
        continue

    rows = []
    failed = []
    for nc_path in solved:
        # {area}_{route}_{start}_{end}: a route never contains an underscore,
        # but an area can (BR_SE), so split from the right.
        area, route, start_date, end_date = nc_path.stem.rsplit("_", 3)
        try:
            n = pypsa.Network()
            n.import_from_netcdf(nc_path)
            run = {"area": area, "route": route,
                   "start_date": start_date, "end_date": end_date}
            run.update(input_variants(scenarios, scenario, run))
            legs = assumptions["transport"]["distance_km"].get(parents.get(area, area))
            summary = extract_summary(
                n, scenario, run, assumptions, legs,
                grid_series.get((area, start_date, end_date)),
                destination_series.get((start_date, end_date)),
            )
            summary["inputs_hash"] = n.meta.get("inputs_hash", "")
            rows.append(summary)
        except Exception as exc:
            failed.append(f"{nc_path.name}: {type(exc).__name__}: {exc}")

    if rows:
        flagged = mark_best_in_country(
            pd.DataFrame(rows), parents, config["report"]["best_zone_by"] or None
        )
        write_report(flagged,
                     out_dir / f"report_{scenario}.csv",
                     out_dir / f".report_{scenario}_diag.csv")

        headline = (
            flagged.loc[:, [c for c in ["area", "route", "lcos_eur_per_t",
                                        "lco_output", "emissions_kg_co2e_per_t_steel"]
                            if c in flagged.columns]]
            .sort_values("area")
        )
        lines.append("")
        lines.append("```")
        lines.append(headline.to_string(index=False, max_rows=400))
        lines.append("```")
    for message in failed:
        lines.append(f"- could not read {message}")
    lines.append("")

(out_dir / "PROGRESS.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines[:6]))
print(f"\nwrote {out_dir}/PROGRESS.md and report_*.csv")
