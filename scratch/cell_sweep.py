"""Per-cell colocated sweep: emit Snakemake inputs, collect results, clean up.

Benchmarks "everything colocated in one cell, brute-forced over every eligible
cell" against the anchor-colocation methodology. Snakemake does the solving —
this script writes the scenario rows and per-cell CF files it needs, then reads
the solved networks back into one table.

Not a Snakemake script: it runs standalone, from the repo root.

    python scratch/cell_sweep.py emit --area VIC1 --stride 1
    snakemake --jobs 6 -- $(python scratch/cell_sweep.py targets --area VIC1)
    python scratch/cell_sweep.py collect --area VIC1
    python scratch/cell_sweep.py compare --area VIC1
    python scratch/cell_sweep.py clean
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow"))

from common._paths import CUTOUTS, REPO_ROOT, RESULTS, SHAPES_RES, TIMESERIES
from common._runs import build_runs_frame, load_scenarios
from scripts.res_cf._helpers_res_cf import (
    eligibility_mask_2d,
    eligibility_weights,
    pick_p95_cell,
)

# Both onshore techs sit in the same cell — that is the whole point of the
# "everything colocated" arm. Offshore wind is excluded: it cannot share a
# land cell, so a colocated triple has no meaning for it.
SWEEP_TECHS = ("wind-onshore", "solar")

# The scenario table is a tracked config file, so `emit` moves the original
# aside under this suffix and `clean` puts it back.
SCENARIO_TABLE = REPO_ROOT / "config" / "scenarios.csv"
SCENARIO_TABLE_BACKUP = SCENARIO_TABLE.with_suffix(".csv.presweep")

SWEEP_DIR = RESULTS / "cell_sweep"

# `-isl` rather than `-islanded` keeps the scenario name short; the mode is
# recovered from the suffix when the results are collected.
MODE_SUFFIXES = {"islanded": "isl", "grid": "grid"}

# The pre-selection arms the sweep is benchmarked against. These are ordinary
# pipeline runs named in config/scenarios.csv, not emitted here: one solved
# network per (mode, route) each. `compare` reads them alongside the swept
# cells so all arms are scored by the same `extract_summary`.
REFERENCE_ARMS = {
    "colo-ref": "anchor-colo-n3 (multi-site)",
    "p95-ref": "bestsite-p95 (single site)",
}

log = logging.getLogger("cell_sweep")


def cell_tag(y: int, x: int) -> str:
    """The variant name for one cutout cell, e.g. `cell-y008x012`."""
    return f"cell-y{y:03d}x{x:03d}"


def load_area_config(area: str) -> tuple[dict, str]:
    """The `res_cf` config block and the region tag the area's geometry is filed under."""
    config = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
    return config, config["areas"][area]["region"]


def emit_sweep_inputs(
    area: str, start_date: str, end_date: str, routes: str, modes: list[str],
    stride: int, limit: int,
) -> None:
    """Write one CF file per eligible cell and append the scenario rows for them."""
    if SCENARIO_TABLE_BACKUP.exists():
        raise FileExistsError(
            f"{SCENARIO_TABLE_BACKUP.name} already exists — an earlier emit was never "
            "cleaned up. Run `clean` first (it restores the scenario table)."
        )

    config, region_tag = load_area_config(area)
    res_cf_config = config["res_cf"]

    cutout_path = CUTOUTS / f"{area}_{start_date}_{end_date}.nc"
    if not cutout_path.exists():
        raise FileNotFoundError(
            f"{cutout_path.relative_to(REPO_ROOT)} is missing. Build it first:\n"
            f"  snakemake --cores 1 -- {cutout_path.relative_to(REPO_ROOT)}\n"
            "At full ERA5 resolution, move any *_backup.nc sibling aside first — the "
            "retrieve rule prefers a backup over CDS and the shipped one is 0.5 degree."
        )
    cutout = atlite.Cutout(str(cutout_path))

    regions = gpd.read_parquet(SHAPES_RES / f"{area}_geo.parquet").to_crs(4326)
    land_geometry = regions.loc[regions["region"] == region_tag].geometry.iloc[0]

    # The same eligibility seam d2_bestsite_p95 and d3_anchor_colo select on, so
    # both arms of the benchmark drop majority-sea coastal cells identically (#41).
    eligible_mask = eligibility_mask_2d(
        cutout,
        land_geometry,
        float(res_cf_config["min_land_fraction"]),
        res_cf_config["eligibility_source"],
    )
    log.info(
        f"{area} {start_date}-{end_date}: {int(eligible_mask.sum())} eligible cells "
        f"of {eligible_mask.size} in the cutout grid"
    )

    log.info("computing per-cell CF grids (both techs, whole grid, one pass)")
    wind_cf_config = res_cf_config["wind_cf"]
    cf_grids = {
        "wind-onshore": cutout.wind(
            turbine=res_cf_config["wind_onshore_turbine"],
            capacity_factor_timeseries=True,
            smooth=wind_cf_config["smooth"],
            add_cutout_windspeed=wind_cf_config["add_cutout_windspeed"],
        ).transpose("time", "y", "x"),
        "solar": cutout.pv(
            panel=res_cf_config["pv_panel"],
            orientation=res_cf_config["pv_orientation"],
            capacity_factor_timeseries=True,
        ).transpose("time", "y", "x"),
    }
    timestamps = pd.DatetimeIndex(cf_grids["solar"].coords["time"].values, name="time")
    latitudes = cutout.data.coords["y"].values
    longitudes = cutout.data.coords["x"].values

    # Which eligible cells actually get solved. `stride` thins the lattice to
    # fit a run budget, but a regular lattice samples the *maximum* badly: in
    # Spain the best-resource cell sits in the Galician northwest, and a stride
    # of 3 picks a best cell ~800 km away and ~14% worse on blended CF. That
    # bias flatters the pre-selection method under test, so three cells are
    # always solved on top of the lattice — the best cell by blended CF, and
    # each tech's P95 anchor, which is the cell d3_anchor_colo fixes its
    # co-location search on. `sampling` in the manifest records why each cell
    # is in, so the unbiased lattice subset stays separable from the forced
    # ones when the results are read.
    lattice_mask = np.zeros_like(eligible_mask)
    lattice_mask[::stride, ::stride] = True
    lattice_cells = eligible_mask & lattice_mask
    sampling_reasons = {
        (int(y), int(x)): ["lattice"] for y, x in zip(*np.where(lattice_cells))
    }

    blended_cf_mean = np.where(
        eligible_mask,
        0.5 * cf_grids["wind-onshore"].values.mean(axis=0)
        + 0.5 * cf_grids["solar"].values.mean(axis=0),
        -np.inf,
    )
    best_y, best_x = np.unravel_index(int(np.argmax(blended_cf_mean)), eligible_mask.shape)
    forced_cells = {(int(best_y), int(best_x)): "proxy-best"}
    anchor_weights = eligibility_weights(
        cutout,
        land_geometry,
        float(res_cf_config["min_land_fraction"]),
        res_cf_config["eligibility_source"],
    )
    for tech, cf_grid in cf_grids.items():
        anchor_y, anchor_x = pick_p95_cell(cf_grid.mean("time"), anchor_weights)
        forced_cells.setdefault((int(anchor_y), int(anchor_x)), f"p95-anchor-{tech}")
    for cell, reason in forced_cells.items():
        sampling_reasons.setdefault(cell, []).append(reason)

    eligible_cells = sorted(sampling_reasons)
    if limit:
        eligible_cells = eligible_cells[:limit]
        sampling_reasons = {cell: sampling_reasons[cell] for cell in eligible_cells}
    log.info(
        f"stride {stride}: {int(lattice_cells.sum())} lattice cells + "
        f"{len(forced_cells)} forced diagnostic cells "
        f"(overlaps merged) = {len(eligible_cells)} to solve"
    )

    TIMESERIES.mkdir(parents=True, exist_ok=True)
    scenario_rows = []
    manifest_records = []
    for y, x in eligible_cells:
        variant = cell_tag(y, x)
        cell_means = {}
        for tech, grid in cf_grids.items():
            series = pd.Series(grid.values[:, y, x], index=timestamps, name=tech)
            cell_means[tech] = float(series.mean())
            series.to_frame().to_parquet(
                TIMESERIES / f"{area}_{tech}_{variant}_{start_date}_{end_date}.parquet",
                index=True,
            )

        for mode in modes:
            scenario = f"{variant}-{MODE_SUFFIXES[mode]}"
            for tech in SWEEP_TECHS:
                scenario_rows.append(
                    (scenario, routes, tech, variant, area, start_date, end_date)
                )
            # Islanded and grid-connected differ by this row alone — no
            # assumptions overlay is involved, so none is written.
            if mode == "grid":
                scenario_rows.append(
                    (scenario, routes, "grid", "dayahead", area, start_date, end_date)
                )

        manifest_records.append(
            {
                "cell_tag": variant,
                "y": y,
                "x": x,
                "latitude": float(latitudes[y]),
                "longitude": float(longitudes[x]),
                "wind_onshore_cf_mean": cell_means["wind-onshore"],
                "solar_cf_mean": cell_means["solar"],
                "sampling": "+".join(sampling_reasons[(y, x)]),
            }
        )

    shutil.copyfile(SCENARIO_TABLE, SCENARIO_TABLE_BACKUP)
    emitted_table = pd.DataFrame(
        scenario_rows,
        columns=["scenario", "route", "tech", "variant", "area", "start_date", "end_date"],
    )
    with open(SCENARIO_TABLE, "a") as table_file:
        table_file.write(f"# --- cell_sweep {area} {start_date}-{end_date} ---\n")
        emitted_table.to_csv(table_file, header=False, index=False)

    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(manifest_records)
    manifest.to_csv(SWEEP_DIR / "manifest.csv", index=False)

    # A row naming an area that the table cannot resolve yields a *runnable*
    # solve job with no CF inputs rather than an error, so the row count per
    # scenario is checked here instead of being discovered 20 hours later.
    validated = load_scenarios(SCENARIO_TABLE, config["areas"])
    expected_rows = {
        "islanded": len(SWEEP_TECHS),
        "grid": len(SWEEP_TECHS) + 1,
    }
    actual_rows = validated["scenario"].value_counts()
    for mode in modes:
        for y, x in eligible_cells:
            scenario = f"{cell_tag(y, x)}-{MODE_SUFFIXES[mode]}"
            if actual_rows.get(scenario, 0) != expected_rows[mode]:
                raise ValueError(
                    f"scenario {scenario} resolved to {actual_rows.get(scenario, 0)} "
                    f"rows, expected {expected_rows[mode]} — the solve would run with "
                    "the wrong inputs"
                )

    all_runs = build_runs_frame(validated, config["areas"])
    run_count = int(all_runs["scenario"].str.startswith("cell-").sum())
    log.info(
        f"emitted {len(eligible_cells)} cells x {len(modes)} modes -> "
        f"{len(scenario_rows)} table rows, {run_count} networks to solve"
    )
    log.info(f"manifest: {(SWEEP_DIR / 'manifest.csv').relative_to(REPO_ROOT)}")
    log.info(f"scenario table backed up to {SCENARIO_TABLE_BACKUP.name}")


def print_sweep_targets(area: str, start_date: str, end_date: str) -> None:
    """Print every network path the emitted rows expect, for a snakemake command line."""
    config, _ = load_area_config(area)
    runs = build_runs_frame(load_scenarios(SCENARIO_TABLE, config["areas"]), config["areas"])
    sweep_runs = runs.loc[
        runs["scenario"].str.startswith("cell-")
        & (runs["area"] == area)
        & (runs["start_date"] == start_date)
        & (runs["end_date"] == end_date)
    ]
    for row in sweep_runs.itertuples():
        print(
            f"results/{row.scenario}/{row.area}_{row.route}_{row.start_date}_{row.end_date}.nc"
        )


def collect_sweep_results(area: str, start_date: str, end_date: str) -> None:
    """Read every solved sweep network into one table using the report's own metrics."""
    # Imported here rather than at module scope: compile_report is a Snakemake
    # script, so importing it pulls the stub shim in and reconfigures logging.
    sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts" / "viz"))
    import compile_report

    config, _ = load_area_config(area)
    scenarios = load_scenarios(SCENARIO_TABLE, config["areas"])
    network_paths = sorted(RESULTS.glob(f"cell-*/{area}_*_{start_date}_{end_date}.nc"))
    if not network_paths:
        raise FileNotFoundError(
            f"no solved networks under {RESULTS.relative_to(REPO_ROOT)}/cell-*/ "
            f"for {area} {start_date}-{end_date}"
        )
    log.info(f"reading {len(network_paths)} solved networks")

    summary_rows = []
    for nc_path in network_paths:
        scenario = nc_path.parent.name
        network_area, route, network_start, network_end = nc_path.stem.split("_")
        network = pypsa.Network()
        network.import_from_netcdf(nc_path)
        run = {
            "area": network_area,
            "route": route,
            "start_date": network_start,
            "end_date": network_end,
        }
        run.update(compile_report.input_variants(scenarios, scenario, run))
        summary = compile_report.extract_summary(network, scenario, run)
        summary["cell_tag"] = scenario.rsplit("-", 1)[0]
        summary["mode"] = scenario.rsplit("-", 1)[1]
        summary_rows.append(summary)

    manifest = pd.read_csv(SWEEP_DIR / "manifest.csv")
    sweep_summary = (
        pd.DataFrame(summary_rows)
        .merge(manifest, on="cell_tag", how="left", validate="many_to_one")
        .sort_values(["mode", "route", "lcos_eur_per_t"])
        .reset_index(drop=True)
    )
    summary_path = SWEEP_DIR / f"sweep_{area}_{start_date}_{end_date}.csv"
    sweep_summary.to_csv(summary_path, index=False)
    log.info(f"wrote {len(sweep_summary)} rows to {summary_path.relative_to(REPO_ROOT)}")

    best_per_arm = (
        sweep_summary
        .loc[sweep_summary.groupby(["mode", "route"])["lcos_eur_per_t"].idxmin()]
        .loc[:, ["mode", "route", "cell_tag", "latitude", "longitude",
                 "wind_onshore_cf_mean", "solar_cf_mean", "lcos_eur_per_t"]]
    )
    log.info(f"cheapest cell per (mode, route):\n{best_per_arm.to_string(index=False)}")


def compare_arms(area: str, start_date: str, end_date: str) -> None:
    """Score the pre-selection arms against the swept per-cell distribution."""
    sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts" / "viz"))
    import compile_report

    summary_path = SWEEP_DIR / f"sweep_{area}_{start_date}_{end_date}.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"{summary_path.relative_to(REPO_ROOT)} is missing — run `collect` first."
        )
    swept = pd.read_csv(summary_path)

    config, _ = load_area_config(area)
    scenarios = load_scenarios(SCENARIO_TABLE, config["areas"])

    reference_rows = []
    for prefix, arm_label in REFERENCE_ARMS.items():
        for mode, suffix in MODE_SUFFIXES.items():
            scenario = f"{prefix}-{suffix}"
            for nc_path in sorted(
                RESULTS.glob(f"{scenario}/{area}_*_{start_date}_{end_date}.nc")
            ):
                network_area, route, network_start, network_end = nc_path.stem.split("_")
                network = pypsa.Network()
                network.import_from_netcdf(nc_path)
                run = {
                    "area": network_area,
                    "route": route,
                    "start_date": network_start,
                    "end_date": network_end,
                }
                run.update(compile_report.input_variants(scenarios, scenario, run))
                summary = compile_report.extract_summary(network, scenario, run)
                summary["arm"] = arm_label
                summary["mode"] = mode
                reference_rows.append(summary)
    if not reference_rows:
        raise FileNotFoundError(
            f"no reference-arm networks under {RESULTS.relative_to(REPO_ROOT)}/"
            f"{{{','.join(REFERENCE_ARMS)}}}-* for {area} {start_date}-{end_date}"
        )
    reference = pd.DataFrame(reference_rows)
    reference.to_csv(
        SWEEP_DIR / f"reference_{area}_{start_date}_{end_date}.csv", index=False
    )

    # The lattice-only best is the honest like-for-like number when the sweep was
    # strided: the forced diagnostic cells were chosen using resource knowledge,
    # which is exactly what the pre-selection arms are allowed to use too.
    comparison_rows = []
    for (mode, route), cells in swept.groupby(["mode", "route"]):
        mode_name = {v: k for k, v in MODE_SUFFIXES.items()}[mode]
        lattice = cells.loc[cells["sampling"] == "lattice", "lcos_eur_per_t"]
        best_overall = cells["lcos_eur_per_t"].min()
        best_cell = cells.loc[cells["lcos_eur_per_t"].idxmin()]
        for _, arm in reference.loc[reference["route"] == route].iterrows():
            if arm["mode"] != mode_name:
                continue
            arm_lcos = arm["lcos_eur_per_t"]
            comparison_rows.append(
                {
                    "mode": mode_name,
                    "route": route,
                    "arm": arm["arm"],
                    "arm_lcos_eur_per_t": arm_lcos,
                    "sweep_best_eur_per_t": best_overall,
                    "sweep_best_cell": best_cell["cell_tag"],
                    "sweep_best_lat": best_cell["latitude"],
                    "sweep_best_lon": best_cell["longitude"],
                    "lattice_best_eur_per_t": lattice.min(),
                    "lattice_median_eur_per_t": lattice.median(),
                    "arm_premium_vs_sweep_pct": 100.0 * (arm_lcos - best_overall) / best_overall,
                    "arm_premium_vs_lattice_pct": 100.0 * (arm_lcos - lattice.min()) / lattice.min(),
                    "arm_percentile_in_lattice": 100.0 * (lattice < arm_lcos).mean(),
                    "cells_swept": len(cells),
                }
            )

    comparison = (
        pd.DataFrame(comparison_rows)
        .sort_values(["mode", "route", "arm"])
        .reset_index(drop=True)
    )
    comparison_path = SWEEP_DIR / f"comparison_{area}_{start_date}_{end_date}.csv"
    comparison.to_csv(comparison_path, index=False)
    log.info(f"wrote {comparison_path.relative_to(REPO_ROOT)}")
    log.info(
        "pre-selection premium over the brute-force optimum:\n"
        + comparison.loc[
            :,
            ["mode", "route", "arm", "arm_lcos_eur_per_t", "sweep_best_eur_per_t",
             "arm_premium_vs_sweep_pct", "arm_percentile_in_lattice"],
        ].to_string(index=False)
    )


def clean_sweep_artifacts(area: str, start_date: str, end_date: str, purge: bool) -> None:
    """Restore the scenario table and delete the per-cell inputs and networks."""
    if SCENARIO_TABLE_BACKUP.exists():
        shutil.move(SCENARIO_TABLE_BACKUP, SCENARIO_TABLE)
        log.info(f"restored {SCENARIO_TABLE.name} from {SCENARIO_TABLE_BACKUP.name}")
    else:
        log.warning(f"no {SCENARIO_TABLE_BACKUP.name} to restore — leaving the table alone")

    cf_files = list(TIMESERIES.glob(f"{area}_*_cell-*_{start_date}_{end_date}.parquet"))
    for cf_file in cf_files:
        cf_file.unlink()
    log.info(f"deleted {len(cf_files)} per-cell CF files")

    result_dirs = [path for path in RESULTS.glob("cell-*") if path.is_dir()]
    for result_dir in result_dirs:
        shutil.rmtree(result_dir)
    log.info(f"deleted {len(result_dirs)} per-cell result directories")

    if purge and SWEEP_DIR.exists():
        shutil.rmtree(SWEEP_DIR)
        log.info(f"purged {SWEEP_DIR.relative_to(REPO_ROOT)} (manifest and summary)")


def main() -> None:
    """Parse the subcommand and run it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=["emit", "targets", "collect", "compare", "clean"]
    )
    parser.add_argument("--area", default="VIC1")
    parser.add_argument("--start-date", default="20250101")
    parser.add_argument("--end-date", default="20251231")
    parser.add_argument(
        "--routes",
        default="h2-dri-eaf|moe-eaf",
        help="route cell for the emitted rows; '|' separates routes",
    )
    parser.add_argument(
        "--modes",
        default="islanded,grid",
        help="comma-separated subset of islanded,grid",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="solve every Nth eligible cell in y and x; 1 sweeps all of them",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="truncate to the first N selected cells, for a smoke run",
    )
    parser.add_argument(
        "--purge",
        action="store_true",
        help="clean: also delete the manifest and collected summary",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s"
    )
    modes = [mode.strip() for mode in args.modes.split(",") if mode.strip()]
    unknown_modes = set(modes) - set(MODE_SUFFIXES)
    if unknown_modes:
        raise ValueError(f"unknown modes {sorted(unknown_modes)}")

    if args.command == "emit":
        if args.stride < 1:
            raise ValueError(f"--stride must be >= 1, got {args.stride}")
        emit_sweep_inputs(
            args.area, args.start_date, args.end_date, args.routes, modes,
            args.stride, args.limit,
        )
    elif args.command == "targets":
        print_sweep_targets(args.area, args.start_date, args.end_date)
    elif args.command == "collect":
        collect_sweep_results(args.area, args.start_date, args.end_date)
    elif args.command == "compare":
        compare_arms(args.area, args.start_date, args.end_date)
    else:
        clean_sweep_artifacts(args.area, args.start_date, args.end_date, args.purge)


if __name__ == "__main__":
    main()
