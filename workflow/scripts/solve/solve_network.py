"""Snakemake entrypoint: assemble inputs, build a PyPSA network, solve it, write netCDF.

The network construction itself lives in `build_network.py` (importable in
isolation). This file is the IO + orchestration shell driven by the
`solve_network` rule in `workflow/rules/solve.smk`.
"""

import json
import logging
import os
from pathlib import Path

# pandas 3.0 defaults to ArrowStringArray for strings; xarray (used by PyPSA
# internally) doesn't support it. Set python-native strings before any import
# that touches pandas string data.
import pandas as pd
pd.options.mode.string_storage = "python"

import geopandas as gpd
import pyarrow.parquet as pq

from build_network import build_network, load_assumptions

from common._logging import configure_logging
from common._provenance import input_manifest

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)


def _area_representative_point(regions_path: Path) -> tuple[float, float]:
    """(lon, lat) of a point guaranteed to lie inside the area's own geometry.

    The plant sits in the middle of the area it is built in. representative_point
    rather than centroid: a centroid can fall outside a concave or multi-part
    country.
    """
    geometry = gpd.read_parquet(regions_path).union_all()
    point = geometry.representative_point()
    return float(point.x), float(point.y)


def _assemble_multisite_cf(
    cf_paths: list[Path], plant_xy: tuple[float, float]
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Build the (site, tech)-keyed CF frame, sites table and demand-site id.

    Candidate CF parquets (the `multi-n{N}` and `anchor-colo-n{N}` variants)
    carry columns named `{tech}@{cell}` and a JSON `site_coords` entry in their
    Arrow schema metadata mapping each column to {lat, lon}. The demand site
    comes from one of two places:

    - a `demand_site` Arrow-metadata key naming one of the file's own columns
      (the anchor-colo case): that column's site doubles as the plant bus, so
      the anchor generator sits directly on the demand bus (its HVDC link is
      skipped by build_network);
    - otherwise `plant_xy`, the area's own representative point, added as
      site_id "plant" (the multi-n* sweep case).

    Returns (cf_timeseries with a (site, tech) MultiIndex on columns, sites
    DataFrame indexed by site_id with x=lon / y=lat, demand_site id).
    """
    demand_site = None

    cf_parts: dict[tuple[str, str], pd.Series] = {}
    coords: dict[str, tuple[float, float]] = {}  # site_id -> (lon, lat)
    for p in cf_paths:
        df = pd.read_parquet(p)
        meta = pq.read_schema(p).metadata or {}
        site_coords = json.loads(meta.get(b"site_coords", b"{}"))
        demand_col = meta.get(b"demand_site", b"").decode() or None
        if demand_col is not None:
            if demand_col not in df.columns:
                raise ValueError(
                    f"demand_site metadata '{demand_col}' in {p} names no column"
                )
            tech, cell = demand_col.split("@", 1)
            demand_from_file = f"{tech}-{cell}"
            if demand_site is not None and demand_site != demand_from_file:
                raise ValueError(
                    f"conflicting demand sites across CF inputs: "
                    f"'{demand_site}' vs '{demand_from_file}' (from {p})"
                )
            demand_site = demand_from_file
        for col in df.columns:
            if "@" not in col:
                raise ValueError(
                    f"multi-site CF column '{col}' in {p} is missing the "
                    "'{tech}@{cell}' key — multi-site scenarios need multi-n* inputs"
                )
            tech, cell = col.split("@", 1)
            # A cell id (c00, c01, …) is only unique within a tech's parquet —
            # solar@c00 and wind-onshore@c00 are different physical cells — so the
            # global site id folds the tech in.
            site = f"{tech}-{cell}"
            key = (site, tech)
            if key in cf_parts:
                raise ValueError(f"duplicate (site, tech) key {key} across CF inputs")
            cc = site_coords.get(col)
            if cc is None:
                raise ValueError(f"no coords in metadata for column '{col}' in {p}")
            cf_parts[key] = df[col]
            coords[site] = (cc["lon"], cc["lat"])

    cf_timeseries = pd.DataFrame(cf_parts)
    cf_timeseries.columns = pd.MultiIndex.from_tuples(cf_timeseries.columns)

    if demand_site is None:
        demand_site = "plant"
        coords["plant"] = plant_xy
    sites = pd.DataFrame(
        {
            "x": {s: c[0] for s, c in coords.items()},
            "y": {s: c[1] for s, c in coords.items()},
        }
    )
    return cf_timeseries, sites, demand_site


def main() -> None:
    """Assemble rule inputs, build and solve the network, export it to netCDF.

    Capacity-factor and grid-price inputs arrive already separated by the rule
    (the CSV's `tech` column does the classifying), so this concatenates the CF
    columns into one multi-tech frame, aligns the optional price series, merges
    base+overlay assumptions, then builds, optimises (HiGHS), and writes the
    network for the one route named by the wildcard.
    """
    scenario = snakemake.wildcards.scenario
    area     = snakemake.wildcards.area
    route    = snakemake.wildcards.route
    run      = f"{scenario}/{area}/{route}"
    out_path = Path(snakemake.output.network)

    cf_paths   = [Path(raw) for raw in snakemake.input.cf_inputs]
    grid_paths = [Path(raw) for raw in snakemake.input.grid_input]

    if len(grid_paths) > 1:
        raise ValueError(f"{run}: multiple grid inputs: {grid_paths}")
    grid_path = grid_paths[0] if grid_paths else None
    price_series = pd.read_parquet(grid_path).iloc[:, 0] if grid_path is not None else None

    # Multi-site mode (one electricity bus per candidate site, distance-costed
    # HVDC links to the demand site) is triggered by the CF data itself: the
    # multi-n* and anchor-colo-n* variants name their columns '{tech}@{cell}',
    # one column per candidate cell. A single-cell input has no '@' and stays
    # on one bus.
    multisite = any(
        "@" in col for p in cf_paths for col in pq.read_schema(p).names
    )

    sites = None
    demand_site = None
    if cf_paths and multisite:
        cf_timeseries, sites, demand_site = _assemble_multisite_cf(
            cf_paths, _area_representative_point(Path(snakemake.input.regions))
        )
    elif cf_paths:
        cf_parts: dict[str, pd.Series] = {}
        for p in cf_paths:
            df = pd.read_parquet(p)
            # Single- and multi-column parquets are uniform: columns are tech keys.
            # d1_area_average.py names single-column outputs by the tech wildcard.
            for col in df.columns:
                if col in cf_parts:
                    raise ValueError(
                        f"{run}: duplicate tech key '{col}' "
                        f"across CF inputs"
                    )
                cf_parts[col] = df[col]
        cf_timeseries = pd.DataFrame(cf_parts)
    elif price_series is not None:
        cf_timeseries = pd.DataFrame(index=price_series.index)
    else:
        raise ValueError(f"{run}: no CF techs and no grid input")

    if price_series is not None:
        price_series = price_series.reindex(cf_timeseries.index)
        if price_series.isna().any():
            raise ValueError(
                f"Price series has {price_series.isna().sum()} missing values after "
                "aligning to CF index. Check that the grid file covers the same period."
            )

    # optional() yields a Namedlist of 0 or 1 paths: present iff a
    # config/assumptions_{scenario}.yaml file exists on disk.
    overlays = list(snakemake.input.assumptions_overlay)
    overlay_path = Path(overlays[0]) if overlays else None
    base_path = Path(snakemake.input.assumptions_base)

    assumptions = load_assumptions(base_path, overlay_path)
    overlay_name = overlay_path.name if overlay_path else "none"
    mode = f"multi-site ({len(sites)} sites)" if sites is not None else "single-site"
    log.info(
        f"building {mode} network for scenario={scenario} area={area} route={route} "
        f"techs={list(cf_timeseries.columns)} (overlay={overlay_name})"
    )
    n = build_network(
        route, assumptions, cf_timeseries, price_series,
        sites=sites, demand_site=demand_site,
    )
    log.info(f"optimising with HiGHS (snapshots={len(n.snapshots)})")
    # HiGHS parallelises across cores by default (not governed by OMP_NUM_THREADS).
    # When running many independent solves in parallel (batch runs), pin each to a
    # single thread via HIGHS_THREADS=1 so N solves use N cores instead of each
    # grabbing all of them and thrashing. Unset => HiGHS default (all cores), which
    # is best for a single interactive solve.
    solver_options = {}
    highs_threads = os.environ.get("HIGHS_THREADS")
    if highs_threads:
        solver_options["threads"] = int(highs_threads)
    # HIGHS_SOLVER selects the algorithm: "ipm" (interior-point / barrier) converges
    # in a small, degeneracy-insensitive iteration count and is far faster than the
    # default simplex on the large degenerate grid+p95 LPs. We only read off the
    # objective and capacities (no basis needed), so crossover to a vertex is skipped
    # unless HIGHS_CROSSOVER=on. Unset => HiGHS default (simplex), best for easy LPs.
    highs_solver = os.environ.get("HIGHS_SOLVER")
    if highs_solver:
        solver_options["solver"] = highs_solver
        if highs_solver == "ipm":
            solver_options["run_crossover"] = os.environ.get("HIGHS_CROSSOVER", "off")
    n.optimize(solver_name="highs", solver_options=solver_options)

    # Every file the rule declared, fingerprinted into the network so the result
    # stays tied to what produced it. compile_report lifts `inputs_hash` into a
    # report column; the per-file map stays here.
    n.meta.update(input_manifest([Path(raw) for raw in snakemake.input]))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(out_path)
    log.info(f"network saved to {out_path} (inputs_hash={n.meta['inputs_hash']})")


if __name__ == "__main__":
    main()
