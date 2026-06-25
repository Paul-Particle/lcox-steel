"""Snakemake entrypoint: assemble inputs, build a PyPSA network, solve it, write netCDF.

The network construction itself lives in `build_network.py` (importable in
isolation). This file is the IO + orchestration shell driven by the
`h2_dri_optimize` rule in `workflow/rules/h2_dri.smk`.
"""

import json
import logging
from pathlib import Path

# pandas 3.0 defaults to ArrowStringArray for strings; xarray (used by PyPSA
# internally) doesn't support it. Set python-native strings before any import
# that touches pandas string data.
import pandas as pd
pd.options.mode.string_storage = "python"

import pyarrow.parquet as pq
import yaml

from build_network import build_network, load_assumptions

from common._logging import configure_logging

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)


def _assemble_multisite_cf(
    cf_paths: list[Path], sites_overlay_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the (site, tech)-keyed CF frame and the sites table for multi-site mode.

    Candidate CF parquets (the `grid-n{N}` variant) carry columns named
    `{tech}@{cell}` and a JSON `site_coords` entry in their Arrow schema metadata
    mapping each column to {lat, lon}. The demand site's coordinates come from the
    sites overlay YAML (`demand_site: {lat, lon}`). Returns (cf_timeseries with a
    (site, tech) MultiIndex on columns, sites DataFrame indexed by site_id with
    x=lon / y=lat). The demand site is added as site_id "plant".
    """
    overlay = yaml.safe_load(Path(sites_overlay_path).read_text()) or {}
    demand = overlay["demand_site"]

    cf_parts: dict[tuple[str, str], pd.Series] = {}
    coords: dict[str, tuple[float, float]] = {}  # site_id -> (lon, lat)
    for p in cf_paths:
        df = pd.read_parquet(p)
        meta = pq.read_schema(p).metadata or {}
        site_coords = json.loads(meta.get(b"site_coords", b"{}"))
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

    coords["plant"] = (float(demand["lon"]), float(demand["lat"]))
    sites = pd.DataFrame(
        {
            "x": {s: c[0] for s, c in coords.items()},
            "y": {s: c[1] for s, c in coords.items()},
        }
    )
    return cf_timeseries, sites


def main() -> None:
    """Assemble rule inputs, build and solve the network, export it to netCDF.

    Classifies each tech-input parquet as a CF series (res_cf) or a grid price
    (entsoe/nem) by its pipeline directory, concatenates CF columns into one
    multi-tech frame, aligns the optional price series, merges base+overlay
    assumptions, then builds, optimises (HiGHS), and writes the network.
    """
    project  = snakemake.wildcards.project
    scenario = snakemake.wildcards.scenario
    out_path = Path(snakemake.output.network)

    # Each tech input is one parquet, classified by its pipeline directory:
    # resources/res_cf/... → capacity-factor series (column names are tech keys);
    # resources/entsoe/... or resources/nem/... → grid price series.
    cf_paths:   list[Path] = []
    grid_paths: list[Path] = []
    for raw in snakemake.input.tech_inputs:
        p = Path(raw)
        (cf_paths if p.parent.name == "res_cf" else grid_paths).append(p)

    if len(grid_paths) > 1:
        raise ValueError(f"{project}/{scenario}: multiple grid inputs: {grid_paths}")
    grid_path = grid_paths[0] if grid_paths else None
    price_series = pd.read_parquet(grid_path).iloc[:, 0] if grid_path is not None else None

    # A config/sites_{project}_{scenario}.yaml overlay (0 or 1 paths via optional())
    # switches the scenario into multi-site mode: one electricity bus per candidate
    # site, distance-costed HVDC links to the demand site. Absent ⇒ single-bus.
    site_overlays = list(snakemake.input.sites_overlay)
    sites_overlay_path = Path(site_overlays[0]) if site_overlays else None
    multisite = sites_overlay_path is not None

    sites = None
    demand_site = None
    if cf_paths and multisite:
        cf_timeseries, sites = _assemble_multisite_cf(cf_paths, sites_overlay_path)
        demand_site = "plant"
    elif cf_paths:
        cf_parts: dict[str, pd.Series] = {}
        for p in cf_paths:
            df = pd.read_parquet(p)
            # Single- and multi-column parquets are uniform: columns are tech keys.
            # 03_build_cf_timeseries.py names single-column outputs by the tech wildcard.
            for col in df.columns:
                if col in cf_parts:
                    raise ValueError(
                        f"{project}/{scenario}: duplicate tech key '{col}' "
                        f"across CF inputs"
                    )
                cf_parts[col] = df[col]
        cf_timeseries = pd.DataFrame(cf_parts)
    elif price_series is not None:
        cf_timeseries = pd.DataFrame(index=price_series.index)
    else:
        raise ValueError(f"{project}/{scenario}: no CF techs and no grid input")

    if price_series is not None:
        price_series = price_series.reindex(cf_timeseries.index)
        if price_series.isna().any():
            raise ValueError(
                f"Price series has {price_series.isna().sum()} missing values after "
                "aligning to CF index. Check that the grid file covers the same period."
            )

    # optional() yields a Namedlist of 0 or 1 paths: present iff a
    # config/assumptions_{project}_{scenario}.yaml file exists on disk.
    overlays = list(snakemake.input.assumptions_overlay)
    overlay_path = Path(overlays[0]) if overlays else None
    base_path = Path(snakemake.input.assumptions_base)

    assumptions = load_assumptions(base_path, overlay_path)
    overlay_name = overlay_path.name if overlay_path else "none"
    mode = f"multi-site ({len(sites)} sites)" if sites is not None else "single-site"
    log.info(
        f"building {mode} network for project={project} scenario={scenario} "
        f"techs={list(cf_timeseries.columns)} (overlay={overlay_name})"
    )
    n = build_network(
        assumptions, cf_timeseries, price_series, sites=sites, demand_site=demand_site
    )
    log.info(f"optimising with HiGHS (snapshots={len(n.snapshots)})")
    n.optimize(solver_name="highs")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(out_path)
    log.info(f"network saved to {out_path}")


if __name__ == "__main__":
    main()
