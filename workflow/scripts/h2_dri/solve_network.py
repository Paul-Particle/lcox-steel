"""
Snakemake entrypoint: assemble inputs, build a PyPSA network, solve it,
write the result to netcdf.

The network construction itself lives in `build_network.py` (importable
in isolation). This file is the IO + orchestration shell driven by the
`h2_dri_optimize` rule in `workflow/rules/h2_dri.smk`.
"""

import logging
from pathlib import Path

# pandas 3.0 defaults to ArrowStringArray for strings; xarray (used by PyPSA
# internally) doesn't support it. Set python-native strings before any import
# that touches pandas string data.
import pandas as pd
pd.options.mode.string_storage = "python"

from build_network import build_network, load_assumptions

from common._logging import configure_logging

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)


def main() -> None:
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

    if cf_paths:
        cf_parts: dict[str, pd.Series] = {}
        for p in cf_paths:
            df = pd.read_parquet(p)
            # Single- and multi-column parquets are uniform: columns are tech keys.
            # build_res_cf_profile.py names single-column outputs by the tech wildcard.
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
    log.info(
        f"building network for project={project} scenario={scenario} "
        f"techs={list(cf_timeseries.columns)} (overlay={overlay_name})"
    )
    n = build_network(assumptions, cf_timeseries, price_series)
    log.info(f"optimising with HiGHS (snapshots={len(n.snapshots)})")
    n.optimize(solver_name="highs")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(out_path)
    log.info(f"network saved to {out_path}")


if __name__ == "__main__":
    main()
