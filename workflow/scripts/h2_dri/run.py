"""
PyPSA DRI-hydrogen optimisation entry point.

Invoked by Snakemake's `script:` directive; one rule fires per
(project, scenario) wildcard pair.
"""

import sys
from pathlib import Path

# pandas 3.0 defaults to ArrowStringArray for strings; xarray (used by PyPSA
# internally) doesn't support it. Set python-native strings before any import
# that touches pandas string data.
import pandas as pd
pd.options.mode.string_storage = "python"

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from network import build_network


def _read_series(path: Path) -> pd.Series:
    """Load a single-column parquet with a time index; return a tz-naive UTC Series."""
    df = pd.read_parquet(path)
    if "time" in df.columns:
        df = df.set_index("time")
    df.index.name = None
    s = df.iloc[:, 0]
    s.index = pd.to_datetime(s.index)
    if s.index.tz is not None:
        s.index = s.index.tz_convert("UTC").tz_localize(None)
    return s


def main() -> None:
    project  = snakemake.wildcards.project
    scenario = snakemake.wildcards.scenario
    techs    = list(snakemake.params.techs)
    out_path = Path(snakemake.output.network)

    with Path(snakemake.input.assumptions).open(encoding="utf-8") as f:
        assumptions = yaml.safe_load(f)

    tech_files  = dict(zip(techs, [Path(p) for p in snakemake.input.tech_inputs]))
    cf_files    = {t: p for t, p in tech_files.items() if t != "grid"}
    grid_path   = tech_files.get("grid")

    price_series = _read_series(grid_path) if grid_path is not None else None

    if cf_files:
        cf_timeseries = pd.DataFrame({t: _read_series(p) for t, p in cf_files.items()})
        cf_timeseries.index = pd.to_datetime(cf_timeseries.index)
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

    n = build_network(assumptions, cf_timeseries, price_series)
    n.optimize(solver_name="highs")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(out_path)
    print(f"Network saved to {out_path}")


if __name__ == "__main__":
    main()
