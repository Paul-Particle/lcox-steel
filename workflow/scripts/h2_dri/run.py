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


PROJECT_NAME      = snakemake.wildcards.project
SCENARIO_NAME     = snakemake.wildcards.scenario
TECH_INPUT_FILES  = list(snakemake.input.tech_inputs)
ASSUMPTIONS_PATH  = Path(snakemake.input.assumptions)
PROJECTS_PATH     = Path(snakemake.input.projects)
OUT_NETWORK       = Path(snakemake.output.network)


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_cf_timeseries(tech_files: dict[str, Path]) -> pd.DataFrame:
    """Load per-tech CF parquets and return a single DataFrame (tech columns)."""
    frames = {}
    for tech, path in tech_files.items():
        if not path.exists():
            raise FileNotFoundError(f"CF file not found: {path}")
        df = pd.read_parquet(path)
        if "time" in df.columns:
            df = df.set_index("time")
        df.index.name = None
        frames[tech] = df.iloc[:, 0]
    result = pd.DataFrame(frames)
    result.index = pd.to_datetime(result.index)
    return result


def load_price_series(path: Path) -> pd.Series:
    """Load a processed grid price parquet; return timezone-naive UTC series."""
    df = pd.read_parquet(path)
    if "time" in df.columns:
        df = df.set_index("time")
    df.index.name = None
    prices = df.iloc[:, 0]
    if prices.index.tz is not None:
        prices.index = prices.index.tz_convert("UTC").tz_localize(None)
    return prices


def run(
    project_name: str,
    scenario_name: str,
    tech_input_files: list[str | Path],
    assumptions: dict,
    projects_df: pd.DataFrame,
    out_network: Path,
) -> None:
    # Techs for this (project, scenario) in projects.csv row order — the same
    # order collect() produced tech_input_files in, so we can zip them.
    rows = projects_df.query(
        "project == @project_name and scenario == @scenario_name"
    )
    techs = rows["tech"].tolist()
    tech_file_map = dict(zip(techs, [Path(f) for f in tech_input_files]))

    cf_files   = {t: f for t, f in tech_file_map.items() if t != "grid"}
    prices_path = tech_file_map.get("grid")

    raw_prices = load_price_series(prices_path) if prices_path is not None else None

    if cf_files:
        cf_timeseries = load_cf_timeseries(cf_files)
    elif raw_prices is not None:
        # Grid-only scenario (no CF techs): the snapshots come from the price series.
        cf_timeseries = pd.DataFrame(index=raw_prices.index)
    else:
        raise ValueError(
            f"{project_name}/{scenario_name}: scenario has no CF techs and no grid input"
        )

    price_series = None
    if raw_prices is not None:
        price_series = raw_prices.reindex(cf_timeseries.index)
        if price_series.isna().any():
            raise ValueError(
                f"Price series has {price_series.isna().sum()} missing values after "
                "aligning to CF index. Check that the entsoe processed file covers "
                "the same period as the CF data."
            )

    n = build_network(assumptions, cf_timeseries, price_series)
    n.optimize(solver_name="highs")

    out_network.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(out_network)


def main() -> None:
    assumptions = load_yaml(ASSUMPTIONS_PATH)
    projects_df = pd.read_csv(PROJECTS_PATH, dtype={"start_date": str, "end_date": str})

    run(
        project_name=PROJECT_NAME,
        scenario_name=SCENARIO_NAME,
        tech_input_files=TECH_INPUT_FILES,
        assumptions=assumptions,
        projects_df=projects_df,
        out_network=OUT_NETWORK,
    )

    print(f"Network saved to {OUT_NETWORK}")


if __name__ == "__main__":
    main()
