"""
PyPSA DRI-hydrogen optimisation entry point.

Invoked by Snakemake's `script:` directive (one rule fires per
(project, scenario) wildcard pair). The hardcoded fallback at the top of the
file lets you also run it standalone for ad-hoc dev work — same convention as
the res_cf scripts.

Standalone invocation:
    python scripts/h2_dri/run.py
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
from costs import compute_lcoh, extract_summary

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._paths import REPO_ROOT


# ── Defaults for standalone runs (DE_2023_baseline / dedicated_res) ──────────
PROJECT_NAME      = "DE_2023_baseline"
SCENARIO_NAME     = "dedicated_res"
TECH_INPUT_FILES  = [
    REPO_ROOT / "resources/res_cf/de_wind_onshore_20230101_20231231.parquet",
    REPO_ROOT / "resources/res_cf/de_solar_20230101_20231231.parquet",
]
ASSUMPTIONS_PATH  = REPO_ROOT / "config/assumptions.yaml"
PROJECTS_PATH     = REPO_ROOT / "config/projects.yaml"
OUT_NETWORK       = REPO_ROOT / "results" / PROJECT_NAME / f"{SCENARIO_NAME}.nc"
OUT_SUMMARY       = REPO_ROOT / "results" / PROJECT_NAME / f"{SCENARIO_NAME}_summary.csv"

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    PROJECT_NAME     = snakemake.wildcards.project
    SCENARIO_NAME    = snakemake.wildcards.scenario
    TECH_INPUT_FILES = list(snakemake.input.tech_inputs)
    ASSUMPTIONS_PATH = Path(snakemake.input.assumptions)
    PROJECTS_PATH    = Path(snakemake.input.projects)
    OUT_NETWORK      = Path(snakemake.output.network)
    OUT_SUMMARY      = Path(snakemake.output.summary)


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
    projects_cfg: dict,
    out_network: Path,
    out_summary: Path,
) -> dict:
    project_cfg  = projects_cfg["projects"][project_name]
    scenario_cfg = project_cfg["scenarios"][scenario_name]
    h2_lhv_kwh_per_kg = assumptions["h2"]["lhv_kwh_per_kg"]

    # Map each tech to its resolved input file (collect preserves list order).
    techs = list(scenario_cfg["techs"])
    tech_file_map = dict(zip(techs, [Path(f) for f in tech_input_files]))

    cf_files   = {t: f for t, f in tech_file_map.items() if t != "grid"}
    prices_path = tech_file_map.get("grid")

    cf_timeseries = load_cf_timeseries(cf_files)

    price_series = None
    if prices_path is not None:
        raw_prices = load_price_series(prices_path)
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

    summary = extract_summary(n, project_name, scenario_name, h2_lhv_kwh_per_kg)
    pd.DataFrame([summary]).to_csv(out_summary, index=False)

    return summary


def main() -> None:
    assumptions  = load_yaml(ASSUMPTIONS_PATH)
    projects_cfg = load_yaml(PROJECTS_PATH)

    summary = run(
        project_name=PROJECT_NAME,
        scenario_name=SCENARIO_NAME,
        tech_input_files=TECH_INPUT_FILES,
        assumptions=assumptions,
        projects_cfg=projects_cfg,
        out_network=OUT_NETWORK,
        out_summary=OUT_SUMMARY,
    )

    label = f"{PROJECT_NAME} / {SCENARIO_NAME}"
    print(f"\n{'=' * 60}\nRan: {label}\n{'=' * 60}")
    print(f"  LCOH:             {summary['lcoh_eur_per_kg']:.3f} €/kg H₂")
    print(f"  Total cost/yr:    {summary['total_annual_cost_eur']:,.0f} €")
    for k, v in summary.items():
        if k.endswith("_mw_opt"):
            print(f"  {k:<30} {v:,.1f} MW")
    if "battery_mwh_opt" in summary:
        print(f"  {'battery_mwh_opt':<30} {summary['battery_mwh_opt']:,.1f} MWh")
    if "h2_buffer_mwh_lhv_opt" in summary:
        print(f"  {'h2_buffer_mwh_lhv_opt':<30} {summary['h2_buffer_mwh_lhv_opt']:,.1f} MWh LHV")
    print(f"\n  Results saved to {OUT_NETWORK} / {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
