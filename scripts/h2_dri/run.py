"""
CLI runner for the DRI-hydrogen PyPSA optimisation.

Usage (from repo root):
    python scripts/h2_dri/run.py --project DE_2023_baseline --scenario dedicated_res
    python scripts/h2_dri/run.py --project DE_2023_baseline  # runs all scenarios
"""

import argparse
import sys
from pathlib import Path

# pandas 3.0 defaults to ArrowStringArray for strings; xarray (used by PyPSA
# internally) doesn't support it. Set python-native strings before any import
# that touches pandas string data.
import pandas as pd
pd.options.mode.string_storage = "python"

import yaml

from network import build_network
from costs import compute_lcoh, extract_summary


REPO_ROOT = Path(__file__).parent.parent.parent
CONFIG_DIR = REPO_ROOT / "config"
RESULTS_DIR = REPO_ROOT / "results"


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_cf_dataframe(country: str, year: int, variant: str, cf_dir: Path) -> pd.DataFrame:
    """
    Load the atlite CF CSV and return a DataFrame with a DatetimeIndex.
    Columns are tech names (wind_onshore, solar, ...) without the '_cf' suffix.
    """
    cc = country.lower()
    filename = f"{cc}_cf_{year}.csv" if variant == "avg" else f"{cc}_cf_{year}_{variant}.csv"
    path = cf_dir / filename
    if not path.exists():
        raise FileNotFoundError(f"CF file not found: {path}")

    df = pd.read_csv(path, parse_dates=["time"], index_col="time")
    df.index.name = None
    df.columns = [c.replace("_cf", "") for c in df.columns]
    return df


def load_price_series(area: str, year: int, processed_path: Path) -> pd.Series:
    """
    Load hourly electricity prices for one area and year from entsoe_processed.feather.
    Returns a Series with DatetimeIndex (timezone-naive, UTC).
    """
    df = pd.read_feather(processed_path)
    prices = df[(area, "price")]
    prices = prices[prices.index.year == year]
    prices.index = prices.index.tz_localize(None) if prices.index.tz is not None else prices.index
    return prices


def run_scenario(
    project_cfg: dict,
    scenario_cfg: dict,
    assumptions: dict,
    data_cfg: dict,
) -> dict:
    """Build, solve, and return summary for one project+scenario."""
    cf_dir = REPO_ROOT / data_cfg["cf_dir"]
    cf_df = load_cf_dataframe(
        country=project_cfg["country"],
        year=project_cfg["year"],
        variant=project_cfg["cf_variant"],
        cf_dir=cf_dir,
    )

    # Keep only the techs requested for this scenario
    techs = scenario_cfg["techs"]
    missing = [t for t in techs if t not in cf_df.columns]
    if missing:
        raise KeyError(f"Techs {missing} not found in CF file columns: {list(cf_df.columns)}")
    cf_df = cf_df[techs]

    # Load prices if grid-connected
    price_series = None
    if scenario_cfg.get("grid_connected", False):
        processed_path = REPO_ROOT / data_cfg["processed_path"]
        raw_prices = load_price_series(
            area=scenario_cfg["grid_price_area"],
            year=project_cfg["year"],
            processed_path=processed_path,
        )
        # Align to CF index
        price_series = raw_prices.reindex(cf_df.index)
        if price_series.isna().any():
            raise ValueError(
                f"Price series has {price_series.isna().sum()} missing values after aligning to CF index. "
                "Check that entsoe_processed.feather covers the same year and hourly resolution as the CF file."
            )

    n = build_network(project_cfg, assumptions, cf_df, price_series)
    n.optimize(solver_name="highs")

    project_name = project_cfg["name"]
    scenario_name = scenario_cfg["name"]

    out_dir = RESULTS_DIR / project_name
    out_dir.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(out_dir / f"{scenario_name}.nc")

    summary = extract_summary(n, project_name, scenario_name)
    pd.DataFrame([summary]).to_csv(out_dir / f"{scenario_name}_summary.csv", index=False)

    return summary


def main():
    parser = argparse.ArgumentParser(description="Run DRI-hydrogen PyPSA optimisation")
    parser.add_argument("--project", help="Project name from projects.yaml (default: first project)")
    parser.add_argument("--scenario", help="Scenario name (default: all scenarios in project)")
    args = parser.parse_args()

    assumptions = load_yaml(CONFIG_DIR / "assumptions.yaml")
    projects_cfg = load_yaml(CONFIG_DIR / "projects.yaml")
    data_cfg = projects_cfg["data"]

    projects = projects_cfg["projects"]
    if args.project:
        projects = [p for p in projects if p["name"] == args.project]
        if not projects:
            print(f"Error: project '{args.project}' not found in projects.yaml", file=sys.stderr)
            sys.exit(1)

    for project_cfg in projects:
        scenarios = project_cfg["scenarios"]
        if args.scenario:
            scenarios = [s for s in scenarios if s["name"] == args.scenario]
            if not scenarios:
                print(f"Error: scenario '{args.scenario}' not in project '{project_cfg['name']}'", file=sys.stderr)
                sys.exit(1)

        for scenario_cfg in scenarios:
            label = f"{project_cfg['name']} / {scenario_cfg['name']}"
            print(f"\n{'='*60}\nRunning: {label}\n{'='*60}")

            summary = run_scenario(project_cfg, scenario_cfg, assumptions, data_cfg)

            print(f"  LCOH:             {summary['lcoh_eur_per_kg']:.3f} €/kg H₂")
            print(f"  Total cost/yr:    {summary['total_annual_cost_eur']:,.0f} €")
            for k, v in summary.items():
                if k.endswith("_mw_opt"):
                    print(f"  {k:<30} {v:,.1f} MW")
            if "battery_mwh_opt" in summary:
                print(f"  {'battery_mwh_opt':<30} {summary['battery_mwh_opt']:,.1f} MWh")
            if "h2_buffer_mwh_lhv_opt" in summary:
                print(f"  {'h2_buffer_mwh_lhv_opt':<30} {summary['h2_buffer_mwh_lhv_opt']:,.1f} MWh LHV")

            print(f"\n  Results saved to results/{project_cfg['name']}/{scenario_cfg['name']}.*")


if __name__ == "__main__":
    main()
