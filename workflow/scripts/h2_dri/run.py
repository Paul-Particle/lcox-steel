"""
PyPSA DRI-hydrogen optimisation entry point.

Invoked by Snakemake's `script:` directive (one rule fires per
(project, scenario) wildcard pair). The hardcoded fallback at the top of the
file lets you also run it standalone for ad-hoc dev work — same convention as
the res_cf scripts.

Standalone invocation:
    python scripts/h2_dri/run.py
"""

from __future__ import annotations

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


# ── Defaults for standalone runs (mirror DE_2023_baseline / dedicated_res) ──
PROJECT_NAME    = "DE_2023_baseline"
SCENARIO_NAME   = "dedicated_res"
CF_PATH         = REPO_ROOT / "resources/res_cf/de_cf_2023.csv"
PRICES_PATH     = REPO_ROOT / "resources/entsoe_processed.feather"
ASSUMPTIONS_PATH = REPO_ROOT / "config/assumptions.yaml"
PROJECTS_PATH   = REPO_ROOT / "config/projects.yaml"
OUT_NETWORK     = REPO_ROOT / "results" / PROJECT_NAME / f"{SCENARIO_NAME}.nc"
OUT_SUMMARY     = REPO_ROOT / "results" / PROJECT_NAME / f"{SCENARIO_NAME}_summary.csv"

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    PROJECT_NAME    = snakemake.wildcards.project
    SCENARIO_NAME   = snakemake.wildcards.scenario
    CF_PATH         = Path(snakemake.input.cf)
    try:
        PRICES_PATH = Path(snakemake.input.prices)
    except AttributeError:
        PRICES_PATH = None
    ASSUMPTIONS_PATH = Path(snakemake.input.assumptions)
    PROJECTS_PATH   = Path(snakemake.input.projects)
    OUT_NETWORK     = Path(snakemake.output.network)
    OUT_SUMMARY     = Path(snakemake.output.summary)


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_cf_timeseries(path: Path) -> pd.DataFrame:
    """
    Load the atlite CF CSV and return a DataFrame with a DatetimeIndex.
    Columns are tech names (wind_onshore, solar, ...) without the '_cf' suffix.
    """
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


def _find(items: list[dict], name: str, kind: str) -> dict:
    for it in items:
        if it["name"] == name:
            return it
    raise KeyError(f"{kind} '{name}' not found")


def run(
    project_name: str,
    scenario_name: str,
    cf_path: Path,
    prices_path: Path | None,
    assumptions: dict,
    projects_cfg: dict,
    out_network: Path,
    out_summary: Path,
) -> dict:
    project_cfg  = _find(projects_cfg["projects"], project_name, "Project")
    scenario_cfg = _find(project_cfg["scenarios"], scenario_name, "Scenario")
    h2_lhv_kwh_per_kg = assumptions["h2"]["lhv_kwh_per_kg"]

    cf_timeseries = load_cf_timeseries(cf_path)

    techs = scenario_cfg["techs"]
    missing = [t for t in techs if t not in cf_timeseries.columns]
    if missing:
        raise KeyError(f"Techs {missing} not found in CF file columns: {list(cf_timeseries.columns)}")
    cf_timeseries = cf_timeseries[techs] if techs else cf_timeseries.iloc[:, :0]

    price_series = None
    if scenario_cfg.get("grid_connected", False):
        if prices_path is None:
            raise ValueError(
                f"Scenario '{scenario_name}' is grid_connected but no prices input was provided."
            )
        if "grid_price_area" not in scenario_cfg:
            raise KeyError(
                f"Scenario '{scenario_name}' is grid_connected but lacks 'grid_price_area' in projects.yaml."
            )
        raw_prices = load_price_series(
            area=scenario_cfg["grid_price_area"],
            year=project_cfg["year"],
            processed_path=prices_path,
        )
        price_series = raw_prices.reindex(cf_timeseries.index)
        if price_series.isna().any():
            raise ValueError(
                f"Price series has {price_series.isna().sum()} missing values after aligning to CF index. "
                "Check that entsoe_processed.feather covers the same year and hourly resolution as the CF file."
            )

    n = build_network(project_cfg, assumptions, cf_timeseries, price_series)
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
        cf_path=CF_PATH,
        prices_path=PRICES_PATH,
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
