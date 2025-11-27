#Loads config.yaml and grabs your first scenario (DE_wind_onshore_baseline).
#Loads data/processed_data.feather.
#Selects the ('DE_LU', 'wind_onshore') column as res_series.
#Normalises that to a per-unit profile (0–1) for the optimisation sweep.
#Infers the time step from the index (e.g. 0.25 h for 15-min data).
#Runs optimise_res_to_el_ratio(...).
#Prints:
#scenario name
#area / RES type
#the best RES:electrolyser ratio
#all metrics at that best point
#Saves a CSV of the full sweep to results/DE_wind_onshore_baseline_res_to_h2_sweep.csv.

import pandas as pd
import numpy as np
from pathlib import Path

import yaml

from res_to_h2_logic import optimise_res_to_el_ratio


CONFIG_PATH = Path("config_hannah.yaml")
PROCESSED_DATA_PATH = Path("data/processed_data.feather")


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    # 1) Load config
    config = load_config(CONFIG_PATH)
    scenario = config["scenarios"][0]  # MVP: single scenario

    area = scenario["area"]          # e.g. "DE_LU"
    res_col = scenario["res_column"] # e.g. "wind_onshore"

    # 2) Load processed data
    df = pd.read_feather(PROCESSED_DATA_PATH)

    # Ensure MultiIndex on columns (should already be the case, but safe)
    if not isinstance(df.columns, pd.MultiIndex):
        df.columns = pd.MultiIndex.from_tuples(df.columns)

    # 3) Extract RES profile for this scenario (DE_LU / wind_onshore)
    try:
        res_series = df[(area, res_col)]
    except KeyError as e:
        raise KeyError(
            f"Could not find column ({area!r}, {res_col!r}) in processed_data.feather"
        ) from e

    # Normalise to per-unit profile (0..1) for the optimisation.
    # If the data is already p.u., you can skip this step.
    max_val = res_series.max()
    if max_val <= 0:
        raise ValueError(f"RES column ({area}, {res_col}) has non-positive max value.")
    res_profile_pu = (res_series / max_val).to_numpy()

    # 4) Determine timestep length [h] from index
    idx = res_series.index
    if len(idx) < 2:
        raise ValueError("Time series must have at least 2 timesteps to infer resolution.")

    timestep_hours = (idx[1] - idx[0]).total_seconds() / 3600.0

    # 5) Run optimisation
    results_df, best_ratio = optimise_res_to_el_ratio(
        res_profile=res_profile_pu,
        config=config,
        timestep_hours=timestep_hours,
    )

    # 6) Print summary
    print("\n=== Dedicated RES → H₂ optimisation ===")
    print(f"Scenario: {scenario['name']}")
    print(f"Area / RES: {area} / {res_col}")
    print(f"Best RES:electrolyser ratio: {best_ratio:.3f}")
    print("\nBest point:")
    print(results_df.loc[best_ratio])

    # 7) Save full sweep results
    out_dir = Path("results")
    out_dir.mkdir(exist_ok=True)

    out_name = f"{scenario['name']}_res_to_h2_sweep.csv"
    out_path = out_dir / out_name
    results_df.to_csv(out_path)

    print(f"\nFull sweep saved to: {out_path}")


if __name__ == "__main__":
    main()
