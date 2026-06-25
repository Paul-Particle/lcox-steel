# =============================================================================
# NOTE – Interpretation of Stage 3B diagnostics plots
#
# These plots DO NOT replay the exact annual Stage 3B optimisation trajectory.
#
# What this script does:
# - Takes the Stage 3B optimal sizing (RES, electrolyser, storage)
# - Re-simulates battery dispatch chronologically on SHORT TIME WINDOWS (e.g. weeks)
# - Uses a RESET initial state of charge (SOC) for each window (empty or full)
#
# Implications:
# - SOC is NOT continuous across the year
# - The plots are NOT exact slices of the annual optimum dispatch
# - Weekly behaviour may differ slightly from the true annual trajectory
#
# What the plots ARE useful for:
# - Understanding intra-week dynamics (surplus, deficit, charging, discharging)
# - Visualising the role of storage (battery → electrolyser contribution)
# - Explaining why storage increases served energy (and thus H2 output / denominator)
#
# What the plots should NOT be used for:
# - Exact validation of annual unserved energy or curtailment
# - Inferring precise SOC levels at specific dates in the year
#
# In short:
# → These are ILLUSTRATIVE DISPATCH DIAGNOSTICS, not exact annual replay.
# =============================================================================


import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import matplotlib.pyplot as plt

CONFIG_PATH = Path("config_hannah.yaml")
DATA_DIR = Path("data/res_cf")
RESULTS_DIR = Path("results")
PLOTS_DIR = RESULTS_DIR / "plots" / "diagnostics"


# ----------------------------
# Helpers
# ----------------------------
def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def scenario_to_name(scenario: dict) -> str:
    return f"{scenario['country']}_{scenario['tech']}_{scenario['variant']}_{scenario['year']}"


def load_atlite_cf_with_time(country: str, tech: str, variant: str, year: int) -> pd.Series:
    cc = country.lower()

    if variant == "avg":
        filename = f"{cc}_cf_{year}.csv"
    elif variant == "bestsite_p95":
        filename = f"{cc}_cf_{year}_bestsite_p95.csv"
    else:
        raise ValueError(f"Unknown variant: {variant}")

    path = DATA_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing CF file: {path}")

    df = pd.read_csv(path, parse_dates=["time"])
    col = f"{tech}_cf"
    if col not in df.columns:
        raise KeyError(f"Column {col} not found in {filename}")

    return pd.Series(df[col].astype(float).to_numpy(), index=pd.to_datetime(df["time"]))


def read_stage3b_best_row(scenario_name: str) -> pd.Series:
    path = RESULTS_DIR / f"{scenario_name}_stage3b_sweep.csv"

    if not path.exists():
        raise FileNotFoundError(f"Missing Stage 3B sweep file: {path}")

    df = pd.read_csv(path)
    df = df[df["storage_feasible"] == True].copy()
    df = df[df["lcoh_with_storage_eur_per_kg"].notna()].copy()

    if df.empty:
        raise ValueError(f"No feasible Stage 3B rows for {scenario_name}")

    best_idx = df["lcoh_with_storage_eur_per_kg"].idxmin()
    best_row = df.loc[best_idx]

    print(
        f"{scenario_name} | Stage3B optimum | "
        f"ratio={best_row['res_to_el_ratio']:.2f} | "
        f"LCOH={best_row['lcoh_with_storage_eur_per_kg']:.2f} | "
        f"storage={best_row['storage_capacity_required_mwh']:.1f} MWh"
    )

    return best_row


# ----------------------------
# Stage 3B time-series simulation
# ----------------------------
def simulate_stage3b_timeseries(
    cf_window: pd.Series,
    electrolyser_mw: float,
    res_capacity_mw: float,
    storage_capacity_mwh: float,
    eta_rt: float,
    initial_soc_mode: str = "empty",
):
    eta_c = np.sqrt(eta_rt)
    eta_d = np.sqrt(eta_rt)

    p_res = cf_window.to_numpy() * float(res_capacity_mw)
    p_el = float(electrolyser_mw)

    direct_res_to_el = np.minimum(p_res, p_el)
    surplus_mw = np.maximum(p_res - p_el, 0.0)
    deficit_mw = np.maximum(p_el - p_res, 0.0)

    n = len(p_res)

    if initial_soc_mode == "empty":
        soc = 0.0
    elif initial_soc_mode == "full":
        soc = float(storage_capacity_mwh)
    else:
        raise ValueError("initial_soc_mode must be 'empty' or 'full'")

    charge_mw = np.zeros(n)
    discharge_mw = np.zeros(n)
    unserved_mw = np.zeros(n)
    curtailed_after_storage_mw = np.zeros(n)
    soc_mwh = np.zeros(n)

    for t in range(n):
        surplus = surplus_mw[t]
        deficit = deficit_mw[t]

        if surplus > 0:
            max_input_due_to_capacity = (storage_capacity_mwh - soc) / eta_c
            max_input_due_to_capacity = max(0.0, max_input_due_to_capacity)

            charge = min(surplus, max_input_due_to_capacity)
            soc += charge * eta_c

            charge_mw[t] = charge
            curtailed_after_storage_mw[t] = surplus - charge

        elif deficit > 0:
            max_output_due_to_soc = soc * eta_d
            discharge = min(deficit, max_output_due_to_soc)

            soc -= discharge / eta_d

            discharge_mw[t] = discharge
            unserved_mw[t] = deficit - discharge

        soc_mwh[t] = soc

    return pd.DataFrame(
        {
            "time": cf_window.index,
            "cf": cf_window.to_numpy(),
            "p_res_mw": p_res,
            "p_el_mw": p_el,
            "direct_res_to_el_mw": direct_res_to_el,
            "surplus_mw": surplus_mw,
            "deficit_mw": deficit_mw,
            "charge_mw": charge_mw,
            "discharge_mw": discharge_mw,
            "unserved_mw": unserved_mw,
            "curtailed_after_storage_mw": curtailed_after_storage_mw,
            "soc_mwh": soc_mwh,
        }
    ).set_index("time")


# ----------------------------
# Plot
# ----------------------------
def plot_window(
    *,
    scenario_name: str,
    window_label: str,
    sim_df: pd.DataFrame,
    ratio: float,
    res_capacity_mw: float,
    electrolyser_mw: float,
    storage_capacity_mwh: float,
    out_path: Path,
):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )

    # ---- top panel: power balance ----
    ax1.plot(sim_df.index, sim_df["p_res_mw"], label="RES generation", linewidth=1.8)
    ax1.axhline(electrolyser_mw, linestyle="--", label="Electrolyser load", linewidth=1.4)
    ax1.plot(sim_df.index, sim_df["direct_res_to_el_mw"], label="Direct RES → EL", linewidth=1.5)
    ax1.plot(sim_df.index, sim_df["discharge_mw"], label="Battery → EL", linewidth=1.5)
    ax1.plot(sim_df.index, sim_df["unserved_mw"], label="Unserved", linewidth=1.5)

    ax1.set_ylabel("Power [MW]")
    ax1.set_title(
        f"{scenario_name} | Stage 3B optimum | {window_label}\n"
        f"ratio={ratio:.2f} | RES={res_capacity_mw:,.0f} MW | "
        f"EL={electrolyser_mw:,.0f} MW | Storage={storage_capacity_mwh:,.0f} MWh"
    )
    ax1.grid(True)
    ax1.legend(loc="upper right", ncol=2)

    # ---- bottom panel: storage state ----
    ax2.plot(sim_df.index, sim_df["soc_mwh"], label="SOC", linewidth=1.8)
    ax2.set_ylabel("SOC [MWh]")
    ax2.set_xlabel("Time")
    ax2.grid(True)
    ax2.legend(loc="upper right")

    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=180)
    plt.close()


# ----------------------------
# Main
# ----------------------------
def main():
    cfg = load_config(CONFIG_PATH)
    defaults = cfg.get("defaults", {})
    scenarios = cfg.get("scenarios", [])

    if not isinstance(scenarios, list) or len(scenarios) == 0:
        raise ValueError("config_hannah.yaml must contain a non-empty list under `scenarios:`")

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    jan_week = ("Jan week", "2023-01-15", "2023-01-21")
    jul_week = ("Jul week", "2023-07-01", "2023-07-07")

    for i, scenario_raw in enumerate(scenarios, start=1):
        scenario = deep_merge(defaults, scenario_raw)
        scenario_name = scenario_to_name(scenario)

        print("\n" + "=" * 80)
        print(f"Plotting scenario {i}/{len(scenarios)}: {scenario_name}")
        print("=" * 80)

        cf = load_atlite_cf_with_time(
            scenario["country"],
            scenario["tech"],
            scenario["variant"],
            int(scenario["year"]),
        )

        best = read_stage3b_best_row(scenario_name)

        ratio = float(best["res_to_el_ratio"])
        electrolyser_mw = float(best["electrolyser_mw"])
        res_capacity_mw = float(best["res_capacity_mw"])
        storage_capacity_mwh = float(best["storage_capacity_required_mwh"])

        storage_cfg = scenario["storage_stage3b"]
        eta_rt = float(storage_cfg["roundtrip_efficiency"])
        initial_soc_mode = str(storage_cfg.get("initial_soc_mode", "empty"))

        # ----------------------------
# FULL-YEAR VALIDATION (Stage 3B)
# ----------------------------
        full_sim_df = simulate_stage3b_timeseries(
            cf_window=cf,  # full year
            electrolyser_mw=electrolyser_mw,
            res_capacity_mw=res_capacity_mw,
            storage_capacity_mwh=storage_capacity_mwh,
            eta_rt=eta_rt,
            initial_soc_mode=initial_soc_mode,
        )

        # annual metrics
        annual_unserved_mwh = full_sim_df["unserved_mw"].sum()  # hourly → MWh

        soc_min = full_sim_df["soc_mwh"].min()
        soc_max = full_sim_df["soc_mwh"].max()

        annual_baseload_mwh = electrolyser_mw * 8760
        allowed_uncovered_mwh = (
            scenario["storage_stage3b"]["allowed_energy_deficit_fraction"]
            * annual_baseload_mwh
        )

        print(
            f"[VALIDATION] SOC range: {soc_min:.1f} → {soc_max:.1f} MWh "
            f"(capacity = {storage_capacity_mwh:.1f})"
        )
        print(
            f"[VALIDATION] Annual unserved: {annual_unserved_mwh:.1f} MWh "
            f"(allowed = {allowed_uncovered_mwh:.1f})"
        )

        windows = [
            (jan_week[0], cf.loc[jan_week[1]:jan_week[2]]),
            (jul_week[0], cf.loc[jul_week[1]:jul_week[2]]),
        ]

        for window_label, cf_window in windows:
            if cf_window.empty:
                print(f"Warning: empty CF window for {scenario_name} | {window_label}. Skipping.")
                continue

            sim_df = full_sim_df.loc[cf_window.index]
            safe_window = "jan_week" if "Jan" in window_label else "jul_week"
            out_name = f"{scenario_name}__stage3b_opt__{safe_window}.png"
            out_path = PLOTS_DIR / out_name

            plot_window(
                scenario_name=scenario_name,
                window_label=window_label,
                sim_df=sim_df,
                ratio=ratio,
                res_capacity_mw=res_capacity_mw,
                electrolyser_mw=electrolyser_mw,
                storage_capacity_mwh=storage_capacity_mwh,
                out_path=out_path,
            )

            print("Saved:", out_path)

    print("\nDone. Plots in:", PLOTS_DIR.resolve())


if __name__ == "__main__":
    main()