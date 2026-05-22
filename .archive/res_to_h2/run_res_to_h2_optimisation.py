import numpy as np
import pandas as pd
from pathlib import Path
import yaml
from datetime import datetime

from res_to_h2_logic import optimise_res_to_el_ratio, dri_to_el_mw, firm_ratio_from_sweep

CONFIG_PATH = Path("config_hannah.yaml")
DATA_DIR = Path("data/res_cf")


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def deep_merge(base: dict, override: dict) -> dict:
    """
    Recursively merges override into base (without mutating inputs).
    Values in override take precedence.
    """
    out = dict(base)  # shallow copy
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_atlite_profile(country: str, tech: str, variant: str, year: int):
    """
    Loads an Atlite-based CF time series.

    country: e.g. "DE"
    tech: e.g. "wind_onshore"
    variant: "avg" or "bestsite_p95"
    year: e.g. 2023
    """
    cc = country.lower()

    if variant == "avg":
        filename = f"{cc}_cf_{year}.csv"
    elif variant == "bestsite_p95":
        filename = f"{cc}_cf_{year}_bestsite_p95.csv"
    else:
        raise ValueError(f"Unknown variant: {variant}")

    path = DATA_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Could not find Atlite CF file: {path}")

    df = pd.read_csv(path, parse_dates=["time"])

    col = f"{tech}_cf"
    if col not in df.columns:
        raise KeyError(f"Column {col} not found in {filename}")

    res_series = df[col].astype(float)

    if not np.isfinite(res_series.to_numpy()).all():
        raise ValueError("CF profile contains NaN or inf values.")
    if res_series.max() <= 0:
        raise ValueError("CF profile has non-positive maximum.")

    timestep_hours = 1.0  # Atlite files are hourly
    return res_series.to_numpy(), timestep_hours


def is_bad_fstring_name(name_from_cfg) -> bool:
    """
    YAML can't evaluate f-strings. Catch common "f\"{...}\"" / braces patterns.
    """
    if not name_from_cfg:
        return True
    s = str(name_from_cfg)
    return ("{" in s) or s.startswith("f\"") or s.startswith("f'")


def main():
    print("Loading config from:", CONFIG_PATH.resolve())
    config = load_config(CONFIG_PATH)

    defaults = config.get("defaults", {})
    scenarios = config.get("scenarios", [])
    if not isinstance(scenarios, list) or len(scenarios) == 0:
        raise ValueError("config_hannah.yaml must contain a non-empty list under `scenarios:`")

    out_dir = Path("results")
    out_dir.mkdir(exist_ok=True)

    for i, scenario_raw in enumerate(scenarios, start=1):
        # 0) Merge defaults + scenario overrides
        scenario = deep_merge(defaults, scenario_raw)

        country = scenario["country"]
        tech = scenario["tech"]
        variant = scenario["variant"]
        year = int(scenario["year"])

        # Scenario name: generate in Python if missing/invalid
        name_from_cfg = scenario.get("name")
        if is_bad_fstring_name(name_from_cfg):
            scenario_name = f"{country}_{tech}_{variant}_{year}"
        else:
            scenario_name = str(name_from_cfg)
        scenario["name"] = scenario_name

        print("\n" + "=" * 80)
        print(f"Running scenario {i}/{len(scenarios)}: {scenario_name}")
        print("=" * 80)

        # 1) Load Atlite CF profile
        res_profile, timestep_hours = load_atlite_profile(
            country=country,
            tech=tech,
            variant=variant,
            year=year,
        )

        print("\n--- CF diagnostics ---")
        print("CF min :", float(res_profile.min()))
        print("CF max :", float(res_profile.max()))
        print("CF mean:", float(res_profile.mean()))
        print("----------------------\n")

        # 2) DRI-scale sizing (from YAML)
        dri_cfg = scenario["dri_sizing"]
        dri_mt_per_year = float(dri_cfg["dri_mt_per_year"])
        h2_intensity = float(dri_cfg["h2_intensity_kg_per_t_dri"])
        availability_target = float(dri_cfg["availability_target"])

        el_eff_kwh_per_kg = float(scenario["electrolyser"]["efficiency_kwh_per_kg"])

        el_mw_dri = dri_to_el_mw(
            dri_mt_per_year=dri_mt_per_year,
            h2_intensity_kg_per_t_dri=h2_intensity,
            efficiency_kwh_per_kg=el_eff_kwh_per_kg,
            availability_target=availability_target,
        )

        print("Optimization grid used:", scenario["optimization"])

        # 3) Run optimisation using DRI-sized electrolyser as the sweep baseline
        results_df, best_ratio = optimise_res_to_el_ratio(
            res_profile=res_profile,
            scenario=scenario,
            finance=config["finance"],
            timestep_hours=timestep_hours,
            electrolyser_mw=el_mw_dri,
        )

        print("Sweep ratio range:", float(results_df.index.min()), "->", float(results_df.index.max()))
        print("Firm reliability at max ratio:", float(results_df.loc[results_df.index.max(), "firm_reliability"]))

        print("\n=== Dedicated RES → H₂ optimisation ===")
        print(f"Scenario: {scenario_name}")
        print(f"Country / Tech: {country} / {tech}")
        print(f"Variant: {variant}")
        print(f"Best RES:electrolyser ratio: {best_ratio:.3f}")

        # DRI-scale RES buildout using optimal ratio
        res_mw_dri = float(best_ratio) * el_mw_dri

        print("\nDRI-scale sizing (illustrative):")
        print(f"DRI output [Mt/y]                : {dri_mt_per_year:,.2f}")
        print(f"Assumed H2 intensity [kg/t]      : {h2_intensity:,.1f}")
        print(f"Availability target [%]          : {availability_target*100:,.1f}")
        print(f"Electrolyser size [MW]           : {el_mw_dri:,.1f}")
        print(f"RES buildout [MW] (ratio-applied): {res_mw_dri:,.1f}")

        # --- Step A: firming baseline (no storage) ---
        # Robust to missing or null `firming:` in YAML
        firm_target = float(scenario["firming"]["firm_power_fraction"])
        firm_info = firm_ratio_from_sweep(results_df, target_reliability=firm_target)

        firm_ratio_interp = firm_info["firm_ratio_interp"]
        firm_ratio_feasible = firm_info["firm_ratio_feasible"]
        firm_row_feasible = firm_info["firm_row_feasible"]
        #if firm_row_feasible is None:
         #   print("DEBUG firm_row_feasible: None (target not achievable in sweep)")
        #else:
         #   print("DEBUG firm_row_feasible annual_storage_cost_eur:", firm_row_feasible.get("annual_storage_cost_eur", "MISSING"))
          #  print("DEBUG firm_row_feasible storage_capacity_required_mwh:", firm_row_feasible.get("storage_capacity_required_mwh", "MISSING"))
           # print("DEBUG firm_row_feasible lcoh_with_storage_eur_per_kg:", firm_row_feasible.get("lcoh_with_storage_eur_per_kg", "MISSING"))

        # Cost-optimal row (intermittent)
        best_row = results_df.loc[best_ratio]

        summary = {
            "scenario": scenario_name,
            "country": country,
            "tech": tech,
            "variant": variant,
            "year": year,

            "electrolyser_mw": float(el_mw_dri),

            # --- Cost-optimal intermittent design ---
            "best_ratio_res_per_el": float(best_ratio),
            "best_res_capacity_mw": float(best_row["res_capacity_mw"]),
            "best_lcoh_eur_per_kg": float(best_row["lcoh_eur_per_kg"]),
            "best_el_utilisation_fraction": float(best_row["el_utilisation"]),
            "best_curtailment_rate_fraction": float(best_row["curtailment_rate"]),
            "best_firm_reliability_fraction": float(best_row["firm_reliability"]),
            "best_energy_served_fraction_no_storage": float(best_row["energy_served_fraction_no_storage"]) if "energy_served_fraction_no_storage" in best_row else np.nan,
            "best_energy_served_fraction_proxy": float(best_row["energy_served_fraction_proxy"]) if "energy_served_fraction_proxy" in best_row else np.nan,
            "best_energy_unserved_mwh_proxy": float(best_row["energy_unserved_mwh_proxy"]) if "energy_unserved_mwh_proxy" in best_row else np.nan,
            "best_lcoh_with_storage_eur_per_kg": float(best_row["lcoh_with_storage_eur_per_kg"]) if "lcoh_with_storage_eur_per_kg" in best_row else np.nan,
            "best_annual_storage_cost_eur_per_year": float(best_row["annual_storage_cost_eur"]) if "annual_storage_cost_eur" in best_row else np.nan,
            "best_storage_capacity_required_mwh": float(best_row["storage_capacity_required_mwh"]) if "storage_capacity_required_mwh" in best_row else np.nan,
            "best_annual_deficit_mwh": float(best_row["annual_deficit_mwh"]) if "annual_deficit_mwh" in best_row else np.nan,
            "best_storage_required_bool": bool(float(best_row.get("annual_storage_cost_eur", 0.0)) > 0.0),
            # --- Firming target ---
            "firm_target_fraction": float(firm_target),
            "firm_ratio_feasible_res_per_el": float(firm_ratio_feasible) if np.isfinite(firm_ratio_feasible) else np.nan,
            "firm_ratio_interp_res_per_el": float(firm_ratio_interp) if np.isfinite(firm_ratio_interp) else np.nan,
        }
    

        if firm_row_feasible is None or not np.isfinite(firm_ratio_feasible):
            summary.update({
                "firm_achievable": False,
                "firm_res_capacity_mw": np.nan,
                "firm_lcoh_eur_per_kg": np.nan,
                "firm_el_utilisation_fraction": np.nan,
                "firm_curtailment_rate_fraction": np.nan,
                "firm_firm_reliability_fraction": np.nan,
                "delta_res_capacity_mw_vs_best": np.nan,
                "delta_lcoh_eur_per_kg_vs_best": np.nan,
                "firm_lcoh_with_storage_eur_per_kg": np.nan,
                "firm_annual_storage_cost_eur_per_year": np.nan,
                "firm_storage_capacity_required_mwh": np.nan,
                "firm_annual_deficit_mwh": np.nan,
                "firm_storage_required_bool": np.nan,
                "firm_energy_served_fraction_no_storage": np.nan,
                "firm_energy_served_fraction_proxy": np.nan,
                "firm_energy_unserved_mwh_proxy": np.nan,

            })
        else:
            summary.update({
                "firm_achievable": True,
                "firm_res_capacity_mw": float(firm_row_feasible["res_capacity_mw"]),
                "firm_lcoh_eur_per_kg": float(firm_row_feasible["lcoh_eur_per_kg"]),
                "firm_el_utilisation_fraction": float(firm_row_feasible["el_utilisation"]),
                "firm_curtailment_rate_fraction": float(firm_row_feasible["curtailment_rate"]),
                "firm_firm_reliability_fraction": float(firm_row_feasible["firm_reliability"]),
                "firm_lcoh_with_storage_eur_per_kg": float(firm_row_feasible["lcoh_with_storage_eur_per_kg"]) if "lcoh_with_storage_eur_per_kg" in firm_row_feasible else np.nan,
                "firm_annual_storage_cost_eur_per_year": float(firm_row_feasible["annual_storage_cost_eur"]) if "annual_storage_cost_eur" in firm_row_feasible else np.nan,
                "firm_storage_capacity_required_mwh": float(firm_row_feasible["storage_capacity_required_mwh"]) if "storage_capacity_required_mwh" in firm_row_feasible else np.nan,
                "firm_annual_deficit_mwh": float(firm_row_feasible["annual_deficit_mwh"]) if "annual_deficit_mwh" in firm_row_feasible else np.nan,
                "firm_storage_required_bool": bool(float(firm_row_feasible.get("annual_storage_cost_eur", 0.0)) > 0.0),
                "delta_res_capacity_mw_vs_best": float((firm_ratio_interp - float(best_ratio)) * float(el_mw_dri)),
                "delta_lcoh_eur_per_kg_vs_best": float(firm_row_feasible["lcoh_eur_per_kg"] - best_row["lcoh_eur_per_kg"]),
                "firm_energy_served_fraction_no_storage": float(firm_row_feasible["energy_served_fraction_no_storage"]) if "energy_served_fraction_no_storage" in firm_row_feasible else np.nan,
                "firm_energy_served_fraction_proxy": float(firm_row_feasible["energy_served_fraction_proxy"]) if "energy_served_fraction_proxy" in firm_row_feasible else np.nan,
                "firm_energy_unserved_mwh_proxy": float(firm_row_feasible["energy_unserved_mwh_proxy"]) if "energy_unserved_mwh_proxy" in firm_row_feasible else np.nan,
            })

        pd.DataFrame([summary]).to_csv(out_dir / f"{scenario_name}_res_to_h2_summary.csv", index=False)       

        # ------------------------------------------------------------
        # Update comparison results table
        # ------------------------------------------------------------

        comparison_path = Path("results/res_to_h2_comparison_table.csv")

        firm_achievable = firm_row_feasible is not None

        new_row = {
            "scenario": scenario_name,
            "country": country,
            "tech": tech,
            "variant": variant,
            "year": year,

            "best_lcoh": float(best_row["lcoh_eur_per_kg"]),
            "best_utilisation": float(best_row["el_utilisation"]),
            "best_curtailment": float(best_row["curtailment_rate"]),

            "firm_achievable": bool(firm_achievable),
            "firm_lcoh": float(firm_row_feasible["lcoh_eur_per_kg"]) if firm_achievable else np.nan,
            "firm_utilisation": float(firm_row_feasible["el_utilisation"]) if firm_achievable else np.nan,
            "firm_curtailment": float(firm_row_feasible["curtailment_rate"]) if firm_achievable else np.nan,

            "run_timestamp_utc": datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        }

        new_df = pd.DataFrame([new_row])

        # If table exists → load and update
        if comparison_path.exists():

            comp = pd.read_csv(comparison_path)

            # remove existing row for same scenario
            comp = comp[comp["scenario"] != scenario_name]

            # append new row
            comp = pd.concat([comp, new_df], ignore_index=True)

        else:

            comp = new_df

        # save updated table
        comp.to_csv(comparison_path, index=False, na_rep="NaN")

        print("\n=== Firming baseline (no storage) ===")
        print(f"Firm reliability target (share of hours P_RES >= P_EL): {firm_target:.2f}")

        if firm_row_feasible is None or not np.isfinite(firm_ratio_feasible):
            print("Target not achievable within sweep range. Increase max_res_per_el_mw in YAML.")
            print("\n=== Comparison: Cost-optimal intermittent design vs Oversizing-only firm design ===")
            print("Cost-optimal intermittent design:")
            print(f"  ratio                 : {float(best_ratio):.3f}")
            print(f"  lcoh_eur_per_kg        : {best_row['lcoh_eur_per_kg']:.3f}")
            print(f"  firm_reliability       : {best_row['firm_reliability']:.3f}")
            print(f"  el_utilisation         : {best_row['el_utilisation']*100:.1f}%")
            print(f"  curtailment_rate       : {best_row['curtailment_rate']*100:.1f}%")
        else:
            # Oversizing deltas relative to best LCOH point (anchor)
            add_res_mw = (firm_ratio_interp - float(best_ratio)) * el_mw_dri

            print(f"First feasible sweep ratio (>= target): {firm_ratio_feasible:.3f}")
            print(f"Interpolated ratio at target:          {firm_ratio_interp:.3f}")
            print(f"Additional RES MW vs best-LCOH point:  {add_res_mw:,.1f} MW")

            print("\nImplied metrics at first feasible sweep point (conservative):")
            print(f"  firm_reliability       : {firm_row_feasible['firm_reliability']:.3f}")
            print(f"  lcoh_eur_per_kg        : {firm_row_feasible['lcoh_eur_per_kg']:.3f}")
            print(f"  el_utilisation         : {firm_row_feasible['el_utilisation']*100:.1f}%")
            print(f"  curtailment_rate       : {firm_row_feasible['curtailment_rate']*100:.1f}%")

            print("\n=== Comparison: Cost-optimal intermittent design vs Oversizing-only firm95 design ===")
            print("Cost-optimal intermittent design:")
            print(f"  ratio                 : {float(best_ratio):.3f}")
            print(f"  lcoh_eur_per_kg        : {best_row['lcoh_eur_per_kg']:.3f}")
            print(f"  firm_reliability       : {best_row['firm_reliability']:.3f}")
            print(f"  el_utilisation         : {best_row['el_utilisation']*100:.1f}%")
            print(f"  curtailment_rate       : {best_row['curtailment_rate']*100:.1f}%")

            print("Oversizing-only firm95 design (first feasible sweep point):")
            print(f"  ratio                 : {firm_ratio_feasible:.3f}")
            print(f"  lcoh_eur_per_kg        : {firm_row_feasible['lcoh_eur_per_kg']:.3f}")
            print(f"  firm_reliability       : {firm_row_feasible['firm_reliability']:.3f}")
            print(f"  el_utilisation         : {firm_row_feasible['el_utilisation']*100:.1f}%")
            print(f"  curtailment_rate       : {firm_row_feasible['curtailment_rate']*100:.1f}%")

        print("\n=== Detailed row: Cost-optimal intermittent design ===")

        def fmt(x):
            if isinstance(x, (int, float, np.floating, np.integer)):
                return f"{float(x):,.2f}"
            return x

        for k, v in best_row.items():
            print(f"{k:30s} : {fmt(v)}")

        if ("el_utilisation" in best_row.index) and ("curtailment_rate" in best_row.index):
            print("\nSystem metrics:")
            print(f"Electrolyser utilisation [%] : {best_row['el_utilisation']*100:.1f}")
            print(f"Curtailment rate [%]        : {best_row['curtailment_rate']*100:.1f}")
        else:
            print("\nSystem metrics: (not available — missing el_utilisation / curtailment_rate)")

        if (
            "lcoe_res_prod_eur_per_mwh" in best_row.index
            and "lcoe_res_used_eur_per_mwh" in best_row.index
        ):
            print("\nDerived RES LCOE metrics:")
            print(f"LCOE (produced) [€/MWh]: {best_row['lcoe_res_prod_eur_per_mwh']:.2f}")
            print(f"LCOE (used)     [€/MWh]: {best_row['lcoe_res_used_eur_per_mwh']:.2f}")

        out_name = f"{scenario_name}_res_to_h2_sweep.csv"
        out_path = out_dir / out_name
        results_df.to_csv(out_path)
        print(f"\nFull sweep saved to: {out_path}")


if __name__ == "__main__":
    main()