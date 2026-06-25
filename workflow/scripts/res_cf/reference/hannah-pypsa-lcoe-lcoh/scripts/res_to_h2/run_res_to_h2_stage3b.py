import numpy as np
import pandas as pd
from pathlib import Path
import yaml
from datetime import datetime, UTC

from res_to_h2_logic import dri_to_el_mw
from run_stage3b_sweep import run_stage3b_sweep

CONFIG_PATH = Path("config_hannah.yaml")
DATA_DIR = Path("data/res_cf")
VALID_RES_TECHS = {"wind_onshore", "solar", "wind_offshore"}
VALID_VARIANTS = {"avg", "bestsite_p95"}


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)


    if "defaults" not in config:
        raise KeyError("Config missing top-level key: 'defaults'")
    if "scenarios" not in config:
        raise KeyError("Config missing top-level key: 'scenarios'")
    if "finance" not in config:
        raise KeyError("Config missing top-level key: 'finance'")

    if not isinstance(config["scenarios"], list) or len(config["scenarios"]) == 0:
        raise ValueError("Config 'scenarios' must be a non-empty list.")

    for i, scenario in enumerate(config["scenarios"]):
        if not isinstance(scenario, dict):
            raise ValueError(f"Scenario at index {i} must be a dict. Got: {type(scenario)}")
        try:
            validate_scenario_config(scenario)
        except Exception as e:
            raise ValueError(f"Error in scenario index {i}: {e}") from e

    return config


def validate_scenario_config(scenario: dict) -> None:
    required_fields = [
        "profile_source",
        "country",
        "tech",
        "variant",
        "year",
        "mix_anchor_tech",
        "res_mix",
    ]

    for field in required_fields:
        if field not in scenario:
            raise ValueError(f"Scenario missing required field '{field}': {scenario}")

    tech = scenario["tech"]
    variant = scenario["variant"]
    mix_anchor_tech = scenario["mix_anchor_tech"]
    res_mix = scenario["res_mix"]

    if tech == "res_mix":
        if variant not in VALID_VARIANTS:
            raise ValueError(
                f"Unsupported variant '{variant}' for res_mix. Allowed: {sorted(VALID_VARIANTS)}"
            )

        if mix_anchor_tech not in VALID_RES_TECHS:
            raise ValueError(
                f"mix_anchor_tech must be one of {sorted(VALID_RES_TECHS)}. "
                f"Got: {mix_anchor_tech}"
            )

        if not isinstance(res_mix, dict) or len(res_mix) == 0:
            raise ValueError(f"res_mix must be a non-empty dict. Got: {res_mix}")

        for mix_tech, weight in res_mix.items():
            if mix_tech not in VALID_RES_TECHS:
                raise ValueError(
                    f"Unsupported res_mix tech '{mix_tech}'. Allowed: {sorted(VALID_RES_TECHS)}"
                )
            if not isinstance(weight, (int, float)):
                raise ValueError(
                    f"Weight for tech '{mix_tech}' must be numeric. Got: {weight}"
                )
            if weight < 0:
                raise ValueError(
                    f"Weight for tech '{mix_tech}' must be >= 0. Got: {weight}"
                )

        if mix_anchor_tech not in res_mix:
            raise ValueError(
                f"mix_anchor_tech '{mix_anchor_tech}' must appear in res_mix keys "
                f"{list(res_mix.keys())}"
            )

        weight_sum = sum(res_mix.values())
        if abs(weight_sum - 1.0) > 1e-6:
            raise ValueError(f"res_mix weights must sum to 1.0. Got: {weight_sum}")

    else:
        if tech not in VALID_RES_TECHS:
            raise ValueError(
                f"Unsupported single-tech '{tech}'. Allowed: {sorted(VALID_RES_TECHS)}"
            )

        if variant not in VALID_VARIANTS:
            raise ValueError(
                f"Unsupported variant '{variant}'. Allowed: {sorted(VALID_VARIANTS)}"
            )

        if mix_anchor_tech is not None:
            raise ValueError(
                f"For single-tech scenario '{tech}', mix_anchor_tech must be null. "
                f"Got: {mix_anchor_tech}"
            )

        if res_mix is not None:
            raise ValueError(
                f"For single-tech scenario '{tech}', res_mix must be null. "
                f"Got: {res_mix}"
            )


def load_atlite_profile(
    country: str,
    tech: str,
    variant: str,
    year: int,
    mix_anchor_tech=None,
    res_mix=None,
):
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

    if tech != "res_mix":
        col = f"{tech}_cf"
        if col not in df.columns:
            raise KeyError(f"{col} not found in {filename}")
        series = df[col].astype(float)
        return series.to_numpy(), 1.0  # hourly

    # res_mix case

    # For res_mix + bestsite_p95, load anchor-specific co-located CF file.
# These files are generated upstream from one anchor-defined location basis:
# - land anchors: land techs use the same land cell
# - offshore remains matched via deterministic counterpart logic upstream
    
    if variant == "bestsite_p95":
        if mix_anchor_tech not in {"wind_onshore", "solar", "wind_offshore"}:
            raise NotImplementedError(
                "res_mix with variant='bestsite_p95' currently supports only "
                "mix_anchor_tech in {'wind_onshore', 'solar', 'wind_offshore'}."
            )

        mix_label = format_res_mix_label(res_mix)

        anchor_filename = (
            f"{cc}_cf_{year}_bestsite_p95_"
            f"anchor-{mix_anchor_tech}_mix-{mix_label}.csv"
        )

        anchor_path = DATA_DIR / anchor_filename

        if not anchor_path.exists():
            raise FileNotFoundError(
                f"Missing scenario-specific anchor-tech CF file: {anchor_path}\n"
                "Run scripts/res_cf/07_make_bestsite_cf_timeseries.py first."
            )

        df_anchor = pd.read_csv(anchor_path, parse_dates=["time"])
        print(f"Loaded CF file: {anchor_filename}")

        cf_mix = None
        for mix_tech in sorted(res_mix.keys()):
            weight = res_mix[mix_tech]
            col = f"{mix_tech}_cf"
            if col not in df_anchor.columns:
                raise KeyError(f"{col} not found in {anchor_filename}")

            series = df_anchor[col].astype(float)
            weighted_series = float(weight) * series

            if cf_mix is None:
                cf_mix = weighted_series
            else:
                cf_mix = cf_mix + weighted_series
       # Guard against tiny floating-point drift outside physical CF bounds         
        cf_mix = cf_mix.clip(lower=0.0, upper=1.0)

        arr = cf_mix.to_numpy()
        if len(arr) != 8760:
            raise ValueError(f"Unexpected time series length: {len(arr)} (expected 8760)")

        return arr, 1.0  # hourly

    cf_mix = None
    for mix_tech in sorted(res_mix.keys()):
        weight = res_mix[mix_tech]
        col = f"{mix_tech}_cf"
        if col not in df.columns:
            raise KeyError(f"{col} not found in {filename}")

        series = df[col].astype(float)
        weighted = weight * series

        if cf_mix is None:
            cf_mix = weighted
        else:
            cf_mix = cf_mix + weighted
            
# Guard against tiny floating-point drift outside physical CF bounds
    cf_mix = cf_mix.clip(lower=0.0, upper=1.0)

    arr = cf_mix.to_numpy()
    if len(arr) != 8760:
        raise ValueError(f"Unexpected time series length: {len(arr)} (expected 8760)")

    return arr, 1.0  # hourly

def format_res_mix_label(res_mix: dict) -> str:
    parts = []
    for tech in sorted(res_mix.keys()):
        weight = res_mix[tech]
        weight_str = str(float(weight)).replace(".", "p")
        parts.append(f"{tech}-{weight_str}")
    return "_".join(parts)

def append_stage3b_comparison_row(
    out_dir: Path,
    name: str,
    scenario: dict,
    best_idx: float,
    best_row: pd.Series,
    electrolyser_mw: float,
) -> None:
    comparison_path = out_dir / "res_to_h2_comparison_table_stage3b.csv"

    if scenario["tech"] == "res_mix":
        res_mix_label = format_res_mix_label(scenario["res_mix"])
        mix_anchor_tech = scenario["mix_anchor_tech"] if scenario["variant"] == "bestsite_p95" else None
    else:
        res_mix_label = None
        mix_anchor_tech = None

    storage_mwh = float(best_row["storage_capacity_required_mwh"])
    storage_duration_h = storage_mwh / float(electrolyser_mw) if electrolyser_mw > 0 else np.nan

    row = {
        "scenario": name,
        "country": scenario["country"],
        "tech": scenario["tech"],
        "variant": scenario["variant"],
        "year": int(scenario["year"]),
        "res_mix": res_mix_label,
        "mix_anchor_tech": mix_anchor_tech,
        "electrolyser_mw": float(electrolyser_mw),
        "best_ratio_stage3b": float(best_idx),
        "best_lcoh_with_storage_eur_per_kg": float(best_row["lcoh_with_storage_eur_per_kg"]),
        "best_lcoh_no_storage_eur_per_kg": float(best_row["lcoh_eur_per_kg"]),
        "storage_capacity_required_mwh": storage_mwh,
        "storage_duration_hours": storage_duration_h,
        "storage_cost_contribution_eur_per_kg": float(best_row["lcoh_storage_component_with_storage_eur_per_kg"]),
        "res_cost_contribution_eur_per_kg": float(best_row["lcoh_res_component_with_storage_eur_per_kg"]),
        "electrolyser_cost_contribution_eur_per_kg": float(best_row["lcoh_electrolyser_component_with_storage_eur_per_kg"]),
        "curtailment_with_storage_mwh": float(best_row["curtailment_with_storage_mwh"]),
        "el_utilisation_at_best": float(best_row["el_utilisation"]),
        "energy_served_fraction_with_storage": float(best_row["energy_served_fraction_with_storage"]),
        "h2_kg_with_storage_per_year": float(best_row["h2_kg_with_storage_per_year"]),
        "h2_kg_no_storage_per_year": float(best_row["h2_kg_no_storage_per_year"]),
        "firm_res_electricity_cost_proxy_eur_per_mwh": float((best_row["annual_res_cost_eur"] + best_row["annual_storage_cost_eur"])/ best_row["elec_served_with_storage_mwh"]), #approximate cost of firm, RES-based electricity at this location/system design -> not yet fully endogeneous elecritcity cost for the whole DRI-EAF plant (is not sized for elec demand from EAF, DRI auxiliary, heat)
        "run_timestamp_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),

    }

    row_df = pd.DataFrame([row])

    if comparison_path.exists():
        existing = pd.read_csv(comparison_path)
        updated = pd.concat([existing, row_df], ignore_index=True)
    else:
        updated = row_df

    COLUMN_ORDER = [
        "scenario",
        "country",
        "tech",
        "variant",
        "year",
        "res_mix",
        "mix_anchor_tech",
        "electrolyser_mw",
        "best_ratio_stage3b",
        "best_lcoh_with_storage_eur_per_kg",
        "best_lcoh_no_storage_eur_per_kg",
        "res_cost_contribution_eur_per_kg",
        "electrolyser_cost_contribution_eur_per_kg",
        "storage_cost_contribution_eur_per_kg",
        "firm_res_electricity_cost_proxy_eur_per_mwh",
        "storage_capacity_required_mwh",
        "storage_duration_hours",
        "curtailment_with_storage_mwh",
        "el_utilisation_at_best",
        "energy_served_fraction_with_storage",
        "h2_kg_with_storage_per_year",
        "h2_kg_no_storage_per_year",
        "run_timestamp_utc",
    ]

    updated = updated[[c for c in COLUMN_ORDER if c in updated.columns]]

    updated.to_csv(comparison_path, index=False)

def main():
    config = load_config(CONFIG_PATH)

    defaults = config["defaults"]
    scenarios = config["scenarios"]
    finance = config["finance"]

    out_dir = Path("results")
    out_dir.mkdir(exist_ok=True)

    for s in scenarios:

        # merge defaults + scenario (simple version)
        scenario = {**defaults, **s}

        

        country = scenario["country"]
        tech = scenario["tech"]
        variant = scenario["variant"]
        year = int(scenario["year"])


        if tech == "res_mix":
            mix_label = format_res_mix_label(scenario["res_mix"])
            anchor_tech = scenario["mix_anchor_tech"]

            if variant == "bestsite_p95":
                name = f"{country}_resmix_anchor-{anchor_tech}_{mix_label}_{variant}_{year}_stage3b"
            else:
                name = f"{country}_resmix_{mix_label}_{variant}_{year}_stage3b"
        else:
            name = f"{country}_{tech}_{variant}_{year}_stage3b"

        print("\n==============================")
        print("Running:", name)
        print("==============================")

        # Load CF
        res_profile, timestep_hours = load_atlite_profile(
            country=country,
            tech=tech,
            variant=variant,
            year=year,
            mix_anchor_tech=scenario["mix_anchor_tech"],
            res_mix=scenario["res_mix"],
        )

        # Electrolyser sizing
        dri_cfg = scenario["dri_sizing"]

        el_mw = dri_to_el_mw(
            dri_mt_per_year=float(dri_cfg["dri_mt_per_year"]),
            h2_intensity_kg_per_t_dri=float(dri_cfg["h2_intensity_kg_per_t_dri"]),
            efficiency_kwh_per_kg=float(scenario["electrolyser"]["efficiency_kwh_per_kg"]),
            availability_target=float(dri_cfg["availability_target"]),
        )

        # Run Stage 3B sweep
        results = run_stage3b_sweep(
            res_profile=res_profile,
            scenario=scenario,
            finance=finance,
            timestep_hours=timestep_hours,
            electrolyser_mw=el_mw,
        )

        # --- Scenario metadata (for traceability) ---
        results["country"] = country
        results["tech"] = tech
        results["variant"] = variant
        results["year"] = year


        if tech == "res_mix":
            mix_label = format_res_mix_label(scenario["res_mix"])
            results["res_mix"] = mix_label

            if variant == "bestsite_p95":
                results["mix_anchor_tech"] = scenario["mix_anchor_tech"]
            else:
                results["mix_anchor_tech"] = None
        else:
            results["res_mix"] = None
            results["mix_anchor_tech"] = None

        # Keep only feasible rows
        feasible = results[results["storage_feasible"] == True]

        if feasible.empty:
            print("No feasible solution found (constraint too strict).")
            best_row = None
        else:
            best_idx = feasible["lcoh_with_storage_eur_per_kg"].idxmin()
            best_row = feasible.loc[best_idx]

            print("\nBest ratio:", best_idx)
            print("LCOH (with storage):", best_row["lcoh_with_storage_eur_per_kg"])
            print("Storage required [MWh]:", best_row["storage_capacity_required_mwh"])

            storage_mwh = best_row["storage_capacity_required_mwh"]
            storage_duration_h = storage_mwh / el_mw

            print(f"Storage duration [hours]: {storage_duration_h:,.1f}")
            print(f"Unserved energy [MWh]: {best_row['unserved_energy_mwh']:,.1f}")
            print(f"Curtailment with storage [MWh]: {best_row['curtailment_with_storage_mwh']:,.1f}")
            print(f"Energy served fraction with storage [-]: {best_row['energy_served_fraction_with_storage']:.4f}")
            print(f"Elec served with storage [MWh]: {best_row['elec_served_with_storage_mwh']:,.1f}")
            print(f"H2 with storage [kg/yr]: {best_row['h2_kg_with_storage_per_year']:,.0f}")
            print(f"H2 no storage [kg/yr]: {best_row['h2_kg_no_storage_per_year']:,.0f}")

            append_stage3b_comparison_row(
                out_dir=out_dir,
                name=name,
                scenario=scenario,
                best_idx=best_idx,
                best_row=best_row,
                electrolyser_mw=el_mw,
            )

        # Save full sweep
        out_path = out_dir / f"{name}_sweep.csv"
        results.to_csv(out_path)
        print("Saved:", out_path)


if __name__ == "__main__":
    main()