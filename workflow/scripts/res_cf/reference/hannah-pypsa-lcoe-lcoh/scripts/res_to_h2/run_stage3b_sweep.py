import numpy as np
import pandas as pd

from res_to_h2_logic import (
    annual_cost_res,
    annual_cost_res_mix,
    annual_cost_electrolyser,
)
from stage3b_storage_dispatch import find_min_storage_capacity


def run_stage3b_sweep(
    res_profile: np.ndarray,
    scenario: dict,
    finance: dict,
    timestep_hours: float,
    electrolyser_mw: float,
) -> pd.DataFrame:
    """
    Stage 3B sweep:
    - sweep across RES:EL ratios
    - simulate direct RES -> electrolyser operation
    - size battery chronologically
    - add storage cost
    - return one row per ratio
    """

    res_cfg = scenario["res"]
    el_cfg = scenario["electrolyser"]
    opt_cfg = scenario["optimization"]
    storage_cfg = scenario["storage_stage3b"]

    wacc = float(finance["default_wacc"])

    min_ratio = float(opt_cfg["min_res_per_el_mw"])
    max_ratio = float(opt_cfg["max_res_per_el_mw"])
    steps = int(opt_cfg["steps"])

    eta_rt = float(storage_cfg["roundtrip_efficiency"])
    allowed_fraction = float(storage_cfg["allowed_energy_deficit_fraction"])
    initial_soc_mode = str(storage_cfg.get("initial_soc_mode", "empty"))

    ratios = np.linspace(min_ratio, max_ratio, steps)
    records = []

    eff_kwh_per_kg = float(el_cfg["efficiency_kwh_per_kg"])

    annual_baseload_mwh = float(electrolyser_mw) * 8760.0
    allowed_uncovered_mwh = allowed_fraction * annual_baseload_mwh

    cf_std = float(np.std(res_profile))
    cf_p5 = float(np.percentile(res_profile, 5))
    cf_p1 = float(np.percentile(res_profile, 1))

    for r in ratios:
        res_capacity_mw = float(r) * float(electrolyser_mw)

        # Hourly RES power
        res_power_mw = res_profile * res_capacity_mw

        # No-storage direct use
        elec_used_mw = np.minimum(res_power_mw, electrolyser_mw)
        curtailed_mw = np.maximum(res_power_mw - electrolyser_mw, 0.0)

        elec_used_mwh = float(elec_used_mw.sum() * timestep_hours)
        curtailed_mwh = float(curtailed_mw.sum() * timestep_hours)

        h2_kg_no_storage = (
            elec_used_mwh * 1000.0 / eff_kwh_per_kg
            if elec_used_mwh > 0
            else 0.0
        )
        # Deficit / surplus for storage model
        deficit_mw = np.maximum(electrolyser_mw - res_power_mw, 0.0)
        surplus_mw = np.maximum(res_power_mw - electrolyser_mw, 0.0)

        deficit_mwh = deficit_mw * timestep_hours
        surplus_mwh = surplus_mw * timestep_hours

        annual_deficit_mwh = float(deficit_mwh.sum())

        # Stage 3B battery sizing
        storage_result = find_min_storage_capacity(
            surplus_mwh=surplus_mwh,
            deficit_mwh=deficit_mwh,
            allowed_uncovered_mwh=allowed_uncovered_mwh,
            eta_rt=eta_rt,
            initial_soc_mode=initial_soc_mode,
            tolerance_mwh=float(storage_cfg.get("tolerance_mwh", 1.0)),
            max_iter=int(storage_cfg.get("max_iter", 60)),
        )

        storage_feasible = bool(storage_result["feasible"])
        storage_capacity_required_mwh = float(storage_result["storage_capacity_required_mwh"]) if storage_feasible else np.nan
        unserved_energy_mwh = float(storage_result["unserved_energy_mwh"])
        curtailment_with_storage_mwh = float(storage_result["curtailment_mwh"])
        max_soc_mwh = float(storage_result["max_soc_mwh"])
        annual_storage_energy_throughput_mwh = float(storage_result["throughput_mwh"])

                # Electricity actually served to electrolyser after storage dispatch
        elec_served_with_storage_mwh = max(
            0.0,
            annual_baseload_mwh - unserved_energy_mwh
        )

        h2_kg_with_storage = (
            elec_served_with_storage_mwh * 1000.0 / eff_kwh_per_kg
            if elec_served_with_storage_mwh > 0
            else 0.0
        )

        # Costs
        if scenario["tech"] == "res_mix":
            annual_res_cost_eur = annual_cost_res_mix(res_capacity_mw, scenario, wacc)
        else:
            annual_res_cost_eur = annual_cost_res(res_capacity_mw, res_cfg, wacc)
        annual_electrolyser_cost_eur = annual_cost_electrolyser(electrolyser_mw, el_cfg, wacc)

        if storage_feasible and np.isfinite(storage_capacity_required_mwh):
            capex_per_kwh = float(storage_cfg["capex_eur_per_kwh"])
            opex_fraction = float(storage_cfg.get("fixed_opex_fraction_of_capex", 0.0))
            lifetime_years = float(storage_cfg["lifetime_years"])

            capex_storage_eur = storage_capacity_required_mwh * 1000.0 * capex_per_kwh
            opex_storage_eur_per_year = capex_storage_eur * opex_fraction

            # same CRF logic as elsewhere
            if wacc == 0:
                crf = 1.0 / lifetime_years
            else:
                crf = (wacc * (1.0 + wacc) ** lifetime_years) / (
                    (1.0 + wacc) ** lifetime_years - 1.0
                )

            annual_storage_cost_eur = capex_storage_eur * crf + opex_storage_eur_per_year
        else:
            annual_storage_cost_eur = np.nan

                # LCOH
        if h2_kg_no_storage > 0:
            lcoh_eur_per_kg = (
                annual_res_cost_eur + annual_electrolyser_cost_eur
            ) / h2_kg_no_storage
        else:
            lcoh_eur_per_kg = np.inf

        if (
            storage_feasible
            and np.isfinite(annual_storage_cost_eur)
            and h2_kg_with_storage > 0
        ):
            lcoh_with_storage_eur_per_kg = (
                annual_res_cost_eur
                + annual_electrolyser_cost_eur
                + annual_storage_cost_eur
            ) / h2_kg_with_storage
        else:
            lcoh_with_storage_eur_per_kg = np.inf

        # Energy-served fractions
        energy_served_fraction_no_storage = float(
            np.clip(elec_used_mwh / annual_baseload_mwh, 0.0, 1.0)
        )

        energy_served_fraction_with_storage = float(
            np.clip(elec_served_with_storage_mwh / annual_baseload_mwh, 0.0, 1.0)
        )

        records.append(
            {
                "res_to_el_ratio": float(r),
                "res_capacity_mw": res_capacity_mw,
                "electrolyser_mw": float(electrolyser_mw),

                # Electricity / hydrogen
                "elec_used_mwh": elec_used_mwh,
                "elec_served_no_storage_mwh": elec_used_mwh,
                "elec_served_with_storage_mwh": elec_served_with_storage_mwh,
                "curtailed_mwh": curtailed_mwh,
                "curtailment_with_storage_mwh": curtailment_with_storage_mwh,
                "h2_kg_per_year": h2_kg_no_storage,
                "h2_kg_no_storage_per_year": h2_kg_no_storage,
                "h2_kg_with_storage_per_year": h2_kg_with_storage,

                # Costs / LCOH
                "annual_res_cost_eur": annual_res_cost_eur,
                "annual_electrolyser_cost_eur": annual_electrolyser_cost_eur,
                "annual_storage_cost_eur": annual_storage_cost_eur,
                "lcoh_eur_per_kg": lcoh_eur_per_kg,
                "lcoh_with_storage_eur_per_kg": lcoh_with_storage_eur_per_kg,

                # Storage diagnostics
                "annual_deficit_mwh": annual_deficit_mwh,
                "allowed_uncovered_mwh": allowed_uncovered_mwh,
                "storage_feasible": storage_feasible,
                "storage_capacity_required_mwh": storage_capacity_required_mwh,
                "unserved_energy_mwh": unserved_energy_mwh,
                "max_soc_mwh": max_soc_mwh,
                "annual_storage_energy_throughput_mwh": annual_storage_energy_throughput_mwh,

                # CF diagnostics
                "cf_std": cf_std,
                "cf_p5": cf_p5,
                "cf_p1": cf_p1,

                # Served-energy fractions
                "energy_served_fraction_no_storage": energy_served_fraction_no_storage,
                "energy_served_fraction_with_storage": energy_served_fraction_with_storage,
            }
        )

    results = pd.DataFrame.from_records(records).set_index("res_to_el_ratio")

    results["res_generation_mwh"] = results["elec_used_mwh"] + results["curtailed_mwh"]
    results["lcoe_res_prod_eur_per_mwh"] = results["annual_res_cost_eur"] / results["res_generation_mwh"]
    results["lcoe_res_used_eur_per_mwh"] = results["annual_res_cost_eur"] / results["elec_used_mwh"]
    results["el_full_load_hours"] = results["elec_used_mwh"] / results["electrolyser_mw"]
    results["el_utilisation"] = results["el_full_load_hours"] / 8760.0
    results["curtailment_rate"] = results["curtailed_mwh"] / results["res_generation_mwh"]
    results["curtailment_rate_with_storage"] = results["curtailment_with_storage_mwh"] / results["res_generation_mwh"]
        # LCOH component contributions using internally consistent denominators
    results["lcoh_res_component_no_storage_eur_per_kg"] = (
        results["annual_res_cost_eur"] / results["h2_kg_no_storage_per_year"]
    )
    results["lcoh_electrolyser_component_no_storage_eur_per_kg"] = (
        results["annual_electrolyser_cost_eur"] / results["h2_kg_no_storage_per_year"]
    )

    results["lcoh_res_component_with_storage_eur_per_kg"] = (
        results["annual_res_cost_eur"] / results["h2_kg_with_storage_per_year"]
    )
    results["lcoh_electrolyser_component_with_storage_eur_per_kg"] = (
        results["annual_electrolyser_cost_eur"] / results["h2_kg_with_storage_per_year"]
    )
    results["lcoh_storage_component_with_storage_eur_per_kg"] = (
        results["annual_storage_cost_eur"] / results["h2_kg_with_storage_per_year"]
    )
    results.replace([np.inf, -np.inf], np.nan, inplace=True)

    return results