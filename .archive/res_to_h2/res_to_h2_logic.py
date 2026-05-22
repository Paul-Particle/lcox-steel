import numpy as np
import pandas as pd


def annuity_factor(wacc: float, lifetime_years: float) -> float:
    if lifetime_years <= 0:
        raise ValueError("lifetime_years must be > 0")

    if wacc == 0:
        return 1.0 / lifetime_years

    return (wacc * (1.0 + wacc) ** lifetime_years) / (
        (1.0 + wacc) ** lifetime_years - 1.0
    )


def annual_cost(
    capex_total: float,
    opex_annual: float,
    wacc: float,
    lifetime_years: float,
) -> float:
    crf = annuity_factor(wacc, lifetime_years)
    return capex_total * crf + opex_annual


def simulate_res_and_electrolyser(
    res_profile: np.ndarray,
    res_capacity_mw: float,
    electrolyser_mw: float,
    delta_hours: float,
) -> tuple[float, float]:

    res_power_mw = res_profile * res_capacity_mw
    elec_used_mw = np.minimum(res_power_mw, electrolyser_mw)
    curtailed_mw = np.maximum(res_power_mw - electrolyser_mw, 0.0)

    elec_used_mwh = elec_used_mw.sum() * delta_hours
    curtailed_mwh = curtailed_mw.sum() * delta_hours

    return elec_used_mwh, curtailed_mwh


def annual_cost_res(
    res_capacity_mw: float,
    res_cfg: dict,
    wacc: float,
) -> float:

    capex = res_cfg["capex_per_mw_eur"] * res_capacity_mw
    opex = res_cfg["opex_per_mw_per_year_eur"] * res_capacity_mw
    lifetime = res_cfg["lifetime_years"]

    return annual_cost(
        capex_total=capex,
        opex_annual=opex,
        wacc=wacc,
        lifetime_years=lifetime,
    )


def annual_cost_electrolyser(
    electrolyser_mw: float,
    el_cfg: dict,
    wacc: float,
) -> float:

    capex = el_cfg["capex_per_mw_eur"] * electrolyser_mw
    opex = el_cfg["opex_per_mw_per_year_eur"] * electrolyser_mw
    lifetime = el_cfg["lifetime_years"]

    return annual_cost(
        capex_total=capex,
        opex_annual=opex,
        wacc=wacc,
        lifetime_years=lifetime,
    )


def dri_to_el_mw(
    dri_mt_per_year: float,
    h2_intensity_kg_per_t_dri: float,
    efficiency_kwh_per_kg: float,
    availability_target: float,
) -> float:

    if dri_mt_per_year <= 0:
        raise ValueError("dri_mt_per_year must be > 0")
    if h2_intensity_kg_per_t_dri <= 0:
        raise ValueError("h2_intensity_kg_per_t_dri must be > 0")
    if efficiency_kwh_per_kg <= 0:
        raise ValueError("efficiency_kwh_per_kg must be > 0")
    if not (0 < availability_target <= 1.0):
        raise ValueError("availability_target must be in (0, 1]")

    dri_t_per_year = dri_mt_per_year * 1_000_000.0
    h2_kg_per_year = dri_t_per_year * h2_intensity_kg_per_t_dri
    elec_mwh_per_year = h2_kg_per_year * efficiency_kwh_per_kg / 1000.0

    hours_available = 8760.0 * availability_target
    el_mw_required = elec_mwh_per_year / hours_available

    return el_mw_required


def optimise_res_to_el_ratio(
    res_profile: np.ndarray,
    scenario: dict,
    finance: dict,
    timestep_hours: float,
    electrolyser_mw: float | None = None,
) -> tuple[pd.DataFrame, float]:

    res_cfg = scenario["res"]
    el_cfg = scenario["electrolyser"]
    opt_cfg = scenario["optimization"]
    storage_cfg = scenario.get("storage_proxy", {})

    wacc = float(finance["default_wacc"])

    min_ratio = float(opt_cfg["min_res_per_el_mw"])
    max_ratio = float(opt_cfg["max_res_per_el_mw"])
    steps = int(opt_cfg["steps"])

    if electrolyser_mw is None:
        electrolyser_mw = float(opt_cfg.get("electrolyser_mw", 10.0))
    else:
        electrolyser_mw = float(electrolyser_mw)

    ratios = np.linspace(min_ratio, max_ratio, steps)
    records: list[dict] = []

    for r in ratios:
        res_capacity_mw = r * electrolyser_mw

        # Hourly power
        res_power_mw = res_profile * res_capacity_mw

        # Firm reliability diagnostic (hour-based)
        firm_reliability = float(np.mean(res_power_mw >= electrolyser_mw))

        # Simulate basic interaction
        elec_used_mwh, curtailed_mwh = simulate_res_and_electrolyser(
            res_profile=res_profile,
            res_capacity_mw=res_capacity_mw,
            electrolyser_mw=electrolyser_mw,
            delta_hours=timestep_hours,
        )

        eff_kwh_per_kg = float(el_cfg["efficiency_kwh_per_kg"])
        h2_kg = elec_used_mwh * 1000.0 / eff_kwh_per_kg

        cost_res = np.nan
        cost_el = np.nan

        annual_storage_cost = 0.0
        storage_capacity_required = 0.0
        annual_deficit_mwh = 0.0

        # NEW: duration-cap diagnostics
        storage_capacity_required_uncapped_mwh = 0.0
        storage_duration_cap_mwh = 0.0
        storage_duration_cap_binding_bool = False
        annual_unserved_mwh_after_storage_cap = 0.0

        energy_served_fraction_no_storage = np.nan
        energy_served_fraction_proxy = np.nan
        energy_unserved_mwh_proxy = np.nan

        if h2_kg <= 0:
            lcoh = np.inf
            lcoh_with_storage = np.inf
        else:
            cost_res = annual_cost_res(res_capacity_mw, res_cfg, wacc)
            cost_el = annual_cost_electrolyser(electrolyser_mw, el_cfg, wacc)
            annual_total_cost = cost_res + cost_el
            lcoh = annual_total_cost / h2_kg

            # -------------------------------
            # Stage 3A: Storage Proxy
            # -------------------------------
            if storage_cfg.get("enabled", False):

                # Hourly deficit & surplus
                deficit_mw = np.maximum(electrolyser_mw - res_power_mw, 0.0)
                surplus_mw = np.maximum(res_power_mw - electrolyser_mw, 0.0)

                deficit_mwh = deficit_mw * timestep_hours
                surplus_mwh = surplus_mw * timestep_hours

                annual_deficit_mwh = float(deficit_mwh.sum())

                annual_baseload_mwh = electrolyser_mw * 8760.0
                allowed_fraction = float(storage_cfg["allowed_energy_deficit_fraction"])
                allowed_uncovered_mwh = allowed_fraction * annual_baseload_mwh

                # Energy-served metrics (annual, energy-based)
                energy_served_fraction_no_storage = float(np.clip(elec_used_mwh / annual_baseload_mwh, 0.0, 1.0))


                if storage_cfg.get("enabled", False):
                    # Proxy assumption: storage covers all deficit energy except the allowed uncovered slice
                    energy_unserved_mwh_proxy = float(min(annual_deficit_mwh, allowed_uncovered_mwh))
                    energy_served_fraction_proxy = float(1.0 - (energy_unserved_mwh_proxy / annual_baseload_mwh))
                else:
                    energy_unserved_mwh_proxy = np.nan
                    energy_served_fraction_proxy = np.nan

                if annual_deficit_mwh > allowed_uncovered_mwh:

                    # Net energy trajectory (proxy)
                    net_balance = surplus_mwh - deficit_mwh
                    cumulative = np.cumsum(net_balance)

                    full_storage_required = float(cumulative.max() - cumulative.min())

                    # Energy that must be covered by shifting (MWh/year)
                    energy_to_cover = annual_deficit_mwh - allowed_uncovered_mwh

                    # Scale storage energy range to the fraction of deficit we require to cover
                    scaling_factor = energy_to_cover / annual_deficit_mwh
                    storage_capacity_required = full_storage_required * scaling_factor

                    # Efficiency adjustment
                    eta = float(storage_cfg["roundtrip_efficiency"])
                    storage_capacity_required = storage_capacity_required / eta

                    # Track uncapped requirement (post-eta) for diagnostics
                    storage_capacity_required_uncapped_mwh = float(storage_capacity_required)

                    # --- Duration cap (battery proxy) ---
                    max_duration = storage_cfg.get("max_storage_duration_hours", None)
                    if max_duration is not None:
                        storage_duration_cap_mwh = float(electrolyser_mw) * float(max_duration)

                        if storage_capacity_required > storage_duration_cap_mwh:
                            storage_duration_cap_binding_bool = True
                            storage_capacity_required = storage_duration_cap_mwh

                        # If cap binds, compute residual unserved energy (lower bound proxy)
                        # Deliverable shifted energy is limited by eta * E_cap
                        deliverable_mwh = eta * storage_capacity_required
                        annual_unserved_mwh_after_storage_cap = float(max(0.0, energy_to_cover - deliverable_mwh))

                    # Costing
                    capex_per_kwh = float(storage_cfg["capex_eur_per_kwh"])
                    lifetime = float(storage_cfg["lifetime_years"])

                    capex_storage = storage_capacity_required * 1000.0 * capex_per_kwh
                    crf = (wacc * (1 + wacc) ** lifetime) / ((1 + wacc) ** lifetime - 1)
                    annual_storage_cost = capex_storage * crf

                lcoh_with_storage = (cost_res + cost_el + annual_storage_cost) / h2_kg

            else:
                lcoh_with_storage = lcoh

        records.append(
            {
                "res_to_el_ratio": r,
                "res_capacity_mw": res_capacity_mw,
                "electrolyser_mw": electrolyser_mw,
                "elec_used_mwh": elec_used_mwh,
                "curtailed_mwh": curtailed_mwh,
                "h2_kg_per_year": h2_kg,
                "lcoh_eur_per_kg": lcoh,
                "lcoh_with_storage_eur_per_kg": lcoh_with_storage,
                "annual_res_cost_eur": cost_res,
                "annual_electrolyser_cost_eur": cost_el,
                "annual_storage_cost_eur": annual_storage_cost,
                "annual_deficit_mwh": annual_deficit_mwh,
                "storage_capacity_required_mwh": storage_capacity_required,
                # NEW: duration-cap diagnostics
                "storage_capacity_required_uncapped_mwh": storage_capacity_required_uncapped_mwh,
                "storage_duration_cap_mwh": storage_duration_cap_mwh,
                "storage_duration_cap_binding_bool": storage_duration_cap_binding_bool,
                "annual_unserved_mwh_after_storage_cap": annual_unserved_mwh_after_storage_cap,
                "firm_reliability": firm_reliability,
                "energy_served_fraction_no_storage": energy_served_fraction_no_storage,
                "energy_served_fraction_proxy": energy_served_fraction_proxy,
                "energy_unserved_mwh_proxy": energy_unserved_mwh_proxy,
            }
        )

    results = pd.DataFrame.from_records(records).set_index("res_to_el_ratio")

    results["res_generation_mwh"] = results["elec_used_mwh"] + results["curtailed_mwh"]
    results["lcoe_res_prod_eur_per_mwh"] = results["annual_res_cost_eur"] / results["res_generation_mwh"]
    results["lcoe_res_used_eur_per_mwh"] = results["annual_res_cost_eur"] / results["elec_used_mwh"]

    results["el_full_load_hours"] = results["elec_used_mwh"] / results["electrolyser_mw"]
    results["el_utilisation"] = results["el_full_load_hours"] / 8760.0
    results["curtailment_rate"] = results["curtailed_mwh"] / results["res_generation_mwh"]

    best_ratio = results["lcoh_eur_per_kg"].idxmin()

    return results, best_ratio


def firm_ratio_from_sweep(
    results: pd.DataFrame,
    target_reliability: float,
) -> dict:

    if "firm_reliability" not in results.columns:
        raise KeyError("results must include column 'firm_reliability'")

    if not (0 < target_reliability <= 1.0):
        raise ValueError("target_reliability must be in (0, 1]")

    s = results["firm_reliability"].sort_index()
    feasible = s[s >= target_reliability]

    if feasible.empty:
        return {
            "firm_ratio_interp": np.nan,
            "firm_ratio_feasible": np.nan,
            "firm_row_feasible": None,
        }

    r_hi = float(feasible.index[0])
    rel_hi = float(s.loc[r_hi])

    idx = list(s.index)
    hi_pos = idx.index(r_hi)

    if hi_pos == 0:
        firm_ratio_interp = r_hi
    else:
        r_lo = float(idx[hi_pos - 1])
        rel_lo = float(s.loc[r_lo])

        if rel_hi == rel_lo:
            firm_ratio_interp = r_hi
        else:
            firm_ratio_interp = r_lo + (target_reliability - rel_lo) * (r_hi - r_lo) / (rel_hi - rel_lo)

    return {
        "firm_ratio_interp": float(firm_ratio_interp),
        "firm_ratio_feasible": r_hi,
        "firm_row_feasible": results.loc[r_hi],
    }
