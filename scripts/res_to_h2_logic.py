import numpy as np
import pandas as pd


def annuity_factor(wacc: float, lifetime_years: float) -> float:
    """
    Capital recovery factor for annualising CAPEX.
    If this is also needed for grid LCOE etc., consider moving it
    into a shared module (e.g. lcoe_core.py) and importing from there.
    """
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
    """
    Annualised cost in €/year for an asset with given CAPEX/OPEX.
    Formula: annual_cost = CAPEX * CRF + OPEX_annual
    """
    crf = annuity_factor(wacc, lifetime_years)
    return capex_total * crf + opex_annual


def simulate_res_and_electrolyser(
    res_profile: np.ndarray,
    res_capacity_mw: float,
    electrolyser_mw: float,
    delta_hours: float,
) -> tuple[float, float]:
    """
    Given:
      - res_profile: per-unit profile over time (e.g. p.u. of installed RES capacity)
      - res_capacity_mw: candidate installed RES capacity [MW] we are testing
      - electrolyser_mw: electrolyser rated power [MW]
      - delta_hours: time step length [h]

    Returns:
      - elec_used_mwh: annual electricity actually consumed by the electrolyser [MWh/year]
      - curtailed_mwh: annual curtailed RES energy [MWh/year]
    """
    # Convert profile in p.u. to MW for this candidate RES capacity
    res_power_mw = res_profile * res_capacity_mw

    # Electrolyser can only take up to its MW rating
    elec_used_mw = np.minimum(res_power_mw, electrolyser_mw)
    curtailed_mw = np.maximum(res_power_mw - electrolyser_mw, 0.0)

    # Integrate over the year
    elec_used_mwh = elec_used_mw.sum() * delta_hours
    curtailed_mwh = curtailed_mw.sum() * delta_hours

    return elec_used_mwh, curtailed_mwh


def annual_cost_res(
    res_capacity_mw: float,
    res_cfg: dict,
    wacc: float,
) -> float:
    """
    Annualised RES cost [€/year] for a given installed RES capacity [MW].
    Uses RES-specific CAPEX/OPEX and lifetime from config.
    """
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
    """
    Annualised electrolyser cost [€/year] for a given electrolyser size [MW].
    Uses electrolyser-specific CAPEX/OPEX and lifetime from config.
    """
    capex = el_cfg["capex_per_mw_eur"] * electrolyser_mw
    opex = el_cfg["opex_per_mw_per_year_eur"] * electrolyser_mw
    lifetime = el_cfg["lifetime_years"]

    return annual_cost(
        capex_total=capex,
        opex_annual=opex,
        wacc=wacc,
        lifetime_years=lifetime,
    )


def optimise_res_to_el_ratio(
    res_profile: np.ndarray,
    config: dict,
    timestep_hours: float,
) -> tuple[pd.DataFrame, float]:
    """
    Sweep RES:electrolyser ratios and compute LCOH for each point.

    For now we assume a single scenario in config["scenarios"][0].
    Returns:
      - results: DataFrame indexed by res_to_el_ratio
      - best_ratio: ratio with minimum LCOH [€/kg]
    """
    scenario = config["scenarios"][0]  # MVP: single scenario
    res_cfg = scenario["res"]
    el_cfg = scenario["electrolyser"]
    opt_cfg = scenario["optimization"]

    # WACC: using global default from config for now
    wacc = config["finance"]["default_wacc"]


    min_ratio = opt_cfg["min_res_per_el_mw"]
    max_ratio = opt_cfg["max_res_per_el_mw"]
    steps = opt_cfg["steps"]

    # Base electrolyser size (MW) for the sweep
    electrolyser_mw = opt_cfg.get("electrolyser_mw", 10.0)

    ratios = np.linspace(min_ratio, max_ratio, steps)
    records: list[dict] = []

    for r in ratios:
        res_capacity_mw = r * electrolyser_mw

        elec_used_mwh, curtailed_mwh = simulate_res_and_electrolyser(
            res_profile=res_profile,
            res_capacity_mw=res_capacity_mw,
            electrolyser_mw=electrolyser_mw,
            delta_hours=timestep_hours,
        )

        # Convert electricity [MWh] to hydrogen [kg]
        eff_kwh_per_kg = el_cfg["efficiency_kwh_per_kg"]
        h2_kg = elec_used_mwh * 1000.0 / eff_kwh_per_kg

        if h2_kg <= 0:
            lcoh = np.inf
        else:
            cost_res = annual_cost_res(res_capacity_mw, res_cfg, wacc)
            cost_el = annual_cost_electrolyser(electrolyser_mw, el_cfg, wacc)
            annual_total_cost = cost_res + cost_el

            # LCOH = annual total cost / annual hydrogen production
            lcoh = annual_total_cost / h2_kg

        records.append(
            {
                "res_to_el_ratio": r,
                "res_capacity_mw": res_capacity_mw,
                "electrolyser_mw": electrolyser_mw,
                "elec_used_mwh": elec_used_mwh,
                "curtailed_mwh": curtailed_mwh,
                "h2_kg_per_year": h2_kg,
                "lcoh_eur_per_kg": lcoh,
            }
        )

    results = pd.DataFrame.from_records(records).set_index("res_to_el_ratio")
    best_ratio = results["lcoh_eur_per_kg"].idxmin()

    return results, best_ratio
