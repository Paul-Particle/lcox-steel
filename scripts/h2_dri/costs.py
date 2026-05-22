"""
Post-solve LCOH and summary metrics from a solved PyPSA network.

Uses explicit cost accounting rather than n.statistics() to be robust
across PyPSA versions and to always include the fixed electrolyser capex.
"""

import pandas as pd
import pypsa

H2_LHV_MWH_PER_KG = 33.33 / 1000  # LHV of hydrogen in MWh/kg


def compute_lcoh(n: pypsa.Network) -> float:
    """LCOH in €/kg H2."""
    return _annual_cost(n) / _h2_produced_kg(n)


def extract_summary(n: pypsa.Network, project_name: str, scenario_name: str) -> dict:
    """Key sizing and cost metrics as a flat dict (suitable for a one-row CSV)."""
    summary = {
        "project": project_name,
        "scenario": scenario_name,
        "lcoh_eur_per_kg": compute_lcoh(n),
        "total_annual_cost_eur": _annual_cost(n),
        "h2_produced_kg": _h2_produced_kg(n),
    }

    # Optimal generator capacities (extendable only)
    for gen in n.generators.index[n.generators.p_nom_extendable]:
        summary[f"{gen}_mw_opt"] = n.generators.at[gen, "p_nom_opt"]

    # Battery
    if "battery" in n.storage_units.index:
        p_opt = n.storage_units.at["battery", "p_nom_opt"]
        summary["battery_mw_opt"] = p_opt
        summary["battery_mwh_opt"] = p_opt * n.storage_units.at["battery", "max_hours"]

    # H2 buffer
    if "h2_buffer" in n.stores.index:
        summary["h2_buffer_mwh_lhv_opt"] = n.stores.at["h2_buffer", "e_nom_opt"]

    # Electrolyser (fixed)
    if "electrolyser" in n.links.index:
        summary["electrolyser_mw"] = n.links.at["electrolyser", "p_nom"]

    return summary


def _annual_cost(n: pypsa.Network) -> float:
    """Sum of all annualized capital costs and variable grid import costs."""
    cost = 0.0

    # Extendable generators (wind, solar)
    mask = n.generators.p_nom_extendable
    cost += (n.generators.loc[mask, "capital_cost"] * n.generators.loc[mask, "p_nom_opt"]).sum()

    # Extendable storage units (battery)
    mask = n.storage_units.p_nom_extendable
    cost += (n.storage_units.loc[mask, "capital_cost"] * n.storage_units.loc[mask, "p_nom_opt"]).sum()

    # Electrolyser — fixed p_nom, always include its capex
    cost += n.links.at["electrolyser", "capital_cost"] * n.links.at["electrolyser", "p_nom"]

    # Extendable stores (H2 buffer)
    mask = n.stores.e_nom_extendable
    cost += (n.stores.loc[mask, "capital_cost"] * n.stores.loc[mask, "e_nom_opt"]).sum()

    # Grid import variable cost — price * dispatch summed over all hours
    if "grid_import" in n.generators.index:
        p = n.generators_t.p.get("grid_import", pd.Series(0.0, index=n.snapshots))
        # marginal_cost may be time-varying (in generators_t) or a scalar (in generators)
        if "grid_import" in n.generators_t.marginal_cost.columns:
            mc = n.generators_t.marginal_cost["grid_import"]
        else:
            mc = n.generators.at["grid_import", "marginal_cost"]
        cost += float((p * mc).sum())

    return float(cost)


def _h2_produced_kg(n: pypsa.Network) -> float:
    """H2 consumed by the DRI load over the simulation period, in kg."""
    h2_mwh_lhv = float(n.loads_t.p["dri_load"].sum())
    return h2_mwh_lhv / H2_LHV_MWH_PER_KG
