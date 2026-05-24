"""
Build a PyPSA network for a single DRI-hydrogen project scenario.

Bus unit convention: MW throughout.
  electricity bus: MW AC
  hydrogen bus:    MW H2 LHV  (1 MWh H2 LHV ≈ 30 kg H2 at LHV ≈ 33.33 kWh/kg)

Electrolyser efficiency is:
  efficiency = h2_lhv_kwh_per_kg / efficiency_kwh_per_kg
             = 33.33 / 55 ≈ 0.606 (MW H2 LHV per MW electricity)
"""

import sys
from pathlib import Path

import pandas as pd
import pypsa

sys.path.insert(0, str(Path(__file__).parent))
from sizing import annuity_factor, dri_to_el_mw


def build_network(
    project_cfg: dict,
    assumptions: dict,
    cf_timeseries: pd.DataFrame,
    price_series: pd.Series | None = None,
) -> pypsa.Network:
    """
    Build (but do not solve) a PyPSA network.

    project_cfg: one entry from projects.yaml (with 'plant' and 'scenarios' resolved)
    assumptions: full assumptions.yaml dict
    cf_timeseries: DataFrame indexed by DatetimeIndex, columns = tech names in assumptions.res
    price_series: optional hourly price Series (€/MWh), same index as cf_timeseries
    """
    wacc = assumptions["finance"]["default_wacc"]
    el_cfg = assumptions["electrolyser"]
    plant = project_cfg["plant"]
    h2_lhv_kwh_per_kg = assumptions["h2"]["lhv_kwh_per_kg"]

    el_mw = dri_to_el_mw(
        dri_mt_per_year=plant["dri_mt_per_year"],
        h2_intensity_kg_per_t_dri=plant["h2_intensity_kg_per_t_dri"],
        efficiency_kwh_per_kg=el_cfg["efficiency_kwh_per_kg"],
        availability_target=plant["availability_target"],
    )
    el_efficiency = h2_lhv_kwh_per_kg / el_cfg["efficiency_kwh_per_kg"]

    n = pypsa.Network()
    n.set_snapshots(cf_timeseries.index)

    _add_buses(n)
    _add_generators(n, cf_timeseries, assumptions["res"], wacc)
    _add_battery(n, assumptions["battery"], wacc)
    _add_electrolyser(n, el_mw, el_efficiency, el_cfg, wacc)
    _add_h2_buffer(n, assumptions["h2_buffer"], wacc)
    _add_dri_load(n, el_mw, el_efficiency)

    if price_series is not None:
        _add_grid_import(n, price_series)

    return n


def _add_buses(n: pypsa.Network) -> None:
    n.add("Bus", "electricity", carrier="AC")
    n.add("Bus", "hydrogen", carrier="H2")


def _add_generators(
    n: pypsa.Network,
    cf_timeseries: pd.DataFrame,
    res_cfg: dict,
    wacc: float,
) -> None:
    for tech in cf_timeseries.columns:
        if tech not in res_cfg:
            raise KeyError(f"No assumptions found for tech '{tech}' — add it to assumptions.yaml")
        cfg = res_cfg[tech]
        cap_cost = (
            annuity_factor(wacc, cfg["lifetime_years"]) * cfg["capex_per_mw_eur"]
            + cfg["opex_per_mw_per_year_eur"]
        )
        n.add(
            "Generator",
            tech,
            bus="electricity",
            carrier=tech,
            p_nom_extendable=True,
            capital_cost=cap_cost,
            marginal_cost=0.0,
            p_max_pu=cf_timeseries[tech],
        )


def _add_battery(n: pypsa.Network, bat_cfg: dict, wacc: float) -> None:
    eta = bat_cfg["efficiency_roundtrip"] ** 0.5
    max_hours = bat_cfg["max_hours"]
    # Fold energy capex into per-MW capital cost (assumes fixed duration = max_hours)
    cap_cost = annuity_factor(wacc, bat_cfg["lifetime_years"]) * (
        bat_cfg["capex_per_mw_eur"] + bat_cfg["capex_per_mwh_eur"] * max_hours
    )
    n.add(
        "StorageUnit",
        "battery",
        bus="electricity",
        carrier="battery",
        p_nom_extendable=True,
        capital_cost=cap_cost,
        marginal_cost=0.0,
        efficiency_store=eta,
        efficiency_dispatch=eta,
        max_hours=max_hours,
        cyclic_state_of_charge=True,
    )


def _add_electrolyser(
    n: pypsa.Network,
    el_mw: float,
    el_efficiency: float,
    el_cfg: dict,
    wacc: float,
) -> None:
    cap_cost = (
        annuity_factor(wacc, el_cfg["lifetime_years"]) * el_cfg["capex_per_mw_eur"]
        + el_cfg["opex_per_mw_per_year_eur"]
    )
    n.add(
        "Link",
        "electrolyser",
        bus0="electricity",
        bus1="hydrogen",
        p_nom=el_mw,
        p_nom_extendable=False,
        efficiency=el_efficiency,
        capital_cost=cap_cost,
        marginal_cost=0.0,
    )


def _add_h2_buffer(n: pypsa.Network, buf_cfg: dict, wacc: float) -> None:
    # Store capital_cost is per MWh of e_nom (H2 LHV energy capacity)
    cap_cost = annuity_factor(wacc, buf_cfg["lifetime_years"]) * buf_cfg["capex_per_mwh_eur"]
    n.add(
        "Store",
        "h2_buffer",
        bus="hydrogen",
        carrier="H2",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=cap_cost,
        marginal_cost=0.0,
    )


def _add_dri_load(n: pypsa.Network, el_mw: float, el_efficiency: float) -> None:
    h2_demand_mw_lhv = el_mw * el_efficiency  # continuous H2 demand in MW LHV
    n.add("Load", "dri_load", bus="hydrogen", carrier="H2", p_set=h2_demand_mw_lhv)


def _add_grid_import(n: pypsa.Network, price_series: pd.Series) -> None:
    n.add(
        "Generator",
        "grid_import",
        bus="electricity",
        carrier="AC",
        p_nom=1e6,           # unconstrained import
        p_nom_extendable=False,
        marginal_cost=price_series,
        capital_cost=0.0,
    )
