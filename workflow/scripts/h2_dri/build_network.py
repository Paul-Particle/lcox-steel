"""
PyPSA network construction for a DRI-hydrogen project scenario.

Pure construction — no IO except YAML loading for assumptions, no snakemake,
no solver call. Importable from a notebook for inspection; the snakemake
entrypoint lives in `solve_network.py`.

Bus unit convention: MW throughout.
  electricity bus: MW AC
  hydrogen bus:    MW H2 LHV  (1 MWh H2 LHV ≈ 30 kg H2 at LHV ≈ 33.33 kWh/kg)

Electrolyser efficiency is:
  efficiency = h2_lhv_kwh_per_kg / efficiency_kwh_per_kg
             = 33.33 / 55 ≈ 0.606 (MW H2 LHV per MW electricity)
"""

import re
from pathlib import Path

import pandas as pd
import pypsa
import yaml

from _helpers import annuity_factor, dri_to_el_mw

from common._constants import H2_LHV_KWH_PER_KG


def _deep_merge(base: dict, overlay: dict) -> dict:
    """Recursively merge `overlay` into `base`. Overlay leaves replace base
    leaves; dict branches are merged key-by-key. Neither input is mutated."""
    out = dict(base)
    for k, v in overlay.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_assumptions(base_path: Path, overlay_path: Path | None) -> dict:
    """Load base assumptions and optionally merge a per-scenario overlay on top.

    `overlay_path == base_path` (or None) means no overlay — caller passes the
    base path twice in the snakemake input when no variant applies."""
    base = yaml.safe_load(Path(base_path).read_text()) or {}
    if overlay_path is None or Path(overlay_path) == Path(base_path):
        return base
    overlay = yaml.safe_load(Path(overlay_path).read_text()) or {}
    return _deep_merge(base, overlay)


def build_network(
    assumptions: dict,
    cf_timeseries: pd.DataFrame,
    price_series: pd.Series | None = None,
) -> pypsa.Network:
    """
    Build (but do not solve) a PyPSA network.

    assumptions: full assumptions dict (from snakemake.config)
    cf_timeseries: DataFrame indexed by DatetimeIndex, columns = tech names in assumptions.res
    price_series: optional hourly price Series (€/MWh), same index as cf_timeseries
    """
    wacc = assumptions["finance"]["default_wacc"]
    el_cfg = assumptions["electrolyser"]
    plant = assumptions["plant"]
    el_mw = dri_to_el_mw(
        dri_mt_per_year=plant["dri_mt_per_year"],
        h2_intensity_kg_per_t_dri=plant["h2_intensity_kg_per_t_dri"],
        efficiency_kwh_per_kg=el_cfg["efficiency_kwh_per_kg"],
        availability_target=plant["availability_target"],
    )
    el_efficiency = H2_LHV_KWH_PER_KG / el_cfg["efficiency_kwh_per_kg"]

    n = pypsa.Network()
    n.set_snapshots(cf_timeseries.index)

    _add_carriers(n, res_techs=list(cf_timeseries.columns))
    _add_buses(n)
    _add_generators(n, cf_timeseries, assumptions["res"], wacc)
    _add_battery(n, assumptions["battery"], wacc)
    _add_electrolyser(n, el_mw, el_efficiency, el_cfg, wacc)
    _add_h2_buffer(n, assumptions["h2_buffer"], wacc)
    _add_dri_load(n, el_mw, el_efficiency, plant["availability_target"])

    if price_series is not None:
        _add_grid_import(n, price_series)

    return n


def _add_carriers(n: pypsa.Network, res_techs: list[str]) -> None:
    # Register every carrier string referenced by a component before that
    # component is added — otherwise PyPSA's consistency check warns
    # ("carriers which are not defined") and n.carriers stays empty,
    # which blocks carrier-aware features (CO2 constraints, grouped stats).
    carriers = list(dict.fromkeys(["AC", "H2", "battery", "electrolyser", *res_techs]))
    n.add("Carrier", carriers)


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
        cfg = res_cfg.get(tech)
        if cfg is None:
            base = re.sub(r"_(east|west)_\d+$|_az\d+$", "", tech)
            cfg = res_cfg.get(base)
        if cfg is None:
            raise KeyError(f"No assumptions found for tech '{tech}' — add it to assumptions.yaml")

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
        carrier="electrolyser",
        p_nom_extendable=True,
        # Floor sized by dri_to_el_mw to meet annual demand at availability_target.
        # Optimiser may grow beyond this when the headroom pays for itself via
        # buffer-mediated arbitrage of cheap RES hours.
        p_nom_min=el_mw,
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


def _add_dri_load(
    n: pypsa.Network, el_mw: float, el_efficiency: float, availability_target: float
) -> None:
    # The plant's hydrogen demand is the annual-average value, not the
    # electrolyser's rated max output. dri_to_el_mw sizes el_mw to deliver the
    # annual quota at the chosen availability, so el_mw * el_efficiency is the
    # rated output and el_mw * el_efficiency * availability_target collapses
    # back to the constant demand implied by dri_mt_per_year * h2_intensity.
    h2_demand_mw_lhv = el_mw * el_efficiency * availability_target
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
