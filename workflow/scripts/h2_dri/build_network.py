"""PyPSA network construction for a DRI-hydrogen project scenario.

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

from _helpers import annuity_factor, dri_to_el_mw, haversine_km

from common._constants import H2_LHV_KWH_PER_KG


def _deep_merge(base: dict, overlay: dict) -> dict:
    """Recursively merge `overlay` into `base` (neither input mutated).

    Overlay leaves replace base leaves; dict branches are merged key-by-key.
    """
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
    sites: pd.DataFrame | None = None,
    demand_site: str | None = None,
) -> pypsa.Network:
    """Build (but do not solve) the PyPSA network for one scenario.

    `assumptions` is the merged base+overlay dict; `price_series` is an optional
    hourly €/MWh grid price on the same index, which adds a grid-import generator
    when present.

    Single-site mode (`sites is None`): `cf_timeseries` has one column per RES tech
    (names matching keys in `assumptions.res`); every generator sits on a single
    `electricity` bus. This is the original, unchanged behaviour.

    Multi-site mode (`sites` given): `cf_timeseries.columns` is a MultiIndex
    (site_id, tech); each site gets its own `electricity_{site_id}` bus carrying its
    generators, and an extendable HVDC link connects every site to the demand site
    (`demand_site`), which hosts the electrolyser, battery, H2 buffer, DRI load and
    any grid import. `sites` is indexed by site_id with columns `x` (lon), `y` (lat).
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

    multisite = sites is not None

    n = pypsa.Network()
    n.set_snapshots(cf_timeseries.index)

    if multisite:
        res_techs = list(dict.fromkeys(tech for _, tech in cf_timeseries.columns))
        elec_bus = f"electricity_{demand_site}"
    else:
        res_techs = list(cf_timeseries.columns)
        elec_bus = "electricity"

    _add_carriers(n, res_techs=res_techs, multisite=multisite)
    if multisite:
        _add_buses_multisite(n, sites, demand_site)
    else:
        _add_buses(n)
    _add_generators(n, cf_timeseries, assumptions["res"], wacc, multisite=multisite)
    _add_battery(n, assumptions["battery"], wacc, bus=elec_bus)
    _add_electrolyser(n, el_mw, el_efficiency, el_cfg, wacc, bus0=elec_bus)
    _add_h2_buffer(n, assumptions["h2_buffer"], wacc)
    _add_dri_load(n, el_mw, el_efficiency, plant["availability_target"])

    if multisite:
        _add_transmission(n, sites, demand_site, assumptions["transmission"], wacc)

    if price_series is not None:
        _add_grid_import(n, price_series, bus=elec_bus)

    return n


def _add_carriers(n: pypsa.Network, res_techs: list[str], multisite: bool = False) -> None:
    """Register every carrier referenced by a component, before those components are added.

    Otherwise PyPSA's consistency check warns ("carriers which are not defined")
    and leaves n.carriers empty, which blocks carrier-aware features (CO2
    constraints, grouped stats). The HVDC carrier is only added in multi-site mode.
    """
    base = ["AC", "H2", "battery", "electrolyser"]
    if multisite:
        base.append("HVDC")
    carriers = list(dict.fromkeys([*base, *res_techs]))
    n.add("Carrier", carriers)


def _add_buses(n: pypsa.Network) -> None:
    """Add the electricity (AC) and hydrogen (H2 LHV) buses."""
    n.add("Bus", "electricity", carrier="AC")
    n.add("Bus", "hydrogen", carrier="H2")


def _add_buses_multisite(
    n: pypsa.Network, sites: pd.DataFrame, demand_site: str
) -> None:
    """Add one AC bus per site (with lon/lat coords) plus the demand-site H2 bus."""
    for site_id, row in sites.iterrows():
        n.add("Bus", f"electricity_{site_id}", carrier="AC", x=row["x"], y=row["y"])
    dem = sites.loc[demand_site]
    n.add("Bus", "hydrogen", carrier="H2", x=dem["x"], y=dem["y"])


def _add_generators(
    n: pypsa.Network,
    cf_timeseries: pd.DataFrame,
    res_cfg: dict,
    wacc: float,
    multisite: bool = False,
) -> None:
    """Add one extendable RES generator per CF column, costed from assumptions.

    Single-site: columns are bare tech keys; each generator is named for its tech
    and sits on the `electricity` bus. Multi-site: columns are a (site, tech)
    MultiIndex; each generator is named `{tech}_{site}` and sits on its site's bus.
    Either way the tech key is looked up in `res_cfg`; orientation-suffixed keys
    (e.g. `solar_az180`, `..._east_30`) fall back to their base tech. Capital cost
    is annuitised CAPEX + fixed OPEX; the CF profile enters as p_max_pu.
    """
    for col in cf_timeseries.columns:
        if multisite:
            site, tech = col
            name = site  # site id already encodes tech+cell, e.g. "solar-c00"
            bus = f"electricity_{site}"
        else:
            tech = col
            name = tech
            bus = "electricity"

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
            name,
            bus=bus,
            carrier=tech,
            p_nom_extendable=True,
            capital_cost=cap_cost,
            marginal_cost=0.0,
            p_max_pu=cf_timeseries[col],
        )


def _add_battery(
    n: pypsa.Network, bat_cfg: dict, wacc: float, bus: str = "electricity"
) -> None:
    """Add an extendable battery; energy CAPEX is folded into the per-MW cost at fixed duration."""
    eta = bat_cfg["efficiency_roundtrip"] ** 0.5
    max_hours = bat_cfg["max_hours"]
    # Fold energy capex into per-MW capital cost (assumes fixed duration = max_hours)
    cap_cost = annuity_factor(wacc, bat_cfg["lifetime_years"]) * (
        bat_cfg["capex_per_mw_eur"] + bat_cfg["capex_per_mwh_eur"] * max_hours
    )
    n.add(
        "StorageUnit",
        "battery",
        bus=bus,
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
    bus0: str = "electricity",
) -> None:
    """Add the extendable electrolyser link (electricity → hydrogen), floored at `el_mw`."""
    cap_cost = (
        annuity_factor(wacc, el_cfg["lifetime_years"]) * el_cfg["capex_per_mw_eur"]
        + el_cfg["opex_per_mw_per_year_eur"]
    )
    n.add(
        "Link",
        "electrolyser",
        bus0=bus0,
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
    """Add the extendable, cyclic H2 storage buffer (a Store on the hydrogen bus)."""
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
    """Add the constant DRI hydrogen demand (annual-average MW H2 LHV) on the hydrogen bus."""
    # The plant's hydrogen demand is the annual-average value, not the
    # electrolyser's rated max output. dri_to_el_mw sizes el_mw to deliver the
    # annual quota at the chosen availability, so el_mw * el_efficiency is the
    # rated output and el_mw * el_efficiency * availability_target collapses
    # back to the constant demand implied by dri_mt_per_year * h2_intensity.
    h2_demand_mw_lhv = el_mw * el_efficiency * availability_target
    n.add("Load", "dri_load", bus="hydrogen", carrier="H2", p_set=h2_demand_mw_lhv)


def _add_grid_import(
    n: pypsa.Network, price_series: pd.Series, bus: str = "electricity"
) -> None:
    """Add an unconstrained grid-import generator priced at the hourly `price_series`."""
    n.add(
        "Generator",
        "grid_import",
        bus=bus,
        carrier="AC",
        p_nom=1e6,           # unconstrained import
        p_nom_extendable=False,
        marginal_cost=price_series,
        capital_cost=0.0,
    )


def _add_transmission(
    n: pypsa.Network,
    sites: pd.DataFrame,
    demand_site: str,
    trans_cfg: dict,
    wacc: float,
) -> None:
    """Add one extendable bidirectional HVDC link from each RES site to the demand site.

    Capital cost and loss both scale with routed distance: capex per MW is the
    annuitised €/MW/km times the great-circle distance inflated by an indirect-route
    factor; efficiency is 1 − losses_pct_per_1000km × km/1000. The optimiser sizes
    each link (zero where the site is uneconomic), so siting emerges from the LP.
    """
    dem = sites.loc[demand_site]
    cap_per_km = (
        annuity_factor(wacc, trans_cfg["lifetime_years"])
        * trans_cfg["cost_per_mw_per_km_eur"]
    )
    loss_per_km = trans_cfg["losses_pct_per_1000km"] / 100.0 / 1000.0
    route = trans_cfg["indirect_route_factor"]

    for site_id, row in sites.iterrows():
        if site_id == demand_site:
            continue
        length_km = haversine_km(dem["x"], dem["y"], row["x"], row["y"]) * route
        n.add(
            "Link",
            f"hvdc_{site_id}",
            bus0=f"electricity_{site_id}",
            bus1=f"electricity_{demand_site}",
            carrier="HVDC",
            p_nom_extendable=True,
            # Unidirectional (default p_min_pu=0): a candidate site is a pure source
            # (generator only, no load/storage), so power only ever flows site→plant
            # and reverse capability would never be used. Verified to give an
            # identical optimum to a bidirectional (p_min_pu=-1) link. NOTE: if
            # per-site batteries/loads are ever added (so a site can import to
            # charge), this must become bidirectional again.
            efficiency=max(0.0, 1.0 - loss_per_km * length_km),
            capital_cost=cap_per_km * length_km,
            marginal_cost=0.0,
        )
