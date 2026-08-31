"""PyPSA network construction for one steelmaking route of one scenario.

Pure construction — no IO except YAML loading for assumptions, no snakemake,
no solver call. Importable from a notebook for inspection; the snakemake
entrypoint lives in `solve_network.py`.

Each network contains exactly one process route, named by the caller (the
`{route}` wildcard, from the `route` column of config/scenarios.csv — never
from the assumptions, which only carry the numbers):
  h2-only:     electricity → electrolyser → H2 → flat H2 load (LCOH model)
  h2-dri-eaf:  electricity → electrolyser → H2 → DRI shaft → iron → EAF → steel
  ng-dri-eaf:  fossil benchmark — gas supply → NG-DRI shaft → iron → EAF → steel
  mix-dri-eaf: H2 and NG DRI shafts side by side on one iron bus; the
               optimiser picks the production split (each shaft pays its own
               capex — a greenfield fuel choice, not co-firing in one shaft)
  moe-eaf:     electricity → molten oxide electrolysis → iron → EAF → steel
  ew-eaf:      electricity → iron electrowinning → iron → EAF → steel

Bus unit convention:
  electricity buses: MW AC
  hydrogen bus:      MW H2 LHV  (1 MWh H2 LHV ≈ 30 kg H2 at LHV ≈ 33.33 kWh/kg)
  gas bus:           MW CH4 LHV (supplied at a flat price incl. optional CO2 cost)
  iron bus:          t/h  (sponge iron / HBI for the DRI routes, electrolytic
                     iron plates for ew-eaf, hot metal for moe-eaf — all
                     treated as freely storable, which hot metal is not; the
                     no-store topology moe-eaf really wants is a separate fix)
  steel bus:         t/h  (liquid steel)

Process steps that consume electricity alongside their bus0 feed (DRI shaft,
EAF) are PyPSA multi-links: bus2 is an electricity bus and efficiency2 is
negative, so p2 = -efficiency2 * p0 is the electricity withdrawal.

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
from common._routes import ROUTES

HOURS_PER_YEAR = 8760.0

# Route groups used when wiring buses/components. Every steel route reaches
# the steel bus through the iron bus and the shared EAF link.
_H2_ROUTES  = ("h2-only", "h2-dri-eaf", "mix-dri-eaf")   # electrolyser + H2 buffer
_GAS_ROUTES = ("ng-dri-eaf", "mix-dri-eaf")              # gas bus + NG-DRI shaft
_H2_DRI_ROUTES = ("h2-dri-eaf", "mix-dri-eaf")           # H2 DRI shaft
_IRON_ROUTES = ("h2-dri-eaf", "ng-dri-eaf", "mix-dri-eaf", "moe-eaf", "ew-eaf")


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
    route: str,
    assumptions: dict,
    cf_timeseries: pd.DataFrame,
    price_series: pd.Series | None = None,
    sites: pd.DataFrame | None = None,
    demand_site: str | None = None,
) -> pypsa.Network:
    """Build (but do not solve) the PyPSA network for one scenario.

    `route` selects the process chain and comes from the `{route}` wildcard;
    `assumptions` is the merged base+overlay dict and carries only numbers.
    `price_series` is an optional hourly
    €/MWh grid price on the same index, which adds a grid-import generator
    (with connection charges from `assumptions["grid"]`) when present.

    Single-site mode (`sites is None`): `cf_timeseries` has one column per RES tech
    (names matching keys in `assumptions.res`); every generator sits on a single
    `electricity` bus.

    Multi-site mode (`sites` given): `cf_timeseries.columns` is a MultiIndex
    (site_id, tech); each site gets its own `electricity_{site_id}` bus carrying its
    generators, and an extendable HVDC link connects every site to the demand site
    (`demand_site`), which hosts the process chain, storage and any grid import.
    `sites` is indexed by site_id with columns `x` (lon), `y` (lat).
    """
    if route not in ROUTES:
        raise ValueError(f"unknown route '{route}' — expected one of {ROUTES}")
    wacc = assumptions["finance"]["default_wacc"]
    plant = assumptions["plant"]

    multisite = sites is not None

    n = pypsa.Network()
    n.set_snapshots(cf_timeseries.index)

    if multisite:
        res_techs = list(dict.fromkeys(tech for _, tech in cf_timeseries.columns))
        elec_bus = f"electricity_{demand_site}"
    else:
        res_techs = list(cf_timeseries.columns)
        elec_bus = "electricity"

    _add_carriers(n, res_techs=res_techs, route=route, multisite=multisite)
    if multisite:
        _add_buses_multisite(n, sites, demand_site, route)
    else:
        _add_buses(n, route)
    _add_generators(n, cf_timeseries, assumptions["res"], wacc, multisite=multisite)
    _add_battery(n, assumptions["battery"], wacc, bus=elec_bus)

    if route in _H2_ROUTES:
        el_cfg = assumptions["electrolyser"]
        el_efficiency = H2_LHV_KWH_PER_KG / el_cfg["efficiency_kwh_per_kg"]
        if route == "h2-only":
            # Floor sized to meet the flat H2 demand at availability_target.
            # Steel routes leave electrolyser sizing entirely to the optimiser —
            # the flat steel load forces enough capacity through the chain anyway.
            el_mw = dri_to_el_mw(
                dri_mt_per_year=plant["dri_mt_per_year"],
                h2_intensity_kg_per_t_dri=plant["h2_intensity_kg_per_t_dri"],
                efficiency_kwh_per_kg=el_cfg["efficiency_kwh_per_kg"],
                availability_target=plant["availability_target"],
            )
        else:
            el_mw = 0.0
        _add_electrolyser(n, el_mw, el_efficiency, el_cfg, wacc, bus0=elec_bus)
        _add_h2_buffer(n, assumptions["h2_buffer"], wacc)

    if route == "h2-only":
        _add_dri_load(n, el_mw, el_efficiency, plant["availability_target"])
    else:
        steel_t_per_h = plant["steel_mt_per_year"] * 1e6 / HOURS_PER_YEAR
        _add_steel_load(n, steel_t_per_h)
        _add_steel_store(n, assumptions["steel_store"], wacc, steel_t_per_h)
        if route in _H2_DRI_ROUTES:
            _add_dri_link(n, plant, assumptions["dri-h2"], wacc, elec_bus)
        if route in _GAS_ROUTES:
            _add_gas_supply(n, assumptions["natural_gas"])
            _add_ng_dri_link(n, assumptions["dri-ng"], wacc, elec_bus)
        if route == "moe-eaf":
            _add_moe_link(n, assumptions["moe"], wacc, elec_bus)
        if route == "ew-eaf":
            _add_ew_link(n, assumptions["ew"], wacc, elec_bus)
        _add_iron_store(n, assumptions["iron_store"], wacc)
        _add_eaf_link(n, assumptions["eaf"], wacc, elec_bus)

    if multisite:
        _add_transmission(n, sites, demand_site, assumptions["transmission"], wacc)

    if price_series is not None:
        _add_grid_import(n, price_series, assumptions["grid"], wacc, bus=elec_bus)

    return n


def _add_carriers(
    n: pypsa.Network, res_techs: list[str], route: str, multisite: bool = False
) -> None:
    """Register every carrier referenced by a component, before those components are added.

    Otherwise PyPSA's consistency check warns ("carriers which are not defined")
    and leaves n.carriers empty, which blocks carrier-aware features (CO2
    constraints, grouped stats). Only the carriers the chosen route actually
    uses are added; the HVDC carrier is only added in multi-site mode.
    """
    base = ["AC", "battery"]
    if route in _H2_ROUTES:
        base += ["H2", "electrolyser"]
    if route in _GAS_ROUTES:
        base += ["gas", "dri-ng"]
    if route in _IRON_ROUTES:
        base += ["iron", "eaf"]
    if route != "h2-only":
        base += ["steel"]
    if route in _H2_DRI_ROUTES:
        base += ["dri-h2"]
    if route == "ew-eaf":
        base += ["ew"]
    if route == "moe-eaf":
        base += ["moe"]
    if multisite:
        base.append("HVDC")
    carriers = list(dict.fromkeys([*base, *res_techs]))
    n.add("Carrier", carriers)


def _route_process_buses(route: str) -> list[tuple[str, str]]:
    """(bus name, carrier) pairs for the process buses the route needs."""
    buses = []
    if route in _H2_ROUTES:
        buses.append(("hydrogen", "H2"))
    if route in _GAS_ROUTES:
        buses.append(("gas", "gas"))
    if route in _IRON_ROUTES:
        buses.append(("iron", "iron"))
    if route != "h2-only":
        buses.append(("steel", "steel"))
    return buses


def _add_buses(n: pypsa.Network, route: str) -> None:
    """Add the electricity (AC) bus plus the route's process buses."""
    n.add("Bus", "electricity", carrier="AC")
    for name, carrier in _route_process_buses(route):
        n.add("Bus", name, carrier=carrier)


def _add_buses_multisite(
    n: pypsa.Network, sites: pd.DataFrame, demand_site: str, route: str
) -> None:
    """Add one AC bus per site (with lon/lat coords) plus the demand-site process buses."""
    for site_id, row in sites.iterrows():
        n.add("Bus", f"electricity_{site_id}", carrier="AC", x=row["x"], y=row["y"])
    dem = sites.loc[demand_site]
    for name, carrier in _route_process_buses(route):
        n.add("Bus", name, carrier=carrier, x=dem["x"], y=dem["y"])


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
        # Floor sized by dri_to_el_mw to meet annual demand at availability_target
        # (h2_only route; steel routes pass 0). Optimiser may grow beyond this when
        # the headroom pays for itself via buffer-mediated arbitrage of cheap RES.
        p_nom_min=el_mw,
        efficiency=el_efficiency,
        capital_cost=cap_cost,
        # Variable opex (water, consumables) per MWh electricity input (p0 side).
        marginal_cost=el_cfg.get("varopex_eur_per_mwh_el", 0.0),
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


def _add_steel_load(n: pypsa.Network, steel_t_per_h: float) -> None:
    """Add the flat steel demand (t/h) on the steel bus."""
    n.add("Load", "steel_load", bus="steel", carrier="steel", p_set=steel_t_per_h)


def _process_capital_cost(cfg: dict, wacc: float, output_t_per_p0_unit: float) -> float:
    """Convert a €/(t/yr of output capacity) capex+opex quote into PyPSA capital_cost.

    Link p_nom is in bus0 units (MW or t/h of input), while process capex is
    quoted per tonne of annual *output* capacity, so scale by output per unit
    of p0 and hours per year. Quotes are taken as nameplate (8760 h) annual
    capacity — at eyeball fidelity the distinction from availability-derated
    quotes is noise.
    """
    per_t_per_year = (
        annuity_factor(wacc, cfg["lifetime_years"]) * cfg["capex_per_t_per_year_eur"]
        + cfg["opex_per_t_per_year_eur"]
    )
    return per_t_per_year * output_t_per_p0_unit * HOURS_PER_YEAR


def _add_dri_link(
    n: pypsa.Network, plant: dict, dri_cfg: dict, wacc: float, elec_bus: str
) -> None:
    """DRI shaft: hydrogen (MW LHV) → sponge iron (t/h), drawing shaft electricity.

    `efficiency` converts MWh H2 into t iron via the plant's H2 intensity;
    `efficiency2` (negative) withdraws shaft electricity in fixed proportion.
    Ore cost enters as marginal_cost, converted from €/t iron to per-MWh-H2 (p0).
    """
    t_iron_per_mwh_h2 = 1000.0 / (
        plant["h2_intensity_kg_per_t_dri"] * H2_LHV_KWH_PER_KG
    )
    n.add(
        "Link",
        "dri-h2",
        bus0="hydrogen",
        bus1="iron",
        bus2=elec_bus,
        carrier="dri-h2",
        p_nom_extendable=True,
        efficiency=t_iron_per_mwh_h2,
        efficiency2=-dri_cfg["el_mwh_per_t"] * t_iron_per_mwh_h2,
        p_min_pu=dri_cfg["p_min_pu"],
        capital_cost=_process_capital_cost(dri_cfg, wacc, t_iron_per_mwh_h2),
        marginal_cost=dri_cfg["ore_eur_per_t"] * t_iron_per_mwh_h2,
    )


def _add_gas_supply(n: pypsa.Network, ng_cfg: dict) -> None:
    """Unlimited natural-gas supply (MW CH4 LHV) at a flat price on the gas bus.

    An extendable zero-capex generator, so the optimiser draws freely; the
    marginal cost is the gas price plus any carbon price on combustion CO2.
    Kept as a generator (not a link marginal cost) so the gas bill lands in
    its own cost group in reports.
    """
    marginal = ng_cfg["price_eur_per_mwh"] + (
        ng_cfg["co2_t_per_mwh"] * ng_cfg["co2_price_eur_per_t"]
    )
    n.add(
        "Generator",
        "gas_supply",
        bus="gas",
        carrier="gas",
        p_nom_extendable=True,
        capital_cost=0.0,
        marginal_cost=marginal,
    )


def _add_ng_dri_link(
    n: pypsa.Network, dri_ng_cfg: dict, wacc: float, elec_bus: str
) -> None:
    """NG-DRI shaft: natural gas (MW LHV) → sponge iron (t/h), drawing shaft electricity.

    Mirrors the H2 shaft (`_add_dri_link`) with gas as the reductant/fuel; the
    gas itself is priced on the gas_supply generator, so only ore enters here
    as marginal cost.
    """
    t_iron_per_mwh_gas = 1.0 / dri_ng_cfg["gas_mwh_per_t"]
    n.add(
        "Link",
        "dri-ng",
        bus0="gas",
        bus1="iron",
        bus2=elec_bus,
        carrier="dri-ng",
        p_nom_extendable=True,
        efficiency=t_iron_per_mwh_gas,
        efficiency2=-dri_ng_cfg["el_mwh_per_t"] * t_iron_per_mwh_gas,
        p_min_pu=dri_ng_cfg["p_min_pu"],
        capital_cost=_process_capital_cost(dri_ng_cfg, wacc, t_iron_per_mwh_gas),
        marginal_cost=dri_ng_cfg["ore_eur_per_t"] * t_iron_per_mwh_gas,
    )


def _add_eaf_link(n: pypsa.Network, eaf_cfg: dict, wacc: float, elec_bus: str) -> None:
    """EAF: iron (t/h) → steel (t/h), drawing melting electricity.

    Per-t-steel quotes (electricity, consumables, capex) are scaled by the
    iron→steel yield onto the link's p0 side (t/h iron).
    """
    t_steel_per_t_iron = 1.0 / eaf_cfg["iron_t_per_t_steel"]
    n.add(
        "Link",
        "eaf",
        bus0="iron",
        bus1="steel",
        bus2=elec_bus,
        carrier="eaf",
        p_nom_extendable=True,
        efficiency=t_steel_per_t_iron,
        efficiency2=-eaf_cfg["el_mwh_per_t"] * t_steel_per_t_iron,
        p_min_pu=eaf_cfg["p_min_pu"],
        capital_cost=_process_capital_cost(eaf_cfg, wacc, t_steel_per_t_iron),
        marginal_cost=eaf_cfg["consumables_eur_per_t"] * t_steel_per_t_iron,
    )


def _add_moe_link(n: pypsa.Network, moe_cfg: dict, wacc: float, elec_bus: str) -> None:
    """Molten oxide electrolysis: electricity (MW) → liquid iron (t/h).

    Ladle metallurgy is the shared EAF link, priced down by a scenario overlay
    where a scenario wants it, rather than a component of its own.
    """
    t_iron_per_mwh = 1.0 / moe_cfg["el_mwh_per_t"]
    n.add(
        "Link",
        "moe",
        bus0=elec_bus,
        bus1="iron",
        carrier="moe",
        p_nom_extendable=True,
        efficiency=t_iron_per_mwh,
        p_min_pu=moe_cfg["p_min_pu"],
        capital_cost=_process_capital_cost(moe_cfg, wacc, t_iron_per_mwh),
        marginal_cost=moe_cfg["ore_eur_per_t"] * t_iron_per_mwh,
    )


def _add_ew_link(n: pypsa.Network, ew_cfg: dict, wacc: float, elec_bus: str) -> None:
    """Iron electrowinning: electricity (MW) → electrolytic iron plates (t/h)."""
    t_iron_per_mwh = 1.0 / ew_cfg["el_mwh_per_t"]
    n.add(
        "Link",
        "ew",
        bus0=elec_bus,
        bus1="iron",
        carrier="ew",
        p_nom_extendable=True,
        efficiency=t_iron_per_mwh,
        p_min_pu=ew_cfg["p_min_pu"],
        capital_cost=_process_capital_cost(ew_cfg, wacc, t_iron_per_mwh),
        marginal_cost=ew_cfg["ore_eur_per_t"] * t_iron_per_mwh,
    )


def _add_iron_store(n: pypsa.Network, store_cfg: dict, wacc: float) -> None:
    """Add the extendable, cyclic iron stockpile (t) on the iron bus.

    Deliberately cheap-but-not-free (see assumptions) so the optimal stockpile
    size is unique and meaningful in reports.
    """
    cap_cost = annuity_factor(wacc, store_cfg["lifetime_years"]) * store_cfg["capex_per_t_eur"]
    n.add(
        "Store",
        "iron_store",
        bus="iron",
        carrier="iron",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=cap_cost,
        marginal_cost=0.0,
    )


def _add_steel_store(
    n: pypsa.Network, store_cfg: dict, wacc: float, steel_t_per_h: float
) -> None:
    """Add the extendable, cyclic finished-steel inventory buffer (t) on the steel bus.

    The supply-side representation of periodic demand: production flexes to fill
    inventory while the flat steel load draws it down, so the plant decouples
    hourly production from steady delivery by carrying stock — demand itself
    stays fixed. `max_weeks` caps the buffer at a realistic inventory horizon so
    the flexibility can't degenerate into whole-year arbitrage; the annuitised
    per-tonne capex keeps the built size unique and reportable.
    """
    cap_cost = annuity_factor(wacc, store_cfg["lifetime_years"]) * store_cfg["capex_per_t_eur"]
    e_nom_max = steel_t_per_h * store_cfg["max_weeks"] * 7 * 24
    n.add(
        "Store",
        "steel_store",
        bus="steel",
        carrier="steel",
        e_nom_extendable=True,
        e_nom_max=e_nom_max,
        e_cyclic=True,
        capital_cost=cap_cost,
        marginal_cost=0.0,
    )


def _add_grid_import(
    n: pypsa.Network,
    price_series: pd.Series,
    grid_cfg: dict,
    wacc: float,
    bus: str = "electricity",
) -> None:
    """Add an extendable grid-import generator with connection charges.

    The optimiser sizes the connection: capacity pays annuitised connection
    capex plus the yearly per-MW fee; imported energy pays the hourly day-ahead
    price plus the volumetric fee. Import-only — selling surplus to the grid is
    deliberately out of scope for now.
    """
    cap_cost = (
        annuity_factor(wacc, grid_cfg["connection_lifetime_years"])
        * grid_cfg["connection_capex_eur_per_mw"]
        + grid_cfg["fee_eur_per_mw_per_year"]
    )
    n.add(
        "Generator",
        "grid_import",
        bus=bus,
        carrier="AC",
        p_nom_extendable=True,
        capital_cost=cap_cost,
        marginal_cost=price_series + grid_cfg["fee_eur_per_mwh"],
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
