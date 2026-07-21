"""Compile per-scenario summaries into a single project-level report CSV.

Invoked by Snakemake's `script:` directive (compile_report rule in viz.smk).
"""

import logging
from functools import lru_cache
from pathlib import Path

import pandas as pd
import pypsa
import yaml

from common._constants import H2_LHV_KWH_PER_KG
from common._logging import configure_logging

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)

# repo root: workflow/scripts/viz/compile_report.py -> parents[3]
_ASSUMPTIONS_PATH = Path(__file__).resolve().parents[3] / "config" / "assumptions.yaml"


@lru_cache(maxsize=1)
def _grid_volumetric_fee_eur_per_mwh() -> float:
    """The volumetric grid fee (€/MWh) added on top of the day-ahead price.

    A global assumption (config/assumptions.yaml), constant across scenarios, so
    the grid import's average priced energy can be split into the market price
    and the fee. Zero if unavailable — the split then collapses to price-only.
    """
    try:
        cfg = yaml.safe_load(_ASSUMPTIONS_PATH.read_text()) or {}
        return float(cfg.get("grid", {}).get("fee_eur_per_mwh", 0.0))
    except (OSError, ValueError, TypeError):
        return 0.0


def _h2_produced_kg(n: pypsa.Network) -> float:
    """Annual H2 produced by the electrolyser, in kg, scaled to 8760 h.

    Read from the electrolyser link's H2-side output rather than the dri_load,
    so the result reflects actual production. The two coincide when the model
    is feasible and the H2 buffer is cyclic, but only the link side stays
    correct if the load formulation changes later (e.g. flexible demand).
    """
    t_hours = len(n.snapshots)
    # PyPSA Link sign convention: p1 < 0 when the link injects power into bus1,
    # so -p1 is the (positive) H2 LHV output on the hydrogen bus.
    h2_mwh_lhv = -float(n.links_t.p1["electrolyser"].sum()) * (8760.0 / t_hours)
    return h2_mwh_lhv / (H2_LHV_KWH_PER_KG / 1000.0)


def _marginal_costs(static: pd.DataFrame, flow_t: pd.DataFrame, mc_t: pd.DataFrame) -> float:
    """Total variable cost over the simulated period for one component class.

    `flow_t` is the priced flow (generators: p; links: p0 — PyPSA applies link
    marginal costs to the input side). Hourly marginal costs (grid import)
    live in `mc_t`; everything else uses its static value.
    """
    total = 0.0
    for name in static.index:
        if name in mc_t.columns:
            mc = mc_t[name]
        else:
            mc = static.at[name, "marginal_cost"]
            if mc == 0.0:
                continue
        if name not in flow_t.columns:
            continue
        total += float((flow_t[name] * mc).sum())
    return total


PROCESS_LINKS = ("dri", "dri_ng", "eaf", "moe", "electrowinning")


def _cost_breakdown(n: pypsa.Network) -> dict[str, float]:
    """Annualised cost per system component group, €/yr; sums to the total.

    Capital costs are already per-year (annualised CAPEX × p_nom_opt); variable
    costs — grid imports, ore, EAF consumables, electrolyser variable opex —
    are scaled from the simulation period up to 8760 h so levelised costs stay
    meaningful on partial-year runs. The reported total annual cost is the sum
    of these groups, so breakdown and total cannot drift apart.
    """
    t_hours = len(n.snapshots)
    annual_scale = 8760.0 / t_hours

    def link_capital(names) -> float:
        idx = [l for l in names if l in n.links.index and n.links.at[l, "p_nom_extendable"]]
        return float((n.links.loc[idx, "capital_cost"] * n.links.loc[idx, "p_nom_opt"]).sum())

    def link_marginal(names) -> float:
        idx = [l for l in names if l in n.links.index]
        return _marginal_costs(n.links.loc[idx], n.links_t.p0, n.links_t.marginal_cost) * annual_scale

    def store_capital(name: str) -> float:
        if name not in n.stores.index or not n.stores.at[name, "e_nom_extendable"]:
            return 0.0
        return float(n.stores.at[name, "capital_cost"] * n.stores.at[name, "e_nom_opt"])

    gens = n.generators
    non_res = ("grid_import", "gas_supply")
    res_idx = gens.index[gens.p_nom_extendable & ~gens.index.isin(non_res)]
    grid_capital = 0.0
    if "grid_import" in gens.index and gens.at["grid_import", "p_nom_extendable"]:
        grid_capital = float(
            gens.at["grid_import", "capital_cost"] * gens.at["grid_import", "p_nom_opt"]
        )
    hvdc = list(n.links.index[n.links.carrier == "HVDC"])

    return {
        "res": float((gens.loc[res_idx, "capital_cost"] * gens.loc[res_idx, "p_nom_opt"]).sum()),
        "battery": float(
            (n.storage_units.capital_cost * n.storage_units.p_nom_opt)[
                n.storage_units.p_nom_extendable
            ].sum()
        ),
        # Generator marginal costs are zero except grid imports (energy price
        # + volumetric fee) and gas, which gets its own group below — so the
        # generic sum minus gas lands in the grid bucket.
        "grid": grid_capital
        + _marginal_costs(
            gens.drop(index="gas_supply", errors="ignore"),
            n.generators_t.p,
            n.generators_t.marginal_cost,
        ) * annual_scale,
        # Gas bill incl. any carbon price (both live on the gas_supply
        # generator's marginal cost).
        "gas": (
            _marginal_costs(
                gens.loc[["gas_supply"]], n.generators_t.p, n.generators_t.marginal_cost
            ) * annual_scale
            if "gas_supply" in gens.index
            else 0.0
        ),
        "electrolyser": link_capital(["electrolyser"]) + link_marginal(["electrolyser"]),
        "h2_buffer": store_capital("h2_buffer"),
        "process": link_capital(PROCESS_LINKS),
        "ore_consumables": link_marginal(PROCESS_LINKS),
        "iron_store": store_capital("iron_store"),
        "steel_store": store_capital("steel_store"),
        "transmission": link_capital(hvdc),
    }


def extract_summary(n: pypsa.Network, project_name: str, scenario_name: str) -> dict:
    """Key sizing and cost metrics as a flat dict (suitable for a one-row CSV).

    The headline levelised cost depends on the network's route: LCOH for the
    pure-H2 model (flat H2 load, no steel chain), LCOS for the steel routes
    (flat steel load). H2 production is reported whenever an electrolyser
    exists, but LCOH is only well-defined when H2 is the end product — on
    steel routes the annual cost covers the whole chain.
    """
    breakdown = _cost_breakdown(n)
    summary = {
        "project": project_name,
        "scenario": scenario_name,
        "total_annual_cost_meur": sum(breakdown.values()) / 1e6,
    }
    # Per-group annual cost columns (zero groups omitted — a scenario without a
    # component shouldn't grow the CSV). plot_lcos_bars stacks these.
    for group, value in breakdown.items():
        if value:
            summary[f"cost_{group}_meur"] = value / 1e6

    total_annual_cost = sum(breakdown.values())

    if "dri_load" in n.loads.index:
        lcoh_eur_per_kg = total_annual_cost / _h2_produced_kg(n)
        summary["lcoh_eur_per_kg"] = lcoh_eur_per_kg
        summary["lcoh_eur_per_mwh_lhv"] = lcoh_eur_per_kg * 1000.0 / H2_LHV_KWH_PER_KG

    if "steel_load" in n.loads.index:
        steel_t_per_year = float(n.loads.at["steel_load", "p_set"]) * 8760.0
        summary["lcos_eur_per_t"] = total_annual_cost / steel_t_per_year
        summary["steel_produced_mt"] = steel_t_per_year / 1e6

    if "electrolyser" in n.links.index:
        summary["h2_produced_kt"] = _h2_produced_kg(n) / 1e6

    # Levelised cost of the underlying energy carriers, €/MWh, so LCOS can be read
    # against the electricity and hydrogen that drive it.
    #   LCOE = electricity-system cost (renewables + battery + grid + transmission)
    #          per MWh of electricity generated (renewable dispatch + grid import).
    #   LCOH = (electrolyser capex/opex + H2 buffer + the electrolyser's electricity
    #          valued at LCOE) per MWh of H2 produced, LHV.
    annual_scale = 8760.0 / len(n.snapshots)
    elec_gens = [g for g in n.generators.index if g != "gas_supply"]
    elec_mwh = (float(n.generators_t.p[elec_gens].sum().sum()) * annual_scale) if elec_gens else 0.0
    elec_cost = sum(breakdown.get(k, 0.0) for k in ("res", "battery", "grid", "transmission"))
    lcoe = elec_cost / elec_mwh if elec_mwh > 0 and elec_cost > 0 else float("nan")
    if lcoe == lcoe:  # not NaN
        summary["lcoe_eur_per_mwh"] = lcoe

    # LCOE decomposition, €/MWh over the same electricity denominator so the parts
    # sum back to LCOE: renewables split by tech, grid split into connection
    # (capacity capital) vs energy (imports). Plus the average grid import price.
    if lcoe == lcoe and elec_mwh > 0:
        res_idx = n.generators.index[
            n.generators.p_nom_extendable & ~n.generators.index.isin(("grid_import", "gas_supply"))
        ]

        def _res_cost(prefix: str) -> float:
            idx = [g for g in res_idx if str(g).startswith(prefix)]
            return float((n.generators.loc[idx, "capital_cost"] * n.generators.loc[idx, "p_nom_opt"]).sum())

        def _res_mwh(prefix: str) -> float:
            idx = [g for g in res_idx if str(g).startswith(prefix)]
            return float(n.generators_t.p[idx].sum().sum()) * annual_scale if idx else 0.0

        grid_conn = grid_energy = grid_mwh = 0.0
        if "grid_import" in n.generators.index:
            gi = n.generators.loc[["grid_import"]]
            if bool(gi.at["grid_import", "p_nom_extendable"]):
                grid_conn = float(gi.at["grid_import", "capital_cost"] * gi.at["grid_import", "p_nom_opt"])
            grid_energy = _marginal_costs(gi, n.generators_t.p, n.generators_t.marginal_cost) * annual_scale
            grid_mwh = float(n.generators_t.p["grid_import"].sum()) * annual_scale

        components = {
            "lcoe_solar": _res_cost("solar"),
            "lcoe_wind_onshore": _res_cost("wind-onshore"),
            "lcoe_wind_offshore": _res_cost("wind-offshore"),
            "lcoe_storage": breakdown.get("battery", 0.0),
            "lcoe_grid_connection": grid_conn,
            "lcoe_grid_energy": grid_energy,
            "lcoe_transmission": breakdown.get("transmission", 0.0),
        }
        res_total = components["lcoe_solar"] + components["lcoe_wind_onshore"] + components["lcoe_wind_offshore"]
        if res_total > 0:
            summary["lcoe_renewables_eur_per_mwh"] = res_total / elec_mwh
        for key, cost in components.items():
            if cost > 0:
                summary[f"{key}_eur_per_mwh"] = cost / elec_mwh

        # Per-technology LCOE (€/MWh over that tech's own generation, not the
        # system total) — a fair unit cost for the renewable itself, distinct from
        # its lcoe_<tech> contribution to the system LCOE above.
        for prefix, tech in (("solar", "solar"), ("wind-onshore", "wind_onshore"),
                             ("wind-offshore", "wind_offshore")):
            cost, mwh = _res_cost(prefix), _res_mwh(prefix)
            if cost > 0 and mwh > 0:
                summary[f"lcoe_{tech}_own_eur_per_mwh"] = cost / mwh

        # Blended renewable own LCOE (over renewable generation only) — the
        # generation-weighted mean of the per-tech own LCOEs, so it sits between
        # them. Distinct from lcoe_renewables (contribution over all electricity,
        # which grid imports drag below the per-tech figures).
        res_mwh = _res_mwh("solar") + _res_mwh("wind-onshore") + _res_mwh("wind-offshore")
        if res_total > 0 and res_mwh > 0:
            summary["lcoe_renewables_own_eur_per_mwh"] = res_total / res_mwh

        # Split the average priced grid energy into the day-ahead market price and
        # the constant volumetric fee, and levelise the connection capex over the
        # same imported-MWh denominator — so all three are €/MWh *imported* and add
        # up to the all-in delivered electricity price (not mixed bases).
        if grid_mwh > 0:
            if grid_energy > 0:
                fee = _grid_volumetric_fee_eur_per_mwh()
                summary["grid_price_eur_per_mwh"] = max(grid_energy / grid_mwh - fee, 0.0)
                summary["grid_fee_eur_per_mwh"] = fee
            if grid_conn > 0:
                summary["grid_connection_eur_per_mwh_imported"] = grid_conn / grid_mwh

    if "electrolyser" in n.links.index and "steel_load" in n.loads.index:
        el_mwh = float(n.links_t.p0["electrolyser"].sum()) * annual_scale
        h2_mwh = _h2_produced_kg(n) * H2_LHV_KWH_PER_KG / 1000.0
        if h2_mwh > 0:
            elec_for_h2 = el_mwh * (lcoe if lcoe == lcoe else 0.0)
            el_cost = breakdown.get("electrolyser", 0.0)
            buf_cost = breakdown.get("h2_buffer", 0.0)
            lcoh = (el_cost + buf_cost + elec_for_h2) / h2_mwh
            summary["lcoh_eur_per_mwh_lhv"] = lcoh
            # LCOH decomposition, €/MWh LHV over the same H2 denominator so the
            # parts sum back to LCOH: the electrolyser plant, its electricity
            # (valued at LCOE) and — where built — the H2 buffer store.
            summary["lcoh_electrolyser_eur_per_mwh_lhv"] = el_cost / h2_mwh
            summary["lcoh_electricity_eur_per_mwh_lhv"] = elec_for_h2 / h2_mwh
            if buf_cost > 0:
                summary["lcoh_h2_storage_eur_per_mwh_lhv"] = buf_cost / h2_mwh

    # Steel-route process links: capacity in output units (t/h of iron or
    # steel — p_nom is input-side, so scale by the link efficiency) plus
    # utilisation, which shows how far each step actually load-follows.
    for link in PROCESS_LINKS:
        if link not in n.links.index:
            continue
        p_nom = n.links.at[link, "p_nom_opt"]
        summary[f"{link}_t_per_h_opt"] = p_nom * n.links.at[link, "efficiency"]
        summary[f"{link}_utilization"] = (
            float(n.links_t.p0[link].mean() / p_nom) if p_nom > 0 else float("nan")
        )

    if "gas_supply" in n.generators.index:
        gas_mwh = float(n.generators_t.p["gas_supply"].sum()) * (8760.0 / len(n.snapshots))
        summary["ng_gwh_lhv"] = gas_mwh / 1e3
    # How much of the iron came from the H2 shaft (production share, not capacity
    # share). Emitted for any DRI route so a pure H2-DRI reads 1.0 and a pure
    # NG-DRI 0.0 — not a missing value that downstream would coerce to 0.
    if "dri" in n.links.index or "dri_ng" in n.links.index:
        iron_h2 = -float(n.links_t.p1["dri"].sum()) if "dri" in n.links.index else 0.0
        iron_ng = -float(n.links_t.p1["dri_ng"].sum()) if "dri_ng" in n.links.index else 0.0
        total = iron_h2 + iron_ng
        summary["iron_from_h2_share"] = iron_h2 / total if total else float("nan")

    if "iron_store" in n.stores.index:
        store_t = n.stores.at["iron_store", "e_nom_opt"]
        summary["iron_store_kt"] = store_t / 1e3
        if "steel_load" in n.loads.index:
            steel_t_per_h = float(n.loads.at["steel_load", "p_set"])
            summary["iron_store_hours_steel"] = (
                store_t / steel_t_per_h if steel_t_per_h else float("nan")
            )

    if "steel_store" in n.stores.index:
        store_t = n.stores.at["steel_store", "e_nom_opt"]
        summary["steel_store_kt"] = store_t / 1e3
        if "steel_load" in n.loads.index:
            steel_t_per_h = float(n.loads.at["steel_load", "p_set"])
            summary["steel_store_hours_steel"] = (
                store_t / steel_t_per_h if steel_t_per_h else float("nan")
            )

    for gen in n.generators.index[n.generators.p_nom_extendable]:
        summary[f"{gen}_gw_opt"] = n.generators.at[gen, "p_nom_opt"] / 1e3

    if "battery" in n.storage_units.index:
        p_opt = n.storage_units.at["battery", "p_nom_opt"]
        summary["battery_gw_opt"] = p_opt / 1e3
        summary["battery_mwh_opt"] = p_opt * n.storage_units.at["battery", "max_hours"]

    if "h2_buffer" in n.stores.index:
        buffer_mwh = n.stores.at["h2_buffer", "e_nom_opt"]
        summary["h2_buffer_gwh"] = buffer_mwh / 1e3
        # Additionally as hours of average H₂ demand: the flat dri_load on the
        # pure-H2 route, or the DRI link's mean H₂ draw on the steel route.
        if "dri_load" in n.loads.index:
            h2_demand_mw = float(n.loads.at["dri_load", "p_set"])
        elif "dri" in n.links.index:
            h2_demand_mw = float(n.links_t.p0["dri"].mean())
        else:
            h2_demand_mw = 0.0
        summary["h2_buffer_hours_dri"] = (
            buffer_mwh / h2_demand_mw if h2_demand_mw else float("nan")
        )

    if "electrolyser" in n.links.index:
        el_cap = n.links.at["electrolyser", "p_nom_opt"]
        summary["electrolyser_gw"] = el_cap / 1e3
        if el_cap > 0 and "electrolyser" in n.links_t.p0.columns:
            summary["electrolyser_utilization"] = float(
                n.links_t.p0["electrolyser"].mean() / el_cap
            )
        else:
            summary["electrolyser_utilization"] = float("nan")

    if "dri_load" in n.loads.index:
        summary["dri_h2_mw_lhv"] = float(n.loads.at["dri_load", "p_set"])

    # Multi-site only (guarded so single-site reports are unchanged): one column
    # per HVDC link capacity, the total annualised transmission cost, and per-tech
    # built-capacity totals summed across candidate sites.
    hvdc = n.links.index[n.links.carrier == "HVDC"]
    if len(hvdc):
        trans_cost = 0.0
        for link in hvdc:
            cap = n.links.at[link, "p_nom_opt"]
            summary[f"{link}_gw_opt"] = cap / 1e3
            trans_cost += n.links.at[link, "capital_cost"] * cap
        summary["transmission_total_annual_cost_meur"] = trans_cost / 1e6

    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    if len(ac_buses) > 1:
        ext = n.generators[n.generators.p_nom_extendable]
        for carrier, grp in ext.groupby("carrier"):
            summary[f"{carrier}_total_gw_opt"] = grp["p_nom_opt"].sum() / 1e3

    return summary


def main() -> None:
    """Load each scenario network for the project and write the combined report CSV.

    Dedupes the netCDF inputs (collect fans out per tech row), extracts one summary
    row per scenario via `extract_summary`, rounds numerics, and writes
    results/report_<project>.csv.
    """
    project_name = snakemake.wildcards.project

    rows = []
    # networks may contain duplicates (collect fans out per tech row); dedupe
    # while preserving order so each scenario appears once.
    network_paths = list(dict.fromkeys(snakemake.input.networks))
    log.info(f"compiling report for project={project_name} ({len(network_paths)} scenarios)")
    for nc_path in network_paths:
        nc_path = Path(nc_path)
        scenario_name = nc_path.stem
        n = pypsa.Network()
        n.import_from_netcdf(nc_path)
        rows.append(extract_summary(n, project_name, scenario_name))

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df[df.select_dtypes("number").columns] = df.select_dtypes("number").round(2)
    df.to_csv(out_path, index=False)
    log.info(f"wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
