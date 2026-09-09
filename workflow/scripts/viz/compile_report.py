"""Compile per-run summaries into the scenario-level report CSV.

Invoked by Snakemake's `script:` directive (compile_report rule in viz.smk).
Writes two files: the report viz reads, and a hidden `_diag.csv` beside it that
keeps the zones the report did not select. Each run carries the `inputs_hash` the
solve stamped into the network, so a number stays tied to its inputs. What the
files look like is common/_report_schema.py's business.
"""

import logging
import re
from pathlib import Path

import pandas as pd
import pypsa
import yaml

from common._constants import H2_LHV_KWH_PER_KG
from common._logging import configure_logging
from common._report_schema import (
    ELECTRICITY_USERS,
    EMISSION_STEPS,
    PROCESS_LINKS,
    RES_TECHS,
    field_stem,
    write_report_file,
)
from common._runs import load_scenarios, zone_parents

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)

def _carrier_key(carrier: str) -> str:
    """The emission table's key for a generator carrier: `wind-onshore` → `wind_onshore`.

    Orientation-suffixed keys fall back to their base tech, the same way
    `build_network._add_generators` resolves them against the res assumptions.
    """
    return field_stem(re.sub(r"_(east|west)_\d+$|_az\d+$", "", carrier))


def _grid_intensity(
    emissions_cfg: dict, area: str, grid_mix: pd.DataFrame | None,
    factors: dict[str, float], snapshots: pd.Index,
) -> tuple[pd.Series, str]:
    """t CO2e per MWh imported, hour by hour, and where the figure came from.

    A grid series solved on `variant: full` carries the area's generation by
    carrier, so the intensity is that mix weighted by the factor table — which
    is the whole point of keying the table to the shared carrier vocabulary. A
    series solved on `dayahead` carries prices only and has no mix in it at all,
    so the area's annual figure stands in and the report says so.

    Production-based either way: a zone that imports coal power from next door
    reads as clean as its own plants. The cross-border columns are in the `full`
    series if that is ever worth fixing.
    """
    if grid_mix is not None:
        carriers = [col for col in grid_mix.columns if col in factors]
        if carriers:
            # ENTSO-E publishes the odd negative generation hour; a negative
            # would otherwise pull the weighted mean the wrong way.
            generation = grid_mix[carriers].clip(lower=0.0)
            generated = generation.sum(axis=1)
            carried = sum(generation[carrier] * factors[carrier] for carrier in carriers)
            hourly = (carried / generated.where(generated > 0)).reindex(snapshots)
            return hourly.ffill().bfill(), "mix"
    return pd.Series(emissions_cfg["grid_t_co2e_per_mwh"][area], index=snapshots), "area_default"


def _electricity_draw(n: pypsa.Network, link: str) -> tuple[pd.Series, str]:
    """A link's electricity draw in MW, and the bus it draws on.

    Some links take electricity as their input (bus0: an electrolyser, a MOE
    cell), others as a by-draw alongside the conversion they exist for (bus2: a
    shaft's auxiliaries, the furnace's melt). Either way PyPSA reads a positive
    port flow as withdrawn from that bus, so both are positive as they stand —
    it is the bus that differs, not the sign.
    """
    if n.buses.at[n.links.at[link, "bus0"], "carrier"] == "AC":
        return n.links_t.p0[link], n.links.at[link, "bus0"]
    return n.links_t.p2[link], n.links.at[link, "bus2"]


def _emissions_breakdown(
    n: pypsa.Network, emissions_cfg: dict, natural_gas_cfg: dict,
    area: str, transport_legs: dict | None, grid_mix: pd.DataFrame | None,
) -> dict[str, object]:
    """What the run emitted in a year, cut four ways.

    `by_step` is annual t CO2e per step and stacks to the total. `electricity_mwh`
    and `electricity_t` are each user's draw and the emissions of that draw alone
    — kept apart from `by_step` because a gas-fired shaft's step total carries its
    combustion too, and dividing that by its MWh would report a furnace as
    buying impossibly dirty power. `sources` splits the total into electricity,
    gas and freight; `grid_source` says where the grid intensity came from.

    Accounting only: none of this reaches the objective, so nothing here can
    move a solve. And it covers the run's energy and freight alone — the process
    steps' own direct emissions are outside the model boundary, which is why the
    result is not a CBAM or an ETS figure. See the `emissions` block in
    config/assumptions.yaml.

    Electricity is attributed hour by hour at the intensity of the system that
    supplied it, so a user that ran when the wind blew carries less than one
    that ran flat out. Two systems can exist: the producing area's, where the
    plant's own renewables mix with whatever it imports, and — on an export
    route — the destination's, which is a grid and nothing else.
    """
    basis = emissions_cfg["basis"]
    factors = {carrier: values[basis]
               for carrier, values in emissions_cfg["electricity_t_co2e_per_mwh"].items()}
    annual = 8760.0 / len(n.snapshots)

    # The destination is found by who feeds it rather than by name, so the bus
    # can be renamed in build_network without breaking the report.
    destination_bus = (n.generators.at["destination_supply", "bus"]
                       if "destination_supply" in n.generators.index else None)
    destination_intensity = pd.Series(
        emissions_cfg["destination_t_co2e_per_mwh"], index=n.snapshots
    )

    grid_source = "none"
    if "grid_import" in n.generators.index:
        grid_hourly, grid_source = _grid_intensity(
            emissions_cfg, area, grid_mix, factors, n.snapshots
        )
    else:
        grid_hourly = pd.Series(0.0, index=n.snapshots)

    # What a MWh on the producing area's electricity buses carried, hour by hour:
    # its own generation at the carriers' factors, its imports at the grid's.
    home_buses = set(n.buses.index[n.buses.carrier == "AC"]) - {destination_bus}
    home_gens = n.generators.index[n.generators.bus.isin(home_buses)]
    generated = n.generators_t.p[home_gens].sum(axis=1)
    carried = pd.Series(0.0, index=n.snapshots)
    for gen in home_gens:
        factor = (grid_hourly if gen == "grid_import"
                  else factors[_carrier_key(n.generators.at[gen, "carrier"])])
        carried += n.generators_t.p[gen] * factor
    home_intensity = (carried / generated.where(generated > 0)).fillna(0.0)

    electricity_mwh = {}
    electricity_by_user = {}
    for user in ELECTRICITY_USERS:
        if user not in n.links.index:
            continue
        draw, bus = _electricity_draw(n, user)
        intensity = destination_intensity if bus == destination_bus else home_intensity
        electricity_mwh[user] = float(draw.sum()) * annual
        electricity_by_user[user] = float((draw * intensity).sum()) * annual
    emissions = dict(electricity_by_user)
    electricity_t = sum(electricity_by_user.values())

    # Gas combustion lands on whichever link burns it, so a blended shaft carries
    # its own gas rather than having it broken out beside it.
    gas_t_per_mwh = (natural_gas_cfg["co2_t_per_mwh"]
                     + emissions_cfg["gas_upstream_t_co2e_per_mwh"][basis])
    gas_t = 0.0
    for link in n.links.index[n.links.bus0 == "gas"]:
        burned = float(n.links_t.p0[link].sum()) * annual * gas_t_per_mwh
        emissions[link] = emissions.get(link, 0.0) + burned
        gas_t += burned

    # Freight over the run's own legs, each mode at its own factor.
    freight = emissions_cfg["freight_kg_co2e_per_t_km"]
    legs = transport_legs or {}
    t_co2e_per_t = sum(freight[mode] * km for mode, km in legs.items()) / 1000.0
    freight_t = 0.0
    for link in ("iron_transport", "steel_transport"):
        if link not in n.links.index:
            continue
        shipped = float(n.links_t.p0[link].sum()) * annual
        emissions[link] = shipped * t_co2e_per_t
        freight_t += emissions[link]

    # Round-trip and line losses are electricity nobody consumed, so they belong
    # to no step. Valued at the year's mean intensity rather than hour by hour:
    # a loss is incurred across the charge and the discharge, and splitting it
    # between them would be a guess dressed up as arithmetic.
    mean_intensity = float(home_intensity.mean())
    if "battery" in n.storage_units.index:
        lost = -float(n.storage_units_t.p["battery"].sum()) * annual
        emissions["battery_losses"] = lost * mean_intensity
        electricity_t += emissions["battery_losses"]
    hvdc = [link for link in n.links.index if link.startswith("hvdc_")]
    if hvdc:
        lost = float((n.links_t.p0[hvdc].sum(axis=1)
                      + n.links_t.p1[hvdc].sum(axis=1)).sum()) * annual
        emissions["transmission_losses"] = lost * mean_intensity
        electricity_t += emissions["transmission_losses"]

    return {
        "by_step": emissions,
        "electricity_mwh": electricity_mwh,
        "electricity_t": electricity_by_user,
        "sources": {"electricity": electricity_t, "gas": gas_t, "freight": freight_t},
        "grid_source": grid_source,
    }


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


def mark_best_in_country(df: pd.DataFrame, parents: dict, metric: str | None) -> pd.DataFrame:
    """Add `country`, and flag each country's best zone for every route.

    A country that supplies its market through zones (Australia through its NEM
    regions) is run once per zone, so several rows describe the same place. Every
    route picks the zone where that route came out cheapest: Australia's moe-eaf
    is reported from whichever NEM region made steel cheapest with moe-eaf, and
    each date range is ranked on its own.

    The losers are flagged, not dropped — this frame is the diagnostic one, where
    what the other zones cost is itself a result and a dropped row would cost a
    re-solve to recover. `write_report` is what narrows it down for viz.

    `metric` names the ranking column, or None to flag everything. A row with no
    value for it (h2-only produces hydrogen, so it has no cost of steel) stays
    flagged rather than being ranked away.
    """
    out = df.copy()
    out.insert(2, "country", out["area"].map(lambda a: parents.get(a, a)))
    if metric is None or metric not in out.columns:
        out["best_in_country"] = True
        return out
    ranked = out.groupby(["country", "route", "start_date", "end_date"])[metric]
    out["best_in_country"] = (out[metric] == ranked.transform("min")) | out[metric].isna()
    return out


def input_variants(scenarios: pd.DataFrame, scenario_name: str, run: dict) -> dict:
    """Which series each tech contributed to one run, as `{tech}_variant` fields.

    A best-site P95 profile and an area average answer different questions, so a
    cost is not interpretable without knowing which produced it. The solve reads
    the parquet, not its name, so the variant survives only in the scenario
    table and is joined back on here. Rows join on the run minus the route —
    every route of a group is built from the same series. The scenario table
    spells a tech `wind-onshore`; the field it lands in is `wind_onshore_variant`.
    """
    group = scenarios[
        (scenarios["scenario"] == scenario_name)
        & (scenarios["area"] == run["area"])
        & (scenarios["start_date"] == run["start_date"])
        & (scenarios["end_date"] == run["end_date"])
    ]
    variants = {f"{field_stem(row.tech)}_variant": row.variant for row in group.itertuples()}
    return variants


def write_report(df: pd.DataFrame, report_path: Path, diagnostic_path: Path) -> None:
    """Write the scenario report and the hidden diagnostic beside it.

    The report is the seam viz reads: which zone represents its country is
    already decided here, so it holds one column per reported place and no
    `best_in_country` row to interpret. The diagnostic keeps every zone and the
    flag, for the question "what would the others have cost?".
    """
    write_report_file(df, diagnostic_path)
    selected = df[df["best_in_country"]].drop(columns="best_in_country")
    write_report_file(selected, report_path)
    log.info(f"wrote {report_path} ({len(selected)} runs) and {diagnostic_path} ({len(df)} runs)")


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
    # Bought power, not built power: these are priced, so they are not part of
    # the renewables the run chose to build.
    non_res = ("grid_import", "gas_supply", "destination_supply")
    res_idx = gens.index[gens.p_nom_extendable & ~gens.index.isin(non_res)]

    def gen_capital(name: str) -> float:
        if name not in gens.index or not gens.at[name, "p_nom_extendable"]:
            return 0.0
        return float(gens.at[name, "capital_cost"] * gens.at[name, "p_nom_opt"])

    def gen_marginal(name: str) -> float:
        if name not in gens.index:
            return 0.0
        return _marginal_costs(
            gens.loc[[name]], n.generators_t.p, n.generators_t.marginal_cost
        ) * annual_scale

    grid_capital = gen_capital("grid_import")
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
            gens.drop(index=["gas_supply", "destination_supply"], errors="ignore"),
            n.generators_t.p,
            n.generators_t.marginal_cost,
        ) * annual_scale,
        # Gas bill incl. any carbon price (both live on the gas_supply
        # generator's marginal cost).
        "gas": gen_marginal("gas_supply"),
        "electrolyser": link_capital(["electrolyser"]) + link_marginal(["electrolyser"]),
        "h2_buffer": store_capital("h2_buffer"),
        "process": link_capital(PROCESS_LINKS),
        "ore_consumables": link_marginal(PROCESS_LINKS),
        "iron_store": store_capital("iron_store"),
        "steel_store": store_capital("steel_store"),
        "transmission": link_capital(hvdc),
        # Freight has no capex — the whole bill is the per-t-km charge. Only
        # one of the two links ever exists: an export route ships its iron, and
        # every other route ships its steel, if it delivers at all.
        "transport": link_marginal(["iron_transport", "steel_transport"]),
        # The destination furnace's own power, kept out of the grid group so an
        # export run shows what it pays at each end.
        "destination_power": (gen_capital("destination_supply")
                              + gen_marginal("destination_supply")),
    }


def extract_summary(
    n: pypsa.Network, scenario_name: str, run: dict, assumptions: dict,
    transport_legs: dict | None = None, grid_mix: pd.DataFrame | None = None,
) -> dict:
    """Key sizing, cost and emission metrics as a flat dict (one row of the CSV).

    `run` identifies the row — area, route, date range and the input variants —
    and leads the columns, so what a number describes reads before the number.

    The headline levelised cost depends on the network's route: LCOH for the
    pure-H2 model (flat H2 load, no steel chain), LCOS for the steel routes
    (flat steel load). H2 production is reported whenever an electrolyser
    exists, but LCOH is only well-defined when H2 is the end product — on
    steel routes the annual cost covers the whole chain.

    `transport_legs` and `grid_mix` are the run's own freight legs and the
    generation mix behind its imports; both only feed the emission fields, and
    both may be absent — a route that ships nothing has no legs, and a grid
    series solved on `dayahead` has no mix.
    """
    breakdown = _cost_breakdown(n)
    summary = {
        "scenario": scenario_name,
        **run,
        "total_annual_cost_meur": sum(breakdown.values()) / 1e6,
    }
    # Per-group annual cost columns, zeros written out: they stack to the total,
    # so a group this route has no component for contributed nothing — which is
    # a result, not a gap. plot_lcos_bars stacks these.
    for group, value in breakdown.items():
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

    # The levelised cost of whatever this run produces — steel for the steel
    # routes, hydrogen for h2-only. One column, so ranking a route's runs against
    # each other needs no per-route special-casing; the unit rides along so a row
    # is readable on its own.
    if "lcos_eur_per_t" in summary:
        summary["lco_output"] = summary["lcos_eur_per_t"]
        summary["lco_output_unit"] = "EUR/t steel"
    elif "lcoh_eur_per_kg" in summary:
        summary["lco_output"] = summary["lcoh_eur_per_kg"]
        summary["lco_output_unit"] = "EUR/kg H2"

    # Levelised cost of the underlying energy carriers, €/MWh, so LCOS can be read
    # against the electricity and hydrogen that drive it.
    #   LCOE = electricity-system cost (renewables + battery + grid + transmission)
    #          per MWh of electricity generated (renewable dispatch + grid import).
    #   LCOH = (electrolyser capex/opex + H2 buffer + the electrolyser's electricity
    #          valued at LCOE) per MWh of H2 produced, LHV.
    annual_scale = 8760.0 / len(n.snapshots)
    elec_gens = [g for g in n.generators.index if g != "gas_supply"]
    elec_mwh = (float(n.generators_t.p[elec_gens].sum().sum()) * annual_scale) if elec_gens else 0.0
    elec_cost = sum(breakdown[k] for k in ("res", "battery", "grid", "transmission"))
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
            **{f"lcoe_{field_stem(tech)}": _res_cost(tech) for tech in RES_TECHS},
            "lcoe_storage": breakdown["battery"],
            "lcoe_grid_connection": grid_conn,
            "lcoe_grid_energy": grid_energy,
            "lcoe_transmission": breakdown["transmission"],
        }
        res_total = sum(components[f"lcoe_{field_stem(tech)}"] for tech in RES_TECHS)
        summary["lcoe_renewables_eur_per_mwh"] = res_total / elec_mwh
        for key, cost in components.items():
            summary[f"{key}_eur_per_mwh"] = cost / elec_mwh

        # Per-technology LCOE (€/MWh over that tech's own generation, not the
        # system total) — a fair unit cost for the renewable itself, distinct from
        # its lcoe_<tech> contribution to the system LCOE above.
        for tech in RES_TECHS:
            cost, mwh = _res_cost(tech), _res_mwh(tech)
            if cost > 0 and mwh > 0:
                summary[f"lcoe_{field_stem(tech)}_own_eur_per_mwh"] = cost / mwh

        # Blended renewable own LCOE (over renewable generation only) — the
        # generation-weighted mean of the per-tech own LCOEs, so it sits between
        # them. Distinct from lcoe_renewables (contribution over all electricity,
        # which grid imports drag below the per-tech figures).
        res_mwh = sum(_res_mwh(tech) for tech in RES_TECHS)
        if res_total > 0 and res_mwh > 0:
            summary["lcoe_renewables_own_eur_per_mwh"] = res_total / res_mwh

        # Capacity factors (annual generation / nameplate × 8760) for the renewables
        # and the grid connection, surfaced next to the built capacities.
        def _res_cap(prefix: str) -> float:
            idx = [g for g in res_idx if str(g).startswith(prefix)]
            return float(n.generators.loc[idx, "p_nom_opt"].sum()) if idx else 0.0

        for tech in RES_TECHS:
            cap = _res_cap(tech)
            if cap > 0:
                summary[f"cf_{field_stem(tech)}"] = _res_mwh(tech) / (cap * 8760.0)
        if grid_mwh > 0:
            grid_p_nom = float(n.generators.at["grid_import", "p_nom_opt"])
            if grid_p_nom > 0:
                summary["cf_grid_connection"] = grid_mwh / (grid_p_nom * 8760.0)

        # Split the average priced grid energy into the day-ahead market price and
        # the constant volumetric fee, and levelise the connection capex over the
        # same imported-MWh denominator — so all three are €/MWh *imported* and add
        # up to the all-in delivered electricity price (not mixed bases).
        if grid_mwh > 0:
            if grid_energy > 0:
                fee = float(assumptions["grid"]["fee_eur_per_mwh"])
                summary["grid_price_eur_per_mwh"] = max(grid_energy / grid_mwh - fee, 0.0)
                summary["grid_fee_eur_per_mwh"] = fee
            if grid_conn > 0:
                summary["grid_connection_eur_per_mwh_imported"] = grid_conn / grid_mwh

    if "electrolyser" in n.links.index and "steel_load" in n.loads.index:
        el_mwh = float(n.links_t.p0["electrolyser"].sum()) * annual_scale
        h2_mwh = _h2_produced_kg(n) * H2_LHV_KWH_PER_KG / 1000.0
        if h2_mwh > 0:
            elec_for_h2 = el_mwh * (lcoe if lcoe == lcoe else 0.0)
            el_cost = breakdown["electrolyser"]
            buf_cost = breakdown["h2_buffer"]
            lcoh = (el_cost + buf_cost + elec_for_h2) / h2_mwh
            summary["lcoh_eur_per_mwh_lhv"] = lcoh
            # LCOH decomposition, €/MWh LHV over the same H2 denominator so the
            # parts sum back to LCOH: the electrolyser plant, its electricity
            # (valued at LCOE) and — where built — the H2 buffer store.
            summary["lcoh_electrolyser_eur_per_mwh_lhv"] = el_cost / h2_mwh
            summary["lcoh_electricity_eur_per_mwh_lhv"] = elec_for_h2 / h2_mwh
            if buf_cost > 0:
                summary["lcoh_h2_storage_eur_per_mwh_lhv"] = buf_cost / h2_mwh

    # Cost detail for the hover: per-plant levelised capex and the ore-vs-consumables
    # split, all €/t steel on the same denominator as the LCOS cost groups — so the
    # per-plant figures sum to the "process" group and ore+consumables to the
    # "ore_consumables" group. Ore is the feedstock priced on the reduction/
    # electrolysis links; consumables is the EAF's marginal cost.
    if "steel_load" in n.loads.index:
        steel_t = float(n.loads.at["steel_load", "p_set"]) * 8760.0
        if steel_t > 0:
            def _link_capex(link_name: str) -> float:
                return (float(n.links.at[link_name, "capital_cost"] * n.links.at[link_name, "p_nom_opt"])
                        if link_name in n.links.index and bool(n.links.at[link_name, "p_nom_extendable"])
                        else 0.0)

            def _link_marginal(link_name: str) -> float:
                if link_name not in n.links.index or link_name not in n.links_t.p0.columns:
                    return 0.0
                return float((n.links_t.p0[link_name] * n.links.at[link_name, "marginal_cost"]).sum()) * annual_scale

            for link in PROCESS_LINKS:
                summary[f"plant_{field_stem(link)}_eur_per_t"] = _link_capex(link) / steel_t
            ore = sum(_link_marginal(l) for l in ("dri-h2", "dri-ng", "moe", "ew"))
            summary["ore_eur_per_t_steel"] = ore / steel_t
            summary["consumables_eur_per_t_steel"] = _link_marginal("eaf") / steel_t

    # Steel-route process links: capacity in output units (t/h of iron or
    # steel — p_nom is input-side, so scale by the link efficiency) plus
    # utilisation, which shows how far each step actually load-follows.
    for link in PROCESS_LINKS:
        if link not in n.links.index:
            continue
        p_nom = n.links.at[link, "p_nom_opt"]
        summary[f"{field_stem(link)}_t_per_h_opt"] = p_nom * n.links.at[link, "efficiency"]
        summary[f"{field_stem(link)}_utilization"] = (
            float(n.links_t.p0[link].mean() / p_nom) if p_nom > 0 else float("nan")
        )

    if "gas_supply" in n.generators.index:
        gas_mwh = float(n.generators_t.p["gas_supply"].sum()) * (8760.0 / len(n.snapshots))
        summary["ng_gwh_lhv"] = gas_mwh / 1e3
    # How much of the iron came from the H2 shaft (production share, not capacity
    # share). Emitted for any DRI route so a pure H2-DRI reads 1.0 and a pure
    # NG-DRI 0.0 — not a missing value that downstream would coerce to 0.
    # On mix-dri-eaf the split happens on the reductant bus, one shaft
    # upstream; on the single-fuel routes it is the shafts themselves.
    h2_link, ng_link = (("reductant-h2", "reductant-ng")
                        if "dri-mix" in n.links.index else ("dri-h2", "dri-ng"))
    if h2_link in n.links.index or ng_link in n.links.index:
        from_h2 = -float(n.links_t.p1[h2_link].sum()) if h2_link in n.links.index else 0.0
        from_ng = -float(n.links_t.p1[ng_link].sum()) if ng_link in n.links.index else 0.0
        total = from_h2 + from_ng
        summary["iron_from_h2_share"] = from_h2 / total if total else float("nan")

    if "iron_store" in n.stores.index:
        store_t = n.stores.at["iron_store", "e_nom_opt"]
        summary["iron_store_kt"] = store_t / 1e3
        if "steel_load" in n.loads.index:
            steel_t_per_h = float(n.loads.at["steel_load", "p_set"])
            summary["iron_store_hours_steel"] = (
                store_t / steel_t_per_h if steel_t_per_h else float("nan")
            )

    annual = 8760.0 / len(n.snapshots)
    for link, field in (("iron_transport", "iron_shipped_kt"),
                        ("steel_transport", "steel_shipped_kt")):
        if link in n.links.index:
            summary[field] = float(n.links_t.p0[link].sum()) * annual / 1e3
            # The distance behind the freight bill, so the report carries it
            # rather than the reader going back to the assumptions.
            summary["transport_km"] = float(n.links.at[link, "length"])

    if "steel_store" in n.stores.index:
        store_t = n.stores.at["steel_store", "e_nom_opt"]
        summary["steel_store_kt"] = store_t / 1e3
        if "steel_load" in n.loads.index:
            steel_t_per_h = float(n.loads.at["steel_load", "p_set"])
            summary["steel_store_hours_steel"] = (
                store_t / steel_t_per_h if steel_t_per_h else float("nan")
            )

    for gen in n.generators.index[n.generators.p_nom_extendable]:
        summary[f"{field_stem(gen)}_gw_opt"] = n.generators.at[gen, "p_nom_opt"] / 1e3

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
        elif "dri-h2" in n.links.index:
            h2_demand_mw = float(n.links_t.p0["dri-h2"].mean())
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
            summary[f"{field_stem(link)}_gw_opt"] = cap / 1e3
            trans_cost += n.links.at[link, "capital_cost"] * cap
        summary["transmission_total_annual_cost_meur"] = trans_cost / 1e6

    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    if len(ac_buses) > 1:
        ext = n.generators[n.generators.p_nom_extendable]
        for carrier, grp in ext.groupby("carrier"):
            summary[f"{field_stem(carrier)}_total_gw_opt"] = grp["p_nom_opt"].sum() / 1e3

    # Emissions, on the basis the assumptions name. Reported in kg so the
    # report's two decimals still say something about a clean route.
    emitted = _emissions_breakdown(
        n, assumptions["emissions"], assumptions["natural_gas"],
        run["area"], transport_legs, grid_mix,
    )
    emissions = emitted["by_step"]
    electricity_mwh = emitted["electricity_mwh"]
    sources = emitted["sources"]
    summary["emissions_basis"] = assumptions["emissions"]["basis"]
    summary["emissions_grid_source"] = emitted["grid_source"]
    emissions_t = sum(emissions.values())
    summary["emissions_kt_co2e_per_year"] = emissions_t / 1e3
    for source, value in sources.items():
        summary[f"emissions_{source}_kt_co2e_per_year"] = value / 1e3

    # Per tonne of steel, and each step's share of it. The steps stack to the
    # total, so a step this route has none of reads 0 rather than blank.
    if "steel_produced_mt" in summary and summary["steel_produced_mt"] > 0:
        steel_t = summary["steel_produced_mt"] * 1e6
        summary["emissions_kg_co2e_per_t_steel"] = emissions_t * 1e3 / steel_t
        for step in EMISSION_STEPS:
            summary[f"emissions_{field_stem(step)}_kg_co2e_per_t_steel"] = (
                emissions.get(step, 0.0) * 1e3 / steel_t
            )
    if emissions_t > 0:
        for step in EMISSION_STEPS:
            summary[f"emissions_{field_stem(step)}_pct"] = (
                emissions.get(step, 0.0) * 100.0 / emissions_t
            )
    # The hydrogen's own emissions over the hydrogen, not the plant's over the
    # hydrogen: this is the figure that compares to the RFNBO ceiling of
    # 3.38 kg CO2e per kg H2, and only the electrolyser's electricity is in the
    # hydrogen's production chain. On a steel route the furnace is downstream of
    # it and belongs to the steel.
    if "electrolyser" in n.links.index:
        h2_kg = _h2_produced_kg(n)
        if h2_kg > 0:
            summary["emissions_kg_co2e_per_kg_h2"] = (
                emitted["electricity_t"].get("electrolyser", 0.0) * 1e3 / h2_kg
            )

    # Electricity by who drew it, and how dirty their own hours were — the number
    # that separates a user chasing cheap renewable hours from one running flat.
    total_el_mwh = sum(electricity_mwh.values())
    if total_el_mwh > 0:
        summary["emissions_kg_co2e_per_mwh_el"] = (
            sources["electricity"] * 1e3 / total_el_mwh
        )
    for user, mwh in electricity_mwh.items():
        summary[f"el_{field_stem(user)}_gwh"] = mwh / 1e3
        if mwh > 0:
            summary[f"emissions_{field_stem(user)}_kg_co2e_per_mwh_el"] = (
                emitted["electricity_t"][user] * 1e3 / mwh
            )

    return summary


def main() -> None:
    """Load every run under the scenario and write the combined report CSV.

    A scenario is an umbrella over runs, so the netCDF name carries the rest of
    the run key — area, route and date range; the scenario table adds the
    variant each tech was solved with. One summary row per run via
    `extract_summary`, then the report and its diagnostic sibling.
    """
    scenario_name = snakemake.wildcards.scenario
    scenarios = load_scenarios(snakemake.input.scenarios, snakemake.config["areas"])
    assumptions = yaml.safe_load(Path(snakemake.input.assumptions).read_text())
    parents = zone_parents(snakemake.config["areas"])

    # The grid series each area was solved against, by area. Only the emission
    # fields read them, and only for the generation mix: a series solved on
    # `dayahead` has none, and an islanded scenario has no series at all.
    grid_paths = {Path(p).stem.split("_")[0]: Path(p)
                  for p in snakemake.input.get("grid_input", [])}

    rows = []
    network_paths = list(dict.fromkeys(snakemake.input.networks))
    log.info(f"compiling report for scenario={scenario_name} ({len(network_paths)} runs)")
    for nc_path in network_paths:
        nc_path = Path(nc_path)
        area, route, start_date, end_date = nc_path.stem.split("_")
        n = pypsa.Network()
        n.import_from_netcdf(nc_path)
        run = {"area": area, "route": route, "start_date": start_date, "end_date": end_date}
        run.update(input_variants(scenarios, scenario_name, run))
        # Freight legs are the country's, not the zone's — a NEM region ships
        # from Australia. Same resolution as solve_network's.
        legs = assumptions["transport"]["distance_km"].get(parents.get(area, area))
        grid_mix = (pd.read_parquet(grid_paths[area]) if area in grid_paths else None)
        summary = extract_summary(n, scenario_name, run, assumptions, legs, grid_mix)
        # Trailing the row, because it identifies the inputs rather than
        # describing them: the per-file map it stands for is in the network.
        summary["inputs_hash"] = n.meta["inputs_hash"]
        rows.append(summary)

    flagged = mark_best_in_country(
        pd.DataFrame(rows), parents, snakemake.params.best_zone_by or None,
    )
    write_report(flagged, Path(snakemake.output.report), Path(snakemake.output.diagnostic))


if __name__ == "__main__":
    main()
