"""Compile one solved run into its report.

Invoked by Snakemake's `script:` directive (compile_report rule in viz.smk).
Writes two files: the report viz reads, and a hidden `_diag.csv` beside it that
also keeps the fields the report holds back. Each run carries the `inputs_hash`
the solve stamped into the network, so a number stays tied to its inputs. What
the files look like is common/_report_schema.py's business.
"""

import calendar
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import yaml

from common._constants import H2_LHV_KWH_PER_KG
from common._logging import configure_logging
from common._report_schema import (
    DIAGNOSTIC_FIELDS,
    ELECTRICITY_USERS,
    IRON_ELECTRICITY_JOBS,
    EMISSION_STEPS,
    FREIGHT_LEGS,
    LEAF_COSTS,
    LEAF_GROUP,
    LEAF_PARENTS,
    LEAF_SUBTOTALS,
    ORE_LINKS,
    ORIGIN_ELECTRICITY_LEAVES,
    REDUCTION_LINKS,
    PROCESS_LINKS,
    RES_TECHS,
    field_stem,
    variant_field,
    write_report_file,
)
from scripts.solve._helpers_solve import annuity_factor, deep_merge, find_transport_legs, read_run_inputs

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)

# The grid carriers that store rather than generate, as both downloaders name
# them. They are not a source of their own, so they are left out of a mix rather
# than counted in it at zero — see `_grid_intensity`.
STORAGE_CARRIERS = ("energy_storage", "pumped_storage")

# ENTSO-E publishes an "Actual Consumption" figure beside the generation of any
# carrier that has one, which the downloader negates and suffixes `_cons`. It is
# a plant's own draw and not a source, whatever the carrier: Germany reports one
# for solar and onshore wind as well as for its storage.
CONSUMPTION_SUFFIX = "_cons"

# A link's second electricity draw (on bus3), and the user it is reported under.
# The H2 shaft preheats its own feed; the blend shaft's preheat is a link of its
# own, `reductant-h2`, and both land in the same row.
SECOND_DRAWS = {"dri-h2": "reductant-h2"}


def _carrier_key(carrier: str) -> str:
    """The emission table's key for a generator carrier: `wind-onshore` → `wind_onshore`.

    Orientation-suffixed keys fall back to their base tech, the same way
    `build_network._add_generators` resolves them against the res assumptions.
    """
    return field_stem(re.sub(r"_(east|west)_\d+$|_az\d+$", "", carrier))


def _hours_per_year(snapshots: pd.DatetimeIndex) -> float:
    """Mean length in hours of the calendar years the snapshots fall in."""
    year_hours = [8784.0 if calendar.isleap(year) else 8760.0 for year in snapshots.year.unique()]
    hours_per_year = float(np.mean(year_hours))
    return hours_per_year


def _grid_intensity(
    area: str, grid_mix: pd.DataFrame | None,
    factors: dict[str, float], snapshots: pd.Index,
) -> pd.Series:
    """t CO2e per MWh imported, hour by hour, from the area's own generation mix.

    Returns None when the series carries no mix at all. A grid series solved on
    `variant: emissions` carries the area's generation by carrier, and those
    column names are the factor table's keys — which is what the downloaders'
    shared vocabulary was for. Not every market publishes one, and nothing stands
    in for it: the report says the intensity is unknown rather than inventing a
    number. A `full` series carries the carriers too, but at native resolution
    and beside load and flow columns, so it is not a mix source either.

    Storage carriers are left out of the mix altogether rather than counted at
    zero: what a reservoir gives back was generated in some earlier hour that is
    already in here, so excluding both its discharge and its charging is what
    makes stored energy carry the average of everything else. Counting the
    discharge at a factor of zero would instead dilute the mix.

    A `_cons` column is left out for a related reason: it is what a carrier's own
    plants drew rather than anything they made, so the mix is over gross
    generation. It is a few MW against tens of GW either way.

    Production-based: a zone importing coal power from next door reads as clean
    as its own plants. The cross-border columns are in the `full` series if that
    is ever worth fixing.
    """
    columns = [col for col in (grid_mix.columns if grid_mix is not None else [])
               if col != "price"
               and not col.endswith(CONSUMPTION_SUFFIX)
               and col not in STORAGE_CARRIERS]
    if not columns:
        # No mix to read, which is the source's limit rather than a mistake in the
        # scenario table: ONS publishes an undivided thermal aggregate that no
        # single factor describes, and AESO and IESO publish no generation at all.
        # The caller substitutes NaN and states the reason in the report.
        return None
    unknown = [col for col in columns if col not in factors]
    if unknown:
        raise ValueError(
            f"{area}: its grid series has generation columns the emission factor "
            f"table has no entry for: {sorted(unknown)}. Either they are carriers "
            f"and belong in the `electricity_t_co2e_per_mwh` table under "
            f"`emissions` in "
            f"the assumptions, or the series is not a `variant: emissions` "
            f"one — a `full` series carries load and flow columns too."
        )
    missing = pd.Index(snapshots).difference(grid_mix.index)
    if len(missing):
        raise ValueError(
            f"{area}: its grid series does not cover {len(missing)} of the run's "
            f"snapshots (first {missing[0]}, last {missing[-1]}). The series and the "
            f"solved network are for different windows, so pairing them would report "
            f"an intensity from the wrong year."
        )
    # ENTSO-E publishes the odd negative generation hour; a negative would
    # otherwise pull the weighted mean the wrong way.
    generation = grid_mix[columns].clip(lower=0.0)
    generated = generation.sum(axis=1)
    carried = sum(generation[carrier] * factors[carrier] for carrier in columns)
    # An hour the zone generated nothing at all is an hour with no mix to read,
    # so it carries the neighbouring hours' rather than a fabricated zero.
    hourly = (carried / generated.where(generated > 0)).reindex(snapshots)
    return hourly.ffill().bfill()


def _storage_carbon(
    storage_flow: pd.Series, closing_level: float, eta_in: float, eta_out: float,
    generation_intensity: pd.Series,
) -> tuple[pd.Series, float]:
    """What storage gave back each hour, in t CO2e, and what it lost over the year.

    `storage_flow` is signed as seen from the electricity bus — positive when the
    battery supplies it, negative when it charges — and `closing_level` is the
    energy left in the store at the final snapshot.

    A battery emits nothing of its own: what comes out carries what went in,
    mixed with whatever was in the tank already. So the tank is followed hour by
    hour — charging adds that hour's carbon and dilutes what was there,
    discharging draws at the blend. The state of charge is the memory, which is
    why there is no averaging window to pick: a tank that turns over daily
    forgets in a day, a seasonal one carries the winter into the spring.

    Both round-trip losses are charged where they happen, at the intensity of
    the energy that suffered them — the charging loss at the hour's generation,
    the discharge loss at the tank's blend — so what the users are charged plus
    what the losses are charged is what the charging drew.

    Two passes, because the state of charge is cyclic: the carbon in the tank at
    the start of the year is whatever the year ends with. The tank mixes, so by
    the second pass it has forgotten where the first one started.
    """
    hourly_intensity = generation_intensity.to_numpy()
    flow = storage_flow.to_numpy()
    carbon = 0.0
    for pass_number in (1, 2):
        opening = carbon
        level = closing_level
        given_back = np.zeros(len(flow))
        lost = 0.0
        for hour in range(len(flow)):
            drawn = max(-flow[hour], 0.0)
            if drawn:
                carbon += drawn * eta_in * hourly_intensity[hour]
                lost += drawn * (1.0 - eta_in) * hourly_intensity[hour]
                level += drawn * eta_in
            injected = max(flow[hour], 0.0)
            if injected:
                taken = injected / eta_out
                blend = carbon / level if level > 0 else 0.0
                carbon -= taken * blend
                level -= taken
                given_back[hour] = injected * blend
                lost += (taken - injected) * blend
        if pass_number == 1:
            # Open the year at the intensity it closed at — the fixed point
            # when the state of charge is cyclic, which is how it is solved.
            carbon = carbon / level * closing_level if level > 0 else 0.0
    # Carbon the tank still holds over what it opened with was paid for and
    # never used, so it is a loss like the others. Nil on a cyclic year,
    # which is what the model solves.
    return pd.Series(given_back, index=storage_flow.index), lost + (carbon - opening)


def _emissions_breakdown(
    n: pypsa.Network, emissions_cfg: dict,
    area: str, transport_legs: dict | None, grid_mix: pd.DataFrame | None,
    destination_area: str, destination_mix: pd.DataFrame | None,
) -> dict[str, object]:
    """What the run emitted in a year, and who to charge it to.

    `by_step` is annual t CO2e per step and stacks to the total. `electricity_mwh`
    and `electricity_t` are each user's draw and the emissions of that draw alone
    — kept apart from `by_step` because a gas-fired shaft's step total carries its
    combustion too, and dividing that by its MWh would report a furnace as
    buying impossibly dirty power. `losses_mwh` is the electricity the round trip
    and the lines took, which no user drew. `sources` splits the total into
    electricity and gas.

    `freight_by_leg` is what delivering the output emitted, and is no part of
    `by_step` or `sources`: the boundary is the plant, and how far a customer
    happens to be is a question asked of a route rather than a property of it.
    It is also the one figure here that does not need a grid mix, so it stands
    when the rest cannot be read.

    Accounting only: nothing here reaches the objective (the optional carbon
    price on gas is set in build_network), so nothing here can move a solve.
    And it covers the run's energy alone — the process steps' own
    direct emissions are outside the model boundary, which is why the result is
    not a CBAM or an ETS figure. See the `emissions` block of the
    assumptions.

    Electricity is attributed hour by hour at the intensity of the system that
    supplied it, so a user that ran when the wind blew carries less than one
    that ran flat out. Two systems can exist: the producing area's, where the
    plant's own renewables mix with whatever it imports, and — on an export
    route — the destination's, which is a grid and nothing else.
    """
    factors = emissions_cfg[emissions_cfg["basis"]]["electricity_t_co2e_per_mwh"]
    annual = _hours_per_year(n.snapshots) / len(n.snapshots)

    # PyPSA's netCDF export drops a time-varying column whose every value is the
    # default, so a component that never ran comes back missing from the
    # dispatch frame rather than zero — an optimum that builds no battery is an
    # ordinary outcome, not a reason for the whole report to fail.
    zero_series = pd.Series(0.0, index=n.snapshots)
    generator_p = n.generators_t.p.reindex(columns=n.generators.index, fill_value=0.0)
    link_p0 = n.links_t.p0.reindex(columns=n.links.index, fill_value=0.0)
    link_p1 = n.links_t.p1.reindex(columns=n.links.index, fill_value=0.0)
    store_e = n.stores_t.e.reindex(columns=n.stores.index, fill_value=0.0)
    # bus2 is optional, so a route whose links all take their electricity on
    # bus0 has neither the column nor the frame.
    link_p2 = (n.links_t["p2"].reindex(columns=n.links.index, fill_value=0.0)
               if "p2" in n.links_t else link_p0 * 0.0)
    link_bus2 = (n.links["bus2"] if "bus2" in n.links.columns
                 else pd.Series("", index=n.links.index))
    # And bus3, which only the H2 shaft uses: the preheat on its H2 feed.
    link_p3 = (n.links_t["p3"].reindex(columns=n.links.index, fill_value=0.0)
               if "p3" in n.links_t else link_p0 * 0.0)
    link_bus3 = (n.links["bus3"] if "bus3" in n.links.columns
                 else pd.Series("", index=n.links.index))

    # The destination is found by who feeds it rather than by name, so the bus
    # can be renamed in build_network without breaking the report.
    destination_bus = (n.generators.at["destination_supply", "bus"]
                       if "destination_supply" in n.generators.index else None)
    if destination_bus is not None and destination_mix is None:
        raise ValueError(
            f"this run melts its iron in {destination_area}, but no market series for "
            f"{destination_area} over the run's own window reached the report, so what "
            f"its furnace emitted is unknown. The run's grid row for its destination "
            f"is what carries it (`grid_input` in rule compile_report)."
        )
    # Its furnace carries the destination market's own hours, read the same way
    # as any other grid's — so an export route's melt is as clean as the country
    # it melts in was in the hours it ran, and no differently.
    # An intensity that cannot be read becomes NaN rather than an exception, so
    # everything downstream of it lands blank while the MWh either side of it
    # stay real. `unavailable` carries the reason into the report.
    unavailable = None
    destination_intensity = (
        _grid_intensity(destination_area, destination_mix, factors, n.snapshots)
        if destination_bus is not None else None
    )
    if destination_bus is not None and destination_intensity is None:
        unavailable = (
            f"{destination_area} publishes no generation mix, so what this route's "
            f"furnace emitted melting its iron there is unknown"
        )
        destination_intensity = float("nan")

    home_intensity_series = None
    if "grid_import" in n.generators.index:
        home_intensity_series = _grid_intensity(area, grid_mix, factors, n.snapshots)
        if home_intensity_series is None:
            unavailable = unavailable or (
                f"{area} publishes prices but no generation mix, so the emission "
                f"intensity of this run's imported electricity is unknown"
            )
            home_intensity_series = float("nan")

    # What a MWh on the producing area's electricity buses carried, hour by hour:
    # its own generation at the carriers' factors, its imports at the grid's.
    home_buses = set(n.buses.index[n.buses.carrier == "AC"]) - {destination_bus}
    home_gens = n.generators.index[n.generators.bus.isin(home_buses)]
    generated = generator_p[home_gens].sum(axis=1)
    carried = pd.Series(0.0, index=n.snapshots)
    for gen in home_gens:
        factor = (home_intensity_series if gen == "grid_import"
                  else factors[_carrier_key(n.generators.at[gen, "carrier"])])
        carried += generator_p[gen] * factor

    # Storage supplies the system like anything else, carrying what charged it
    # (see `_storage_carbon`). Without it an hour supplied out of the battery
    # has nothing generating in it, and everything drawn in that hour comes out
    # free: a solar plant that runs its furnace through the night reports about
    # a third light.
    generation_intensity = (carried / generated.where(generated > 0)).fillna(0.0)
    # Signed as the electricity bus sees it: the charger's p0 is a withdrawal,
    # and the discharger's p1 is already negative leaving the link, so negating
    # it lands positive. A network built without a battery has neither link.
    has_battery = "battery_charger" in n.links.index
    charged = link_p0["battery_charger"] if has_battery else zero_series
    discharged = -link_p1["battery_discharger"] if has_battery else zero_series
    stored_carried, storage_lost_t = _storage_carbon(
        discharged - charged,
        float(store_e.at[n.snapshots[-1], "battery"]) if has_battery else 0.0,
        float(n.links.at["battery_charger", "efficiency"]) if has_battery else 1.0,
        float(n.links.at["battery_discharger", "efficiency"]) if has_battery else 1.0,
        generation_intensity,
    )
    supplied = generated + discharged
    home_intensity = ((carried + stored_carried)
                      / supplied.where(supplied > 0)).fillna(0.0)

    # Which links draw electricity is read off the network — on bus0 as their
    # input (an electrolyser, a MOE cell) or on bus2 alongside the conversion
    # they exist for (a shaft's auxiliaries, the furnace's melt) — and either way
    # PyPSA reads a positive port flow as withdrawn from that bus, so both are
    # positive as they stand. The schema's list is what the report has columns
    # for, so a new drawing link is an error here rather than a silent hole in
    # the total.
    electricity_mwh = {}
    electricity_by_user = {}
    for link in n.links.index:
        if n.links.at[link, "carrier"] == "HVDC":
            continue    # a line rather than a user; its loss is charged below
        if n.links.at[link, "carrier"] == "battery":
            continue    # storage moves a draw between hours rather than making
            # one; what it keeps is charged as battery_losses
        on_bus0 = n.buses.at[n.links.at[link, "bus0"], "carrier"] == "AC"
        bus2 = link_bus2[link]
        if not on_bus0 and not (bus2 and n.buses.at[bus2, "carrier"] == "AC"):
            continue
        if link not in ELECTRICITY_USERS:
            raise ValueError(
                f"{link} draws electricity but the report has no field for it. Add it "
                f"to ELECTRICITY_USERS in common/_report_schema.py, or its draw goes "
                f"missing from the run's total while every share still reads 100 %."
            )
        draw, bus = ((link_p0[link], n.links.at[link, "bus0"]) if on_bus0
                     else (link_p2[link], bus2))
        intensity = destination_intensity if bus == destination_bus else home_intensity
        electricity_mwh[link] = float(draw.sum()) * annual
        electricity_by_user[link] = float((draw * intensity).sum()) * annual
        # A second draw on bus3 is reported under the user it is for, not the
        # link it is on: the H2 shaft's preheat sits with the blend shaft's, so
        # the preheat has one row whichever shaft drew it.
        bus3 = link_bus3[link]
        if bus3 and n.buses.at[bus3, "carrier"] == "AC":
            if link not in SECOND_DRAWS:
                raise ValueError(
                    f"{link} draws electricity on bus3 but the report does not know "
                    f"what for. Name the user it belongs to in SECOND_DRAWS."
                )
            user = SECOND_DRAWS[link]
            intensity = destination_intensity if bus3 == destination_bus else home_intensity
            electricity_mwh[user] = electricity_mwh.get(user, 0.0) + float(link_p3[link].sum()) * annual
            electricity_by_user[user] = (electricity_by_user.get(user, 0.0)
                                         + float((link_p3[link] * intensity).sum()) * annual)
    emissions = dict(electricity_by_user)
    electricity_t = sum(electricity_by_user.values())

    # Gas combustion lands on whichever link burns it, so a blended shaft carries
    # its own gas rather than having it broken out beside it. Not scope 2: it is
    # what gives the NG-DRI comparison case a figure to hold the grid runs against.
    gas_factors = emissions_cfg["natural_gas_reductant_t_co2e_per_mwh"]
    gas_t_per_mwh = gas_factors["process"] + gas_factors["upstream"]
    gas_t = 0.0
    for link in n.links.index[n.links.bus0 == "gas"]:
        burned = float(link_p0[link].sum()) * annual * gas_t_per_mwh
        emissions[link] = emissions.get(link, 0.0) + burned
        gas_t += burned

    # Freight over the run's own legs, each mode at its own factor. A diagnostic:
    # kept out of `emissions` so that no total counts it.
    freight = emissions_cfg["freight_kg_co2e_per_t_km"]
    legs = transport_legs or {}
    t_co2e_per_t = sum(freight[mode] * km for mode, km in legs.items()) / 1000.0
    freight_by_leg = {}
    for link in FREIGHT_LEGS:
        if link not in n.links.index:
            continue
        shipped = float(link_p0[link].sum()) * annual
        freight_by_leg[link] = shipped * t_co2e_per_t

    # Round-trip and line losses are electricity nobody consumed, so they belong
    # to no step. The battery's is charged inside the tank, where the intensity
    # of what suffered it is known; a line's is hour by hour, on the system that
    # fed it that hour. Between them and the storage tank, what the users are
    # charged plus what the losses are charged is what the generation carried.
    losses_mwh = 0.0
    if has_battery:
        emissions["battery_losses"] = storage_lost_t * annual
        electricity_t += emissions["battery_losses"]
        losses_mwh += float((charged - discharged).sum()) * annual
    hvdc = list(n.links.index[n.links.carrier == "HVDC"])
    if hvdc:
        lost_hourly = link_p0[hvdc].sum(axis=1) + link_p1[hvdc].sum(axis=1)
        emissions["transmission_losses"] = (
            float((lost_hourly * home_intensity).sum()) * annual
        )
        electricity_t += emissions["transmission_losses"]
        losses_mwh += float(lost_hourly.sum()) * annual

    if unavailable:
        # What each user drew, and what the losses took, are known whatever the
        # mix was — only what it emitted is not. So the MWh stay and every t of
        # CO2e inside the boundary goes blank, the gas included: it stacks into
        # one total with the electricity, and a total that is part known and part
        # not is the number most likely to be quoted as if it were whole. The
        # freight is a distance and a factor, and no total holds it, so it stands.
        nan = float("nan")
        emissions = {step: nan for step in emissions}
        electricity_by_user = {user: nan for user in electricity_by_user}
        electricity_t = gas_t = nan

    return {
        "by_step": emissions,
        "electricity_mwh": electricity_mwh,
        "electricity_t": electricity_by_user,
        "losses_mwh": losses_mwh,
        "sources": {"electricity": electricity_t, "gas": gas_t},
        "freight_by_leg": freight_by_leg,
        "unavailable": unavailable,
    }


def _h2_produced_kg(n: pypsa.Network) -> float:
    """Annual H2 produced by the electrolyser, in kg, scaled to a year.

    Read from the electrolyser link's H2-side output, so it is what was made
    rather than what any one consumer drew.
    """
    t_hours = len(n.snapshots)
    # PyPSA Link sign convention: p1 < 0 when the link injects power into bus1,
    # so -p1 is the (positive) H2 LHV output on the hydrogen bus.
    h2_mwh_lhv = -float(n.links_t.p1["electrolyser"].sum()) * (_hours_per_year(n.snapshots) / t_hours)
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


def _cost_breakdown(n: pypsa.Network) -> dict[str, float]:
    """Annualised cost per system component group, €/yr; sums to the total.

    Capital costs are already per-year (annualised CAPEX × p_nom_opt); variable
    costs — grid imports, ore, EAF consumables, electrolyser variable opex —
    are scaled from the simulation period up to a year so levelised costs stay
    meaningful on partial-year runs. The reported total annual cost is the sum
    of these groups, so breakdown and total cannot drift apart.
    """
    t_hours = len(n.snapshots)
    annual_scale = _hours_per_year(n.snapshots) / t_hours

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
        # Energy on the store, power on the charger — the discharger carries no
        # capex, since the pair shares one inverter rating (see solve_network).
        "battery": store_capital("battery") + link_capital(["battery_charger"]),
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
        "destination_electricity": (gen_capital("destination_supply")
                              + gen_marginal("destination_supply")),
    }


def _fixed_om_share(wacc: float, cfg: dict, capex_key: str, opex_key: str) -> float:
    """Fixed O&M's share of one component's annual cost, as build_network priced it.

    A `capital_cost` in the network is `annuity x capex + fixed opex` per unit of
    capacity, so the share carries over unchanged to the annual cost of whatever
    the solve went on to build.
    """
    annual_capex = annuity_factor(wacc, cfg["lifetime_years"]) * cfg[capex_key]
    total = annual_capex + cfg[opex_key]
    return cfg[opex_key] / total if total else 0.0


def _leaf_breakdown(
    n: pypsa.Network, assumptions: dict, breakdown: dict[str, float], steel_t: float,
    electricity: dict,
) -> dict[str, float]:
    """The cost groups cut as finely as the model allows, as report fields.

    `cost_*_eur_per_t_steel` comes out, one column per priced thing plus what each
    parent group of them comes to, stacking to `lcos_eur_per_t`. So does
    `cost_electricity_*_eur_per_t_steel`, which divides one of those groups — the electricity — by
    the job each euro of it paid for. Both are checked against the total they
    stack to rather than assumed to close.

    Every leaf is one component's own annual cost off the solved network —
    `capital_cost x p_nom_opt` for what was built, `marginal_cost x dispatch` for
    what was consumed. Nothing is apportioned by a share of something else: a
    battery's power and energy halves, a grid connection against its energy, and
    each renewable against the others are each the cost the solve incurred, not a
    ratio of config quotes or of levelised contributions.

    The one place the quotes are needed is the capital/fixed-O&M line
    inside a single component, because `capital_cost` is their sum and the
    network keeps no record of the two parts. `assumptions` must therefore be the
    merged base+overlay the solve itself was given.

    `electricity` carries what `_emissions_breakdown` read off the network: each
    user's annual draw in `by_user`, the round trip and the lines in
    `losses_mwh`, and the generation both are supplied out of in
    `generation_mwh`.
    """
    annual_scale = _hours_per_year(n.snapshots) / len(n.snapshots)
    # The furnace is financed where it stands, everything before it where the
    # iron is made, exactly as build_network financed them. It stands at the
    # destination when it buys that market's power.
    origin_wacc = assumptions["finance"]["origin_wacc"]
    furnace_wacc = (assumptions["finance"]["destination_wacc"]
                    if "destination_supply" in n.generators.index else origin_wacc)
    wacc_of = lambda link: furnace_wacc if link == "eaf" else origin_wacc
    leaves = dict.fromkeys(LEAF_COSTS, 0.0)

    def link_capital(name: str) -> float:
        if name not in n.links.index or not bool(n.links.at[name, "p_nom_extendable"]):
            return 0.0
        return float(n.links.at[name, "capital_cost"] * n.links.at[name, "p_nom_opt"])

    def link_marginal(name: str) -> float:
        if name not in n.links.index or name not in n.links_t.p0.columns:
            return 0.0
        priced = n.links_t.p0[name] * n.links.at[name, "marginal_cost"]
        return float(priced.sum()) * annual_scale

    def store_capital(name: str) -> float:
        if name not in n.stores.index or not bool(n.stores.at[name, "e_nom_extendable"]):
            return 0.0
        return float(n.stores.at[name, "capital_cost"] * n.stores.at[name, "e_nom_opt"])

    def generator_capital(name: str) -> float:
        gens = n.generators
        if name not in gens.index or not bool(gens.at[name, "p_nom_extendable"]):
            return 0.0
        return float(gens.at[name, "capital_cost"] * gens.at[name, "p_nom_opt"])

    # -- feedstock: the ore quote on whichever reduction step was built, and the
    # furnace's consumables. Read off the links that carry them rather than from
    # a list of the routes that have them, so a shaft reducing with a blend pays
    # for its ore here like either single-fuel shaft does.
    leaves["ore"] = sum(link_marginal(link) for link in ORE_LINKS)
    leaves["consumables"] = link_marginal("eaf")

    # -- process plants, each cut into annualised capital and fixed O&M
    for link in PROCESS_LINKS:
        annual = link_capital(link)
        if not annual:
            continue
        share = _fixed_om_share(wacc_of(link), assumptions[link],
                                "capex_per_t_per_year_eur", "opex_per_t_per_year_eur")
        leaves[f"{field_stem(link)}_fom"] = annual * share
        leaves[f"{field_stem(link)}_capex"] = annual * (1.0 - share)

    # -- natural gas: the fuel bill, and separately any carbon price on it. Both
    # ride on the gas generator's marginal cost, so the carbon comes out at the
    # rate it was charged and the fuel is what is left of the group.
    emissions_cfg = assumptions["emissions"]
    gas_mwh = (float(n.generators_t.p["gas_supply"].sum()) * annual_scale
               if "gas_supply" in n.generators.index else 0.0)
    carbon_per_mwh = (emissions_cfg["co2_price_eur_per_t"]
                      * emissions_cfg["natural_gas_reductant_t_co2e_per_mwh"]["process"])
    leaves["gas_carbon"] = gas_mwh * carbon_per_mwh
    leaves["gas_fuel"] = breakdown["gas"] - leaves["gas_carbon"]

    # -- hydrogen: the electrolyser's water and consumables are its marginal
    # cost, so the network separates them from its capital directly, without
    # going through the hydrogen it made and the efficiency it nominally made it at.
    leaves["electrolyser_water"] = link_marginal("electrolyser")
    electrolyser_capital = link_capital("electrolyser")
    if electrolyser_capital:
        share = _fixed_om_share(origin_wacc, assumptions["electrolyser"],
                                "capex_per_mw_eur", "opex_per_mw_per_year_eur")
        leaves["electrolyser_fom"] = electrolyser_capital * share
        leaves["electrolyser_capex"] = electrolyser_capital * (1.0 - share)
    leaves["h2_buffer"] = breakdown["h2_buffer"]

    # -- the electricity system, by technology and by component
    res_idx = n.generators.index[
        n.generators.p_nom_extendable
        & ~n.generators.index.isin(("grid_import", "gas_supply", "destination_supply"))
    ]

    def res_annual(names) -> float:
        gens = n.generators.loc[list(names)]
        return float((gens["capital_cost"] * gens["p_nom_opt"]).sum()) if len(gens) else 0.0

    named = set()
    for tech in RES_TECHS:
        # A multi-site run names a generator per candidate site, so the tech is a
        # prefix rather than the whole name.
        idx = [gen for gen in res_idx if str(gen).startswith(tech)]
        named.update(idx)
        annual = res_annual(idx)
        if not annual:
            continue
        share = _fixed_om_share(origin_wacc, assumptions["res"][tech],
                                "capex_per_mw_eur", "opex_per_mw_per_year_eur")
        leaves[f"res_{field_stem(tech)}_fom"] = annual * share
        leaves[f"res_{field_stem(tech)}_capex"] = annual * (1.0 - share)
    # A generator none of the three techs names has no capex and O&M split to
    # put it in, and no row of its own: an error, rather than a cost parked
    # somewhere undivided or dropped from the stack.
    unnamed = sorted(str(gen) for gen in res_idx if gen not in named)
    if unnamed:
        raise ValueError(
            f"the run built renewable generators the report has no technology for: "
            f"{unnamed}. Add the technology to RES_TECHS in "
            f"common/_report_schema.py, which gives it its own capex and fixed O&M rows."
        )

    # The inverter rating and the energy are two components the solve sized
    # separately, so what each cost needs no assumption about duration.
    leaves["battery_power_capex"] = link_capital("battery_charger")
    leaves["battery_energy_capex"] = store_capital("battery")

    # The connection pays annuitised capex plus a yearly charge on the same built
    # MW; the imported energy pays the market price plus a flat charge on every
    # MWh. So each divides at the rate it was charged at, and the market price —
    # the one part that is neither a quote nor a capacity — is the remainder.
    # The same for both connections: the origin's, and on an export route the
    # destination furnace's, each financed where it stands.
    grid_cfg = assumptions["grid"]

    def split_connection(generator: str, wacc: float, group_total: float, prefix: str) -> None:
        connection = generator_capital(generator)
        if connection:
            annual_capex = (annuity_factor(wacc, grid_cfg["connection_lifetime_years"])
                            * grid_cfg["connection_capex_eur_per_mw"])
            yearly_fee = grid_cfg["fee_eur_per_mw_per_year"]
            fee_share = yearly_fee / (annual_capex + yearly_fee) if (annual_capex + yearly_fee) else 0.0
            leaves[f"{prefix}_capacity_fee"] = connection * fee_share
            leaves[f"{prefix}_connection_capex"] = connection * (1.0 - fee_share)
        mwh = (float(n.generators_t.p[generator].sum()) * annual_scale
               if generator in n.generators.index else 0.0)
        leaves[f"{prefix}_fee"] = mwh * float(grid_cfg["fee_eur_per_mwh"])
        leaves[f"{prefix}_market"] = group_total - connection - leaves[f"{prefix}_fee"]

    split_connection("grid_import", origin_wacc, breakdown["grid"], "grid")
    split_connection("destination_supply", assumptions["finance"]["destination_wacc"],
                     breakdown["destination_electricity"], "destination_grid")

    # -- what is already as fine as the model cuts it
    for group in ("transmission", "iron_store", "steel_store", "transport"):
        leaves[group] = breakdown[group]

    # Every group is either split into leaves or carried into one, so the leaves
    # come back to the total annual cost. A group that grows a component with no
    # leaf to put it in would otherwise go missing from the stack while every
    # share on the chart still read 100 % of LCOS.
    total = sum(breakdown.values())
    leaf_total = sum(leaves.values())
    if abs(leaf_total - total) > max(1.0, 1e-9 * abs(total)):
        raise ValueError(
            f"the cost leaves come to {leaf_total:,.0f} EUR/yr against a total annual "
            f"cost of {total:,.0f} EUR/yr. Every cost_*_meur group has to be split "
            f"into leaves or carried into one here; LEAF_COSTS in "
            f"common/_report_schema.py is the list of them."
        )

    fields = {f"cost_{leaf}_eur_per_t_steel": value / steel_t for leaf, value in leaves.items()}
    for parent in LEAF_PARENTS:
        fields[f"cost_{parent}_eur_per_t_steel"] = sum(
            value for leaf, value in leaves.items() if LEAF_GROUP[leaf] == parent
        ) / steel_t
    # A few components' leaves added up, beside them rather than among them.
    for name, parts in LEAF_SUBTOTALS.items():
        fields[f"cost_{name}_eur_per_t_steel"] = sum(leaves[leaf] for leaf in parts) / steel_t

    # -- the electricity bill, divided by the job the electricity did. Two systems
    # can pay it: the plant's own, and — on an export route — the market its
    # furnace stands in, which supplies nothing but that furnace. Each system's
    # whole cost is the megawatt-hours it delivered, so pricing every draw at its
    # own system's rate makes the jobs below close on the bill exactly, with no
    # remainder to take up the slack.
    draws = electricity["by_user"]
    destination_mwh = (float(n.generators_t.p["destination_supply"].sum())
                       * annual_scale if "destination_supply" in n.generators.index
                       else 0.0)
    home_cost = sum(breakdown[group] for group in
                    ("res", "grid", "battery", "transmission"))
    home_mwh = electricity["generation_mwh"] - destination_mwh
    home_rate = home_cost / home_mwh if home_mwh > 0 else 0.0
    # The furnace melts on whichever system it stands in; everything else is at
    # home whatever the route. `draws` is keyed by link id as the network spells
    # it — `dri-h2`, not `dri_h2` — which is what REDUCTION_LINKS names too.
    melts_abroad = destination_mwh > 0
    home_draws = {link: mwh for link, mwh in draws.items()
                  if not (melts_abroad and link == "eaf")}
    reduction_mwh = sum(mwh for link, mwh in home_draws.items()
                        if link in REDUCTION_LINKS)
    electrolyser_mwh = home_draws.get("electrolyser", 0.0)
    melt_mwh = 0.0 if melts_abroad else home_draws.get("eaf", 0.0)
    # The press, by name rather than as whatever is left: a link that fits none
    # of the jobs then fails the check below instead of hiding in this one.
    briquetting_mwh = home_draws.get("briquetting", 0.0)
    # And the electricity nobody drew at all: the battery's round trip and the
    # lines'. Its own job rather than part of another, as no step drew it.
    losses_mwh = electricity["losses_mwh"]
    # Every megawatt-hour the home system generated is in exactly one of those
    # five, or the bands below divide the bill by shares that do not add up to
    # it. Checked rather than left as a remainder: a remainder closes on the
    # total whatever it has absorbed, so a step missing from the five would read
    # as part of the leftovers instead of as an error.
    attributed = (reduction_mwh + electrolyser_mwh + melt_mwh
                  + briquetting_mwh + losses_mwh)
    if home_mwh > 0 and abs(attributed - home_mwh) > max(1.0, 1e-6 * home_mwh):
        raise ValueError(
            f"the electricity jobs account for {attributed:,.0f} MWh of the "
            f"{home_mwh:,.0f} MWh this plant generated. Every drawing link has to be "
            f"in REDUCTION_LINKS or be the furnace, the electrolyser or the "
            f"briquetting press; a new one needs a job of its own in "
            f"ELECTRICITY_JOBS — and `draws` is keyed the way the network spells a "
            f"link id (common/_report_schema.py names the lists)."
        )
    jobs = {
        "reduction": reduction_mwh * home_rate,
        "eaf": breakdown["destination_electricity"] if melts_abroad else melt_mwh * home_rate,
        "electrolyser": electrolyser_mwh * home_rate,
        "briquetting": briquetting_mwh * home_rate,
        "losses": losses_mwh * home_rate,
    }
    # The same bill the `electricity` leaves price, divided by what it was for
    # rather than by what was bought, so the two cuts have to land on the same
    # number. A job that swallowed another's megawatt-hours would still close on
    # the bill, which is what the megawatt-hour check above is separately for.
    electricity_bill = home_cost + breakdown["destination_electricity"]
    job_total = sum(jobs.values())
    if abs(job_total - electricity_bill) > max(1.0, 1e-9 * abs(electricity_bill)):
        raise ValueError(
            f"the electricity jobs come to {job_total:,.0f} EUR/yr against an "
            f"electricity bill of {electricity_bill:,.0f} EUR/yr. The jobs divide the "
            f"bill by the megawatt-hours each drew, so they only close while every "
            f"drawing link is in REDUCTION_LINKS or is the furnace, the "
            f"electrolyser or the briquetting press "
            f"(common/_report_schema.py names the lists)."
        )
    fields.update({f"cost_electricity_{job}_eur_per_t_steel": value / steel_t
                   for job, value in jobs.items()})

    # -- ironmaking's share of the origin's electricity, per tonne of iron. Every
    # MWh at the system's average cost, as the jobs above are priced, so the share
    # is the ironmaking steps' draw over everything drawn at the origin — the
    # losses go with it in proportion. On an export route the furnace melts
    # abroad, draws nothing here, and the share is the whole origin system.
    iron_mwh = reduction_mwh + electrolyser_mwh + briquetting_mwh
    iron_share = iron_mwh / (iron_mwh + melt_mwh) if iron_mwh + melt_mwh > 0 else 0.0
    iron_t = _iron_produced_t(n) * annual_scale
    if iron_t > 0:
        per_t_iron = {leaf: leaves[leaf] * iron_share / iron_t for leaf in ORIGIN_ELECTRICITY_LEAVES}
        fields["cost_electricity_eur_per_t_iron"] = sum(per_t_iron.values())
        fields.update({f"cost_electricity_{leaf}_eur_per_t_iron": value
                       for leaf, value in per_t_iron.items()})
        fields.update({f"cost_electricity_{name}_eur_per_t_iron": sum(per_t_iron[leaf] for leaf in parts)
                       for name, parts in LEAF_SUBTOTALS.items()
                       if all(leaf in per_t_iron for leaf in parts)})
        # And by job: each ironmaking job's MWh at the average rate, the losses
        # at ironmaking's share of them, so the four close on the same total.
        iron_jobs = {
            "reduction": reduction_mwh * home_rate,
            "electrolyser": electrolyser_mwh * home_rate,
            "briquetting": briquetting_mwh * home_rate,
            "losses": losses_mwh * iron_share * home_rate,
        }
        assert set(iron_jobs) == set(IRON_ELECTRICITY_JOBS)
        fields.update({f"cost_electricity_{job}_eur_per_t_iron": value / iron_t
                       for job, value in iron_jobs.items()})
        # The hydrogen's own costs per tonne of iron: all of them, since all the
        # hydrogen goes into ironmaking. The total, which also needs the
        # electricity at the LCOH's rate, is written beside the LCOH.
        electrolyser = {part: leaves[f"electrolyser_{part}"] for part in ("capex", "fom", "water")}
        fields["cost_electrolyser_eur_per_t_iron"] = sum(electrolyser.values()) / iron_t
        fields.update({f"cost_electrolyser_{part}_eur_per_t_iron": value / iron_t
                       for part, value in electrolyser.items()})
        fields["cost_h2_buffer_eur_per_t_iron"] = leaves["h2_buffer"] / iron_t

    # The same capital/upkeep split on the carriers, over each carrier's own
    # denominator so the parts stack to the reported part they divide.
    generation_mwh = electricity["generation_mwh"]
    if generation_mwh > 0:
        for part in ("capex", "fom"):
            fields[f"lcoe_res_{part}_eur_per_mwh"] = sum(
                leaves[f"res_{field_stem(tech)}_{part}"] for tech in RES_TECHS
            ) / generation_mwh
    for part in ("capex", "fom", "water"):
        fields[f"cost_electrolyser_{part}_meur"] = leaves[f"electrolyser_{part}"] / 1e6
    h2_kg = _h2_produced_kg(n) if "electrolyser" in n.links.index else 0.0
    if h2_kg > 0:
        for part in ("capex", "fom", "water"):
            fields[f"lcoh_electrolyser_{part}_eur_per_kg"] = leaves[f"electrolyser_{part}"] / h2_kg

    # -- what a tonne of steel took, measured off the run
    fields["el_eaf_mwh_per_t_steel"] = draws.get("eaf", 0.0) / steel_t
    fields["el_reduction_mwh_per_t_steel"] = reduction_mwh / steel_t
    fields["el_electrolyser_mwh_per_t_steel"] = draws.get("electrolyser", 0.0) / steel_t
    fields["gas_mwh_per_t_steel"] = gas_mwh / steel_t
    if "electrolyser" in n.links.index:
        fields["h2_produced_kg_per_t_steel"] = _h2_produced_kg(n) / steel_t
    if "battery" in n.stores.index and "battery_charger" in n.links.index:
        power_mw = float(n.links.at["battery_charger", "p_nom_opt"])
        energy_mwh = float(n.stores.at["battery", "e_nom_opt"])
        if power_mw > 0:
            fields["battery_duration_hours"] = energy_mwh / power_mw
        if energy_mwh > 0:
            fields["max_battery_c_rate"] = power_mw / energy_mwh

    # -- the inputs behind the split, so a cost can be checked against them
    for link in PROCESS_LINKS:
        fields[f"annuity_factor_{field_stem(link)}"] = annuity_factor(
            wacc_of(link), assumptions[link]["lifetime_years"]
        )
    for link in ORE_LINKS:
        fields[f"ore_quote_{field_stem(link)}_eur_per_t_iron"] = assumptions[link]["ore_eur_per_t"]
    if "eaf" in n.links.index:
        # The furnace's efficiency is tonnes of steel per tonne of iron, so its
        # reciprocal is the charge this route feeds it — off the network, rather
        # than back through the route's own charge table.
        fields["eaf_iron_t_per_t_steel"] = 1.0 / float(n.links.at["eaf", "efficiency"])
    return fields


def _res_energy(n: pypsa.Network, hours_per_year: float) -> dict[str, float]:
    """What each renewable delivered in a year, and what it was switched off for, GWh.

    Summed over every generator of the tech — a candidate site or a panel
    orientation is a generator of its own, named after the tech — so the pair
    describes the tech whatever the siting. Generated is the denominator of
    `lcoe_{tech}_eur_per_mwh_{tech}`; curtailed is what the built capacity
    could have produced on top of it.
    """
    available = n.get_switchable_as_dense("Generator", "p_max_pu") * n.generators.p_nom_opt
    annual_scale = hours_per_year / len(n.snapshots)
    built = n.generators.index[
        n.generators.p_nom_extendable
        & ~n.generators.index.isin(("grid_import", "gas_supply", "destination_supply"))
    ]
    fields = {}
    for tech in RES_TECHS:
        idx = [g for g in built if str(g).startswith(tech)]
        generated = float(n.generators_t.p[idx].sum().sum()) * annual_scale if idx else 0.0
        potential = float(available[idx].sum().sum()) * annual_scale if idx else 0.0
        fields[f"{field_stem(tech)}_generated_gwh"] = generated / 1e3
        # Clipped at zero: the solver's tolerance can put dispatch a hair above
        # what was available, which is not negative curtailment.
        fields[f"{field_stem(tech)}_curtailed_gwh"] = max(potential - generated, 0.0) / 1e3
    return fields


def _iron_produced_t(n: pypsa.Network) -> float:
    """Tonnes of iron the ironmaking steps put out over the snapshots.

    The route's own iron product: what leaves the shaft or the cells, or on a
    route that briquettes it, the briquettes — the press is part of making the
    iron, and HBI is what it hands on. Not annualised.
    """
    if "briquetting" in n.links.index:
        return float(-n.links_t.p1["briquetting"].sum())
    return float(sum(-n.links_t.p1[link].sum() for link in ORE_LINKS if link in n.links.index))


def _h2_buffer_discharged_mwh(n: pypsa.Network, hours_per_year: float) -> float:
    """Hydrogen the buffer gave back to the hydrogen bus in a year, MWh LHV.

    A Store's p is positive when it feeds its bus, so the positive hours are its
    discharge. Lossless and cyclic, so the same again went in.
    """
    discharge = n.stores_t.p["h2_buffer"].clip(lower=0.0)
    return float(discharge.sum()) * hours_per_year / len(n.snapshots)


def _battery_discharged_mwh(n: pypsa.Network, hours_per_year: float) -> float:
    """What the battery gave back to the plant in a year, MWh, after the discharge loss.

    The discharger's output side, so a MWh lost on the way out is not counted as
    delivered. The denominator of `lcoe_storage_eur_per_mwh_discharged`.
    """
    drawn = float(n.links_t.p0["battery_discharger"].sum())
    return drawn * float(n.links.at["battery_discharger", "efficiency"]) * hours_per_year / len(n.snapshots)


def extract_summary(
    n: pypsa.Network, scenario_name: str, run: dict, assumptions: dict,
    transport_legs: dict | None = None, grid_mix: pd.DataFrame | None = None,
    destination_mix: pd.DataFrame | None = None,
) -> dict:
    """Key sizing, cost and emission metrics as a flat dict (one row of the CSV).

    `run` identifies the row — origin, destination, route, date range and the input variants —
    and leads the columns, so what a number describes reads before the number.

    The headline is the levelised cost of steel (LCOS) over the flat steel
    load; hydrogen's output and cost are reported wherever an electrolyser runs.

    `transport_legs` and `grid_mix` are the run's own freight legs and the
    generation mix behind its imports; both only feed the emission fields, and
    both may be absent — a route that ships nothing has no legs, and a grid
    series solved on `dayahead` has no mix. `destination_mix` is the same thing
    for the market an `-export` route melts in, and is absent for every route
    that melts its iron at home.
    """
    breakdown = _cost_breakdown(n)
    hours_per_year = _hours_per_year(n.snapshots)
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

    steel_t_per_h = float(n.loads.at["steel_load", "p_set"])
    steel_t_per_year = steel_t_per_h * hours_per_year
    summary["lcos_eur_per_t"] = total_annual_cost / steel_t_per_year
    summary["steel_produced_mt"] = steel_t_per_year / 1e6

    if "electrolyser" in n.links.index:
        summary["h2_produced_kt"] = _h2_produced_kg(n) / 1e6

    # Levelised cost of the underlying energy carriers, €/MWh, so LCOS can be read
    # against the electricity and hydrogen that drive it.
    #   LCOE = electricity-system cost (renewables + battery + grid + transmission
    #          + the destination furnace's own power on an export route) per MWh of
    #          electricity generated (renewable dispatch + grid import + that same
    #          destination supply). Route-wide: every megawatt-hour the route draws
    #          is in the denominator, so every one of them has to be paid for in the
    #          numerator: `elec_gens` picks up the destination supply, so `elec_cost`
    #          carries `destination_electricity` alongside it. Numerator and denominator
    #          cover the same generators, or an export route reads cheaper than its
    #          domestic twin by the MWh nobody was charged for.
    #   LCOH = (electrolyser capex/opex + H2 buffer + the electrolyser's electricity
    #          valued at LCOE) per MWh of H2 produced, LHV.
    annual_scale = hours_per_year / len(n.snapshots)
    elec_gens = [g for g in n.generators.index if g != "gas_supply"]
    elec_mwh = (float(n.generators_t.p[elec_gens].sum().sum()) * annual_scale) if elec_gens else 0.0
    elec_cost = sum(breakdown[k]
                    for k in ("res", "battery", "grid", "transmission", "destination_electricity"))
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
            # Zero on a domestic route, which builds no destination supply. It is a
            # part like the others so the decomposition still closes on LCOE.
            "lcoe_destination_electricity": breakdown["destination_electricity"],
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
                summary[f"lcoe_{field_stem(tech)}_eur_per_mwh_{field_stem(tech)}"] = cost / mwh

        # Blended renewable own LCOE (over renewable generation only) — the
        # generation-weighted mean of the per-tech own LCOEs, so it sits between
        # them. Distinct from lcoe_renewables (contribution over all electricity,
        # which grid imports drag below the per-tech figures).
        res_mwh = sum(_res_mwh(tech) for tech in RES_TECHS)
        if res_total > 0 and res_mwh > 0:
            summary["lcoe_renewables_eur_per_mwh_renewables"] = res_total / res_mwh

        # Capacity factors (annual generation / nameplate × hours per year) for the renewables
        # and the grid connection, surfaced next to the built capacities.
        def _res_cap(prefix: str) -> float:
            idx = [g for g in res_idx if str(g).startswith(prefix)]
            return float(n.generators.loc[idx, "p_nom_opt"].sum()) if idx else 0.0

        for tech in RES_TECHS:
            cap = _res_cap(tech)
            if cap > 0:
                summary[f"cf_{field_stem(tech)}"] = _res_mwh(tech) / (cap * hours_per_year)
        if grid_mwh > 0:
            grid_p_nom = float(n.generators.at["grid_import", "p_nom_opt"])
            if grid_p_nom > 0:
                summary["cf_grid_connection"] = grid_mwh / (grid_p_nom * hours_per_year)

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

    # What the destination furnace paid for its power, market price only and on
    # the same footing as `grid_price_eur_per_mwh`: its bill over the MWh it
    # drew, less the volumetric fee. Draw-weighted rather than the year's mean,
    # because the furnace picks its hours — which is the whole reason its series
    # is hourly. Not floored at zero: a furnace that ran in the hours the price
    # was negative really was paid to melt, and rounding that up to nothing
    # would report a bill it never had.
    if "destination_supply" in n.generators.index:
        destination_mwh = float(n.generators_t.p["destination_supply"].sum())
        if destination_mwh > 0:
            destination_energy = _marginal_costs(
                n.generators.loc[["destination_supply"]],
                n.generators_t.p, n.generators_t.marginal_cost,
            )
            summary["destination_price_eur_per_mwh"] = (
                destination_energy / destination_mwh
                - float(assumptions["grid"]["fee_eur_per_mwh"])
            )

    if "electrolyser" in n.links.index:
        el_mwh = float(n.links_t.p0["electrolyser"].sum()) * annual_scale
        h2_mwh = _h2_produced_kg(n) * H2_LHV_KWH_PER_KG / 1000.0
        if h2_mwh > 0:
            elec_for_h2 = el_mwh * (lcoe if lcoe == lcoe else 0.0)
            el_cost = breakdown["electrolyser"]
            buf_cost = breakdown["h2_buffer"]
            lcoh = (el_cost + buf_cost + elec_for_h2) / h2_mwh
            summary["lcoh_eur_per_mwh_lhv"] = lcoh
            # And per kg, which is the unit the cost-breakdown page charts
            # hydrogen in.
            summary["lcoh_eur_per_kg"] = lcoh * H2_LHV_KWH_PER_KG / 1000.0
            # LCOH decomposition, €/kg over the same H2 denominator so the parts
            # sum back to `lcoh_eur_per_kg`: the electrolyser plant, its
            # electricity (valued at LCOE) and — where built — the H2 buffer store.
            h2_kg = h2_mwh * 1000.0 / H2_LHV_KWH_PER_KG
            summary["lcoh_electrolyser_eur_per_kg"] = el_cost / h2_kg
            summary["lcoh_electricity_eur_per_kg"] = elec_for_h2 / h2_kg
            if buf_cost > 0:
                summary["lcoh_h2_storage_eur_per_kg"] = buf_cost / h2_kg
            # The same three as yearly M€. Not cost groups: the electricity is
            # part of the electricity groups already, so these add to each
            # other and to `cost_electrolyser_meur` / `cost_h2_buffer_meur` only.
            summary["cost_electricity_electrolyser_meur"] = elec_for_h2 / 1e6
            summary["cost_h2_production_meur"] = (el_cost + buf_cost + elec_for_h2) / 1e6
            summary["cost_h2_production_eur_per_t_steel"] = (
                (el_cost + buf_cost + elec_for_h2) / steel_t_per_year)
            iron_t = _iron_produced_t(n) * annual_scale
            if iron_t > 0:
                summary["cost_h2_production_eur_per_t_iron"] = (
                    (el_cost + buf_cost + elec_for_h2) / iron_t)

    # Cost detail for the hover: per-plant levelised capex and the ore-vs-consumables
    # split, all €/t steel on the same denominator as the LCOS cost groups — so the
    # per-plant figures sum to the "process" group and ore+consumables to the
    # "ore_consumables" group. Ore is the feedstock priced on the reduction/
    # electrolysis links; consumables is the EAF's marginal cost.
    def _link_capex(link_name: str) -> float:
        return (float(n.links.at[link_name, "capital_cost"] * n.links.at[link_name, "p_nom_opt"])
                if link_name in n.links.index and bool(n.links.at[link_name, "p_nom_extendable"])
                else 0.0)

    def _link_marginal(link_name: str) -> float:
        if link_name not in n.links.index or link_name not in n.links_t.p0.columns:
            return 0.0
        return float((n.links_t.p0[link_name] * n.links.at[link_name, "marginal_cost"]).sum()) * annual_scale

    for link in PROCESS_LINKS:
        summary[f"plant_{field_stem(link)}_eur_per_t_steel"] = _link_capex(link) / steel_t_per_year
    # Every link that buys ore — the blended shaft included, where the
    # ore is over a third of the cost of the steel — so this column and
    # `consumables` add up to the `ore_consumables` group and the leaf
    # split under it reaches the same total.
    ore = sum(_link_marginal(link) for link in ORE_LINKS)
    summary["ore_eur_per_t_steel"] = ore / steel_t_per_year
    summary["consumables_eur_per_t_steel"] = _link_marginal("eaf") / steel_t_per_year

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
        gas_mwh = float(n.generators_t.p["gas_supply"].sum()) * (hours_per_year / len(n.snapshots))
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
        summary["reductant_h2_share"] = from_h2 / total if total else float("nan")
        # The same split on what went in rather than what it reduced: each link's
        # input side is the fuel itself, MWh LHV of hydrogen or of gas.
        h2_fed = float(n.links_t.p0[h2_link].sum()) if h2_link in n.links.index else 0.0
        ng_fed = float(n.links_t.p0[ng_link].sum()) if ng_link in n.links.index else 0.0
        fed = h2_fed + ng_fed
        summary["reductant_h2_share_energy"] = h2_fed / fed if fed > 0 else float("nan")

    if "iron_store" in n.stores.index:
        store_t = n.stores.at["iron_store", "e_nom_opt"]
        summary["iron_store_kt"] = store_t / 1e3
        summary["iron_store_hours_steel"] = store_t / steel_t_per_h

    annual = hours_per_year / len(n.snapshots)
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
        summary["steel_store_hours_steel"] = store_t / steel_t_per_h

    for gen in n.generators.index[n.generators.p_nom_extendable]:
        summary[f"{field_stem(gen)}_gw_opt"] = n.generators.at[gen, "p_nom_opt"] / 1e3

    summary.update(_res_energy(n, hours_per_year))
    if "grid_import" in n.generators.index:
        summary["grid_imported_gwh"] = (float(n.generators_t.p["grid_import"].sum())
                                        * hours_per_year / len(n.snapshots) / 1e3)

    if "battery" in n.stores.index:
        # Power is the shared inverter rating, energy is the store — now two
        # separate answers from the solve rather than one and a fixed ratio.
        summary["battery_gw_opt"] = n.links.at["battery_charger", "p_nom_opt"] / 1e3
        summary["battery_mwh_opt"] = n.stores.at["battery", "e_nom_opt"]
        discharged_mwh = _battery_discharged_mwh(n, hours_per_year)
        summary["battery_discharged_gwh"] = discharged_mwh / 1e3
        if summary["battery_mwh_opt"] > 0:
            summary["battery_cycles_per_year"] = discharged_mwh / summary["battery_mwh_opt"]
        if discharged_mwh > 0 and breakdown["battery"] > 0:
            summary["lcoe_storage_eur_per_mwh_discharged"] = breakdown["battery"] / discharged_mwh

    if "h2_buffer" in n.stores.index:
        buffer_mwh = n.stores.at["h2_buffer", "e_nom_opt"]
        summary["h2_buffer_gwh"] = buffer_mwh / 1e3
        summary["h2_buffer_kt"] = buffer_mwh * 1000.0 / H2_LHV_KWH_PER_KG / 1e6
        # Additionally as hours of the shaft's mean H₂ draw: what the H2 shaft
        # takes, or on the blend shaft what its hydrogen feed takes.
        h2_demand_mw = sum(float(n.links_t.p0[link].mean())
                           for link in ("dri-h2", "reductant-h2") if link in n.links.index)
        summary["h2_buffer_duration_hours"] = (
            buffer_mwh / h2_demand_mw if buffer_mwh > 0 and h2_demand_mw > 0 else float("nan")
        )
        discharged_mwh = _h2_buffer_discharged_mwh(n, hours_per_year)
        discharged_kg = discharged_mwh * 1000.0 / H2_LHV_KWH_PER_KG
        summary["h2_buffer_discharged_kt"] = discharged_kg / 1e6
        if buffer_mwh > 0:
            summary["h2_buffer_cycles_per_year"] = discharged_mwh / buffer_mwh
        if "electrolyser" in n.links.index and _h2_produced_kg(n) > 0:
            summary["h2_buffer_share_of_h2_pct"] = discharged_kg * 100.0 / _h2_produced_kg(n)

    if "electrolyser" in n.links.index:
        el_cap = n.links.at["electrolyser", "p_nom_opt"]
        summary["electrolyser_gw"] = el_cap / 1e3
        if el_cap > 0 and "electrolyser" in n.links_t.p0.columns:
            summary["electrolyser_utilization"] = float(
                n.links_t.p0["electrolyser"].mean() / el_cap
            )
        else:
            summary["electrolyser_utilization"] = float("nan")

    # Multi-site only: one column per HVDC link capacity, and per-tech built
    # capacity summed across candidate sites. What the links cost is
    # `cost_transmission_meur`, with the other cost groups.
    for link in n.links.index[n.links.carrier == "HVDC"]:
        summary[f"{field_stem(link)}_gw_opt"] = n.links.at[link, "p_nom_opt"] / 1e3

    # A run offered several candidate sites names a generator per site, so no
    # single `{tech}_gw_opt` says what that tech built and the total is reported
    # beside them. Only the sited techs: a grid connection and a gas supply are
    # one generator whatever the siting, and on an export route the destination's
    # supply carries the same `AC` carrier as the home connection — summing those
    # two adds a connection in one country to a connection in another.
    extendable = n.generators[n.generators.p_nom_extendable]
    sited = extendable[extendable["carrier"].isin(RES_TECHS)]
    multi_site = len(sited) > sited["carrier"].nunique()
    if multi_site:
        for carrier, at_sites in sited.groupby("carrier"):
            summary[f"{field_stem(carrier)}_total_gw_opt"] = (
                at_sites["p_nom_opt"].sum() / 1e3
            )

    # Emissions, on the basis the assumptions name. Reported in kg so the
    # report's two decimals still say something about a clean route.
    emitted = _emissions_breakdown(
        n, assumptions["emissions"],
        run["origin"], transport_legs, grid_mix,
        run["destination"], destination_mix,
    )
    emissions = emitted["by_step"]
    electricity_mwh = emitted["electricity_mwh"]
    sources = emitted["sources"]
    summary["emissions_basis"] = assumptions["emissions"]["basis"]
    if emitted["unavailable"]:
        summary["emissions_unavailable_reason"] = emitted["unavailable"]
    emissions_t = sum(emissions.values())
    freight_by_leg = emitted["freight_by_leg"]
    summary["emissions_kt_co2e_per_year"] = emissions_t / 1e3
    for source, value in sources.items():
        summary[f"emissions_{source}_kt_co2e_per_year"] = value / 1e3
    summary["emissions_freight_kt_co2e_per_year"] = sum(freight_by_leg.values()) / 1e3

    # Per tonne of steel, and each step's share of it. The steps stack to the
    # total, so a step this route has none of reads 0 rather than blank. The
    # freight sits beside them rather than in them: it is what delivering the
    # tonne added, and no total above counts it.
    summary["emissions_kg_co2e_per_t_steel"] = emissions_t * 1e3 / steel_t_per_year
    for step in EMISSION_STEPS:
        summary[f"emissions_{field_stem(step)}_kg_co2e_per_t_steel"] = (
            emissions.get(step, 0.0) * 1e3 / steel_t_per_year
        )
    for leg in FREIGHT_LEGS:
        summary[f"emissions_{field_stem(leg)}_kg_co2e_per_t_steel"] = (
            freight_by_leg.get(leg, 0.0) * 1e3 / steel_t_per_year
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
    #
    # The losses ride along in proportion to what the electrolyser drew: the
    # battery cycle and the line that delivered its power were incurred to
    # deliver it, and a run whose every MWh exists to make hydrogen would
    # otherwise read a third light against a threshold it is being held to.
    drawn_t = sum(emitted["electricity_t"].values())
    if "electrolyser" in n.links.index:
        h2_kg = _h2_produced_kg(n)
        if h2_kg > 0 and drawn_t > 0:
            summary["emissions_kg_co2e_per_kg_h2"] = (
                emitted["electricity_t"].get("electrolyser", 0.0)
                * (sources["electricity"] / drawn_t) * 1e3 / h2_kg
            )
        elif h2_kg > 0 and drawn_t == 0:
            # Every MWh came from carriers this basis rates at zero, so the
            # hydrogen carries nothing. That is the figure the RFNBO ceiling is
            # read against, so it is reported rather than left blank. A run whose
            # mix is unknown has a NaN here instead and stays blank.
            summary["emissions_kg_co2e_per_kg_h2"] = 0.0

    # Electricity by who drew it, and how dirty their own hours were — the number
    # that separates a user chasing cheap renewable hours from one running flat.
    # The losses are on both sides of the system average: they emitted, and they
    # were MWh the system drew, so a user can sit either side of it.
    total_el_mwh = sum(electricity_mwh.values()) + emitted["losses_mwh"]
    if total_el_mwh > 0:
        summary["emissions_kg_co2e_per_mwh_el"] = (
            sources["electricity"] * 1e3 / total_el_mwh
        )
    summary["el_users_gwh"] = sum(electricity_mwh.values()) / 1e3
    for user, mwh in electricity_mwh.items():
        summary[f"el_{field_stem(user)}_gwh"] = mwh / 1e3
        if mwh > 0:
            summary[f"emissions_{field_stem(user)}_kg_co2e_per_mwh_el"] = (
                emitted["electricity_t"][user] * 1e3 / mwh
            )
    # And what nobody drew, so the draws above account for the whole of what was
    # generated rather than most of it.
    summary["el_losses_gwh"] = emitted["losses_mwh"] / 1e3

    # The cost groups cut as finely as the model allows — the two taxonomies the
    # cost-breakdown page plots. Last, because it reads each user's electricity
    # draw off the emissions breakdown above rather than reading the ports a
    # second time.
    summary.update(_leaf_breakdown(
        n, assumptions, breakdown, steel_t_per_year,
        {"by_user": electricity_mwh, "losses_mwh": emitted["losses_mwh"],
         "generation_mwh": elec_mwh},
    ))

    return summary


def input_variants(rows: pd.DataFrame, origin: str) -> dict:
    """Which series each tech contributed to one run, as its `variant_field`s.

    A best-site P95 profile and an area average answer different questions, so a
    cost is not interpretable without knowing which produced it. The table spells
    a tech `wind-onshore`; the field it lands in is `wind_onshore_variant`.

    The grid's is the origin's own connection, so a blank `grid_variant_origin`
    means the plant was islanded. An export run's destination has a grid row too, for
    the furnace that buys its power there; it is always `emissions`, and the
    `destination` field already names the market, so it gets no field of its own
    and must not overwrite the origin's.
    """
    own = rows[(rows["tech"] != "grid") | (rows["area"] == origin)]
    return {variant_field(row.tech): row.variant for row in own.itertuples()}


def main() -> None:
    """Summarise one solved run into its report and the hidden diagnostic beside it.

    The run's rows in its scenario table name its origin, destination and route,
    and the variant each tech was solved with; the netCDF name carries none of
    them. One summary row via `extract_summary`, then the two files.
    """
    scenario_name = snakemake.wildcards.scenario
    run_name = snakemake.wildcards.run
    run_inputs = read_run_inputs(Path(snakemake.input.scenario_table), scenario_name,
                                 run_name, list(snakemake.input.grid_input))
    rows = run_inputs["rows"]
    # The same merge the solve did, so the cost leaves are split on the quotes
    # that priced the run rather than on the base file's.
    base_text = Path(snakemake.input.assumptions_base).read_text()
    overlays = list(snakemake.input.assumptions_overlay)
    overlay_text = Path(overlays[0]).read_text() if overlays else None
    assumptions = deep_merge(yaml.safe_load(base_text), yaml.safe_load(overlay_text or "") or {})
    # The run's own copy of what priced it: the base, the overlay and the merge
    # of them, so a reader does not have to merge the two in their head. The
    # overlay is written even when the scenario has none, because "this changed
    # nothing" is a fact worth recording and an absent file cannot say it.
    Path(snakemake.output.base_assumptions).parent.mkdir(parents=True, exist_ok=True)
    Path(snakemake.output.base_assumptions).write_text(base_text)
    Path(snakemake.output.overlay).write_text(
        overlay_text or f"# {scenario_name} has no overlay: it ran on the base assumptions unchanged.\n"
    )
    Path(snakemake.output.merged_assumptions).write_text(
        "# The base with this scenario's overlay merged in: what actually priced the\n"
        "# run, rather than the two files a reader would have to merge in their head.\n"
        + yaml.safe_dump(assumptions, sort_keys=False)
    )

    run = {
        "run": run_name,
        "origin": run_inputs["origin"],
        "destination": run_inputs["destination"],
        "route": run_inputs["route"],
        "start_date": rows["start_date"].iloc[0],
        "end_date": rows["end_date"].iloc[0],
        **input_variants(rows, run_inputs["origin"]),
    }
    n = pypsa.Network()
    n.import_from_netcdf(Path(snakemake.input.network))
    read_mix = lambda path: pd.read_parquet(path) if path is not None else None
    summary = extract_summary(
        n, scenario_name, run, assumptions,
        find_transport_legs(assumptions["transport"]["distance_km"], run_inputs["origin"],
                            run_inputs["destination"], snakemake.config["areas"]),
        read_mix(run_inputs["grid_path"]), read_mix(run_inputs["destination_path"]),
    )
    # Trailing the row, because it identifies the inputs rather than describing
    # them: the per-file map it stands for is in the network.
    summary["inputs_hash"] = n.meta["inputs_hash"]
    summary["code_hash"] = n.meta.get("code_hash", "")

    # The report is the seam viz reads; the diagnostic keeps every field,
    # including the ones the report is not the place for.
    report = pd.DataFrame([summary])
    write_report_file(report, Path(snakemake.output.diagnostic))
    write_report_file(report, Path(snakemake.output.report), hold_back=DIAGNOSTIC_FIELDS)
    log.info(f"wrote {snakemake.output.report} and {snakemake.output.diagnostic}")


if __name__ == "__main__":
    main()
