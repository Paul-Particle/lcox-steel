"""Unit tests for the report's two fine cuts of the cost of steel.

`cost_*_eur_per_t_steel` (by what was bought) and `cost_electricity_*_eur_per_t_steel` (by what it was
for) each have to stack back to LCOS. The cost-breakdown page prints the
*reported* total on top of the stack it draws, so anything the split does not
reach goes missing from the bars while every share still reads 100 %.

Built on real networks from `build_network`, so the capital and marginal costs
are the ones that price a solve, with a solved state written in by hand rather
than optimised: closure is arithmetic and holds whatever the optimiser would
have picked, so the tests need no solver and no cutouts.
"""

import copy

import numpy as np
import pandas as pd
import pytest
import yaml

from common._constants import H2_LHV_KWH_PER_KG
from common._workbook_translation import VALUES_OF_KEYS_NOT_REQUIRED_BY_THE_WORKBOOK_SCENARIO_CONFIGURATION
from common.generate_inputs import set_nested
from scripts.viz import compile_report
from scripts.solve._helpers_solve import annuity_factor
from scripts.solve.build_network import build_network
from common._report_schema import (
    ELECTRICITY_USERS,
    LEAF_COSTS,
    LEAF_GROUP,
    LEAF_PARENTS,
    LEAF_SUBTOTALS,
    ORE_LINKS,
    plant_stem,
    IRON_ELECTRICITY_JOBS,
    ORIGIN_ELECTRICITY_LEAVES,
    PROCESS_LINKS,
    ELECTRICITY_JOBS,
    field_stem,
)
from conftest import REPO_ROOT

HOURS = 24
# Every route that makes steel and can be built without a multi-site overlay.
# The blended shaft and the export twins are the ones the split has to reach.
ROUTES = ["h2-dri-eaf", "mix-dri-eaf", "ng-dri-eaf", "moe-eaf", "ew-eaf",
          "h2-dri-eaf-export", "moe-eaf-export", "ew-eaf-export",
          "mix-dri-eaf-export", "ng-dri-eaf-export"]


@pytest.fixture(scope="module")
def assumptions() -> dict:
    """The real assumptions, so the capital costs are composed as a solve's are."""
    return yaml.safe_load((REPO_ROOT / "config" / ".template.yaml").read_text())


def _solved(route: str, assumptions: dict):
    """One route's network with a solved state written in, no optimiser involved.

    Capacities and flows are arbitrary but non-zero on every component, which is
    the point: each leaf has to pick up its own share of whatever was built.
    """
    index = pd.date_range("2025-01-01", periods=HOURS, freq="h")
    profile = np.clip(np.sin(np.arange(HOURS) / 12 * np.pi), 0.05, None)
    cf = pd.DataFrame({"solar": profile, "wind-onshore": profile * 0.8}, index=index)
    legs = {"sea": 8945} if route.endswith("-export") else None
    price = pd.Series(np.linspace(20.0, 60.0, HOURS), index=index)
    n = build_network(
        route, assumptions, cf, price_series=price, transport_legs=legs,
        destination_price=price * 1.1 if route.endswith("-export") else None,
    )
    for component, capacity, flow in (("generators", "p_nom", "p"),
                                      ("links", "p_nom", "p0"),
                                      ("stores", "e_nom", "e")):
        static = getattr(n, component)
        extendable = static.index[static[f"{capacity}_extendable"]]
        static.loc[extendable, f"{capacity}_opt"] = 100.0
        frame = pd.DataFrame(50.0, index=n.snapshots, columns=static.index)
        getattr(n, component + "_t")[flow] = frame
    # A link's output side, which the iron-from-H2 share and the hydrogen read.
    n.links_t.p1 = pd.DataFrame(-40.0, index=n.snapshots, columns=n.links.index)
    return n


def _fields(n, assumptions: dict) -> tuple[dict, dict, float]:
    """(leaf fields, cost groups, tonnes of steel) for one written-in network."""
    breakdown = compile_report._cost_breakdown(n)
    steel_t = float(n.loads.at["steel_load", "p_set"]) * 8760.0
    # The drawing links and the generation, as extract_summary hands them over.
    annual = 8760.0 / len(n.snapshots)
    # Keyed by link id, exactly as `_emissions_breakdown` hands them over: the
    # network spells a shaft `dri-h2`, and a stemmed key here would let a
    # mis-keyed lookup in the split pass the test and fail on a real report.
    drawn = {link: float(n.links_t.p0[link].sum()) * annual
             for link in ELECTRICITY_USERS if link in n.links.index}
    # Generation balances the draws and the losses, because that is the one
    # thing the electricity bands need to be true: they divide the bill by the
    # megawatt-hours each job took, so the jobs have to be the generation. The
    # written-in dispatch puts the same flow on every link, which by itself
    # draws more than the generators make.
    losses_mwh = 0.02 * sum(drawn.values())
    elec_mwh = sum(drawn.values()) + losses_mwh
    fields = compile_report._leaf_breakdown(
        n, assumptions, breakdown, steel_t,
        {"by_user": drawn, "losses_mwh": losses_mwh, "generation_mwh": elec_mwh},
    )
    return fields, breakdown, steel_t


@pytest.mark.parametrize("route", ROUTES)
def test_the_leaves_stack_to_the_cost_of_steel(route, assumptions):
    """Every group is split into leaves or carried into one, so nothing is lost."""
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)

    leaves = sum(fields[f"cost_{leaf}_eur_per_t_steel"] for leaf in LEAF_COSTS)
    assert leaves == pytest.approx(sum(breakdown.values()) / steel_t, rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_the_electricity_jobs_stack_to_the_electricity_bill(route, assumptions):
    """The four jobs divide one bill, with none of them left as a remainder.

    They only close while every drawing link is accounted for — which is what
    `_leaf_breakdown` raises about rather than flooring a residual. Closing on
    the bill rather than on LCOS is the point: the jobs re-cut one cost group,
    so they must land on what that group's own leaves come to.
    """
    n = _solved(route, assumptions)
    fields, _, _ = _fields(n, assumptions)

    jobs = sum(fields[f"cost_electricity_{job}_eur_per_t_steel"] for job in ELECTRICITY_JOBS)
    assert jobs == pytest.approx(fields["cost_electricity_eur_per_t_steel"], rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_each_parent_group_is_its_own_leaves(route, assumptions):
    """The coarse reading beside the leaves is those leaves and nothing else."""
    n = _solved(route, assumptions)
    fields, _, _ = _fields(n, assumptions)

    for parent in LEAF_PARENTS:
        own = sum(fields[f"cost_{leaf}_eur_per_t_steel"] for leaf in LEAF_COSTS
                  if LEAF_GROUP[leaf] == parent)
        assert fields[f"cost_{parent}_eur_per_t_steel"] == pytest.approx(own, rel=1e-9)


@pytest.mark.parametrize("route", ["h2-dri-eaf", "mix-dri-eaf", "moe-eaf",
                                  "h2-dri-eaf-export"])
def test_a_plants_two_halves_are_its_whole_annual_cost(route, assumptions):
    """Capital and fixed O&M divide the plant cost the report already published.

    `capital_cost` in the network is their sum, so this is the one split that
    still needs the quotes — and the one most easily got wrong.
    """
    n = _solved(route, assumptions)
    fields, _, steel_t = _fields(n, assumptions)

    for link in PROCESS_LINKS:
        if link not in n.links.index:
            continue
        # The route's one ironmaking step reports as the ironmaking plant.
        stem = plant_stem(link)
        whole = float(n.links.at[link, "capital_cost"]
                      * n.links.at[link, "p_nom_opt"]) / steel_t
        halves = (fields[f"cost_{stem}_capex_eur_per_t_steel"]
                  + fields[f"cost_{stem}_fom_eur_per_t_steel"])
        assert halves == pytest.approx(whole, rel=1e-9)
        # And neither half is the whole of it: a fixed-opex quote of about an
        # eighth of capex puts the upkeep well inside these bounds.
        assert 0.0 < fields[f"cost_{stem}_fom_eur_per_t_steel"] < whole


@pytest.mark.parametrize("route", ["h2-dri-eaf", "mix-dri-eaf", "moe-eaf", "ew-eaf"])
def test_the_power_that_made_the_iron_is_its_own_band(route, assumptions):
    """The reduction step is named rather than swept in with the leftovers.

    It was a remainder before, under a label that called it the rest of the
    plant — and on the electrolytic routes it is the largest single electricity
    item there is. A remainder also hides a mis-keyed draw: a lookup that missed
    every shaft still closed on the total, with the whole reduction step filed
    as something nothing drew.
    """
    n = _solved(route, assumptions)
    fields, _, _ = _fields(n, assumptions)

    assert fields["cost_electricity_reduction_eur_per_t_steel"] > 0
    assert fields["el_reduction_mwh_per_t_steel"] > 0
    assert (fields["cost_electricity_reduction_eur_per_t_steel"]
            > fields["cost_electricity_briquetting_eur_per_t_steel"] + fields["cost_electricity_losses_eur_per_t_steel"])


@pytest.mark.parametrize("route", ["mix-dri-eaf", "mix-dri-eaf-export"])
def test_a_blended_shaft_pays_for_its_ore(route, assumptions):
    """The ore leaf names every link that buys ore, the blended shaft included.

    `ORE_LINKS` has to name the blended shaft too, or on this route both the
    column and the leaf under it read zero while the run paid for the ore.
    """
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)

    assert fields["cost_ore_eur_per_t_steel"] > 0
    assert fields["cost_ore_eur_per_t_steel"] == pytest.approx(breakdown["ore"] / steel_t, rel=1e-9)
    assert fields["cost_eaf_consumables_eur_per_t_steel"] == pytest.approx(
        breakdown["eaf_consumables"] / steel_t, rel=1e-9)


def test_the_electrolyser_pays_its_water_per_mwh_of_hydrogen_out(assumptions):
    """Quoted per MWh of H2, so a better efficiency does not change the water bill per MWh of H2."""
    electrolyser = _solved("h2-dri-eaf", assumptions).links.loc["electrolyser"]

    assert electrolyser["marginal_cost"] / electrolyser["efficiency"] == pytest.approx(
        assumptions["electrolyser"]["varopex_eur_per_mwh_h2"])


def test_the_hydrogen_buffer_pays_its_fixed_om_on_top_of_its_annuitised_capex(assumptions):
    """Per MWh of storage a year, as every other element's capital cost carries its fixed O&M."""
    buffer = _solved("h2-dri-eaf", assumptions).stores.loc["h2_buffer"]
    buffer_cfg = assumptions["h2_buffer"]
    annuitised_capex = (annuity_factor(assumptions["finance"]["origin_wacc"], buffer_cfg["financing_period"])
                        * buffer_cfg["capex_per_mwh_eur"])

    assert buffer["capital_cost"] == pytest.approx(annuitised_capex + buffer_cfg["opex_per_mwh_per_year_eur"])


def test_the_battery_pays_fixed_om_on_its_energy_and_its_inverter(assumptions):
    """Each capacity carries its own fixed O&M on top of its annuitised capex."""
    n = _solved("h2-dri-eaf", assumptions)
    battery_cfg = assumptions["battery"]
    annuity = annuity_factor(assumptions["finance"]["origin_wacc"], battery_cfg["financing_period"])

    assert n.stores.at["battery", "capital_cost"] == pytest.approx(
        annuity * battery_cfg["capex_per_mwh_eur"] + battery_cfg["opex_per_mwh_per_year_eur"])
    assert n.links.at["battery_charger", "capital_cost"] == pytest.approx(
        annuity * battery_cfg["capex_per_mw_eur"] + battery_cfg["opex_per_mw_per_year_eur"])


@pytest.mark.parametrize("route", ["h2-dri-eaf", "moe-eaf", "ew-eaf",
                                   "h2-dri-eaf-export", "moe-eaf-export", "ew-eaf-export"])
def test_the_workbooks_routes_build_and_cost_without_the_keys_its_scenarios_do_not_require(route, assumptions):
    """NaN where no workbook scenario requires a value stops only the runs that read it, and none of its routes do."""
    without_unrequired = copy.deepcopy(assumptions)
    for path, value in VALUES_OF_KEYS_NOT_REQUIRED_BY_THE_WORKBOOK_SCENARIO_CONFIGURATION.items():
        set_nested(without_unrequired, path, value)
    n = _solved(route, without_unrequired)
    fields, breakdown, steel_t = _fields(n, without_unrequired)

    assert all(np.isfinite(value) for value in fields.values())
    leaves = sum(fields[f"cost_{leaf}_eur_per_t_steel"] for leaf in LEAF_COSTS)
    assert leaves == pytest.approx(sum(breakdown.values()) / steel_t, rel=1e-9)


def summarise_written_in_run_with_generation_matching_the_draws(route: str, assumptions: dict) -> dict:
    """The report row of a written-in run whose home generation is exactly its draws and losses.

    The written-in flows draw more than the generators make, which the report's
    electricity jobs refuse, so the home generators are scaled to what was drawn.
    """
    n = _solved(route, assumptions)
    # What the stores gave back, which the report's stockpile fields read.
    n.stores_t.p = pd.DataFrame(0.0, index=n.snapshots, columns=n.stores.index)
    is_export = route.endswith("-export")
    run = {"run": "test", "origin": "ESP", "destination": "DEU" if is_export else "ESP"}
    legs = {"sea": 8945} if is_export else None
    grid_mix = pd.DataFrame({"gas": 100.0, "solar": 100.0}, index=n.snapshots)
    emitted = compile_report._emissions_breakdown(n, assumptions["emissions"], run["origin"], legs, grid_mix,
                                                  run["destination"], grid_mix)
    home_draws_mwh = sum(mwh for link, mwh in emitted["electricity_mwh"].items()
                         if not (is_export and link == "eaf"))
    home_generators = [g for g in n.generators.index if g not in ("gas_supply", "destination_supply")]
    annual = 8760.0 / len(n.snapshots)
    generated_mwh = float(n.generators_t.p[home_generators].sum().sum()) * annual
    n.generators_t.p[home_generators] *= (home_draws_mwh + emitted["losses_mwh"]) / generated_mwh
    summary = compile_report.extract_summary(n, "test", run, assumptions, legs, grid_mix, grid_mix)
    return summary


def test_hydrogen_is_priced_at_the_origins_electricity_on_an_export_route_too(assumptions):
    """The furnace's power abroad made no hydrogen, so LCOH reads the home rate, as the electrolyser's electricity job does."""
    domestic = summarise_written_in_run_with_generation_matching_the_draws("h2-dri-eaf", assumptions)
    export = summarise_written_in_run_with_generation_matching_the_draws("h2-dri-eaf-export", assumptions)

    for summary in (domestic, export):
        steel_t = summary["steel_produced_mt"] * 1e6
        assert summary["cost_electricity_electrolyser_meur"] * 1e6 == pytest.approx(
            summary["cost_electricity_electrolyser_eur_per_t_steel"] * steel_t, rel=1e-9)
    # Both twins build the same plant at home, so they make hydrogen at the same cost.
    assert export["lcoh_eur_per_mwh_lhv"] == pytest.approx(domestic["lcoh_eur_per_mwh_lhv"], rel=1e-9)
    assert export["lcoe_eur_per_mwh"] != pytest.approx(domestic["lcoe_eur_per_mwh"], rel=1e-6)


@pytest.mark.parametrize("group", ["ladle_metallurgy", "total_plant"])
def test_a_cost_group_with_nowhere_to_go_is_refused(group, assumptions):
    """A cost group no leaf accounts for stops the report.

    Silence here is the failure that matters: the stack would simply draw short
    of the total it prints on top of itself. Both shapes of it are refused — a
    group added to the breakdown and never split, and a group that grew a
    component the split does not reach.

    Two leaves are deliberately residuals of their group, so they absorb an
    inflated one instead: the gas bill less its carbon price, and the grid's
    energy less its connection and its volumetric fee. Nothing can be added to
    either group without the market price or the fuel bill being what it is.
    """
    n = _solved("h2-dri-eaf", assumptions)
    breakdown = compile_report._cost_breakdown(n)
    breakdown[group] = breakdown.get(group, 0.0) + 1e6

    with pytest.raises(ValueError, match="LEAF_COSTS"):
        compile_report._leaf_breakdown(
            n, assumptions, breakdown,
            float(n.loads.at["steel_load", "p_set"]) * 8760.0,
            {"by_user": {"eaf": 1e5}, "losses_mwh": 0.0, "generation_mwh": 1e6},
        )


def test_a_renewable_the_report_has_no_technology_for_is_refused(assumptions):
    """A generator none of RES_TECHS names stops the report rather than being parked.

    Its cost has no capex and O&M split to go into and no row of its own, so the
    only alternatives to stopping are a catch-all row or a stack that draws short.
    """
    n = _solved("h2-dri-eaf", assumptions)
    n.add("Generator", "geothermal", bus="electricity", p_nom_extendable=True,
          capital_cost=1e5)
    n.generators.loc["geothermal", "p_nom_opt"] = 100.0

    with pytest.raises(ValueError, match="geothermal"):
        compile_report._leaf_breakdown(
            n, assumptions, compile_report._cost_breakdown(n),
            float(n.loads.at["steel_load", "p_set"]) * 8760.0,
            {"by_user": {"eaf": 1e5}, "losses_mwh": 0.0, "generation_mwh": 1e6},
        )


@pytest.mark.parametrize("route", [r for r in ROUTES if r.endswith("-export")])
def test_the_destination_s_connection_splits_like_the_origin_s(route, assumptions):
    """Four leaves, as the origin's grid has, and their total reported beside them.

    The total is what the destination's cost group is, so splitting it moved
    nothing; and it is not a leaf, so the stack does not count it twice.
    """
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)
    parts = [f"cost_destination_grid_{part}_eur_per_t_steel"
             for part in ("connection_capex", "capacity_fee", "market", "fee")]

    assert all(fields[part] >= 0 for part in parts)
    assert fields["cost_destination_grid_connection_capex_eur_per_t_steel"] > 0
    assert sum(fields[part] for part in parts) == pytest.approx(
        fields["cost_destination_electricity_eur_per_t_steel"], rel=1e-9)
    assert fields["cost_destination_electricity_eur_per_t_steel"] == pytest.approx(
        breakdown["destination_electricity"] / steel_t, rel=1e-9)
    assert "destination_electricity" not in LEAF_COSTS



@pytest.mark.parametrize("route", ROUTES)
def test_a_subtotal_is_its_leaves_and_its_cost_group(route, assumptions):
    """Each subtotal adds up its own leaves, and where a cost group covers exactly
    those leaves — renewables, battery, grid, destination — it is that group."""
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)
    for name, parts in LEAF_SUBTOTALS.items():
        assert fields[f"cost_{name}_eur_per_t_steel"] == pytest.approx(
            sum(fields[f"cost_{leaf}_eur_per_t_steel"] for leaf in parts), rel=1e-9)
    for name in ("res", "battery", "grid", "destination_electricity"):
        assert fields[f"cost_{name}_eur_per_t_steel"] == pytest.approx(
            breakdown[name] / steel_t, rel=1e-9, abs=1e-9)


def test_a_drawing_link_with_no_job_of_its_own_is_refused(assumptions):
    """Briquetting is the press and only the press, so a new user cannot hide in it.

    Its megawatt-hours go to no job, the jobs then fall short of what the plant
    generated, and the report stops rather than pricing the bill on the rest.
    """
    n = _solved("h2-dri-eaf", assumptions)
    with pytest.raises(ValueError, match="ELECTRICITY_JOBS"):
        compile_report._leaf_breakdown(
            n, assumptions, compile_report._cost_breakdown(n),
            float(n.loads.at["steel_load", "p_set"]) * 8760.0,
            {"by_user": {"eaf": 1e5, "compressor": 5e4}, "losses_mwh": 0.0,
             "generation_mwh": 1.5e5},
        )



@pytest.mark.parametrize("route", ROUTES)
def test_the_iron_s_electricity_is_its_share_of_the_steel_s(route, assumptions):
    """The iron's electricity rows stack to their total, each subtotal is its parts,
    and the iron's share plus the furnace's comes back to the steel's whole bill.

    Every MWh is at the system's average cost, so the iron carries its own draw
    and the losses in proportion to it, and the furnace the rest — or, melting
    abroad, its own destination bill and nothing at the origin.
    """
    n = _solved(route, assumptions)
    fields, _, steel_t = _fields(n, assumptions)
    iron_t = compile_report._iron_produced_t(n) * 8760.0 / len(n.snapshots)
    total = fields["cost_electricity_eur_per_t_iron"]

    assert total > 0
    assert sum(fields[f"cost_electricity_{leaf}_eur_per_t_iron"]
               for leaf in ORIGIN_ELECTRICITY_LEAVES) == pytest.approx(total, rel=1e-9)
    for name, parts in LEAF_SUBTOTALS.items():
        if set(parts) <= set(ORIGIN_ELECTRICITY_LEAVES):
            assert fields[f"cost_electricity_{name}_eur_per_t_iron"] == pytest.approx(
                sum(fields[f"cost_electricity_{leaf}_eur_per_t_iron"] for leaf in parts), rel=1e-9)
    assert pd.isna(fields.get("cost_electricity_battery_fom_eur_per_t_iron", float("nan")))

    # Ironmaking stops at the DRI: the press and the furnace after it carry the
    # rest of the origin's draw, and with it the rest of the losses.
    job = lambda name: fields[f"cost_electricity_{name}_eur_per_t_steel"]
    iron_jobs = job("reduction") + job("electrolyser")
    furnace_home = 0.0 if route.endswith("-export") else job("eaf")
    iron_share = iron_jobs / (iron_jobs + job("briquetting") + furnace_home)
    after_iron = job("eaf") + job("briquetting") + (1 - iron_share) * job("losses")
    assert total * iron_t / steel_t + after_iron == pytest.approx(
        fields["cost_electricity_eur_per_t_steel"], rel=1e-9)
    # The press per tonne of the DRI it pressed is its steel job over the iron.
    assert fields["cost_electricity_briquetting_eur_per_t_iron"] == pytest.approx(
        job("briquetting") * steel_t / iron_t, rel=1e-9)



@pytest.mark.parametrize("route", ROUTES)
def test_the_iron_s_electricity_by_job_closes_on_the_same_total(route, assumptions):
    """Ironmaking's three jobs, losses at its share, stack to its electricity total;
    the press is beside them, not one of them."""
    n = _solved(route, assumptions)
    fields, _, _ = _fields(n, assumptions)
    jobs = sum(fields[f"cost_electricity_{job}_eur_per_t_iron"] for job in IRON_ELECTRICITY_JOBS)
    assert jobs == pytest.approx(fields["cost_electricity_eur_per_t_iron"], rel=1e-9)
    assert "eaf" not in IRON_ELECTRICITY_JOBS
    assert "briquetting" not in IRON_ELECTRICITY_JOBS


def _iron_and_after(route, assumptions):
    """(fields, yearly € of each leaf, iron t, steel t, what the iron store holds)."""
    n = _solved(route, assumptions)
    fields, _, steel_t = _fields(n, assumptions)
    iron_t = compile_report._iron_produced_t(n) * 8760.0 / len(n.snapshots)
    yearly = {leaf: fields[f"cost_{leaf}_eur_per_t_steel"] * steel_t for leaf in LEAF_COSTS}
    holds = n.stores.at["iron_store", "bus"] if "iron_store" in n.stores.index else None
    return n, fields, yearly, iron_t, steel_t, holds


@pytest.mark.parametrize("route", ROUTES)
def test_the_iron_and_everything_after_it_are_the_whole_steel(route, assumptions):
    """The plant gate splits the cost of steel in two with nothing lost or doubled:
    the iron, and the press, furnace, stores, freight and the rest of the power."""
    _, fields, yearly, iron_t, _, holds = _iron_and_after(route, assumptions)
    origin_power = sum(yearly[leaf] for leaf in ORIGIN_ELECTRICITY_LEAVES)
    iron_power = fields["cost_electricity_eur_per_t_iron"] * iron_t
    after = (yearly["eaf_consumables"]
             + sum(yearly[f"{stem}_{half}"] for stem in ("steelmaking_plant", "briquetting_plant")
                   for half in ("capex", "fom"))
             + yearly["steel_store"] + yearly["transport"]
             + (yearly["iron_or_hbi_store"] if holds != "iron" else 0.0)
             + sum(v for leaf, v in yearly.items() if leaf.startswith("destination_grid_"))
             + origin_power - iron_power)
    assert fields["lcoi_eur_per_t_iron"] * iron_t + after == pytest.approx(
        sum(yearly.values()), rel=1e-9)


@pytest.mark.parametrize("route", [r for r in ROUTES if r.endswith("-export") and "dri" in r])
def test_the_briquettes_are_the_iron_plus_the_press_and_their_stockpile(route, assumptions):
    n, fields, yearly, iron_t, steel_t, holds = _iron_and_after(route, assumptions)
    hbi_t = float(-n.links_t.p1["briquetting"].sum()) * 8760.0 / len(n.snapshots)
    job = lambda name: fields[f"cost_electricity_{name}_eur_per_t_steel"] * steel_t
    press_share = job("briquetting") / (job("reduction") + job("electrolyser") + job("briquetting"))
    press = (yearly["briquetting_plant_capex"] + yearly["briquetting_plant_fom"]
             + job("briquetting") + job("losses") * press_share)

    assert holds == "hbi"
    # Its breakdown: the iron, the press's two halves, its power, the stockpile.
    parts = ("ironmaking", "briquetting_plant_capex", "briquetting_plant_fom",
             "electricity_briquetting", "iron_or_hbi_store")
    assert sum(fields[f"cost_{part}_eur_per_t_hbi"] for part in parts) == pytest.approx(
        fields["lcohbi_eur_per_t_hbi"], rel=1e-9)
    assert fields["cost_briquetting_plant_eur_per_t_hbi"] == pytest.approx(
        fields["cost_briquetting_plant_capex_eur_per_t_hbi"]
        + fields["cost_briquetting_plant_fom_eur_per_t_hbi"], rel=1e-9)
    # And per tonne of steel: the same euros, and the steel-basis rows that are
    # its parts add up to it.
    assert fields["cost_hbi_eur_per_t_steel"] * steel_t == pytest.approx(
        fields["lcohbi_eur_per_t_hbi"] * hbi_t, rel=1e-9)
    steel_parts = ("ironmaking", "briquetting_plant", "electricity_briquetting_incl_losses",
                   "iron_or_hbi_store")
    assert sum(fields[f"cost_{part}_eur_per_t_steel"] for part in steel_parts) == pytest.approx(
        fields["cost_hbi_eur_per_t_steel"], rel=1e-9)
    assert fields["lcohbi_eur_per_t_hbi"] * hbi_t == pytest.approx(
        fields["lcoi_eur_per_t_iron"] * iron_t + press + yearly["iron_or_hbi_store"], rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_the_cost_of_iron_is_the_sum_of_its_parts(route, assumptions):
    """The breakdown stacks to the levelised cost of iron, the hydrogen without its
    power (that is in the electricity) and the carbon apart from the gas."""
    _, fields, _, iron_t, steel_t, _ = _iron_and_after(route, assumptions)
    # The plant per tonne of steel is the same yearly euros as per tonne of iron.
    assert fields["cost_ironmaking_plant_eur_per_t_steel"] * steel_t == pytest.approx(
        fields["cost_ironmaking_plant_eur_per_t_iron"] * iron_t, rel=1e-9)
    parts = ("ore", "ironmaking_plant_capex", "ironmaking_plant_fom", "gas_purchase", "gas_carbon",
             "electrolyser", "h2_buffer", "electricity", "iron_store")
    assert sum(fields[f"cost_{part}_eur_per_t_iron"] for part in parts) == pytest.approx(
        fields["lcoi_eur_per_t_iron"], rel=1e-9)
    assert fields["cost_ironmaking_plant_eur_per_t_iron"] == pytest.approx(
        fields["cost_ironmaking_plant_capex_eur_per_t_iron"]
        + fields["cost_ironmaking_plant_fom_eur_per_t_iron"], rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_the_cost_of_iron_per_tonne_of_steel_is_its_parts_on_that_basis(route, assumptions):
    """The iron per tonne of steel is the same euros as per tonne of iron, and the
    steel-basis rows that are its parts add up to it."""
    _, fields, _, iron_t, steel_t, _ = _iron_and_after(route, assumptions)
    total = fields["cost_ironmaking_eur_per_t_steel"]

    assert total * steel_t == pytest.approx(fields["lcoi_eur_per_t_iron"] * iron_t, rel=1e-9)
    parts = ("ore", "ironmaking_plant", "gas_purchase", "gas_carbon", "electrolyser", "h2_buffer",
             "electricity_ironmaking", "ironmaking_store")
    assert sum(fields[f"cost_{part}_eur_per_t_steel"] for part in parts) == pytest.approx(
        total, rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_the_steel_s_electricity_is_ironmaking_with_its_losses_and_the_rest(route, assumptions):
    """Ironmaking (its losses in), the press, the furnace and the losses left over
    close on the steel's whole electricity bill."""
    _, fields, _, _, _, _ = _iron_and_after(route, assumptions)
    job = lambda name: fields[f"cost_electricity_{name}_eur_per_t_steel"]
    assert (job("ironmaking") + job("briquetting") + job("eaf")
            + job("losses_excl_ironmaking")) == pytest.approx(
        fields["cost_electricity_eur_per_t_steel"], rel=1e-9)



@pytest.mark.parametrize("route", ROUTES)
def test_the_iron_delivered_is_the_gate_plus_its_freight(route, assumptions):
    """At the furnace: the gate's product and its freight, over what arrived. A
    home-made route ships nothing first, so there it is the gate's own cost."""
    n, fields, yearly, iron_t, _, _ = _iron_and_after(route, assumptions)
    total = fields["lcoi_delivered_eur_per_t_iron_or_hbi"]
    assert total == pytest.approx(fields["cost_plant_gate_eur_per_t_iron_or_hbi"]
                                  + fields["cost_transport_eur_per_t_iron_or_hbi"], rel=1e-9)
    if "iron_transport" not in n.links.index:
        assert fields["cost_transport_eur_per_t_iron_or_hbi"] == 0.0
        assert total == pytest.approx(fields["lcoi_eur_per_t_iron"], rel=1e-9)
    else:
        arrived_t = float(-n.links_t.p1["iron_transport"].sum()) * 8760.0 / len(n.snapshots)
        assert fields["cost_transport_eur_per_t_iron_or_hbi"] * arrived_t == pytest.approx(
            yearly["transport"], rel=1e-9)
