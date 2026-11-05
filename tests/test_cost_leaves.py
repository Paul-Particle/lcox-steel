"""Unit tests for the report's two fine cuts of the cost of steel.

`cost_*_eur_per_t` (by what was bought) and `el_*_eur_per_t` (by what it was
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

    leaves = sum(fields[f"cost_{leaf}_eur_per_t"] for leaf in LEAF_COSTS)
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

    jobs = sum(fields[f"el_{job}_eur_per_t"] for job in ELECTRICITY_JOBS)
    assert jobs == pytest.approx(fields["cost_electricity_eur_per_t"], rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_each_parent_group_is_its_own_leaves(route, assumptions):
    """The coarse reading beside the leaves is those leaves and nothing else."""
    n = _solved(route, assumptions)
    fields, _, _ = _fields(n, assumptions)

    for parent in LEAF_PARENTS:
        own = sum(fields[f"cost_{leaf}_eur_per_t"] for leaf in LEAF_COSTS
                  if LEAF_GROUP[leaf] == parent)
        assert fields[f"cost_{parent}_eur_per_t"] == pytest.approx(own, rel=1e-9)


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
        stem = field_stem(link)
        whole = float(n.links.at[link, "capital_cost"]
                      * n.links.at[link, "p_nom_opt"]) / steel_t
        halves = (fields[f"cost_{stem}_capex_eur_per_t"]
                  + fields[f"cost_{stem}_fom_eur_per_t"])
        assert halves == pytest.approx(whole, rel=1e-9)
        # And neither half is the whole of it: a fixed-opex quote of about an
        # eighth of capex puts the upkeep well inside these bounds.
        assert 0.0 < fields[f"cost_{stem}_fom_eur_per_t"] < whole


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

    assert fields["el_reduction_eur_per_t"] > 0
    assert fields["reduction_el_mwh_per_t_steel"] > 0
    assert (fields["el_reduction_eur_per_t"]
            > fields["el_handling_eur_per_t"] + fields["el_losses_eur_per_t"])


@pytest.mark.parametrize("route", ["mix-dri-eaf", "mix-dri-eaf-export"])
def test_a_blended_shaft_pays_for_its_ore(route, assumptions):
    """The ore leaf names every link that buys ore, the blended shaft included.

    `ORE_LINKS` has to name the blended shaft too, or on this route both the
    column and the leaf under it read zero while the run paid for the ore.
    """
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)

    assert fields["cost_ore_eur_per_t"] > 0
    assert (fields["cost_ore_eur_per_t"] + fields["cost_consumables_eur_per_t"]
            == pytest.approx(breakdown["ore_consumables"] / steel_t, rel=1e-9))


def test_the_electrolyser_pays_its_water_per_mwh_of_hydrogen_out(assumptions):
    """Quoted per MWh of H2, so a better efficiency does not change the water bill per MWh of H2."""
    electrolyser = _solved("h2-dri-eaf", assumptions).links.loc["electrolyser"]

    assert electrolyser["marginal_cost"] / electrolyser["efficiency"] == pytest.approx(
        assumptions["electrolyser"]["varopex_eur_per_mwh_h2"])


def test_the_hydrogen_buffer_pays_its_fixed_om_on_top_of_its_annuitised_capex(assumptions):
    """Per MWh of storage a year, as every other element's capital cost carries its fixed O&M."""
    buffer = _solved("h2-dri-eaf", assumptions).stores.loc["h2_buffer"]
    buffer_cfg = assumptions["h2_buffer"]
    annuitised_capex = (annuity_factor(assumptions["finance"]["origin_wacc"], buffer_cfg["lifetime_years"])
                        * buffer_cfg["capex_per_mwh_eur"])

    assert buffer["capital_cost"] == pytest.approx(annuitised_capex + buffer_cfg["opex_per_mwh_per_year_eur"])


def test_the_battery_pays_fixed_om_on_its_energy_and_its_inverter(assumptions):
    """Each capacity carries its own fixed O&M on top of its annuitised capex."""
    n = _solved("h2-dri-eaf", assumptions)
    battery_cfg = assumptions["battery"]
    annuity = annuity_factor(assumptions["finance"]["origin_wacc"], battery_cfg["lifetime_years"])

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
    leaves = sum(fields[f"cost_{leaf}_eur_per_t"] for leaf in LEAF_COSTS)
    assert leaves == pytest.approx(sum(breakdown.values()) / steel_t, rel=1e-9)


def summarise_written_in_run_with_generation_matching_the_draws(route: str, assumptions: dict) -> dict:
    """The report row of a written-in run whose home generation is exactly its draws and losses.

    The written-in flows draw more than the generators make, which the report's
    electricity jobs refuse, so the home generators are scaled to what was drawn.
    """
    n = _solved(route, assumptions)
    is_export = route.endswith("-export")
    run = {"run": "test", "area": "ESP", "destination": "DEU" if is_export else "ESP"}
    legs = {"sea": 8945} if is_export else None
    grid_mix = pd.DataFrame({"gas": 100.0, "solar": 100.0}, index=n.snapshots)
    emitted = compile_report._emissions_breakdown(n, assumptions["emissions"], run["area"], legs, grid_mix,
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
    """The furnace's power abroad made no hydrogen, so LCOH reads the home rate, as `el_hydrogen` does."""
    domestic = summarise_written_in_run_with_generation_matching_the_draws("h2-dri-eaf", assumptions)
    export = summarise_written_in_run_with_generation_matching_the_draws("h2-dri-eaf-export", assumptions)

    for summary in (domestic, export):
        h2_mwh = summary["h2_produced_kt"] * 1e6 * H2_LHV_KWH_PER_KG / 1000.0
        steel_t = summary["steel_produced_mt"] * 1e6
        assert summary["lcoh_electricity_eur_per_mwh_lhv"] * h2_mwh == pytest.approx(
            summary["el_hydrogen_eur_per_t"] * steel_t, rel=1e-9)
    # Both twins build the same plant at home, so they make hydrogen at the same cost.
    assert export["lcoh_eur_per_mwh_lhv"] == pytest.approx(domestic["lcoh_eur_per_mwh_lhv"], rel=1e-9)
    assert export["lcoe_eur_per_mwh"] != pytest.approx(domestic["lcoe_eur_per_mwh"], rel=1e-6)


@pytest.mark.parametrize("group", ["ladle_metallurgy", "process"])
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
