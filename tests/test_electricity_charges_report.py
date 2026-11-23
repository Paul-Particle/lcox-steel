"""The report charges the solved plant what its areas' workbooks bill: the grid
bill at the regime's rates, the connection, and what is charged on
consumption, with a grid or without.

On the same two-day networks and charges as test_electricity_charges_solve.
"""

import numpy as np
import pandas as pd
import pytest

from scripts.grid import _tariff
from scripts.solve.build_network import build_network
from scripts.solve.solve_network import _tie_battery_inverter
from scripts.viz import compile_report
from test_electricity_charges_solve import (CONNECTION, GRID, INDEX, TEMPLATE, _annuity, _charges,
                                            _regime, _solve)

MIX = pd.DataFrame({"gas": 400.0, "hard_coal": 100.0, "solar": np.where(INDEX.hour < 12, 500.0, 0.0)},
                   index=INDEX)
RUN = {"run": "DEU", "origin": "DEU", "destination": "DEU", "route": "h2-dri-eaf",
       "start_date": "20250301", "end_date": "20250302"}


# -- the report: the grid bill and the consumption charge reach the stack and the table


def _billed(consumption):
    """Charges with a component per kind of grid charge."""
    regime = _regime("full", 50_000, 10.0, fixed=1_000.0, components={
        "leistungspreis": {"capacity_eur_per_mw_per_year": 50_000},
        "arbeitspreis": {"energy_eur_per_mwh": 10.0},
        "messstellenbetrieb": {"fixed_eur_per_year": 1_000.0}})
    return _charges("DEU", regime, consumption=consumption)


def _report(consumption):
    charges = _billed(consumption)
    n, _ = _solve(origin=charges)
    s = compile_report.extract_summary(n, "tariff", RUN, TEMPLATE, None, MIX, None,
                                       {"origin": charges, "destination": None})
    return n, s, charges


def test_the_grid_bill_is_the_capacity_fee_and_the_fee():
    _, s, _ = _report(consumption=-20.0)
    assert s["grid_tariff_regime"] == "full"
    assert (s["cost_grid_capacity_fee_eur_per_t_steel"] + s["cost_grid_fee_eur_per_t_steel"]
            == pytest.approx(s["grid_charges_eur_per_t_steel"]))
    # The consumption charge is a line of its own, and the leaves still close on LCOS.
    assert s["cost_consumption_charges_eur_per_t_steel"] < 0
    leaves = sum(s[f"cost_{leaf}_eur_per_t_steel"] for leaf in compile_report.LEAF_COSTS)
    assert leaves == pytest.approx(s["lcos_eur_per_t_steel"])


def test_the_subsidy_lowers_the_lcos_by_what_it_paid_on_consumption():
    n, with_subsidy, _ = _report(consumption=-20.0)
    _, without, _ = _report(consumption=0.0)
    steel_t = with_subsidy["steel_produced_mt"] * 1e6
    consumed = with_subsidy["el_users_gwh"] * 1e3
    assert without["lcos_eur_per_t_steel"] - with_subsidy["lcos_eur_per_t_steel"] == pytest.approx(
        20.0 * consumed / steel_t)
    assert with_subsidy["cost_consumption_charges_eur_per_t_steel"] == pytest.approx(
        -20.0 * consumed / steel_t)


def test_the_subsidy_is_paid_on_what_the_plant_drew_not_on_what_it_imported():
    """The processes' draw, as the report counts it: the battery's losses are not
    consumption, so a plant that cycles its battery imports more than it consumed."""
    n, s, _ = _report(consumption=-20.0)
    consumed = compile_report._consumption_mwh(n, at_destination=False)
    assert consumed <= n.meta["grid_tariff"]["meters"]["grid_import"]["energy_mwh_per_year"] * (1 + 1e-9)
    assert consumed == pytest.approx(s["el_users_gwh"] * 1e3, rel=1e-6)


def test_a_run_is_charged_as_the_consumer_it_asks_for():
    """One area's charges for two consumers, as build_electricity_charges writes
    them: a run asked for `standard` is solved and reported on standard's rates,
    with no subsidy, and says so."""
    def own(capacity, energy, consumption):
        return {"grid": {"regimes": [_regime("full", capacity, energy, fixed=1_000.0)]},
                "consumption": {"eur_per_mwh": consumption,
                                "components": {"strompreiskompensation": consumption} if consumption else {}},
                "connection": None}
    area = {"area": "DEU", "source": "test", "component_groups": {},
            "consumers": {"iron_steel": own(50_000, 10.0, -20.0), "standard": own(60_000, 30.0, 0.0)}}
    reports = {}
    for consumer in ("iron_steel", "standard"):
        charges = _tariff.for_consumer(area, consumer)
        n, _ = _solve(origin=charges)
        reports[consumer] = compile_report.extract_summary(
            n, "tariff", RUN, TEMPLATE, None, MIX, None, {"origin": charges, "destination": None})

    steel, standard = reports["iron_steel"], reports["standard"]
    assert (steel["electricity_consumer"], standard["electricity_consumer"]) == ("iron_steel", "standard")
    assert standard["cost_consumption_charges_eur_per_t_steel"] == 0.0
    assert steel["cost_consumption_charges_eur_per_t_steel"] < 0
    assert standard["lcos_eur_per_t_steel"] > steel["lcos_eur_per_t_steel"]


def test_the_charges_table_lists_every_component_and_closes_on_each_total():
    n, _, charges = _report(consumption=-20.0)
    table = compile_report.electricity_charges(n, "tariff", "DEU", {"origin": charges})
    # Only the consumer the run was charged as.
    assert set(table["consumer"]) == {"iron_steel"}
    own = table
    for base, components in (("grid", ["leistungspreis", "arbeitspreis", "messstellenbetrieb"]),
                             ("consumption", ["strompreiskompensation"])):
        rows = own[own["base"] == base]
        parts = rows[rows["component"] != "total"]
        assert list(parts["component"]) == components
        assert parts["eur_per_year"].sum() == pytest.approx(
            rows.loc[rows["component"] == "total", "eur_per_year"].item())


def test_an_export_charges_each_side_on_its_own_consumption():
    """The furnace abroad is the destination's consumption, everything else the
    origin's; together they are every user's draw."""
    origin = _charges("AUS", _regime("a", 50_000, 10.0))
    destination = _charges("DEU", _regime("d", 50_000, 10.0), consumption=-20.0)
    n, _ = _solve("moe-eaf-export", origin, destination)
    run = {**RUN, "run": "AUS-DEU", "origin": "AUS", "route": "moe-eaf-export"}
    s = compile_report.extract_summary(n, "tariff", run, TEMPLATE, {"sea": 17_965}, MIX, MIX,
                                       {"origin": origin, "destination": destination})

    at_origin = compile_report._consumption_mwh(n, at_destination=False)
    at_destination = compile_report._consumption_mwh(n, at_destination=True)
    assert at_origin + at_destination == pytest.approx(s["el_users_gwh"] * 1e3, rel=1e-6)
    assert at_destination == pytest.approx(s["el_eaf_gwh"] * 1e3, rel=1e-6)
    # Only the destination's area charges consumption, so only its side pays.
    steel_t = s["steel_produced_mt"] * 1e6
    assert s["cost_consumption_charges_eur_per_t_steel"] == 0.0
    assert s["cost_destination_consumption_charges_eur_per_t_steel"] == pytest.approx(
        -20.0 * at_destination / steel_t)


def test_a_plant_on_its_own_renewables_is_charged_on_its_consumption_too():
    """No grid, so no tariff regime: the area's consumption charge still applies."""
    sun = pd.DataFrame({"solar": np.where((INDEX.hour >= 6) & (INDEX.hour < 18), 0.7, 0.0)},
                       index=INDEX)
    n = build_network("h2-dri-eaf", TEMPLATE, sun)
    status, condition = n.optimize(solver_name="highs", extra_functionality=_tie_battery_inverter)
    assert condition == "optimal"
    charges = _charges("DEU", _regime("full", 50_000, 10.0), consumption=-20.0)
    s = compile_report.extract_summary(n, "res-only", RUN, TEMPLATE, None, None, None,
                                       {"origin": charges, "destination": None})

    steel_t = s["steel_produced_mt"] * 1e6
    assert "grid_import" not in n.generators.index
    assert s["cost_consumption_charges_eur_per_t_steel"] == pytest.approx(
        -20.0 * s["el_users_gwh"] * 1e3 / steel_t)
    leaves = sum(s[f"cost_{leaf}_eur_per_t_steel"] for leaf in compile_report.LEAF_COSTS)
    assert leaves == pytest.approx(s["lcos_eur_per_t_steel"])


def test_the_connection_capex_leaf_holds_the_fixed_sum_and_the_leaves_close():
    with_sum = _charges("DEU", _regime("full", 50_000, 10.0), connection=CONNECTION)
    without_sum = _charges("DEU", _regime("full", 50_000, 10.0),
                           connection={**CONNECTION, "capex_eur": 0.0,
                                       "components": {name: {**parts, "capex_eur": 0.0}
                                                      for name, parts in CONNECTION["components"].items()}})
    summaries = {}
    for name, charges in (("with", with_sum), ("without", without_sum)):
        n, _ = _solve(origin=charges)
        summaries[name] = compile_report.extract_summary(n, "tariff", RUN, TEMPLATE, None, MIX, None,
                                                         {"origin": charges, "destination": None})
    s = summaries["with"]
    steel_t = s["steel_produced_mt"] * 1e6
    lump_sum = _annuity(TEMPLATE["finance"]["origin_wacc"]) * 2_000_000.0

    # The same plant (the fixed sum moves no decision), its capex leaf higher by the sum.
    assert (s["cost_grid_connection_capex_eur_per_t_steel"]
            - summaries["without"]["cost_grid_connection_capex_eur_per_t_steel"]) == pytest.approx(
        lump_sum / steel_t)
    assert s["cost_grid_market_eur_per_t_steel"] == pytest.approx(
        summaries["without"]["cost_grid_market_eur_per_t_steel"])
    leaves = sum(s[f"cost_{leaf}_eur_per_t_steel"] for leaf in compile_report.LEAF_COSTS)
    assert leaves == pytest.approx(s["lcos_eur_per_t_steel"])
    # The per-MWh halves of the connection still add up to the whole.
    assert (s["grid_connection_capex_eur_per_mwh_imported"] + s["grid_capacity_fee_eur_per_mwh_imported"]
            == pytest.approx(s["grid_connection_eur_per_mwh_imported"]))


def test_the_charges_table_lists_the_connection_by_component():
    charges = _charges("DEU", _regime("full", 50_000, 10.0), connection=CONNECTION)
    n, result = _solve(origin=charges)
    table = compile_report.electricity_charges(n, "tariff", "DEU", {"origin": charges})
    rows = table[table["base"] == "connection"]
    parts = rows[rows["component"] != "total"]

    assert list(parts["component"]) == ["netzanschlusskosten", "baukostenzuschuss"]
    assert set(parts["group"]) == {"grid_connection"}
    assert parts["eur_per_year"].sum() == pytest.approx(
        rows.loc[rows["component"] == "total", "eur_per_year"].item())
    connection = result["meters"]["grid_import"]["connection"]
    assert rows.loc[rows["component"] == "total", "eur_per_year"].item() == pytest.approx(
        connection["capex_eur_per_year"] + connection["lump_sum_eur_per_year"])
