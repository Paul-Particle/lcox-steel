"""Unit tests for the report's finest cut of the cost of steel.

The leaves and the alternative bands each have to stack back to LCOS. The
cost-breakdown page prints the *reported* total on top of the stack it draws, so
a group that grows a component with no leaf to put it in goes missing from the
bars while every share still reads 100 % — which is what happened to a blended
shaft's ore, a third of the cost of that route's steel.

Built on real networks from `build_network`, so the capital and marginal costs
are the ones that price a solve, with a solved state written in by hand rather
than optimised: closure is arithmetic and holds whatever the optimiser would
have picked, so the tests need no solver and no cutouts.
"""

import numpy as np
import pandas as pd
import pytest
import yaml

import compile_report  # sys.path set by conftest
from build_network import build_network
from common._report_schema import (
    ALT_LCOS_PARTS,
    LEAF_COSTS,
    LEAF_GROUP,
    LEAF_PARENTS,
    PROCESS_LINKS,
    field_stem,
)
from conftest import REPO_ROOT

HOURS = 24
# Every route that makes steel and can be built without a multi-site overlay.
# The blended shaft and the export twins are the ones that used to fall through.
ROUTES = ["h2-dri-eaf", "mix-dri-eaf", "ng-dri-eaf", "moe-eaf", "ew-eaf",
          "h2-dri-eaf-export", "moe-eaf-export", "ew-eaf-export",
          "mix-dri-eaf-export", "ng-dri-eaf-export"]


@pytest.fixture(scope="module")
def assumptions() -> dict:
    """The real assumptions, so the capital costs are composed as a solve's are."""
    return yaml.safe_load((REPO_ROOT / "config" / "assumptions.yaml").read_text())


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
    drawn = {link: float(n.links_t.p0[link].sum()) * annual
             for link in ("eaf", "electrolyser") if link in n.links.index}
    elec_mwh = float(
        n.generators_t.p[[g for g in n.generators.index if g != "gas_supply"]].sum().sum()
    ) * annual
    fields = compile_report._leaf_breakdown(
        n, assumptions, breakdown, steel_t, 55.0, elec_mwh, drawn,
    )
    return fields, breakdown, steel_t


@pytest.mark.parametrize("route", ROUTES)
def test_the_leaves_stack_to_the_cost_of_steel(route, assumptions):
    """Every group is split into leaves or carried into one, so nothing is lost."""
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)

    leaves = sum(fields[f"leaf_{leaf}_eur_per_t"] for leaf in LEAF_COSTS)
    assert leaves == pytest.approx(sum(breakdown.values()) / steel_t, rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_the_alternative_bands_stack_to_the_cost_of_steel(route, assumptions):
    """The other taxonomy closes on the same total, cut a different way."""
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)

    bands = sum(fields[f"alt_lcos_{part}_eur_per_t"] for part in ALT_LCOS_PARTS)
    assert bands == pytest.approx(sum(breakdown.values()) / steel_t, rel=1e-9)


@pytest.mark.parametrize("route", ROUTES)
def test_each_parent_group_is_its_own_leaves(route, assumptions):
    """The coarse reading beside the leaves is those leaves and nothing else."""
    n = _solved(route, assumptions)
    fields, _, _ = _fields(n, assumptions)

    for parent in LEAF_PARENTS:
        own = sum(value for key, value in fields.items()
                  if key.startswith("leaf_")
                  and LEAF_GROUP[key[len("leaf_"):-len("_eur_per_t")]] == parent)
        assert fields[f"group_{parent}_eur_per_t"] == pytest.approx(own, rel=1e-9)


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
        halves = (fields[f"leaf_{stem}_capex_eur_per_t"]
                  + fields[f"leaf_{stem}_fom_eur_per_t"])
        assert halves == pytest.approx(whole, rel=1e-9)
        # And neither half is the whole of it: a fixed-opex quote of about an
        # eighth of capex puts the upkeep well inside these bounds.
        assert 0.0 < fields[f"leaf_{stem}_fom_eur_per_t"] < whole


@pytest.mark.parametrize("route", ["mix-dri-eaf", "mix-dri-eaf-export"])
def test_a_blended_shaft_pays_for_its_ore(route, assumptions):
    """The ore leaf names every link that buys ore, the blended shaft included.

    `ore_eur_per_t_steel` used to name four of the five, so on this route both
    it and the leaf under it read zero while the run had paid for the ore.
    """
    n = _solved(route, assumptions)
    fields, breakdown, steel_t = _fields(n, assumptions)

    assert fields["leaf_ore_eur_per_t"] > 0
    assert (fields["leaf_ore_eur_per_t"] + fields["leaf_consumables_eur_per_t"]
            == pytest.approx(breakdown["ore_consumables"] / steel_t, rel=1e-9))


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
            55.0, 1e6, {"eaf": 1e5},
        )
