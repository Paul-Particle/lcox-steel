"""Unit tests for how a run's emissions are attributed.

The interesting part is not the multiplication, it is who gets charged: a link
that draws electricity as its input against one that draws it alongside another
conversion, a furnace on the far side of an ocean against one at the plant, and
a shaft that burns gas as well as drawing power. These pin that on hand-built
networks with the dispatch written in — no cutouts, no solves, so they run
anywhere and a change in `config/assumptions.yaml` cannot move them.
"""

import pandas as pd
import pypsa
import pytest

import compile_report  # sys.path set by conftest
from common._report_schema import EMISSION_STEPS

# A deliberately round factor table: coal at 1, gas at 0.5, wind free, so every
# expected number below can be read off by hand.
EMISSIONS = {
    "basis": "test",
    "electricity_t_co2e_per_mwh": {
        "hard_coal":    {"test": 1.0},
        "gas":          {"test": 0.5},
        "wind_onshore": {"test": 0.0},
        "solar":        {"test": 0.2},
    },
    "grid_t_co2e_per_mwh": {"VIC1": 0.8},
    "destination_t_co2e_per_mwh": 0.4,
    "gas_upstream_t_co2e_per_mwh": {"test": 0.1},
    "freight_kg_co2e_per_t_km": {"sea": 0.004, "rail": 0.01},
}
NATURAL_GAS = {"co2_t_per_mwh": 0.2}
# Four snapshots stand for the year, so an hourly MW is 2190 MWh annualised.
SCALE = 8760 / 4


def _network() -> pypsa.Network:
    """A two-bus skeleton: an electricity bus feeding an iron bus through an EAF."""
    n = pypsa.Network()
    n.set_snapshots(range(4))
    n.add("Bus", "electricity", carrier="AC")
    n.add("Bus", "iron", carrier="iron")
    n.add("Bus", "steel", carrier="steel")
    return n


def _dispatch(n: pypsa.Network, component: str, attr: str, values: dict) -> None:
    """Write a solved dispatch straight onto the network, no optimiser involved."""
    frame = pd.DataFrame(values, index=n.snapshots)
    getattr(n, component + "_t")[attr] = frame


def _grid_mix(n: pypsa.Network) -> pd.DataFrame:
    """A `variant: emissions` series that works out to 0.8 t/MWh: four parts coal
    to one part wind, and the price column that rides along with it."""
    return pd.DataFrame(
        {"hard_coal": [40.0] * 4, "wind_onshore": [10.0] * 4, "price": [50.0] * 4},
        index=n.snapshots,
    )


def test_a_bus0_draw_and_a_bus2_draw_are_both_charged():
    """An electrolyser takes power as its input, a furnace as a by-draw. Both count."""
    n = _network()
    n.add("Bus", "hydrogen", carrier="H2")
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "electrolyser", bus0="electricity", bus1="hydrogen")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [30.0] * 4})
    _dispatch(n, "links", "p0", {"electrolyser": [10.0] * 4, "eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"electrolyser": [0.0] * 4, "eaf": [20.0] * 4})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, _grid_mix(n)
    )

    assert emitted["electricity_mwh"]["electrolyser"] == pytest.approx(10.0 * 4 * SCALE)
    assert emitted["electricity_mwh"]["eaf"] == pytest.approx(20.0 * 4 * SCALE)
    # Both on the grid's 0.8 t/MWh, because that is all that supplied them.
    assert emitted["by_step"]["eaf"] == pytest.approx(20.0 * 4 * SCALE * 0.8)


def test_a_flat_user_and_a_flexible_one_carry_different_intensities():
    """Whoever runs when the wind blows carries less. This is the point of the hourly cut."""
    n = _network()
    n.add("Bus", "hydrogen", carrier="H2")
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "electrolyser", bus0="electricity", bus1="hydrogen")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    # Two windy hours then two still ones; the electrolyser only runs in the wind.
    _dispatch(n, "generators", "p", {
        "wind-onshore": [10.0, 10.0, 0.0, 0.0],
        "grid_import":  [0.0, 0.0, 10.0, 10.0],
    })
    _dispatch(n, "links", "p0", {"electrolyser": [5.0, 5.0, 0.0, 0.0], "eaf": [1.0] * 4})
    _dispatch(n, "links", "p2", {"electrolyser": [0.0] * 4, "eaf": [5.0] * 4})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, _grid_mix(n)
    )

    # The electrolyser drew only in free hours; the furnace drew half and half.
    assert emitted["by_step"]["electrolyser"] == pytest.approx(0.0)
    assert emitted["by_step"]["eaf"] == pytest.approx(2 * 5.0 * SCALE * 0.8)


def test_a_gas_shaft_carries_its_combustion_but_not_in_its_power_intensity():
    """Its step total includes the gas; its per-MWh figure must not, or it reads as
    having bought impossibly dirty electricity."""
    n = _network()
    n.add("Bus", "gas", carrier="gas")
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "dri-ng", bus0="gas", bus1="iron", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [2.0] * 4})
    _dispatch(n, "links", "p0", {"dri-ng": [10.0] * 4})
    _dispatch(n, "links", "p2", {"dri-ng": [2.0] * 4})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, _grid_mix(n)
    )

    gas_t = 10.0 * 4 * SCALE * (0.2 + 0.1)
    power_t = 2.0 * 4 * SCALE * 0.8
    assert emitted["by_step"]["dri-ng"] == pytest.approx(gas_t + power_t)
    assert emitted["electricity_t"]["dri-ng"] == pytest.approx(power_t)
    assert emitted["sources"]["gas"] == pytest.approx(gas_t)


def test_a_destination_furnace_is_charged_to_the_destination_grid():
    """An export route's furnace runs where the iron lands, not where it was made."""
    n = _network()
    n.add("Bus", "electricity_destination", carrier="AC")
    n.add("Bus", "iron_destination", carrier="iron")
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Generator", "destination_supply", bus="electricity_destination", carrier="AC")
    n.add("Link", "eaf", bus0="iron_destination", bus1="steel",
          bus2="electricity_destination")
    _dispatch(n, "generators", "p", {
        "wind-onshore": [10.0] * 4, "destination_supply": [6.0] * 4,
    })
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [6.0] * 4})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, None
    )

    # 0.4 t/MWh at the destination, not the 0 of the wind that made the iron.
    assert emitted["by_step"]["eaf"] == pytest.approx(6.0 * 4 * SCALE * 0.4)


def test_freight_is_charged_per_tonne_over_the_run_s_own_legs():
    n = _network()
    n.add("Bus", "iron_destination", carrier="iron")
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Link", "iron_transport", bus0="iron", bus1="iron_destination")
    _dispatch(n, "generators", "p", {"wind-onshore": [0.0] * 4})
    _dispatch(n, "links", "p0", {"iron_transport": [1.0] * 4})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", {"sea": 10_000, "rail": 300}, None
    )

    per_t = (0.004 * 10_000 + 0.01 * 300) / 1000.0
    assert emitted["by_step"]["iron_transport"] == pytest.approx(1.0 * 4 * SCALE * per_t)
    assert emitted["sources"]["freight"] == pytest.approx(emitted["by_step"]["iron_transport"])


def test_the_grid_s_intensity_comes_from_its_own_generation_mix():
    """The carrier columns of a `variant: emissions` series, weighted by the table.
    A price column rides along in that series and must not be mistaken for one."""
    n = _network()
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    # Three quarters coal, one quarter wind, plus the price column to ignore.
    mix = pd.DataFrame(
        {"hard_coal": [30.0] * 4, "wind_onshore": [10.0] * 4, "price": [50.0] * 4},
        index=n.snapshots,
    )
    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, mix
    )
    assert emitted["by_step"]["eaf"] == pytest.approx(10.0 * 4 * SCALE * 0.75)


def test_a_grid_run_without_a_mix_is_an_error_not_a_default():
    """There is nothing to stand in for a missing mix, so say so and stop. The
    message has to name the fix, because the fix is one cell of scenarios.csv."""
    n = _network()
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    with pytest.raises(ValueError, match="variant: emissions"):
        compile_report._emissions_breakdown(n, EMISSIONS, NATURAL_GAS, "VIC1", None, None)

    # A price-only series is the same case: no carrier columns, no intensity.
    prices_only = pd.DataFrame({"price": [50.0] * 4}, index=n.snapshots)
    with pytest.raises(ValueError, match="no generation mix"):
        compile_report._emissions_breakdown(
            n, EMISSIONS, NATURAL_GAS, "VIC1", None, prices_only
        )


def test_an_islanded_run_needs_no_mix_at_all():
    """It imports nothing, so there is no grid intensity to be missing."""
    n = _network()
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"wind-onshore": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, None
    )
    assert emitted["by_step"]["eaf"] == pytest.approx(0.0)


def test_the_steps_stack_to_the_total():
    """Every emitting thing lands in exactly one declared step, so the shares mean
    something. A step this network has none of contributes nothing."""
    n = _network()
    n.add("Bus", "gas", carrier="gas")
    n.add("Bus", "hydrogen", carrier="H2")
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "electrolyser", bus0="electricity", bus1="hydrogen")
    n.add("Link", "dri-ng", bus0="gas", bus1="iron", bus2="electricity")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    n.add("StorageUnit", "battery", bus="electricity", carrier="battery")
    _dispatch(n, "generators", "p", {"grid_import": [40.0] * 4})
    _dispatch(n, "links", "p0", {
        "electrolyser": [10.0] * 4, "dri-ng": [8.0] * 4, "eaf": [5.0] * 4,
    })
    _dispatch(n, "links", "p2", {
        "electrolyser": [0.0] * 4, "dri-ng": [3.0] * 4, "eaf": [12.0] * 4,
    })
    # Charges twice as hard as it discharges: the difference is the round-trip loss.
    _dispatch(n, "storage_units", "p", {"battery": [-2.0, -2.0, 1.0, 1.0]})

    emitted = compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", None, _grid_mix(n)
    )

    assert set(emitted["by_step"]) <= set(EMISSION_STEPS)
    assert emitted["by_step"]["battery_losses"] == pytest.approx(2.0 * SCALE * 0.8)
    total = sum(emitted["by_step"].values())
    assert total == pytest.approx(sum(emitted["sources"].values()))
