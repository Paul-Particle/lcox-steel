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
    "gas_upstream_t_co2e_per_mwh": {"test": 0.1},
    "freight_kg_co2e_per_t_km": {"sea": {"test": 0.004}, "rail": {"test": 0.01}},
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
    """Write a solved dispatch straight onto the network, no optimiser involved.

    Merges, so the battery's links and a route's own links can be written
    separately without one erasing the other.
    """
    frame = getattr(n, component + "_t")[attr].copy()
    for name, series in values.items():
        frame[name] = pd.Series(list(series), index=n.snapshots)
    getattr(n, component + "_t")[attr] = frame


def _battery(n: pypsa.Network, flow: list[float] | None = None) -> None:
    """Add the battery the way build_network does, dispatched from a bus-side flow.

    `flow` is signed as the electricity bus sees it — positive when the battery
    supplies it, negative when it charges. Pass None to leave the dispatch
    frames empty, which is the shape a battery that never ran comes back as.
    """
    n.add("Bus", "electricity_battery", carrier="battery")
    n.add("Store", "battery", bus="electricity_battery", carrier="battery")
    n.add("Link", "battery_charger", bus0="electricity", bus1="electricity_battery",
          carrier="battery")
    n.add("Link", "battery_discharger", bus0="electricity_battery", bus1="electricity",
          carrier="battery")
    if flow is None:
        return
    signed = pd.Series(flow, index=n.snapshots)
    _dispatch(n, "links", "p0", {
        "battery_charger": (-signed).clip(lower=0.0),
        "battery_discharger": signed.clip(lower=0.0),
    })
    _dispatch(n, "links", "p1", {"battery_discharger": -signed.clip(lower=0.0)})


def _grid_mix(n: pypsa.Network) -> pd.DataFrame:
    """A `variant: emissions` series that works out to 0.8 t/MWh: four parts coal
    to one part wind, and the price column that rides along with it."""
    return pd.DataFrame(
        {"hard_coal": [40.0] * 4, "wind_onshore": [10.0] * 4, "price": [50.0] * 4},
        index=n.snapshots,
    )


def _breakdown(n: pypsa.Network, legs: dict | None = None,
               grid_mix: pd.DataFrame | None = None,
               destination_mix: pd.DataFrame | None = None) -> dict:
    """The function under test, with only the inputs a case varies named at the call.

    The factor table, the gas figures and the two place names are the same in
    every case, so they sit here rather than in eighteen argument lists."""
    return compile_report._emissions_breakdown(
        n, EMISSIONS, NATURAL_GAS, "VIC1", legs, grid_mix, "DEU", destination_mix
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

    emitted = _breakdown(n, grid_mix=_grid_mix(n))

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

    emitted = _breakdown(n, grid_mix=_grid_mix(n))

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

    emitted = _breakdown(n, grid_mix=_grid_mix(n))

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

    # The destination's own mix: half coal, half wind, so 0.5 t/MWh.
    destination_mix = pd.DataFrame(
        {"hard_coal": [5.0] * 4, "wind_onshore": [5.0] * 4, "price": [90.0] * 4},
        index=n.snapshots,
    )
    emitted = _breakdown(n, destination_mix=destination_mix)

    # The destination market's 0.5 t/MWh, not the 0 of the wind that made the iron.
    assert emitted["by_step"]["eaf"] == pytest.approx(6.0 * 4 * SCALE * 0.5)


def test_the_destination_furnace_carries_the_hours_it_melted_in():
    """The destination is a market with hours of its own, so a furnace that melts
    in its clean ones carries less than its year — the same cut as at home."""
    n = _network()
    n.add("Bus", "electricity_destination", carrier="AC")
    n.add("Bus", "iron_destination", carrier="iron")
    n.add("Generator", "destination_supply", bus="electricity_destination", carrier="AC")
    n.add("Link", "eaf", bus0="iron_destination", bus1="steel",
          bus2="electricity_destination")
    # Two dirty hours then two clean ones; the furnace only melts in the clean ones.
    _dispatch(n, "generators", "p", {"destination_supply": [0.0, 0.0, 6.0, 6.0]})
    _dispatch(n, "links", "p0", {"eaf": [0.0, 0.0, 5.0, 5.0]})
    _dispatch(n, "links", "p2", {"eaf": [0.0, 0.0, 6.0, 6.0]})
    destination_mix = pd.DataFrame(
        {"hard_coal": [10.0, 10.0, 0.0, 0.0],
         "wind_onshore": [0.0, 0.0, 10.0, 10.0],
         "price": [120.0, 120.0, 20.0, 20.0]},
        index=n.snapshots,
    )

    emitted = _breakdown(n, destination_mix=destination_mix)

    # Nothing, where the destination's flat annual average of 0.5 t/MWh would
    # have charged it for coal it never bought.
    assert emitted["by_step"]["eaf"] == pytest.approx(0.0)


def test_a_destination_furnace_without_its_market_s_series_is_an_error():
    """The furnace's power is most of an export route's carbon, so a missing
    destination series has to stop the report rather than read as clean."""
    n = _network()
    n.add("Bus", "electricity_destination", carrier="AC")
    n.add("Bus", "iron_destination", carrier="iron")
    n.add("Generator", "destination_supply", bus="electricity_destination", carrier="AC")
    n.add("Link", "eaf", bus0="iron_destination", bus1="steel",
          bus2="electricity_destination")
    _dispatch(n, "generators", "p", {"destination_supply": [6.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [6.0] * 4})

    with pytest.raises(ValueError, match="DEU"):
        _breakdown(n)


def test_freight_is_charged_per_tonne_over_the_run_s_own_legs():
    n = _network()
    n.add("Bus", "iron_destination", carrier="iron")
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Link", "iron_transport", bus0="iron", bus1="iron_destination")
    _dispatch(n, "generators", "p", {"wind-onshore": [0.0] * 4})
    _dispatch(n, "links", "p0", {"iron_transport": [1.0] * 4})

    emitted = _breakdown(n, legs={"sea": 10_000, "rail": 300})

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
    emitted = _breakdown(n, grid_mix=mix)
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
        _breakdown(n)

    # A price-only series is the same case: no carrier columns, no intensity.
    prices_only = pd.DataFrame({"price": [50.0] * 4}, index=n.snapshots)
    with pytest.raises(ValueError, match="no generation mix"):
        _breakdown(n, grid_mix=prices_only)


def test_an_islanded_run_needs_no_mix_at_all():
    """It imports nothing, so there is no grid intensity to be missing."""
    n = _network()
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"wind-onshore": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    emitted = _breakdown(n)
    assert emitted["by_step"]["eaf"] == pytest.approx(0.0)


def test_what_comes_out_of_storage_is_not_free():
    """An hour supplied out of the battery has nothing generating in it, so the
    draw in that hour would otherwise be charged at nothing. The run's total has
    to come to what its generation carried, whichever hour it was used in."""
    n = _network()
    n.add("Generator", "solar", bus="electricity", carrier="solar")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    # Sun in the odd hours only; the furnace runs through the dark ones on the
    # battery, which gives back 4 of every 5 MWh it takes.
    _battery(n, [-5.0, 4.0, -5.0, 4.0])
    _dispatch(n, "generators", "p", {"solar": [10.0, 0.0, 10.0, 0.0]})
    _dispatch(n, "links", "p0", {"eaf": [1.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [5.0, 4.0, 5.0, 4.0]})

    emitted = _breakdown(n)

    # 20 MWh of solar at 0.2 is all this system ever emitted, however it moved.
    assert sum(emitted["by_step"].values()) == pytest.approx(20.0 * SCALE * 0.2)
    assert emitted["by_step"]["eaf"] == pytest.approx(18.0 * SCALE * 0.2)
    assert emitted["by_step"]["battery_losses"] == pytest.approx(2.0 * SCALE * 0.2)
    assert emitted["losses_mwh"] == pytest.approx(2.0 * SCALE)


def test_charging_clean_dilutes_what_the_battery_gives_back():
    """The tank mixes, so its intensity moves as it is charged: coal in the
    morning comes back out as coal, and wind on top of it dilutes what is left.
    One average over the year would split this the wrong way between the two
    users — same total, different tonnes each."""
    n = _network()
    n.add("Bus", "hydrogen", carrier="H2")
    n.add("Generator", "coal", bus="electricity", carrier="hard_coal")
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Link", "electrolyser", bus0="electricity", bus1="hydrogen")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    # Fill on coal, give a third of it back, top up on wind, empty the rest.
    _battery(n, [-10.0, 5.0, -10.0, 15.0])
    _dispatch(n, "generators", "p", {
        "coal":         [10.0, 0.0, 0.0, 0.0],
        "wind-onshore": [0.0, 0.0, 10.0, 0.0],
    })
    _dispatch(n, "links", "p0", {"electrolyser": [0.0, 5.0, 0.0, 0.0], "eaf": [1.0] * 4})
    _dispatch(n, "links", "p2", {"electrolyser": [0.0] * 4, "eaf": [0.0, 0.0, 0.0, 15.0]})

    emitted = _breakdown(n)

    # The electrolyser drew before the wind arrived, so it got undiluted coal;
    # the furnace drew from a tank that was two thirds wind by then. A single
    # charge-weighted average would have said 2.5 and 7.5.
    assert emitted["by_step"]["electrolyser"] == pytest.approx(5.0 * SCALE)
    assert emitted["by_step"]["eaf"] == pytest.approx(5.0 * SCALE)
    assert emitted["by_step"].get("battery_losses", 0.0) == pytest.approx(0.0)
    assert sum(emitted["by_step"].values()) == pytest.approx(10.0 * SCALE)


def test_stored_energy_carries_the_hours_it_charged_in():
    """Not the year's average, which is a number the run does not contain. A
    battery filled on wind gives back wind, however dirty the rest of the year
    was — and its round-trip loss is charged at the same rate."""
    n = _network()
    n.add("Generator", "wind-onshore", bus="electricity", carrier="wind-onshore")
    n.add("Generator", "solar", bus="electricity", carrier="solar")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    # Hour 1 is all wind and fills the battery; hour 2 runs off it; hour 3 is
    # all solar. The year's mean intensity would be 0.1, halfway between them.
    _battery(n, [-5.0, 4.0, 0.0, 0.0])
    _dispatch(n, "generators", "p", {
        "wind-onshore": [10.0, 0.0, 0.0, 0.0],
        "solar":        [0.0, 0.0, 10.0, 0.0],
    })
    _dispatch(n, "links", "p0", {"eaf": [1.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [5.0, 4.0, 10.0, 0.0]})

    emitted = _breakdown(n)

    # The furnace's dark hour is free because the wind that filled the battery
    # was; only the solar hour costs anything.
    assert emitted["by_step"]["eaf"] == pytest.approx(10.0 * SCALE * 0.2)
    assert emitted["by_step"]["battery_losses"] == pytest.approx(0.0)
    assert sum(emitted["by_step"].values()) == pytest.approx(10.0 * SCALE * 0.2)


def test_a_battery_that_never_ran_is_not_a_missing_column():
    """PyPSA's netCDF export drops an all-zero dispatch column, so a run that
    built no battery comes back with one in the index and none in the frame. An
    optimum like that is ordinary, and used to take the whole report down."""
    n = _network()
    n.add("Generator", "solar", bus="electricity", carrier="solar")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _battery(n)
    _dispatch(n, "generators", "p", {"solar": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [1.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [5.0] * 4})
    # Neither link appears in the dispatch — the shape the round trip leaves behind.
    assert "battery_charger" not in n.links_t.p0.columns
    assert "battery_discharger" not in n.links_t.p1.columns

    emitted = _breakdown(n)
    assert emitted["by_step"]["battery_losses"] == pytest.approx(0.0)
    assert emitted["by_step"]["eaf"] == pytest.approx(5.0 * 4 * SCALE * 0.2)


def test_a_mix_from_another_window_is_an_error_not_a_clean_hour():
    """Reindexing a series onto snapshots it does not cover gives NaN, and a NaN
    intensity used to become a zero one — a provenance error reading as
    zero-carbon electricity."""
    n = _network()
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    other_year = pd.DataFrame(
        {"hard_coal": [40.0] * 4, "wind_onshore": [10.0] * 4, "price": [50.0] * 4},
        index=range(100, 104),
    )
    with pytest.raises(ValueError, match="different windows"):
        _breakdown(n, grid_mix=other_year)


def test_a_carrier_the_factor_table_does_not_know_is_an_error():
    """Dropping it would renormalise the mix over the rest and read as though the
    unknown carrier generated nothing — the config says it is an error."""
    n = _network()
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    mix = _grid_mix(n).assign(fusion=[100.0] * 4)
    with pytest.raises(ValueError, match="no entry for"):
        _breakdown(n, grid_mix=mix)


def test_storage_is_left_out_of_the_grid_mix_rather_than_counted_at_zero():
    """What a reservoir gives back was generated in some earlier hour that is
    already in the mix, so counting it at a factor of zero would dilute it."""
    n = _network()
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    # All the generation is coal; a third of what leaves the zone came out of
    # pumped storage, and its charging reads as negative generation.
    mix = pd.DataFrame(
        {"hard_coal": [30.0] * 4, "pumped_storage": [10.0] * 4,
         "pumped_storage_cons": [-12.0] * 4, "price": [50.0] * 4},
        index=n.snapshots,
    )
    emitted = _breakdown(n, grid_mix=mix)
    assert emitted["by_step"]["eaf"] == pytest.approx(10.0 * 4 * SCALE * 1.0)


def test_a_carrier_s_own_consumption_is_not_a_source_of_its_own():
    """ENTSO-E reports an "Actual Consumption" figure for any carrier that has
    one, not only for storage — Germany publishes one for solar and onshore wind.
    It is what those plants drew, so it is not in the mix, and it is certainly
    not a carrier the factor table should have to know."""
    n = _network()
    n.add("Generator", "grid_import", bus="electricity", carrier="AC")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    _dispatch(n, "generators", "p", {"grid_import": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [5.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [10.0] * 4})

    # Three quarters coal to one quarter wind, and both renewables draw a little.
    mix = pd.DataFrame(
        {"hard_coal": [30.0] * 4, "wind_onshore": [10.0] * 4,
         "wind_onshore_cons": [-2.0] * 4, "solar": [0.0] * 4,
         "solar_cons": [-1.0] * 4, "price": [50.0] * 4},
        index=n.snapshots,
    )
    emitted = _breakdown(n, grid_mix=mix)
    assert emitted["by_step"]["eaf"] == pytest.approx(10.0 * 4 * SCALE * 0.75)


def test_a_link_drawing_power_the_report_cannot_name_is_an_error():
    """The report has one column per declared user, and the run's total is the sum
    over them — so an undeclared drawing link would go missing from the total
    while every share still stacked to a tidy 100 %."""
    n = _network()
    n.add("Bus", "hydrogen", carrier="H2")
    n.add("Generator", "solar", bus="electricity", carrier="solar")
    n.add("Link", "eaf", bus0="iron", bus1="steel", bus2="electricity")
    n.add("Link", "compressor", bus0="electricity", bus1="hydrogen")
    _dispatch(n, "generators", "p", {"solar": [10.0] * 4})
    _dispatch(n, "links", "p0", {"eaf": [1.0] * 4, "compressor": [2.0] * 4})
    _dispatch(n, "links", "p2", {"eaf": [5.0] * 4, "compressor": [0.0] * 4})

    with pytest.raises(ValueError, match="ELECTRICITY_USERS"):
        _breakdown(n)


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
    _battery(n, [-2.0, -2.0, 1.0, 1.0])
    _dispatch(n, "generators", "p", {"grid_import": [40.0] * 4})
    _dispatch(n, "links", "p0", {
        "electrolyser": [10.0] * 4, "dri-ng": [8.0] * 4, "eaf": [5.0] * 4,
    })
    _dispatch(n, "links", "p2", {
        "electrolyser": [0.0] * 4, "dri-ng": [3.0] * 4, "eaf": [12.0] * 4,
    })
    # Charges twice as hard as it discharges: the difference is the round-trip loss.

    emitted = _breakdown(n, grid_mix=_grid_mix(n))

    assert set(emitted["by_step"]) <= set(EMISSION_STEPS)
    assert emitted["by_step"]["battery_losses"] == pytest.approx(2.0 * SCALE * 0.8)
    total = sum(emitted["by_step"].values())
    assert total == pytest.approx(sum(emitted["sources"].values()))
