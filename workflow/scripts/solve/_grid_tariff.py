"""Price a run's grid imports by tariff regime and pick the cheapest.

A grid connection is a "meter": the origin's grid import and, on an export
route, the destination furnace's supply. Each is billed on its own area's
tariff, which `build_electricity_charges` writes as regimes (see
scripts/grid/_tariff.py): within one, the bill is linear in the connection's
peak (p_nom) and its annual energy E, and holds only inside an (h, E) box,
h = E / p_nom. The solve runs once per regime, or per pair of regimes on an
export, with the box as linear constraints, and keeps the cheapest. An area
without a workbook has one regime, `config`: the `grid` block of the
assumptions with no box, which is what the model priced before tariffs.

E is annualised with the objective weightings, which build_network sets so the
window's hours count as one average year, as the tariffs' bands are per year.

The connection is paid once when it is built, the same whichever regime bills
it: a part per MW, annualised into the meter's capital cost so the LP sizes the
connection with it, and a fixed sum. An area without a workbook pays the
config's `grid.connection_capex_eur_per_mw` and no fixed sum.

Two costs are outside the LP (`outside_lp_eur_per_year`): the regime's fixed
€/yr (metering, the fixed term of the tiered surcharges) and the connection's
fixed sum. Both are owed for building a connection at all, a yes/no decision
an LP could only price with a binary variable (a MILP). But every regime of a
workbook area holds from 1 GWh of import a year (see scripts/grid/_tariff.py),
so a grid run always builds a connection, and both are constants: added to
every regime's cost for the comparison, which is the same as having them in
the objective, as a constant moves no decision. An area without a workbook has
neither. A plant with no connection at all is the RES-only scenario of the
same route, technology and site.

An area's charges here are those of the one consumer the run is charged as
(`scripts.grid._tariff.for_consumer`, the config's
`electricity_charges.consumer`), and only their `grid` and `connection` parts
are read. What the workbook charges on consumption (the Industriestrompreis,
the Strompreiskompensation) is the report's business: it applies to every run
in the area, with a grid or without, and moves no decision the LP takes.
"""

from __future__ import annotations

import itertools
import math

import pypsa
import xarray as xr

from scripts.solve._helpers_solve import annuity_factor

# The meters, by the generator build_network adds for each.
GRID_IMPORT = "grid_import"
DESTINATION_SUPPLY = "destination_supply"


def config_regime(grid_cfg: dict) -> dict:
    """The one regime of an area without a tariff workbook."""
    return {
        "name": "config",
        "h_min": None, "h_max": None, "e_min_mwh": None, "e_max_mwh": None,
        "capacity_eur_per_mw_per_year": grid_cfg["fee_eur_per_mw_per_year"],
        "energy_eur_per_mwh": grid_cfg["fee_eur_per_mwh"],
        "fixed_eur_per_year": 0.0,
        # Named so the report's list of charges has the config's two fees as lines.
        "components": {
            "config_capacity_fee": {"capacity_eur_per_mw_per_year": grid_cfg["fee_eur_per_mw_per_year"]},
            "config_volumetric_fee": {"energy_eur_per_mwh": grid_cfg["fee_eur_per_mwh"]},
        },
    }


def _from_workbook(charges: dict | None) -> bool:
    """Whether an area's charges come from a workbook: an area without one has no regimes."""
    return bool(charges) and (charges.get("grid") or {}).get("regimes") is not None


def regimes_of(charges: dict | None, grid_cfg: dict) -> list[dict]:
    """An area's grid regimes: its workbook's, or the config's one."""
    if not _from_workbook(charges):
        return [config_regime(grid_cfg)]
    return charges["grid"]["regimes"]


def connection_of(charges: dict | None, grid_cfg: dict) -> dict:
    """An area's one-off connection costs: its workbook's, or the config's €/MW."""
    if charges and charges.get("connection"):
        return charges["connection"]
    return {
        "capex_eur_per_mw": grid_cfg["connection_capex_eur_per_mw"],
        "capex_eur": 0.0,
        "components": {"config_connection_capex": {"capex_eur_per_mw": grid_cfg["connection_capex_eur_per_mw"],
                                                   "capex_eur": 0.0}},
    }


def _connection_annuity(grid_cfg: dict, wacc: float) -> float:
    return annuity_factor(wacc, grid_cfg["connection_financing_period"])


def lump_sum_eur_per_year(connection: dict, grid_cfg: dict, wacc: float) -> float:
    """The connection's fixed one-off cost, annualised: a constant the LP leaves out.

    Owed for building a connection at all, whatever its size, and a grid run of
    a workbook area always builds one, as its regimes hold from 1 GWh of import
    a year. Zero for an area without a workbook.
    """
    return _connection_annuity(grid_cfg, wacc) * float(connection["capex_eur"])


def outside_lp_eur_per_year(regime: dict, connection: dict, grid_cfg: dict, wacc: float) -> float:
    """What a meter costs a year that the LP does not price: the regime's fixed
    €/yr and the connection's lump sum, both constants of the regime."""
    return float(regime["fixed_eur_per_year"]) + lump_sum_eur_per_year(connection, grid_cfg, wacc)


def regime_bill(regime: dict, p_nom: float, energy: float) -> float:
    """What a regime charges a connection of p_nom MW that imports `energy` MWh a year."""
    return (regime["capacity_eur_per_mw_per_year"] * p_nom
            + regime["energy_eur_per_mwh"] * energy + regime["fixed_eur_per_year"])


def meters(n: pypsa.Network) -> list[str]:
    """The grid connections the network has: origin import, destination supply, or both."""
    return [name for name in (GRID_IMPORT, DESTINATION_SUPPLY) if name in n.generators.index]


def candidates(regimes_by_meter: dict[str, list[dict]]) -> list[dict[str, dict]]:
    """Every combination of one regime per meter.

    A domestic DEU run has one meter and 5 combinations; an AUS→DEU export, 1
    regime at the origin (no workbook) times 5 at the destination.
    """
    names = list(regimes_by_meter)
    return [dict(zip(names, combo)) for combo in itertools.product(*regimes_by_meter.values())]


def apply_regime(n: pypsa.Network, meter: str, regime: dict, price, grid_cfg: dict, wacc: float,
                 connection: dict) -> None:
    """Set a meter's capacity and energy charges to the regime's, its capital
    cost to the connection's annualised €/MW plus the regime's capacity charge.

    The same two prices build_network's _add_grid_import sets from the config's
    fees: `capital_cost` per MW of the connection a year, `marginal_cost` per
    MWh imported, the hourly market price plus the regime's energy rate.
    """
    n.generators.at[meter, "capital_cost"] = (
        _connection_annuity(grid_cfg, wacc) * float(connection["capex_eur_per_mw"])
        + regime["capacity_eur_per_mw_per_year"]
    )
    n.generators_t.marginal_cost[meter] = price + regime["energy_eur_per_mwh"]


def _bound(value, default: float) -> float:
    """A regime's bound from the YAML, where an open one is null: `default` then."""
    return default if value is None else float(value)


def add_regime_constraints(n: pypsa.Network, meter: str, regime: dict) -> None:
    """Keep the meter's (h, E) inside the regime's box. Call from extra_functionality.

    E is the meter's import over one average year: the hourly import weighted
    with the objective weightings, which build_network sets so a window of any
    length counts as one year, as the tariff's bands do. h = E / p_nom is a
    ratio of two variables, so its bounds go in multiplied out
    (E ≥ h_min · p_nom, E ≤ h_max · p_nom), which keeps them linear. A bound
    left as None in the tariff adds no constraint, and a regime without any —
    the config's — leaves the LP as it was.
    """
    h_min, h_max = _bound(regime["h_min"], 0.0), _bound(regime["h_max"], math.inf)
    e_min, e_max = _bound(regime["e_min_mwh"], 0.0), _bound(regime["e_max_mwh"], math.inf)
    if h_min == 0.0 and h_max == math.inf and e_min == 0.0 and e_max == math.inf:
        return

    # PyPSA's names for the LP's variables: the hourly import and the size.
    model = n.model
    weights = xr.DataArray(
        n.snapshot_weightings["objective"].values, coords={"snapshot": n.snapshots}, dims="snapshot"
    )
    energy = (model.variables["Generator-p"].loc[:, meter] * weights).sum("snapshot")
    p_nom = model.variables["Generator-p_nom"].loc[meter]
    tag = f"{meter}_tariff"
    if h_min > 0.0:
        model.add_constraints(energy - h_min * p_nom >= 0, name=f"{tag}_h_min")
    if h_max < math.inf:
        model.add_constraints(energy - h_max * p_nom <= 0, name=f"{tag}_h_max")
    if e_min > 0.0:
        model.add_constraints(energy >= e_min, name=f"{tag}_e_min")
    if e_max < math.inf:
        model.add_constraints(energy <= e_max, name=f"{tag}_e_max")


def annual_energy_mwh(n: pypsa.Network, meter: str) -> float:
    """The meter's solved import over one average year."""
    return float((n.generators_t.p[meter] * n.snapshot_weightings["objective"]).sum())


def connection_outcome(connection: dict, grid_cfg: dict, wacc: float, p_nom: float,
                       groups: dict) -> dict:
    """What building the connection came to per year, for n.meta and the report.

    The €/MW part is in the meter's capital cost, which the LP priced; the fixed
    sum is a constant the LP leaves out, charged here. The lines split both by
    component and add up to their total.
    """
    annuity = _connection_annuity(grid_cfg, wacc)
    lines = [{"component": component,
              "group": groups.get(component),
              "eur_per_year": annuity * (float(parts.get("capex_eur_per_mw", 0.0)) * p_nom
                                         + float(parts.get("capex_eur", 0.0)))}
             for component, parts in connection["components"].items()]
    return {
        "capex_eur_per_mw": float(connection["capex_eur_per_mw"]),
        "capex_eur": float(connection["capex_eur"]),
        "annuity_factor": annuity,
        "capex_eur_per_year": annuity * float(connection["capex_eur_per_mw"]) * p_nom,
        "lump_sum_eur_per_year": lump_sum_eur_per_year(connection, grid_cfg, wacc),
        "charges": lines,
    }


def outcome(n: pypsa.Network, meter: str, regime: dict, charges: dict | None,
            grid_cfg: dict, wacc: float) -> dict:
    """What the meter's grid bill and its connection came to under its regime,
    for n.meta and the report."""
    energy = annual_energy_mwh(n, meter)
    p_nom = float(n.generators.at[meter, "p_nom_opt"])
    peak = float(n.generators_t.p[meter].max())
    from_workbook = _from_workbook(charges)
    groups = (charges or {}).get("component_groups") or {}

    # The bill line by line, each component at its own rates: a tiered
    # surcharge's fixed offset stays with that surcharge. The lines add up to
    # the regime's bill.
    lines = []
    for component, parts in regime["components"].items():
        amount = (parts.get("capacity_eur_per_mw_per_year", 0.0) * p_nom
                  + parts.get("energy_eur_per_mwh", 0.0) * energy
                  + parts.get("fixed_eur_per_year", 0.0))
        lines.append({"component": component,
                      "group": groups.get(component) if from_workbook else "config",
                      "eur_per_year": amount})
    bill = regime_bill(regime, p_nom, energy)

    return {
        "area": charges["area"] if charges else None,
        "source": charges["source"] if charges else "config",
        "consumer": charges.get("consumer") if charges else None,
        "regime": regime["name"],
        "h_min": regime["h_min"], "h_max": regime["h_max"],
        "e_min_mwh": regime["e_min_mwh"], "e_max_mwh": regime["e_max_mwh"],
        "capacity_eur_per_mw_per_year": regime["capacity_eur_per_mw_per_year"],
        "energy_eur_per_mwh": regime["energy_eur_per_mwh"],
        "fixed_eur_per_year": regime["fixed_eur_per_year"],
        "components": regime["components"],
        "p_nom_mw": p_nom,
        "peak_mw": peak,
        "energy_mwh_per_year": energy,
        "utilisation_hours": energy / p_nom if p_nom > 0 else None,
        "charges": lines,
        "bill_eur_per_year": bill,
        # Not part of the tariff bill: what building the connection cost a year.
        "connection": connection_outcome(connection_of(charges, grid_cfg), grid_cfg, wacc, p_nom,
                                         groups if from_workbook else {}),
    }
