"""Price a run's grid imports by tariff regime and pick the cheapest.

`build_grid_tariff` writes each area's tariff as regimes: within one, the bill
is linear in the connection's peak (p_nom) and its annual energy E, and holds
only inside an (h, E) box, h = E / p_nom. The solve runs once per regime, with
that box as linear constraints, and keeps the solve whose objective plus the
regime's fixed €/yr is lowest. An export run has two meters, the origin's grid
import and the destination's supply, so it runs every pair of their regimes.

An area without a workbook has one regime, `config`: the `grid` block of the
assumptions with no box, which is what the model priced before tariffs.

E is annualised with the objective weightings, which build_network sets so the
window's hours count as one average year — the tariffs' bands are per year.

The subsidy rate of a tariff stays out of the LP (see scripts/grid/_tariff.py)
and is reported per meter after the solve.
"""

from __future__ import annotations

import itertools
import math

import pypsa
import xarray as xr

from scripts.solve._helpers_solve import annuity_factor

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
        "components": {},
    }


def regimes_of(tariff: dict | None, grid_cfg: dict) -> list[dict]:
    """An area's regimes: its workbook's, or the config's one."""
    if tariff is None or tariff.get("regimes") is None:
        return [config_regime(grid_cfg)]
    return tariff["regimes"]


def subsidy_of(tariff: dict | None) -> tuple[float, list[str]]:
    """The €/MWh paid back after the solve, and what it is made of."""
    if tariff is None or tariff.get("regimes") is None:
        return 0.0, []
    post = tariff.get("post_solve") or {}
    return float(post.get("subsidy_eur_per_mwh") or 0.0), list(post.get("subsidy_components") or [])


def meters(n: pypsa.Network) -> list[str]:
    """The grid connections the network has: origin import, destination supply, or both."""
    return [name for name in (GRID_IMPORT, DESTINATION_SUPPLY) if name in n.generators.index]


def candidates(regimes_by_meter: dict[str, list[dict]]) -> list[dict[str, dict]]:
    """Every combination of one regime per meter."""
    names = list(regimes_by_meter)
    return [dict(zip(names, combo)) for combo in itertools.product(*regimes_by_meter.values())]


def apply_regime(n: pypsa.Network, meter: str, regime: dict, price, grid_cfg: dict, wacc: float) -> None:
    """Set a meter's capacity and energy charges to the regime's."""
    n.generators.at[meter, "capital_cost"] = (
        annuity_factor(wacc, grid_cfg["connection_financing_period"])
        * grid_cfg["connection_capex_eur_per_mw"]
        + regime["capacity_eur_per_mw_per_year"]
    )
    n.generators_t.marginal_cost[meter] = price + regime["energy_eur_per_mwh"]


def _bound(value, default: float) -> float:
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


def outcome(n: pypsa.Network, meter: str, regime: dict, tariff: dict | None) -> dict:
    """What the meter's bill came to under its regime, for n.meta and the report."""
    energy = annual_energy_mwh(n, meter)
    p_nom = float(n.generators.at[meter, "p_nom_opt"])
    peak = float(n.generators_t.p[meter].max())
    subsidy, subsidy_components = subsidy_of(tariff)
    return {
        "area": tariff["area"] if tariff else None,
        "source": tariff["source"] if tariff else "config",
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
        "subsidy_eur_per_mwh": subsidy,
        "subsidy_components": subsidy_components,
        "subsidy_eur_per_year": subsidy * energy,
    }
