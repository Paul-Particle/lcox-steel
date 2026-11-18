"""Turn one country's electricity-charges workbook into what the model prices.

The workbook (data/electricity_charges/{area}.xlsx) lists every charge on top
of the wholesale price, and its `model_variable` says what each is charged on:

* **the grid connection** — `fee_eur_per_mw_per_year`, `fee_eur_per_mwh`,
  `fee_eur_per_year` and `multiplier` rows. They exist only where the run buys
  grid power, and they become the regimes below, which the solve tries.
* **consumption** — `consumption_eur_per_mwh` rows, on every MWh the plant's
  processes draw, own renewables included (the Industriestrompreis, the
  Strompreiskompensation; negative is a subsidy, positive a levy). The report
  applies them to every run in the area, with a grid or without, after the
  solve: the consumption is set by the steel the plant makes, so they move no
  decision, and as a negative €/MWh on the grid import in the LP they would
  pay the plant to import energy and burn it in storage losses. Rows sharing a
  `one_of` param exclude each other; the cheaper for the plant is kept.

A grid bill depends on the plant's utilisation hours h = E / p_nom and its
annual energy E. Tariffs define h as annual import over peak import; p_nom of
the grid import stands in for the peak because every MW of it is charged, so
the optimum builds no more than its highest hourly import. Billing peaks are
15-minute averages, so the hourly model reads h slightly high. Neither is known before the solve, and h is a ratio of two
decision variables, so the LP cannot price it directly. Instead the tariff is
cut into regimes: within each, the bill is linear,

    bill = capacity · p_nom + energy · E + fixed,

and valid only while h_min ≤ h ≤ h_max and e_min ≤ E ≤ e_max. The solve runs
once per regime with those bounds as linear constraints and keeps the cheapest.

Where the regimes come from, in the tariff sheet's terms:

* **h bands** — rate rows with `h_min`/`h_max` (the Leistungspreis and
  Arbeitspreis below and above 2500 h). Their boundaries cut [0, ∞) into base
  bands; a rate row without h bounds applies in every band.
* **tiers** — rate rows with `tier_from_mwh`/`tier_to_mwh`, where only the
  energy inside the tier pays that rate (KWKG, offshore, §19 surcharges). Within
  tier k the bill is r_k·E + Σ_{i<k} r_i·width_i − r_k·from_k, so each tier
  boundary cuts E into segments with their own energy rate and fixed offset.
* **multipliers** — discount rows (§19 StromNEV) with their own h range and
  `e_min_mwh`, scaling the components named in `applies_to`. Each adds regimes
  next to the base bands rather than replacing them. A discount may overlap its
  base band because the cheapest regime wins; that is why multipliers above 1
  are refused.

What stays out of the regimes:

* **connection costs** (`connection_capex_eur`, `connection_capex_eur_per_mw`)
  are not read yet. The config's `grid.connection_*` still prices the connection.

The workbook's own sheets compute the values; pandas reads the values Excel
cached on its last save, so the workbook has to be saved by Excel.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

CAPACITY = "fee_eur_per_mw_per_year"
ENERGY = "fee_eur_per_mwh"
FIXED = "fee_eur_per_year"
MULTIPLIER = "multiplier"
CONSUMPTION = "consumption_eur_per_mwh"
CONNECTION = ("connection_capex_eur", "connection_capex_eur_per_mw")

UNITS = {
    CAPACITY: "EUR/MW/yr",
    ENERGY: "EUR/MWh",
    FIXED: "EUR/yr",
    MULTIPLIER: "factor",
    CONSUMPTION: "EUR/MWh",
    "connection_capex_eur": "EUR",
    "connection_capex_eur_per_mw": "EUR/MW",
}
# Columns whose meaning the regimes do not cover yet; FRA and ESP will need them.
UNSUPPORTED_COLUMNS = ("period", "offpeak_min", "offpeak_max")
KNOWN_PARAMS = ("one_of",)

INF = math.inf


@dataclass
class Regime:
    """One linear piece of the bill and the (h, E) box it holds in."""

    name: str
    h_min: float
    h_max: float
    e_min_mwh: float
    e_max_mwh: float
    capacity_eur_per_mw_per_year: float = 0.0
    energy_eur_per_mwh: float = 0.0
    fixed_eur_per_year: float = 0.0
    components: dict = field(default_factory=dict)

    def add(self, component: str, capacity=0.0, energy=0.0, fixed=0.0) -> None:
        self.capacity_eur_per_mw_per_year += capacity
        self.energy_eur_per_mwh += energy
        self.fixed_eur_per_year += fixed
        entry = self.components.setdefault(
            component,
            {"capacity_eur_per_mw_per_year": 0.0, "energy_eur_per_mwh": 0.0, "fixed_eur_per_year": 0.0},
        )
        entry["capacity_eur_per_mw_per_year"] += capacity
        entry["energy_eur_per_mwh"] += energy
        entry["fixed_eur_per_year"] += fixed

    def to_dict(self) -> dict:
        # Drop components that add nothing, so the YAML shows what is billed.
        components = {
            name: {key: _num(value) for key, value in parts.items() if value != 0.0}
            for name, parts in self.components.items()
            if any(value != 0.0 for value in parts.values())
        }
        return {
            "name": self.name,
            "h_min": _num(self.h_min),
            "h_max": _num(self.h_max),
            "e_min_mwh": _num(self.e_min_mwh),
            "e_max_mwh": _num(self.e_max_mwh),
            "capacity_eur_per_mw_per_year": _num(self.capacity_eur_per_mw_per_year),
            "energy_eur_per_mwh": _num(self.energy_eur_per_mwh),
            "fixed_eur_per_year": _num(self.fixed_eur_per_year),
            "components": components,
        }


def _num(value: float):
    """A float for the YAML, with an open bound as null."""
    if value == INF:
        return None
    value = float(value)
    return int(value) if value.is_integer() else round(value, 6)


def _bound(value, default: float) -> float:
    return default if pd.isna(value) else float(value)


def _parse_params(text) -> dict:
    if pd.isna(text) or not str(text).strip():
        return {}
    params = {}
    for item in str(text).split(";"):
        key, _, value = item.partition("=")
        key, value = key.strip(), value.strip()
        if key not in KNOWN_PARAMS:
            raise ValueError(f"unknown param {key!r} in {text!r}; known: {KNOWN_PARAMS}")
        params[key] = value
    return params


def _applies_to(text) -> list[str]:
    if pd.isna(text):
        return []
    return [name.strip() for name in str(text).split(",") if name.strip()]


def read_workbook(path: Path) -> tuple[pd.DataFrame, dict]:
    """The tariff sheet and the settings sheet as name → value."""
    sheets = pd.read_excel(path, sheet_name=["tariff", "settings"])
    tariff = sheets["tariff"].dropna(how="all")
    settings = dict(zip(sheets["settings"]["name"], sheets["settings"]["value"]))
    return tariff, settings


def _select_rows(tariff: pd.DataFrame, consumer: str, source: str) -> pd.DataFrame:
    """The rows that bill this consumer, checked for what the regimes can express."""
    missing_values = tariff[tariff["value"].isna()]
    if not missing_values.empty:
        raise ValueError(
            f"{source}: rows {list(missing_values['component'])} have no cached value — "
            f"open the workbook in Excel and save it so its formulas are evaluated"
        )
    rows = tariff[tariff["consumer"].isna() | (tariff["consumer"] == consumer)].copy()

    years = rows["year"].dropna().unique()
    if len(years) > 1:
        raise ValueError(f"{source}: one tariff year per workbook, found {sorted(years)}")

    unknown = set(rows["model_variable"]) - set(UNITS)
    if unknown:
        raise ValueError(f"{source}: unknown model_variable {sorted(unknown)}; known: {sorted(UNITS)}")
    for row in rows.itertuples():
        if row.unit != UNITS[row.model_variable]:
            raise ValueError(
                f"{source}: {row.component} is a {row.model_variable} in {row.unit!r}, "
                f"expected {UNITS[row.model_variable]!r}"
            )
    for column in UNSUPPORTED_COLUMNS:
        if column in rows and rows[column].notna().any():
            used = list(rows.loc[rows[column].notna(), "component"])
            raise NotImplementedError(f"{source}: column {column!r} (used by {used}) is not modelled yet")

    rows["params"] = rows["params"].map(_parse_params) if "params" in rows else [{}] * len(rows)
    return rows


def _consumption(rows: pd.DataFrame) -> tuple[float, dict[str, float]]:
    """The €/MWh charged on consumption, and the rate of each component it came from."""
    charged = rows[rows["model_variable"] == CONSUMPTION]
    used = {}
    exclusive = {}
    for row in charged.itertuples():
        group = row.params.get("one_of")
        if group is None:
            used[row.component] = used.get(row.component, 0.0) + row.value
        else:
            exclusive.setdefault(group, []).append(row)
    for group_rows in exclusive.values():
        # The rows exclude each other and apply to the same MWh: the plant takes
        # whichever is cheaper for it, the most negative.
        best = min(group_rows, key=lambda row: row.value)
        used[best.component] = used.get(best.component, 0.0) + best.value
    return sum(used.values()), used


def _covers(row, lo: float, hi: float, lo_col: str, hi_col: str) -> bool:
    return _bound(getattr(row, lo_col), 0.0) <= lo and _bound(getattr(row, hi_col), INF) >= hi


def _breakpoints(rows: pd.DataFrame, lo_col: str, hi_col: str) -> list[tuple[float, float]]:
    points = {0.0, INF}
    points |= {float(v) for v in rows[lo_col].dropna()}
    points |= {float(v) for v in rows[hi_col].dropna()}
    ordered = sorted(points)
    return list(zip(ordered[:-1], ordered[1:]))


def _component_terms(rows: pd.DataFrame, band, segment, source: str) -> dict:
    """capacity/energy/fixed of each rate component within one h band and E segment."""
    terms = {}
    for component, comp_rows in rows.groupby("component", sort=False):
        in_band = [row for row in comp_rows.itertuples() if _covers(row, *band, "h_min", "h_max")]
        if not in_band:
            raise ValueError(
                f"{source}: {component} has no rate for {band[0]:g}–{band[1]:g} h; "
                f"its h ranges must cover every hour band"
            )
        tiered = [row for row in in_band if pd.notna(row.tier_from_mwh) or pd.notna(row.tier_to_mwh)]
        flat = [row for row in in_band if row not in tiered]
        if tiered and flat or len(flat) > 1:
            raise ValueError(f"{source}: {component} has more than one rate for {band[0]:g}–{band[1]:g} h")
        if flat:
            row = flat[0]
            key = "capacity" if row.model_variable == CAPACITY else "energy"
            terms[component] = {key: row.value}
            continue

        if any(row.model_variable != ENERGY for row in tiered):
            raise ValueError(f"{source}: {component} is tiered, which only an energy rate can be")
        tiers = sorted(tiered, key=lambda row: _bound(row.tier_from_mwh, 0.0))
        edge = 0.0
        for row in tiers:
            if _bound(row.tier_from_mwh, 0.0) != edge:
                raise ValueError(f"{source}: {component}'s tiers leave a gap or overlap at {edge:g} MWh")
            edge = _bound(row.tier_to_mwh, INF)
        if edge != INF:
            raise ValueError(f"{source}: {component}'s last tier must be open-ended")

        current = next(row for row in tiers if _covers(row, *segment, "tier_from_mwh", "tier_to_mwh"))
        start = _bound(current.tier_from_mwh, 0.0)
        below = sum(
            row.value * (_bound(row.tier_to_mwh, INF) - _bound(row.tier_from_mwh, 0.0))
            for row in tiers if _bound(row.tier_to_mwh, INF) <= start
        )
        terms[component] = {"energy": current.value, "fixed": below - current.value * start}
    return terms


def _label(prefix: str, lo: float, hi: float) -> str:
    return f"{prefix}{lo:g}-" + ("" if hi == INF else f"{hi:g}")


def build_regimes(rows: pd.DataFrame, source: str) -> list[Regime]:
    """Every regime of the selected rows. Unconditioned €/yr rows sit in each one's fixed term."""
    rates = rows[rows["model_variable"].isin([CAPACITY, ENERGY])]
    fixed_rows = rows[rows["model_variable"] == FIXED]
    multipliers = rows[rows["model_variable"] == MULTIPLIER]

    conditioned_fixed = fixed_rows[fixed_rows[["h_min", "h_max", "e_min_mwh", "tier_from_mwh", "tier_to_mwh"]].notna().any(axis=1)]
    if not conditioned_fixed.empty:
        raise NotImplementedError(f"{source}: conditions on €/yr rows {list(conditioned_fixed['component'])} are not modelled yet")
    if rates["e_min_mwh"].notna().any():
        raise NotImplementedError(
            f"{source}: an e_min_mwh condition on a rate row ({list(rates.loc[rates['e_min_mwh'].notna(), 'component'])}) "
            f"is not modelled yet — only multipliers carry one"
        )
    for row in multipliers.itertuples():
        if not 0.0 <= row.value <= 1.0:
            raise ValueError(
                f"{source}: {row.component} multiplies by {row.value}; only discounts (0–1) are allowed, "
                f"because the cheapest regime is the one kept"
            )
        names = _applies_to(row.applies_to)
        if not names:
            raise ValueError(f"{source}: {row.component} has no applies_to")
        unknown = set(names) - set(rates["component"])
        if unknown:
            raise ValueError(f"{source}: {row.component} applies to unknown components {sorted(unknown)}")

    bands = _breakpoints(rates, "h_min", "h_max")
    segments = _breakpoints(rates, "tier_from_mwh", "tier_to_mwh")

    def regime(name, h_lo, h_hi, e_lo, e_hi, band, segment, scale=None):
        piece = Regime(name, h_lo, h_hi, e_lo, e_hi)
        for row in fixed_rows.itertuples():
            piece.add(row.component, fixed=row.value)
        for component, parts in _component_terms(rates, band, segment, source).items():
            factor = (scale or {}).get(component, 1.0)
            piece.add(
                component,
                capacity=parts.get("capacity", 0.0) * factor,
                energy=parts.get("energy", 0.0) * factor,
                fixed=parts.get("fixed", 0.0) * factor,
            )
        return piece

    many_segments = len(segments) > 1
    regimes = []
    for band in bands:
        for segment in segments:
            name = _label("h", *band) + ("_" + _label("e", *segment) if many_segments else "")
            regimes.append(regime(name, *band, *segment, band, segment))

    for row in multipliers.itertuples():
        lo, hi = _bound(row.h_min, 0.0), _bound(row.h_max, INF)
        e_min = _bound(row.e_min_mwh, 0.0)
        scale = {name: row.value for name in _applies_to(row.applies_to)}
        for band in bands:
            h_lo, h_hi = max(lo, band[0]), min(hi, band[1])
            if h_lo >= h_hi:
                continue
            for segment in segments:
                e_lo = max(e_min, segment[0])
                if e_lo >= segment[1]:
                    continue
                name = f"{row.component}_{_label('h', h_lo, h_hi)}" + (
                    "_" + _label("e", e_lo, segment[1]) if many_segments or e_lo > 0 else "")
                regimes.append(regime(name, h_lo, h_hi, e_lo, segment[1], band, segment, scale))

    return regimes


def _consumer_charges(tariff: pd.DataFrame, consumer: str, source: str) -> tuple[list, dict]:
    """One consumer's grid regimes and consumption charge."""
    rows = _select_rows(tariff, consumer, source)
    rate, components = _consumption(rows)
    consumption = {"eur_per_mwh": _num(rate),
                   "components": {name: _num(value) for name, value in components.items()}}
    return [piece.to_dict() for piece in build_regimes(rows, source)], consumption


def build_charges(path: Path, area: str) -> dict:
    """One area's charges as the dict written to resources/electricity_charges/{area}.yaml.

    `grid` holds the regimes the solve tries where the run buys grid power;
    `consumption` the charge the report puts on every MWh the plant consumes.
    Both are the settings sheet's `consumer`'s. Each also carries `reference`,
    the same for the sheet's `reference_consumer` if it names one: the report
    bills the plant's solved profile on it too, to show what the plant's
    privileges are worth against a consumer without them.
    """
    source = str(path)
    tariff, settings = read_workbook(path)
    consumer = settings.get("consumer")
    if not consumer:
        raise ValueError(f"{source}: the settings sheet names no consumer")
    regimes, consumption = _consumer_charges(tariff, consumer, source)

    reference_consumer = settings.get("reference_consumer")
    if pd.isna(reference_consumer) or not str(reference_consumer).strip():
        reference_consumer, reference_regimes, reference_consumption = None, None, None
    else:
        reference_consumer = str(reference_consumer).strip()
        reference_regimes, reference_consumption = _consumer_charges(tariff, reference_consumer, source)

    rows = tariff[tariff["model_variable"].isin(UNITS)]
    ignored = sorted(set(rows.loc[rows["model_variable"].isin(CONNECTION), "component"]))
    years = rows["year"].dropna().unique()
    # What each component is, for the report's list of charges: the workbook's
    # own `group` column, carried through unread by the regimes.
    groups = {
        row.component: (None if pd.isna(row.group) else row.group)
        for row in rows.itertuples() if row.model_variable not in CONNECTION
    }

    return {
        "area": area,
        "source": source.replace("\\", "/"),
        "year": int(years[0]) if len(years) else None,
        "consumer": consumer,
        "reference_consumer": reference_consumer,
        "component_groups": groups,
        "grid": {"regimes": regimes, "reference": reference_regimes},
        "consumption": {**consumption, "reference": reference_consumption},
        "not_read_yet": ignored,
    }


def default_charges(area: str) -> dict:
    """An area without a workbook: its grid is priced from the config's `grid`
    block, and nothing is charged on consumption."""
    return {
        "area": area, "source": "config", "year": None,
        "consumer": None, "reference_consumer": None, "component_groups": {},
        "grid": {"regimes": None, "reference": None},
        "consumption": {"eur_per_mwh": 0, "components": {}, "reference": None},
        "not_read_yet": [],
    }
