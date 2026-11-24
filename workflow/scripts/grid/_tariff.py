"""Turn one country's electricity-charges workbook into what the model prices.

The workbook (data/electricity_charges/{area}.xlsx) lists every charge on top
of the wholesale price, and its `model_variable` says what each is charged on:

* **the grid connection** — `fee_eur_per_mw_per_year`, `fee_eur_per_mwh`,
  `fee_eur_per_year` and `multiplier` rows. They exist only where the run buys
  grid power, and they become the regimes below, which the solve tries.
* **consumption** — `consumption_eur_per_mwh` rows, on every MWh the plant's
  processes draw, own renewables included (the Industriestrompreis, the
  Strompreiskompensation; they are negative because they are subsidies). The
  report applies them to every run in the area, with a grid or without, after
  the solve. The reason: the battery is currently set up so that it can charge
  and discharge at the same time, which cannot happen in reality. If the
  wholesale price were below the subsidy, the plant would earn money just by
  charging and discharging the battery at the same time, which again does not
  happen in reality. Therefore the subsidy on consumption is only applied
  afterwards. Rows sharing a `one_of` param exclude each other; the cheaper for
  the plant is kept.

  TODO: fix the battery so that it cannot charge and discharge at the same
  time, and then include the charges on consumption in the optimisation.

A grid bill depends on the plant's utilisation hours h = E / p_nom and its
annual energy E. Tariffs define h as annual import over peak import; p_nom of
the grid import stands in for the peak because every MW of it is charged, so
the optimum builds no more than its highest hourly import. Billing peaks are
15-minute averages, so the hourly model reads h slightly high. Neither h nor E
is known before the solve, and h is a ratio of two decision variables, so the
LP cannot price it directly. Instead the tariff is cut into regimes: within
each, the bill is linear,

    bill = capacity · p_nom + energy · E + fixed,

and valid only while h_min ≤ h ≤ h_max and e_min ≤ E ≤ e_max. The solve runs
once per regime with those bounds as linear constraints and keeps the cheapest.
This leads to several runs for one scenario and more computational effort.

Where the regimes come from, in the tariff sheet's terms:

* **h bands** — rate rows with `h_min`/`h_max` (the Leistungspreis and
  Arbeitspreis below and above 2500 h). Their boundaries cut [0, ∞) into base
  bands; a rate row without h bounds applies in every band.
* **tiers** — rate rows with `tier_from_mwh`/`tier_to_mwh`, where only the
  energy inside the tier pays that rate (KWKG, offshore, §19 surcharges). Within
  tier k the bill is r_k·E + Σ_{i<k} r_i·width_i − r_k·from_k, so each tier
  boundary cuts E into segments with their own energy rate and fixed offset.
  A grid run is assumed to build a connection and import at least 1 GWh a
  year through it (`ASSUMED_MIN_ANNUAL_IMPORT_MWH`), so every regime holds
  from 1 GWh: the segments below it get no regime of their own, which reduces
  the number of runs per scenario and therefore the computational effort, and
  the first GWh's higher rate is exact in the fixed term of the others. The
  case of a plant on its own renewables with no connection at all is not lost:
  it is the RES-only scenario of the same route, technology and site.
* **multipliers** — discount rows (§19 StromNEV, the French TURPE reduction)
  with their own h range and `e_min_mwh`, scaling the components named in
  `applies_to`: rates and €/yr rows alike (the TURPE reduction also lowers the
  metering and customer management fees). The workbook lists every component
  a discount scales; a group name is not enough. Each adds regimes next to the
  base bands rather than replacing them. A discount may overlap its base band
  because the cheapest regime wins; that is why multipliers above 1 are
  refused.
* **shares** — `percent_of` rows, a percentage of the components named in
  `applies_to`, as they stand in each regime after its discount (the French
  CTA: 5 % of the metering, customer management and second-feed fees after
  the TURPE reduction). A share takes the shape of what it is a share of: of a
  €/yr fee it is €/yr, of a rate it is a rate.

What stays out of the regimes:

* **connection costs** — `connection_capex_eur_per_mw` and `connection_capex_eur`
  rows, paid once when the connection is built and the same whichever regime
  bills it. They replace the config's `grid.connection_capex_eur_per_mw`; the
  solve annualises both over the config's `grid.connection_financing_period`
  (see scripts/solve/_grid_tariff.py). The part per MW is in the LP and sizes
  the connection. The fixed one-off sum is owed for building a connection at
  all, a yes/no decision an LP could only price with a binary variable (a
  MILP); as every grid run of a workbook area builds one (see the tiers
  above), it is a constant instead, which the solve adds to every regime's
  cost.

The workbook's own sheets compute the values; pandas reads the values Excel
cached on its last save, so the workbook has to be saved by Excel.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# The values of the tariff sheet's `model_variable` column: what a row is
# charged on.
CAPACITY = "fee_eur_per_mw_per_year"  # per MW of the connection's peak (Leistungspreis)
ENERGY = "fee_eur_per_mwh"  # per MWh imported (Arbeitspreis, surcharges, taxes)
FIXED = "fee_eur_per_year"  # per year, whatever the plant does (metering)
MULTIPLIER = "multiplier"  # a factor on other rows (the §19 discount, the TURPE reduction)
PERCENT_OF = "percent_of"  # a percentage of other rows, after any discount (the French CTA)
CONSUMPTION = "consumption_eur_per_mwh"  # per MWh consumed, own renewables included
CONNECTION = ("connection_capex_eur", "connection_capex_eur_per_mw")  # paid once for the connection

# The unit the `unit` column must give for each model_variable. A check, not a
# conversion: the workbook converts to these units itself, and a row in any
# other unit (ct/kWh where EUR/MWh is expected) stops the rule rather than be
# read ten times too small. A model_variable missing here is refused as unknown.
UNITS = {
    CAPACITY: "EUR/MW/yr",
    ENERGY: "EUR/MWh",
    FIXED: "EUR/yr",
    MULTIPLIER: "factor",
    PERCENT_OF: "percent",
    CONSUMPTION: "EUR/MWh",
    "connection_capex_eur": "EUR",
    "connection_capex_eur_per_mw": "EUR/MW",
}
# Columns whose meaning the regimes do not cover yet; FRA and ESP will need them.
# A workbook may have them, empty: only a row that fills one is refused, so it
# is not read as if its time window were not there.
UNSUPPORTED_COLUMNS = ("period", "offpeak_min", "offpeak_max")
# The keys the `params` column may hold (see _parse_params).
KNOWN_PARAMS = ("one_of",)

INF = math.inf
# What a grid run is assumed to import a year at least, through a connection it
# always builds. Every regime holds from here (see build_regimes).
ASSUMED_MIN_ANNUAL_IMPORT_MWH = 1_000.0


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
        """Add one component's rates to the regime.

        The totals are what the solve prices; `components` keeps the same
        numbers by component, for the report's list of charges. Both grow
        together here, so the components always add up to the totals.
        """
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
    """A bound from a workbook cell, `default` where the cell is empty.

    The caller says what empty means: every lower bound (h_min, e_min_mwh,
    tier_from_mwh) is read with 0 and every upper bound with INF, so an empty
    cell is no limit.
    """
    return default if pd.isna(value) else float(value)


def _parse_params(text) -> dict:
    """The `params` cell as a dict: "one_of=ets_or_isp" → {"one_of": "ets_or_isp"}.

    Several params are separated by ";"; an unknown key stops the rule.
    """
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
    """The `applies_to` cell as names: "leistungspreis, arbeitspreis" → ["leistungspreis", "arbeitspreis"]."""
    if pd.isna(text):
        return []
    return [name.strip() for name in str(text).split(",") if name.strip()]


def read_workbook(path: Path) -> pd.DataFrame:
    """The tariff sheet, the only one the model reads; the others compute its values."""
    return pd.read_excel(path, sheet_name="tariff").dropna(how="all")


def _select_rows(tariff: pd.DataFrame, consumer: str, source: str) -> pd.DataFrame:
    """The rows that bill this consumer, checked for what the regimes can express."""
    missing_values = tariff[tariff["value"].isna()]
    if not missing_values.empty:
        raise ValueError(
            f"{source}: rows {list(missing_values['component'])} have no cached value — "
            f"open the workbook in Excel and save it so its formulas are evaluated"
        )
    # A row with no consumer bills everyone; a row that names one bills only it.
    # For iron_steel: every unnamed row, its own Stromsteuer (0) and §19 surcharge
    # (0.25); standard's Stromsteuer and §19 surcharge are left out.
    rows = tariff[tariff["consumer"].isna() | (tariff["consumer"] == consumer)].copy()

    # Each check stops the rule with what is wrong, rather than let a row be
    # priced wrong or dropped without a word.
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
    """Whether a row's range (its lo_col..hi_col cells) holds the whole of lo..hi.

    The Leistungspreis row of 2500 h and up covers the band 2500–∞ but not
    0–2500. The same for tiers, with the tier columns.
    """
    return _bound(getattr(row, lo_col), 0.0) <= lo and _bound(getattr(row, hi_col), INF) >= hi


def _breakpoints(rows: pd.DataFrame, lo_col: str, hi_col: str) -> list[tuple[float, float]]:
    """Every boundary the rows' ranges name, as neighbouring pairs from 0 to ∞.

    DEU's h bounds 0 and 2500 give [(0, 2500), (2500, ∞)]; its tier bounds,
    [(0, 1000), (1000, ∞)]. A third band in the workbook, at 5000 h, would
    give three pairs without a change here.
    """
    points = {0.0, INF}
    points |= {float(v) for v in rows[lo_col].dropna()}
    points |= {float(v) for v in rows[hi_col].dropna()}
    ordered = sorted(points)
    return list(zip(ordered[:-1], ordered[1:]))


def _component_terms(rows: pd.DataFrame, band, segment, source: str) -> dict:
    """capacity/energy/fixed of each rate component within one h band and E segment."""
    terms = {}
    for component, comp_rows in rows.groupby("component", sort=False):
        # The component's rows that hold in this band: the Leistungspreis keeps
        # one of its two, the KWKG surcharge (no h bounds) both its tiers.
        in_band = [row for row in comp_rows.itertuples() if _covers(row, *band, "h_min", "h_max")]
        if not in_band:
            raise ValueError(
                f"{source}: {component} has no rate for {band[0]:g}–{band[1]:g} h; "
                f"its h ranges must cover every hour band"
            )
        # A component is either one flat rate in the band, or a set of tiers.
        # Two flat rates in the same band would leave it unclear which applies.
        tiered = [row for row in in_band if pd.notna(row.tier_from_mwh) or pd.notna(row.tier_to_mwh)]
        flat = [row for row in in_band if row not in tiered]
        if tiered and flat or len(flat) > 1:
            raise ValueError(f"{source}: {component} has more than one rate for {band[0]:g}–{band[1]:g} h")
        if flat:
            row = flat[0]
            key = "capacity" if row.model_variable == CAPACITY else "energy"
            terms[component] = {key: row.value}
            continue

        # The tiers, checked to run from 0 to ∞ without a gap or an overlap.
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

        # The LP can price one €/MWh per regime, but a tiered charge has several:
        # the first MWh pay one rate, the ones above it another. Within the
        # segment, every MWh falls in the same tier, so the bill can be written
        # as that tier's rate on all the MWh plus a constant, which corrects for
        # the MWh below the tier having been priced at the wrong rate:
        #
        #     bill = rate · E + (what the earlier tiers cost in full − rate · start)
        #
        # KWKG above 1 GWh: 4.46 €/MWh on the first 1000 MWh, 0.669 above, so
        # 0.669 · E + (4.46 · 1000 − 0.669 · 1000) = 0.669 · E + 3,791 €/yr,
        # which at E = 5000 MWh is 7,136 €, as 4.46 · 1000 + 0.669 · 4000 is.
        # `current` is the tier the segment lies in, `start` where it begins, and
        # `below` what every tier before it charges in full.
        current = next(row for row in tiers if _covers(row, *segment, "tier_from_mwh", "tier_to_mwh"))
        start = _bound(current.tier_from_mwh, 0.0)
        below = sum(
            row.value * (_bound(row.tier_to_mwh, INF) - _bound(row.tier_from_mwh, 0.0))
            for row in tiers if _bound(row.tier_to_mwh, INF) <= start
        )
        terms[component] = {"energy": current.value, "fixed": below - current.value * start}
    return terms


def _label(prefix: str, lo: float, hi: float) -> str:
    """A range for a regime's name: ("h", 2500, ∞) → "h2500-", ("e", 0, 1000) → "e0-1000"."""
    return f"{prefix}{lo:g}-" + ("" if hi == INF else f"{hi:g}")


def build_regimes(rows: pd.DataFrame, source: str) -> list[Regime]:
    """Every regime of the selected rows. Unconditioned €/yr rows sit in each one's fixed term.

    A base regime for each h band (and energy segment, where a tier cuts
    one), then one more for each discount and band it overlaps. DEU: h0-2500,
    h2500- and the three §19 discounts.
    """
    # Four piles by what a row is charged on. The rates make the bands and the
    # bill; the €/yr rows go into every regime's fixed term; the multipliers add
    # the discount regimes; the shares are added last, on what the others came
    # to. Consumption and connection rows are in none of them: they are not part
    # of the grid bill (see _consumption, _connection).
    rates = rows[rows["model_variable"].isin([CAPACITY, ENERGY])]
    fixed_rows = rows[rows["model_variable"] == FIXED]
    multipliers = rows[rows["model_variable"] == MULTIPLIER]
    shares = rows[rows["model_variable"] == PERCENT_OF]

    conditions = ["h_min", "h_max", "e_min_mwh", "tier_from_mwh", "tier_to_mwh"]
    conditioned_fixed = fixed_rows[fixed_rows[conditions].notna().any(axis=1)]
    if not conditioned_fixed.empty:
        raise NotImplementedError(f"{source}: conditions on €/yr rows {list(conditioned_fixed['component'])} are not modelled yet")
    # A share follows its components into every regime; it has no box of its own.
    conditioned_shares = shares[shares[conditions].notna().any(axis=1)]
    if not conditioned_shares.empty:
        raise NotImplementedError(f"{source}: conditions on percent_of rows {list(conditioned_shares['component'])} are not modelled yet")
    if rates["e_min_mwh"].notna().any():
        raise NotImplementedError(
            f"{source}: an e_min_mwh condition on a rate row ({list(rates.loc[rates['e_min_mwh'].notna(), 'component'])}) "
            f"is not modelled yet — only multipliers carry one"
        )
    # What a discount or a share may name: the components that make up the bill
    # before any share. A share of a share is not allowed, so the order in which
    # shares are added cannot change the bill.
    billed = set(rates["component"]) | set(fixed_rows["component"])
    for row in multipliers.itertuples():
        if not 0.0 <= row.value <= 1.0:
            raise ValueError(
                f"{source}: {row.component} multiplies by {row.value}; only discounts (0–1) are allowed, "
                f"because the cheapest regime is the one kept"
            )
    for row in pd.concat([multipliers, shares]).itertuples():
        names = _applies_to(row.applies_to)
        if not names:
            raise ValueError(f"{source}: {row.component} has no applies_to")
        unknown = set(names) - billed
        if unknown:
            raise ValueError(f"{source}: {row.component} applies to unknown components {sorted(unknown)}")

    bands = _breakpoints(rates, "h_min", "h_max")

    # Each energy segment is priced at its own tier rates (`segment`) and holds
    # for a range of annual import (`window`, the regime's e_min..e_max box).
    # The 1 GWh boundary in DEU is a tier: the tiered surcharges charge the
    # first GWh a year at a higher rate than every GWh after it. We assume a
    # grid run builds a connection and imports at least
    # ASSUMED_MIN_ANNUAL_IMPORT_MWH (1 GWh) a year through it, so no window
    # starts below it:
    # - the segments below 1 GWh get no regimes of their own; each would be one
    #   more solve that a plant consuming TWh never wins;
    # - the first segment above it holds from 1 GWh, and the first GWh's higher
    #   rate sits exactly in its fixed term, as every regime's import is at
    #   least that;
    # - so every grid run builds a connection, which makes the connection's
    #   fixed costs (metering, the fixed one-off sum) certain: constants the
    #   solve adds to every regime's cost rather than yes/no decisions.
    # A plant on its own renewables with no connection at all is not lost: it
    # is the RES-only scenario of the same route, technology and site, which
    # has no grid. A grid run that would be cheaper without a connection shows
    # as dearer than that scenario.
    segments = [s for s in _breakpoints(rates, "tier_from_mwh", "tier_to_mwh")
                if s[1] > ASSUMED_MIN_ANNUAL_IMPORT_MWH]
    windows = [(max(lo, ASSUMED_MIN_ANNUAL_IMPORT_MWH), hi) for lo, hi in segments]

    # One regime: its box (h_lo..h_hi, e_lo..e_hi), every €/yr row, and each
    # rate component's terms in its band and segment, scaled by `scale` where a
    # discount applies to the component (1, no discount, where it does not).
    # Then the shares, each a percentage of its components as they stand now,
    # after the discount: the CTA is 5 % of the TURPE fees the plant actually
    # pays, so in the TURPE reduction's regime it is 5 % of the reduced fees.
    def regime(name, h_lo, h_hi, e_lo, e_hi, band, segment, scale=None):
        scale = scale or {}
        piece = Regime(name, h_lo, h_hi, e_lo, e_hi)
        for row in fixed_rows.itertuples():
            piece.add(row.component, fixed=row.value * scale.get(row.component, 1.0))
        for component, parts in _component_terms(rates, band, segment, source).items():
            factor = scale.get(component, 1.0)
            piece.add(
                component,
                capacity=parts.get("capacity", 0.0) * factor,
                energy=parts.get("energy", 0.0) * factor,
                fixed=parts.get("fixed", 0.0) * factor,
            )
        # A share names no other share (see `billed`), so the order they are
        # added in does not matter.
        for row in shares.itertuples():
            of = [piece.components[name] for name in _applies_to(row.applies_to)]
            share = row.value / 100.0
            piece.add(
                row.component,
                capacity=share * sum(parts["capacity_eur_per_mw_per_year"] for parts in of),
                energy=share * sum(parts["energy_eur_per_mwh"] for parts in of),
                fixed=share * sum(parts["fixed_eur_per_year"] for parts in of),
            )
        return piece

    # The base regimes: every band with every segment. A tariff with one segment
    # names its regimes by their band alone.
    many_segments = len(segments) > 1
    regimes = []
    for band in bands:
        for segment, window in zip(segments, windows):
            name = _label("h", *band) + ("_" + _label("e", *window) if many_segments else "")
            regimes.append(regime(name, *band, *window, band, segment))

    # The discount regimes, next to the base ones rather than in their place. A
    # discount holds where its own range overlaps a band: §19's 7000–7500 h
    # lies in the band 2500–∞, whose rates it scales, and misses 0–2500. Its
    # energy condition (≥ 10 GWh) narrows the segment's window the same way. A
    # plant inside the discount's box also fits the base regime of its band;
    # both are solved and the cheaper, the discount, wins. That is why only
    # factors up to 1 are allowed: one above would never win and never apply.
    for row in multipliers.itertuples():
        lo, hi = _bound(row.h_min, 0.0), _bound(row.h_max, INF)
        e_min = _bound(row.e_min_mwh, 0.0)
        scale = {name: row.value for name in _applies_to(row.applies_to)}
        for band in bands:
            h_lo, h_hi = max(lo, band[0]), min(hi, band[1])
            if h_lo >= h_hi:
                continue
            for segment, window in zip(segments, windows):
                e_lo = max(e_min, window[0])
                if e_lo >= window[1]:
                    continue
                name = f"{row.component}_{_label('h', h_lo, h_hi)}" + (
                    "_" + _label("e", e_lo, window[1]) if many_segments or e_lo > 0 else "")
                regimes.append(regime(name, h_lo, h_hi, e_lo, window[1], band, segment, scale))

    return regimes


def _connection(rows: pd.DataFrame) -> dict:
    """The one-off cost of building the connection: per MW, and in one sum.

    Each in total and by component, e.g. the Netzanschlusskosten with a part per
    MW and a fixed part, and the Baukostenzuschuss per MW.
    """
    components = {}
    for row in rows[rows["model_variable"].isin(CONNECTION)].itertuples():
        entry = components.setdefault(row.component, {"capex_eur_per_mw": 0.0, "capex_eur": 0.0})
        key = "capex_eur" if row.model_variable == "connection_capex_eur" else "capex_eur_per_mw"
        entry[key] += row.value
    return {
        "capex_eur_per_mw": _num(sum(entry["capex_eur_per_mw"] for entry in components.values())),
        "capex_eur": _num(sum(entry["capex_eur"] for entry in components.values())),
        "components": {name: {key: _num(value) for key, value in entry.items()}
                       for name, entry in components.items()},
    }


def _consumer_charges(tariff: pd.DataFrame, consumer: str, source: str) -> dict:
    """One consumer's grid regimes, consumption charge and connection costs."""
    rows = _select_rows(tariff, consumer, source)
    rate, components = _consumption(rows)
    return {
        "grid": {"regimes": [piece.to_dict() for piece in build_regimes(rows, source)]},
        "consumption": {"eur_per_mwh": _num(rate),
                        "components": {name: _num(value) for name, value in components.items()}},
        "connection": _connection(rows),
    }


def build_charges(path: Path, area: str) -> dict:
    """One area's charges as the dict written to resources/electricity_charges/{area}.yaml.

    One entry under `consumers` per consumer the tariff sheet's `consumer`
    column names (DEU: iron_steel and standard), each with what a plant charged as
    that consumer pays: `grid`, the regimes the solve tries where the run buys
    grid power; `consumption`, the charge the report puts on every MWh the
    plant consumes; `connection`, the one-off cost of building the grid
    connection. A row with no consumer bills every one. Which consumer a run is
    charged as is the config's `electricity_charges.consumer`, not the
    workbook's: iron_steel, unless a run asks for another (see `for_consumer`).
    """
    source = str(path)
    tariff = read_workbook(path)
    # Every consumer the workbook names, each priced in full; the run picks one.
    consumers = sorted(tariff["consumer"].dropna().astype(str).str.strip().unique())
    if not consumers:
        raise ValueError(f"{source}: the tariff sheet's consumer column names no consumer")

    rows = tariff[tariff["model_variable"].isin(UNITS)]
    # What each component is, for the report's list of charges: the workbook's
    # own `group` column, carried through unread by the regimes.
    groups = {row.component: (None if pd.isna(row.group) else row.group) for row in rows.itertuples()}

    return {
        "area": area,
        "source": source.replace("\\", "/"),
        "component_groups": groups,
        "consumers": {consumer: _consumer_charges(tariff, consumer, source) for consumer in consumers},
    }


def default_charges(area: str) -> dict:
    """An area without a workbook: its grid and its connection are priced from
    the config's `grid` block, and nothing is charged on consumption, whichever
    consumer the run is charged as."""
    return {"area": area, "source": "config", "component_groups": {}, "consumers": None}


def for_consumer(charges: dict | None, consumer: str) -> dict | None:
    """An area's charges as one consumer pays them, as the solve and the report read them.

    `charges` is one resources/electricity_charges/{area}.yaml; `consumer` the
    config's `electricity_charges.consumer`. An area without a workbook charges
    every consumer the same: the config's grid fees and nothing on consumption.
    """
    if charges is None:
        return None
    flat = {"area": charges["area"], "source": charges["source"], "consumer": consumer,
            "component_groups": charges.get("component_groups") or {}}
    if charges.get("consumers") is None:
        return {**flat, "grid": {"regimes": None},
                "consumption": {"eur_per_mwh": 0.0, "components": {}}, "connection": None}
    if consumer not in charges["consumers"]:
        raise ValueError(f"{charges['source']}: no charges for consumer {consumer!r}; "
                         f"the workbook names {sorted(charges['consumers'])}")
    return {**flat, **charges["consumers"][consumer]}
