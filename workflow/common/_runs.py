"""What a run is, and how the scenario rows expand into runs.

A scenario name is an umbrella label and carries no meaning for the workflow.
What identifies one network is the run key — scenario, area, date range, route —
and every part of it but the scenario comes from the rows filed under that name.

Two cells expand rather than naming one value: `route` accepts `all-routes` or a
`|`-separated list, and `area` accepts `all-areas`. Both are resolved here, at
DAG-build time, so neither ever reaches a wildcard or a filename.

Deliberately free of heavy imports: the Snakefile reads this on every DAG build.
"""

import pandas as pd

ROUTES = ("h2-only", "h2-dri-eaf", "ng-dri-eaf", "mix-dri-eaf", "moe-eaf", "ew-eaf",
          "moe-eaf-export", "ew-eaf-export")

# Reserved: an `-export` route ships its iron to the EU and melts it there. The
# DRI routes cannot yet, because sponge iron has to be briquetted before it will
# survive the trip and the briquetting step is not built. The two routes above
# make cold solid iron to begin with and needed no such step.
PLANNED_ROUTES = ("h2-dri-eaf-export", "ng-dri-eaf-export", "mix-dri-eaf-export")

# How the iron reaches the furnace, and which step made it. The first decides
# how much electricity the EAF needs — melting from cold is most of an EAF's
# bill, and iron that arrives hot or liquid has already been paid for upstream.
# The second decides how much iron a t of steel takes, which is gangue: DR
# pellets bring some, and the electrolytic routes almost none.
#
# An export route always charges cold. Its iron went across an ocean.
EAF_CHARGE = {
    "h2-dri-eaf":     ("hot", "dri-h2"),
    "ng-dri-eaf":     ("hot", "dri-ng"),
    "mix-dri-eaf":    ("hot", "dri-h2"),    # both shafts run on DR pellets
    "moe-eaf":        ("liquid", "moe"),
    "ew-eaf":         ("cold", "ew"),
    "moe-eaf-export": ("cold", "moe"),
    "ew-eaf-export":  ("cold", "ew"),
}

ALL_ROUTES = "all-routes"
ALL_AREAS = "all-areas"

# Several routes in one cell are separated by this.
ROUTE_SEPARATOR = "|"

# Longest first: `moe-eaf` is a prefix of `moe-eaf-export`, and a regex
# alternation takes the first branch that matches.
ROUTE_PATTERN = "|".join(sorted(ROUTES, key=len, reverse=True))

# The columns that, with the scenario, identify one network.
RUN_KEY = ["scenario", "area", "start_date", "end_date"]


def expand_route_cell(cell: str) -> list[str]:
    """One `route` cell -> the concrete route ids it names."""
    if cell == ALL_ROUTES:
        return list(ROUTES)
    routes = cell.split(ROUTE_SEPARATOR)
    unknown = set(routes) - set(ROUTES)
    if unknown:
        planned = unknown & set(PLANNED_ROUTES)
        hint = " (reserved but not implemented yet)" if planned else ""
        raise ValueError(
            f"config/scenarios.csv names unknown route(s) {sorted(unknown)}{hint} — "
            f"expected {ALL_ROUTES!r} or one of {ROUTES}"
        )
    return routes


def area_aliases(areas: dict) -> dict:
    """Every code that names an area -> the area. Includes each market code.

    A bidding zone or NEM region is how a market knows an area, so writing
    DE_LU where DEU is meant resolves rather than failing.
    """
    aliases = {area: area for area in areas}
    for area, cfg in areas.items():
        market_area = cfg.get("market_area")
        if market_area:
            aliases.setdefault(market_area, area)
    return aliases


def zone_parents(areas: dict) -> dict:
    """Zone -> the area that contains it. An area that is nobody's zone is its own.

    Reporting ranks a country's zones against each other, so it needs to know
    which runs are alternative sites for the same place rather than different
    places.
    """
    parents = {area: area for area in areas}
    for area, cfg in areas.items():
        for zone in cfg.get("zones", []):
            parents[zone] = area
    return parents


def top_level_areas(areas: dict) -> list[str]:
    """The areas `all-areas` starts from — every area that is not someone's zone."""
    nested = {zone for cfg in areas.values() for zone in cfg.get("zones", [])}
    return sorted(set(areas) - nested)


def resolve_market_areas(area: str, areas: dict) -> list[str]:
    """`area` if it trades in a market itself, else the zones through which it does."""
    if areas[area].get("market"):
        return [area]
    return sorted(areas[area].get("zones", []))


def scenario_areas(scenario: str, rows: pd.DataFrame, areas: dict) -> list[str]:
    """The areas one scenario covers, resolving aliases, `all-areas` and zones.

    `all-areas` starts from every top-level area. A scenario with a `tech=grid`
    row then needs prices, so an area that has no market of its own is replaced
    by the zones through which it trades — that is how naming Australia in a
    grid scenario yields its NEM regions separately, while an islanded scenario
    keeps the whole country. An area with neither drops out.

    Resolving per scenario rather than per row is what keeps a scenario's
    renewables and its prices in the same place.
    """
    aliases = area_aliases(areas)
    written = set(rows["area"]) - {ALL_AREAS}
    unknown = written - set(aliases)
    if unknown:
        raise ValueError(
            f"scenario {scenario!r} names unknown area(s) {sorted(unknown)} — "
            f"config.yaml's `areas` registry accepts {sorted(aliases)}"
        )
    named = {aliases[a] for a in written}

    if ALL_AREAS in set(rows["area"]):
        if named:
            raise ValueError(
                f"scenario {scenario!r} mixes {ALL_AREAS!r} with named areas "
                f"{sorted(named)}; use one or the other."
            )
        named = set(top_level_areas(areas))

    if not (rows["tech"] == "grid").any():
        return sorted(named)

    resolved = sorted({z for area in named for z in resolve_market_areas(area, areas)})
    if not resolved:
        raise ValueError(
            f"scenario {scenario!r} has a tech=grid row, but none of {sorted(named)} "
            f"trades in a wholesale market or names zones that do. Model it islanded, "
            f"or give the area a `market` or `zones` entry in config.yaml."
        )
    return resolved


def build_runs_frame(scenarios: pd.DataFrame, areas: dict) -> pd.DataFrame:
    """(scenario, area, start_date, end_date, route) — one row per network solved."""
    frames = []
    for scenario, rows in scenarios.groupby("scenario"):
        resolved = scenario_areas(scenario, rows, areas)
        dates = rows[["start_date", "end_date"]].drop_duplicates()
        routes = sorted({r for cell in rows["route"] for r in expand_route_cell(cell)})
        frames.append(
            dates.merge(pd.Series(resolved, name="area"), how="cross")
            .merge(pd.Series(routes, name="route"), how="cross")
            .assign(scenario=scenario)
        )
    runs = pd.concat(frames)[RUN_KEY + ["route"]]
    return runs.sort_values(RUN_KEY + ["route"]).reset_index(drop=True)


def expand_scenario_rows(scenarios: pd.DataFrame, areas: dict) -> pd.DataFrame:
    """The input rows with `all-areas` resolved, so each row names one real area.

    This is what the solve rule's `lookup` reads: it joins on a concrete area, so
    the expansion has to have happened before the DAG is built.
    """
    frames = []
    for scenario, rows in scenarios.groupby("scenario"):
        resolved = scenario_areas(scenario, rows, areas)
        # Every row is re-pointed at the resolved areas, not just the `all-areas`
        # ones: naming AUS in a grid scenario means its zones for the CF rows too,
        # or the renewables would sit in a different place from the prices.
        frames.append(
            rows.drop(columns="area")
            .merge(pd.Series(resolved, name="area"), how="cross")[scenarios.columns]
        )
    return pd.concat(frames).reset_index(drop=True)


def load_scenarios(path, areas: dict) -> pd.DataFrame:
    """Read the scenario table into validated rows, one per input timeseries.

    `#` rows are planned scenarios parked next to their live siblings (see the
    export routes) and are skipped.
    """
    scenarios = pd.read_csv(
        path, comment="#", dtype={"start_date": str, "end_date": str}
    )
    expanded = expand_scenario_rows(scenarios, areas)
    check_one_series_per_tech(expanded)
    check_run_coverage(expanded)
    return expanded


def check_one_series_per_tech(scenarios: pd.DataFrame) -> None:
    """Reject a run that would be handed two timeseries for the same tech.

    One network is one (scenario, area, date range) group, so a tech may appear
    once per group — twice and the solve gets two series for one tech.
    """
    key = RUN_KEY + ["tech"]
    if scenarios.duplicated(subset=key).any():
        offenders = scenarios[scenarios.duplicated(subset=key, keep=False)][key]
        raise ValueError(
            f"the scenario table has duplicate rows for one run and tech:\n{offenders}"
        )


def check_run_coverage(scenarios: pd.DataFrame) -> None:
    """Reject a scenario that pairs renewables with prices from somewhere else.

    Rows join into a network by run key, so a grid row filed under a different
    area or date range from the CF rows never meets them: it splits off a
    generation-free network of its own and leaves the renewables unpriced. Both
    halves are well-formed on their own, which is what makes the mistake quiet —
    the signature is a scenario holding a CF-only run *and* a grid-only run at
    the same time.

    A scenario made only of grid rows (a market-supplied benchmark) is fine, and
    so is one that runs some areas islanded and others grid-connected.
    """
    for scenario, rows in scenarios.groupby("scenario"):
        cf_runs = set(map(tuple, rows.loc[rows["tech"] != "grid", RUN_KEY].values))
        grid_runs = set(map(tuple, rows.loc[rows["tech"] == "grid", RUN_KEY].values))
        unpriced = cf_runs - grid_runs
        generation_free = grid_runs - cf_runs
        if unpriced and generation_free:
            raise ValueError(
                f"scenario {scenario!r} has renewables and grid prices that never meet:\n"
                f"  renewables with no price: {sorted(r[1:] for r in unpriced)}\n"
                f"  price with no renewables: {sorted(r[1:] for r in generation_free)}\n"
                f"Rows join by (area, start_date, end_date). Site the renewables in the "
                f"area whose prices they pay, or split the two into separate scenarios."
            )
