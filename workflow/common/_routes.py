"""The steelmaking-route registry, shared by the workflow and the network builder.

One place defines which routes exist, so `all-routes` in config/scenarios.csv and
the `route` wildcard can never drift apart. Deliberately free of heavy imports —
the Snakefile reads this on every DAG build.

Route ids are hyphenated, matching the `route` column in config/scenarios.csv and
the process-block names in config/assumptions.yaml.
"""

import pandas as pd

ROUTES = ("h2-only", "h2-dri-eaf", "ng-dri-eaf", "mix-dri-eaf", "moe-eaf", "ew-eaf")

# Reserved for the trade scenarios (produce abroad, ship the product to the EU).
# Not implemented: `*-export` appears in config/scenarios.csv as commented rows
# only, and is excluded from ALL_ROUTES until a route actually builds a network.
PLANNED_ROUTES = ("h2-dri-eaf-export", "ng-dri-eaf-export", "mix-dri-eaf-export",
                  "moe-eaf-export", "ew-eaf-export")

# Value in the CSV's `route` column meaning "every implemented route". It never
# reaches a filename — build_routes_frame expands it here, and the `route`
# wildcard only ever matches a concrete id.
ALL_ROUTES = "all-routes"

# Several routes in one cell are separated by this.
ROUTE_SEPARATOR = "|"

ROUTE_PATTERN = "|".join(ROUTES)


def expand_route_cell(cell: str) -> list[str]:
    """One `route` cell from scenarios.csv -> the concrete route ids it names."""
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


def build_routes_frame(scenarios: pd.DataFrame) -> pd.DataFrame:
    """(scenario, route) frame — one row per network that will be solved.

    A scenario's route set is the union over its rows, so `route` can vary
    between the rows of one scenario without changing what gets built.
    """
    expanded = (
        scenarios[["scenario", "route"]]
        .drop_duplicates()
        .assign(route=lambda df: df["route"].map(expand_route_cell))
        .explode("route")
        .drop_duplicates()
        .sort_values(["scenario", "route"])
        .reset_index(drop=True)
    )
    return expanded
