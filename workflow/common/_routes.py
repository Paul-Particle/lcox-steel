"""The route vocabulary, and the expansion from scenario rows to runs.

A scenario name is an umbrella label and carries no meaning for the workflow.
What identifies one network is the run key — scenario, area, date range, route —
and every part of it but the scenario comes from the rows filed under that name.

Deliberately free of heavy imports: the Snakefile reads this on every DAG build.
"""

import pandas as pd

ROUTES = ("h2-only", "h2-dri-eaf", "ng-dri-eaf", "mix-dri-eaf", "moe-eaf", "ew-eaf")

# Reserved for the trade scenarios (produce abroad, ship the product to the EU).
# Not implemented: `*-export` appears in config/scenarios.csv as commented rows
# only, and is excluded from ROUTES until a route actually builds a network.
PLANNED_ROUTES = ("h2-dri-eaf-export", "ng-dri-eaf-export", "mix-dri-eaf-export",
                  "moe-eaf-export", "ew-eaf-export")

# Value in the CSV's `route` column meaning "every implemented route". It never
# reaches a filename — build_runs_frame expands it here, and the `route`
# wildcard only ever matches a concrete id.
ALL_ROUTES = "all-routes"

# Several routes in one cell are separated by this.
ROUTE_SEPARATOR = "|"

ROUTE_PATTERN = "|".join(ROUTES)

# The columns that, with the scenario, identify one network.
RUN_KEY = ["scenario", "area", "start_date", "end_date"]


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


def build_runs_frame(scenarios: pd.DataFrame) -> pd.DataFrame:
    """(scenario, area, start_date, end_date, route) — one row per network solved.

    A scenario's route set is the union over its rows, so `route` can vary
    between the rows of one run without changing what gets built.
    """
    runs = (
        scenarios[RUN_KEY + ["route"]]
        .drop_duplicates()
        .assign(route=lambda df: df["route"].map(expand_route_cell))
        .explode("route")
        .drop_duplicates()
        .sort_values(RUN_KEY + ["route"])
        .reset_index(drop=True)
    )
    return runs
