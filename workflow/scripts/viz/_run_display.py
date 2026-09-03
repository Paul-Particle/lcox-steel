"""How a run is labelled in a chart.

A scenario is an umbrella over runs, so a chart of one scenario holds several
areas, date ranges and routes. Both plot scripts label them the same way, and
neither reads snakemake at import time, so this lives here where a test can
reach it.
"""

import pandas as pd


def run_label(row: pd.Series) -> str:
    """A bar's label: the run key minus the scenario, which the title carries.

    The route alone no longer identifies a row — a scenario spanning areas would
    render every one of them under the same name. The date range joins in only
    when the chart actually spans more than one year, so the common
    single-year case stays short.
    """
    label = f"{row['area']} {row['route']}"
    years = {str(row["start_date"])[:4], str(row["end_date"])[:4]}
    return label if len(years) == 1 else f"{label} {'-'.join(sorted(years))}"
