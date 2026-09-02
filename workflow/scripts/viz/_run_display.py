"""How a run is labelled in a chart, and which runs a chart shows.

A scenario is an umbrella over runs, so a chart of one scenario holds several
areas, date ranges and routes. Both plot scripts need the same answers, and
neither reads snakemake at import time, so these live here where a test can
reach them.
"""

import logging

import pandas as pd

log = logging.getLogger(__name__)


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


def best_zones_only(df: pd.DataFrame) -> pd.DataFrame:
    """Drop the zones the report did not flag as their country's cheapest.

    Australia is solved once per NEM region; the report ranks them and flags one
    (see compile_report.mark_best_in_country). Charting all of them would show
    the same country several times over.
    """
    if "best_in_country" not in df.columns:
        return df
    kept = df[df["best_in_country"].astype(bool)]
    hidden = len(df) - len(kept)
    if hidden:
        log.info(f"hiding {hidden} zone(s) beaten by a cheaper one in the same country")
    return kept
