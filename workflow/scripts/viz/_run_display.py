"""How a run is labelled in a chart.

A scenario is an umbrella over runs, so a chart of one scenario holds several
areas, date ranges and routes. Both plot scripts label them the same way, and
neither reads snakemake at import time, so this lives here where a test can
reach it.
"""

import pandas as pd


def run_labels(runs: pd.DataFrame) -> pd.Series:
    """A label per run: the run key minus the scenario, which the title carries.

    Labelled against the frame as a whole, because whether the date range is
    needed to tell two bars apart is a property of the chart and not of either
    bar: the range joins in when `runs` holds more than one window, and a window
    that crosses new year names both years.

    Callers that plot a subset pass the whole frame and select afterwards, so a
    run keeps one label across every chart built from the same report.
    """
    label = runs["area"].astype(str) + " " + runs["route"].astype(str)
    windows = runs[["start_date", "end_date"]].drop_duplicates()
    if len(windows) == 1:
        return label

    start_year = runs["start_date"].astype(str).str[:4]
    end_year = runs["end_date"].astype(str).str[:4]
    span = start_year.where(start_year == end_year, start_year + "-" + end_year)
    return label + " " + span
