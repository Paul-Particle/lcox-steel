"""Give a scenario the files of the run that stands for it.

Invoked by Snakemake's `script:` directive (pick_best_run rule in viz.smk). A
scenario is usually one run. An area on the grid that trades through market
zones is run once per zone, and the run with the lowest `report.best_run_by`
stands for the scenario: its report, diagnostic, network and assumptions are
copied up into the scenario's own folder. Every run keeps its files under the
scenario's `.runs/`, beside `runs.csv`, which lists them and says which one was
picked — so `.runs/` is the complete record of what was solved.
"""

import logging
import shutil
from pathlib import Path

import pandas as pd

from common._logging import configure_logging
from common._report_schema import read_report

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)


def main() -> None:
    """Rank the scenario's runs, copy the best one's files up and list them all."""
    metric = snakemake.params.best_run_by
    # One row per run, in the order the rule collected their files.
    runs = pd.concat([read_report(path) for path in snakemake.input.diagnostic], ignore_index=True)
    # A run with no value for the metric is never picked over one that has it.
    best = runs[metric].sort_values(na_position="last").index[0]
    # Every output but the run list has an input of the same name, one path per run.
    for name in ("report", "diagnostic", "network", "overlay", "base_assumptions", "merged_assumptions"):
        Path(getattr(snakemake.output, name)).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(getattr(snakemake.input, name)[best], getattr(snakemake.output, name))
    (
        runs[["run", metric]]
        .assign(picked=runs.index == best)
        .to_csv(snakemake.output.runs, index=False)
    )
    log.info(f"{snakemake.wildcards.scenario}: run {runs.at[best, 'run']} stands for the scenario "
             f"({metric} {runs.at[best, metric]}, {len(runs)} run(s))")


if __name__ == "__main__":
    main()
