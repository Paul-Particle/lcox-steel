"""One area's electricity charges, from its workbook or the config default.

An area with a workbook in data/electricity_charges/ gets its grid regimes, its
consumption charge and its connection costs from it, for every consumer the
workbook names (see scripts/grid/_tariff.py); the run picks one. An area
without one gets a file that leaves its grid and its connection on the `grid`
block of the assumptions, as before, so a scenario overlay on those fees still
applies, and charges nothing on consumption.
"""

import logging
from pathlib import Path

import yaml

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.grid import _tariff

configure_logging(snakemake)
log = logging.getLogger(__name__)


def main() -> None:
    area = snakemake.wildcards.area
    # The rule's input is optional: one path when the area has a workbook,
    # none when it has not.
    workbook = list(snakemake.input.workbook)

    if workbook:
        charges = _tariff.build_charges(Path(workbook[0]), area)
        # What was read, for checking a changed workbook in the log rather than
        # in the YAML: per consumer, each regime's rates, the charge on
        # consumption and the connection.
        log.info(f"{area}: charges for {', '.join(charges['consumers'])} from {workbook[0]}")
        for consumer, own in charges["consumers"].items():
            grid, consumption, connection = own["grid"], own["consumption"], own["connection"]
            log.info(f"{area}/{consumer}: {len(grid['regimes'])} grid regimes")
            for regime in grid["regimes"]:
                log.info(
                    f"  {regime['name']}: {regime['capacity_eur_per_mw_per_year']} EUR/MW/yr, "
                    f"{regime['energy_eur_per_mwh']} EUR/MWh, {regime['fixed_eur_per_year']} EUR/yr"
                )
            log.info(f"{area}/{consumer}: {consumption['eur_per_mwh']} EUR/MWh on consumption "
                     f"({', '.join(consumption['components']) or 'none'})")
            log.info(f"{area}/{consumer}: connection {connection['capex_eur_per_mw']} EUR/MW plus "
                     f"{connection['capex_eur']} EUR once ({', '.join(connection['components']) or 'none'})")
    else:
        charges = _tariff.default_charges(area)
        log.info(f"{area}: no electricity-charges workbook, the config's grid fees apply")

    # In the order build_charges gives it, with ä, ü and € written as they are.
    out = Path(snakemake.output.charges)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(charges, sort_keys=False, allow_unicode=True), encoding="utf-8")


if __name__ == "__main__":
    main()
