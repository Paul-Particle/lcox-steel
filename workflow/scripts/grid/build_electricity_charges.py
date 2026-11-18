"""One area's electricity charges, from its workbook or the config default.

An area with a workbook in data/electricity_charges/ gets its grid regimes and
its consumption charge from it (see scripts/grid/_tariff.py). An area without
one gets a file that leaves its grid on the `grid` block of the assumptions, as
before, so a scenario overlay on those fees still applies, and charges nothing
on consumption.
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
    workbook = list(snakemake.input.workbook)

    if workbook:
        charges = _tariff.build_charges(Path(workbook[0]), area)
        grid, consumption = charges["grid"], charges["consumption"]
        log.info(f"{area}: {len(grid['regimes'])} grid regimes from {workbook[0]}")
        for regime in grid["regimes"]:
            log.info(
                f"  {regime['name']}: {regime['capacity_eur_per_mw_per_year']} EUR/MW/yr, "
                f"{regime['energy_eur_per_mwh']} EUR/MWh, {regime['fixed_eur_per_year']} EUR/yr"
            )
        log.info(f"{area}: {consumption['eur_per_mwh']} EUR/MWh on consumption "
                 f"({', '.join(consumption['components']) or 'none'})")
        if charges["reference_consumer"]:
            log.info(f"{area}: reference consumer {charges['reference_consumer']}, "
                     f"{len(grid['reference'])} grid regimes, "
                     f"{charges['consumption']['reference']['eur_per_mwh']} EUR/MWh on consumption")
        if charges["not_read_yet"]:
            log.info(f"{area}: not read yet, the config prices these: {charges['not_read_yet']}")
    else:
        charges = _tariff.default_charges(area)
        log.info(f"{area}: no electricity-charges workbook, the config's grid fees apply")

    out = Path(snakemake.output.charges)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(charges, sort_keys=False, allow_unicode=True), encoding="utf-8")


if __name__ == "__main__":
    main()
