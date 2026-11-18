"""One area's grid tariff as LP regimes, from its workbook or the config default.

An area with a workbook in data/grid_tariffs/ gets its regimes from it (see
scripts/grid/_tariff.py). An area without one gets a marker that tells the
solve to price its grid import from the `grid` block of the assumptions, as
before, so a scenario overlay on those fees still applies.
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
        tariff = _tariff.build_tariff(Path(workbook[0]), area)
        log.info(
            f"{area}: {len(tariff['regimes'])} regimes from {workbook[0]}, "
            f"subsidy after the solve {tariff['post_solve']['subsidy_eur_per_mwh']} EUR/MWh "
            f"({', '.join(tariff['post_solve']['subsidy_components']) or 'none'})"
        )
        for regime in tariff["regimes"]:
            log.info(
                f"  {regime['name']}: {regime['capacity_eur_per_mw_per_year']} EUR/MW/yr, "
                f"{regime['energy_eur_per_mwh']} EUR/MWh, {regime['fixed_eur_per_year']} EUR/yr"
            )
        if tariff["not_read_yet"]:
            log.info(f"{area}: not read yet, the config prices these: {tariff['not_read_yet']}")
    else:
        tariff = _tariff.default_tariff(area)
        log.info(f"{area}: no tariff workbook, the config's grid fees apply")

    out = Path(snakemake.output.tariff)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(tariff, sort_keys=False, allow_unicode=True), encoding="utf-8")


if __name__ == "__main__":
    main()
