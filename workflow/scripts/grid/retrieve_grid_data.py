"""Grid market data for one (area, variant, date-range) slice.

Single entry point for every `tech=grid` row in config/scenarios.csv. The area
is a market-zone code, and the zone decides which source implementation runs:
AEMO NEM regions go to `_nem`, everything else to ENTSO-E via `_entsoe`. Both
modules keep their own processed cache and are unaware of this dispatcher.

Areas in the scenario CSV come from two namespaces (see config/config.yaml):
market zones like DE_LU / VIC1, and three-letter national codes like AUS used
by the RES pipeline where a whole-country cutout is wanted. A country has no
single wholesale price, so a national code on a grid row is rejected here —
this is the only place that knows the difference.
"""

import logging
import re
from pathlib import Path

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
import _entsoe
import _nem

configure_logging(snakemake)
log = logging.getLogger(__name__)

# AEMO NEM regions. Every other market zone is assumed to be an ENTSO-E bidding
# zone; _entsoe validates it against data/entsoe_cache/entsoe_bidding_zones.csv.
NEM_AREAS = ("NSW1", "QLD1", "SA1", "TAS1", "VIC1")

# Three uppercase letters = an ISO 3166-1 alpha-3 national code (see the area
# convention in config/config.yaml). Market zones never match this shape.
NATIONAL_AREA = re.compile(r"[A-Z]{3}")


def main() -> None:
    """Route the request to the NEM or ENTSO-E implementation for this area."""
    area = snakemake.wildcards.area

    if NATIONAL_AREA.fullmatch(area):
        raise ValueError(
            f"{area!r} is a national area code, which has no wholesale price series. "
            f"A tech=grid row needs a market zone — an ENTSO-E bidding zone (DE_LU, "
            f"FR, ES, ...) or an AEMO NEM region ({', '.join(NEM_AREAS)}). National "
            f"codes are for RES rows only, where they select a whole-country cutout."
        )

    source = "NEM" if area in NEM_AREAS else "ENTSO-E"
    log.info(f"{area}: retrieving {snakemake.wildcards.variant} from {source}")
    retrieve = _nem.retrieve if area in NEM_AREAS else _entsoe.retrieve
    retrieve(snakemake)


if __name__ == "__main__":
    main()
