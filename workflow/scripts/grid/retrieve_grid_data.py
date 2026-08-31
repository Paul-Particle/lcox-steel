"""Grid market data for one (area, variant, date-range) slice.

Single entry point for every `tech=grid` row in config/scenarios.csv. The area
registry (config/config.yaml `areas`) says which wholesale market an area trades
in, and that decides which source implementation runs. The sources themselves
are unaware of this dispatcher.

Areas are a flat namespace — nothing infers anything from the shape of an area
code. An area with no `market` entry has no price series, so a tech=grid row on
it is an error; that is the only thing separating an area you can grid-connect
from one you can only model islanded.
"""

import logging

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
import _entsoe
import _nem

configure_logging(snakemake)
log = logging.getLogger(__name__)

SOURCES = {"entsoe": _entsoe.retrieve, "nem": _nem.retrieve}


def main() -> None:
    """Route the request to the source that carries this area's prices."""
    area = snakemake.wildcards.area
    market = snakemake.params.market

    if not market:
        priced = sorted(a for a, cfg in snakemake.config["areas"].items() if cfg.get("market"))
        raise ValueError(
            f"area {area!r} has no `market` entry in config.yaml, so no wholesale "
            f"price series exists for it — a tech=grid row needs one of {priced}. "
            f"Model {area!r} islanded, or add its market and price source."
        )
    if market not in SOURCES:
        raise ValueError(
            f"area {area!r} declares market {market!r}; expected one of {sorted(SOURCES)}"
        )

    log.info(f"{area}: retrieving {snakemake.wildcards.variant} from {market}")
    SOURCES[market](snakemake)


if __name__ == "__main__":
    main()
