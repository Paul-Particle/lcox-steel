"""The one carrier vocabulary the two markets and the factor table share.

This is what makes a single emission factor table serve ENTSO-E and the NEM: the
two downloaders were deliberately given the same short carrier names, and those
names are the keys of `emissions.electricity_t_co2e_per_mwh`. Nothing translates
between them anywhere, so the three have to be held to each other here — a
carrier one downloader can emit and the table has no factor for makes the report
refuse the run, and a factor for a carrier nobody emits is dead weight that
invites a typo.
"""

from pathlib import Path

import yaml

import download_entsoe
import download_nem  # sys.path set by conftest

REPO_ROOT = Path(__file__).resolve().parents[1]
FACTORS = yaml.safe_load((REPO_ROOT / "config" / "assumptions.yaml").read_text())[
    "emissions"
]["electricity_t_co2e_per_mwh"]


def test_every_carrier_either_downloader_can_emit_has_a_factor():
    """A missing factor is a run the report refuses, so it cannot be found later."""
    emitted = set(download_entsoe.CARRIER_NAMES.values()) | set(download_nem.CARRIER_PATTERNS)
    assert emitted - set(FACTORS) == set()


def test_the_factor_table_has_no_carrier_neither_downloader_emits():
    """The other direction: a factor nothing can ever look up is dead weight."""
    emitted = set(download_entsoe.CARRIER_NAMES.values()) | set(download_nem.CARRIER_PATTERNS)
    assert set(FACTORS) - emitted == set()


def test_the_two_markets_differ_only_where_they_are_documented_to():
    """ENTSO-E splits hydro into river and reservoir; the NEM reports it whole.
    That one key is the whole difference, and it is why the table carries all
    three."""
    entsoe = set(download_entsoe.CARRIER_NAMES.values())
    nem = set(download_nem.CARRIER_PATTERNS)
    assert nem - entsoe == {"hydro"}


def test_every_factor_carries_all_three_bases():
    """A factor keyed on one basis only would be read on whichever basis the run
    was on, silently — which is what the basis label exists to prevent."""
    for carrier, values in FACTORS.items():
        assert set(values) == {"combustion", "delegated_act", "lifecycle"}, carrier
