"""Unit tests for the input fingerprint a solve stamps into its network.

A report row carries `inputs_hash`; the network carries the per-file map behind
it. What has to hold is that the fingerprint tracks *content*, not paths or
mtimes — assumptions get edited and timeseries get rebuilt under the same name,
and a fingerprint that missed that would connect a number to the wrong inputs.
"""

from common._provenance import input_manifest  # sys.path set by conftest


def _write(directory, name, text):
    path = directory / name
    path.write_text(text)
    return path


def test_manifest_lists_every_input_by_content(tmp_path):
    paths = [
        _write(tmp_path, "assumptions.yaml", "wacc: 0.07"),
        _write(tmp_path, "solar.parquet", "profile"),
    ]
    manifest = input_manifest(paths)

    assert set(manifest["inputs"]) == {p.as_posix() for p in paths}
    assert len(manifest["inputs_hash"]) == 12


def test_the_same_inputs_fingerprint_the_same(tmp_path):
    """Order must not matter — the rule hands the files over in whatever order."""
    first = _write(tmp_path, "a.yaml", "wacc: 0.07")
    second = _write(tmp_path, "b.parquet", "profile")

    assert input_manifest([first, second]) == input_manifest([second, first])


def test_an_edited_input_changes_the_fingerprint(tmp_path):
    """The case this exists for: a file edited in place, same name, same run key."""
    assumptions = _write(tmp_path, "assumptions.yaml", "wacc: 0.07")
    before = input_manifest([assumptions])

    assumptions.write_text("wacc: 0.09")
    after = input_manifest([assumptions])

    assert before["inputs_hash"] != after["inputs_hash"]
    assert before["inputs"][assumptions.as_posix()] != after["inputs"][assumptions.as_posix()]


def test_a_dropped_input_changes_the_fingerprint(tmp_path):
    """A run solved without its assumptions overlay is not the same run."""
    base = _write(tmp_path, "assumptions.yaml", "wacc: 0.07")
    overlay = _write(tmp_path, "assumptions_x.yaml", "wacc: 0.09")

    assert (
        input_manifest([base])["inputs_hash"]
        != input_manifest([base, overlay])["inputs_hash"]
    )
