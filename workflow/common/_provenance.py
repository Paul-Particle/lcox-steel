"""Tying an output back to the exact inputs that produced it."""

import hashlib
import json
from pathlib import Path


def input_manifest(input_paths: list[Path]) -> dict:
    """Content digests of a set of input files, and one digest over them all.

    A result is only recoverable if you can say which inputs produced it, and a
    path does not say that — assumptions get edited and timeseries get rebuilt
    under the same name. So each input is hashed by content. `inputs_hash` alone
    is enough to tell two runs apart or match them up; the per-file map says
    which input moved.
    """
    files = {}
    for path in sorted(input_paths):
        files[path.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
    fingerprint = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    manifest = {"inputs": files, "inputs_hash": fingerprint[:12]}
    return manifest
