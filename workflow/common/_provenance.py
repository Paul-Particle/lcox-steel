"""Tying an output back to the exact inputs and code that produced it."""

import hashlib
import json
import sys
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
        files[path.as_posix()] = _digest(path)
    fingerprint = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    manifest = {"inputs": files, "inputs_hash": fingerprint[:12]}
    return manifest


def code_manifest(repo_root: Path, extra_paths: tuple[Path, ...] = ()) -> dict:
    """Content digests of the repo's own modules that were loaded, and one over them all.

    `inputs_hash` says which data made a result; this says which code did, so a number
    that moved because the model changed can be told from one that moved because an
    assumption did. Read off `sys.modules` rather than a written-down list, because a
    list goes stale the first time an import changes and nothing notices.

    Only modules under `repo_root` are counted: what the environment ships is the
    environment's business, and pinning it is what the lockfile is for.

    `extra_paths` exists for one case. Snakemake runs a `script:` entry point from a
    copy under `.snakemake/scripts/`, so the entry script is in `sys.modules` under a
    temporary path and never under its own. A caller that wants itself counted passes
    its repo path here.
    """
    files = {}
    for module in list(sys.modules.values()):
        module_file = getattr(module, "__file__", None)
        if module_file is None:
            continue
        resolved = Path(module_file).resolve()
        if repo_root in resolved.parents:
            files[resolved.relative_to(repo_root).as_posix()] = _digest(resolved)
    for path in extra_paths:
        resolved = Path(path).resolve()
        if resolved.exists():
            files[resolved.relative_to(repo_root).as_posix()] = _digest(resolved)
    fingerprint = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    return {"code": files, "code_hash": fingerprint[:12]}


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]
