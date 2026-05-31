# Local stand-ins for snakemake helpers that don't ship yet.
#
# Delete this file (and its `include:` line in the Snakefile) once the upstream
# equivalents land — there is in-progress work toward a proper `optional()`
# helper on snakemake itself.

from pathlib import Path


def optional(template: str):
    """Defer wildcard substitution + existence check until job-evaluation time.

    Returns an input function that, given `wildcards`, substitutes them into
    `template` and yields the resolved path if the file exists on disk, else `[]`.

    Use exactly like a snakemake helper:

        input:
            cfg = optional("config/assumptions_{project}_{scenario}.yaml"),

    The wrapped attribute is always a Namedlist with 0 or 1 entries.
    """
    def _resolve(wildcards):
        path = template.format(**wildcards)
        return [path] if Path(path).exists() else []
    return _resolve
