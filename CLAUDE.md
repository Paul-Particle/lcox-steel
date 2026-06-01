# Project conventions

## Logging
- f-strings are fine in log calls (`log.info(f"...")`). The old PyPSA-Eur-style
  "%-style only" rule was dropped — lazy-eval/exception-safety wins don't
  bite at this codebase's scale.
- Don't rewrite existing %-style calls just to switch styles.

## Snakefile / `.smk` style
Keep raw Python in `Snakefile` and `workflow/rules/*.smk` to a minimum.
Inside rule bodies (`input:`, `output:`, `params:`, `log:`, …) and at module
scope, **do not** use:

- f-strings (use Snakemake's wildcard-template strings: `"path/{wildcard}/..."`),
- lambdas,
- list comprehensions,
- input functions (def-based callables passed as `input=...`).

Reach for Snakemake's declarative helpers wherever possible:

- `collect(...)` / `expand(...)` for fan-out over wildcards or dataframe rows.
- `lookup(query=..., within=...)` / `lookup(dpath=..., within=config)` for
  pulling values out of a dataframe or the config dict.

Imports and top-of-Snakefile dataframe loads are fine — the constraint is
about per-rule Python, not the file header.
