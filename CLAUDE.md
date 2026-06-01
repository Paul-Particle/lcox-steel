# Project conventions

## Logging
- f-strings are fine in log calls (`log.info(f"...")`). The old PyPSA-Eur-style
  "%-style only" rule was dropped — lazy-eval/exception-safety wins don't
  bite at this codebase's scale.
- The `format:` template in `config/config.yaml` stays %-style — that's
  `logging.Formatter`'s own substitution syntax (`%(asctime)s`, `%(message)s`),
  not the same thing as call-site formatting.

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

## Two distinct script patterns — don't conflate them

Every Snakemake-invoked script in this repo carries the **stub-fallback
block** so linters/IDEs can resolve module-level `snakemake.*` references:

```python
if "snakemake" not in globals():
    from common._stubs import snakemake
```

Keep this in **every** script (`grid/`, `h2_dri/`, `res_cf/`, `viz/`).
It's a linter shim — independent of whether the script can actually run
standalone.

The **hardcoded-default block** (`_VAR = "de"` etc. with an
`if "snakemake" in globals() and hasattr(snakemake, "wildcards"): override`
guard) is a *separate* pattern that lets a script run without
Snakemake. Only `res_cf/` scripts use it by design. Don't add it to
`grid/`, `h2_dri/`, or `viz/` — those are Snakemake-only.
