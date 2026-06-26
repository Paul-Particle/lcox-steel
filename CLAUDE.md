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

## VM NUL-byte corruption (persistent issue)

### What happens
The host VM's filesystem occasionally writes **trailing NUL (`\x00`) padding
bytes** when a file is flushed to disk (a sparse-write / page-size rounding
artefact at the hypervisor level). This can affect any file the VM writes —
Python scripts, `.smk` rule files, the Snakefile itself, etc.

**The git consequence is severe**: if a corrupted file is staged, git detects
the NUL bytes and stores the object as a *binary blob*. This means:
- `git diff` produces no textual output for that file (just "binary files differ")
- The corruption is silently baked into every future commit until manually fixed
- `git fsck` may report `badRefName` / integrity errors

Commit `a6cd7b9` is a historical example where `workflow/Snakefile` was
committed as a binary blob with 29 trailing NULs.

### Three layers of defence now in place

**Layer 1 — `.gitattributes`** (repo root)
Declares all source/config file types (`Snakefile`, `*.smk`, `*.py`, `*.yaml`,
`*.csv`, …) as `text`. Git will never classify them as binary blobs regardless
of stray NULs. True binary formats (`.parquet`, `.feather`, `.nc`, …) are
explicitly marked `binary`.

**Layer 2 — `.githooks/pre-commit`**
Runs automatically on every `git commit` (wired via `core.hookspath=.githooks`
already set in the repo config — no extra setup needed). For each staged text
file it reads the *staged blob*, strips any NUL bytes, re-stages the clean
version, and prints a warning. The commit still proceeds — it self-heals.

**Layer 3 — Snakemake `onstart` block** (`workflow/Snakefile`)
Runtime last-resort. At the start of every Snakemake run, scans
`scripts/**/*.py`, `rules/**/*.smk`, and the Snakefile itself for NUL bytes
and strips them in-place before any rules execute.

### If you see NUL-byte corruption on a file not covered above
1. Strip manually: `python3 -c "p=open('path/to/file','rb+'); d=p.read(); p.seek(0); p.write(d.replace(b'\x00',b'')); p.truncate()"`
2. Re-stage: `git add path/to/file`
3. Check the staged blob is clean: `git cat-file blob $(git ls-files --stage path/to/file | awk '{print $2}') | python3 -c "import sys; d=sys.stdin.buffer.read(); print(d.count(b'\\x00'), 'NULs')"`
4. Consider extending `.gitattributes` or the `onstart` glob to cover that file type.
