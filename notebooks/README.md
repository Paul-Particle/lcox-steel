# Notebooks

Exploratory / analysis notebooks. These are kept **small in git**: cell outputs
are stripped from the stored version via [nbstripout](https://github.com/kynan/nbstripout),
so only source is committed. Your local copy keeps its rendered outputs — only
what git stores is stripped.

## One-time setup (per clone)

The `.gitattributes` entry (`*.ipynb filter=nbstripout`) is committed, but the
filter command itself is per-clone git config. After cloning, run once:

```bash
nbstripout --install
```

`nbstripout` ships in `environment.yaml`. Without this step, notebooks still
work — the filter entry is simply a no-op, so outputs could sneak into a commit.
`git add`/`git commit` handle the stripping automatically once it's installed.

## Notebooks

- `res_cf_coastal_cf_exploration.ipynb` — coastal onshore-wind capacity-factor
  exploration (land/sea contamination of ERA5 cells). A fully-rendered copy with
  outputs lives on the `pr8-review-analysis` branch.
