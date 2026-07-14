# Tests

## Grid pipeline output completeness

`test_window_completeness.py` pins the behaviour of `retrieve_entsoe`'s
`_assert_window_complete` guard, which runs on every pipeline output and fails
loudly when the produced window has holes. Truncated raw-cache months (a fetch
that stopped mid-month) and partial downloads are otherwise silent — the slice
just has fewer rows, and the hole only surfaces far downstream as NaN after a
reindex (e.g. a solve aligning prices to a full capacity-factor year). The
energy-charts gate below does **not** catch this: it masks NaN before comparing
and tolerates hundreds of missing hours.

- `dayahead` is resampled to a clean hourly grid, so every hour of the window
  must be present and non-null.
- `full` spans a mixed 15-min/hourly resolution (ENTSO-E switched DE_LU
  day-ahead to 15-min in Oct 2025), so instead of a fixed grid the guard flags
  any gap larger than the ffill tolerance — which only occurs on truncation.

These tests use synthetic frames, so they need no raw cache and run anywhere.

## Grid pipeline vs. energy-charts.info (issue #20)

`test_entsoe_vs_energycharts.py` cross-checks the ENTSO-E grid pipeline against
[energy-charts.info](https://energy-charts.info) (Fraunhofer ISE), which
publishes the same German day-ahead price, generation-by-type, and load. Both
ultimately draw on ENTSO-E, so they should agree closely — this turns the
"eyeball the dashboard" check into something automatic.

**What it does**

- Rebuilds the `full`-variant frame for **DE_LU / 2023** from the local ENTSO-E
  raw monthly cache, through the pipeline's own `_process_full_month` (so it
  exercises the real processing code).
- Compares it, on the hourly overlap, to a pinned energy-charts reference
  (`data/energycharts_de_2023_hourly.csv`).
- Asserts per-metric correlation and relative-MAE thresholds (`price`, `solar`,
  `wind_onshore`, `wind_offshore`, `load`, and the `wind`/`res` aggregates).
  Thresholds carry ~2x headroom over the agreement measured on 2023 data. A
  small positive bias (our value slightly above energy-charts, from their
  revisions/gap-filling) is expected and does not fail the gate.

**Running**

```bash
pytest                        # from the repo root
pytest tests/test_entsoe_vs_energycharts.py -v
```

The pipeline side needs the ENTSO-E raw monthly cache
(`data/entsoe_cache/DE_LU/2023-*/`) present locally; if it isn't, the test
**skips** with a message rather than failing. The energy-charts side is the
committed CSV, so no network access or API key is needed to run the gate.

**Refreshing the reference**

The reference CSV is pinned for determinism. To regenerate it (e.g. for a new
year), run:

```bash
python tests/data/refresh_energycharts_reference.py
```
