# TODO / Notes

## ENTSO-E bidding zone list

`data/entsoe_cache/entsoe_bidding_zones.csv` is currently hand-crafted (60 zones).
The `entsoe` Python library already enumerates all bidding zones via its `entsoe.Area`
enum. Replace the CSV with a file generated from that enum, or validate directly against
it in `retrieve_entsoe.py`, so the list stays up to date automatically without manual edits.

## Cutout cache (like the grid data cache)

Grid data uses a persistent per-variant processed cache (`resources/entsoe/{variant}.parquet`,
`resources/nem/{variant}.parquet`) that accumulates area/month slices and avoids re-downloading
covered periods. Cutouts have no equivalent: each `download_cutout` rule invocation is a
one-shot ERA5 pull for a fixed `(cf_area, start_date, end_date)` triple.

**Current stopgap**: `download_cutout.py` checks for a sibling `cutouts/{name}_backup.nc`
file at the start of its run and copies from it instead of hitting CDS. Renaming an
existing cutout to `*_backup.nc` thus pins it across code-trigger reruns. `protected()` is
no longer used on the rule output — the backup IS the safety net.

A proper cutout cache would: store ERA5 data by area in a persistent location (analogous to
`data/entsoe_cache/`), check spatial and temporal coverage before requesting, and support
partial fills. This would also make the `geo → cutout → timeseries` chain more natural:
`build_regions`/`build_offshore_regions` would define the area, the cache layer would ensure
coverage, and `build_res_cf_profile` would slice from the cache. When this lands, the backup
hack and its README/HANDOFF callouts can come out.

## CDS download monitoring

Atlite ERA5 cutout downloads go through the CDS API. With `monthly_requests=True`
(set in `download_cutout.py`), atlite submits one CDS job per month — up to 12 jobs
for a full year. Each job queues independently; global queue depths of 5000+ are
common and waits of several hours are normal.

### Checking status from the terminal

```python
import cdsapi
c = cdsapi.Client(quiet=True)
for j in c.client.get_jobs()._json_dict['jobs']:
    qos  = j['metadata']['qos']['status']
    user = qos.get('user', [{}])[0]
    print(
        j['status'], j['processID'],
        '| user queue:', user.get('queued', 0),
        'running:', user.get('running', 0),
    )
```

Key fields:
- `status` — `accepted` (queued), `running`, `successful`, `failed`
- `user.queued` / `user.running` — your own per-user position (limit: 1 running at a time)
- Global queue (`qos.limit`) gives context on overall wait times

### Future: automatic status updates

Wire this into Snakemake logging or a small polling script that prints progress
as monthly jobs complete, so long cutout runs don't require manual checks.

## viz/_helpers.py cleanup

`workflow/scripts/viz/_helpers.py` was ported from a standalone notebook and has several
issues to resolve before it's considered production-ready:

- `darkmode()` silently switches plot style based on wall-clock hour (≤18:00 = light,
  >18:00 = dark). Breaks reproducibility. Either remove or make explicit with a `force=`
  parameter.
- Global `save` flag + `toggle_save()` / `set_save()` / `save_is_on()` setters are
  redundant — `save_fig()` already takes `save_condition`. Pass `save=` directly instead.
- `OUTPUTPATH = Path('./results/')` and `DATAPATH = Path('./data/')` use cwd-relative
  paths; should use `common._paths.RESULTS` / `common._paths.DATA` for consistency.
- The 258-line inline colormap at the bottom (`cm_data`) could move to `viz/_cmap.py`
  so the rest of the file stays readable.
- No module docstring (now added as a WIP marker).
