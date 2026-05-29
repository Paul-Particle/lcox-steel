# TODO / Notes

## ENTSO-E bidding zone list

`data/entsoe_cache/entsoe_bidding_zones.csv` is currently hand-crafted (60 zones).
The `entsoe` Python library already enumerates all bidding zones via its `entsoe.Area`
enum. Replace the CSV with a file generated from that enum, or validate directly against
it in `retrieve_entsoe.py`, so the list stays up to date automatically without manual edits.

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
