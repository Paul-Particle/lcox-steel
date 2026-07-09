"""Background CDS job-status logging for long ERA5 downloads.

`atlite.Cutout.prepare()` blocks while ERA5 requests sit in the CDS queue —
often for hours, since the per-user limit is one running job at a time and the
global queue routinely holds thousands. With third-party logging at WARNING the
run looks hung. `cds_progress_logger` wraps `prepare()` in a daemon thread that
polls the CDS jobs endpoint and logs progress.

Log discipline (matters for humans *and* agents tailing the file):
- A line is emitted only when the job-status counts **change**, so the log is a
  timeline of real events, not hundreds of identical snapshots.
- While nothing changes, a single "still working" heartbeat is emitted every
  `heartbeat_s` seconds, stating how long the state has held and that a normal
  ERA5 job runs ~20-25 min. This is the fix for the classic misread where one
  point-in-time snapshot ("running=1 queued=4") is mistaken for a stall.
- `done_this_run` counts completions since monitoring began (0 -> N as the
  cutout's feature requests finish), giving an at-a-glance progress signal that
  the account-wide `successful` count does not.

Polling still happens every `interval_s` so changes are noticed promptly; only
the *logging* is throttled. The status reflects all jobs on the account, not
only this run's — for the single cutout downloads this pipeline issues, that is
effectively this run. All CDS access is best-effort: any polling error is logged
at DEBUG and swallowed so the monitor can never break the download it observes.
"""

from __future__ import annotations

import logging
import threading
import time
from collections import Counter
from contextlib import contextmanager

log = logging.getLogger(__name__)


def _poll() -> tuple[dict, dict]:
    """Return (status_counts, your_queue) from CDS (raises on API error)."""
    import cdsapi

    client = cdsapi.Client(quiet=True)
    # `._json_dict` is cdsapi-internal; guarded by the caller's try/except.
    jobs = client.client.get_jobs()._json_dict["jobs"]
    counts = Counter(job["status"] for job in jobs)

    user = {}
    for job in jobs:
        status = job.get("metadata", {}).get("qos", {}).get("status", {})
        if status.get("user"):
            user = status["user"][0]
            break

    return dict(counts), user


def _mins(seconds: float) -> str:
    return f"{int(seconds // 60)}m"


@contextmanager
def cds_progress_logger(interval_s: float = 60.0, heartbeat_s: float = 300.0):
    """Log CDS job status while the block runs (see module docstring)."""
    stop = threading.Event()
    start = time.monotonic()

    def loop() -> None:
        last_key = None          # (running, queued, successful, failed)
        last_change = start
        last_heartbeat = start
        base_successful = None   # successful count when monitoring began

        while not stop.is_set():
            now = time.monotonic()
            try:
                counts, _user = _poll()
                key = (
                    counts.get("running", 0),
                    counts.get("accepted", 0),
                    counts.get("successful", 0),
                    counts.get("failed", 0),
                )
                if base_successful is None:
                    base_successful = key[2]

                summary = (
                    f"running={key[0]} queued={key[1]} "
                    f"done_this_run={max(0, key[2] - base_successful)} "
                    f"failed={key[3]}"
                )

                if key != last_key:
                    if last_key is not None and key[3] > last_key[3]:
                        log.warning(f"[+{_mins(now - start)}] CDS: {summary}  (a CDS job FAILED)")
                    else:
                        note = ""
                        if last_key is not None and key[2] > last_key[2]:
                            note = "  (a CDS job completed)"
                        log.info(f"[+{_mins(now - start)}] CDS: {summary}{note}")
                    last_key, last_change, last_heartbeat = key, now, now
                elif now - last_heartbeat >= heartbeat_s:
                    if key[0] > 0:
                        why = "a job is running (ERA5 jobs take ~20-25 min) — expected, not a stall"
                    elif key[1] > 0:
                        why = (
                            "jobs queued but none started yet — waiting on the CDS global "
                            "queue, which can be slow at peak times; not a stall"
                        )
                    else:
                        why = "no active CDS jobs — atlite is likely processing/writing locally"
                    log.info(
                        f"[+{_mins(now - start)}] CDS: still working — {summary} "
                        f"(unchanged {_mins(now - last_change)}; {why})"
                    )
                    last_heartbeat = now
            except Exception as exc:
                log.debug(f"CDS status poll failed: {exc!r}")
            stop.wait(interval_s)

    thread = threading.Thread(target=loop, name="cds-progress", daemon=True)
    log.info(
        f"CDS: monitoring queue (status logged on change; heartbeat every {_mins(heartbeat_s)})"
    )
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join(timeout=2.0)
        log.info(f"CDS: monitoring stopped after {_mins(time.monotonic() - start)}")
