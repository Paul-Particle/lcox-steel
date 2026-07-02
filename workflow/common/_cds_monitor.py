"""Background CDS job-status logging for long ERA5 downloads.

`atlite.Cutout.prepare()` blocks while ERA5 requests sit in the CDS queue —
often for hours, since the per-user limit is one running job at a time and the
global queue routinely holds thousands. With third-party logging at WARNING the
run looks hung. `cds_progress_logger` wraps `prepare()` in a daemon thread that
polls the CDS jobs endpoint and logs one concise status line per interval, so
the wait is observable.

The status reflects *all* jobs on the account, not only this run's — for the
single cutout downloads this pipeline issues, that is effectively this run.

All CDS access is best-effort: any polling error is logged at DEBUG and swallowed
so the monitor can never break or slow the download it is observing.
"""

from __future__ import annotations

import logging
import threading
from collections import Counter
from contextlib import contextmanager

log = logging.getLogger(__name__)


def _poll_once() -> str:
    """Query CDS and return a one-line status summary (raises on API error)."""
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

    return (
        f"CDS jobs: running={counts.get('running', 0)} "
        f"queued={counts.get('accepted', 0)} "
        f"successful={counts.get('successful', 0)} "
        f"failed={counts.get('failed', 0)}  "
        f"(your queue: {user.get('queued', '?')} queued / "
        f"{user.get('running', '?')} running)"
    )


@contextmanager
def cds_progress_logger(interval_s: float = 30.0):
    """Log CDS job status every `interval_s` seconds until the block exits."""
    stop = threading.Event()

    def loop() -> None:
        while not stop.is_set():
            try:
                log.info(_poll_once())
            except Exception as exc:
                log.debug(f"CDS status poll failed: {exc!r}")
            stop.wait(interval_s)

    thread = threading.Thread(target=loop, name="cds-progress", daemon=True)
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join(timeout=2.0)
