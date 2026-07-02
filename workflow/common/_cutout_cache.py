"""Deterministic content-addressed cache for atlite ERA5 cutouts.

Each ERA5 download is expensive (hours of CDS queue, up to GBs), so a cutout
should be re-fetched only when the request itself changes. The cache keys on the
*actual* request parameters — module, bounding box, resolution, time range — not
on the ``(cf_area, start, end)`` filename triple. That closes the stale-bounds
hole of the old ``_backup.nc`` stopgap, where a ``mainland_bbox`` or
``offshore_max_distance_km`` edit changed the real bbox but reused a cutout
cached under the same filename.

Layout (under a gitignored ``cutouts/cache/``):
    <cf_area>_<start>_<end>_<key>.nc     the cached cutout
    <cf_area>_<start>_<end>_<key>.json   its request parameters (for inspection)

Entries are shared with the rule output via hardlink (copy fallback), so the
cache costs no extra disk and survives Snakemake deleting a rule output.

Scope: exact-parameter match only. Coverage-aware reuse (slice a sub-request out
of a larger cached cutout; fill only missing months) is a deferred follow-up —
see TODO.md "Cutout cache".
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
from pathlib import Path

from common._paths import CUTOUTS

log = logging.getLogger(__name__)

CACHE_DIR = CUTOUTS / "cache"


def cache_params(
    module: str,
    x: slice,
    y: slice,
    dx: float | None,
    dy: float | None,
    start_date: str,
    end_date: str,
) -> dict:
    """Canonical, JSON-serialisable description of a cutout request.

    Coordinates are rounded to 5 dp so float noise in the bbox does not spuriously
    change the key.
    """

    def r(value):
        return None if value is None else round(float(value), 5)

    return {
        "module": module,
        "x0": r(x.start),
        "x1": r(x.stop),
        "y0": r(y.start),
        "y1": r(y.stop),
        "dx": r(dx),
        "dy": r(dy),
        "start_date": start_date,
        "end_date": end_date,
    }


def cache_key(params: dict) -> str:
    blob = json.dumps(params, sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()[:12]


def cache_paths(cf_area: str, params: dict) -> tuple[Path, Path]:
    """Return (cutout_path, params_json_path) for the given request."""
    key = cache_key(params)
    stem = f"{cf_area}_{params['start_date']}_{params['end_date']}_{key}"
    return CACHE_DIR / f"{stem}.nc", CACHE_DIR / f"{stem}.json"


def link_or_copy(src: Path, dst: Path) -> None:
    """Materialise `dst` from `src`, preferring a hardlink over a full copy."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copyfile(src, dst)


def store_in_cache(source: Path, cf_area: str, params: dict) -> Path:
    """Add a freshly downloaded cutout to the cache; return the cache path."""
    cutout_path, params_path = cache_paths(cf_area, params)
    link_or_copy(source, cutout_path)
    params_path.write_text(json.dumps(params, indent=2, sort_keys=True))
    log.info(f"stored cutout in cache: {cutout_path.name}")
    return cutout_path
