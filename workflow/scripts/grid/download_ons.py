"""ONS download primitives — imported by retrieve_ons.py.

ONS publishes one parquet per dataset per calendar year on a public S3 bucket
(no API key, no rate limit, no pagination), so the raw cache here is a plain
per-(dataset, year) file rather than the per-month fetches ENTSO-E and NEM need.

Both datasets come in long format keyed on (id_subsistema, din_instante); the
pivot to the (area, metric) wide shape used by the rest of the grid pipeline
happens in retrieve_ons.py.
"""

import logging
from pathlib import Path

import pandas as pd
import requests

# Module-level logger only — the rule script (retrieve_ons.py) installs handlers.
log = logging.getLogger(__name__)

S3_ROOT = "https://ons-aws-prod-opendata.s3.amazonaws.com/dataset"

# dataset key -> (bucket sub-path, filename stem). Paths are lowercase, filenames upper.
DATASETS = {
    "cmo": ("cmo_tm", "CMO_SEMIHORARIO"),
    "balance": ("balanco_energia_subsistema_ho", "BALANCO_ENERGIA_SUBSISTEMA"),
}

# Subsystem codes as they appear in id_subsistema. "SIN" (national total) also
# appears in the balance dataset but is not a market area, so it is not listed.
SUBSYSTEMS = ["SE", "S", "NE", "N"]


def download_year(dataset: str, year: int, cache_dir: Path) -> Path:
    """Fetch one dataset-year parquet into the raw cache, returning its path.

    A present file is left alone — ONS revises recent months in place, so refresh
    by deleting the year file rather than by passing a flag.
    """
    sub_path, stem = DATASETS[dataset]
    cache_path = cache_dir / f"{stem}_{year}.parquet"
    if cache_path.exists():
        return cache_path

    url = f"{S3_ROOT}/{sub_path}/{stem}_{year}.parquet"
    log.info(f"{dataset}/{year}: fetching {url}")
    response = requests.get(url, timeout=300)
    response.raise_for_status()

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(response.content)
    log.info(f"{dataset}/{year}: cached ({len(response.content) / 1e6:.1f} MB)")
    return cache_path


def read_area_year(dataset: str, area: str, year: int, cache_dir: Path) -> pd.DataFrame:
    """Return one subsystem's rows for one year, indexed by din_instante (Brasilia-naive).

    The `val_*` columns are coerced to float because ONS is not consistent about
    it across year files — CMO 2025 arrives as float64 and CMO 2026 as str, with
    no change to the published data dictionary.
    """
    raw = pd.read_parquet(download_year(dataset, year, cache_dir))
    area_rows = raw[raw["id_subsistema"] == area].copy()
    if area_rows.empty:
        raise ValueError(
            f"{dataset} {year}: no rows for subsystem {area!r}. "
            f"Expected one of {SUBSYSTEMS}."
        )
    area_rows = area_rows.set_index("din_instante").sort_index()
    area_rows = area_rows.drop(columns=["id_subsistema", "nom_subsistema"])
    return area_rows.astype(float)
