#!/usr/bin/env bash
#
# Batch-download ERA5 cutouts through the QC-gated, cached `download_cutout` rule,
# one target at a time (CDS runs one job per user), continuing past failures and
# QC-summarising at the end.
#
# Crash-resilient: cutouts already on disk / in the cache are skipped, so simply
# re-running resumes where an interrupted run (or a VM crash) left off — the
# expensive downloads live in `cutouts/cache/`, not in this script.
#
# Usage:
#   bash workflow/scripts/res_cf/download_cutouts_batch.sh                 # default target set
#   bash workflow/scripts/res_cf/download_cutouts_batch.sh cutouts/de_20240101_20241231.nc ...
#   nohup bash workflow/scripts/res_cf/download_cutouts_batch.sh &         # detached, survives logout
#
# Env overrides: SMK / PY (default to `snakemake` / `python` on PATH — activate
# the project env first, or point these at its binaries).
set -o pipefail

# Repo root = three levels up from workflow/scripts/res_cf/.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT" || exit 1

SMK="${SMK:-snakemake}"
PY="${PY:-python}"
LOG="logs/overnight_downloads.log"

# Default target set: 2024 (breadth-first), then 2025, then the 2023 baseline
# fill. de_2023 / aus_2025 already ship / exist, so they are not listed.
DEFAULT_TARGETS=(
  cutouts/de_20240101_20241231.nc
  cutouts/fr_20240101_20241231.nc
  cutouts/es_20240101_20241231.nc
  cutouts/aus_20240101_20241231.nc
  cutouts/bra_20240101_20241231.nc
  cutouts/de_20250101_20251231.nc
  cutouts/fr_20250101_20251231.nc
  cutouts/es_20250101_20251231.nc
  cutouts/bra_20250101_20251231.nc
  cutouts/fr_20230101_20231231.nc
  cutouts/es_20230101_20231231.nc
  cutouts/aus_20230101_20231231.nc
  cutouts/bra_20230101_20231231.nc
)
# Targets from CLI args if given, else the default set.
if [ "$#" -gt 0 ]; then TARGETS=("$@"); else TARGETS=("${DEFAULT_TARGETS[@]}"); fi

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" | tee -a "$LOG"; }

log "================ DOWNLOAD RUN START (resume-safe) ================"
log "targets: ${#TARGETS[@]}"
for t in "${TARGETS[@]}"; do
  if [ -f "$t" ]; then log "SKIP (exists): $t"; continue; fi
  log ">>> START $t"
  start=$(date +%s)
  "$SMK" --snakefile workflow/Snakefile --workflow-profile none \
      --cores 4 --keep-going --rerun-incomplete --nolock \
      --rerun-triggers mtime -- "$t" >> "$LOG" 2>&1
  rc=$?
  dur=$(( ($(date +%s) - start) / 60 ))
  if [ $rc -eq 0 ] && [ -f "$t" ]; then
    log "<<< DONE  $t  (${dur} min)  size=$(du -h "$t" | cut -f1)"
  else
    log "<<< FAIL  $t  rc=$rc (${dur} min) — continuing"
  fi
done

log "================ FINAL QC SUMMARY ================"
PYTHONPATH=workflow "$PY" - <<'PYEOF' 2>&1 | grep -viE "SerializationWarning|FutureWarning|warnings.warn" | tee -a "$LOG"
import re
from pathlib import Path
from common._cutout_qc import check_cutout
for f in sorted(Path("cutouts").glob("*.nc")):
    if "backup" in f.name or "spot" in f.name:
        continue
    m = re.match(r"[a-z]+_(\d{8})_(\d{8})\.nc$", f.name)
    if not m:
        continue
    try:
        r = check_cutout(f, m.group(1), m.group(2))
        print(f"{f.name}: {'OK' if r.ok else 'FAILED -> ' + '; '.join(r.problems)}")
    except Exception as e:
        print(f"{f.name}: ERROR opening -> {e!r}")
PYEOF

log "================ DOWNLOAD RUN END ================"
