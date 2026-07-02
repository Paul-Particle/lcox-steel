"""Structural QC for atlite ERA5 cutouts.

A cutout that downloaded partially, or that mixed final ERA5 with preliminary
ERA5T data, is the historical failure mode of the download stage: it either
carries an ``expver`` dimension/coordinate (ERA5 vs ERA5T) or has missing /
all-NaN timesteps. Neither raises on its own during ``cutout.prepare()``, so a
broken cutout can flow silently into the CF pipeline.

``check_cutout`` reads a finished cutout and returns a report of every problem
found; ``validate_cutout`` is the pipeline gate — it raises ``CutoutQCError`` so
the ``download_cutout`` rule fails and Snakemake discards the bad output.

The same function backs the standalone spot-check driver, so ad-hoc probes and
the pipeline apply identical criteria.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Physical fields that ERA5 defines everywhere over the cutout domain (land and
# sea, day and night). A fully-NaN timestep in any of these means the hour did
# not download. `height` is static (no time dim) so it is not time-checked.
KEY_TIME_VARS = ("influx_direct", "influx_diffuse", "wnd100m", "temperature", "runoff")


class CutoutQCError(ValueError):
    """Raised by validate_cutout when a cutout fails structural QC."""


@dataclass
class CutoutQCReport:
    path: Path
    dims: dict = field(default_factory=dict)
    time_start: pd.Timestamp | None = None
    time_end: pd.Timestamp | None = None
    n_steps: int = 0
    expected_steps: int | None = None
    n_duplicate_steps: int = 0
    n_nonhourly_steps: int = 0
    nan_fraction: dict = field(default_factory=dict)
    n_allnan_steps: dict = field(default_factory=dict)
    problems: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems

    def summary(self) -> str:
        head = f"{self.path.name}: {'OK' if self.ok else 'FAILED'}"
        lines = [head, f"  dims={self.dims}"]
        if self.time_start is not None:
            lines.append(
                f"  time {self.time_start} -> {self.time_end}  "
                f"steps={self.n_steps} expected={self.expected_steps} "
                f"dupes={self.n_duplicate_steps} non-1h={self.n_nonhourly_steps}"
            )
        for var, frac in self.nan_fraction.items():
            note = ""
            if self.n_allnan_steps.get(var):
                note = f"  <-- {self.n_allnan_steps[var]} fully-NaN timesteps"
            lines.append(f"  {var:16s} NaN {frac * 100:6.2f}%{note}")
        for p in self.problems:
            lines.append(f"  PROBLEM: {p}")
        return "\n".join(lines)


def _expected_hourly_steps(start_yyyymmdd: str, end_yyyymmdd: str) -> int:
    """Hourly step count for the inclusive [start 00:00, end 23:00] range."""
    start = pd.Timestamp(f"{start_yyyymmdd[:4]}-{start_yyyymmdd[4:6]}-{start_yyyymmdd[6:8]} 00:00")
    end = pd.Timestamp(f"{end_yyyymmdd[:4]}-{end_yyyymmdd[4:6]}-{end_yyyymmdd[6:8]} 23:00")
    return int((end - start) / pd.Timedelta(hours=1)) + 1


def check_cutout(
    path: str | Path,
    start_date: str | None = None,
    end_date: str | None = None,
    max_nan_fraction: float = 0.0,
) -> CutoutQCReport:
    """Inspect a cutout and return a CutoutQCReport (never raises on bad data).

    If ``start_date``/``end_date`` (YYYYMMDD) are given, the time axis is checked
    for exact coverage of [start 00:00, end 23:00]; otherwise only internal
    continuity is checked. ``max_nan_fraction`` is the tolerated NaN share per
    key field (fully-NaN timesteps always fail regardless).
    """
    path = Path(path)
    report = CutoutQCReport(path=path)
    ds = xr.open_dataset(path)
    try:
        report.dims = dict(ds.sizes)

        # ERA5 / ERA5T mixing marker — atlite strips expver on a clean download,
        # so its presence means preliminary data leaked through.
        for marker in ("expver", "number"):
            if marker in ds.sizes or marker in ds.coords:
                report.problems.append(
                    f"'{marker}' present (dim/coord) — ERA5/ERA5T mix or unprocessed download"
                )

        if "time" in ds.coords:
            t = pd.DatetimeIndex(pd.to_datetime(ds["time"].values)).sort_values()
            report.time_start, report.time_end = t[0], t[-1]
            report.n_steps = len(t)
            report.n_duplicate_steps = int(t.duplicated().sum())
            diffs = t.to_series().diff().dropna()
            report.n_nonhourly_steps = int((diffs != pd.Timedelta(hours=1)).sum())

            if report.n_duplicate_steps:
                report.problems.append(f"{report.n_duplicate_steps} duplicate timesteps")
            if report.n_nonhourly_steps:
                report.problems.append(f"{report.n_nonhourly_steps} non-hourly steps (gaps)")

            if start_date and end_date:
                report.expected_steps = _expected_hourly_steps(start_date, end_date)
                expected_start = pd.Timestamp(
                    f"{start_date[:4]}-{start_date[4:6]}-{start_date[6:8]} 00:00"
                )
                expected_end = pd.Timestamp(
                    f"{end_date[:4]}-{end_date[4:6]}-{end_date[6:8]} 23:00"
                )
                if report.n_steps != report.expected_steps:
                    report.problems.append(
                        f"{report.n_steps} timesteps, expected {report.expected_steps}"
                    )
                if report.time_start != expected_start:
                    report.problems.append(
                        f"starts at {report.time_start}, expected {expected_start}"
                    )
                if report.time_end != expected_end:
                    report.problems.append(
                        f"ends at {report.time_end}, expected {expected_end}"
                    )
        else:
            report.problems.append("no 'time' coordinate")

        for var in KEY_TIME_VARS:
            if var not in ds.data_vars:
                report.problems.append(f"missing expected variable '{var}'")
                continue
            da = ds[var]
            report.nan_fraction[var] = float(np.isnan(da.values).mean())
            if "time" in da.dims:
                other_dims = [d for d in da.dims if d != "time"]
                allnan = int(np.isnan(da).all(dim=other_dims).sum())
                report.n_allnan_steps[var] = allnan
                if allnan:
                    report.problems.append(f"'{var}' has {allnan} fully-NaN timesteps")
            if report.nan_fraction[var] > max_nan_fraction:
                report.problems.append(
                    f"'{var}' NaN fraction {report.nan_fraction[var]:.4f} > {max_nan_fraction}"
                )
    finally:
        ds.close()
    return report


def validate_cutout(
    path: str | Path,
    start_date: str | None = None,
    end_date: str | None = None,
    max_nan_fraction: float = 0.0,
) -> CutoutQCReport:
    """Pipeline gate: run check_cutout and raise CutoutQCError if it failed."""
    report = check_cutout(path, start_date, end_date, max_nan_fraction)
    if not report.ok:
        raise CutoutQCError(f"Cutout failed QC:\n{report.summary()}")
    return report
