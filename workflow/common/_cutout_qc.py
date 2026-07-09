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


@dataclass(frozen=True)
class NaNAllowance:
    """How much NaN a single variable may contain before QC fails.

    max_total        total NaN cells tolerated across the whole field.
    max_consecutive  longest tolerated run of consecutive timesteps in which the
                     variable contains *any* NaN cell.

    A *fully*-NaN timestep (an entire missing hour) always fails regardless of
    this allowance — see `check_cutout`.
    """

    max_total: int = 0
    max_consecutive: int = 0


# Per-variable NaN tolerance applied by `check_cutout`. Empty by default: any
# variable not listed here is held to a hard zero (NaNAllowance() == 0 total,
# 0 consecutive), so a single NaN anywhere fails QC. This strictness is
# deliberate — every cutout downloaded so far is flawless and we do not yet know
# whether *any* NaN is normal for ERA5 over these domains.
#
# If a source of slightly-degraded-but-still-acceptable data ever appears (e.g.
# a variable that legitimately carries a few isolated NaNs), relax it explicitly
# here rather than loosening the global check, e.g.:
#
#     DEFAULT_NAN_ALLOWANCE = {
#         "runoff": NaNAllowance(max_total=48, max_consecutive=3),
#     }
#
# A fully-NaN hour still fails even for a variable listed here.
DEFAULT_NAN_ALLOWANCE: dict[str, NaNAllowance] = {}


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
    n_nan: dict = field(default_factory=dict)
    n_allnan_steps: dict = field(default_factory=dict)
    max_consecutive_nan_steps: dict = field(default_factory=dict)
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
            parts = [f"NaN {frac * 100:6.2f}%", f"n={self.n_nan.get(var, 0)}"]
            if self.n_allnan_steps.get(var):
                parts.append(f"{self.n_allnan_steps[var]} fully-NaN hrs")
            if self.max_consecutive_nan_steps.get(var):
                parts.append(f"max {self.max_consecutive_nan_steps[var]} consec")
            lines.append(f"  {var:16s} " + "  ".join(parts))
        for p in self.problems:
            lines.append(f"  PROBLEM: {p}")
        return "\n".join(lines)


def _expected_hourly_steps(start_yyyymmdd: str, end_yyyymmdd: str) -> int:
    """Hourly step count for the inclusive [start 00:00, end 23:00] range."""
    start = pd.Timestamp(f"{start_yyyymmdd[:4]}-{start_yyyymmdd[4:6]}-{start_yyyymmdd[6:8]} 00:00")
    end = pd.Timestamp(f"{end_yyyymmdd[:4]}-{end_yyyymmdd[4:6]}-{end_yyyymmdd[6:8]} 23:00")
    return int((end - start) / pd.Timedelta(hours=1)) + 1


def _longest_true_run(mask: np.ndarray) -> int:
    """Longest run of consecutive True values in a 1-D boolean array."""
    best = run = 0
    for flag in mask:
        run = run + 1 if flag else 0
        if run > best:
            best = run
    return best


def check_cutout(
    path: str | Path,
    start_date: str | None = None,
    end_date: str | None = None,
    nan_allowance: dict[str, NaNAllowance] | None = None,
) -> CutoutQCReport:
    """Inspect a cutout and return a CutoutQCReport (never raises on bad data).

    If ``start_date``/``end_date`` (YYYYMMDD) are given, the time axis is checked
    for exact coverage of [start 00:00, end 23:00]; otherwise only internal
    continuity is checked.

    ``nan_allowance`` maps a variable name to how much NaN it may carry (see
    ``NaNAllowance``); it defaults to ``DEFAULT_NAN_ALLOWANCE`` (empty), which
    holds every key variable to zero NaN. A fully-NaN timestep (a missing hour)
    always fails regardless of the allowance.
    """
    path = Path(path)
    allowances = DEFAULT_NAN_ALLOWANCE if nan_allowance is None else nan_allowance
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
            allowance = allowances.get(var, NaNAllowance())

            total_nan = int(np.isnan(da.values).sum())
            report.n_nan[var] = total_nan
            report.nan_fraction[var] = total_nan / da.size if da.size else 0.0

            if "time" in da.dims:
                other_dims = [d for d in da.dims if d != "time"]
                isnan = np.isnan(da)
                allnan_mask = isnan.all(dim=other_dims).values
                anynan_mask = isnan.any(dim=other_dims).values
                report.n_allnan_steps[var] = int(allnan_mask.sum())
                report.max_consecutive_nan_steps[var] = _longest_true_run(anynan_mask)

                # A whole missing hour always fails, regardless of the allowance.
                if report.n_allnan_steps[var]:
                    report.problems.append(
                        f"'{var}' has {report.n_allnan_steps[var]} fully-NaN timesteps (missing hours)"
                    )
                if report.max_consecutive_nan_steps[var] > allowance.max_consecutive:
                    report.problems.append(
                        f"'{var}' has {report.max_consecutive_nan_steps[var]} consecutive NaN "
                        f"timesteps (allowed {allowance.max_consecutive})"
                    )

            if total_nan > allowance.max_total:
                report.problems.append(
                    f"'{var}' has {total_nan} NaN values (allowed {allowance.max_total})"
                )
    finally:
        ds.close()
    return report


def validate_cutout(
    path: str | Path,
    start_date: str | None = None,
    end_date: str | None = None,
    nan_allowance: dict[str, NaNAllowance] | None = None,
) -> CutoutQCReport:
    """Pipeline gate: run check_cutout and raise CutoutQCError if it failed."""
    report = check_cutout(path, start_date, end_date, nan_allowance)
    if not report.ok:
        raise CutoutQCError(f"Cutout failed QC:\n{report.summary()}")
    return report
