"""Where each sweep's ordering changes, read off the solved runs.

Linear interpolation between the two bracketing sweep points. The point of doing
it from the reports rather than from the cost equation is that the equation holds
the rest of the system fixed, and the optimiser does not: it re-sizes everything
around the input that moved.
"""
import re
from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parents[1] / "results"

SWEEPS = {
    "ew-capex": ("ew-eaf", "moe-eaf", "EW capex, EUR/(t.yr)"),
    "gas": ("ng-dri-eaf", "h2-dri-eaf", "gas price, EUR/MWh"),
}

for prefix, (route_a, route_b, label) in SWEEPS.items():
    points = {}
    for path in sorted(RESULTS.glob(f"report_{prefix}-*.csv")):
        value = int(re.search(rf"{prefix}-(\d+)", path.stem).group(1))
        report = pd.read_csv(path, index_col=0).T
        report["lcos"] = pd.to_numeric(report["lcos_eur_per_t"], errors="coerce")
        for area in report["area"].unique():
            a = report[(report["area"] == area) & (report["route"] == route_a)]["lcos"]
            b = report[(report["area"] == area) & (report["route"] == route_b)]["lcos"]
            if len(a) and len(b):
                points.setdefault(area, {})[value] = float(a.iloc[0]) - float(b.iloc[0])

    print(f"=== {prefix}: {route_a} minus {route_b}, EUR/t ===")
    print(f"{'area':6s} " + " ".join(f"{v:>8d}" for v in sorted(next(iter(points.values()))))
          + "   crossing")
    for area, series in sorted(points.items()):
        values = sorted(series)
        row = " ".join(f"{series[v]:>8.1f}" for v in values)
        crossing = "—"
        for lo, hi in zip(values, values[1:]):
            if (series[lo] < 0) != (series[hi] < 0):
                frac = -series[lo] / (series[hi] - series[lo])
                crossing = f"{lo + frac * (hi - lo):.0f}"
                break
        print(f"{area:6s} {row}   {crossing}")
    print(f"       ({label}; negative means {route_a} is cheaper)\n")
