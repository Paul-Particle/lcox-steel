"""Validate the ENTSO-E grid pipeline against the energy-charts.info reference.

energy-charts.info (Fraunhofer ISE) publishes the same German day-ahead price,
generation-by-type, and load that our pipeline derives from ENTSO-E, so it is an
independent cross-check for what had until now only been eyeballed on the
dashboard (issue #20). Both draw on ENTSO-E, so we expect very close agreement:
the thresholds below carry roughly 2x headroom over the agreement measured on
2023 data.

The comparison is on the hourly overlap. energy-charts applies its own revisions
and gap-filling, so a small positive bias (our value slightly above
energy-charts) is expected and does not fail the gate.
"""

from dataclasses import dataclass

import pandas as pd
import pytest

# (min correlation, max relative MAE %) per metric.
# Calibrated on DE_LU 2023; measured agreement in parentheses.
THRESHOLDS = {
    "price":         (0.99, 5.0),    # measured: corr 1.000, relMAE 0.0%
    "solar":         (0.99, 8.0),    # measured: corr 1.000, relMAE 3.8%
    "wind_onshore":  (0.99, 8.0),    # measured: corr 1.000, relMAE 3.1%
    "wind_offshore": (0.99, 10.0),   # measured: corr 1.000, relMAE 0.1%
    "load":          (0.995, 4.0),   # measured: corr 1.000, relMAE 1.0%
    "wind":          (0.99, 8.0),    # onshore + offshore
    "res":           (0.99, 8.0),    # wind + solar
}


@dataclass
class Agreement:
    n: int
    mae: float
    bias: float
    corr: float
    rel_mae_pct: float


def _agreement(ours: pd.Series, ref: pd.Series) -> Agreement:
    """MAE / bias / correlation / relative MAE on the hourly-mean overlap."""
    ours_h = ours.resample("1h").mean()
    ref_h = ref.resample("1h").mean()
    common = ours_h.index.intersection(ref_h.index)
    a, b = ours_h.loc[common], ref_h.loc[common]
    mask = a.notna() & b.notna()
    a, b = a[mask], b[mask]
    mae = (a - b).abs().mean()
    return Agreement(
        n=len(a),
        mae=mae,
        bias=(a - b).mean(),
        corr=a.corr(b),
        rel_mae_pct=100 * mae / b.abs().mean(),
    )


@pytest.mark.parametrize("metric", list(THRESHOLDS))
def test_pipeline_matches_energycharts(metric, pipeline_full_de_2023, energycharts_ref):
    # The reference stores primitives; derive the aggregates the same way the
    # pipeline does so the comparison is apples-to-apples.
    ref = energycharts_ref.assign(
        wind=lambda d: d["wind_onshore"] + d["wind_offshore"],
        res=lambda d: d["wind_onshore"] + d["wind_offshore"] + d["solar"],
    )

    agreement = _agreement(pipeline_full_de_2023[metric], ref[metric])
    min_corr, max_rel_mae = THRESHOLDS[metric]

    assert agreement.n > 8000, f"{metric}: only {agreement.n} overlapping hours"
    assert agreement.corr >= min_corr, (
        f"{metric}: corr {agreement.corr:.4f} < {min_corr} "
        f"(bias {agreement.bias:.1f}, relMAE {agreement.rel_mae_pct:.1f}%)"
    )
    assert agreement.rel_mae_pct <= max_rel_mae, (
        f"{metric}: relMAE {agreement.rel_mae_pct:.1f}% > {max_rel_mae}% "
        f"(corr {agreement.corr:.4f}, bias {agreement.bias:.1f})"
    )
