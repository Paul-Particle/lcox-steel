"""Financial and sizing utilities for the DRI-hydrogen model."""


def annuity_factor(wacc: float, lifetime_years: float) -> float:
    """Capital recovery factor: annualises a lump-sum capex over a project life."""
    if lifetime_years <= 0:
        raise ValueError("lifetime_years must be > 0")
    if wacc == 0:
        return 1.0 / lifetime_years
    return (wacc * (1.0 + wacc) ** lifetime_years) / (
        (1.0 + wacc) ** lifetime_years - 1.0
    )


def dri_to_el_mw(
    dri_mt_per_year: float,
    h2_intensity_kg_per_t_dri: float,
    efficiency_kwh_per_kg: float,
    availability_target: float,
) -> float:
    """Size the electrolyser (MW electricity input) to supply a DRI plant continuously.

    `efficiency_kwh_per_kg` already encodes MWh-electricity per kg-H2, so the H2
    LHV doesn't appear here. `availability_target` is the fraction of 8760 h the
    plant runs (e.g. 1.0). Returns the nameplate MW of electricity input needed.
    """
    if dri_mt_per_year <= 0:
        raise ValueError("dri_mt_per_year must be > 0")
    if h2_intensity_kg_per_t_dri <= 0:
        raise ValueError("h2_intensity_kg_per_t_dri must be > 0")
    if efficiency_kwh_per_kg <= 0:
        raise ValueError("efficiency_kwh_per_kg must be > 0")
    if not (0 < availability_target <= 1.0):
        raise ValueError("availability_target must be in (0, 1]")

    h2_kg_per_year = dri_mt_per_year * 1_000_000.0 * h2_intensity_kg_per_t_dri
    elec_mwh_per_year = h2_kg_per_year * efficiency_kwh_per_kg / 1000.0
    return elec_mwh_per_year / (8760.0 * availability_target)
