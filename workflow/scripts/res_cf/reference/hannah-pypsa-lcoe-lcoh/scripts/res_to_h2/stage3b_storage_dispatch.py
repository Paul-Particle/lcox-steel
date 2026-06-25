import numpy as np


def simulate_battery_dispatch(
    surplus_mwh: np.ndarray,
    deficit_mwh: np.ndarray,
    storage_capacity_mwh: float,
    eta_rt: float,
    initial_soc_mode: str = "empty",
):
    """
    Chronological battery dispatch (no power limits, no optimisation)

    Returns:
        dict with:
            unserved_energy_mwh
            curtailment_mwh
            max_soc_mwh
            throughput_mwh (discharged energy)
    """

    # efficiencies
    eta_c = np.sqrt(eta_rt)
    eta_d = np.sqrt(eta_rt)

    n = len(surplus_mwh)

    # initial SOC
    if initial_soc_mode == "empty":
        soc = 0.0
    elif initial_soc_mode == "full":
        soc = storage_capacity_mwh
    else:
        raise ValueError("initial_soc_mode must be 'empty' or 'full'")

    max_soc = soc
    total_unserved = 0.0
    total_curtailment = 0.0
    total_throughput = 0.0  # discharged energy

    for t in range(n):

        surplus = surplus_mwh[t]
        deficit = deficit_mwh[t]

        # -----------------
        # Case 1: surplus → charge
        # -----------------
        if surplus > 0:

            # max input energy limited by SOC headroom
            max_input_due_to_capacity = (storage_capacity_mwh - soc) / eta_c

            charge_input = min(surplus, max_input_due_to_capacity)

            soc += charge_input * eta_c

            curtailment = surplus - charge_input
            total_curtailment += curtailment

        # -----------------
        # Case 2: deficit → discharge
        # -----------------
        elif deficit > 0:

            # max deliverable energy from SOC
            max_output_due_to_soc = soc * eta_d

            discharge_output = min(deficit, max_output_due_to_soc)

            soc -= discharge_output / eta_d

            unserved = deficit - discharge_output
            total_unserved += unserved

            total_throughput += discharge_output

        # update max SOC
        if soc > max_soc:
            max_soc = soc

    return {
        "unserved_energy_mwh": total_unserved,
        "curtailment_mwh": total_curtailment,
        "max_soc_mwh": max_soc,
        "throughput_mwh": total_throughput,
    }

import numpy as np


def find_min_storage_capacity(
    surplus_mwh: np.ndarray,
    deficit_mwh: np.ndarray,
    allowed_uncovered_mwh: float,
    eta_rt: float,
    initial_soc_mode: str = "empty",
    tolerance_mwh: float = 1.0,
    max_iter: int = 60,
):
    """
    Find the minimum battery energy capacity [MWh] such that:

        annual_unserved_energy_mwh <= allowed_uncovered_mwh

    Uses binary search.

    Returns:
        dict with:
            feasible
            storage_capacity_required_mwh
            unserved_energy_mwh
            curtailment_mwh
            max_soc_mwh
            throughput_mwh
    """

    if allowed_uncovered_mwh < 0:
        raise ValueError("allowed_uncovered_mwh must be >= 0")
    if tolerance_mwh <= 0:
        raise ValueError("tolerance_mwh must be > 0")
    if max_iter <= 0:
        raise ValueError("max_iter must be > 0")

    # Lower bound: no storage
    lo = 0.0

    # Upper bound: safe first guess
    hi = float(np.sum(surplus_mwh))

    # Edge case: no surplus at all
    if hi <= 0:
        test = simulate_battery_dispatch(
            surplus_mwh=surplus_mwh,
            deficit_mwh=deficit_mwh,
            storage_capacity_mwh=0.0,
            eta_rt=eta_rt,
            initial_soc_mode=initial_soc_mode,
        )
        return {
            "feasible": test["unserved_energy_mwh"] <= allowed_uncovered_mwh,
            "storage_capacity_required_mwh": 0.0,
            **test,
        }

    # First test: even with very large storage, is target achievable?
    hi_test = simulate_battery_dispatch(
        surplus_mwh=surplus_mwh,
        deficit_mwh=deficit_mwh,
        storage_capacity_mwh=hi,
        eta_rt=eta_rt,
        initial_soc_mode=initial_soc_mode,
    )

    best_test = hi_test

    if hi_test["unserved_energy_mwh"] > allowed_uncovered_mwh:
        return {
            "feasible": False,
            "storage_capacity_required_mwh": np.nan,
            **hi_test,
        }

    # Binary search
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)

        mid_test = simulate_battery_dispatch(
            surplus_mwh=surplus_mwh,
            deficit_mwh=deficit_mwh,
            storage_capacity_mwh=mid,
            eta_rt=eta_rt,
            initial_soc_mode=initial_soc_mode,
        )

        if mid_test["unserved_energy_mwh"] <= allowed_uncovered_mwh:
            hi = mid
            best_test = mid_test
        else:
            lo = mid

        if (hi - lo) <= tolerance_mwh:
            break

    return {
        "feasible": True,
        "storage_capacity_required_mwh": hi,
        **best_test,
    }