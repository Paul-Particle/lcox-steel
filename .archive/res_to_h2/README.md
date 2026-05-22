# RES → H2 optimisation (Dedicated intermittent system)

This folder contains a simplified techno-economic model that computes **levelised cost of hydrogen (LCOH)** from dedicated renewable energy (RES) supply using **hourly capacity factor (CF) time series** (e.g. from the Atlite pipeline).

The model is designed for:

- Electrolyser sizing  
- Cross-country comparison  
- Green DRI modelling  
- System-level RES → H2 analysis  

All scenarios are defined in `config_hannah.yaml`.

All scenarios use identical modelling logic to ensure cross-scenario comparability.

---

# Outputs (authoritative)

For each scenario defined under `scenarios:` in the YAML file, the script produces:

### Sweep file

`results/<scenario_name>_res_to_h2_sweep.csv`

This file contains the **full optimisation sweep across RES-to-electrolyser ratios**.

### Summary file

`results/<scenario_name>_res_to_h2_summary.csv`

This file contains **design-point metrics** for:

- cost-optimal intermittent design
- oversizing-only firm reliability design (if achievable)

---

# Sweep output columns

Each row corresponds to one tested **RES-to-electrolyser ratio**.

Core system variables

- `res_to_el_ratio`
- `res_capacity_mw`
- `electrolyser_mw`
- `res_generation_mwh`
- `elec_used_mwh`
- `curtailed_mwh`

Hydrogen production

- `h2_kg_per_year`

Cost metrics

- `annual_res_cost_eur`
- `annual_electrolyser_cost_eur`
- `annual_storage_cost_eur`
- `lcoh_eur_per_kg`
- `lcoh_with_storage_eur_per_kg`

System performance

- `el_full_load_hours`
- `el_utilisation`
- `curtailment_rate`
- `firm_reliability`

Renewable electricity cost diagnostics

- `lcoe_res_prod_eur_per_mwh`
- `lcoe_res_used_eur_per_mwh`

Storage proxy diagnostics

- `annual_deficit_mwh`
- `storage_capacity_required_mwh`
- `storage_capacity_required_uncapped_mwh`
- `storage_duration_cap_mwh`
- `storage_duration_cap_binding_bool`
- `annual_unserved_mwh_after_storage_cap`

Energy service diagnostics

- `energy_served_fraction_no_storage`
- `energy_served_fraction_proxy`
- `energy_unserved_mwh_proxy`

---

# Method summary

## 1. Hourly dispatch logic

The system assumes:

- Dedicated RES supply  
- No grid imports  
- Perfect electrolyser flexibility within nameplate limits  

At each hourly timestep `t`:

P_res(t) = CF(t) × P_res_installed  
P_used(t) = min(P_res(t), P_el)  
P_curtailed(t) = max(P_res(t) − P_el, 0)

Where:

- `CF(t)` = hourly capacity factor  
- `P_res_installed` = installed renewable capacity  
- `P_el` = electrolyser capacity

Annual electricity:

E_used = Σ P_used(t) × Δt  
E_curtailed = Σ P_curtailed(t) × Δt  
E_res_generation = Σ P_res(t) × Δt  

with Δt = 1 hour.

---

## 2. Hydrogen production

H2_kg = (E_used × 1000) / efficiency_kwh_per_kg

Where:

- `H2_kg` = annual hydrogen production  
- `efficiency_kwh_per_kg` = electrolyser electricity intensity

---

## 3. Annualised cost calculation

Capital recovery factor:

CRF = w × (1 + w)^n / ((1 + w)^n − 1)

Where:

- `w` = weighted average cost of capital  
- `n` = lifetime

Annual cost components:

C_res = annualised renewable cost  
C_el  = annualised electrolyser cost  

Total system cost:

C_total = C_res + C_el

---

## 4. Levelised cost of hydrogen

LCOH = C_total / H2_kg

Units: €/kg

---

# System performance metrics

Electrolyser utilisation:

EL_utilisation = E_used / (P_el × 8760)

Curtailment rate:

Curtailment_rate = E_curtailed / E_res_generation

Electrolyser full load hours:

FLH = E_used / P_el

---

# Firming diagnostics (oversizing-only baseline)

The model reports a **firm reliability metric** defined as:

firm_reliability = share of hours where P_res(t) ≥ P_el

This measures how often the electrolyser can operate without deficit.

The optimisation sweep identifies the first RES-to-electrolyser ratio satisfying:

firm_reliability ≥ firm_power_fraction

Default target:

firm_power_fraction = 0.95

If the target cannot be reached within the sweep range, the model reports that firming is not achievable within the tested ratios.

---

# Stage 3A storage proxy (energy-based firming)

A simplified storage proxy estimates storage required to achieve a target **annual energy service level**.

This stage **does not simulate chronological storage dispatch**.

Instead it approximates storage requirements using cumulative surplus and deficit energy.

---

## Hourly deficit and surplus

Deficit(t) = max(P_el − P_res(t), 0)  
Surplus(t) = max(P_res(t) − P_el, 0)

Energy equivalents:

Deficit_energy(t) = Deficit(t) × Δt  
Surplus_energy(t) = Surplus(t) × Δt  

Annual deficit energy:

Annual_deficit = Σ Deficit_energy(t)

---

## Allowed annual deficit

A configurable fraction of annual energy demand may remain uncovered.

allowed_energy_deficit_fraction = f

Default:

f = 0.05

Annual baseload energy:

E_baseload = P_el × 8760

Allowed uncovered energy:

E_allowed_uncovered = f × E_baseload

Energy that must be covered by storage:

E_to_cover = Annual_deficit − E_allowed_uncovered

If Annual_deficit ≤ E_allowed_uncovered, no storage is required.

---

## Storage energy requirement proxy

Net system energy balance:

Net(t) = Surplus_energy(t) − Deficit_energy(t)

Cumulative balance:

Cumulative(t) = cumulative_sum(Net(t))

Proxy storage capacity requirement:

Storage_required_uncapped = max(Cumulative) − min(Cumulative)

---

## Scaling to required coverage

Scaling_factor = E_to_cover / Annual_deficit

Storage_required_scaled = Storage_required_uncapped × Scaling_factor

---

## Efficiency adjustment

Storage_required = Storage_required_scaled / η

Where η is storage round-trip efficiency.

---

## Storage duration cap (battery proxy)

To avoid unrealistic seasonal storage assumptions a maximum duration constraint is applied.

max_storage_energy = P_el × max_storage_duration_hours

Final storage capacity:

Storage_capacity = min(Storage_required, max_storage_energy)

The model reports whether the cap binds.

---

## Residual unmet energy

If the duration cap binds the model reports remaining unmet energy:

annual_unserved_mwh_after_storage_cap

This indicates the battery duration assumption is insufficient to meet the energy service target.

---

# Energy service diagnostics

Energy served without storage:

energy_served_fraction_no_storage = E_used / (P_el × 8760)

Energy served with proxy storage:

energy_served_fraction_proxy = 1 − allowed_energy_deficit_fraction

Remaining uncovered energy:

energy_unserved_mwh_proxy = allowed_energy_deficit_fraction × P_el × 8760

---

# DRI-scale sizing

Electrolyser capacity derived from DRI output:

H2_annual = dri_mt_per_year × 1,000,000 × h2_intensity_kg_per_t_dri

Electricity_annual = H2_annual × efficiency_kwh_per_kg / 1000

MW_el = Electricity_annual / (8760 × availability_target)

Where:

- `dri_mt_per_year` = DRI output  
- `h2_intensity_kg_per_t_dri` = hydrogen demand per tonne DRI  
- `availability_target` = design supply fraction

Important distinction:

availability_target determines installed electrolyser capacity.

Actual utilisation is determined by RES variability and oversizing.

---

# Diagnostics script

`scripts/res_to_h2/diagnostics_res_to_h2.py`

This script visualises hourly system behaviour for selected design points.

For each scenario it plots:

- cost-optimal intermittent design
- firm reliability design (if achievable)

Time windows analysed:

- January representative week  
- January representative day  
- July representative week  
- July representative day  

Plots produced:

1. Power time series  
   P_RES(t), P_used(t), and P_EL

2. Deficit duration curve  
   Sorted hourly deficits

3. Cumulative deficit energy curve  
   Compared against allowed deficit threshold

Plots are written to:

results/plots/

---

# Scripts (pipeline overview)

`scripts/res_to_h2/res_to_h2_logic.py`

Contains dispatch simulation, cost calculations, optimisation sweep, and storage proxy logic.

`scripts/res_to_h2/run_res_to_h2_optimisation.py`

Main execution script:

- loads configuration
- loads Atlite CF data
- computes DRI-scale electrolyser size
- runs optimisation
- writes output CSV files

Run with:

python scripts/res_to_h2/run_res_to_h2_optimisation.py

---

# Assumptions

- Linear CAPEX scaling  
- Perfect electrolyser flexibility  
- No minimum load constraint  
- No degradation  
- Dedicated RES supply  
- No grid imports  
- No hybrid dispatch  

Storage proxy assumptions:

- energy-only proxy  
- no chronological dispatch  
- configurable annual deficit allowance  
- battery-style duration cap

---

# Future extensions

- Stage 3B chronological storage dispatch model  
- Hybrid RES + grid supply  
- Electrolyser minimum load constraint  
- Seasonal storage technologies  
- Degradation modelling 