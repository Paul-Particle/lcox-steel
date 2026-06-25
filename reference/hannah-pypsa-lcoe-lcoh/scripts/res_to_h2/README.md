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

`results/<scenario_name>_stage3b_sweep.csv`
`results/<scenario_name>_res_to_h2_sweep.csv`

This file contains the **full optimisation sweep across RES-to-electrolyser ratios**.
Stage 3B outputs (`*_stage3b_sweep.csv`) are the **authoritative results** for all analyses involving storage and firm hydrogen supply.

### Summary file (legacy / optional)

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
- `elec_served_no_storage_mwh`
- `elec_served_with_storage_mwh`

Hydrogen production

- `h2_kg_per_year`
- `h2_kg_no_storage_per_year`
- `h2_kg_with_storage_per_year`

Cost metrics

- `annual_res_cost_eur`
- `annual_electrolyser_cost_eur`
- `annual_storage_cost_eur`
- `lcoh_eur_per_kg`→ **intermittent system (no storage)**
- `lcoh_with_storage_eur_per_kg`→ **Stage 3B result (authoritative)**

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

Without storage:

H2_no_storage = (E_used × 1000) / efficiency_kwh_per_kg

With Stage 3B storage:

H2_with_storage = (E_served_stage3b × 1000) / efficiency_kwh_per_kg

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

Note:

Stage 3B diagnostics plots are based on short-window re-simulation with reset initial SOC.
They are illustrative of dispatch behaviour but do not represent the exact annual SOC trajectory.

---

# Scripts (pipeline overview)

`scripts/res_to_h2/res_to_h2_logic.py`

Contains shared techno-economic logic including:

- annualised cost functions
- DRI-scale electrolyser sizing
- intermittent RES → electrolyser sweep logic
- Stage 3A storage proxy logic

`scripts/res_to_h2/stage3b_storage_dispatch.py`

Contains Stage 3B chronological battery dispatch logic, including:

- hourly charge / discharge simulation
- annual unserved-energy calculation
- binary search for minimum storage capacity satisfying the allowed annual deficit constraint

`scripts/res_to_h2/run_stage3b_sweep.py`

Runs the Stage 3B optimisation sweep across RES-to-electrolyser ratios for a given scenario and returns:

- no-storage metrics
- Stage 3B storage sizing outputs
- with-storage LCOH outputs
- Stage 3B cost-component outputs

`scripts/res_to_h2/run_res_to_h2_stage3b.py`

Main execution script for Stage 3B:

- loads configuration
- loads Atlite CF data
- computes DRI-scale electrolyser size
- runs the Stage 3B sweep
- writes `*_stage3b_sweep.csv` outputs

Run with:

python scripts/res_to_h2/run_res_to_h2_stage3b.py

`scripts/res_to_h2/run_res_to_h2_optimisation.py`

Legacy / baseline optimisation script for:

- intermittent no-storage case
- oversizing-only firming baseline
- Stage 3A storage proxy outputs

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

# Stage 3B storage model (chronological dispatch-based firming)

Stage 3B extends the Stage 3A proxy by introducing a **chronological storage dispatch simulation**.  
This allows explicit modelling of **hourly charge/discharge behaviour**, instead of relying on cumulative energy approximations.

This stage provides the chronologically consistent and cost-consistent LCOH including storage,
using electricity actually served to the electrolyser after storage dispatch.

- `lcoh_with_storage_eur_per_kg`
- Stage 3B sweep outputs (`*_stage3b_sweep.csv`)
- LCOH breakdown plots (RES / electrolyser / storage)

---

## Conceptual difference vs Stage 3A

| Feature                     | Stage 3A (proxy)                    | Stage 3B (implemented)              |
|----------------------------|------------------------------------|-------------------------------------|
| Time resolution            | Aggregated (cumulative)            | Hourly chronological                |
| Storage dispatch           | No                                 | Yes                                 |
| Energy shifting realism    | Approximate                        | Explicit                            |
| Power constraints          | Not represented                    | Implicit via timestep balance       |
| Cost consistency           | Partial                            | Full (used in LCOH)                 |

---

## 1. Hourly storage dispatch logic

At each timestep `t`:

Residual balance:

Net(t) = P_res(t) − P_el

Storage state evolution:

If Net(t) > 0 (surplus):
- Charge storage (subject to capacity limit)

If Net(t) < 0 (deficit):
- Discharge storage (subject to available energy)

State of charge (SOC):

Charging efficiency: η_c = √η_rt  
Discharging efficiency: η_d = √η_rt  

SOC(t) = SOC(t−1) + charge(t) × η_c − discharge(t) / η_d

Where:

- `η` = storage round-trip efficiency (applied on charging side)
- SOC is bounded:

0 ≤ SOC(t) ≤ Storage_capacity

---

## 2. Storage sizing

Storage capacity is determined via binary search on chronological dispatch.

The model finds the minimum storage capacity such that:

annual_unserved_energy_mwh ≤ allowed_uncovered_mwh

This ensures the annual energy service constraint is satisfied under full hourly dispatch.

## 3. Energy served with storage

The dispatch simulation tracks:

- energy served to electrolyser
- unmet demand after storage

Residual unmet energy:

unserved(t) = max(deficit(t) − discharge(t), 0)

Annual unserved energy:

E_unserved_stage3b = Σ unserved(t)

---

## 4. Feasibility criterion

A scenario is considered **storage-feasible** if:

E_unserved_stage3b ≤ allowed_energy_deficit_fraction × (P_el × 8760)

This corresponds to the same energy service target as Stage 3A, but now validated chronologically.

Output column:

- `storage_feasible` (boolean)

Only feasible points are used for:

- `lcoh_with_storage_eur_per_kg`
- optimal system identification

---

## 5. Storage cost integration

Annual storage cost:

C_storage = annualised storage CAPEX + fixed OPEX

This is included in total system cost:

C_total_with_storage = C_res + C_el + C_storage

Resulting LCOH:

With storage, hydrogen output is based on electricity actually served to the electrolyser after storage dispatch.

Annual served electricity:

E_served_stage3b = (P_el × 8760) − E_unserved_stage3b

Hydrogen production:

H2_stage3b = (E_served_stage3b × 1000) / efficiency_kwh_per_kg

Levelised cost:

LCOH_with_storage = (C_res + C_el + C_storage) / H2_stage3b

## 6. Key outputs added in Stage 3B

Additional columns in sweep output:

- `storage_feasible`
- `lcoh_with_storage_eur_per_kg`

These represent the **fully costed and chronologically validated system**.

---

## 7. Interpretation

Stage 3B reveals important system effects not captured in Stage 3A:

- Chronological mismatch between surplus and deficit
- Inefficient utilisation of storage due to timing
- Increased storage requirements vs proxy in some cases
- Potential infeasibility despite sufficient annual energy balance

---

## 8. Important modelling implications

- Best-site RES does **not necessarily reduce LCOH**:
  - Higher CF → higher utilisation
  - But also → higher curtailment or storage needs depending on profile shape

- Storage costs can dominate total LCOH in:
  - low-CF regions
  - highly intermittent solar-dominated systems

- The model remains:
  - energy-constrained (not power-constrained storage)
  - battery-like (duration cap enforced)
  - deterministic (single-year CF input)

# RES mix scenarios (multi-technology systems)

The model supports combined renewable systems via:

- `tech: "res_mix"`
- `res_mix: {tech: weight}`
- `mix_anchor_tech` (for best-site scenarios)

## RES mix definition

The hourly RES profile is constructed as:

CF_mix(t) = Σ (w_i × CF_i(t))

Where:
- w_i = technology weights (sum to 1)
- CF_i(t) = hourly capacity factor of technology i

## Co-location logic (best-site scenarios)

For `variant: "bestsite_p95"`:

- A single anchor technology defines the location
- All technologies are evaluated at the same location basis

This uses anchor-specific CF files:
- `<cc>_cf_2023_bestsite_p95_anchor-wind_onshore.csv`
- `<cc>_cf_2023_bestsite_p95_anchor-solar.csv`

Important:
- Land-only mixes → exact same grid cell
- Offshore technologies are mapped via deterministic matching logic

## Interpretation

RES mix scenarios represent **co-located hybrid systems**, not spatial portfolios.
---

# Future extensions


- Hybrid RES + grid supply  
- Electrolyser minimum load constraint  
- Seasonal storage technologies  
- Degradation modelling 