"""
09b_battery_cost_sweep.py

Diagnostic script to find the battery cost threshold where offshore wind
switches on/off in the PyPSA optimisation.
"""

import pypsa
import pandas as pd
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent.parent
CF_DIR       = PROJECT_ROOT / "data" / "res_cf"

# ── Load CF data ───────────────────────────────────────────────────────────────
cf_file = CF_DIR / "de_cf_2023_bestsite_p95_anchor-wind_onshore.csv"
df = pd.read_csv(cf_file, parse_dates=["time"])
df = df.set_index("time")

wind_onshore_cf  = df["wind_onshore_cf"]
wind_offshore_cf = df["wind_offshore_cf"]
solar_cf         = df["solar_cf"]

# ── Fixed parameters ───────────────────────────────────────────────────────────
electrolyser_mw   = 628.0
battery_max_hours = 6
electrolyser_load = pd.Series(electrolyser_mw, index=df.index)

cost_wind_onshore  = 80_000
cost_wind_offshore = 130_000
cost_solar         = 40_000

# ── Battery cost sweep ─────────────────────────────────────────────────────────
battery_costs_to_test = [5_000, 8_000, 10_000, 12_000, 15_000,
                          20_000, 25_000, 30_000, 40_000, 50_000]

print("\n--- Battery Cost Sweep ---")
print(f"{'Battery €/MWh':>15} {'Offshore MW':>12} {'Onshore MW':>12} {'Solar MW':>12} {'Battery MW':>12} {'Total Cost €M':>15}")

for bat_cost in battery_costs_to_test:
    n = pypsa.Network()
    n.set_snapshots(df.index)
    n.add("Bus", "electricity_bus")
    n.add("Generator", "wind_onshore",  bus="electricity_bus",
          p_nom_extendable=True, p_max_pu=wind_onshore_cf,
          carrier="wind_onshore",  capital_cost=cost_wind_onshore)
    n.add("Generator", "wind_offshore", bus="electricity_bus",
          p_nom_extendable=True, p_max_pu=wind_offshore_cf,
          carrier="wind_offshore", capital_cost=cost_wind_offshore)
    n.add("Generator", "solar",         bus="electricity_bus",
          p_nom_extendable=True, p_max_pu=solar_cf,
          carrier="solar",        capital_cost=cost_solar)
    n.add("StorageUnit", "battery",     bus="electricity_bus",
          p_nom_extendable=True, max_hours=battery_max_hours,
          efficiency_store=0.95, efficiency_dispatch=0.95,
          cyclic_state_of_charge=True,
          capital_cost=bat_cost * battery_max_hours)
    n.add("Load", "electrolyser", bus="electricity_bus",
          p_set=electrolyser_load)

    n.optimize(solver_name="highs", solver_options={"output_flag": False})

    offshore = n.generators.loc["wind_offshore", "p_nom_opt"]
    onshore  = n.generators.loc["wind_onshore",  "p_nom_opt"]
    solar    = n.generators.loc["solar",          "p_nom_opt"]
    battery  = n.storage_units.loc["battery",     "p_nom_opt"]
    cost     = n.objective / 1e6

    print(f"{bat_cost:>15,} {offshore:>12.1f} {onshore:>12.1f} {solar:>12.1f} {battery:>12.1f} {cost:>15.1f}")