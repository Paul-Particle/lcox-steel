"""
09_pypsa_optimise.py

Purpose
-------
Runs PyPSA optimisation for each bestsite anchor scenario.
Takes hourly CF time series from script 07 outputs and finds
the optimal RES + battery capacity mix to minimise system cost.

Inputs
------
data/res_cf/*_cf_*_bestsite_p95_anchor-*.csv  (excluding mix files)

Outputs
-------
results/pypsa/pypsa_results.csv  (one row per run, appended)
"""

import pypsa
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, UTC

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent.parent
CF_DIR       = PROJECT_ROOT / "data" / "res_cf"
RESULTS_DIR  = PROJECT_ROOT / "results" / "pypsa"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Fixed parameters ───────────────────────────────────────────────────────────
electrolyser_mw   = 628.0    # derived from DRI sizing in config (hardcoded for now)
battery_max_hours = 6        # max storage duration in hours — also used for cost conversion

# ── Placeholder costs (€/MW/year annualised) ──────────────────────────────────
# TODO: replace with real costs from config_hannah.yaml
cost_wind_onshore  = 80_000    # €/MW/year
cost_wind_offshore = 130_000   # €/MW/year
cost_solar         = 40_000    # €/MW/year
cost_battery_mwh   = 15_000    # €/MWh/year (energy capacity)

# ── Find bestsite anchor files ─────────────────────────────────────────────────
# Exclude mix files (files with "mix-" in the name)
bestsite_files = sorted([
    f for f in CF_DIR.glob("*_cf_*_bestsite_p95_anchor-*.csv")
    if "mix-" not in f.name
])

print(f"Found {len(bestsite_files)} bestsite anchor files to process\n")

# ── Loop over all anchor files ─────────────────────────────────────────────────
for cf_file in bestsite_files:

    print(f"Processing: {cf_file.name}")

    # ── Load CF data ───────────────────────────────────────────────────────────
    df = pd.read_csv(cf_file, parse_dates=["time"])
    df = df.set_index("time")

    wind_onshore_cf  = df["wind_onshore_cf"]
    wind_offshore_cf = df["wind_offshore_cf"]
    solar_cf         = df["solar_cf"]

    # ── Derive metadata from filename ──────────────────────────────────────────
    filename    = cf_file.stem
    country     = filename.split("_")[0].upper()
    scenario    = "bestsite"
    anchor_tech = filename.split("anchor-")[1]
    year        = int(filename.split("_cf_")[1][:4])

    # ── Electrolyser load ──────────────────────────────────────────────────────
    electrolyser_load = pd.Series(electrolyser_mw, index=df.index)

    # ── Build PyPSA network ────────────────────────────────────────────────────
    network = pypsa.Network()
    network.set_snapshots(df.index)

    # Single bus — connection point where all generation and demand meet
    network.add("Bus", "electricity_bus")

    # RES generators
    # p_max_pu = hourly CF — max fraction of installed capacity that can generate
    # p_nom_extendable = True — PyPSA decides how many MW to build
    network.add("Generator",
        "wind_onshore",
        bus="electricity_bus",
        p_nom_extendable=True,
        p_max_pu=wind_onshore_cf,
        carrier="wind_onshore",
        capital_cost=cost_wind_onshore,
    )
    network.add("Generator",
        "wind_offshore",
        bus="electricity_bus",
        p_nom_extendable=True,
        p_max_pu=wind_offshore_cf,
        carrier="wind_offshore",
        capital_cost=cost_wind_offshore,
    )
    network.add("Generator",
        "solar",
        bus="electricity_bus",
        p_nom_extendable=True,
        p_max_pu=solar_cf,
        carrier="solar",
        capital_cost=cost_solar,
    )

    # Battery storage
    # cyclic_state_of_charge = battery state at end of year = start of year
    # capital_cost = €/MWh × max_hours converts energy cost to per-MW basis
    network.add("StorageUnit",
        "battery",
        bus="electricity_bus",
        p_nom_extendable=True,
        max_hours=battery_max_hours,
        efficiency_store=0.95,
        efficiency_dispatch=0.95,
        cyclic_state_of_charge=True,
        capital_cost=cost_battery_mwh * battery_max_hours,
    )

    # Electrolyser — flat baseload demand every hour
    network.add("Load",
        "electrolyser",
        bus="electricity_bus",
        p_set=electrolyser_load,
    )

    # ── Solve ──────────────────────────────────────────────────────────────────
    network.optimize(
        solver_name="highs",
        solver_options={"output_flag": False},
    )

    # ── Extract results ────────────────────────────────────────────────────────
    wind_onshore_gen_mwh  = network.generators_t.p["wind_onshore"].sum()
    wind_offshore_gen_mwh = network.generators_t.p["wind_offshore"].sum()
    solar_gen_mwh         = network.generators_t.p["solar"].sum()

    # Curtailment = max possible generation minus actual generation
    curtailment      = (network.generators_t.p_max_pu * network.generators.p_nom_opt
                        - network.generators_t.p)
    annual_curtailment = curtailment.sum()
    total_possible     = (network.generators_t.p_max_pu * network.generators.p_nom_opt).sum()
    curtailment_rate   = annual_curtailment / total_possible * 100

    def safe_curtailment_pct(tech):
        # Returns 0 if tech not built, avoids division by zero
        cap = network.generators.loc[tech, "p_nom_opt"]
        if cap > 0:
            return float(curtailment_rate[tech])
        return 0.0

    battery_charge_mwh    = network.storage_units_t.p_store["battery"].sum()
    battery_discharge_mwh = network.storage_units_t.p_dispatch["battery"].sum()

    lcoe_eur_mwh = network.objective / electrolyser_load.sum()

    # ── Assemble results row ───────────────────────────────────────────────────
    results_row = {
        "country":                       country,
        "scenario":                      scenario,
        "anchor_tech":                   anchor_tech,
        "year":                          year,
        "wind_onshore_mw":               round(network.generators.loc["wind_onshore",  "p_nom_opt"], 1),
        "wind_offshore_mw":              round(network.generators.loc["wind_offshore", "p_nom_opt"], 1),
        "solar_mw":                      round(network.generators.loc["solar",         "p_nom_opt"], 1),
        "battery_mw":                    round(network.storage_units.loc["battery", "p_nom_opt"], 1),
        "battery_mwh":                   round(network.storage_units.loc["battery", "p_nom_opt"] * battery_max_hours, 1),
        "electrolyser_mw":               electrolyser_mw,
        "total_annual_cost_eur":         round(network.objective),
        "lcoe_eur_mwh":                  round(lcoe_eur_mwh, 2),
        "res_oversizing_ratio":          round((network.generators.loc["wind_onshore",  "p_nom_opt"] +
                                                network.generators.loc["wind_offshore", "p_nom_opt"] +
                                                network.generators.loc["solar",         "p_nom_opt"]) / electrolyser_mw, 2),  # total RES MW built per MW of electrolyser capacity
        "battery_oversizing_ratio":      round(network.storage_units.loc["battery", "p_nom_opt"] / electrolyser_mw, 2),       # battery MW built per MW of electrolyser capacity
        "res_to_demand_ratio":           round((wind_onshore_gen_mwh + wind_offshore_gen_mwh + solar_gen_mwh +
                                                annual_curtailment["wind_onshore"] +
                                                annual_curtailment["wind_offshore"] +
                                                annual_curtailment["solar"]) / electrolyser_load.sum(), 2),                    # max possible generation vs annual electrolyser demand
        "wind_onshore_gen_mwh":          round(wind_onshore_gen_mwh),
        "wind_offshore_gen_mwh":         round(wind_offshore_gen_mwh),
        "solar_gen_mwh":                 round(solar_gen_mwh),
        "wind_onshore_curtailment_mwh":  round(annual_curtailment["wind_onshore"]),
        "wind_offshore_curtailment_mwh": round(annual_curtailment["wind_offshore"]),
        "solar_curtailment_mwh":         round(annual_curtailment["solar"]),
        "wind_onshore_curtailment_pct":  round(safe_curtailment_pct("wind_onshore"), 1),
        "wind_offshore_curtailment_pct": round(safe_curtailment_pct("wind_offshore"), 1),
        "solar_curtailment_pct":         round(safe_curtailment_pct("solar"), 1),
        "battery_charge_mwh":            round(battery_charge_mwh),
        "battery_discharge_mwh":         round(battery_discharge_mwh),
        "run_timestamp_utc":             datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    # ── Save to CSV ────────────────────────────────────────────────────────────
    results_df = pd.DataFrame([results_row])
    out_path   = RESULTS_DIR / "pypsa_results.csv"

    if out_path.exists():
        existing = pd.read_csv(out_path)
        updated  = pd.concat([existing, results_df], ignore_index=True)
    else:
        updated = results_df

    updated.to_csv(out_path, index=False)
    print(f"  → saved: {country} | {scenario} | {anchor_tech} | LCOE: €{lcoe_eur_mwh:.1f}/MWh\n")

print("All files processed.")