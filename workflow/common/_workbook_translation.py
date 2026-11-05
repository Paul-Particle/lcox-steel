"""The input workbook's field names, translated into the model's keys.

What the scenario tab's terms mean for a scenario is the generator's business
(generate_inputs.py); this file only says which row feeds which key.
"""

CARRIERS = ("biomass", "brown_coal", "coal_gas", "energy_storage", "gas", "geothermal", "hard_coal",
            "hydro", "hydro_reservoir", "hydro_river", "marine", "nuclear", "oil", "oil_shale", "other",
            "other_re", "peat", "pumped_storage", "solar", "waste", "wind_offshore", "wind_onshore")

# Model key -> the workbook row it is taken from.
RENAMED_SOURCES = {
    "plant.h2_intensity_kg_per_t_dri": "dri_h2.kg_h2_per_t_dri",
    "res.wind-onshore.capex_eur_per_mw": "res.wind-onshore.capex_eur_per_mw_base",
    "res.wind-onshore.fom_eur_per_mw_per_year": "res.wind-onshore.opex_eur_per_mw_per_year_base",
    "res.solar.capex_eur_per_mw": "res.solar.capex_eur_per_mw_base",
    "res.solar.fom_eur_per_mw_per_year": "res.solar.opex_eur_per_mw_per_year_base",
    "res.wind-offshore.capex_eur_per_mw": "res.wind-offshore.capex_eur_per_mw_base",
    "res.wind-offshore.fom_eur_per_mw_per_year": "res.wind-offshore.opex_eur_per_mw_per_year_base",
    "electrolyser.capex_eur_per_mw": "electrolyser.capex_eur_per_mw_base",
    "electrolyser.fom_eur_per_mw_per_year": "electrolyser.opex_eur_per_mw_per_year_base",
    "electrolyser.varopex_eur_per_mwh_h2": "electrolyser.varopex_eur_per_mwh_h2_base",
    "electrolyser.efficiency_kwh_per_kg": "electrolyser.efficiency_kwh_per_kg",
    "battery.efficiency_roundtrip": "battery.efficiency_roundtrip",
    "h2_buffer.capex_eur_per_mwh": "h2_buffer.capex_eur_per_mwh",
    "h2_buffer.fom_eur_per_mwh": "h2_buffer.fom_eur_per_mwh",
    "dri-ng.gas_mwh_per_t": "dri_ng.gas_mwh_per_t",
    "dri-ng.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "eaf.capex_eur_per_t": "eaf.capex_eur_per_t",
    "moe.capex_eur_per_t": "moe.capex_eur_per_t",
    "moe.el_mwh_per_t": "moe.el_mwh_per_t",
    "moe.ore_eur_per_t": "moe.ore_eur_per_t_liquid_iron",
    "ew.capex_eur_per_t": "ewin.capex_eur_per_t",
    "ew.el_mwh_per_t": "ewin.el_mwh_per_t",
    "ew.ore_eur_per_t": "ewin.ore_eur_per_t_electrolytic_iron",
    "briquetting.el_mwh_per_t": "briquetting.el_mwh_per_t_hbi",
    **{f"res.{tech}.financing_period": "res.financing_period" for tech in ("wind-onshore", "solar", "wind-offshore")},
    "electrolyser.financing_period": "electrolyser.financing_period",
    "battery.financing_period": "battery.financing_period",
    "h2_buffer.financing_period": "h2_buffer.financing_period",
    **{f"{plant}.financing_period": "IS.financing_period" for plant in ("dri-h2", "dri-mix", "dri-ng", "eaf", "moe", "ew")},
    **{f"emissions.electricity_t_co2e_per_mwh.{carrier}.combustion":
       f"emissions.electricity_t_co2e_per_mwh.{carrier}.combustion" for carrier in CARRIERS},
}
CONVERTED_SOURCES = {
    **{f"transport.distance_km.{origin}.{destination}.sea":
       f"shipping.distance_km_{origin.lower()}_{destination.lower()}"
       for origin in ("AUS", "BRA") for destination in ("ESP", "DEU", "FRA")},
}
SHARED_SOURCES = {
    "dri-h2.capex_eur_per_t": "dri.capex_eur_per_t",
    "dri-mix.capex_eur_per_t": "dri.capex_eur_per_t",
    "dri-ng.capex_eur_per_t": "dri.capex_eur_per_t",
    "dri-h2.ore_eur_per_t": "dri.ore_eur_per_t",
    "dri-mix.ore_eur_per_t": "dri.ore_eur_per_t",
    "dri-ng.ore_eur_per_t": "dri.ore_eur_per_t",
    "dri-h2.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    "dri-mix.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    "dri-ng.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    "dri-h2.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "dri-mix.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "iron_store.capex_eur_per_t": "iron_steel_store.capex_eur_per_t",
    "steel_store.capex_eur_per_t": "iron_steel_store.capex_eur_per_t",
    "iron_store.financing_period": "iron_steel_store.financing_period",
    "steel_store.financing_period": "iron_steel_store.financing_period",
    "transport.sea.iron.cost_eur_per_t": "shipping.cost_eur_per_t",
    "transport.sea.steel.cost_eur_per_t": "shipping.cost_eur_per_t",
    "transport.sea.iron.cost_eur_per_t_per_km": "shipping.cost_eur_per_t_per_km",
    "transport.sea.steel.cost_eur_per_t_per_km": "shipping.cost_eur_per_t_per_km",
}
MODEL_NAMED_SOURCES = {
    "eaf.charge.briquettes.el_mwh_per_t": "eaf_hbi.el_mwh_per_t",
    "eaf.charge.plates.el_mwh_per_t": "eaf_ewin.el_mwh_per_t",
    "eaf.charge.hot.el_mwh_per_t": "eaf_dri.el_mwh_per_t",
    "eaf.charge.liquid.el_mwh_per_t": "eaf_moe.el_mwh_per_t",
    "ew.iron_t_per_t_steel": "eaf_ewin.iron_t_per_t_steel",
    "moe.iron_t_per_t_steel": "eaf_moe.iron_t_per_t_steel",
    "dri-h2.h2_preheat_el_mwh_per_t": "dri_h2.el_mwh_per_t",
    "dri-mix.h2_preheat_el_mwh_per_t": "dri_h2.el_mwh_per_t",
}
WORKBOOK_SOURCES = RENAMED_SOURCES | CONVERTED_SOURCES | SHARED_SOURCES | MODEL_NAMED_SOURCES
# Workbook row -> why the generator does not read it. A row the generator neither reads
# nor finds here stops it.
UNREAD_ROWS = {
    # Held at 0 in both value and unit: not provided yet, so the template's value stands
    # in and the generator warns. Once given, each is read under the model key named.
    **{row: f"not provided: {row}" for row in (
        "steel_store.max_weeks",
        "grid.connection_capex_eur_per_mw", "grid.connection_lifetime_years",
        "grid.fee_eur_per_mw_per_year", "grid.fee_eur_per_mwh",
        "transmission.cost_per_mw_per_km_eur", "transmission.lifetime_years",
        "transmission.losses_pct_per_1000km", "transmission.indirect_route_factor")},
    # Held at 0 in both value and unit; the template's minimum loads stand in where the
    # model has one (dri-h2, dri-mix, dri-ng, moe).
    **{row: "to be discussed: is a minimum load the right approach?"
       for row in ("dri.p_min_pu", "moe.p_min_pu", "eaf.p_min_pu", "ewin.p_min_pu")},
    "finance.res.eu_wacc": "to be discussed: a WACC of its own for the renewables",
    "finance.res.noneu_wacc": "to be discussed: a WACC of its own for the renewables",
    "dri_h2.ng_mwh_per_t": "to be discussed: natural gas into the H2 shaft?",
    "eaf_hbi.iron_t_per_t_steel": "to be discussed: the iron ratio of the HBI charge, which repeats the DRI's",
    "battery.max_hours": "the model does not use it: the optimiser sizes the battery's power and energy apart",
}
# Keys the model needs and the workbook has no row for; the template's value stands in,
# and the generator warns until the workbook gives one:
#   "anchor.min_capacity_mw"  (the anchor floor: https://github.com/Paul-Particle/lcox-steel/issues/47)
#   "battery.fom_eur_per_mw_per_year"   (fixed O&M; the template's is eyeballed)
#   "battery.fom_eur_per_mwh_per_year"
#   "briquetting.capex_eur_per_t"
#   "briquetting.fom_eur_per_t_capacity_per_year"
#   "briquetting.financing_period"
#   "briquetting.yield_t_per_t"
#   "eaf.consumables_eur_per_t"
#   "plant.steel_mt_per_year"
#   "natural_gas.price_eur_per_mwh"  (the workbook gives only its low and high)

# Model key -> its value where no scenario the workbook configures requires the key. NaN
# stops any run that reads it; the carbon price is zero because no scenario prices carbon.
VALUES_OF_KEYS_NOT_REQUIRED_BY_THE_WORKBOOK_SCENARIO_CONFIGURATION = {
    **{f"transport.distance_km.{origin}.ESP.rail": float("nan") for origin in ("DEU", "FRA", "ESP", "AB", "ONT")},
    **{f"transport.distance_km.{origin}.ESP.sea": float("nan") for origin in ("AB", "ONT")},
    "transport.distance_km.DEMO.ES_SYN.sea": float("nan"),
    **{f"transport.rail.{cargo}.{price}": float("nan")
       for cargo in ("iron", "steel") for price in ("cost_eur_per_t", "cost_eur_per_t_per_km")},
    # Diagnostic only: no total counts freight emissions.
    "emissions.freight_kg_co2e_per_t_km.sea": float("nan"),
    "emissions.freight_kg_co2e_per_t_km.rail": float("nan"),
    "emissions.natural_gas_reductant_t_co2e_per_mwh.process": float("nan"),
    "emissions.natural_gas_reductant_t_co2e_per_mwh.upstream": float("nan"),
    "dri-mix.reductant_store_hours": float("nan"),
    "emissions.co2_price_eur_per_t": 0.0,
}
