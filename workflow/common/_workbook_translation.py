"""The input workbook's field names, translated into the model's keys.

What the scenario tab's terms mean for a scenario is the generator's business
(generate_inputs.py); this file only says which row feeds which key.
"""

RES_TECHS = ("wind-onshore", "solar", "wind-offshore")
CARRIERS = ("biomass", "brown_coal", "coal_gas", "energy_storage", "gas", "geothermal", "hard_coal",
            "hydro", "hydro_reservoir", "hydro_river", "marine", "nuclear", "oil", "oil_shale", "other",
            "other_re", "peat", "pumped_storage", "solar", "waste", "wind_offshore", "wind_onshore")

# Model key -> the workbook row it is taken from.
RENAMED_SOURCES = {
    "plant.h2_intensity_kg_per_t_dri": "dri_h2.kg_h2_per_t_dri",
    "res.wind-onshore.capex_per_mw_eur": "res.wind-onshore.capex_eur_per_mw_base",
    "res.wind-onshore.opex_per_mw_per_year_eur": "res.wind-onshore.opex_eur_per_mw_per_year_base",
    "res.solar.capex_per_mw_eur": "res.solar.capex_eur_per_mw_base",
    "res.solar.opex_per_mw_per_year_eur": "res.solar.opex_eur_per_mw_per_year_base",
    "res.wind-offshore.capex_per_mw_eur": "res.wind-offshore.capex_eur_per_mw_base",
    "res.wind-offshore.opex_per_mw_per_year_eur": "res.wind-offshore.opex_eur_per_mw_per_year_base",
    "electrolyser.capex_per_mw_eur": "electrolyser.capex_eur_per_mw_base",
    "electrolyser.opex_per_mw_per_year_eur": "electrolyser.opex_eur_per_mw_per_year_base",
    "electrolyser.efficiency_kwh_per_kg": "electrolyser.efficiency_kwh_per_kg",
    "battery.efficiency_roundtrip": "battery.efficiency_roundtrip",
    "h2_buffer.capex_per_mwh_eur": "h2_buffer.capex_eur_per_mwh",
    "dri-ng.gas_mwh_per_t": "dri_ng.gas_mwh_per_t",
    "dri-ng.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "eaf.capex_per_t_per_year_eur": "eaf.capex_eur_per_t",
    "moe.capex_per_t_per_year_eur": "moe.capex_eur_per_t",
    "moe.el_mwh_per_t": "moe.el_mwh_per_t",
    "moe.ore_eur_per_t": "moe.ore_eur_per_t_liquid_iron",
    "ew.capex_per_t_per_year_eur": "ewin.capex_eur_per_t",
    "ew.el_mwh_per_t": "ewin.el_mwh_per_t",
    "ew.ore_eur_per_t": "ewin.ore_eur_per_t_electrolytic_iron",
    "briquetting.el_mwh_per_t": "briquetting.el_mwh_per_t_hbi",
    **{f"res.{tech}.lifetime_years": "res.financing_period" for tech in RES_TECHS},
    "electrolyser.lifetime_years": "electrolyser.financing_period",
    "battery.lifetime_years": "battery.financing_period",
    "h2_buffer.lifetime_years": "h2_buffer.financing_period",
    **{f"{plant}.lifetime_years": "IS.financing_period" for plant in ("dri-h2", "dri-mix", "dri-ng", "eaf", "moe", "ew")},
    **{f"emissions.default_grid_emission_factors.electricity_t_co2e_per_mwh.{carrier}":
       f"emissions.electricity_t_co2e_per_mwh.{carrier}.combustion" for carrier in CARRIERS},
}
CONVERTED_SOURCES = {
    **{f"transport.distance_km.{origin}.{destination}.sea":
       f"shipping.distance_km_{origin.lower()}_{destination.lower()}"
       for origin in ("AUS", "BRA") for destination in ("ESP", "DEU", "FRA")},
}
SHARED_SOURCES = {
    "dri-h2.capex_per_t_per_year_eur": "dri.capex_eur_per_t",
    "dri-mix.capex_per_t_per_year_eur": "dri.capex_eur_per_t",
    "dri-ng.capex_per_t_per_year_eur": "dri.capex_eur_per_t",
    "dri-h2.ore_eur_per_t": "dri.ore_eur_per_t",
    "dri-mix.ore_eur_per_t": "dri.ore_eur_per_t",
    "dri-ng.ore_eur_per_t": "dri.ore_eur_per_t",
    "dri-h2.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    "dri-mix.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    "dri-ng.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    "dri-h2.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "dri-mix.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "iron_store.capex_per_t_eur": "iron_steel_store.capex_eur_per_t",
    "steel_store.capex_per_t_eur": "iron_steel_store.capex_eur_per_t",
    "iron_store.lifetime_years": "iron_steel_store.financing_period",
    "steel_store.lifetime_years": "iron_steel_store.financing_period",
    "transport.sea.iron.eur_per_t": "shipping.cost_eur_per_t",
    "transport.sea.steel.eur_per_t": "shipping.cost_eur_per_t",
    "transport.sea.iron.eur_per_t_km": "shipping.cost_eur_per_t_per_km",
    "transport.sea.steel.eur_per_t_km": "shipping.cost_eur_per_t_per_km",
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
# Keys the model needs and the workbook gives no value for; the template's value
# stands in, and the generator warns until the workbook gives one. Where the workbook
# has a row, it holds 0 in both value and unit and is left unread.
#   "dri-h2.p_min_pu": "dri.p_min_pu",
#   "dri-mix.p_min_pu": "dri.p_min_pu",
#   "dri-ng.p_min_pu": "dri.p_min_pu",
#   "moe.p_min_pu": "moe.p_min_pu",
#   "steel_store.max_weeks": "steel_store.max_weeks",
#   "grid.connection_capex_eur_per_mw": "grid.connection_capex_eur_per_mw",
#   "grid.connection_lifetime_years": "grid.connection_lifetime_years",
#   "grid.fee_eur_per_mw_per_year": "grid.fee_eur_per_mw_per_year",
#   "grid.fee_eur_per_mwh": "grid.fee_eur_per_mwh",
#   "transmission.cost_per_mw_per_km_eur": "transmission.cost_per_mw_per_km_eur",
#   "transmission.lifetime_years": "transmission.lifetime_years",
#   "transmission.losses_pct_per_1000km": "transmission.losses_pct_per_1000km",
#   "transmission.indirect_route_factor": "transmission.indirect_route_factor",
#   "anchor.min_capacity_mw"  (the anchor floor: https://github.com/Paul-Particle/lcox-steel/issues/47)
#   "briquetting.capex_per_t_per_year_eur"
#   "briquetting.opex_per_t_per_year_eur"
#   "briquetting.lifetime_years"
#   "briquetting.yield_t_per_t"
#   "eaf.consumables_eur_per_t"
#   "plant.steel_mt_per_year"
#   "natural_gas.price_eur_per_mwh"  (the workbook gives only its low and high)
#   "electrolyser.varopex_eur_per_mwh_el"  (the workbook gives it per MWh of H2: a model change)
# Rows held at 0 for keys the model does not have (a minimum load: to be discussed):
#   "eaf.p_min_pu", "ewin.p_min_pu"

# Workbook override keys that stand for several model keys.
OVERRIDE_KEYS = {"finance.default_wacc": ("finance.origin_wacc", "finance.destination_wacc")}


def translate_override(key: str, value, *source_tables: dict) -> list[tuple[str, object]]:
    """An override written in the workbook's keys, as (model key, value or row name) pairs.

    The workbook names a key by its stem (`res.solar.capex_eur_per_mw`) and one of its
    named rows (a sensitivity level or a region) by its full key
    (`res.solar.capex_eur_per_mw_low`), which becomes that row's name (`low`). One
    workbook row can feed several model keys, and the override then reaches all of
    them. A model key, a number or a `±n%` change passes through unchanged.
    """
    paths_and_names_by_row = {}
    for source_table in source_tables:
        for path, spec in source_table.items():
            named_rows = spec.items() if isinstance(spec, dict) else [(None, spec)]
            for name, row in named_rows:
                if isinstance(row, str):
                    paths_and_names_by_row.setdefault(row, []).append((path, name))
    matches = paths_and_names_by_row.get(key) or paths_and_names_by_row.get(f"{key}_base") or []
    paths = OVERRIDE_KEYS.get(key) or list(dict.fromkeys(path for path, _ in matches)) or [key]
    row_name = next((name for _, name in paths_and_names_by_row.get(value, []) if name), None)
    translated = [(path, row_name or value) for path in paths]
    return translated
