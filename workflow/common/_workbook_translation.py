"""The input workbook's field names, translated into the model's keys.

What the scenario tab's terms mean for a scenario is the generator's business
(generate_inputs.py); this file only says which row feeds which key.
"""

RES_TECHS = ("wind-onshore", "solar", "wind-offshore")
CARRIERS = ("biomass", "brown_coal", "coal_gas", "energy_storage", "gas", "geothermal", "hard_coal",
            "hydro", "hydro_reservoir", "hydro_river", "marine", "nuclear", "oil", "oil_shale", "other",
            "other_re", "peat", "pumped_storage", "solar", "waste", "wind_offshore", "wind_onshore")

# Model key -> the workbook row(s) it is taken from: one key; {"sum" | "mean": [keys]};
# or named alternatives, the first being the default.
RENAMED_SOURCES = {
    "plant.h2_intensity_kg_per_t_dri": "dri_h2.kg_h2_per_t_dri",
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
PER_SCENARIO_SOURCES = {
    **{path: source for tech in RES_TECHS for path, source in {
        f"res.{tech}.capex_per_mw_eur": {v: f"res.{tech}.capex_eur_per_mw_{v}" for v in ("base", "low", "high")},
        f"res.{tech}.opex_per_mw_per_year_eur": {v: f"res.{tech}.opex_eur_per_mw_per_year_{v}"
                                                 for v in ("base", "low", "high")},
    }.items()},
    "natural_gas.price_eur_per_mwh": {"default": {"mean": ["natural_gas.price_eur_per_mwh_high",
                                                           "natural_gas.price_eur_per_mwh_low"]},
                                      "high": "natural_gas.price_eur_per_mwh_high",
                                      "low": "natural_gas.price_eur_per_mwh_low"},
}
WORKBOOK_SOURCES = (RENAMED_SOURCES | CONVERTED_SOURCES | SHARED_SOURCES | MODEL_NAMED_SOURCES
                    | PER_SCENARIO_SOURCES)
# Rows the workbook still holds at 0, value and unit:
#   "dri-h2.p_min_pu": "dri.p_min_pu",
#   "dri-mix.p_min_pu": "dri.p_min_pu",
#   "dri-ng.p_min_pu": "dri.p_min_pu",
#   "eaf.p_min_pu": "eaf.p_min_pu",
#   "moe.p_min_pu": "moe.p_min_pu",
#   "ew.p_min_pu": "ewin.p_min_pu",
#   "steel_store.max_weeks": "steel_store.max_weeks",
#   "grid.connection_capex_eur_per_mw": "grid.connection_capex_eur_per_mw",
#   "grid.connection_lifetime_years": "grid.connection_lifetime_years",
#   "grid.fee_eur_per_mw_per_year": "grid.fee_eur_per_mw_per_year",
#   "grid.fee_eur_per_mwh": "grid.fee_eur_per_mwh",
#   "transmission.cost_per_mw_per_km_eur": "transmission.cost_per_mw_per_km_eur",
#   "transmission.lifetime_years": "transmission.lifetime_years",
#   "transmission.losses_pct_per_1000km": "transmission.losses_pct_per_1000km",
#   "transmission.indirect_route_factor": "transmission.indirect_route_factor",

# Workbook override keys that stand for several model keys.
OVERRIDE_KEYS = {"finance.default_wacc": ("finance.origin_wacc", "finance.destination_wacc")}


def translate_override(key: str, value, sources: dict) -> list[tuple[str, object]]:
    """An override written in the workbook's keys, as (model key, value or alternative name) pairs.

    The workbook names a key by its stem (`electrolyser.capex_eur_per_mw`) and an
    alternative by its full key (`electrolyser.capex_eur_per_mw_low`). One workbook
    row can feed several model keys, and the override then reaches all of them. A
    model key, a number or a `±n%` change passes through unchanged.
    """
    alternatives_by_workbook_key = {}
    for path, spec in sources.items():
        for name, source in ({"default": spec} if isinstance(spec, str) else spec).items():
            if isinstance(source, str):
                alternatives_by_workbook_key.setdefault(source, []).append((path, name))
    matches = alternatives_by_workbook_key.get(key) or alternatives_by_workbook_key.get(f"{key}_base") or []
    paths = OVERRIDE_KEYS.get(key) or list(dict.fromkeys(path for path, _ in matches)) or [key]
    translated_value = alternatives_by_workbook_key[value][0][1] if value in alternatives_by_workbook_key else value
    translated = [(path, translated_value) for path in paths]
    return translated
