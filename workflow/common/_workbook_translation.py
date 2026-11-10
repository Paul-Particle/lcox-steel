"""The input workbook's field names, translated into the model's keys.

A workbook row named as a model key is read as that key and is not listed here;
this file says which row feeds which key where the names differ. What the scenario
tab's terms mean for a scenario, and which rows give a key its sensitivity levels
and regional values, is the generator's business (generate_inputs.py).
"""

# Model key -> the workbook row it is taken from, where the row names the same thing
# in other words.
RENAMED_SOURCES = {
    "h2-dri.kg_h2_per_t_dri": "dri_h2.kg_h2_per_t_dri",
    "ng-dri.gas_mwh_per_t": "dri_ng.gas_mwh_per_t",
    "ng-dri.el_mwh_per_t": "dri_ng.el_mwh_per_t",
}
# Model key -> the workbook row it is taken from, where the model names it for what it
# is in the model.
MODEL_NAMED_SOURCES = {
    "eaf.charge.briquettes.el_mwh_per_t": "eaf_hbi.el_mwh_per_t",
    "eaf.charge.plates.el_mwh_per_t": "eaf_ewin.el_mwh_per_t",
    "eaf.charge.hot.el_mwh_per_t": "eaf_dri.el_mwh_per_t",
    "eaf.charge.liquid.el_mwh_per_t": "eaf_moe.el_mwh_per_t",
    "ewin.iron_t_per_t_steel": "eaf_ewin.iron_t_per_t_steel",
    "moe.iron_t_per_t_steel": "eaf_moe.iron_t_per_t_steel",
    "h2-dri.h2_preheat_el_mwh_per_t": "dri_h2.el_mwh_per_t",
}
# Model key -> the workbook row it is converted from.
CONVERTED_SOURCES = {
    **{f"transport.distance_km.{origin}.{destination}.sea":
       f"shipping.distance_km_{origin.lower()}_{destination.lower()}"
       for origin in ("AUS", "BRA") for destination in ("ESP", "DEU", "FRA")},
}
# Model key -> the workbook row it is taken from, where one row feeds several keys. A
# row renamed or model-named above feeds its own key there and the other keys here.
SHARED_SOURCES = {
    # Renamed above, shared here.
    "mix-dri.kg_h2_per_t_dri": "dri_h2.kg_h2_per_t_dri",
    "h2-dri.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "mix-dri.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    # Model-named above, shared here.
    "mix-dri.h2_preheat_el_mwh_per_t": "dri_h2.el_mwh_per_t",
    # No one key of the model is the row's own.
    **{f"res.{tech}.financing_period": "res.financing_period" for tech in ("wind-onshore", "solar", "wind-offshore")},
    **{f"{plant}.financing_period": "IS.financing_period" for plant in ("h2-dri", "mix-dri", "ng-dri", "eaf", "moe", "ewin")},
    **{f"{shaft}.capex_eur_per_t": "dri.capex_eur_per_t" for shaft in ("h2-dri", "mix-dri", "ng-dri")},
    **{f"{shaft}.ore_eur_per_t_dri": "dri.ore_eur_per_t" for shaft in ("h2-dri", "mix-dri", "ng-dri")},
    **{f"{shaft}.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel" for shaft in ("h2-dri", "mix-dri", "ng-dri")},
    **{f"{store}.{key}": f"iron_steel_store.{key}"
       for store in ("iron_store", "steel_store") for key in ("capex_eur_per_t", "financing_period")},
    **{f"transport.sea.{cargo}.{key}": f"shipping.{key}"
       for cargo in ("iron", "steel") for key in ("cost_eur_per_t", "cost_eur_per_t_per_km")},
}
WORKBOOK_SOURCES = RENAMED_SOURCES | MODEL_NAMED_SOURCES | CONVERTED_SOURCES | SHARED_SOURCES
# Held at 0 in both value and unit: not provided yet, so the template's value stands in.
# Each row already carries the model key it will be read under once given.
ROWS_NOT_PROVIDED_YET = (
    "steel_store.max_weeks",
    "grid.connection_capex_eur_per_mw", "grid.connection_financing_period",
    "grid.fee_eur_per_mw_per_year", "grid.fee_eur_per_mwh",
    "transmission.capex_eur_per_mw_per_km", "transmission.financing_period",
    "transmission.losses_pct_per_1000km", "transmission.indirect_route_factor",
)
# Workbook row -> why the generator does not read it, a row named as a model key
# included. A row the generator neither reads nor finds here stops it.
UNREAD_ROWS = {
    **{row: "not provided yet" for row in ROWS_NOT_PROVIDED_YET},
    # Held at 0 in both value and unit; the template's minimum loads stand in where the
    # model has one (h2-dri, mix-dri, ng-dri, moe).
    **{row: "to be discussed: is a minimum load the right approach?"
       for row in ("dri.p_min_pu", "moe.p_min_pu", "eaf.p_min_pu", "ewin.p_min_pu")},
    "finance.res.eu_wacc": "to be discussed: a WACC of its own for the renewables",
    "finance.res.noneu_wacc": "to be discussed: a WACC of its own for the renewables",
    "dri_h2.ng_mwh_per_t": "to be discussed: natural gas into the H2 shaft?",
    "eaf_hbi.iron_t_per_t_steel": "to be discussed: the iron ratio of the HBI charge, which repeats the DRI's",
    "battery.max_hours": "the model does not use it: the optimiser sizes the battery's power and energy apart",
}
# Model key -> why the template's value stands in for the workbook's scenarios. The
# generator warns with these reasons, and stops on a key it takes from the template that
# is not listed here, or one listed that the workbook or the scenarios' overlays now set.
TEMPLATE_BACKSTOPPED_KEYS = {
    "anchor.min_capacity_mw": "no workbook row: the anchor floor, https://github.com/Paul-Particle/lcox-steel/issues/47",
    "battery.fom_eur_per_mw_per_year": "no workbook row; the template's is eyeballed",
    "battery.fom_eur_per_mwh_per_year": "no workbook row; the template's is eyeballed",
    "briquetting.capex_eur_per_t": "no workbook row",
    "briquetting.fom_eur_per_t_capacity_per_year": "no workbook row",
    "briquetting.financing_period": "no workbook row",
    "briquetting.yield_t_per_t": "no workbook row",
    "eaf.consumables_eur_per_t": "no workbook row",
    "plant.steel_mt_per_year": "no workbook row",
    "natural_gas.price_eur_per_mwh": "the workbook gives only its low and high",
    **{key: "its workbook row is not provided yet" for key in ROWS_NOT_PROVIDED_YET},
    **{f"{plant}.p_min_pu": "to be discussed: is a minimum load the right approach?"
       for plant in ("h2-dri", "mix-dri", "ng-dri", "moe")},
}

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
    "mix-dri.reductant_store_hours": float("nan"),
    "emissions.co2_price_eur_per_t": 0.0,
}
