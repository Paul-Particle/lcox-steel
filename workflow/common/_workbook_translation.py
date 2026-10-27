"""The input workbook's keys and scenario vocabulary, translated into the model's.

Everything that depends on how data/assumptions/Assumptions_yaml_inputs.xlsx
names things lives here, so generate_inputs.py only ever sees the model's names.
"""

import pandas as pd

RES_TECHS = ("wind-onshore", "solar", "wind-offshore")
IRON_OPEX = {"default": {"mean": ["iron.fom_usd_per_t_capacity_per_year_aus",
                                  "iron.fom_usd_per_t_capacity_per_year_bra",
                                  "iron.fom_usd_per_t_capacity_per_year_eu"]},
             "aus": "iron.fom_usd_per_t_capacity_per_year_aus",
             "bra": "iron.fom_usd_per_t_capacity_per_year_bra",
             "eu": "iron.fom_usd_per_t_capacity_per_year_eu"}
CARRIERS = ("biomass", "brown_coal", "coal_gas", "energy_storage", "gas", "geothermal", "hard_coal",
            "hydro", "hydro_reservoir", "hydro_river", "marine", "nuclear", "oil", "oil_shale", "other",
            "other_re", "peat", "pumped_storage", "solar", "waste", "wind_offshore", "wind_onshore")

# Model key -> the workbook row(s) it is taken from: one key; {"sum" | "mean": [keys]};
# or named alternatives, the first being the default.
WORKBOOK_SOURCES = {
    "finance.origin_wacc": "finance.other.eu_wacc",
    "finance.destination_wacc": "finance.other.eu_wacc",
    "plant.h2_intensity_kg_per_t_dri": "dri_h2.kg_h2_per_t_dri",
    **{path: source for tech in RES_TECHS for path, source in {
        f"res.{tech}.capex_per_mw_eur": {v: f"res.{tech}.capex_usd_per_mw_{v}" for v in ("base", "low", "high")},
        f"res.{tech}.opex_per_mw_per_year_eur": {v: f"res.{tech}.opex_usd_per_mw_per_year_{v}"
                                                 for v in ("base", "low", "high")},
        f"res.{tech}.lifetime_years": "res.financing_period",
    }.items()},
    "electrolyser.capex_per_mw_eur": "electrolyser.capex_usd_per_mw_base",
    "electrolyser.opex_per_mw_per_year_eur": "electrolyser.opex_usd_per_mw_per_year_base",
    "electrolyser.lifetime_years": "electrolyser.financing_period",
    "electrolyser.efficiency_kwh_per_kg": "electrolyser.efficiency_kwh_per_kg",
    "battery.capex_per_mw_eur": {"global": "battery.capex_usd_per_mw_global", "europe": "battery.capex_usd_per_mw_eu"},
    "battery.capex_per_mwh_eur": {"global": "battery.capex_usd_per_mwh_global", "europe": "battery.capex_usd_per_mwh_eu"},
    "battery.lifetime_years": "battery.financing_period",
    "battery.efficiency_roundtrip": "battery.efficiency_roundtrip",
    "h2_buffer.capex_per_mwh_eur": "h2_buffer.capex_usd_per_mwh",
    "h2_buffer.lifetime_years": "h2_buffer.financing_period",
    **{path: source for shaft in ("dri-h2", "dri-mix", "dri-ng") for path, source in {
        f"{shaft}.capex_per_t_per_year_eur": "dri.capex_usd_per_t",
        f"{shaft}.opex_per_t_per_year_eur": IRON_OPEX,
        f"{shaft}.lifetime_years": "IS.financing_period",
        f"{shaft}.ore_eur_per_t": "dri.ore_usd_per_t",
        f"{shaft}.iron_t_per_t_steel": "eaf_dri.iron_t_per_t_steel",
    }.items()},
    "dri-h2.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "dri-h2.h2_preheat_el_mwh_per_t": "dri_h2.el_mwh_per_t",
    "dri-mix.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "dri-mix.h2_preheat_el_mwh_per_t": "dri_h2.el_mwh_per_t",
    "dri-ng.gas_mwh_per_t": "dri_ng.gas_mwh_per_t",
    "dri-ng.el_mwh_per_t": "dri_ng.el_mwh_per_t",
    "eaf.capex_per_t_per_year_eur": "eaf.capex_usd_per_t",
    "eaf.opex_per_t_per_year_eur": {"default": {"mean": ["eaf.fom_usd_per_t_capacity_per_year_aus",
                                                         "eaf.fom_usd_per_t_capacity_per_year_bra",
                                                         "eaf.fom_usd_per_t_capacity_per_year_eu"]},
                                    "aus": "eaf.fom_usd_per_t_capacity_per_year_aus",
                                    "bra": "eaf.fom_usd_per_t_capacity_per_year_bra",
                                    "eu": "eaf.fom_usd_per_t_capacity_per_year_eu"},
    "eaf.lifetime_years": "IS.financing_period",
    "eaf.charge.briquettes.el_mwh_per_t": "eaf_hbi.el_mwh_per_t",
    "eaf.charge.plates.el_mwh_per_t": "eaf_ewin.el_mwh_per_t",
    "eaf.charge.hot.el_mwh_per_t": "eaf_dri.el_mwh_per_t",
    "eaf.charge.liquid.el_mwh_per_t": "eaf_moe.el_mwh_per_t",
    "moe.capex_per_t_per_year_eur": "moe.capex_usd_per_t",
    "moe.opex_per_t_per_year_eur": IRON_OPEX,
    "moe.lifetime_years": "IS.financing_period",
    "moe.el_mwh_per_t": "moe.el_mwh_per_t",
    "moe.ore_eur_per_t": "moe.ore_usd_per_t_liquid_iron",
    "ew.capex_per_t_per_year_eur": "ewin.capex_usd_per_t",
    "ew.opex_per_t_per_year_eur": IRON_OPEX,
    "ew.lifetime_years": "IS.financing_period",
    "ew.el_mwh_per_t": "ewin.el_mwh_per_t",
    "ew.ore_eur_per_t": "ewin.ore_usd_per_t_electrolytic_iron",
    "natural_gas.price_eur_per_mwh": {"default": {"mean": ["natural_gas.price_usd_per_mwh_high",
                                                           "natural_gas.price_usd_per_mwh_low"]},
                                      "high": "natural_gas.price_usd_per_mwh_high",
                                      "low": "natural_gas.price_usd_per_mwh_low"},
    "briquetting.el_mwh_per_t": "briquetting.el_mwh_per_t_hbi",
    "iron_store.capex_per_t_eur": "iron_steel_store.capex_usd_per_t",
    "iron_store.lifetime_years": "iron_steel_store.financing_period",
    "transport.sea.iron.eur_per_t": "shipping.cost_usd_per_t",
    "transport.sea.iron.eur_per_t_km": "shipping.cost_usd_per_t_per_km",
    "transport.sea.steel.eur_per_t": "shipping.cost_usd_per_t",
    "transport.sea.steel.eur_per_t_km": "shipping.cost_usd_per_t_per_km",
    **{f"transport.distance_km.{origin}.{destination}.sea":
       f"shipping.distance_km_{origin.lower()}_{destination.lower()}"
       for origin in ("AUS", "BRA") for destination in ("ESP", "DEU", "FRA")},
    **{f"emissions.default_grid_emission_factors.electricity_t_co2e_per_mwh.{carrier}":
       f"emissions.electricity_t_co2e_per_mwh.{carrier}.combustion" for carrier in CARRIERS},
}

# Values the generator picks per scenario rather than writing into the base:
# `wacc.{region}` by the region of an origin or destination.
LOOKUP_SOURCES = {
    "wacc.eu": "finance.other.eu_wacc",
    "wacc.noneu": "finance.other.noneu_wacc",
}
# Which WACC row of the workbook finances a country.
WACC_REGIONS = {"DEU": "eu", "ESP": "eu", "FRA": "eu", "AUS": "noneu", "BRA": "noneu"}

# Workbook override keys that stand for several model keys.
OVERRIDE_KEYS = {"finance.default_wacc": ("finance.origin_wacc", "finance.destination_wacc")}

# The scenario tab's vocabulary, where it differs from the model's ids.
AREAS = {"GER": "DEU"}
ROUTES = {"ewin-eaf": "ew-eaf"}
SITE_VARIANTS = {"national-average": "area-average"}
POWER_SUPPLIES = {"islanded_res": "islanded", "res_grid": "grid"}
# A route whose steel is made somewhere else than its iron is the model's
# `-export` twin of the route the workbook names.
EXPORT_SUFFIX = "-export"


def translate_scenario_rows(scenario_rows: pd.DataFrame) -> pd.DataFrame:
    """The scenario tab's rows with areas, routes, site variants and power supply in the model's ids."""
    translated_rows = scenario_rows.assign(
        iron_location=lambda df: df["iron_location"].replace(AREAS),
        steel_location=lambda df: df["steel_location"].replace(AREAS),
        route=lambda df: df["route"].replace(ROUTES),
        site_variant=lambda df: df["site_variant"].replace(SITE_VARIANTS),
        power_supply=lambda df: df["power_supply"].replace(POWER_SUPPLIES),
    )
    return translated_rows


def translate_override(key: str, value) -> list[tuple[str, object]]:
    """An override written in the workbook's keys, as (model key, value or alternative name) pairs.

    The workbook names a key by its stem (`electrolyser.capex_usd_per_mw`) and an
    alternative by its full key (`electrolyser.capex_usd_per_mw_low`). One workbook
    row can feed several model keys, and the override then reaches all of them. A
    model key, a number or a `±n%` change passes through unchanged.
    """
    alternatives_by_workbook_key = {}
    for path, spec in WORKBOOK_SOURCES.items():
        for name, source in ({"default": spec} if isinstance(spec, str) else spec).items():
            if isinstance(source, str):
                alternatives_by_workbook_key.setdefault(source, []).append((path, name))
    matches = alternatives_by_workbook_key.get(key) or alternatives_by_workbook_key.get(f"{key}_base") or []
    paths = OVERRIDE_KEYS.get(key) or list(dict.fromkeys(path for path, _ in matches)) or [key]
    translated_value = alternatives_by_workbook_key[value][0][1] if value in alternatives_by_workbook_key else value
    translated = [(path, translated_value) for path in paths]
    return translated
