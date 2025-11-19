import numpy as np
import pandas as pd
import sys

sys.path.insert(1, "../")

# config
# TODO: Move these into model_config.yml
operation_and_maintenance_estimate = 0.05  # fraction of initial capex, per year
WACC = 0.07  # post-tax, real terms
financing_period = 15  # years
capital_recovery_factor = WACC / (1 - ((1 + WACC) ** (-financing_period)))

electrolyzer_efficiency = 0.7  # estimate

# fixed value per MWh for grid cost. Based on highest consumption band.
# Should be actual calculation based on demand profile etc., because it has significant impact and average annual eurostat data is probably mostly based on large consumers with constant demand.
grid_cost = (
    29.7  # EUR/MWh Source: Eurostat. For Germany only (temp), will be a lookup later.
)

# emission_threshold = 28.2 # g/MJ H2 LHV; from delegated act on low carbon fuels
emission_threshold = 0.101  # t / MWh

# emission factors by tech (Source: Agorameter, except for nuclear, which is from the DA) t_CO2e/MWh
emissions_factors = pd.Series(
    {
        "biomass": 0,
        "energy_storage": 0,
        "brown_coal": 1.12,
        "coal_gas": 1.57,
        "gas": 0.38,
        "hard_coal": 0.82,
        "oil": 1.57,
        "oil_shale": 1.57,
        "peat": 1.57,
        "geothermal": 0,
        "pumped_storage": 0,
        "hydro_river": 0,
        "hydro_reservoir": 0,
        "marine": 0,
        "nuclear": 0.005454,
        "other": 1.57,
        "other_re": 0,
        "solar": 0,
        "waste": 1.57,
        "wind_offshore": 0,
        "wind_onshore": 0,
    }
)

CAPEX = pd.Series([750, 1500, 3000])  # EUR/kW
CAPEX.index = pd.Index(CAPEX.astype(str))
annutitized_capital_cost = CAPEX * capital_recovery_factor
annual_operations_and_maintenance = CAPEX * operation_and_maintenance_estimate
fixed_costs = (
    (annual_operations_and_maintenance + annutitized_capital_cost)
    / electrolyzer_efficiency
    * 1000
)  # EUR/MWh H2 LHV

# Get data and process
df = pd.read_feather("../data/processed_data.feather")
# some data from the 28th in 23 and 24 is missing for DE_LU, unclear why.
df = df.dropna()


# all of the following will be a loop over areas with one dfCountry for each area
def reorder_columns(dfCountry):
    macro_data = [
        "price",
        "demand_forecast",
        "wind_onshore_forecast",
        "wind_offshore_forecast",
        "wind_forecast",
        "solar_forecast",
        "vre_forecast",
        "residual",
        "demand",
    ]
    imports = [
        c for c in dfCountry.columns if c.startswith("from") and c != "from_sum"
    ] + ["from_sum"]
    exports = [c for c in dfCountry.columns if c.startswith("to") and c != "to_sum"] + [
        "to_sum"
    ]
    consumption = [c for c in dfCountry.columns if c.endswith("_cons")]
    time_data = ["dt", "hour", "doy", "doy_season"]

    new_columns = (
        time_data + macro_data + sorted(imports) + sorted(exports) + sorted(consumption)
    )
    production = set(dfCountry.columns) - set(new_columns)
    new_columns += sorted(production)

    return dfCountry.loc[:, new_columns]


dfC = reorder_columns(df.DE_LU)


# Calculate LCOH
def add_levelized_capex_cols(df_in, factors):
    levelized_capex_df = df_in.capacity_factor.apply(
        lambda x: factors / (len(df_in) * x)
    )
    return df_in.join(levelized_capex_df)


def add_total_cost_cols(df_in, cols):
    total_cost_df = (df_in.loc[:, cols].apply(lambda x: x + df_in.average_cost)).rename(
        columns=lambda c: f"{c}_total"
    )
    return df_in.join(total_cost_df)


def calculate_LCOH(dfCountry, year, fixed_costs):
    dfLCOH = (
        dfCountry.copy()
        .loc[year, ["price"]]
        .resample("1h")
        .first()
        .dropna()
        .sort_values("price")
        .assign(
            variable_cost=lambda df: (df.price + grid_cost) / electrolyzer_efficiency
        )  # grid cost also needs to account for inefficiency
        .assign(
            average_cost=lambda df: df.variable_cost.cumsum()
            / np.arange(1, len(df) + 1)
        )  # sorted price -> summing as we use more and more hours of the year + dividing by number of hours we use -> average for that capacity factor/hours of the year used
        .assign(capacity_factor=lambda df: np.arange(1, len(df) + 1) / len(df))
        .pipe(add_levelized_capex_cols, factors=fixed_costs)
        .pipe(add_total_cost_cols, cols=fixed_costs.index)
    )
    return dfLCOH


dfLCOH23 = calculate_LCOH(dfC, "2023", fixed_costs)
dfLCOH24 = calculate_LCOH(dfC, "2024", fixed_costs)


# results
def get_perfect_foresight_results(dfLCOH):
    m = dfLCOH.loc[:, [f"{c}_total" for c in fixed_costs.index]].idxmin()
    results = dfLCOH.loc[m]
    results.index = m.index
    optima = {}
    for r in results.iterrows():
        index = r[0]
        data = r[1]
        optimum = data
        optimum = optimum[
            [
                "price_max",
                "variable_cost_max",
                "average_cost",
                "capacity_factor",
                index.split("_")[0],
                index,
            ]
        ]
        optima[f"capital_cost_{index.split('_')[0]}"] = optimum.rename(
            {index.split("_")[0]: "fixed_cost", index: "total_cost"}
        )
    results = pd.DataFrame.from_dict(optima)
    return results


results23 = get_perfect_foresight_results(dfLCOH23)
results24 = get_perfect_foresight_results(dfLCOH24)


# heuristic results: Use previous year's optima times a factor to adjust for decreasing electricity costs as the limit over which no production takes place. Assumes electrolyzers are incentivized to minimize LCOH, does not reflect behavior expected for operation under a marginal hydrogen price.
def get_heuristic_results(
    previous_year_dfLCOH, dfLCOH, cost_decrease_expectation_factor=0.9
):
    m = previous_year_dfLCOH.loc[:, [f"{c}_total" for c in fixed_costs.index]].idxmin()
    previous_optima = previous_year_dfLCOH.loc[m]
    previous_optima.index = m.index
    heuristic_optima = {}
    for r in previous_optima.iterrows():
        index = r[0]
        data = r[1]
        m = dfLCOH.loc[:, "price"] <= data.price * cost_decrease_expectation_factor
        optimum = dfLCOH.loc[dfLCOH.loc[m].price.idxmax()]
        optimum = optimum[
            [
                "price_max",
                "variable_cost_max",
                "average_cost",
                "capacity_factor",
                index.split("_")[0],
                index,
            ]
        ]
        heuristic_optima[f"capital_cost_{index.split('_')[0]}"] = optimum.rename(
            {index.split("_")[0]: "fixed_cost", index: "total_cost"}
        )
    results_heuristic = pd.DataFrame.from_dict(heuristic_optima)
    return results_heuristic


results24_heuristic = get_heuristic_results(dfLCOH23, dfLCOH24)
delta_results = results24_heuristic - results24


# Note: not using a cf of 1 would be more realistic but still would implicitly assume one can somehow time downtime to fall perfectly into the hours of high prices.
def get_baseload_results(dfLCOH, baseload_cf=1):
    m = dfLCOH.query("capacity_factor >= @baseload_cf").capacity_factor.idxmax()
    columns = [f"capital_cost_{c}" for c in fixed_costs.index]
    common_values = pd.concat(
        [
            dfLCOH.loc[m].loc[
                ["price_max", "variable_cost_max", "average_cost", "capacity_factor"]
            ]
        ]
        * 3,
        axis=1,
    )
    common_values.columns = columns
    fixed_cost_row = dfLCOH.loc[m].loc[[f"{c}" for c in fixed_costs.index]].to_frame().T
    fixed_cost_row.columns = columns
    fixed_cost_row.index = ["fixed_costs"]
    total_cost_row = (
        dfLCOH.loc[m].loc[[f"{c}_total" for c in fixed_costs.index]].to_frame().T
    )
    total_cost_row.columns = columns
    total_cost_row.index = ["total_costs"]

    baseload_results = pd.concat([common_values, fixed_cost_row, total_cost_row])
    return baseload_results


results24_baseload = get_baseload_results(dfLCOH24)


def calculate_emissions_limited_LCOH(dfCountry):
    pass


# temporary output
print("Results in EUR/MWh H2 LHV")
print("perfect_foresight: 2023")
print(results23, "\n")
print("perfect_foresight: 2024")
print(results24, "\n")
print("heuristic: 2024")
print(results24_heuristic, "\n")
print("delta heuristic vs perfect foresight: 2024")
print(delta_results, "\n")
print("baseload: 2024")
print(results24_baseload)
