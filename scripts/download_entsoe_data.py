import argparse
import pandas as pd
import os
import entsoe
from dotenv import load_dotenv
from pathlib import Path

def get_enabled_areas():
    """Reads the areas.csv file and returns a list of enabled area codes."""
    areas_df = pd.read_csv("areas.csv")
    return areas_df[areas_df["enabled"]]["area_code"].tolist()

AREAS = get_enabled_areas()


# Keeping the api key out of the code is good practice. Instead:
def get_entsoe_client():
    """Initializes and returns the EntsoePandasClient."""
    load_dotenv()
    API_KEY = os.environ.get("ENTSOE_API_KEY")
    if not API_KEY:
        raise ValueError(
            "ENTSOE_API_KEY environment variable not set. Please set it in your .env file."
        )

    client = entsoe.EntsoePandasClient(api_key=API_KEY)
    return client


def download_price_data(client, areas, start, end, output_file):
    """Downloads day-ahead prices for a list of areas and saves to a feather file."""
    print("Downloading day-ahead prices...")
    dfs_price = []
    for country_code in areas:
        print(country_code, end=" ")
        try:
            data = client.query_day_ahead_prices(country_code, start=start, end=end)
        except Exception as e:
            print('ERROR: ', repr(e))
        print("ok", end=" | ")
        data.name = country_code  # returns a series, so .name not .columns
        dfs_price.append(data)

    df_price = pd.concat(dfs_price, axis=1)

    # Connection errors can cause duplicated columns, save the later one
    df_price = df_price.loc[:, ~df_price.columns.duplicated(keep="last")]


    df_price.to_feather(output_file)
    print(f"Successfully saved price data to {output_file}")


def download_load_forecast_data(client, areas, start, end, output_file):
    """Downloads load forecast data for a list of areas and saves to a feather file."""
    print("Downloading load forecast data...")
    dfs_power_fc = []
    for country_code in areas:
        print(country_code, end=" ")
        try:
            data = client.query_load_forecast(country_code, start=start, end=end)
        except Exception as e:
            print('ERRROR', repr(e))
        print("ok", end=" | ")
        data.columns = [country_code]  # returns a dataframe, so .columns not .name
        dfs_power_fc.append(data)

    df_power_fc = pd.concat(dfs_power_fc, axis=1)

    # Connection errors can cause duplicated columns, save the later one
    df_power = df_power_fc.loc[:, ~df_power_fc.columns.duplicated(keep="last")]

    df_power.to_feather(output_file)
    print(f"Successfully saved load forecast data to {output_file}")


def download_load_data(client, areas, start, end, output_file):
    """Downloads actual load data for a list of areas and saves to a feather file."""
    print("Downloading actual load data...")
    dfs_power = []
    for country_code in areas:
        print(country_code, end=" ")
        try:
            data = client.query_load(country_code, start=start, end=end)
        except Exception as e:
            print('ERRROR', repr(e))
        print("ok", end=" | ")
        data.columns = [country_code]  # returns a dataframe, so .columns not .name
        dfs_power.append(data)

    df_power = pd.concat(dfs_power, axis=1)

    # Connection errors can cause duplicated columns, save the later one
    df_power = df_power.loc[:, ~df_power.columns.duplicated(keep="last")]

    df_power.to_feather(output_file)
    print(f"Successfully saved actual load data to {output_file}")


def download_vre_forecast_data(client, areas, start, end, output_file):
    print("Downloading wind and solar forecast data...")

    vre_names = {"Solar", "Wind Onshore", "Wind Offshore"}
    vre_new_names = {
        "Solar": "solar_forecast",
        "Wind Onshore": "wind_onshore_forecast",
        "Wind Offshore": "wind_offshore_forecast",
    }
    dfs_vre_tup = []
    for country_code in areas:
        print(country_code, end=" ")
        try:
            data = client.query_wind_and_solar_forecast(
                country_code, start=start, end=end
            )
        except Exception as e:
            print(repr(e))
        print("ok", end=" | ")
        dfs_vre_tup.append((data, country_code))

    # take care of duplicted countries
    seen_codes = set()
    dfs_vre_tup_temp = []
    for data, country_code in reversed(dfs_vre_tup):
        if country_code in seen_codes:
            continue
        seen_codes.add(country_code)
        dfs_vre_tup_temp.append((data, country_code))
    dfs_vre_tup = list(reversed(dfs_vre_tup_temp))
        

    dfs_vre = []
    for data, country_code in dfs_vre_tup:
        missing = vre_names - set(data.columns)
        data = data.copy()
        for m in missing:
            data[m] = 0.0
        index = pd.MultiIndex.from_tuples(
            [(country_code, vre_new_names[vre]) for vre in data.columns]
        )  # list of tuples
        data.columns = index
        dfs_vre.append(data)
    df_vre = pd.concat(dfs_vre, axis=1)
    df_vre.to_feather(output_file)
    print(f"Successfully saved VRE data to {output_file}")

def download_generation_data(client, areas, start, end, output_file):
    print("Downloading actual generation data...")

    gen_names = {
        "Biomass": "biomass",
        "Energy storage": "energy_storage",
        "Fossil Brown coal/Lignite": "brown_coal",
        "Fossil Coal-derived gas": "coal_gas",
        "Fossil Gas": "gas",
        "Fossil Hard coal": "hard_coal",
        "Fossil Oil": "oil",
        "Fossil Oil shale": "oil_shale",
        "Fossil Peat": "peat",
        "Geothermal": "geothermal",
        "Hydro Pumped Storage": "pumped_storage",
        "Hydro Run-of-river and poundage": "hydro_river",
        "Hydro Water Reservoir": "hydro_reservoir",
        "Marine": "marine",
        "Nuclear": "nuclear",
        "Other": "other",
        "Other renewable": "other_re",
        "Solar": "solar",
        "Waste": "waste",
        "Wind Offshore": "wind_offshore",
        "Wind Onshore": "wind_onshore",
    }
    gen_names_production = {k + "_Actual Aggregated": v for k, v in gen_names.items()}
    gen_names_consumption = {k + "_Actual Consumption": v + "_cons" for k, v in gen_names.items()}
    gen_names = gen_names_production | gen_names_consumption

    dfs_gen_tup = []
    for country_code in areas:
        print(country_code, end=" ")
        try:
            data = client.query_generation(
                country_code, start=start, end=end
            )
        except Exception as e:
            print(repr(e))
        if data.columns.nlevels == 2:
            data.columns = ['_'.join(col) for col in data.columns.values]
        else:
            data.columns = ['_'.join([col, 'Actual Aggregated']) for col in data.columns]
        cols_to_change = data.filter(regex='_cons$', axis=1).columns
        data.loc[:, cols_to_change] = data.loc[:, cols_to_change] * -1
        dfs_gen_tup.append((data, country_code))
        print("ok", end=" | ")


    # take care of duplicted countries
    seen_codes = set()
    dfs_gen_tup_temp = []
    for data, country_code in reversed(dfs_gen_tup):
        if country_code in seen_codes:
            continue
        seen_codes.add(country_code)
        dfs_gen_tup_temp.append((data, country_code))
    dfs_gen_tup = list(reversed(dfs_gen_tup_temp))
        

    dfs_gen = []
    for data, country_code in dfs_gen_tup:
        missing = gen_names.keys() - set(data.columns)
        data = data.copy()
        for m in missing:
            data[m] = 0.0
        index = pd.MultiIndex.from_tuples(
            [(country_code, gen_names[gen]) for gen in data.columns]
        )  # list of tuples
        data.columns = index
        dfs_gen.append(data)
    df_gen = pd.concat(dfs_gen, axis=1)
    df_gen.to_feather(output_file)
    print(f"Successfully saved actual generation data to {output_file}")




if __name__ == "__main__":
    # Load environment variables from .env file
    load_dotenv()

    parser = argparse.ArgumentParser(description="Download data from the ENTSO-E API.")
    parser.add_argument(
        "data_type",
        choices=["prices", "load_forecast", "load_actual", "vre", "generation"],
        help="The type of data to download.",
    )
    parser.add_argument("output_file", type=Path, help="The path to the output file.")
    parser.add_argument("--start", required=True, help="Start date in YYYYMMDD format.")
    parser.add_argument("--end", required=True, help="End date in YYYYMMDD format.")
    parser.add_argument(
        "--areas",
        nargs="+",
        default=AREAS,
        help="List of area codes to download data for.",
    )

    args = parser.parse_args()

    client = get_entsoe_client()
    start_ts = pd.Timestamp(args.start, tz="Europe/Brussels")
    end_ts = pd.Timestamp(args.end, tz="Europe/Brussels")

    if args.data_type == "prices":
        download_price_data(client, args.areas, start_ts, end_ts, args.output_file)
    elif args.data_type == "load_forecast":
        download_load_forecast_data(client, args.areas, start_ts, end_ts, args.output_file)
    elif args.data_type == "load_actual":
        download_load_data(client, args.areas, start_ts, end_ts, args.output_file)
    elif args.data_type == "vre":
        download_vre_forecast_data(client, args.areas, start_ts, end_ts, args.output_file)
    elif args.data_type == "generation":
        download_generation_data(client, args.areas, start_ts, end_ts, args.output_file)