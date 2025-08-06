import argparse
import pandas as pd
import os
import entsoe
from dotenv import load_dotenv
from pathlib import Path
import country_converter as coco

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
            print(repr(e))
        print("ok", end=" | ")
        data.name = country_code  # returns a series, so .name not .columns
        dfs_price.append(data)

    df_price = pd.concat(dfs_price, axis=1)

    df_price.to_feather(output_file)
    print(f"Successfully saved price data to {output_file}")


def download_load_data(client, areas, start, end, output_file):
    """Downloads load forecast data for a list of areas and saves to a feather file."""
    print("Downloading load forecast data...")
    dfs_power = []
    for country_code in areas:
        print(country_code, end=" ")
        try:
            # continue
            data = client.query_load_forecast(country_code, start=start, end=end)
        except Exception as e:
            print(repr(e))
            # continue
        print("ok", end=" | ")
        data.columns = [country_code]  # returns a dataframe, so .columns not .name
        dfs_power.append(data)

    df_power = pd.concat(dfs_power, axis=1)
    df_power.to_feather(output_file)
    print(f"Successfully saved load data to {output_file}")


def download_vre_data(client, areas, start, end, output_file):
    """Downloads"""
    print("Downloading wind and solar forecast data...")

    vre_names = {"Solar", "Wind Onshore", "Wind Offshore"}
    vre_new_names = {
        "Solar": "solar",
        "Wind Onshore": "wind_onshore",
        "Wind Offshore": "wind_offshore",
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

    dfs_vre = []
    for data, country_code in dfs_vre_tup:
        missing = vre_names - set(data.columns)
        print(data.columns)
        print(country_code, missing)
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


if __name__ == "__main__":
    # Load environment variables from .env file
    load_dotenv()

    parser = argparse.ArgumentParser(description="Download data from the ENTSO-E API.")
    parser.add_argument(
        "data_type",
        choices=["prices", "load", "vre"],
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
    elif args.data_type == "load":
        download_load_data(client, args.areas, start_ts, end_ts, args.output_file)
    elif args.data_type == "vre":
        download_vre_data(client, args.areas, start_ts, end_ts, args.output_file)