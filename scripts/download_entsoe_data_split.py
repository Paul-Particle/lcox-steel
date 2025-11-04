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

def get_entsoe_client():
    """Initializes and returns the EntsoePandasClient."""
    load_dotenv()
    API_KEY = os.environ.get("ENTSOE_API_KEY")
    if not API_KEY:
        raise ValueError(
            "ENTSOE_API_KEY environment variable not set. Please set it in your .env file."
        )

    client = entsoe.EntsoePandasClient(api_key=API_KEY) # pyright: ignore[reportPrivateImportUsage]
    return client


def download_price_data(client, area, start, end, output_file):
    """Downloads day-ahead prices for a single area and saves to a feather file."""
    print(f"Downloading day-ahead prices for {area}...")
    try:
        data = client.query_day_ahead_prices(area, start=start, end=end)
    except Exception as e:
        print('ERROR: ', repr(e))
    else:
        print("ok")
        data.name = area
        df = data.to_frame()
        df.to_feather(output_file)
        print(f"Successfully saved price data to {output_file}")


def download_load_forecast_data(client, area, start, end, output_file):
    """Downloads load forecast data for a single area and saves to a feather file."""
    print(f"Downloading load forecast data for {area}...")
    try:
        data = client.query_load_forecast(area, start=start, end=end)
    except Exception as e:
        print('ERROR', repr(e))
    else:
        print("ok")
        data.columns = [area]
        data.to_feather(output_file)
        print(f"Successfully saved load forecast data to {output_file}")


def download_load_data(client, area, start, end, output_file):
    """Downloads actual load data for a single area and saves to a feather file."""
    print(f"Downloading actual load data for {area}...")
    try:
        data = client.query_load(area, start=start, end=end)
    except Exception as e:
        print('ERROR', repr(e))
    else:
        print("ok")
        data.columns = [area]
        data.to_feather(output_file)
        print(f"Successfully saved actual load data to {output_file}")


def download_vre_forecast_data(client, area, start, end, output_file):
    print(f"Downloading wind and solar forecast data for {area}...")

    vre_new_names = {
        "Solar": "solar_forecast",
        "Wind Onshore": "wind_onshore_forecast",
        "Wind Offshore": "wind_offshore_forecast",
    }
    try:
        data = client.query_wind_and_solar_forecast(
            area, start=start, end=end
        )
    except Exception as e:
        print(repr(e))
    else:
        print("ok")
        index = pd.MultiIndex.from_tuples(
            [(area, vre_new_names[vre]) for vre in data.columns]
        )
        data.columns = index
        data.to_feather(output_file)
        print(f"Successfully saved VRE data to {output_file}")

def download_generation_data(client, area, start, end, output_file):
    print(f"Downloading actual generation data for {area}...")

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

    try:
        data = client.query_generation(
            area, start=start, end=end
        )
    except Exception as e:
        print(repr(e))
    else:
        if data.columns.nlevels == 2:
            data.columns = ['_'.join(col) for col in data.columns.values]
        else:
            data.columns = ['_'.join([col, 'Actual Aggregated']) for col in data.columns]
        cols_to_change = data.filter(regex='Consumption$', axis=1).columns
        data.loc[:, cols_to_change] = data.loc[:, cols_to_change] * -1
        
        index = pd.MultiIndex.from_tuples(
            [(area, gen_names[gen]) for gen in data.columns]
        )
        data.columns = index
        data.to_feather(output_file)
        print(f"Successfully saved actual generation data to {output_file}")

def download_export_import_data(client, area, start, end, output_file):
    print(f"Downloading crossborder export-import flow data for {area}...")
    dfs_crossborder = []
    
    print(f"{area} in")
    try: 
        data_in = client.query_physical_crossborder_allborders(area, start=start, end=end, export=False, per_hour=False)
    except Exception as e:
        print(repr(e))
    else:
        print("ok")
        data_in.columns = ['from_' + col for col in data_in.columns]
        dfs_crossborder.append(data_in)

    print(f"{area} out")
    try: 
        data_out = client.query_physical_crossborder_allborders(area, start=start, end=end, export=True, per_hour=False)
    except Exception as e:
        print(repr(e))
    else:
        print("ok")
        data_out.columns = ['to_' + col for col in data_out.columns]
        dfs_crossborder.append(data_out * -1)

    if dfs_crossborder:
        df_crossborder = pd.concat(dfs_crossborder, axis=1)
        index = pd.MultiIndex.from_tuples(
            [(area, col) for col in df_crossborder.columns]
        )
        df_crossborder.columns = index
        df_crossborder.to_feather(output_file)
        print(f"Successfully saved crossborder data to {output_file}")


if __name__ == "__main__":
    # Load environment variables from .env file
    load_dotenv()

    parser = argparse.ArgumentParser(description="Download data from the ENTSO-E API.")
    parser.add_argument(
        "data_type",
        choices=["prices", "load_forecast", "load_actual", "vre", "generation", "crossborder"],
        help="The type of data to download.",
    )
    parser.add_argument("output_dir", type=Path, help="The path to the output directory.")
    parser.add_argument("area", help="Area code to download data for.")
    parser.add_argument("year", type=int, help="Year to download data for.")
    parser.add_argument("month", type=int, help="Month to download data for.")


    args = parser.parse_args()

    client = get_entsoe_client()
    start_ts = pd.Timestamp(f"{args.year}-{args.month}-01", tz="Europe/Brussels")
    end_ts = start_ts + pd.offsets.MonthEnd(0)

    output_path = args.output_dir / args.area / f"{args.year}-{args.month:02d}"
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / f"{args.data_type}.feather"


    print(f"Downloading data from {start_ts} to {end_ts}")
    print(f"Area: {args.area}")

    if args.data_type == "prices":
        download_price_data(client, args.area, start_ts, end_ts, output_file)
    elif args.data_type == "load_forecast":
        download_load_forecast_data(client, args.area, start_ts, end_ts, output_file)
    elif args.data_type == "load_actual":
        download_load_data(client, args.area, start_ts, end_ts, output_file)
    elif args.data_type == "vre":
        download_vre_forecast_data(client, args.area, start_ts, end_ts, output_file)
    elif args.data_type == "generation":
        download_generation_data(client, args.area, start_ts, end_ts, output_file)
    elif args.data_type == "crossborder":
        download_export_import_data(client, args.area, start_ts, end_ts, output_file)