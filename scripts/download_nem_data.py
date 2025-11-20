import pandas as pd
import os
from datetime import datetime
from nemosis import dynamic_data_compiler
from pathlib import Path

# This block is for linters and IDEs. It will not be executed by Snakemake.
if "snakemake" not in globals():
    from _stubs import snakemake


def download_price_data(start_time, end_time, cache_dir, rebuild):
    """Downloads price data or gets it from the cached feather files or csv files if rebuild=True."""
    print("Fetching prices...")
    prices = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCHPRICE",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "REGIONID", "RRP"], # regional reference price in Australian Dollar
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    prices["SETTLEMENTDATE"] = pd.to_datetime(prices["SETTLEMENTDATE"])
    prices = prices.pivot_table(
        index="SETTLEMENTDATE", columns="REGIONID", values="RRP"
    )
    prices.columns = pd.MultiIndex.from_tuples([(col, 'price') for col in prices.columns])
    prices.index.name = None
    return prices

def download_generation_data(start_time, end_time, cache_dir, rebuild):
    """Downloads and processes generation data."""
    generator_file_path = cache_dir / "NEM Registration and Exemption List.xlsx"

    if not os.path.exists(generator_file_path):
        print(f"Please ensure the generator excel file is at: {generator_file_path}")
        raise FileNotFoundError(
            "Missing manual excel config file. See readme for download instructions."
        )

    generator_info = pd.read_excel(
        generator_file_path, sheet_name="PU and Scheduled Loads", engine="openpyxl"
    )
    generator_info.columns = generator_info.columns.str.strip().str.replace("\n", " ")

    
    def determine_gen_type(df):
        dfi = df['gen_info']
        type_regex = {
            'hard_coal' : 'black coal',
            'brown_coal' : 'brown coal',
            'coal_gas' : 'coal seam methane|coal mine gas',
            'waste' : 'methane',
            'oil' : 'diesel|ethane|kerosene',
            'biomass' : 'bagasse|biomass|biogas',
            'gas' : 'natural gas|natrual gas',
            'solar' : 'solar',
            'energy_storage' : 'battery',
            'pumped_storage' : 'pump storage|^- -$',
            'other_re' : 'sewerage',
            'wind_onshore' : 'wind - onshore|wind',
            'hydro_river' : 'run of river',
            'hydro' : 'hydro',
        }
        for t, r in type_regex.items():
            m = dfi.str.contains(r)
            df.loc[m, 'gen_type'] = t
            df.loc[m, 'gen_info_temp'] = df.loc[m, 'gen_info']
            df.loc[m, 'gen_info'] = ''
        df = df.drop(columns=['gen_info'])
        df.columns = ['gen_type', 'gen_info']
        return df.loc[:,['gen_info', 'gen_type']]


    nem_gen_names = dict(generator_info.copy()
        .loc[:, ["DUID", "Fuel Source - Descriptor", "Technology Type - Descriptor"]]
        .dropna()
        .pivot_table(
            index=["Fuel Source - Descriptor"],
            columns=["Technology Type - Descriptor"],
            values="DUID",
            aggfunc="count",
            fill_value=0,
        )
        .stack() # getting all combinations of Fuel Source and Technology Type
        .loc[lambda s: s>0] # filter for existing combinations
        .reset_index() # put the columns back into the df
        .assign(gen_info = lambda df: df["Fuel Source - Descriptor"].str.lower() + ' ' + df["Technology Type - Descriptor"].str.lower())
        .drop(columns=["Technology Type - Descriptor", "Fuel Source - Descriptor", 0]) # 0 column is number of plants for each combination, i.e. not relevant
        .pipe(determine_gen_type) # apply regex-based categorization
        .values
    )
    # Ensure that there are no nan values in the dict, to catch changes to future versions of the excel.
    # This can happen if a new 'Fuel Source - Descriptor' and 'Technology Type - Descriptor' combination
    # is not covered by the regex in `determine_gen_type`.
    assert not any(pd.isna(v) for v in nem_gen_names.values()), "NaN values found in nem_gen_names mapping. Please update the regex in 'determine_gen_type'."
    

    # final mapping of generator ID to generator type
    map_df = (generator_info.copy()
        .loc[:, ["DUID", "Region", "Fuel Source - Descriptor", "Technology Type - Descriptor"]]
        .assign(gen_info = lambda df: df["Fuel Source - Descriptor"].str.lower() + ' ' + df["Technology Type - Descriptor"].str.lower())
        .drop(columns=["Technology Type - Descriptor", "Fuel Source - Descriptor"])
        .assign(gen_type = lambda df: df["gen_info"].map(nem_gen_names))
    )


    print("Fetching generation...")
    scada = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCH_UNIT_SCADA",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "DUID", "SCADAVALUE"], # SCADA = generation in MW
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    scada["SETTLEMENTDATE"] = pd.to_datetime(scada["SETTLEMENTDATE"])

    generation = (
        scada.merge(map_df.dropna(subset=['gen_type']), on="DUID", how="left")
        .groupby(["SETTLEMENTDATE", "Region", "gen_type"])["SCADAVALUE"]
        .sum()
        .reset_index()
        .pivot_table(
            index="SETTLEMENTDATE",
            columns=["Region", "gen_type"],
            values="SCADAVALUE",
        )
        .fillna(0)
    )
    generation = generation.clip(lower=0)
    generation.index.name = None
    return generation


def download_load_data(start_time, end_time, cache_dir, rebuild):
    """Downloads and processes load data."""
    print("Fetching load...")
    load = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCHREGIONSUM",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "REGIONID", "TOTALDEMAND"],
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    load["SETTLEMENTDATE"] = pd.to_datetime(load["SETTLEMENTDATE"])

    # Pivot load
    load = load.pivot_table(
        index="SETTLEMENTDATE", columns="REGIONID", values="TOTALDEMAND"
    )
    load.columns = pd.MultiIndex.from_tuples([(col, 'load') for col in load.columns])
    load.index.name = None
    return load


def download_data(snakemake):
    cache_dir = Path(snakemake.params.nemosis_cache_dir)
    rebuild = snakemake.params.get("rebuild", False) # pyright: ignore[reportGeneralTypeIssues]
    # Parse the YYYYMMDD date from config and format for nemosis as YYYY/MM/DD HH:MM:SS
    start_time = datetime.strptime(snakemake.params.start_date, "%Y%m%d")
    start_time = start_time.strftime("%Y/%m/%d") + " 00:00:00"
    end_time = datetime.strptime(snakemake.params.end_date, "%Y%m%d")
    end_time = end_time.strftime("%Y/%m/%d") + " 23:59:59"
    
    prices = download_price_data(start_time, end_time, cache_dir, rebuild)
    generation = download_generation_data(start_time, end_time, cache_dir, rebuild)
    load = download_load_data(start_time, end_time, cache_dir, rebuild)

    df_combined = pd.concat([prices, load, generation], axis=1)

    # Calculate Region-Specific Variables using a dictionary comprehension
    # This is faster and cleaner than the manual list append loop
    areas = ["NSW1", "VIC1", "QLD1", "SA1", "TAS1"]
    frames = {}

    for a in areas:
        if (a, 'load') in df_combined.columns:
            # Extract Demand and Generation
            d = df_combined[(a, 'load')]
            # Use .get() to handle cases where a region might miss a fuel type (e.g., no solar in TAS1 historically)
            w = df_combined.get((a, "wind_onshore"), pd.Series(0, index=df_combined.index))
            s = df_combined.get((a, "solar"), pd.Series(0, index=df_combined.index))

            # Create a DataFrame for this region
            region_df = pd.DataFrame(index=generation.index)
            region_df["wind"] = w 
            region_df["residual"] = d - (w + s)

            # Combine with other data for the region
            frames[a] = pd.concat([df_combined[a], region_df], axis=1)

    # Keys become the top level (Region), Columns become the second level (Variable)
    df_final = pd.concat(frames.values(), axis=1, keys=frames.keys())
    
    df_final.index = df_final.index.tz_localize(None) # pyright: ignore[reportAttributeAccessIssue]
    df_final.to_feather(snakemake.output[0])

    print(f"Successfully saved NEM data to {snakemake.output[0]}")


if __name__ == "__main__":
    download_data(snakemake)
