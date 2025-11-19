import pandas as pd
import os
from datetime import datetime
from nemosis import dynamic_data_compiler
from pathlib import Path

# This block is for linters and IDEs. It will not be executed by Snakemake.
if "snakemake" not in globals():
    from _stubs import snakemake


def download_data(snakemake):
    cache_dir = Path(snakemake.params.nemosis_cache_dir)
    generator_file_path = cache_dir / "NEM Registration and Exemption List.xlsx"

    if not os.path.exists(generator_file_path):
        print(f"Please ensure the generator excel file is at: {generator_file_path}")
        raise FileNotFoundError(
            "Missing manual config file. See readme for download instructions."
        )

    print("Reading generator map...")
    gen_info = pd.read_excel(
        generator_file_path, sheet_name="PU and Scheduled Loads", engine="openpyxl"
    )
    gen_info.columns = gen_info.columns.str.strip().str.replace("\n", " ")

    cols_to_keep = ["DUID", "Region", "Fuel Source - Primary"]
    map_df = gen_info[cols_to_keep].copy()

    print(f"Successfully mapped {len(map_df)} generator units.")

    # Parse the YYYYMMDD date from config and format for nemosis
    start_dt = datetime.strptime(snakemake.params.start_date, "%Y%m%d")
    end_dt = datetime.strptime(snakemake.params.end_date, "%Y%m%d")

    # Format it to what nemosis expects: YYYY/MM/DD HH:MM:SS
    start_time = start_dt.strftime("%Y/%m/%d") + " 00:00:00"
    end_time = end_dt.strftime("%Y/%m/%d") + " 23:59:59"

    # SCADA = Supervisory Control and Data Acquisition = actual physical generation data
    print("Fetching generation...")
    scada = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCH_UNIT_SCADA",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "DUID", "SCADAVALUE"],
        fformat="feather",
        keep_csv=True,
        rebuild=False,
    )

    print("Fetching load...")
    demand = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCHREGIONSUM",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "REGIONID", "TOTALDEMAND"],
        fformat="feather",
        keep_csv=True,
        rebuild=False,
    )

    # Ensure Datetimes
    scada["SETTLEMENTDATE"] = pd.to_datetime(scada["SETTLEMENTDATE"])
    demand["SETTLEMENTDATE"] = pd.to_datetime(demand["SETTLEMENTDATE"])

    # Pivot Generation
    gen_pivot = (
        scada.merge(map_df, on="DUID", how="left")
        .groupby(["SETTLEMENTDATE", "Region", "Fuel Source - Primary"])["SCADAVALUE"]
        .sum()
        .reset_index()
        .pivot_table(
            index="SETTLEMENTDATE",
            columns=["Region", "Fuel Source - Primary"],
            values="SCADAVALUE",
        )
        .fillna(0)
    )

    # Pivot Demand
    dem_pivot = demand.pivot_table(
        index="SETTLEMENTDATE", columns="REGIONID", values="TOTALDEMAND"
    )

    gen_pivot = gen_pivot.clip(lower=0)

    # Resample to Hourly (Mean)
    # NEM data is 5-min/30-min; ENTSO-E is hourly.
    # generation = gen_pivot.resample('1h').mean()
    # load = dem_pivot.resample('1h').mean()
    generation = gen_pivot
    load = dem_pivot

    # Calculate Region-Specific Variables using a dictionary comprehension
    # This is faster and cleaner than the manual list append loop
    regions = ["NSW1", "VIC1", "QLD1", "SA1", "TAS1"]
    frames = {}

    for r in regions:
        if r in load.columns:
            # Extract Demand and Generation
            d = load[r]
            # Use .get() to handle cases where a region might miss a fuel type (e.g., no solar in TAS1 historically)
            w = generation.get((r, "Wind"), pd.Series(0, index=generation.index))
            s = generation.get((r, "Solar"), pd.Series(0, index=generation.index))

            # Create a DataFrame for this region
            region_df = pd.DataFrame(index=generation.index)
            region_df["wind_onshore"] = w
            region_df["solar"] = s  # User 'vre' equivalent
            region_df["residual"] = d - (w + s)

            # Add to dictionary
            frames[r] = region_df

    # Keys become the top level (Region), Columns become the second level (Variable)
    df_final = pd.concat(frames, axis=1, names=["Region", "Variable"])

    df_final.to_feather(snakemake.output[0])

    print(f"Successfully saved NEM data to {snakemake.output[0]}")


if __name__ == "__main__":
    download_data(snakemake)
