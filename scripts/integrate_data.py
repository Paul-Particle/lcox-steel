import argparse
import pandas as pd
from pathlib import Path

def integrate_data(input_dir, output_dir, data_types, areas, year_months):
    """
    Integrates downloaded data from a structured directory into monolithic feather files.
    """
    print("Integrating data...")
    output_dir.mkdir(parents=True, exist_ok=True)

    for data_type in data_types:
        print(f"Processing data type: {data_type}")
        dfs = []
        for area in areas:
            for year_month in year_months:
                file_path = input_dir / area / year_month / f"{data_type}.feather"
                if file_path.exists():
                    df = pd.read_feather(file_path)
                    dfs.append(df)
        
        if dfs:
            # Concatenate all dataframes for the current data type
            df_integrated = pd.concat(dfs, axis=1)
            
            # Connection errors can cause duplicated columns, save the later one
            df_integrated = df_integrated.loc[:, ~df_integrated.columns.duplicated(keep="last")]

            output_file = output_dir / f"{data_type}.feather"
            df_integrated.to_feather(output_file)
            print(f"Successfully saved integrated data to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrate ENTSO-E data.")
    parser.add_argument("input_dir", type=Path, help="Path to the input directory.")
    parser.add_argument("output_dir", type=Path, help="Path to the output directory.")
    parser.add_argument("--data_types", nargs="+", required=True, help="List of data types to integrate.")
    parser.add_argument("--areas", nargs="+", required=True, help="List of area codes to integrate.")
    parser.add_argument("--year_months", nargs="+", required=True, help="List of year-month combinations to integrate (e.g., 2025-10).")
    
    args = parser.parse_args()
    
    integrate_data(args.input_dir, args.output_dir, args.data_types, args.areas, args.year_months)
