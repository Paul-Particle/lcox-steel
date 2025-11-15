import pandas as pd
from pathlib import Path
from collections import defaultdict

if "snakemake" not in globals():
    from _stubs import snakemake

def integrate_data(snakemake):
    """
    Integrates downloaded data from a list of input files into monolithic feather files
    using the snakemake object.
    """
    print("Integrating data...")
    output_dir = Path(snakemake.output[0]).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # Group input files by area and then by data_type from the file path
    files_by_area_and_type = defaultdict(lambda: defaultdict(list))
    for file_path in snakemake.input:
        p = Path(file_path)
        area = p.parent.parent.name  # e.g., data/raw/DE/2023-01/prices.feather -> DE
        data_type = p.stem
        files_by_area_and_type[data_type][area].append(p)

    for data_type, areas in files_by_area_and_type.items():
        print(f"Processing data type: {data_type}")
        all_area_dfs = []
        for area, file_paths in areas.items():
            # For each area, concatenate the time-split files
            area_dfs = []
            for file_path in sorted(file_paths): # sort to be chronological
                df = pd.read_feather(file_path)
                area_dfs.append(df)
            
            if area_dfs:
                # Concatenate along rows (axis=0) for time series data
                df_area = pd.concat(area_dfs, axis=0)
                # Remove duplicate index entries, keeping the last one
                df_area = df_area[~df_area.index.duplicated(keep='last')]
                all_area_dfs.append(df_area)

        if all_area_dfs:
            # Concatenate all areas horizontally
            df_integrated = pd.concat(all_area_dfs, axis=1)
            df_integrated = df_integrated.sort_index()

            output_file = output_dir / f"{data_type}.feather"
            df_integrated.to_feather(output_file)
            print(f"Successfully saved integrated data to {output_file}")

if __name__ == "__main__":
    integrate_data(snakemake)
