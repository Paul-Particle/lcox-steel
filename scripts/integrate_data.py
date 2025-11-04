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

    # Group input files by data_type from the file path
    files_by_type = defaultdict(list)
    for file_path in snakemake.input:
        p = Path(file_path)
        data_type = p.stem
        files_by_type[data_type].append(p)

    for data_type, file_paths in files_by_type.items():
        print(f"Processing data type: {data_type}")
        dfs = []
        for file_path in file_paths:
            df = pd.read_feather(file_path)
            dfs.append(df)
        
        if dfs:
            df_integrated = pd.concat(dfs, axis=1)
            df_integrated = df_integrated.loc[:, ~df_integrated.columns.duplicated(keep="last")]

            output_file = output_dir / f"{data_type}.feather"
            df_integrated.to_feather(output_file)
            print(f"Successfully saved integrated data to {output_file}")

if __name__ == "__main__":
    integrate_data(snakemake)
