import pandas as pd
from pathlib import Path
import country_converter as coco

if "snakemake" not in globals():
    from _stubs import snakemake

def get_enabled_areas(areas_file_path):
    """Reads the areas.csv file and returns a list of enabled area codes."""
    areas_df = pd.read_csv(areas_file_path)
    return areas_df[areas_df["enabled"]]["area_code"].tolist()

def process_data(snakemake):
    """
    Processes data into a single dataframe.
    """
    print("Processing data...")
    
    # Use the input from snakemake to get the enabled areas
    AREAS = get_enabled_areas(snakemake.input.areas_config)

    # Access inputs by name, not by index, for robustness
    df_price = pd.read_feather(snakemake.input.prices)
    df_power_fc = pd.read_feather(snakemake.input.load_forecast)
    df_power_ac = pd.read_feather(snakemake.input.load_actual)
    df_vre = pd.read_feather(snakemake.input.vre)
    df_gen = pd.read_feather(snakemake.input.generation)
    df_xim = pd.read_feather(snakemake.input.crossborder)

    # forward fill interpolation, as some data might have higher resolution
    df_price = df_price.ffill(limit=3).fillna(0.0)
    df_power_fc = df_power_fc.ffill(limit=3).fillna(0.0)
    df_power_ac = df_power_ac.ffill(limit=3).fillna(0.0)
    df_vre = df_vre.ffill(limit=3).fillna(0.0)
    df_gen = df_gen.ffill(limit=3).fillna(0.0)
    df_xim = df_xim.ffill(limit=3).fillna(0.0)

    # Create MultiIndex columns for merging
    df_p = df_price.copy()
    index_p = pd.MultiIndex.from_tuples([(col, 'price') for col in df_p.columns])
    df_p.columns = index_p

    df_d_fc = df_power_fc.copy()
    index_d = pd.MultiIndex.from_tuples([(col, 'demand_forecast') for col in df_d_fc.columns])
    df_d_fc.columns = index_d

    df_d_ac = df_power_ac.copy()
    index_d = pd.MultiIndex.from_tuples([(col, 'demand') for col in df_d_ac.columns])
    df_d_ac.columns = index_d

    df_v = df_vre.copy()

    df_g = df_gen.copy()

    df_x = df_xim.copy()

    # Combine the dataframes
    df_EU = pd.concat([df_p, df_d_fc, df_d_ac, df_v, df_g, df_x], axis=1)

    # Process each area
    dfs_EU = {}
    for area in [col for col in df_EU.columns.get_level_values(0).unique() if col in AREAS]:
        # Select data for the current area, resulting in a DataFrame with single-level columns
        df_area_data = df_EU[area].copy()

        # Perform the assign operations
        df_area_processed = (df_area_data
            .sort_index(axis=1)
            .assign(hour          = lambda x: x.index.hour)
            .assign(doy           = lambda x: x.index.dayofyear)
            .assign(doy_season    = lambda x: ((x.doy + 8*30+4)%365)+1)
            .assign(dt            = lambda x: x.index)
            .assign(wind_forecast = lambda x: x.get('wind_onshore_foreacast', 0) + x.get('wind_offshore_forecast', 0))
            .assign(vre_forecast  = lambda x: x.get('wind_forecast', 0) + x.get('solar_forecast', 0))
            .assign(residual      = lambda x: x.get('demand_forecast', 0) - x.get('vre_forecast', 0))
            .ffill(limit=3)
        )
        
        # Create a new MultiIndex for the columns of the processed area data
        # The first level is the area, the second level are the new and old metrics
        new_columns = pd.MultiIndex.from_product([[area], df_area_processed.columns])
        df_area_processed.columns = new_columns
        
        dfs_EU[area] = df_area_processed

    # Combine all areas into a single dataframe horizontally
    df_EU_all = pd.concat(dfs_EU.values(), axis=1)
    
    # Save the processed data
    df_EU_all.to_feather(snakemake.output[0])
    print(f"Successfully saved processed data to {snakemake.output[0]}")


if __name__ == "__main__":
    process_data(snakemake)
