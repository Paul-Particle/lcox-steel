import pandas as pd
import plotly.graph_objects as go
from pathlib import Path

# This block is for linters and IDEs. It will not be executed by Snakemake.
if "snakemake" not in globals():
    from _stubs import snakemake

def verify_data_availability(snakemake):
    """
    Analyzes processed data to visualize data availability for each area and data type.
    """
    print("Verifying data availability...")
    processed_data_file = Path(snakemake.input[0])
    output_file = Path(snakemake.output[0])
    output_dir = output_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_feather(processed_data_file)
    areas = df.columns.get_level_values(0).unique()
    data_types = df.columns.get_level_values(1).unique()

    fig = go.Figure()

    # Define a color map for data types
    color_map = {
        'prices': 'blue',
        'load_forecast': 'green',
        'load_actual': 'red',
        'vre': 'orange',
        'generation': 'purple',
        'crossborder': 'brown'
    }

    y_offset = 0
    for area in areas:
        area_df = df[area]
        for i, data_type in enumerate(data_types):
            if data_type in area_df.columns:
                start_date = area_df[data_type].first_valid_index()
                end_date = area_df[data_type].last_valid_index()
                if pd.notna(start_date) and pd.notna(end_date):
                    fig.add_trace(go.Scatter(
                        x=[start_date, end_date],
                        y=[y_offset + i, y_offset + i],
                        mode='lines',
                        name=f"{area} - {data_type}",
                        line=dict(color=color_map.get(data_type, 'black'), width=10),
                        showlegend=True
                    ))
        y_offset += len(data_types)

    fig.update_layout(
        title="Data Availability by Area and Data Type",
        xaxis_title="Date",
        yaxis_title="Area / Data Type",
        yaxis=dict(
            tickmode='array',
            tickvals=[i * len(data_types) + len(data_types)/2 for i in range(len(areas))],
            ticktext=areas
        ),
        showlegend=True
    )

    fig.write_html(output_file)
    print(f"Successfully saved data availability plot to {output_file}")

if __name__ == "__main__":
    verify_data_availability(snakemake)