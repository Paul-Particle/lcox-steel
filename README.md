# LCOX Steel: Green Hydrogen Cost Analysis

This project provides a data analysis pipeline to determine the optimal Levelized Cost of Hydrogen (LCOH) based on electricity market data or other supply scenarios. It uses a Snakemake workflow to download, process, and analyze data from two primary sources:

*   **ENTSO-E:** The European Network of Transmission System Operators for Electricity, for European market data.
*   **NEM:** The National Electricity Market of Australia, for Australian market data.

The ultimate goal is to model the costs associated with producing steel with green hydrogen.

## Project Structure

```
├── Snakefile           # Snakemake workflow definition
├── config.yaml         # Configuration for the Snakemake pipeline
├── environment.yml     # Conda environment definition
├── data/               # Raw and processed data
├── notebooks/          # Jupyter notebooks for exploration and analysis
├── results/            # Output files, such as plots and tables
└── scripts/            # Python scripts for data downloading, processing, and plotting
```

## Setup and Usage
0.  **Get API key:**

    See [here](https://transparencyplatform.zendesk.com/hc/en-us/articles/12845911031188-How-to-get-security-token). Send an email to transparency@entsoe.eu with “Restful API access” in the subject line. Indicate the email address you entered during registration in the email body. 

1.  **Create a .env file**
    ```
    touch .env
    ```
    Add the line ```ENTSOE_API_KEY=<your-api-key-here>```

2.  **Create the Conda environment:**

    ```bash
    conda env create -f environment.yaml
    ```

2.  **Activate the environment:**

    ```bash
    conda activate lcox-env
    ```

3.  **Configure the pipeline:**

    Edit the `config.yaml` file to specify the data download period and other parameters.

4.  **Run the Snakemake pipeline:**

    ```bash
    snakemake --cores 1
    ```

    This will execute the entire workflow, from data downloading to plot generation.

## Workflow

The `Snakefile` defines the following steps:

1.  **Download Data:** The pipeline fetches data from two separate sources.
    *   `download_entsoe_data.py`: Downloads data from the ENTSO-E API for specified European regions. It downloads data for each data type (prices, load, generation, etc.) into separate files for each month to manage API request sizes.
    *   `download_nem_data.py`: Downloads generation and demand data from Australia's National Electricity Market (NEM). It processes also processes this data directly to match ENTSO-E data (apart from the data resolution, which is 5 min) and calculates VRE (Variable Renewable Energy) and residual load for each Australian region. It relies on a static excel file containing information about generators which you can download [here](https://www.aemo.com.au/-/media/files/electricity/nem/participant_information/nem-registration-and-exemption-list.xlsx). It is also included in this repo, but may change in the future. Note that NEM data has no forecasts.

2.  **Integrate ENTSO-E Data:**
    *   `integrate_entsoe_data.py`: This script takes the monthly raw files downloaded from ENTSO-E and combines them into a single file for each data type (e.g., `prices.feather`, `generation.feather`). This prepares the data for processing.

3.  **Process ENTSO-E Data:**
    *   `process_entsoe_data.py`: Merges the various integrated ENTSO-E datasets (price, load, generation, etc.) into a single, wide-format DataFrame. It uses a multi-level column index where the top level is the country code. It also calculates additional useful metrics like `vre_forecast` (sum of wind and solar forecasts) and `residual_forecast` (load forecast minus VRE forecast).

4.  **Analyze and Optimize (LCOH Calculation):**
    *   `optimize_grid_capacity_factor.py`: This is the first core analysis script of the project. It is not currently part of the Snakemake workflow but can be run after the pipeline completes. It reads the processed data and performs the following:
        *   Calculates the Levelized Cost of Hydrogen (LCOH) based on electricity prices and configurable costs for electrolyzers (CAPEX, O&M, efficiency).
        *   Determines the optimal production schedule to minimize LCOH under a "perfect foresight" model, where production only occurs during the cheapest electricity price hours of the year.
        *   Compares this to a "heuristic" model (based on previous year's prices) and a simple "baseload" model (continuous operation).
        *   The results, including optimal capacity factor, max electricity price, and total LCOH, are printed to the console.

5.  **Plotting:**
    *   `plotting_utils.py`: This script is not directly executed yet but contains helper functions and style definitions used by other scripts (like those in the `notebooks/` directory) to create standardized plots for analysis.

## Other notes

```.gitattributes``` contains settings to filter out large ```.ipynb``` notebooks to store the using ```git-lfs```.