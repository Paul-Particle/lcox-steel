import os
import sys
from pathlib import Path

import pandas as pd
from entsoe import EntsoePandasClient

try:
    from dotenv import load_dotenv
except ImportError:
    load_dotenv = None


# --- config paths ---------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
PROCESSED_PATH = ROOT / "data" / "processed_data.feather"
OUT_PATH = ROOT / "data" / "installed_capacities_2024.csv"


# --- ENTSO-E client ------------------------------------------------------

def get_entsoe_client() -> EntsoePandasClient:
    if load_dotenv is not None:
        # load .env in project root if present
        env_path = ROOT / ".env"
        if env_path.exists():
            load_dotenv(env_path)

    api_key = os.getenv("ENTSOE_API_KEY")
    if not api_key:
        raise RuntimeError(
            "ENTSOE_API_KEY not set. Put it in .env or as an environment variable."
        )
    return EntsoePandasClient(api_key=api_key)


# --- main logic ----------------------------------------------------------

GEN_NAMES = {
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


def main():
    if not PROCESSED_PATH.exists():
        print(f"Cannot find {PROCESSED_PATH}. Run Snakemake first.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading areas from {PROCESSED_PATH}...")
    df = pd.read_feather(PROCESSED_PATH)
    areas = ["DE_LU"] #For Debugging 
    #areas = df.columns.get_level_values(0).unique()

    client = get_entsoe_client()

    #start = pd.Timestamp(year=2024, month=1, day=1, hour=0, tz="UTC")
    #end = pd.Timestamp(year=2025, month=1, day=1, hour=0, tz="UTC")
    # (not used; we query capacities year-by-year below)

    records = []

    for area in areas:
        area_cc = area
        if area_cc.startswith("SE"):
            area_cc = "SE"

        print(f"Querying installed capacity for {area_cc}...")

        dfTs = []

        for y in range(2015, 2026):
            y_start = pd.Timestamp(year=y, month=1, day=1, hour=0, tz="UTC")
            y_end = pd.Timestamp(year=y + 1, month=1, day=1, hour=0, tz="UTC")

            try:
                df_year = client.query_installed_generation_capacity(
                    country_code=area_cc, start=y_start, end=y_end
                )
            except Exception:
                continue

    # keep only if at least one real number exists
            if df_year is not None and df_year.notna().any().any():
                dfTs.append(df_year)

        if not dfTs:
            print(f"  No capacity data found for {area_cc}", file=sys.stderr)
            continue

        dfT = pd.concat(dfTs).sort_index()

# Use a recent window to avoid historical weirdness; fall back if empty
        df_recent = dfT.loc[pd.Timestamp("2023-01-01", tz="UTC") :]
        if df_recent.empty:
            df_recent = dfT

        ser_mean = df_recent.mean(numeric_only=True)
        ser_mean = ser_mean.rename(GEN_NAMES)
        ser_mean = ser_mean[ser_mean.index.isin(GEN_NAMES.values())]


        row = ser_mean.to_dict()
        row["area"] = area  # original MultiIndex area name
        records.append(row)

    if not records:
        print("No capacity records collected. Aborting.", file=sys.stderr)
        sys.exit(1)

    df_out = pd.DataFrame.from_records(records).set_index("area").fillna(0.0)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(OUT_PATH)

    print(f"Saved installed capacities to {OUT_PATH}")


if __name__ == "__main__":
    main()
