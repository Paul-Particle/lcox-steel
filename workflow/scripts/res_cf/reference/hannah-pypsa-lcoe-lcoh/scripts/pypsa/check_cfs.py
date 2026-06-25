import pandas as pd

files_to_check = [
    "data/res_cf/bra_cf_2023_bestsite_p95_anchor-wind_onshore.csv",
    "data/res_cf/de_cf_2023_bestsite_p95_anchor-solar.csv",
    "data/res_cf/es_cf_2023_bestsite_p95_anchor-solar.csv",
    "data/res_cf/fr_cf_2023_bestsite_p95_anchor-solar.csv",
]

for path in files_to_check:
    df = pd.read_csv(path)
    cc = path.split("/")[-1].split("_")[0].upper()
    anchor = path.split("anchor-")[1].replace(".csv", "")
    print(f"\n{cc} | {anchor}")
    print(f"  wind_onshore:  {df['wind_onshore_cf'].mean():.3f}")
    print(f"  wind_offshore: {df['wind_offshore_cf'].mean():.3f}")
    print(f"  solar:         {df['solar_cf'].mean():.3f}")