import pandas as pd
from common._paths import RES_CF

def monthly_means(path):
    s = pd.read_parquet(path)["cf"]
    m = s.resample("MS").mean()
    print("\n===", path, "(monthly mean CF) ===")
    print(m.to_string())

monthly_means(str(RES_CF / "de_wind_onshore_cf_2023.parquet"))
monthly_means(str(RES_CF / "de_solar_cf_2023.parquet"))
