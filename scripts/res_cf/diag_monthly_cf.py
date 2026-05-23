import pandas as pd

def monthly_means(path):
    df = pd.read_csv(path)
    df["time"] = pd.to_datetime(df["time"])
    s = df.set_index("time")["cf"]
    m = s.resample("MS").mean()
    print("\n===", path, "(monthly mean CF) ===")
    print(m.to_string())

monthly_means(str(RES_CF / "de_wind_onshore_cf_2023.csv"))
monthly_means(str(RES_CF / "de_solar_cf_2023.csv"))
