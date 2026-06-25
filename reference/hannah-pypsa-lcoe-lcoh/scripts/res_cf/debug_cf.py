import pandas as pd
from pathlib import Path

files = sorted(Path("data/res_cf").glob("*_cf_2023_bestsite_p95.csv"))

for f in files:
    df = pd.read_csv(f)

    print("\n" + "="*60)
    print(f.name)

    for col in ["wind_onshore_cf", "wind_offshore_cf"]:
        s = df[col]

        mean = s.mean()
        maxv = s.max()
        p95 = s.quantile(0.95)
        share_high = (s > 0.95).mean()
        share_sat = (s > 0.99).mean()

        print(f"\n{col}:")
        print(f"  mean        = {mean:.3f}")
        print(f"  max         = {maxv:.3f}")
        print(f"  p95         = {p95:.3f}")
        print(f"  >0.95 share = {share_high:.3%}")
        print(f"  >0.99 share = {share_sat:.3%}")