"""Regenerate the committed energy-charts.info reference CSV.

The validation test (tests/test_entsoe_vs_energycharts.py) compares the grid
pipeline against a *pinned* hourly reference so the gate is deterministic and
runs offline. This script is how that reference is produced — run it to refresh
it (e.g. for a new year) or to see exactly where the numbers come from.

    python tests/data/refresh_energycharts_reference.py

Source: https://api.energy-charts.info (Fraunhofer ISE), endpoints /price and
/public_power for the German bidding zone. No API key required.
"""

from pathlib import Path

import pandas as pd
import requests

BASE = "https://api.energy-charts.info"
BZN = "DE-LU"       # /price bidding zone
COUNTRY = "de"      # /public_power country
START, END = "2023-01-01", "2023-12-31"
OUT = Path(__file__).with_name("energycharts_de_2023_hourly.csv")

# energy-charts production-type name -> reference column
PRODUCTION_TYPES = {
    "Solar": "solar",
    "Wind onshore": "wind_onshore",
    "Wind offshore": "wind_offshore",
    "Load": "load",
}


def _index(unix_seconds) -> pd.DatetimeIndex:
    return pd.to_datetime(pd.Series(unix_seconds), unit="s", utc=True).dt.tz_localize(None)


def main() -> None:
    price_json = requests.get(
        f"{BASE}/price", params={"bzn": BZN, "start": START, "end": END}, timeout=180
    ).json()
    power_json = requests.get(
        f"{BASE}/public_power", params={"country": COUNTRY, "start": START, "end": END}, timeout=180
    ).json()

    price = pd.Series(price_json["price"], index=_index(price_json["unix_seconds"]), name="price")

    power_index = _index(power_json["unix_seconds"])
    series = {"price": price}
    by_name = {p["name"]: p["data"] for p in power_json["production_types"]}
    for api_name, column in PRODUCTION_TYPES.items():
        series[column] = pd.Series(by_name[api_name], index=power_index)

    reference = (
        pd.DataFrame(series)
        .resample("1h")
        .mean()
        .loc[f"{START}":f"{END} 23:00"]
        .round(3)
    )
    reference.index.name = "time"
    reference.to_csv(OUT)
    print(f"wrote {OUT} ({len(reference)} hourly rows)")


if __name__ == "__main__":
    main()
