import os
from dotenv import load_dotenv
from pathlib import Path
import entsoe
import pandas as pd

# 1) Load API key
load_dotenv()
api_key = os.getenv("ENTSOE_API_KEY")
if not api_key:
    raise RuntimeError("ENTSOE_API_KEY not found")

# 2) Create raw client
raw = entsoe.EntsoeRawClient(api_key=api_key)

# 3) Define query
area_cc = "DE"
year = 2020
start = pd.Timestamp(year=year, month=1, day=1, hour=0, tz="UTC")
end = pd.Timestamp(year=year + 1, month=1, day=1, hour=0, tz="UTC")

# 4) Query XML
xml = raw.query_installed_generation_capacity(
    country_code=area_cc,
    start=start,
    end=end
)

print("XML length:", len(xml))
print(xml[:1200])

# 5) Save XML to disk
out = Path("data/debug")
out.mkdir(parents=True, exist_ok=True)

xml_path = out / "entsoe_capacity_DE_2020.xml"
xml_path.write_text(xml, encoding="utf-8")

print(f"Saved XML to {xml_path.resolve()}")

