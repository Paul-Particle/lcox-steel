import os
from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd
import entsoe
from dotenv import load_dotenv


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data"
OUT_PATH = OUT_DIR / "installed_capacities_latest.csv"

PROCESSED_PATH = ROOT / "data" / "processed_data.feather"

def infer_start_year(default=2023) -> int:
    if PROCESSED_PATH.exists():
        gen = pd.read_feather(PROCESSED_PATH)
        return int(gen.index.max().year)
    return default

AREAS = ["DE_LU", "FR", "ES"]

START_YEAR = infer_start_year(default=2023)
print(f"[fetch_entsoe_capacities] Using START_YEAR={START_YEAR} (inferred from {PROCESSED_PATH.name} if present)")

MIN_YEAR = 2018

# Map for this endpoint only (country_code expected)
AREA_TO_COUNTRY_CODE = {
    "DE_LU": "DE",
    "FR": "FR",
    "ES": "ES",
}

PSR_CODE_TO_NAME = {
    "B01": "Biomass",
    "B02": "Fossil Brown coal/Lignite",
    "B03": "Fossil Coal-derived gas",
    "B04": "Fossil Gas",
    "B05": "Fossil Hard coal",
    "B06": "Fossil Oil",
    "B07": "Fossil Oil shale",
    "B08": "Fossil Peat",
    "B09": "Geothermal",
    "B10": "Hydro Pumped Storage",
    "B11": "Hydro Run-of-river and poundage",
    "B12": "Hydro Water Reservoir",
    "B13": "Marine",
    "B14": "Nuclear",
    "B15": "Other",
    "B16": "Solar",
    "B17": "Waste",
    "B18": "Wind Offshore",
    "B19": "Wind Onshore",
    "B20": "Other renewable",
}

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


def get_api_key() -> str:
    load_dotenv(ROOT / ".env")
    api_key = os.getenv("ENTSOE_API_KEY")
    if not api_key:
        raise RuntimeError("ENTSOE_API_KEY not set. Put it in .env or environment.")
    return api_key


def _detect_ns(root: ET.Element) -> dict:
    """Extract namespace from root tag like '{urn:...}GL_MarketDocument'."""
    if root.tag.startswith("{") and "}" in root.tag:
        uri = root.tag.split("}")[0].strip("{")
        return {"ns": uri}
    # fallback: no namespace
    return {"ns": ""}


def parse_installed_capacity_xml(xml: str) -> dict:
    """
    Return dict {psr_code: quantity_float}.
    If multiple Points exist for a TimeSeries, take max(quantity).
    """
    root = ET.fromstring(xml)
    NS = _detect_ns(root)

    # If namespace is empty, ElementTree expects plain tags
    def q(tag: str) -> str:
        return f"ns:{tag}" if NS["ns"] else tag

    out = {}

    ts_path = q("TimeSeries")
    for ts in root.findall(ts_path, NS if NS["ns"] else None):
        psr = ts.findtext(
            f"{q('MktPSRType')}/{q('psrType')}",
            default=None,
            namespaces=NS if NS["ns"] else None,
        )
        if psr is None:
            continue

        # Collect all Point quantities
        quantities = []
        for pt in ts.findall(
            f"{q('Period')}/{q('Point')}",
            NS if NS["ns"] else None,
        ):
            qty_txt = pt.findtext(
                q("quantity"),
                default=None,
                namespaces=NS if NS["ns"] else None,
            )
            if qty_txt is None:
                continue
            try:
                quantities.append(float(qty_txt))
            except ValueError:
                continue

        if not quantities:
            continue

        out[psr] = max(quantities)

    return out


def fetch_capacity_for_year(raw_client, area: str, year: int) -> dict:
    """
    Query raw XML for a given area/year, parse it, return dict {internal_tech: capacity_mw}.
    """
    start = pd.Timestamp(year=year, month=1, day=1, hour=0, tz="UTC")
    end = pd.Timestamp(year=year + 1, month=1, day=1, hour=0, tz="UTC")

    country_code = AREA_TO_COUNTRY_CODE.get(area, area)

    xml = raw_client.query_installed_generation_capacity(
        country_code=country_code, start=start, end=end
    )

    psr_to_qty = parse_installed_capacity_xml(xml)

    out = {}
    for psr_code, qty in psr_to_qty.items():
        entsoe_name = PSR_CODE_TO_NAME.get(psr_code)
        if not entsoe_name:
            continue
        internal = GEN_NAMES.get(entsoe_name)
        if not internal:
            continue
        out[internal] = out.get(internal, 0.0) + qty

    return out


def most_recent_capacity(raw_client, area: str) -> tuple[int | None, dict]:
    for year in range(START_YEAR, MIN_YEAR - 1, -1):
        try:
            caps = fetch_capacity_for_year(raw_client, area, year)
        except Exception as e:
            print(f"  Warning: {area} {year} fetch failed: {type(e).__name__}: {e}")
            continue

        # Consider "non-empty" as having any positive values
        if any(v > 0 for v in caps.values()):
            return year, caps

    return None, {}


def main():
    api_key = get_api_key()
    raw_client = entsoe.EntsoeRawClient(api_key=api_key)

    records = []
    for area in AREAS:
        print(f"Fetching installed capacity for {area} (trying {START_YEAR}↓{MIN_YEAR}) ...")
        year, caps = most_recent_capacity(raw_client, area)

        if year is None:
            print(f"  No non-zero capacity found for {area} in {START_YEAR}..{MIN_YEAR}")
            record = {"area": area, "year_used": None}
        else:
            print(f"  Using year {year} for {area}")
            record = {"area": area, "year_used": year}
            record.update(caps)

        records.append(record)

    df_out = pd.DataFrame(records).fillna(0.0)

    base_cols = ["area", "year_used"]
    tech_cols = sorted([c for c in df_out.columns if c not in base_cols])
    df_out = df_out[base_cols + tech_cols]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(OUT_PATH, index=False)

    print(f"\nSaved capacities to: {OUT_PATH}")


if __name__ == "__main__":
    main()