"""Cut the four ONS submarkets from Natural Earth Admin 1 and give them cutouts.

The submarkets are groups of states, so the land geometry is a dissolve over the
states each one contains. The mapping is electrical, not geographic, and three
entries are worth knowing about: Acre and Rondonia sit in Sudeste/Centro-Oeste
rather than Norte, which is where the Madeira HVDC links put them; Maranhao sits
in Norte rather than Nordeste; and Roraima joined the SIN only during 2025.

Each submarket takes the whole-Brazil cutout rather than a slice of it. A
submarket's cells are a strict subset of the national grid — same ERA5 lattice,
same year — and the CF scripts select cells by the region geometry, not by the
cutout's extent, so the result is identical either way. Slicing would only save
compute, and it does not: the Sudeste box alone spans Brazil from Acre to
Espirito Santo, so it is half the national file, while a whole-Brazil CF build
takes about six minutes. Writing four multi-gigabyte slices to save four minutes
each is the wrong trade, so the cutout is hardlinked and the geometry does the
work.

Run from the repo root with the project env.
"""
import os
from pathlib import Path

import geopandas as gpd

REPO_ROOT = Path(__file__).resolve().parents[1]

NATURAL_EARTH_ZIP = (REPO_ROOT / "data/shapes/ne_10m_admin_1_states_provinces"
                     / "ne_10m_admin_1_states_provinces.zip")
EEZ_ZIP = REPO_ROOT / "data/shapes/offshore_zones/eez_v12.zip"
SHAPES_OUT = REPO_ROOT / "resources/shapes"
CUTOUTS = REPO_ROOT / "cutouts"
BRAZIL_CUTOUT = CUTOUTS / "BRA_20250101_20251231.nc"
START_DATE, END_DATE = "20250101", "20251231"
OFFSHORE_MAX_KM = 200.0

# ONS submarket -> the states it covers, as ISO 3166-2 subdivision codes.
SUBMARKET_STATES = {
    "BR_SE": ["SP", "RJ", "MG", "ES", "GO", "DF", "MT", "MS", "AC", "RO"],
    "BR_S":  ["PR", "SC", "RS"],
    "BR_NE": ["BA", "SE", "AL", "PE", "PB", "RN", "CE", "PI"],
    "BR_N":  ["AM", "PA", "AP", "TO", "MA", "RR"],
}

# IBGE 2022 state areas, km2 — summed per submarket to check the dissolve. They
# total 8,512,692 against Brazil's 8,510,417, so the four lists together cover
# the country once and only once.
IBGE_STATE_AREA_KM2 = {
    "SP": 248_219, "RJ": 43_750, "MG": 586_513, "ES": 46_074, "GO": 340_111,
    "DF": 5_760, "MT": 903_207, "MS": 357_145, "AC": 164_123, "RO": 237_765,
    "PR": 199_307, "SC": 95_730, "RS": 281_730,
    "BA": 564_733, "SE": 21_925, "AL": 27_848, "PE": 98_148, "PB": 56_585,
    "RN": 52_811, "CE": 148_921, "PI": 251_577,
    "AM": 1_559_168, "PA": 1_245_870, "AP": 142_470, "TO": 277_621,
    "MA": 331_937, "RR": 223_644,
}

states = gpd.read_file(str(NATURAL_EARTH_ZIP)).to_crs(4326)
brazilian_states = states.loc[states["adm0_a3"] == "BRA"]
eez = gpd.read_file(str(EEZ_ZIP)).to_crs(4326)
brazilian_eez = eez.loc[(eez["ISO_TER1"] == "BRA") & (eez["POL_TYPE"] == "200NM")].union_all()

for area, state_codes in SUBMARKET_STATES.items():
    wanted = [f"BR-{code}" for code in state_codes]
    selected = brazilian_states.loc[brazilian_states["iso_3166_2"].isin(wanted)]
    missing = sorted(set(wanted) - set(selected["iso_3166_2"]))
    if missing:
        raise ValueError(f"{area}: Natural Earth has no state for {missing}")

    land = gpd.GeoDataFrame(
        {"region": [area], "geometry": [selected.geometry.union_all()]}, crs=4326
    )
    land["geometry"] = land["geometry"].buffer(0)
    land_km2 = land.to_crs(6933)["geometry"].area.iloc[0] / 1e6
    official_km2 = sum(IBGE_STATE_AREA_KM2[code] for code in state_codes)
    land.to_parquet(SHAPES_OUT / f"{area}_geo.parquet")

    # Same construction as b_make_offshore_geometry: the national EEZ minus this
    # submarket's land, kept within OFFSHORE_MAX_KM of its own coast — which is
    # what confines a national EEZ to one submarket's stretch of it.
    land_buffer = land.to_crs(6933).copy()
    land_buffer["geometry"] = land_buffer.buffer(OFFSHORE_MAX_KM * 1000)
    buffer_geometry = land_buffer.to_crs(4326)["geometry"].iloc[0]
    offshore_geometry = (
        brazilian_eez.difference(land["geometry"].iloc[0]).intersection(buffer_geometry)
    )
    offshore = gpd.GeoDataFrame({"region": [area], "geometry": [offshore_geometry]}, crs=4326)
    offshore["geometry"] = offshore["geometry"].buffer(0)
    offshore.to_parquet(SHAPES_OUT / f"{area}_offshore_geo.parquet")
    offshore_km2 = offshore.to_crs(6933)["geometry"].area.iloc[0] / 1e6

    cutout_path = CUTOUTS / f"{area}_{START_DATE}_{END_DATE}.nc"
    if not cutout_path.exists():
        os.link(BRAZIL_CUTOUT, cutout_path)

    print(f"{area:6s} land {land_km2:>10,.0f} km2 vs IBGE {official_km2:>10,.0f} "
          f"({100 * (land_km2 / official_km2 - 1):+.1f}%)  "
          f"offshore {offshore_km2:>9,.0f} km2  cutout -> {BRAZIL_CUTOUT.name}")
