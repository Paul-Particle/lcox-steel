"""Cut Alberta and Ontario from Natural Earth Admin 1 into region parquets.

The same hand-cut route the NEM zone geometries took: make_area_geometry reads
Admin 0 and refuses for an area that is someone's zone, so a province cannot
come from the workflow until issue #64 lands. Natural Earth 10m Admin 1 is a
manual download into data/shapes/ne_10m_admin_1_states_provinces/.
"""

from pathlib import Path
import geopandas as gpd

REPO_ROOT = Path(__file__).resolve().parents[1]
NATURAL_EARTH_ZIP = (
    REPO_ROOT
    / "data/shapes/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.zip"
)
SHAPES_OUT = REPO_ROOT / "resources/shapes"

# Statistics Canada total area including water — the check that the cut is sane.
OFFICIAL_AREA_KM2 = {"AB": 661_848, "ON": 1_076_395}

provinces = gpd.read_file(str(NATURAL_EARTH_ZIP)).to_crs(4326)
canadian_provinces = provinces.loc[provinces["adm0_a3"] == "CAN"]

for area, iso_3166_2, region in [("AB", "CA-AB", "AB"), ("ONT", "CA-ON", "ON")]:
    selected = canadian_provinces.loc[canadian_provinces["iso_3166_2"] == iso_3166_2]
    province_geometry = gpd.GeoDataFrame(
        {"region": [region], "geometry": [selected.geometry.union_all()]}, crs=4326
    )
    province_geometry["geometry"] = province_geometry["geometry"].buffer(0)

    area_km2 = province_geometry.to_crs(6933).area.sum() / 1e6
    deviation_pct = 100 * (area_km2 / OFFICIAL_AREA_KM2[region] - 1)
    print(
        f"{area}: {area_km2:,.0f} km2 vs official {OFFICIAL_AREA_KM2[region]:,} "
        f"({deviation_pct:+.1f}%)"
    )

    province_geometry.to_parquet(SHAPES_OUT / f"{area}_geo.parquet")
