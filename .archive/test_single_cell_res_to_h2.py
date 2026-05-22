"""
Quick single-cell smoke test for the res_to_h2 optimisation pipeline.

Extracts wind CF for one inland northern-Germany grid cell (~52.5N, 10.5E,
near Gifhorn/Lower Saxony) from the existing 2-week Atlite cutout, then runs
the RES:electrolyser sweep with a small step count.

Run from the project root:
    python scripts/test_single_cell_res_to_h2.py

Note: 2-week CF window means LCOH numbers are not representative of annual economics.
"""

import sys
from pathlib import Path

import atlite
import geopandas as gpd
import numpy as np
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).parent / "res_to_h2"))
from res_to_h2_logic import dri_to_el_mw, optimise_res_to_el_ratio

# ---------------------------------------------------------------------------
# Target cell: Gifhorn, Lower Saxony — inland, northern Germany
# ---------------------------------------------------------------------------
LON, LAT = 10.5, 52.5
CELL_HALF = 0.15  # ~0.3-degree box; captures roughly one ERA5 cell

CUTOUT_PATH = Path("data/cutouts/de_2023_jan2w.nc")


def extract_cell_wind_cf(cutout_path: Path, lon: float, lat: float, half: float) -> np.ndarray:
    cutout = atlite.Cutout(str(cutout_path))

    cell_gdf = gpd.GeoDataFrame(
        {"region": ["test_cell"]},
        geometry=[box(lon - half, lat - half, lon + half, lat + half)],
        crs=4326,
    )

    matrix = cutout.indicatormatrix(cell_gdf)

    wind_xr = cutout.wind(
        matrix=matrix,
        turbine="Vestas_V112_3MW",
        capacity_factor=False,
        per_unit=True,
    )

    obj = wind_xr.to_pandas()
    s = obj.iloc[:, 0] if isinstance(obj, type(obj)) and obj.ndim == 2 else obj
    return s.clip(0, 1).to_numpy()


def main():
    print(f"Cutout : {CUTOUT_PATH}")
    print(f"Cell   : {LAT}°N, {LON}°E (Gifhorn, Lower Saxony)\n")

    cf = extract_cell_wind_cf(CUTOUT_PATH, LON, LAT, CELL_HALF)
    print(f"CF array : {len(cf)} timesteps")
    print(f"  mean={cf.mean():.3f}  min={cf.min():.3f}  max={cf.max():.3f}\n")

    el_mw = dri_to_el_mw(
        dri_mt_per_year=2.0,
        h2_intensity_kg_per_t_dri=50.0,
        efficiency_kwh_per_kg=55.0,
        availability_target=1.0,
    )
    print(f"Electrolyser size : {el_mw:.1f} MW\n")

    scenario = {
        "res": {
            "capex_per_mw_eur": 1_000_000,
            "opex_per_mw_per_year_eur": 30_000,
            "lifetime_years": 25,
        },
        "electrolyser": {
            "capex_per_mw_eur": 700_000,
            "opex_per_mw_per_year_eur": 20_000,
            "lifetime_years": 20,
            "efficiency_kwh_per_kg": 55,
        },
        "optimization": {
            "min_res_per_el_mw": 1.0,
            "max_res_per_el_mw": 10.0,
            "steps": 20,
        },
        "storage_proxy": {"enabled": False},
    }

    finance = {"default_wacc": 0.07}

    print("Running sweep (20 steps, storage proxy disabled)...")
    results_df, best_ratio = optimise_res_to_el_ratio(
        res_profile=cf,
        scenario=scenario,
        finance=finance,
        timestep_hours=1.0,
        electrolyser_mw=el_mw,
    )

    best = results_df.loc[best_ratio]
    print("\n=== Results ===")
    print(f"Best RES:EL ratio  : {best_ratio:.2f}")
    print(f"LCOH               : {best['lcoh_eur_per_kg']:.3f} €/kg")
    print(f"EL utilisation     : {best['el_utilisation']*100:.1f}%")
    print(f"Curtailment rate   : {best['curtailment_rate']*100:.1f}%")
    print(f"Firm reliability   : {best['firm_reliability']*100:.1f}%")


if __name__ == "__main__":
    main()
