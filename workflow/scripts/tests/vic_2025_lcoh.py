"""
End-to-end LCOH test for Victoria, Australia (NEM zone VIC1), 2025 weather + prices.

Picks 4 cells inside an approximate VIC bounding box (wind P95, solar P95, plus
2 randomly chosen cells with a fixed seed), runs the PyPSA model for each cell
under 3 scenarios (grid_only, captive, hybrid) — 12 optimisations total — and
emits:
  results/vic_2025/summary.csv     LCOH + sizing + VIC-wide CF averages per scenario
  results/vic_2025/capacities.csv  per-(location, scenario) wind/solar capacities + max accepted grid price
  results/vic_2025/locations.csv   chosen cells with lat/lon
  results/vic_2025/cf_grids.nc     per-cell annual-mean CF grids + VIC mask
  results/vic_2025/price_series.csv VIC1 hourly EUR price series used (for the price-distribution viz)
  results/vic_2025/{loc}_{scen}.nc 12 solved PyPSA networks

Currency: NEM prices are AUD-native. Converted to EUR inline using the
top-level `fx.eur_per_aud` rate in config/projects.yaml. Bump that if you
need real accuracy.
"""

from __future__ import annotations

# pandas 3.0 defaults to ArrowStringArray for strings; xarray (used by PyPSA)
# doesn't support it. Set python-native strings before any import that touches
# pandas string data.
import pandas as pd
pd.options.mode.string_storage = "python"

import importlib.util
import sys
from pathlib import Path

import atlite
import geopandas as gpd
import numpy as np
import xarray as xr
import yaml
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common._paths import REPO_ROOT, CUTOUTS, RESULTS, SHAPES_RES  # noqa: E402

SEED = 42
NEM_REGION = "VIC1"


def load_bestsite_module():
    """Load make_bestsite_cf.py via importlib (same pattern as complementarity.py)."""
    spec = importlib.util.spec_from_file_location(
        "bestsite", REPO_ROOT / "scripts/res_cf/make_bestsite_cf.py"
    )
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def build_vic_geom(cfg: dict):
    """AUS land polygon ∩ vic_bbox → VIC-approximating geometry."""
    aus = cfg["res_cf"]["countries"]["aus"]
    lon_min, lon_max, lat_min, lat_max = aus["vic_bbox"]
    regions = gpd.read_parquet(SHAPES_RES / "regions.parquet")
    aus_geom = regions.loc[regions["region"] == "AUS", "geometry"].iloc[0]
    return aus_geom.intersection(box(lon_min, lat_min, lon_max, lat_max))


def pick_cells(wind_cf_year: xr.DataArray, solar_cf_year: xr.DataArray, vic_geom,
               find_p95_cell, mask_cells_inside) -> list[tuple[str, int, int]]:
    """Return 4 (name, y_idx, x_idx) tuples."""
    wind_p95  = find_p95_cell(wind_cf_year, vic_geom)
    solar_p95 = find_p95_cell(solar_cf_year, vic_geom)

    wind_mean  = wind_cf_year.mean("time").values
    solar_mean = solar_cf_year.mean("time").values
    in_vic     = mask_cells_inside(wind_cf_year.mean("time"), vic_geom)

    valid = in_vic & np.isfinite(wind_mean) & np.isfinite(solar_mean)
    ys, xs = np.where(valid)
    taken = {wind_p95, solar_p95}
    candidates = [(int(y), int(x)) for y, x in zip(ys, xs) if (int(y), int(x)) not in taken]
    if len(candidates) < 2:
        raise RuntimeError(f"Only {len(candidates)} candidate VIC cells available; cannot pick 2.")

    rng = np.random.default_rng(SEED)
    picks = rng.choice(len(candidates), size=2, replace=False)
    r1, r2 = candidates[picks[0]], candidates[picks[1]]

    return [
        ("wind_p95",  wind_p95[0],  wind_p95[1]),
        ("solar_p95", solar_p95[0], solar_p95[1]),
        ("random_1",  r1[0], r1[1]),
        ("random_2",  r2[0], r2[1]),
    ]


def vic_area_weighted_mean(cf_mean: np.ndarray, vic_mask: np.ndarray, lat_1d: np.ndarray) -> float:
    """cos(lat) area-weighted mean of cf_mean over cells where vic_mask is True."""
    lat_2d = np.broadcast_to(lat_1d[:, None], cf_mean.shape)
    w = np.cos(np.deg2rad(lat_2d)) * vic_mask
    return float((cf_mean * w).sum() / w.sum())


def main() -> None:
    cfg         = yaml.safe_load(open(REPO_ROOT / "config/config.yaml"))
    assumptions = yaml.safe_load(open(REPO_ROOT / "config/assumptions.yaml"))
    projects    = yaml.safe_load(open(REPO_ROOT / "config/projects.yaml"))
    year        = int(cfg["res_cf"]["year"])
    eur_per_aud = float(projects["fx"]["eur_per_aud"])

    bestsite = load_bestsite_module()
    # Override standalone defaults in case they were picked up before our config was applied
    bestsite.YEAR = year
    bestsite.COUNTRIES = ["aus"]

    # h2_dri imports (after pandas string_storage shim above)
    sys.path.insert(0, str(REPO_ROOT / "scripts/h2_dri"))
    from network import build_network          # noqa: E402
    from costs import extract_summary          # noqa: E402

    # ── VIC geometry & per-cell CF grids ──
    vic_geom = build_vic_geom(cfg)
    print("Building wind CF grid for AUS 2025 (atlite, 4 quarters)…")
    wind_cf_year  = bestsite.build_cf_year("AUS", "wind_onshore")
    print("Building solar CF grid for AUS 2025 (atlite, 4 quarters)…")
    solar_cf_year = bestsite.build_cf_year("AUS", "solar")

    # ── Cell selection ──
    locations = pick_cells(
        wind_cf_year, solar_cf_year, vic_geom,
        bestsite._find_p95_cell, bestsite.mask_cells_inside,
    )
    lon_1d = wind_cf_year.x.values
    lat_1d = wind_cf_year.y.values
    print("Locations:")
    for name, y, x in locations:
        print(f"  {name:9s}: y={y} x={x} lon={lon_1d[x]:.2f} lat={lat_1d[y]:.2f}")

    # ── Per-cell hourly CF series ──
    cell_cfs = {}
    for name, y, x in locations:
        wind_ts  = bestsite.extract_cell_timeseries(wind_cf_year,  y, x, "wind_onshore")
        solar_ts = bestsite.extract_cell_timeseries(solar_cf_year, y, x, "solar")
        # extract_cell_timeseries returns a pd.Series; align indices on the common one
        df = pd.DataFrame({"wind_onshore": wind_ts, "solar": solar_ts})
        cell_cfs[name] = df

    # ── Price series ──
    nem        = pd.read_parquet(REPO_ROOT / "resources/nem_processed.parquet")
    price_aud  = nem[(NEM_REGION, "price")]
    price_aud  = price_aud[price_aud.index.year == year]
    price_eur  = (price_aud * eur_per_aud).reindex(cell_cfs["wind_p95"].index)
    if price_eur.isna().any():
        n_nan = int(price_eur.isna().sum())
        print(f"WARNING: filling {n_nan} NaN prices (post-reindex) with the median.")
        price_eur = price_eur.fillna(price_eur.median())

    # ── 12 PyPSA runs ──
    scenarios = {
        "grid_only": dict(techs=[], grid_connected=True),
        "captive":   dict(techs=["wind_onshore", "solar"], grid_connected=False),
        "hybrid":    dict(techs=["wind_onshore", "solar"], grid_connected=True),
    }
    plant_cfg = projects["projects"][0]["plant"]  # reuse DE_2023_baseline plant spec
    project_cfg = {"plant": plant_cfg, "country": "AUS", "year": year}

    out_dir = RESULTS / "vic_2025"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for loc_name, y_idx, x_idx in locations:
        for sc_name, sc in scenarios.items():
            cf_ts = (
                cell_cfs[loc_name][sc["techs"]].copy()
                if sc["techs"]
                else pd.DataFrame(index=cell_cfs[loc_name].index)
            )
            price = price_eur if sc["grid_connected"] else None

            print(f"Solving {loc_name} / {sc_name}…")
            n = build_network(project_cfg, assumptions, cf_ts, price_series=price)
            n.optimize(solver_name="highs")

            row = extract_summary(n, "vic_2025", f"{loc_name}_{sc_name}",
                                  assumptions["h2"]["lhv_kwh_per_kg"])
            row["location"] = loc_name
            row["scenario"] = sc_name

            if sc["grid_connected"] and "grid_import" in n.generators_t.p.columns:
                grid_p   = n.generators_t.p["grid_import"]
                accepted = price_eur.loc[grid_p > 1e-3]
                row["max_grid_import_price_eur_mwh"] = float(accepted.max()) if len(accepted) else float("nan")
            else:
                row["max_grid_import_price_eur_mwh"] = float("nan")

            n.export_to_netcdf(out_dir / f"{loc_name}_{sc_name}.nc")
            rows.append(row)

    summary = pd.DataFrame(rows)

    # ── VIC-wide CF averages (area-weighted by cos(lat)) ──
    wind_mean  = wind_cf_year.mean("time").values
    solar_mean = solar_cf_year.mean("time").values
    vic_mask   = bestsite.mask_cells_inside(wind_cf_year.mean("time"), vic_geom)

    wind_cf_vic_avg  = vic_area_weighted_mean(wind_mean,  vic_mask, lat_1d)
    solar_cf_vic_avg = vic_area_weighted_mean(solar_mean, vic_mask, lat_1d)
    print(f"\nVIC-wide CF (cos-lat weighted): wind {wind_cf_vic_avg:.3f}, solar {solar_cf_vic_avg:.3f}")

    def mix(row):
        w = row.get("wind_onshore_mw_opt") or 0.0
        s = row.get("solar_mw_opt") or 0.0
        if w + s <= 0:
            return float("nan")
        return (w * wind_cf_vic_avg + s * solar_cf_vic_avg) / (w + s)

    summary["mix_cf_vic_avg"]   = summary.apply(mix, axis=1)
    summary["wind_cf_vic_avg"]  = wind_cf_vic_avg
    summary["solar_cf_vic_avg"] = solar_cf_vic_avg

    # ── Outputs ──
    summary_cols = [
        "location", "scenario",
        "lcoh_eur_per_kg",
        "electrolyser_mw",
        "wind_onshore_mw_opt", "solar_mw_opt",
        "battery_mw_opt", "battery_mwh_opt",
        "h2_buffer_mwh_lhv_opt",
        "max_grid_import_price_eur_mwh",
        "mix_cf_vic_avg", "wind_cf_vic_avg", "solar_cf_vic_avg",
    ]
    for c in summary_cols:
        if c not in summary.columns:
            summary[c] = float("nan")
    summary = summary[summary_cols]
    summary.to_csv(out_dir / "summary.csv", index=False)
    print(f"\nWrote {out_dir / 'summary.csv'}")
    print(summary[["location", "scenario", "lcoh_eur_per_kg",
                   "wind_onshore_mw_opt", "solar_mw_opt"]].to_string(index=False))

    # capacities + locations + price series + cf grids for viz
    summary[["location", "scenario",
             "wind_onshore_mw_opt", "solar_mw_opt",
             "max_grid_import_price_eur_mwh"]].to_csv(out_dir / "capacities.csv", index=False)

    pd.DataFrame([
        {"location": name, "y_idx": int(y), "x_idx": int(x),
         "lon": float(lon_1d[x]), "lat": float(lat_1d[y])}
        for name, y, x in locations
    ]).to_csv(out_dir / "locations.csv", index=False)

    price_eur.rename("price_eur_per_mwh").to_csv(out_dir / "price_series.csv")

    xr.Dataset(
        {
            "wind_cf_annual_mean":  (("y", "x"), wind_mean),
            "solar_cf_annual_mean": (("y", "x"), solar_mean),
            "in_vic":               (("y", "x"), vic_mask),
        },
        coords={"y": lat_1d, "x": lon_1d},
    ).to_netcdf(out_dir / "cf_grids.nc")

    print(f"Wrote {out_dir / 'cf_grids.nc'}, capacities.csv, locations.csv, price_series.csv")


if __name__ == "__main__":
    main()
