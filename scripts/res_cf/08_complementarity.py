"""
08_complementarity.py

Purpose
-------
Identifies the best combination of RES cells (onshore wind, offshore wind, solar)
for each country by screening for complementarity — i.e. combinations whose hourly
profiles interlock well to maximise electrolyser utilisation.

This script produces two outputs per country:
  (a) Top N complementary cell triplets for the best-site RES mix scenario
  (b) National mean profiles for the country-average scenario

Both outputs feed directly into the PyPSA optimisation (scripts/h2_dri/run.py).

Method
------
For each country:
  1. Load per-cell CF grids from Atlite cutouts via build_cf_year() (reused from 07)
  2. Pre-filter cells to quality_floor percentile within correct geometry
  3. Pre-compute actual lon/lat coordinates for all candidate cells
  4. Screen all valid triplets (onshore_i, offshore_j, solar_k):
     - spatial filter: max pairwise distance <= max_radius_km
     - complementarity score:
         coincidence  = fraction of hours where combined CF > threshold
         mean_corr    = mean of 3 pairwise Pearson correlations
         score        = w1 * coincidence - w2 * mean_corr
  5. Rank triplets by score, save top N
  6. Save national mean profiles for average scenario

Outputs
-------
data/res_cf/<cc>_complementarity_top{N}_2023.csv
data/res_cf/<cc>_average_profiles_2023.csv
"""

from __future__ import annotations

import importlib.util
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# ── Import reusable functions from script 07 ─────────────────────────────────
_spec = importlib.util.spec_from_file_location(
    "bestsite",
    Path(__file__).parent / "07_make_bestsite_cf_timeseries.py"
)
_bestsite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bestsite)

build_cf_year           = _bestsite.build_cf_year
extract_cell_timeseries = _bestsite.extract_cell_timeseries
geometry_for_tech       = _bestsite.geometry_for_tech
haversine_distance_km   = _bestsite.haversine_distance_km
mask_cells_inside       = _bestsite.mask_cells_inside

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent.parent
CONFIG_PATH  = PROJECT_ROOT / "config_hannah.yaml"
OUTDIR       = PROJECT_ROOT / "data" / "res_cf"
CF_DIR       = PROJECT_ROOT / "data" / "res_cf"

YEAR  = 2023
TECHS = ["wind_onshore", "wind_offshore", "solar"]

MAX_TRIPLETS_BRUTE_FORCE = 500_000

# Module-level country variable set in main loop
_current_country: str = ""


# ── Config ────────────────────────────────────────────────────────────────────

def load_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_pypsa_config(config: dict) -> dict:
    pypsa = config.get("defaults", {}).get("pypsa", {})
    spatial = config.get("defaults", {}).get("spatial_matching_res_mix", {})
    return {
        "countries":             pypsa.get("countries", ["DE", "FR", "ES", "AUS", "BRA"]),
        "year":                  pypsa.get("year", YEAR),
        "top_n":                 pypsa.get("top_n_combinations", 10),
        "coincidence_threshold": pypsa.get("coincidence_threshold", 0.20),
        "w_coincidence":         pypsa.get("complementarity_weights", {}).get("w_coincidence", 0.6),
        "w_correlation":         pypsa.get("complementarity_weights", {}).get("w_correlation", 0.4),
        "max_radius_km":         float(spatial.get("max_radius_km", 300.0)),
        "quality_floor":         float(spatial.get("quality_floor_fraction", 0.90)),
    }


# ── Coordinate helpers ────────────────────────────────────────────────────────

def get_candidate_coords(cf, ys: np.ndarray, xs: np.ndarray):
    """Return (lons, lats) arrays for all candidate cells."""
    lons = np.array([float(cf.x.values[xi]) for xi in xs])
    lats = np.array([float(cf.y.values[yi]) for yi in ys])
    return lons, lats


def max_pairwise_distance_km_coords(
    lon_on:  float, lat_on:  float,
    lon_off: float, lat_off: float,
    lon_sol: float, lat_sol: float,
) -> float:
    """Maximum great-circle distance between any two cells in the triplet."""
    d1 = haversine_distance_km(lon_on,  lat_on,
                                np.array([lon_off]), np.array([lat_off]))[0]
    d2 = haversine_distance_km(lon_on,  lat_on,
                                np.array([lon_sol]), np.array([lat_sol]))[0]
    d3 = haversine_distance_km(lon_off, lat_off,
                                np.array([lon_sol]), np.array([lat_sol]))[0]
    return float(max(d1, d2, d3))


# ── Complementarity scoring ───────────────────────────────────────────────────

def score_triplet(
    ts_on:  np.ndarray,
    ts_off: np.ndarray,
    ts_sol: np.ndarray,
    threshold:    float,
    w_coincidence: float,
    w_correlation: float,
) -> tuple[float, float, float]:
    """
    Returns (score, coincidence, mean_pairwise_corr).
    Higher score = better complementarity.
    """
    combined    = (ts_on + ts_off + ts_sol) / 3.0
    coincidence = float(np.mean(combined > threshold))

    corr_on_off  = float(np.corrcoef(ts_on,  ts_off)[0, 1])
    corr_on_sol  = float(np.corrcoef(ts_on,  ts_sol)[0, 1])
    corr_off_sol = float(np.corrcoef(ts_off, ts_sol)[0, 1])
    mean_corr    = (corr_on_off + corr_on_sol + corr_off_sol) / 3.0

    score = w_coincidence * coincidence - w_correlation * mean_corr
    return score, coincidence, mean_corr


# ── Brute-force screen ────────────────────────────────────────────────────────

def brute_force_screen(
    ys_on,  xs_on,  lons_on,  lats_on,
    ys_off, xs_off, lons_off, lats_off,
    ys_sol, xs_sol, lons_sol, lats_sol,
    ts_on_all:  np.ndarray,
    ts_off_all: np.ndarray,
    ts_sol_all: np.ndarray,
    cf_on, cf_off, cf_sol,
    max_radius_km:  float,
    threshold:      float,
    w_coincidence:  float,
    w_correlation:  float,
    top_n:          int,
) -> list[dict]:
    records = []

    for i, j, k in product(range(len(ys_on)),
                             range(len(ys_off)),
                             range(len(ys_sol))):

        dist = max_pairwise_distance_km_coords(
            lons_on[i],  lats_on[i],
            lons_off[j], lats_off[j],
            lons_sol[k], lats_sol[k],
        )
        if dist > max_radius_km:
            continue

        score, coincidence, mean_corr = score_triplet(
            ts_on_all[:, i], ts_off_all[:, j], ts_sol_all[:, k],
            threshold, w_coincidence, w_correlation,
        )

        records.append({
            "score":         score,
            "coincidence":   coincidence,
            "mean_corr":     mean_corr,
            "dist_km":       dist,
            "onshore_y_idx":  int(ys_on[i]),
            "onshore_x_idx":  int(xs_on[i]),
            "onshore_x":      float(lons_on[i]),
            "onshore_y":      float(lats_on[i]),
            "offshore_y_idx": int(ys_off[j]),
            "offshore_x_idx": int(xs_off[j]),
            "offshore_x":     float(lons_off[j]),
            "offshore_y":     float(lats_off[j]),
            "solar_y_idx":    int(ys_sol[k]),
            "solar_x_idx":    int(xs_sol[k]),
            "solar_x":        float(lons_sol[k]),
            "solar_y":        float(lats_sol[k]),
        })

    records.sort(key=lambda r: r["score"], reverse=True)
    return records[:top_n]


# ── Greedy search + neighbourhood top-N ──────────────────────────────────────

def greedy_screen(
    ys_on,  xs_on,  lons_on,  lats_on,
    ys_off, xs_off, lons_off, lats_off,
    ys_sol, xs_sol, lons_sol, lats_sol,
    ts_on_all:  np.ndarray,
    ts_off_all: np.ndarray,
    ts_sol_all: np.ndarray,
    cf_on, cf_off, cf_sol,
    max_radius_km:  float,
    threshold:      float,
    w_coincidence:  float,
    w_correlation:  float,
    top_n:          int,
) -> list[dict]:
    """
    Greedy sequential selection tried from all three anchor techs,
    followed by neighbourhood search to return top N diverse results.
    """
    print("  → Using greedy search (candidate space too large for brute force)")

    def _find_best_triplet(anchor: str):
        if anchor == "wind_onshore":
            best_i = int(np.argmax(ts_on_all.mean(axis=0)))

            best_j, best_score = None, -np.inf
            for j in range(len(ys_off)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i], lats_on[best_i],
                    lons_off[j],     lats_off[j],
                    lons_on[best_i], lats_on[best_i],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, j], ts_on_all[:, best_i],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_j = j
            if best_j is None:
                return None

            best_k, best_score = None, -np.inf
            for k in range(len(ys_sol)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i],  lats_on[best_i],
                    lons_off[best_j], lats_off[best_j],
                    lons_sol[k],      lats_sol[k],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, best_j], ts_sol_all[:, k],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_k = k
            if best_k is None:
                return None

        elif anchor == "wind_offshore":
            best_j = int(np.argmax(ts_off_all.mean(axis=0)))

            best_i, best_score = None, -np.inf
            for i in range(len(ys_on)):
                d = max_pairwise_distance_km_coords(
                    lons_on[i],       lats_on[i],
                    lons_off[best_j], lats_off[best_j],
                    lons_on[i],       lats_on[i],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, i], ts_off_all[:, best_j], ts_on_all[:, i],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_i = i
            if best_i is None:
                return None

            best_k, best_score = None, -np.inf
            for k in range(len(ys_sol)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i],  lats_on[best_i],
                    lons_off[best_j], lats_off[best_j],
                    lons_sol[k],      lats_sol[k],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, best_j], ts_sol_all[:, k],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_k = k
            if best_k is None:
                return None

        else:  # anchor == "solar"
            best_k = int(np.argmax(ts_sol_all.mean(axis=0)))

            best_i, best_score = None, -np.inf
            for i in range(len(ys_on)):
                d = max_pairwise_distance_km_coords(
                    lons_on[i],       lats_on[i],
                    lons_sol[best_k], lats_sol[best_k],
                    lons_on[i],       lats_on[i],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, i], ts_sol_all[:, best_k], ts_on_all[:, i],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_i = i
            if best_i is None:
                return None

            best_j, best_score = None, -np.inf
            for j in range(len(ys_off)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i],  lats_on[best_i],
                    lons_off[j],      lats_off[j],
                    lons_sol[best_k], lats_sol[best_k],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, j], ts_sol_all[:, best_k],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_j = j
            if best_j is None:
                return None

        # Swap passes until convergence
        best_score, _, _ = score_triplet(
            ts_on_all[:, best_i], ts_off_all[:, best_j], ts_sol_all[:, best_k],
            threshold, w_coincidence, w_correlation,
        )
        improved = True
        while improved:
            improved = False
            for i in range(len(ys_on)):
                d = max_pairwise_distance_km_coords(
                    lons_on[i],       lats_on[i],
                    lons_off[best_j], lats_off[best_j],
                    lons_sol[best_k], lats_sol[best_k],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, i], ts_off_all[:, best_j], ts_sol_all[:, best_k],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_i = i
                    improved = True
            for j in range(len(ys_off)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i], lats_on[best_i],
                    lons_off[j],     lats_off[j],
                    lons_sol[best_k], lats_sol[best_k],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, j], ts_sol_all[:, best_k],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_j = j
                    improved = True
            for k in range(len(ys_sol)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i],  lats_on[best_i],
                    lons_off[best_j], lats_off[best_j],
                    lons_sol[k],      lats_sol[k],
                )
                if d > max_radius_km:
                    continue
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, best_j], ts_sol_all[:, k],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_k = k
                    improved = True

        return best_i, best_j, best_k, best_score

    # Run greedy from all three anchors
    anchor_results = []
    for anchor in ["wind_onshore", "wind_offshore", "solar"]:
        print(f"    trying anchor: {anchor}...", end=" ", flush=True)
        result = _find_best_triplet(anchor)
        if result is None:
            print("no valid triplet found")
            continue
        bi, bj, bk, bs = result
        anchor_results.append((bi, bj, bk, bs))
        print(f"score={bs:.4f}")

    if not anchor_results:
        return []

    # Neighbourhood search around all anchor results
    records = {}

    def _add(i, j, k):
        key = (i, j, k)
        if key in records:
            return
        d = max_pairwise_distance_km_coords(
            lons_on[i],  lats_on[i],
            lons_off[j], lats_off[j],
            lons_sol[k], lats_sol[k],
        )
        if d > max_radius_km:
            return
        score, coincidence, mean_corr = score_triplet(
            ts_on_all[:, i], ts_off_all[:, j], ts_sol_all[:, k],
            threshold, w_coincidence, w_correlation,
        )
        records[key] = {
            "score":          score,
            "coincidence":    coincidence,
            "mean_corr":      mean_corr,
            "dist_km":        d,
            "onshore_y_idx":  int(ys_on[i]),
            "onshore_x_idx":  int(xs_on[i]),
            "onshore_x":      float(lons_on[i]),
            "onshore_y":      float(lats_on[i]),
            "offshore_y_idx": int(ys_off[j]),
            "offshore_x_idx": int(xs_off[j]),
            "offshore_x":     float(lons_off[j]),
            "offshore_y":     float(lats_off[j]),
            "solar_y_idx":    int(ys_sol[k]),
            "solar_x_idx":    int(xs_sol[k]),
            "solar_x":        float(lons_sol[k]),
            "solar_y":        float(lats_sol[k]),
        }

    for best_i, best_j, best_k, _ in anchor_results:
        _add(best_i, best_j, best_k)
        for i in range(len(ys_on)):
            _add(i, best_j, best_k)
        for j in range(len(ys_off)):
            _add(best_i, j, best_k)
        for k in range(len(ys_sol)):
            _add(best_i, best_j, k)

    sorted_records = sorted(records.values(), key=lambda r: r["score"], reverse=True)
    return sorted_records[:top_n]


# ── Average profiles ──────────────────────────────────────────────────────────

def save_average_profiles(cc: str, year: int) -> None:
    src = CF_DIR / f"{cc.lower()}_cf_{year}.csv"
    dst = CF_DIR / f"{cc.lower()}_average_profiles_{year}.csv"
    if not src.exists():
        raise FileNotFoundError(f"National mean CF file not found: {src}")
    df = pd.read_csv(src, parse_dates=["time"])
    df.to_csv(dst, index=False)
    print(f"  Saved average profiles → {dst.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    global _current_country

    config = load_config()
    cfg    = get_pypsa_config(config)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    print("Countries from config:", cfg["countries"])
    for country in cfg["countries"]:
        _current_country = country.upper()
        cc = country.lower()

        print(f"\n{'='*60}")
        print(f"Country: {_current_country}")
        print(f"{'='*60}")

        save_average_profiles(cc, cfg["year"])

        print("  Building CF grids (this takes a few minutes per tech)...")
        cf_grids = {}
        for tech in TECHS:
            print(f"    {tech}...", end=" ", flush=True)
            cf_grids[tech] = build_cf_year(_current_country, tech)
            print("done")

        print("  Identifying candidate cells...")
        candidates   = {}
        ts_matrices  = {}
        coord_arrays = {}

        for tech in TECHS:
            cf   = cf_grids[tech]
            geom = geometry_for_tech(_current_country, tech)

            cell_mean = cf.mean("time")
            inside    = mask_cells_inside(cell_mean, geom)
            mean_vals = cell_mean.values
            valid     = inside & np.isfinite(mean_vals) & (mean_vals > 0)

            ys, xs = np.where(valid)
            n      = len(ys)
            print(f"    {tech}: {n} candidate cells")

            n_timesteps = int(cf.sizes["time"])
            ts_matrix   = np.zeros((n_timesteps, n), dtype=np.float32)
            for idx in range(n):
                ts = extract_cell_timeseries(cf, int(ys[idx]), int(xs[idx]), tech)
                ts_matrix[:, idx] = ts.values

            lons, lats = get_candidate_coords(cf, ys, xs)

            candidates[tech]   = (ys, xs, mean_vals[valid])
            ts_matrices[tech]  = ts_matrix
            coord_arrays[tech] = (lons, lats)

        ys_on,  xs_on,  _ = candidates["wind_onshore"]
        ys_off, xs_off, _ = candidates["wind_offshore"]
        ys_sol, xs_sol, _ = candidates["solar"]

        lons_on,  lats_on  = coord_arrays["wind_onshore"]
        lons_off, lats_off = coord_arrays["wind_offshore"]
        lons_sol, lats_sol = coord_arrays["solar"]

        n_triplets = len(ys_on) * len(ys_off) * len(ys_sol)
        print(f"  Total candidate triplets: {n_triplets:,}")

        screen_kwargs = dict(
            ys_on=ys_on,   xs_on=xs_on,   lons_on=lons_on,   lats_on=lats_on,
            ys_off=ys_off, xs_off=xs_off, lons_off=lons_off, lats_off=lats_off,
            ys_sol=ys_sol, xs_sol=xs_sol, lons_sol=lons_sol, lats_sol=lats_sol,
            ts_on_all=ts_matrices["wind_onshore"],
            ts_off_all=ts_matrices["wind_offshore"],
            ts_sol_all=ts_matrices["solar"],
            cf_on=cf_grids["wind_onshore"],
            cf_off=cf_grids["wind_offshore"],
            cf_sol=cf_grids["solar"],
            max_radius_km=cfg["max_radius_km"],
            threshold=cfg["coincidence_threshold"],
            w_coincidence=cfg["w_coincidence"],
            w_correlation=cfg["w_correlation"],
            top_n=cfg["top_n"],
        )

        if n_triplets <= MAX_TRIPLETS_BRUTE_FORCE:
            top_records = brute_force_screen(**screen_kwargs)
        else:
            top_records = greedy_screen(**screen_kwargs)

        if not top_records:
            print(f"  WARNING: No valid triplets found for {_current_country}. "
                  f"Try increasing max_radius_km (currently {cfg['max_radius_km']} km).")
            continue

        df = pd.DataFrame(top_records)
        df.insert(0, "rank",    range(1, len(df) + 1))
        df.insert(1, "country", _current_country)

        out_path = OUTDIR / f"{cc}_complementarity_top{cfg['top_n']}_{cfg['year']}.csv"
        df.to_csv(out_path, index=False)

        print(f"\n  Top {len(df)} complementary triplets → {out_path.name}")
        print(f"  Best score: {df['score'].iloc[0]:.4f}  "
              f"coincidence: {df['coincidence'].iloc[0]:.3f}  "
              f"mean_corr: {df['mean_corr'].iloc[0]:.3f}  "
              f"dist_km: {df['dist_km'].iloc[0]:.1f}")


if __name__ == "__main__":
    main()
