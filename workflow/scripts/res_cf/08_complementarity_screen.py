"""
08_complementarity_screen.py

Purpose
-------
Identifies the best combination of RES cells (onshore wind, offshore wind, solar)
for each country by screening for complementarity — i.e. combinations whose hourly
profiles interlock well to maximise electrolyser utilisation.

This script produces two outputs per country:
  (a) Top N complementary cell triplets for the best-site RES mix scenario
  (b) National mean profiles for the country-average scenario

Both outputs feed directly into the PyPSA optimisation

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
resources/res_cf/<cc>_complementarity_top{N}_2023.parquet
resources/res_cf/<cc>_average_profiles_2023.parquet
"""

from __future__ import annotations

import importlib.util
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import logging
from common._logging import configure_logging, progress
from common._paths import REPO_ROOT, RES_CF
from scripts.res_cf._helpers import annual_cutout_path, haversine_distance_km

if "snakemake" not in globals():
    from common._stubs import snakemake
configure_logging(snakemake)
log = logging.getLogger(__name__)

# ── Import reusable functions from script 07 ─────────────────────────────────
_spec = importlib.util.spec_from_file_location(
    "bestsite",
    Path(__file__).parent / "07_make_bestsite_cf_timeseries.py"
)
_bestsite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bestsite)

# In Snakemake mode, override the imported module's path globals so that
# geometry_for_tech reads from the rule's declared input files, not 07's defaults.
if "snakemake" in globals() and hasattr(snakemake, "input"):
    _bestsite.REGIONS_PATH = Path(snakemake.input.regions)
    _bestsite.OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)

SM_CUTOUT_PATH: Path | None = None
SM_SCREEN_OUT:  Path | None = None
SM_CF_OUT:      Path | None = None

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    SM_CUTOUT_PATH = Path(snakemake.input.cutout)
    SM_SCREEN_OUT  = Path(snakemake.output.screen)
    SM_CF_OUT      = Path(snakemake.output.cf)

build_cf_year           = _bestsite.build_cf_year
extract_cell_timeseries = _bestsite.extract_cell_timeseries
geometry_for_tech       = _bestsite.geometry_for_tech
# haversine_distance_km   = _bestsite.haversine_distance_km #imported from helpers now
mask_cells_inside       = _bestsite.mask_cells_inside

# ── Paths ─────────────────────────────────────────────────────────────────────


OUTDIR = RES_CF
CF_DIR = RES_CF


TECHS = ["wind_onshore", "wind_offshore", "solar"]

MAX_TRIPLETS_BRUTE_FORCE = 500_000


# ── Config ────────────────────────────────────────────────────────────────────

def load_config() -> dict:
    with open(REPO_ROOT / "config/config.yaml") as f:
        c = yaml.safe_load(f)


# just one config here, no def get_pypsa_config(config: dict) -> dict:
    rc = c["res_cf"]
    comp = rc.get("complementarity", {})
    return {
        "countries":             [x.upper() for x in rc["countries"]],
        "year":                  2023,
        "top_n":                 comp.get("top_n", 10),
        "coincidence_threshold": comp.get("coincidence_threshold", 0.20),
        "w_coincidence":         comp.get("w_coincidence", 0.6),
        "w_correlation":         comp.get("w_correlation", 0.4),
        "max_radius_km":         float(comp.get("max_radius_km", 300.0)),
        "quality_floor":         float(comp.get("quality_floor", 0.90)),
        "max_triplets_brute_force": comp.get("max_triplets_brute_force", 500_000),
    }


# ── Coordinate helpers ────────────────────────────────────────────────────────

def get_candidate_coords(cf, ys: np.ndarray, xs: np.ndarray):
    """
    Return (lons, lats) arrays for all candidate cells.
    ys, xs are grid indices into the DataArray dimensions.
    """
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
    n_total = len(ys_on) * len(ys_off) * len(ys_sol)
    triplets = product(range(len(ys_on)), range(len(ys_off)), range(len(ys_sol)))
    for i, j, k in progress(triplets, desc="brute-force triplets", total=n_total, unit="trip"):

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

MAX_TRIPLETS_BRUTE_FORCE = load_config()['max_triplets_brute_force']


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
    # WIP NOTE: The three `if anchor == ...` blocks share ~80% of their body and
    # are the highest-value refactor target in this file. Leave until the broader
    # cutout-cache refactor stabilises the API this script depends on.
    log.info("using greedy search (candidate space too large for brute force)")

    def _find_best_triplet(anchor: str):
        """
        Run greedy from a given anchor tech.
        Returns (best_i, best_j, best_k, best_score) or None if no valid triplet found.
        """
        if anchor == "wind_onshore":
            best_i = int(np.argmax(ts_on_all.mean(axis=0)))

            # Step 2: best offshore complement
            best_j, best_score = None, -np.inf
            for j in range(len(ys_off)):
                d = max_pairwise_distance_km_coords(
                    lons_on[best_i], lats_on[best_i],
                    lons_off[j],     lats_off[j],
                    lons_on[best_i], lats_on[best_i],
                )
                if d > max_radius_km:
                    continue
                # WIP NOTE: ts_on_all[:, best_i] is passed as the solar placeholder while
                # solar hasn't been chosen yet. This is an intentional bootstrap (using the
                # anchor as a stand-in) or a copy-paste bug — needs verification.
                s, _, _ = score_triplet(
                    ts_on_all[:, best_i], ts_off_all[:, j], ts_on_all[:, best_i],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_j = j
            if best_j is None:
                return None

            # Step 3: best solar complement
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

            # Step 2: best onshore complement
            best_i, best_score = None, -np.inf
            for i in range(len(ys_on)):
                d = max_pairwise_distance_km_coords(
                    lons_on[i],       lats_on[i],
                    lons_off[best_j], lats_off[best_j],
                    lons_on[i],       lats_on[i],
                )
                if d > max_radius_km:
                    continue
                # WIP NOTE: ts_on_all[:, i] is used as solar placeholder — see comment above.
                s, _, _ = score_triplet(
                    ts_on_all[:, i], ts_off_all[:, best_j], ts_on_all[:, i],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_i = i
            if best_i is None:
                return None

            # Step 3: best solar complement
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

            # Step 2: best onshore complement
            best_i, best_score = None, -np.inf
            for i in range(len(ys_on)):
                d = max_pairwise_distance_km_coords(
                    lons_on[i],       lats_on[i],
                    lons_sol[best_k], lats_sol[best_k],
                    lons_on[i],       lats_on[i],
                )
                if d > max_radius_km:
                    continue
                # WIP NOTE: ts_on_all[:, i] is used as offshore placeholder — see comment above.
                s, _, _ = score_triplet(
                    ts_on_all[:, i], ts_sol_all[:, best_k], ts_on_all[:, i],
                    threshold, w_coincidence, w_correlation,
                )
                if s > best_score:
                    best_score = s
                    best_i = i
            if best_i is None:
                return None

            # Step 3: best offshore complement
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

        # Step 4: swap passes until no improvement
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

    # ── Run greedy from all three anchors ─────────────────────────────────────
    anchor_results = []
    for anchor in ["wind_onshore", "wind_offshore", "solar"]:
        log.info(f"trying anchor: {anchor} ...")
        result = _find_best_triplet(anchor)
        if result is None:
            log.info(f"anchor={anchor}: no valid triplet found")
            continue
        bi, bj, bk, bs = result
        anchor_results.append((bi, bj, bk, bs))
        log.info(f"anchor={anchor}: score={bs:.4f}")

    if not anchor_results:
        return []

    # ── Neighbourhood search around all anchor results ────────────────────────
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

    # Add all anchor results and their neighbourhoods
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
    """Combine per-tech country-average parquets into a single diagnostic file.

    Source files are the outputs of `build_res_cf_profile`; they are written once
    per tech, not per year, so we glob by the area+tech prefix.
    """
    frames = {}
    for tech in TECHS:
        pattern = f"{cc.lower()}_{tech.replace('_', '-')}_country-average_*.parquet"
        matches = sorted(CF_DIR.glob(pattern))
        if not matches:
            log.warning(f"no country-average parquet for {tech} ({cc}) — skipping")
            continue
        src = matches[-1]  # most recent if multiple
        df_tech = pd.read_parquet(src)
        # The country-average parquet has a single CF column; rename it to the tech.
        cf_col = [c for c in df_tech.columns if c != "time"]
        frames[tech.replace("_", "-")] = df_tech.set_index("time")[cf_col[0]]

    if not frames:
        log.warning(f"no country-average sources found for {cc} — skipping average profiles")
        return

    dst = CF_DIR / f"{cc.lower()}_average_profiles_{year}.parquet"
    pd.DataFrame(frames).to_parquet(dst, index=True)
    log.info(f"saved average profiles → {dst.name}")


def _write_sm_outputs(
    top_records: list[dict],
    ts_matrices: dict[str, np.ndarray],
    candidates: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> None:
    """Write Snakemake outputs for the top-1 complementarity triplet.

    screen: top-N metadata parquet (all triplets, not just top-1)
    cf:     3-column timeseries parquet (wind-onshore, wind-offshore, solar)
            for the top-1 triplet; columns match the naming convention used
            by build_anchored_cf so solve_network.py can consume it unchanged.
    """
    if "snakemake" not in globals():
        return

    if SM_SCREEN_OUT is not None:
        df_screen = pd.DataFrame(top_records)
        df_screen.insert(0, "rank", range(1, len(df_screen) + 1))
        df_screen.to_parquet(SM_SCREEN_OUT, index=False)
        log.info(f"wrote screen output → {SM_SCREEN_OUT.name}")

    if SM_CF_OUT is not None and top_records:
        best = top_records[0]
        ys_on,  xs_on,  _ = candidates["wind_onshore"]
        ys_off, xs_off, _ = candidates["wind_offshore"]
        ys_sol, xs_sol, _ = candidates["solar"]

        i = int(np.where(
            (ys_on == best["onshore_y_idx"]) & (xs_on == best["onshore_x_idx"])
        )[0][0])
        j = int(np.where(
            (ys_off == best["offshore_y_idx"]) & (xs_off == best["offshore_x_idx"])
        )[0][0])
        k = int(np.where(
            (ys_sol == best["solar_y_idx"]) & (xs_sol == best["solar_x_idx"])
        )[0][0])

        pd.DataFrame({
            "wind-onshore":  ts_matrices["wind_onshore"][:, i],
            "wind-offshore": ts_matrices["wind_offshore"][:, j],
            "solar":         ts_matrices["solar"][:, k],
        }).to_parquet(SM_CF_OUT, index=False)
        log.info(f"wrote top-1 CF output → {SM_CF_OUT.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    sm_country = None
    if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
        sm_country = snakemake.wildcards.cf_area.upper()

    config = load_config()
    # no second config

    OUTDIR.mkdir(parents=True, exist_ok=True)
    countries = [sm_country] if sm_country is not None else [c.upper() for c in config["countries"]]
    log.info(f"countries: {countries}")
    for country in countries:
        _current_country = country.upper()
        cc = country.lower()

        log.info(f"country={_current_country} start")



        # 1. Save average profiles
        save_average_profiles(cc, config["year"])

        # 2. Build full CF grids
        log.info("building CF grids (this takes a few minutes per tech)")
        cutout = SM_CUTOUT_PATH if SM_CUTOUT_PATH is not None else annual_cutout_path(cc, config["year"])
        cf_grids = {}
        for tech in TECHS:

            cf_grids[tech] = build_cf_year(cutout, tech)
            log.info(f"tech={tech} CF grid built")

        # 3. Get candidate cells per tech
        log.info("identifying candidate cells")
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
            log.info(f"tech={tech}: {n} candidate cells")

            # Extract time series matrix (8760, n)
            n_timesteps = int(cf.sizes["time"])
            ts_matrix   = np.zeros((n_timesteps, n), dtype=np.float32)
            for idx in progress(range(n), desc=f"extract {tech} cells", total=n, unit="cell"):
                ts = extract_cell_timeseries(cf, int(ys[idx]), int(xs[idx]))
                ts_matrix[:, idx] = ts.values

            # Pre-compute actual coordinates
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
        log.info(f"total candidate triplets: {n_triplets}")

        # 4. Screen triplets
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
            max_radius_km=config["max_radius_km"],
            threshold=config["coincidence_threshold"],
            w_coincidence=config["w_coincidence"],
            w_correlation=config["w_correlation"],
            top_n=config["top_n"],
        )

        if n_triplets <= config["max_triplets_brute_force"]:
            top_records = brute_force_screen(**screen_kwargs)
        else:
            top_records = greedy_screen(**screen_kwargs)

        if not top_records:
            log.warning(
                f"no valid triplets found for {_current_country}. "
                f"Try increasing max_radius_km (currently {config['max_radius_km']} km).")
            continue

        # 5. Save results
        df = pd.DataFrame(top_records)
        df.insert(0, "rank",    range(1, len(df) + 1))
        df.insert(1, "country", _current_country)

        out_path = OUTDIR / f"{cc}_complementarity_top{config['top_n']}_{config['year']}.parquet"
        df.to_parquet(out_path, index=False)

        log.info(f"top {len(df)} complementary triplets → {out_path.name}")
        best = df.iloc[0]
        log.info(
            f"best score={best['score']:.4f} coincidence={best['coincidence']:.3f} "
            f"mean_corr={best['mean_corr']:.3f} dist_km={best['dist_km']:.1f}"
        )

        _write_sm_outputs(top_records, ts_matrices, candidates)


if __name__ == "__main__":
    main()