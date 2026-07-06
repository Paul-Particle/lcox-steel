"""
07b_make_anchor_colocated_cf_timeseries.py

Purpose
-------
For each anchor technology (wind_onshore, wind_offshore, solar), builds a
co-located project scenario: the anchor cell is fixed (P95, from 07's logic),
and the other two technologies each independently propose their own top-N
best complementary candidate cells nearby.

Anchor-based vs free triplet search (08)
-----------------------------------------
07b answers: given this location (P95), what are the best candidate sites to
feed into an optimized mix?" — one cell is fixed first (the anchor), then
candidate partners are shortlisted around it. PyPSA (09) then decides the
actual cost-minimizing capacities across those candidates. This is cheap
and gives a realistic "flagship project" story.

08's brute-force/greedy triplet screen answers a different question: "ignoring
any single anchor, what's the single best-complementing site combo in the
whole country?" More thorough, more expensive, no anchor story. 07b is not a
replacement for 08 — different scenario type.

Method
------
- Anchor cell: fixed P95 cell for the anchor tech (reused from 07)
- Candidate cells (other 2 techs): cells within max_radius_km of the anchor
  (checked first, cheaply) that are also valid (inside geometry, real CF
  data). Scored pairwise (anchor vs candidate) by:
      score = w_coincidence * coincidence - w_correlation * correlation
  Top `n_candidates` (default 3, configurable) kept per tech.
- No distance term in the score — PyPSA prices transmission by distance
  separately, so this would double-count. Radius is a pre-filter/hard feasibility
  cutoff only.
- Offshore: if no valid cell within radius, CF series is all zeros (no
  fallback to nearest-overall — see design doc decision, no exceptions).

Outputs
-------
resources/res_cf/<cc>_cf_2023_bestsite_p95_anchor-<anchor>_candidates.parquet

One file per anchor per country. Wide format: one CF column per candidate
(anchor + up to n_candidates per other tech), metadata columns (x, y,
dist_km, score, coincidence, correlation, n_candidates_found) repeated on
every row — same convention as 07's add_location_metadata.
"""

from pathlib import Path
import importlib.util
import logging

import numpy as np
import pandas as pd
import xarray as xr

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from common._paths import RES_CF, SHAPES_RES
from scripts.res_cf._helpers import annual_cutout_path, haversine_distance_km, load_res_cf_cfg

configure_logging(snakemake)
log = logging.getLogger(__name__)

# ── Import reusable functions from script 07 (same importlib pattern as 08) ──
_spec = importlib.util.spec_from_file_location(
    "bestsite",
    Path(__file__).parent / "07_make_bestsite_cf_timeseries.py"
)
_bestsite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bestsite)

build_cf_year     = _bestsite.build_cf_year
geometry_for_tech = _bestsite.geometry_for_tech
mask_cells_inside = _bestsite.mask_cells_inside
find_p95_cell      = _bestsite.find_p95_cell
get_cell_coords    = _bestsite.get_cell_coords
extract_cell_timeseries = _bestsite.extract_cell_timeseries

if "snakemake" in globals() and hasattr(snakemake, "input"):
    _bestsite.REGIONS_PATH = Path(snakemake.input.regions)
    _bestsite.OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)

_RES_CF_CFG = load_res_cf_cfg()
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _RES_CF_CFG = snakemake.config["res_cf"]

_ANCHOR_CFG = _RES_CF_CFG.get("anchor_colocation", {})

MAX_RADIUS_KM         = float(_ANCHOR_CFG.get("max_radius_km", 100.0))
N_CANDIDATES          = int(_ANCHOR_CFG.get("n_candidates", 3))
COINCIDENCE_THRESHOLD = float(_ANCHOR_CFG.get("coincidence_threshold", 0.20))
W_COINCIDENCE         = float(_ANCHOR_CFG.get("w_coincidence", 0.6))
W_CORRELATION         = float(_ANCHOR_CFG.get("w_correlation", 0.4))

YEAR = 2023
OUTDIR = RES_CF
COUNTRIES = ["de"]  # standalone default

CUTOUT_DIR = None
if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    COUNTRIES = [snakemake.wildcards.cf_area.lower()]
    CUTOUT_DIR = Path(snakemake.input.cutout)

TECHS = ["wind_onshore", "wind_offshore", "solar"]

def score_candidate(ts_anchor: np.ndarray, ts_candidate: np.ndarray,
                     threshold: float, w_coincidence: float,
                     w_correlation: float) -> tuple[float, float, float]:
    """Pairwise complementarity score between anchor and one candidate cell.

    Returns (score, coincidence, correlation).
    """
    combined = (ts_anchor + ts_candidate) / 2.0
    coincidence = float(np.mean(combined > threshold))

    correlation = float(np.corrcoef(ts_anchor, ts_candidate)[0, 1])

    score = w_coincidence * coincidence - w_correlation * correlation
    return score, coincidence, correlation


def find_top_candidates(anchor_x: float, anchor_y: float, cf_year: xr.DataArray,
                         geom, max_radius_km: float, n_candidates: int,
                         ts_anchor: np.ndarray, threshold: float,
                         w_coincidence: float, w_correlation: float) -> list[dict]:
    """Find up to `n_candidates` cells within max_radius_km, ranked by score.

    Order of operations (radius first, since it's cheap and shrinks the grid
    before the more expensive geometry/validity check runs):
    1. Compute distance from anchor to every grid cell
    2. Keep only cells within max_radius_km
    3. Of those, keep valid cells (inside geometry, real CF data, mean CF > 0)
    4. Score survivors against the anchor, sort, keep top n_candidates

    Returns a list of dicts (one per candidate), each with:
    y_idx, x_idx, x, y, dist_km, score, coincidence, correlation.
    Empty list if none found within radius.
    """
    xs_all = cf_year.x.values
    ys_all = cf_year.y.values
    xx, yy = np.meshgrid(xs_all, ys_all)

    dist_km = haversine_distance_km(anchor_x, anchor_y, xx, yy)
    within_radius = dist_km <= max_radius_km

    if not np.any(within_radius):
        return []

    cell_mean = cf_year.mean("time")
    inside = mask_cells_inside(cell_mean, geom)
    mean_vals = cell_mean.values
    valid = within_radius & inside & np.isfinite(mean_vals) & (mean_vals > 0)

    ys, xs = np.where(valid)

    if len(ys) == 0:
        return []

    candidates = []
    for y_idx, x_idx in zip(ys, xs):
        d = float(dist_km[y_idx, x_idx])
        ts_candidate = extract_cell_timeseries(cf_year, int(y_idx), int(x_idx)).values
        score, coincidence, correlation = score_candidate(
            ts_anchor, ts_candidate, threshold, w_coincidence, w_correlation
        )
        x, y = get_cell_coords(cf_year, int(y_idx), int(x_idx))
        candidates.append({
            "y_idx": int(y_idx), "x_idx": int(x_idx),
            "x": x, "y": y, "dist_km": d,
            "score": score, "coincidence": coincidence, "correlation": correlation,
        })

    candidates.sort(key=lambda c: c["score"], reverse=True)
    return candidates[:n_candidates]


def build_anchor_scenario(cutout_path: Path, country_upper: str, anchor_tech: str,
                           max_radius_km: float, n_candidates: int,
                           threshold: float, w_coincidence: float,
                           w_correlation: float) -> tuple[pd.Series, dict[str, list[dict]], dict[str, list[pd.Series]]]:
    """Build one anchor's full co-located scenario: anchor series + candidate
    lists + candidate series for the other two techs.

    Offshore special case: if no valid offshore cell exists within radius,
    a single zero-CF placeholder candidate is used instead of leaving the
    tech empty (no fallback to nearest-overall — see design decision).
    """
    anchor_geom = geometry_for_tech(country_upper, anchor_tech)
    anchor_cf_year = build_cf_year(cutout_path, anchor_tech)
    anchor_y_idx, anchor_x_idx = find_p95_cell(anchor_cf_year, anchor_geom)
    anchor_series = extract_cell_timeseries(anchor_cf_year, anchor_y_idx, anchor_x_idx)
    anchor_x, anchor_y = get_cell_coords(anchor_cf_year, anchor_y_idx, anchor_x_idx)

    other_techs = [t for t in TECHS if t != anchor_tech]

    candidates_meta = {}
    candidates_series = {}

    for tech in other_techs:
        tech_cf_year = build_cf_year(cutout_path, tech)
        tech_geom = geometry_for_tech(country_upper, tech)

        # TODO(perf): build_cf_year() rebuilds each tech's full-country CF grid from
        # scratch here, and this function runs once per anchor (3x per country) —
        # so the same tech's grid gets recomputed up to 3x per country. Consider
        # caching cf_year per (country, tech) across anchor runs if this gets slow.

        found = find_top_candidates(
            anchor_x, anchor_y, tech_cf_year, tech_geom,
            max_radius_km, n_candidates,
            anchor_series.values, threshold, w_coincidence, w_correlation,
        )

        if not found and tech == "wind_offshore":
            zero_series = pd.Series(0.0, index=anchor_series.index, name="cf")
            candidates_meta[tech] = [{
                "y_idx": None, "x_idx": None, "x": None, "y": None,
                "dist_km": None, "score": None,
                "coincidence": None, "correlation": None,
            }]
            candidates_series[tech] = [zero_series]
            continue

        meta_list = []
        series_list = []
        for cand in found:
            ts = extract_cell_timeseries(tech_cf_year, cand["y_idx"], cand["x_idx"])
            meta_list.append(cand)
            series_list.append(ts)

        candidates_meta[tech] = meta_list
        candidates_series[tech] = series_list

    return anchor_series, candidates_meta, candidates_series

def assemble_output_dataframe(anchor_tech: str, anchor_series: pd.Series,
                               candidates_meta: dict[str, list[dict]],
                               candidates_series: dict[str, list[pd.Series]],
                               max_radius_km: float, n_candidates: int,
                               threshold: float, w_coincidence: float,
                               w_correlation: float) -> pd.DataFrame:
    """Assemble the wide-format output: time, anchor_tech, <anchor>_cf,
    <tech>_cand{n}_cf columns, plus repeated metadata columns
    (_x, _y, _dist_km, _score, _coincidence, _correlation), one
    _n_candidates_found count column per other tech, and the run's config
    settings (so results are traceable if config changes later).

    Note: the offshore zero-CF placeholder (see build_anchor_scenario) is
    not a real candidate, so it's excluded from n_candidates_found (reported
    as 0) even though its CF column is still written.
    """
    df = pd.DataFrame({"time": anchor_series.index})
    df["anchor_tech"] = anchor_tech
    df[f"{anchor_tech}_cf"] = anchor_series.values

    df["cfg_max_radius_km"] = max_radius_km
    df["cfg_n_candidates"] = n_candidates
    df["cfg_coincidence_threshold"] = threshold
    df["cfg_w_coincidence"] = w_coincidence
    df["cfg_w_correlation"] = w_correlation
    df["run_timestamp_utc"] = pd.Timestamp.utcnow().isoformat()


    for tech, series_list in candidates_series.items():
        meta_list = candidates_meta[tech]

        is_placeholder = len(meta_list) == 1 and meta_list[0]["dist_km"] is None
        df[f"{tech}_n_candidates_found"] = 0 if is_placeholder else len(series_list)

        for i, (series, meta) in enumerate(zip(series_list, meta_list), start=1):
            df[f"{tech}_cand{i}_cf"] = series.values
            for key in ("x", "y", "dist_km", "score", "coincidence", "correlation"):
                df[f"{tech}_cand{i}_{key}"] = meta[key]

    return df


def main() -> None:
    for cc in COUNTRIES:
        country_upper = cc.upper()
        cutout = CUTOUT_DIR or annual_cutout_path(cc, YEAR)

        for anchor_tech in TECHS:
            anchor_series, candidates_meta, candidates_series = build_anchor_scenario(
                cutout, country_upper, anchor_tech,
                MAX_RADIUS_KM, N_CANDIDATES,
                COINCIDENCE_THRESHOLD, W_COINCIDENCE, W_CORRELATION,
            )

            df = assemble_output_dataframe(
                anchor_tech, anchor_series, candidates_meta, candidates_series,
                MAX_RADIUS_KM, N_CANDIDATES,
                COINCIDENCE_THRESHOLD, W_COINCIDENCE, W_CORRELATION,
            )

            out_path = OUTDIR / f"{cc}_cf_{YEAR}_bestsite_p95_anchor-{anchor_tech}_candidates.parquet"
            df.to_parquet(out_path, index=False)
            log.info(f"{country_upper}: wrote {out_path.name}")


if __name__ == "__main__":
    main()
