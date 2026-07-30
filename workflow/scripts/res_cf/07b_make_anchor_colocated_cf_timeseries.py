"""
07b_make_anchor_colocated_cf_timeseries.py

Purpose
-------
For one anchor technology (wind-onshore, wind-offshore or solar), builds a
co-located project scenario: the anchor cell is fixed (P95, from 07's logic),
and the other two technologies each independently propose their own top-N
best complementary candidate cells nearby.

07b answers: "given this location (P95), what are the best candidate sites to
feed into an optimized mix?" — one cell is fixed first (the anchor), then
candidate partners are shortlisted around it. PyPSA (solve_network) then
decides the actual cost-minimizing capacities across those candidates. This is
cheap and gives a realistic "flagship project" story. (A free triplet search
over the whole country — the removed script 08 — asked a different question
and was retired in favour of this anchor-based scenario type.)

Method
------
- Anchor cell: fixed P95 cell for the anchor tech (reused from 07)
- Candidate cells (other 2 techs): cells within max_radius_km of the anchor
  (checked first, cheaply) that are also valid (inside geometry, real CF
  data). Scored pairwise (anchor vs candidate) by:
      score = w_coincidence * coincidence - w_correlation * correlation
  Top `n_candidates` (from the variant wildcard, e.g. anchor-colo-n3) kept
  per tech.
- No distance term in the score — PyPSA prices transmission by distance
  separately (candidate → anchor HVDC links in multi-site mode), so this
  would double-count. Radius is a pre-filter/hard feasibility cutoff only.
- Offshore: if no valid cell within radius, the tech is simply absent from
  the output (no zero-CF placeholder, no fallback to nearest-overall — see
  design doc decision). Its absence is recorded in the file metadata.

Outputs
-------
resources/res_cf/{cf_area}_{tech}_anchor-colo-n{N}_{start_date}_{end_date}.parquet

One file per (anchor tech, country, period), in the multi-site CF contract
consumed by solve_network (same as 03c's multi-n* outputs):

- table: DatetimeIndex named `time`; one CF column per site, named
  `{tech}@anchor` (the fixed anchor cell) and `{tech}@c00`, `{tech}@c01`, …
  (candidates, best score first). Tech keys are hyphenated to match
  assumptions.yaml (`wind-onshore`, not `wind_onshore`).
- Arrow schema metadata:
  - `site_coords`: JSON {column: {lat, lon}} — required by solve_network's
    multi-site assembly.
  - `demand_site`: the anchor column name. Its presence tells solve_network
    to enter multi-site mode with the anchor cell as the plant/demand bus
    (no config/sites_*.yaml overlay needed — the anchor IS the plant site).
  - `anchor_colocation`: JSON with per-candidate diagnostics (score,
    coincidence, correlation, dist_km), candidate counts, the run's config
    and a timestamp — the traceability data that used to be repeated
    row-wise as table columns.
"""

from pathlib import Path
import importlib.util
import json
import logging

import numpy as np
import pandas as pd
import xarray as xr
import atlite
import pyarrow as pa
import pyarrow.parquet as pq

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from common._paths import CUTOUTS, RES_CF, SHAPES_RES
from scripts.res_cf._helpers import (
    eligibility_mask_2d,
    eligibility_weights,
    haversine_distance_km,
    load_res_cf_cfg,
    mask_cells_inside,
    pick_p95_cell,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

# ── Import reusable functions from script 07 (importlib, not a package import) ──
_spec = importlib.util.spec_from_file_location(
    "bestsite",
    Path(__file__).parent / "07_make_bestsite_cf_timeseries.py"
)
_bestsite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bestsite)

build_cf_year     = _bestsite.build_cf_year
geometry_for_tech = _bestsite.geometry_for_tech
get_cell_coords    = _bestsite.get_cell_coords
extract_cell_timeseries = _bestsite.extract_cell_timeseries

# Standalone defaults (res_cf hardcoded-default pattern)
_CF_AREA = "de"
_ANCHOR_TECH = "wind-onshore"   # hyphenated, as in the tech wildcard
_START_DATE = "20230101"
_END_DATE = "20231231"
_CUTOUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_REGIONS_PATH = SHAPES_RES / "de_geo.parquet"
_OFFSHORE_REGIONS_PATH = SHAPES_RES / "de_offshore_geo.parquet"
_ANCHOR_CFG = load_res_cf_cfg().get("anchor_colocation", {})
_MIN_LAND_FRACTION = float(load_res_cf_cfg().get("min_land_fraction", 0.0))
_ELIGIBILITY_SOURCE = load_res_cf_cfg().get("eligibility_source", "indicatormatrix")
_N_CANDIDATES = int(_ANCHOR_CFG.get("n_candidates", 3))
_OUT = RES_CF / f"{_CF_AREA}_{_ANCHOR_TECH}_anchor-colo-n{_N_CANDIDATES}_{_START_DATE}_{_END_DATE}.parquet"

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _CF_AREA = snakemake.wildcards.cf_area.lower()
    _ANCHOR_TECH = snakemake.wildcards.tech
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _CUTOUT_PATH = Path(snakemake.input.cutout)
    _REGIONS_PATH = Path(snakemake.input.regions)
    _OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    _ANCHOR_CFG = snakemake.params.anchor_colocation
    _MIN_LAND_FRACTION = float(snakemake.config["res_cf"].get("min_land_fraction", 0.0))
    _ELIGIBILITY_SOURCE = snakemake.config["res_cf"].get("eligibility_source", "indicatormatrix")
    # The variant wildcard (anchor-colo-n{N}) pins n_candidates per scenario;
    # it overrides the config default so the filename always tells the truth.
    _N_CANDIDATES = int(snakemake.wildcards.variant.rsplit("-n", 1)[1])
    _OUT = Path(snakemake.output[0])

_bestsite.REGIONS_PATH = _REGIONS_PATH
_bestsite.OFFSHORE_REGIONS_PATH = _OFFSHORE_REGIONS_PATH

MAX_RADIUS_KM         = float(_ANCHOR_CFG.get("max_radius_km", 100.0))
COINCIDENCE_THRESHOLD = float(_ANCHOR_CFG.get("coincidence_threshold", 0.20))
W_COINCIDENCE         = float(_ANCHOR_CFG.get("w_coincidence", 0.6))
W_CORRELATION         = float(_ANCHOR_CFG.get("w_correlation", 0.4))

TECHS = ["wind_onshore", "wind_offshore", "solar"]  # snake_case: 07's convention


def _hyphen(tech: str) -> str:
    """Internal snake_case tech name → hyphenated key used in assumptions/columns."""
    return tech.replace("_", "-")


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
                         w_coincidence: float, w_correlation: float,
                         eligible: np.ndarray | None = None) -> list[dict]:
    """Find up to `n_candidates` cells within max_radius_km, ranked by score.

    Order of operations (radius first, since it's cheap and shrinks the grid
    before the more expensive geometry/validity check runs):
    1. Compute distance from anchor to every grid cell
    2. Keep only cells within max_radius_km
    3. Of those, keep valid cells (inside geometry, real CF data, mean CF > 0,
       and — for onshore techs — land-sea eligible per `eligible`)
    4. Score survivors against the anchor, sort, keep top n_candidates

    `eligible` is an optional (y, x) boolean land-sea mask (#41). When given,
    sea-contaminated coastal border cells are dropped from candidate selection;
    pass None (or an all-True mask) for offshore, which is sited on sea.

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
    if eligible is not None:
        valid = valid & eligible

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
                           w_correlation: float,
                           min_land_fraction: float = 0.0,
                           eligibility_source: str = "indicatormatrix") -> tuple[pd.Series, tuple[float, float], dict[str, list[dict]], dict[str, list[pd.Series]]]:
    """Build one anchor's full co-located scenario: anchor series + anchor
    (x, y) + candidate lists + candidate series for the other two techs.

    A tech with no valid cell within radius (typically offshore for inland
    anchors) gets an empty list — no placeholder series. The output assembly
    simply omits that tech and records the zero count in the file metadata.

    `min_land_fraction` is the #41 land-sea eligibility cutoff, applied to the
    onshore anchor P95 selection and to onshore candidate cells; offshore is
    never thresholded (sited on sea by design).
    """
    anchor_geom = geometry_for_tech(country_upper, anchor_tech)
    anchor_cf_year = build_cf_year(cutout_path, anchor_tech)

    co = atlite.Cutout(path=str(cutout_path))
    anchor_min_lf = min_land_fraction if anchor_tech != "wind_offshore" else 0.0
    anchor_weights = eligibility_weights(co, anchor_geom, anchor_min_lf, eligibility_source)

    anchor_y_idx, anchor_x_idx = pick_p95_cell(anchor_cf_year.mean("time"), anchor_weights)
    anchor_series = extract_cell_timeseries(anchor_cf_year, anchor_y_idx, anchor_x_idx)
    anchor_x, anchor_y = get_cell_coords(anchor_cf_year, anchor_y_idx, anchor_x_idx)

    other_techs = [t for t in TECHS if t != anchor_tech]

    candidates_meta = {}
    candidates_series = {}

    for tech in other_techs:
        tech_cf_year = build_cf_year(cutout_path, tech)
        tech_geom = geometry_for_tech(country_upper, tech)

        tech_min_lf = min_land_fraction if tech != "wind_offshore" else 0.0
        eligible = eligibility_mask_2d(co, tech_geom, tech_min_lf, eligibility_source)

        found = find_top_candidates(
            anchor_x, anchor_y, tech_cf_year, tech_geom,
            max_radius_km, n_candidates,
            anchor_series.values, threshold, w_coincidence, w_correlation,
            eligible=eligible,
        )

        if not found:
            log.warning(
                f"{country_upper}: no valid {tech} cell within "
                f"{max_radius_km:.0f} km of the {anchor_tech} anchor — "
                "tech omitted from output"
            )

        meta_list = []
        series_list = []
        for cand in found:
            ts = extract_cell_timeseries(tech_cf_year, cand["y_idx"], cand["x_idx"])
            meta_list.append(cand)
            series_list.append(ts)

        candidates_meta[tech] = meta_list
        candidates_series[tech] = series_list

    return anchor_series, (anchor_x, anchor_y), candidates_meta, candidates_series


def assemble_output(anchor_tech: str, anchor_series: pd.Series,
                     anchor_xy: tuple[float, float],
                     candidates_meta: dict[str, list[dict]],
                     candidates_series: dict[str, list[pd.Series]],
                     max_radius_km: float, n_candidates: int,
                     threshold: float, w_coincidence: float,
                     w_correlation: float) -> tuple[pd.DataFrame, dict, dict]:
    """Assemble the multi-site CF table and its two metadata payloads.

    Returns (df, site_coords, anchor_meta):
    - df: DatetimeIndex `time`, columns `{tech}@anchor` / `{tech}@cNN`
      (hyphenated tech keys, candidates in score order).
    - site_coords: {column: {lat, lon}} for every column.
    - anchor_meta: diagnostics + config, written under the
      `anchor_colocation` metadata key.
    """
    anchor_x, anchor_y = anchor_xy
    anchor_col = f"{_hyphen(anchor_tech)}@anchor"

    columns: dict[str, pd.Series] = {anchor_col: anchor_series}
    site_coords: dict[str, dict[str, float]] = {
        anchor_col: {"lat": anchor_y, "lon": anchor_x}
    }

    candidates_out: dict[str, list[dict]] = {}
    n_found: dict[str, int] = {}
    for tech, series_list in candidates_series.items():
        tech_key = _hyphen(tech)
        meta_list = candidates_meta[tech]
        n_found[tech_key] = len(series_list)
        candidates_out[tech_key] = []
        for i, (series, meta) in enumerate(zip(series_list, meta_list)):
            col = f"{tech_key}@c{i:02d}"
            columns[col] = series
            site_coords[col] = {"lat": meta["y"], "lon": meta["x"]}
            candidates_out[tech_key].append({
                "column": col,
                "lat": meta["y"], "lon": meta["x"],
                "dist_km": meta["dist_km"],
                "score": meta["score"],
                "coincidence": meta["coincidence"],
                "correlation": meta["correlation"],
            })

    df = pd.DataFrame(columns)
    df.index.name = "time"

    anchor_meta = {
        "anchor_tech": _hyphen(anchor_tech),
        "anchor": {"lat": anchor_y, "lon": anchor_x},
        "candidates": candidates_out,
        "n_candidates_found": n_found,
        "config": {
            "max_radius_km": max_radius_km,
            "n_candidates": n_candidates,
            "coincidence_threshold": threshold,
            "w_coincidence": w_coincidence,
            "w_correlation": w_correlation,
        },
        "run_timestamp_utc": pd.Timestamp.now("UTC").isoformat(),
    }
    return df, site_coords, anchor_meta


def write_parquet(df: pd.DataFrame, site_coords: dict, anchor_meta: dict,
                   anchor_col: str, out_path: Path) -> None:
    """Write the CF table with site_coords / demand_site / anchor_colocation
    attached as Arrow schema metadata (self-describing multi-site contract)."""
    table = pa.Table.from_pandas(df, preserve_index=True)
    meta = dict(table.schema.metadata or {})
    meta[b"site_coords"] = json.dumps(site_coords).encode()
    meta[b"demand_site"] = anchor_col.encode()
    meta[b"anchor_colocation"] = json.dumps(anchor_meta).encode()
    table = table.replace_schema_metadata(meta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, str(out_path))


def main() -> None:
    country_upper = _CF_AREA.upper()
    anchor_tech = _ANCHOR_TECH.replace("-", "_")  # 07's functions use snake_case
    if anchor_tech not in TECHS:
        raise ValueError(f"Unknown anchor tech: {_ANCHOR_TECH!r}")

    anchor_series, anchor_xy, candidates_meta, candidates_series = build_anchor_scenario(
        _CUTOUT_PATH, country_upper, anchor_tech,
        MAX_RADIUS_KM, _N_CANDIDATES,
        COINCIDENCE_THRESHOLD, W_COINCIDENCE, W_CORRELATION,
        _MIN_LAND_FRACTION, _ELIGIBILITY_SOURCE,
    )

    df, site_coords, anchor_meta = assemble_output(
        anchor_tech, anchor_series, anchor_xy, candidates_meta, candidates_series,
        MAX_RADIUS_KM, _N_CANDIDATES,
        COINCIDENCE_THRESHOLD, W_COINCIDENCE, W_CORRELATION,
    )

    anchor_col = f"{_hyphen(anchor_tech)}@anchor"
    write_parquet(df, site_coords, anchor_meta, anchor_col, _OUT)
    log.info(
        f"{country_upper}: wrote {_OUT.name} "
        f"({len(df)} rows, {len(df.columns)} site columns, "
        f"anchor={_hyphen(anchor_tech)} at "
        f"lat={anchor_xy[1]:.2f} lon={anchor_xy[0]:.2f})"
    )


if __name__ == "__main__":
    main()
