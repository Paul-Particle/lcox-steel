"""Build per-cell candidate CF profiles for multi-site plant-location optimization.

Where `build_res_cf_profile.py` collapses a country to one CF series via the
indicator matrix, this script keeps the cell dimension: it sub-samples cutout
cells inside the region on a coarse regular lattice (targeting ~`n_cells`), and
writes one CF column per surviving cell. The downstream h2_dri network builder
turns each column into its own electricity bus + extendable generator, and the
LP chooses how much to build where (siting emerges from the optimization).

Output is a single multi-column parquet:
  - columns named `{tech}@c{ii}` (e.g. `solar@c00`) — the part before `@` is the
    bare tech key looked up in assumptions, the part after is the site id;
  - each column's `{lat, lon}` is stored in the parquet's Arrow schema metadata
    under `site_coords` (JSON), so the file is self-describing and the network
    builder needs no separate coordinates input.

Cell extraction reuses the per-cell-grid pattern from build_solar_tilt_mix_p95.py
and plot_cf_map.py (`capacity_factor_timeseries=True`), not the indicatormatrix
aggregation.
"""

import json
import logging
from pathlib import Path

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from common._paths import CUTOUTS, RES_CF, SHAPES_RES

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging

configure_logging(snakemake)
log = logging.getLogger(__name__)

# Standalone defaults
_CF_AREA = "de"
_TECH = "wind-onshore"
_N_CELLS = 10
_START_DATE = "20230101"
_END_DATE = "20231231"
_CUTOUT_PATH = CUTOUTS / "de_20230101_20231231.nc"
_REGIONS_PATH = SHAPES_RES / "de_geo.parquet"
_OFFSHORE_REGIONS_PATH = SHAPES_RES / "de_offshore_geo.parquet"
_OUT = RES_CF / "de_wind-onshore_grid-n10_20230101_20231231.parquet"
_REGION = "DE"
_WIND_ONSHORE_TURBINE = "Vestas_V112_3MW"
_WIND_OFFSHORE_TURBINE = "NREL_ReferenceTurbine_5MW_offshore"
_PV_PANEL = "CSi"
_PV_ORIENTATION = "latitude_optimal"
_WIND_CF = {"smooth": True, "add_cutout_windspeed": True}

if "snakemake" in globals() and hasattr(snakemake, "wildcards"):
    _CF_AREA = snakemake.wildcards.cf_area
    _TECH = snakemake.wildcards.tech
    _N_CELLS = int(snakemake.wildcards.n_cells)
    _START_DATE = snakemake.wildcards.start_date
    _END_DATE = snakemake.wildcards.end_date
    _CUTOUT_PATH = Path(snakemake.input.cutout)
    _REGIONS_PATH = Path(snakemake.input.regions)
    _OFFSHORE_REGIONS_PATH = Path(snakemake.input.offshore_regions)
    _OUT = Path(snakemake.output[0])
    _REGION = snakemake.params.region
    _WIND_ONSHORE_TURBINE = snakemake.params.wind_onshore_turbine
    _WIND_OFFSHORE_TURBINE = snakemake.params.wind_offshore_turbine
    _PV_PANEL = snakemake.params.pv_panel
    _PV_ORIENTATION = snakemake.params.pv_orientation
    _WIND_CF = snakemake.params.wind_cf


_WIND_SMOOTH = _WIND_CF.get("smooth", True)
_WIND_ADD_CUTOUT_WS = _WIND_CF.get("add_cutout_windspeed", True)


def _get_region_geometry(path: Path, region: str):
    """Return the geometry of `region` from the GeoParquet at `path` (raises if missing)."""
    gdf = gpd.read_parquet(path).to_crs(4326)
    row = gdf.loc[gdf["region"] == region]
    if row.empty:
        raise ValueError(f"Region '{region}' not found in {path}")
    return row.geometry.iloc[0]


def _mask_cells_inside(cell_mean, geom) -> np.ndarray:
    """Return a boolean grid marking cutout cells whose centre lies within `geom`.

    Copied from build_solar_tilt_mix_p95.py to keep this script self-contained.
    """
    xs = cell_mean.coords["x"].values
    ys = cell_mean.coords["y"].values
    xx, yy = np.meshgrid(xs, ys)
    points = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel()), crs=4326)
    inside = points.within(geom) | points.touches(geom)
    return inside.values.reshape(cell_mean.shape)


def _select_lattice_cells(
    cell_mean, inside: np.ndarray, n_cells: int
) -> list[tuple[int, int]]:
    """Pick ~`n_cells` in-region cells on a coarse regular lattice.

    Strides the (y, x) index grid so the regular sublattice intersected with the
    region yields about `n_cells` cells; if that overshoots, keeps the highest
    annual-mean-CF cells (restoring row-major order for stable column naming).
    Deterministic — no randomness.
    """
    yy, xx = np.where(inside)
    n_inside = len(yy)
    if n_inside == 0:
        raise ValueError("no cutout cells fall inside the region")
    if n_inside <= n_cells:
        return sorted(zip(yy.tolist(), xx.tolist()))

    stride = max(1, int(round((n_inside / n_cells) ** 0.5)))
    keep = (yy % stride == 0) & (xx % stride == 0)
    chosen = list(zip(yy[keep].tolist(), xx[keep].tolist()))

    # The lattice can over- or under-shoot n_cells; reconcile using annual-mean CF
    # as the tie-breaker so we always return ~n_cells cells (never empty).
    if len(chosen) > n_cells:
        cf_vals = np.array([cell_mean.values[y, x] for y, x in chosen])
        chosen = [chosen[i] for i in np.argsort(-cf_vals)[:n_cells]]
    elif len(chosen) < n_cells:
        chosen_set = set(chosen)
        for k in np.argsort(-cell_mean.values[yy, xx]):
            c = (int(yy[k]), int(xx[k]))
            if c not in chosen_set:
                chosen.append(c)
                chosen_set.add(c)
                if len(chosen) >= n_cells:
                    break

    return sorted(chosen)  # row-major order → stable c00, c01, … naming


def _cf_grid(cutout):
    """Return the per-cell CF time-series grid (dims time, y, x) for the tech."""
    if _TECH == "solar":
        return cutout.pv(
            panel=_PV_PANEL,
            orientation=_PV_ORIENTATION,
            capacity_factor_timeseries=True,
        )
    if _TECH == "wind-onshore":
        return cutout.wind(
            turbine=_WIND_ONSHORE_TURBINE,
            capacity_factor_timeseries=True,
            smooth=_WIND_SMOOTH,
            add_cutout_windspeed=_WIND_ADD_CUTOUT_WS,
        )
    if _TECH == "wind-offshore":
        return cutout.wind(
            turbine=_WIND_OFFSHORE_TURBINE,
            capacity_factor_timeseries=True,
            smooth=_WIND_SMOOTH,
            add_cutout_windspeed=_WIND_ADD_CUTOUT_WS,
        )
    raise ValueError(f"Unknown tech: {_TECH!r}")


def main() -> None:
    """Sub-sample in-region cells on a regular lattice and write per-cell CF columns."""
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    cutout = atlite.Cutout(str(_CUTOUT_PATH))

    regions_path = (
        _OFFSHORE_REGIONS_PATH if _TECH == "wind-offshore" else _REGIONS_PATH
    )
    geom = _get_region_geometry(regions_path, _REGION)

    log.info(f"computing {_TECH} CF grid for {_CF_AREA}")
    cf_grid = _cf_grid(cutout)
    cell_mean = cf_grid.mean("time")

    inside = _mask_cells_inside(cell_mean, geom)
    cells = _select_lattice_cells(cell_mean, inside, _N_CELLS)
    log.info(f"selected {len(cells)} candidate cells (target {_N_CELLS})")

    results: dict[str, pd.Series] = {}
    site_coords: dict[str, dict[str, float]] = {}
    for ii, (y_idx, x_idx) in enumerate(cells):
        ts = cf_grid.isel(y=y_idx, x=x_idx).to_pandas()
        ts.index = pd.to_datetime(ts.index)
        ts = ts.clip(0.0, 1.0)
        lat = float(cutout.data.coords["y"].isel(y=y_idx))
        lon = float(cutout.data.coords["x"].isel(x=x_idx))
        col = f"{_TECH}@c{ii:02d}"
        results[col] = ts
        site_coords[col] = {"lat": lat, "lon": lon}
        log.info(
            f"  {col}: lat={lat:.2f} lon={lon:.2f} annual_mean_cf={float(ts.mean()):.3f}"
        )

    df = pd.DataFrame(results)
    df.index.name = "time"

    # Attach per-column coordinates to the parquet schema metadata so the network
    # builder reads them straight off the file (no separate coords input).
    table = pa.Table.from_pandas(df, preserve_index=True)
    meta = dict(table.schema.metadata or {})
    meta[b"site_coords"] = json.dumps(site_coords).encode()
    table = table.replace_schema_metadata(meta)
    pq.write_table(table, str(_OUT))
    log.info(f"wrote {_OUT} ({len(df)} rows, {len(results)} candidate cells)")


if __name__ == "__main__":
    main()
