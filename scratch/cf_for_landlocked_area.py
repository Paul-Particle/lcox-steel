"""Build best-site P95 capacity factors for an area that has no coastline.

`d2_bestsite_p95` computes all three technologies on every call, whichever one the
rule asked for, and picks each one's P95 cell over the cells its geometry covers.
A landlocked area has no offshore geometry, so offshore wind has no eligible cell
and `pick_p95_cell` fails on an all-NaN distance array — taking the onshore and
solar outputs down with it, which are perfectly well defined.

This produces those two outputs directly, using the same functions and the same
one-column-per-tech contract, so the files are what the rule would have written.
It is a stopgap: the rule itself should skip a technology whose geometry covers no
cells and say so, which is the fix recorded in TODO.md. Doing it there would
change the script that every capacity-factor output depends on, and re-running
those would in turn invalidate every solved network, so it waits for a quiet moment.

Usage: cf_for_landlocked_area.py AREA REGION START END
"""
import importlib.util
import sys
from pathlib import Path

import atlite

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflow"))
sys.path.insert(0, str(REPO / "workflow/scripts"))

from res_cf._helpers_res_cf import eligibility_weights, load_res_cf_cfg, pick_p95_cell

area, region, start_date, end_date = sys.argv[1:5]

spec = importlib.util.spec_from_file_location(
    "bestsite", REPO / "workflow/scripts/res_cf/d2_bestsite_p95.py")
bestsite = importlib.util.module_from_spec(spec)
sys.modules["bestsite"] = bestsite
spec.loader.exec_module(bestsite)

bestsite.REGIONS_PATH = REPO / f"resources/shapes/{area}_geo.parquet"
bestsite.OFFSHORE_REGIONS_PATH = REPO / f"resources/shapes/{area}_offshore_geo.parquet"

config = load_res_cf_cfg()
min_land_fraction = float(config.get("min_land_fraction", 0.0))
eligibility_source = config.get("eligibility_source", "indicatormatrix")

cutout_path = REPO / f"cutouts/{area}_{start_date}_{end_date}.nc"
cutout = atlite.Cutout(path=str(cutout_path))

for tech in ["wind_onshore", "solar"]:
    cf_year = bestsite.build_cf_year(cutout_path, tech)
    geometry = bestsite.geometry_for_tech(region, tech)
    weights = eligibility_weights(cutout, geometry, min_land_fraction, eligibility_source)
    y_idx, x_idx = pick_p95_cell(cf_year.mean("time").values, weights)
    series = bestsite.extract_cell_timeseries(cf_year, y_idx, x_idx)

    column = tech.replace("_", "-")
    out_path = (REPO / "resources/timeseries"
                / f"{area}_{column}_bestsite-p95_{start_date}_{end_date}.parquet")
    series.rename(column).to_frame().to_parquet(out_path, index=True)
    print(f"{area} {column}: best-site mean CF {float(series.mean()):.3f} "
          f"over {int((weights > 0).sum())} eligible cells -> {out_path.name}")
