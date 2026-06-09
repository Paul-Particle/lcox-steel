"""Compile per-scenario summaries into a single project-level report CSV.

Invoked by Snakemake's `script:` directive (compile_report rule in viz.smk).
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from common._constants import H2_LHV_KWH_PER_KG
from common._logging import configure_logging

if "snakemake" not in globals():
    from common._stubs import snakemake

configure_logging(snakemake)
log = logging.getLogger(__name__)


def _h2_produced_kg(n: pypsa.Network) -> float:
    """Annual H2 produced by the electrolyser, in kg, scaled to 8760 h.

    Read from the electrolyser link's H2-side output rather than the dri_load,
    so the result reflects actual production. The two coincide when the model
    is feasible and the H2 buffer is cyclic, but only the link side stays
    correct if the load formulation changes later (e.g. flexible demand).
    """
    t_hours = len(n.snapshots)
    # PyPSA Link sign convention: p1 < 0 when the link injects power into bus1,
    # so -p1 is the (positive) H2 LHV output on the hydrogen bus.
    h2_mwh_lhv = -float(n.links_t.p1["electrolyser"].sum()) * (8760.0 / t_hours)
    return h2_mwh_lhv / (H2_LHV_KWH_PER_KG / 1000.0)


def _annual_cost(n: pypsa.Network) -> float:
    """Annualised capital costs + variable grid-import costs (scaled to 8760 h).

    Capital costs are already per-year (annualised CAPEX × p_nom_opt). Grid-import
    variable costs are scaled from the simulation period up to 8760 h so LCOH stays
    meaningful on partial-year runs.
    """
    t_hours = len(n.snapshots)
    annual_scale = 8760.0 / t_hours

    cost = 0.0

    mask = n.generators.p_nom_extendable
    cost += (n.generators.loc[mask, "capital_cost"] * n.generators.loc[mask, "p_nom_opt"]).sum()

    mask = n.storage_units.p_nom_extendable
    cost += (n.storage_units.loc[mask, "capital_cost"] * n.storage_units.loc[mask, "p_nom_opt"]).sum()

    # All extendable links: the electrolyser (always) and, in multi-site runs, the
    # HVDC transmission links. Summing the mask captures both without double count.
    mask = n.links.p_nom_extendable
    cost += (n.links.loc[mask, "capital_cost"] * n.links.loc[mask, "p_nom_opt"]).sum()

    mask = n.stores.e_nom_extendable
    cost += (n.stores.loc[mask, "capital_cost"] * n.stores.loc[mask, "e_nom_opt"]).sum()

    if "grid_import" in n.generators.index:
        p = n.generators_t.p.get("grid_import", pd.Series(0.0, index=n.snapshots))
        if "grid_import" in n.generators_t.marginal_cost.columns:
            mc = n.generators_t.marginal_cost["grid_import"]
        else:
            mc = n.generators.at["grid_import", "marginal_cost"]
        cost += float((p * mc).sum()) * annual_scale

    return float(cost)


def _compute_lcoh(n: pypsa.Network) -> float:
    return _annual_cost(n) / _h2_produced_kg(n)


def extract_summary(n: pypsa.Network, project_name: str, scenario_name: str) -> dict:
    """Key sizing and cost metrics as a flat dict (suitable for a one-row CSV)."""
    lcoh_eur_per_kg = _compute_lcoh(n)
    summary = {
        "project": project_name,
        "scenario": scenario_name,
        "lcoh_eur_per_kg": lcoh_eur_per_kg,
        "lcoh_eur_per_mwh_lhv": lcoh_eur_per_kg * 1000.0 / H2_LHV_KWH_PER_KG,
        "total_annual_cost_meur": _annual_cost(n) / 1e6,
        "h2_produced_kt": _h2_produced_kg(n) / 1e6,
    }

    for gen in n.generators.index[n.generators.p_nom_extendable]:
        summary[f"{gen}_gw_opt"] = n.generators.at[gen, "p_nom_opt"] / 1e3

    if "battery" in n.storage_units.index:
        p_opt = n.storage_units.at["battery", "p_nom_opt"]
        summary["battery_gw_opt"] = p_opt / 1e3
        summary["battery_mwh_opt"] = p_opt * n.storage_units.at["battery", "max_hours"]

    if "h2_buffer" in n.stores.index:
        # Reported as hours of DRI H₂ demand (buffer MWh LHV / load MW LHV).
        buffer_mwh = n.stores.at["h2_buffer", "e_nom_opt"]
        if "dri_load" in n.loads.index:
            dri_mw = float(n.loads.at["dri_load", "p_set"])
            summary["h2_buffer_hours_dri"] = buffer_mwh / dri_mw if dri_mw else float("nan")
        else:
            summary["h2_buffer_hours_dri"] = float("nan")

    if "electrolyser" in n.links.index:
        el_cap = n.links.at["electrolyser", "p_nom_opt"]
        summary["electrolyser_gw"] = el_cap / 1e3
        if el_cap > 0 and "electrolyser" in n.links_t.p0.columns:
            summary["electrolyser_utilization"] = float(
                n.links_t.p0["electrolyser"].mean() / el_cap
            )
        else:
            summary["electrolyser_utilization"] = float("nan")

    if "dri_load" in n.loads.index:
        summary["dri_h2_mw_lhv"] = float(n.loads.at["dri_load", "p_set"])

    # Multi-site only (guarded so single-site reports are unchanged): one column
    # per HVDC link capacity, the total annualised transmission cost, and per-tech
    # built-capacity totals summed across candidate sites.
    hvdc = n.links.index[n.links.carrier == "HVDC"]
    if len(hvdc):
        trans_cost = 0.0
        for link in hvdc:
            cap = n.links.at[link, "p_nom_opt"]
            summary[f"{link}_gw_opt"] = cap / 1e3
            trans_cost += n.links.at[link, "capital_cost"] * cap
        summary["transmission_total_annual_cost_meur"] = trans_cost / 1e6

    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    if len(ac_buses) > 1:
        ext = n.generators[n.generators.p_nom_extendable]
        for carrier, grp in ext.groupby("carrier"):
            summary[f"{carrier}_total_gw_opt"] = grp["p_nom_opt"].sum() / 1e3

    return summary


def main() -> None:
    """Load each scenario network for the project and write the combined report CSV.

    Dedupes the netCDF inputs (collect fans out per tech row), extracts one summary
    row per scenario via `extract_summary`, rounds numerics, and writes
    results/report_<project>.csv.
    """
    project_name = snakemake.wildcards.project

    rows = []
    # networks may contain duplicates (collect fans out per tech row); dedupe
    # while preserving order so each scenario appears once.
    network_paths = list(dict.fromkeys(snakemake.input.networks))
    log.info(f"compiling report for project={project_name} ({len(network_paths)} scenarios)")
    for nc_path in network_paths:
        nc_path = Path(nc_path)
        scenario_name = nc_path.stem
        n = pypsa.Network()
        n.import_from_netcdf(nc_path)
        rows.append(extract_summary(n, project_name, scenario_name))

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df[df.select_dtypes("number").columns] = df.select_dtypes("number").round(2)
    df.to_csv(out_path, index=False)
    log.info(f"wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
