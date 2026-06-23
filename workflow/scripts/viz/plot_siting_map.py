"""Geographic siting map for one multi-site h2_dri scenario.

The multi-site analogue of plot_cf_map: instead of a CF heatmap it reads a solved
network (results/{project}/{scenario}.nc) and shows, on a real map, where the LP
chose to build. Candidate sites are markers sized by built capacity and coloured
by tech; the demand/electrolyser site is starred; HVDC links are drawn from each
built site to the demand site with width proportional to optimal link capacity
(unbuilt candidate links are faint grey). Everything is derived from the network's
bus coordinates, so no cutout/region inputs are needed.

Snakemake rule: plot_siting_map (in viz.smk).
"""

import logging
from math import asin, cos, radians, sin, sqrt
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import pypsa

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.viz.style import (
    apply_header,
    blue_black,
    dark_gray,
    fca_blue,
    fca_template,
    light_gray,
    magenta_red,
    sand_yellow,
    save_figure,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

_BUILT_MW = 1.0  # capacity below this counts as "not built"
_TECH_COLOR = {"solar": sand_yellow, "wind-onshore": fca_blue, "wind-offshore": blue_black}


def _haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    dlon, dlat = radians(lon2 - lon1), radians(lat2 - lat1)
    a = sin(dlat / 2) ** 2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2) ** 2
    return 6371.0 * 2 * asin(sqrt(a))


def demand_site_bus(n: pypsa.Network) -> str:
    """The electricity bus hosting the electrolyser (= the demand/plant site)."""
    return n.links.at["electrolyser", "bus0"]


def build_site_table(n: pypsa.Network) -> pd.DataFrame:
    """One row per candidate electricity bus: coords, tech, built gen MW, link MW, CF, km.

    The demand bus (electrolyser bus0) is excluded — it carries no candidate
    generator. CF is the mean p_max_pu of the site's generator; km is the routed
    great-circle distance to the demand site (matching the network's link length).
    """
    demand = demand_site_bus(n)
    dx, dy = n.buses.at[demand, "x"], n.buses.at[demand, "y"]

    gen_by_bus = {row.bus: name for name, row in n.generators.iterrows()}
    link_by_src = {n.links.at[k, "bus0"]: k for k in n.links.index[n.links.carrier == "HVDC"]}

    rows = []
    for bus in n.buses.index[n.buses.carrier == "AC"]:
        if bus == demand:
            continue
        gen = gen_by_bus.get(bus)
        link = link_by_src.get(bus)
        lon, lat = n.buses.at[bus, "x"], n.buses.at[bus, "y"]
        rows.append({
            "site": bus.removeprefix("electricity_"),
            "lon": lon,
            "lat": lat,
            "tech": n.generators.at[gen, "carrier"] if gen else "?",
            "gen_mw": float(n.generators.at[gen, "p_nom_opt"]) if gen else 0.0,
            "link_mw": float(n.links.at[link, "p_nom_opt"]) if link else 0.0,
            "cf": float(n.generators_t.p_max_pu[gen].mean()) if gen else float("nan"),
            "km": _haversine_km(dx, dy, lon, lat),
        })
    return pd.DataFrame(rows).set_index("site")


def _marker_sizes(mw: pd.Series, mw_max: float) -> list[float]:
    """sqrt-scaled marker sizes; sub-threshold (unbuilt) sites get a small dot."""
    out = []
    for v in mw:
        out.append(8.0 + 34.0 * sqrt(v / mw_max) if v >= _BUILT_MW and mw_max > 0 else 5.0)
    return out


def plot(n: pypsa.Network, sites: pd.DataFrame, out: Path, project: str, scenario: str) -> None:
    """Render the siting map (markers + HVDC links + demand site) to PNG + HTML."""
    demand = demand_site_bus(n)
    dx, dy = n.buses.at[demand, "x"], n.buses.at[demand, "y"]
    gen_max = max(sites["gen_mw"].max(), 1.0)
    link_max = max(sites["link_mw"].max(), 1.0)

    fig = go.Figure()

    # HVDC links: faint grey for unbuilt candidate connections, coloured + width-
    # scaled for built ones (one trace each so width can vary; shared legend group).
    unbuilt_lon, unbuilt_lat = [], []
    built_shown = False
    for site, r in sites.iterrows():
        if r["link_mw"] < _BUILT_MW:
            unbuilt_lon += [r["lon"], dx, None]
            unbuilt_lat += [r["lat"], dy, None]
            continue
        fig.add_trace(go.Scattergeo(
            lon=[r["lon"], dx], lat=[r["lat"], dy], mode="lines",
            line=dict(width=1.5 + 6.0 * r["link_mw"] / link_max, color=dark_gray),
            opacity=0.75, legendgroup="hvdc", showlegend=not built_shown,
            name="HVDC link (built)",
            hovertemplate=f"{site} → plant<br>link {r['link_mw']:.0f} MW<extra></extra>",
        ))
        built_shown = True
    if unbuilt_lon:
        fig.add_trace(go.Scattergeo(
            lon=unbuilt_lon, lat=unbuilt_lat, mode="lines",
            line=dict(width=0.6, color=light_gray), opacity=0.5,
            name="candidate link (not built)", hoverinfo="skip",
        ))

    # Candidate site markers, one trace per tech (size ∝ built MW, faint if zero).
    for tech, grp in sites.groupby("tech"):
        fig.add_trace(go.Scattergeo(
            lon=grp["lon"], lat=grp["lat"], mode="markers",
            marker=dict(
                size=_marker_sizes(grp["gen_mw"], gen_max),
                color=_TECH_COLOR.get(tech, dark_gray),
                line=dict(width=1.0, color="white"),
                opacity=[0.95 if v >= _BUILT_MW else 0.45 for v in grp["gen_mw"]],
            ),
            name=tech,
            customdata=list(zip(grp.index, grp["gen_mw"], grp["cf"] * 100, grp["km"])),
            hovertemplate=("%{customdata[0]}<br>built %{customdata[1]:.0f} MW<br>"
                           "CF %{customdata[2]:.1f}%<br>%{customdata[3]:.0f} km to plant"
                           "<extra></extra>"),
        ))

    # Demand / electrolyser site.
    fig.add_trace(go.Scattergeo(
        lon=[dx], lat=[dy], mode="markers",
        marker=dict(symbol="star", size=18, color=magenta_red, line=dict(width=1.2, color="white")),
        name="demand site (electrolyser)",
        hovertemplate=f"demand site<br>lon {dx:.2f}, lat {dy:.2f}<extra></extra>",
    ))

    built = sites[sites["gen_mw"] >= _BUILT_MW]
    fig.update_layout(
        template=fca_template,
        geo=dict(scope="europe", resolution=50, fitbounds="locations",
                 showcountries=True, countrycolor=dark_gray,
                 showland=True, landcolor="rgb(243,243,243)",
                 showocean=True, oceancolor="rgb(228,241,247)",
                 lataxis=dict(showgrid=True, gridcolor=light_gray),
                 lonaxis=dict(showgrid=True, gridcolor=light_gray)),
        # Legend bottom-left so it clears the bottom-right brand logo.
        legend=dict(x=0.01, y=0.01, xanchor="left", yanchor="bottom",
                    bgcolor="rgba(255,255,255,0.75)"),
    )

    apply_header(
        fig,
        title=f"{project} / {scenario} — chosen siting",
        subtitle=(f"{len(built)} of {len(sites)} candidate sites built · "
                  "marker size ∝ capacity · line width ∝ HVDC capacity"),
        fig_width=900, fig_height=760,
        margin_l=20, margin_r=20, margin_t=90, margin_b=50,
    )
    saved = save_figure(fig, out.parent, out.stem)
    log.info(f"saved {' + '.join(saved)}")


def main() -> None:
    """Load the solved network and render its siting map."""
    project = snakemake.wildcards.project
    scenario = snakemake.wildcards.scenario
    n = pypsa.Network()
    n.import_from_netcdf(snakemake.input.network)
    if int((n.buses.carrier == "AC").sum()) <= 1:
        raise ValueError(
            f"{project}/{scenario} is a single-site network — plot_siting_map is "
            "only meaningful for multi-site scenarios (those with a sites overlay)."
        )
    sites = build_site_table(n)
    plot(n, sites, Path(snakemake.output.png), project, scenario)


if __name__ == "__main__":
    main()
