"""Per-site built capacity + transmission for one multi-site h2_dri scenario.

The multi-site analogue of plot_capacity_bars: where that chart puts scenarios on
the x-axis, this puts the candidate *sites* of a single multi-site scenario on the
x-axis (only those the LP actually built), with grouped bars for generation MW
(coloured by tech) and the site's HVDC link MW. Read straight from the solved
network (results/{project}/{scenario}.nc).

Snakemake rule: plot_site_capacity_bars (in viz.smk).
"""

import logging
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import pypsa

if "snakemake" not in globals():
    from common._stubs import snakemake

from common._logging import configure_logging
from scripts.viz.style import (
    apply_header,
    dark_gray,
    fca_blue,
    fca_template,
    sand_yellow,
    save_figure,
)

configure_logging(snakemake)
log = logging.getLogger(__name__)

_BUILT_MW = 1.0
_WIND = ("wind-onshore", "wind-offshore")


def build_site_table(n: pypsa.Network) -> pd.DataFrame:
    """One row per candidate site (excludes the demand bus): tech, gen MW, link MW."""
    demand = n.links.at["electrolyser", "bus0"]
    gen_by_bus = {row.bus: name for name, row in n.generators.iterrows()}
    link_by_src = {n.links.at[k, "bus0"]: k for k in n.links.index[n.links.carrier == "HVDC"]}

    rows = []
    for bus in n.buses.index[n.buses.carrier == "AC"]:
        if bus == demand:
            continue
        gen = gen_by_bus.get(bus)
        link = link_by_src.get(bus)
        rows.append({
            "site": bus.removeprefix("electricity_"),
            "tech": n.generators.at[gen, "carrier"] if gen else "?",
            "gen_mw": float(n.generators.at[gen, "p_nom_opt"]) if gen else 0.0,
            "link_mw": float(n.links.at[link, "p_nom_opt"]) if link else 0.0,
        })
    return pd.DataFrame(rows).set_index("site")


def plot(sites: pd.DataFrame, out: Path, project: str, scenario: str) -> None:
    """Grouped bars (solar gen / wind gen / HVDC link MW) per built site → PNG + HTML."""
    built = sites[(sites["gen_mw"] >= _BUILT_MW) | (sites["link_mw"] >= _BUILT_MW)]
    built = built.sort_values("gen_mw", ascending=False)
    x = list(built.index)

    solar_y = [r["gen_mw"] if r["tech"] == "solar" else 0.0 for _, r in built.iterrows()]
    wind_y  = [r["gen_mw"] if r["tech"] in _WIND else 0.0 for _, r in built.iterrows()]
    link_y  = list(built["link_mw"])

    fig = go.Figure()
    fig.add_trace(go.Bar(x=x, y=solar_y, name="Solar (MW)", marker_color=sand_yellow,
                         hovertemplate="%{x}<br>solar %{y:.0f} MW<extra></extra>"))
    fig.add_trace(go.Bar(x=x, y=wind_y, name="Wind (MW)", marker_color=fca_blue,
                         hovertemplate="%{x}<br>wind %{y:.0f} MW<extra></extra>"))
    fig.add_trace(go.Bar(x=x, y=link_y, name="HVDC link (MW)", marker_color=dark_gray,
                         hovertemplate="%{x}<br>link %{y:.0f} MW<extra></extra>"))

    fig.update_layout(
        template=fca_template,
        barmode="group", bargap=0.25, bargroupgap=0.08,
        # No rotated y-axis title — the unit (MW) lives in the subtitle.
        yaxis_title=None, xaxis_title=None,
        legend=dict(x=0.99, y=0.99, xanchor="right", yanchor="top",
                    bgcolor="rgba(255,255,255,0.65)"),
    )
    fig.update_xaxes(type="category", tickangle=-45)
    fig.update_yaxes(rangemode="tozero")

    apply_header(
        fig,
        title=f"{project} / {scenario} — built capacity by site",
        subtitle=f"MW · {len(built)} sites built · generation (by tech) and HVDC link capacity",
        fig_width=max(720, 90 * len(x) + 320), fig_height=560,
        margin_l=80, margin_r=40, margin_t=100, margin_b=120,
    )
    saved = save_figure(fig, out.parent, out.stem)
    log.info(f"saved {' + '.join(saved)}")


def main() -> None:
    """Load the solved network and render its per-site capacity bars."""
    project = snakemake.wildcards.project
    scenario = snakemake.wildcards.scenario
    n = pypsa.Network()
    n.import_from_netcdf(snakemake.input.network)
    if int((n.buses.carrier == "AC").sum()) <= 1:
        raise ValueError(
            f"{project}/{scenario} is a single-site network — plot_site_capacity_bars "
            "is only meaningful for multi-site scenarios (those with a sites overlay)."
        )
    plot(build_site_table(n), Path(snakemake.output.png), project, scenario)


if __name__ == "__main__":
    main()
