"""Append two one-input sweeps to the scenario table, with their overlays.

Both headline orderings in the September run turn on a single input sitting near
one end of its own stated range, so each sweep moves that one input across the
range and reports where the ordering changes.

  ew-capex-*   electrowinning capex, 500-3000 EUR/(t/yr) — the span the
               assumptions file already names. MOE retakes first place somewhere
               between 852 and 1110 on the grid cases, which is below the value
               the model used a day earlier.
  gas-*        natural gas, 30-70 EUR/MWh against the run's 21.9. The fossil
               advantage disappears between 45 and 67 depending on geography.

Three areas each, chosen to bracket the predicted crossings rather than to cover
the map: the area cell is a scenario-level set, so a scenario can name at most as
many areas as it has rows, and three rows is what a grid scenario takes.
"""
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WINDOW = ("20250101", "20251231")

# (scenario prefix, overlay dict template, routes, values, areas)
SWEEPS = [
    ("ew-capex", "ew", "capex_per_t_per_year_eur", "ew-eaf|moe-eaf",
     [500, 900, 1200, 1800, 3000],
     # break-even predicted at ESP 852, DEU 966, SA1 1110 — one below, one inside,
     # one above the value the model used before this run.
     ["ESP", "DEU", "SA1"]),
    ("gas", "natural_gas", "price_eur_per_mwh", "ng-dri-eaf|h2-dri-eaf",
     [30, 40, 50, 60, 70],
     # predicted sign change at SA1 45.8, DEU 63.4, ESP 65.7.
     ["SA1", "DEU", "ESP"]),
]

lines = []
for prefix, block, key, routes, values, areas in SWEEPS:
    for value in values:
        scenario = f"{prefix}-{value}"
        overlay = REPO / f"config/assumptions_{scenario}.yaml"
        overlay.write_text(
            f"# One input moved, everything else as the base run. See the\n"
            f"# `{prefix}` sweep in scratch/add_sweep_scenarios.py for why this range.\n"
            f"{block}:\n  {key}: {value}\n"
        )
        # Grid rows: the crossings being looked for are all in grid cases.
        techs = [("wind-onshore", "bestsite-p95"), ("solar", "bestsite-p95"),
                 ("grid", "emissions")]
        for (tech, variant), area in zip(techs, areas):
            lines.append(f"{scenario},{routes},{tech},{variant},{area},{WINDOW[0]},{WINDOW[1]}")

table = REPO / "config/scenarios.csv"
text = table.read_text()
marker = "# September 2026 test run"
head, _, tail = text.partition(marker)
table.write_text(head + "\n".join(lines) + "\n" + marker + tail)
print(f"added {len(lines)} rows across {sum(len(s[4]) for s in SWEEPS)} scenarios")
