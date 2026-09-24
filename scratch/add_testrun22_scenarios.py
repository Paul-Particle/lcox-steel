"""Append the testrun 2.2 flexibility set to the scenario table, with its overlays.

Three questions, all Spain, 2025, best-site wind and solar, islanded and grid:

  esp-base, esp-ramp-*   the ramp set of testrun 2.1 re-solved on this code, cavern
                         storage, so both halves of the storage comparison share one
                         code_hash.
  esp-tank-*             the same ramp cases on tank-priced hydrogen storage, where
                         the shaft can no longer lean on a months-deep buffer.
  esp-*dri-pmin-0*       the shaft's minimum load and its ramp limit apart and
                         together, at both storage prices.
  esp-moe-*              the MOE cell's minimum load at 0.7 (one pot's figure) and 0
                         (pots switched off), each with no ramp limit, 20 %/h and 5 %/h.

A scenario already in the table is left alone, overlay included, so the testrun 2.1
ramp cases keep their own files and re-running this script adds nothing twice.
"""
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
WINDOW = ("20250101", "20251231")
TANK_EUR_PER_MWH = 60050

ramp = lambda rate: {"ramp_limit_up": rate, "ramp_limit_down": rate}
tank = {"h2_buffer": {"capex_per_mwh_eur": TANK_EUR_PER_MWH}}

# (scenario, route, overlay, one-line description); None means no overlay file.
SCENARIOS = [
    ("esp-base", "h2-dri-eaf", None, ""),
    ("esp-ramp-20", "h2-dri-eaf",
     {"electrolyser": ramp(0.20), "dri-h2": ramp(0.20)},
     "Electrolyser and H2 shaft each held to 0.20 of nameplate per hour."),
    ("esp-ramp-05", "h2-dri-eaf",
     {"electrolyser": ramp(0.05), "dri-h2": ramp(0.05)},
     "Electrolyser and H2 shaft each held to 0.05 of nameplate per hour."),
    ("esp-ramp-05-elyzr", "h2-dri-eaf", {"electrolyser": ramp(0.05)},
     "Electrolyser alone held to 0.05 of nameplate per hour."),
    ("esp-ramp-05-dri", "h2-dri-eaf", {"dri-h2": ramp(0.05)},
     "H2 shaft alone held to 0.05 of nameplate per hour."),
    ("esp-tank", "h2-dri-eaf", tank,
     "Hydrogen storage priced as tank with compressor, no ramp limit."),
    ("esp-tank-ramp-20", "h2-dri-eaf",
     {**tank, "electrolyser": ramp(0.20), "dri-h2": ramp(0.20)},
     "Tank storage; electrolyser and H2 shaft each held to 0.20 of nameplate per hour."),
    ("esp-tank-ramp-05", "h2-dri-eaf",
     {**tank, "electrolyser": ramp(0.05), "dri-h2": ramp(0.05)},
     "Tank storage; electrolyser and H2 shaft each held to 0.05 of nameplate per hour."),
    ("esp-tank-ramp-05-elyzr", "h2-dri-eaf", {**tank, "electrolyser": ramp(0.05)},
     "Tank storage; electrolyser alone held to 0.05 of nameplate per hour."),
    ("esp-tank-ramp-05-dri", "h2-dri-eaf", {**tank, "dri-h2": ramp(0.05)},
     "Tank storage; H2 shaft alone held to 0.05 of nameplate per hour."),
    ("esp-dri-pmin-0", "h2-dri-eaf", {"dri-h2": {"p_min_pu": 0.0}},
     "H2 shaft with no minimum load and no ramp limit."),
    ("esp-dri-pmin-0-ramp-05", "h2-dri-eaf", {"dri-h2": {"p_min_pu": 0.0, **ramp(0.05)}},
     "H2 shaft with no minimum load, held to 0.05 of nameplate per hour."),
    ("esp-tank-dri-pmin-0", "h2-dri-eaf", {**tank, "dri-h2": {"p_min_pu": 0.0}},
     "Tank storage; H2 shaft with no minimum load and no ramp limit."),
    ("esp-tank-dri-pmin-0-ramp-05", "h2-dri-eaf",
     {**tank, "dri-h2": {"p_min_pu": 0.0, **ramp(0.05)}},
     "Tank storage; H2 shaft with no minimum load, held to 0.05 of nameplate per hour."),
]
for pmin in (0.7, 0.0):
    label = f"{round(pmin * 100):02d}" if pmin else "0"
    for rate in (None, 0.20, 0.05):
        suffix = f"-ramp-{round(rate * 100):02d}" if rate else ""
        block = {"p_min_pu": pmin, **(ramp(rate) if rate else {})}
        limit = f"held to {rate:.2f} of nameplate per hour" if rate else "no ramp limit"
        SCENARIOS.append((
            f"esp-moe-pmin-{label}{suffix}", "moe-eaf", {"moe": block},
            f"MOE cell with minimum load {pmin}, {limit}.",
        ))

islanded_techs = [("wind-onshore", "bestsite-p95"), ("solar", "bestsite-p95")]
grid_techs = islanded_techs + [("grid", "emissions")]

table = REPO / "config/scenarios.csv"
existing = {line.split(",")[0] for line in table.read_text().splitlines()}
lines = []
for base_name, route, overlay, description in SCENARIOS:
    for scenario, techs in ((base_name, islanded_techs), (f"{base_name}-grid", grid_techs)):
        if scenario in existing:
            continue
        if overlay is not None:
            header = (
                f"# Overlay for the `{scenario}` scenario, deep-merged over assumptions.yaml.\n"
                f"#\n# {description}\n#\n"
                f"# Part of testrun 2.2; see scratch/add_testrun22_scenarios.py.\n"
            )
            body = yaml.safe_dump(overlay, sort_keys=False)
            (REPO / f"config/overlays/{scenario}.yaml").write_text(header + body)
        lines += [f"{scenario},{route},{tech},{variant},ESP,{WINDOW[0]},{WINDOW[1]}"
                  for tech, variant in techs]

marker = "# Testrun 2.2: hydrogen storage price, minimum load and ramp limits, ESP only."
text = table.read_text().rstrip("\n") + "\n"
if marker not in text:
    text += marker + "\n"
table.write_text(text + "".join(line + "\n" for line in lines))
print(f"added {len(lines)} rows; {len(SCENARIOS) * 2} scenarios in the set")
