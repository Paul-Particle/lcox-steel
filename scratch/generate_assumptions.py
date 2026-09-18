#!/usr/bin/env python3
"""Generate config/assumptions.yaml from the clean assumptions input workbook.

Expected workbook
-----------------
File: data/assumptions/Assumptions_yaml_inputs.xlsx
Sheet: assumptions_inputs
Columns:
  A  YAML key          required; dotted path used in assumptions.yaml
  B  Value             required; written to YAML as-is (no unit conversion)
  C  Unit              documentation only
  D  Old YAML key      documentation only
  E  Confirmed row ID  provenance only

Only columns A and B drive the generated YAML. Columns C-E may be retained for
human documentation and code-rewiring/provenance work.

Run from the repository root:
    python scratch/generate_assumptions.py

Optional custom paths:
    python scratch/generate_assumptions.py INPUT.xlsx OUTPUT.yaml
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import openpyxl
import yaml

SHEET = "assumptions_inputs"
DEFAULT_INPUT = Path("data/assumptions/Assumptions_yaml_inputs.xlsx")
DEFAULT_OUTPUT = Path("config/assumptions.yaml")
EXPECTED_HEADERS = ["YAML key", "Value", "Unit", "Old YAML key", "Confirmed row ID"]

# Inputs deliberately represented in the model interface but not yet researched.
# They are not rows in the clean Excel input and are therefore added as null.
NULL_KEYS = (
    "dri.p_min_pu",
    "eaf.p_min_pu",
    "moe.p_min_pu",
    "ewin.p_min_pu",
    "steel_store.max_weeks",
    "grid.connection_capex_eur_per_mw",
    "grid.connection_lifetime_years",
    "grid.fee_eur_per_mw_per_year",
    "grid.fee_eur_per_mwh",
    "transmission.cost_per_mw_per_km_eur",
    "transmission.lifetime_years",
    "transmission.losses_pct_per_1000km",
    "transmission.indirect_route_factor",
)

# Readability only; YAML semantics do not depend on top-level ordering.
TOP_LEVEL_ORDER = (
    "finance",
    "res", "electrolyser", "battery", "h2_buffer",
    "IS",
    # Ironmaking and iron handling together.
    "dri", "dri_h2", "dri_ng", "moe", "ewin", "iron", "briquetting",
    # Steelmaking and storage.
    "eaf", "eaf_dri", "eaf_hbi", "eaf_ewin", "eaf_moe",
    "iron_steel_store", "steel_store",
    # External inputs / infrastructure.
    "natural_gas", "grid", "transmission", "shipping",
)

# Presentation rounding agreed for assumptions.yaml. The Excel remains the
# source of the unrounded values.
ROUND_2_PREFIXES = (
    "dri_h2.",
    "battery.capex_",
    "dri_ng.",
    "eaf.fom_",
    "eaf_dri.el_mwh_per_t",
    "eaf_hbi.el_mwh_per_t",
    "eaf_ewin.el_mwh_per_t",
    "eaf_moe.el_mwh_per_t",
)
ROUND_2_KEYS = {
    "electrolyser.efficiency_kwh_per_kg",
    "electrolyser.varopex_usd_per_mwh_h2_base",
    "electrolyser.varopex_usd_per_mwh_h2_low",
    "electrolyser.varopex_usd_per_mwh_h2_high",
    "h2_buffer.fom_usd_per_mwh",
    "h2_buffer.capex_usd_per_mwh",
    "dri.ore_usd_per_t",
    "moe.el_mwh_per_t",
    "moe.ore_usd_per_t_liquid_iron",
    "ewin.el_mwh_per_t",
    "ewin.ore_usd_per_t_electrolytic_iron",
    "natural_gas.price_usd_per_mwh_high",
    "natural_gas.price_usd_per_mwh_low",
    "iron_steel_store.capex_usd_per_t",
}


def clean_scalar(value: Any) -> Any:
    """Convert Excel scalar types to clean YAML-safe Python values."""
    if value is None:
        return None
    if isinstance(value, float):
        if math.isnan(value):
            return None
        if value.is_integer():
            return int(value)
    return value


def value_for_yaml(key: str, value: Any) -> Any:
    value = clean_scalar(value)
    if isinstance(value, (int, float)) and (
        key in ROUND_2_KEYS or key.startswith(ROUND_2_PREFIXES)
    ):
        return round(value, 2)
    return value


def nested_set(root: dict[str, Any], dotted_key: str, value: Any) -> None:
    """Insert a value using a dotted YAML path, e.g. res.solar.capex..."""
    parts = dotted_key.split(".")
    node = root
    for part in parts[:-1]:
        existing = node.get(part)
        if existing is None:
            node[part] = {}
        elif not isinstance(existing, dict):
            raise ValueError(
                f"YAML hierarchy collision at {part!r} while inserting {dotted_key!r}"
            )
        node = node[part]

    leaf = parts[-1]
    if leaf in node:
        raise ValueError(f"Duplicate YAML key: {dotted_key}")
    node[leaf] = value


def reorder(data: dict[str, Any]) -> dict[str, Any]:
    ordered = {key: data[key] for key in TOP_LEVEL_ORDER if key in data}
    ordered.update({key: value for key, value in data.items() if key not in ordered})

    # financing_period applies to all renewable technologies, so make that
    # visually explicit by keeping it first under res:.
    res = ordered.get("res")
    if isinstance(res, dict) and "financing_period" in res:
        ordered["res"] = {
            "financing_period": res["financing_period"],
            **{key: value for key, value in res.items() if key != "financing_period"},
        }
    return ordered


def build_assumptions(xlsx_path: Path) -> dict[str, Any]:
    wb = openpyxl.load_workbook(xlsx_path, data_only=True, read_only=True)
    if SHEET not in wb.sheetnames:
        raise KeyError(f"Workbook must contain sheet {SHEET!r}; found {wb.sheetnames!r}")

    ws = wb[SHEET]
    headers = [ws.cell(1, col).value for col in range(1, 6)]
    if headers != EXPECTED_HEADERS:
        raise ValueError(
            "Unexpected assumptions_inputs layout. Expected first five headers "
            f"{EXPECTED_HEADERS!r}; found {headers!r}."
        )

    result: dict[str, Any] = {}
    seen: set[str] = set()

    for excel_row in range(2, ws.max_row + 1):
        raw_key = ws.cell(excel_row, 1).value
        value = ws.cell(excel_row, 2).value

        if raw_key in (None, ""):
            continue
        key = str(raw_key).strip()
        if not key:
            continue
        if key in seen:
            raise ValueError(f"Duplicate YAML key {key!r} at Excel row {excel_row}")
        seen.add(key)

        # These keys are deliberately unresolved. They may remain in the clean
        # input sheet for documentation/provenance, but their placeholder Value
        # is ignored until they are removed from NULL_KEYS.
        if key in NULL_KEYS:
            nested_set(result, key, None)
            continue

        if value is None:
            raise ValueError(
                f"Excel row {excel_row} ({key}) has a blank Value. "
                "Populate it, or add the key to NULL_KEYS if it is deliberately unresolved."
            )

        nested_set(result, key, value_for_yaml(key, value))

    # Also add unresolved interface keys that are not represented as rows in Excel.
    for key in NULL_KEYS:
        if key not in seen:
            nested_set(result, key, None)

    return reorder(result)


def write_yaml(data: dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# AUTO-GENERATED by scratch/generate_assumptions.py.\n"
        "# Source values: data/assumptions/Assumptions_yaml_inputs.xlsx, sheet assumptions_inputs.\n"
        "# Excel Value is written directly; no unit conversion is performed.\n"
        "# null = deliberately retained model input whose value is still being researched.\n\n"
    )
    body = yaml.safe_dump(
        data,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
        width=1000,
    )
    output_path.write_text(header + body, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("output", nargs="?", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = build_assumptions(args.input)
    write_yaml(data, args.output)
    print(f"Generated {args.output} from {args.input}")


if __name__ == "__main__":
    main()
