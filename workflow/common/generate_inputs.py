#!/usr/bin/env python3
"""Write the model's input files from the team's workbook: config/assumptions.yaml from
config/.assumptions.yaml and data/assumptions/Assumptions_yaml_inputs.xlsx.

Run: python -m common.generate_inputs [WORKBOOK] [--verbose]
"""
import argparse
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from common._paths import REPO_ROOT

log = logging.getLogger(__name__)

TEMPLATE_PATH = REPO_ROOT / "config" / ".assumptions.yaml"
OUTPUT_PATH = REPO_ROOT / "config" / "assumptions.yaml"
DEFAULT_WORKBOOK = REPO_ROOT / "data" / "assumptions" / "Assumptions_yaml_inputs.xlsx"
INPUTS_SHEET = "assumptions_inputs"
INPUTS_COLUMNS = {"Confirmed row ID": "id", "YAML key": "workbook_key", "Value": "value", "Unit": "unit"}
# ECB reference rate, 2025 annual average (EXR.A.USD.EUR.SP00.A): the rate the
# master register converts with (its Conversions tab, EUR_to_USD_2025).
USD_PER_EUR = 1.13
COMBINATIONS = {"mean": np.mean, "sum": np.sum}
GENERATED_HEADER = f"""\
# GENERATED FILE: edit config/.assumptions.yaml instead, then run
#   python -m common.generate_inputs
# Every value tagged `!register` there is filled in here from the input
# workbook, converted from USD to EUR at {USD_PER_EUR} USD/EUR, with the register
# row ID, the workbook's key, value and unit in the comment above it.

"""


@dataclass
class RegisterReference:
    spec: int | dict
    line: int
    start_column: int
    end_column: int
    path: str = ""
    values: dict[str, float] = field(default_factory=dict)
    comments: list[str] = field(default_factory=list)


class TemplateLoader(yaml.SafeLoader):
    pass


def construct_register_reference(loader: TemplateLoader, node: yaml.Node) -> RegisterReference:
    if node.start_mark.line != node.end_mark.line:
        raise ValueError(f"line {node.start_mark.line + 1}: a !register reference must fit on one line")
    spec = (int(loader.construct_scalar(node)) if isinstance(node, yaml.ScalarNode)
            else loader.construct_mapping(node, deep=True))
    reference = RegisterReference(spec, node.start_mark.line, node.start_mark.column, node.end_mark.column)
    return reference


TemplateLoader.add_constructor("!register", construct_register_reference)


def read_workbook_inputs(workbook_path: Path) -> pd.DataFrame:
    """The workbook's input rows indexed by register row ID; rows without one are left out."""
    workbook_inputs = (
        pd.read_excel(workbook_path, sheet_name=INPUTS_SHEET)
        .loc[:, list(INPUTS_COLUMNS)]
        .rename(columns=INPUTS_COLUMNS)
        .dropna(subset=["id"])
        .astype({"id": int})
        .set_index("id")
    )
    if not workbook_inputs.index.is_unique:
        raise ValueError(f"Confirmed row IDs used twice: {workbook_inputs.index[workbook_inputs.index.duplicated()].tolist()}")
    return workbook_inputs


def read_template_references(template_text: str) -> tuple[list[RegisterReference], list[str]]:
    """Every `!register` value in the template in file order, and the paths of the backstop values."""
    pending = [("", yaml.load(template_text, Loader=TemplateLoader))]
    references, backstop_paths = [], []
    while pending:
        path, node = pending.pop()
        if isinstance(node, RegisterReference):
            node.path = path
            references.append(node)
        elif isinstance(node, dict):
            pending.extend((f"{path}.{key}".lstrip("."), value) for key, value in node.items())
        else:
            backstop_paths.append(path)
    ordered_references = sorted(references, key=lambda r: (r.line, r.start_column))
    return ordered_references, sorted(backstop_paths)


def resolve_reference(reference: RegisterReference, workbook_inputs: pd.DataFrame) -> None:
    """Fill in the reference's value per alternative (the first is the default) and its comments."""
    is_alternatives = isinstance(reference.spec, dict) and not set(reference.spec) <= set(COMBINATIONS)
    alternatives = reference.spec if is_alternatives else {"default": reference.spec}
    row_ids = {name: [entry] if isinstance(entry, int) else next(iter(entry.values()))
               for name, entry in alternatives.items()}
    unknown_ids = sorted({i for ids in row_ids.values() for i in ids} - set(workbook_inputs.index))
    if unknown_ids:
        raise KeyError(f"{reference.path}: no row with Confirmed row ID {unknown_ids} in the workbook")
    rows = workbook_inputs.loc[[i for ids in row_ids.values() for i in ids]]
    if rows["value"].isna().any():
        raise ValueError(f"{reference.path}: workbook rows {rows.index[rows['value'].isna()].tolist()} have no value")
    if rows["unit"].nunique() > 1:
        raise ValueError(f"{reference.path}: workbook rows are in different units {rows['unit'].unique().tolist()}")
    unit = rows["unit"].iloc[0]
    if ("eur" in reference.path.rsplit(".", 1)[-1].split("_")) != ("USD" in unit or "EUR" in unit):
        raise ValueError(f"{reference.path}: workbook unit '{unit}' does not match the key's currency")

    row_labels = {}
    for name, entry in alternatives.items():
        combine = COMBINATIONS[next(iter(entry))] if isinstance(entry, dict) else np.sum
        combined = float(combine(workbook_inputs.loc[row_ids[name], "value"]))
        reference.values[name] = combined / USD_PER_EUR if "USD" in unit else combined
        operation = f"{next(iter(entry))} of " if isinstance(entry, dict) else ""
        row_labels[name] = operation + ", ".join(map(str, row_ids[name]))
    default_name = next(iter(alternatives))
    default_label = f"{default_name} = {row_labels[default_name]}" if is_alternatives else row_labels[default_name]
    reference.comments = [f"register {default_label}: " + ", ".join(
        f"{row.workbook_key} {format_number(row.value)}" for row in rows.drop_duplicates().itertuples()
    ) + f" {unit}"]
    if is_alternatives:
        reference.comments.append("alternatives in EUR: " + ", ".join(
            f"{name} = {format_number(reference.values[name])}" for name in list(alternatives)[1:]))


def format_number(value: float) -> str:
    """Six significant digits, never in exponent notation (YAML 1.1 would read that as a string)."""
    formatted_number = np.format_float_positional(float(value), precision=6, unique=False, fractional=False, trim="-")
    return formatted_number


def render_assumptions(template_text: str, references: list[RegisterReference]) -> str:
    lines = template_text.splitlines()
    for reference in reversed(references):
        line = lines[reference.line]
        default_value = format_number(next(iter(reference.values.values())))
        lines[reference.line] = line[:reference.start_column] + default_value + line[reference.end_column:]
    for line_number in sorted({r.line for r in references}, reverse=True):
        indentation = re.match(r"\s*", lines[line_number]).group()
        comments = dict.fromkeys(c for r in references if r.line == line_number for c in r.comments)
        lines[line_number:line_number] = [f"{indentation}# {comment}" for comment in comments]
    # The template's opening comment block explains the template, not the output.
    body_start = next(n for n, line in enumerate(lines) if not line.startswith("#"))
    rendered_text = GENERATED_HEADER + "\n".join(lines[body_start:]).lstrip("\n") + "\n"
    return rendered_text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", nargs="?", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--verbose", "-v", action="store_true", help="also list the backstop values")
    arguments = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log.setLevel(logging.DEBUG if arguments.verbose else logging.INFO)

    workbook_inputs = read_workbook_inputs(arguments.workbook)
    template_text = TEMPLATE_PATH.read_text()
    references, backstop_paths = read_template_references(template_text)
    for reference in references:
        resolve_reference(reference, workbook_inputs)
    rendered_text = render_assumptions(template_text, references)

    unchanged = OUTPUT_PATH.exists() and OUTPUT_PATH.read_text() == rendered_text
    if not unchanged:
        OUTPUT_PATH.write_text(rendered_text)
    log.info(f"{len(references)} values from {arguments.workbook.name}, "
             f"{len(backstop_paths)} backstop values from {TEMPLATE_PATH.name}: "
             f"{OUTPUT_PATH.relative_to(REPO_ROOT)} {'unchanged' if unchanged else 'written'}")
    log.debug("backstop values:\n  " + "\n  ".join(backstop_paths))


if __name__ == "__main__":
    main()
