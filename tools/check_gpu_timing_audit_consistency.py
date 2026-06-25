#!/usr/bin/env python3
"""Verify published timing notes match the validation evidence ledger."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
AUDIT_PATH = Path("docs/gpu_port_coverage_audit.md")

CASE_LABELS = {
    "les_core_channel": "channel flow",
    "adm_disk": "ADM disk",
    "atm_line": "ATM line",
}


def load_cases() -> dict[str, dict]:
    data = json.loads(EVIDENCE_PATH.read_text(encoding="utf-8"))
    return {
        case["id"]: case
        for case in data.get("cases", [])
        if isinstance(case, dict) and isinstance(case.get("id"), str)
    }


def parse_timing_table() -> tuple[dict[str, dict[str, float]], dict[str, float] | None]:
    rows: dict[str, dict[str, float]] = {}
    total: dict[str, float] | None = None
    headers: list[str] | None = None

    def split_row(line: str) -> list[str]:
        return [cell.strip() for cell in line.strip().strip("|").split("|")]

    def numeric(cell: str, suffix: str | None = None) -> float:
        value = cell.strip()
        if suffix and value.endswith(suffix):
            value = value[: -len(suffix)].strip()
        return float(value)

    for line in AUDIT_PATH.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            continue
        cells = split_row(line)
        if not cells or cells[0] == "---" or cells[0].startswith("---"):
            continue
        if cells[0] == "Case":
            headers = cells
            continue
        if headers is None or len(cells) != len(headers):
            continue
        row = dict(zip(headers, cells, strict=True))
        label = row["Case"]
        if label == "three-case total":
            total = {
                "cumulative_seconds": numeric(row["Cumulative wall"], "s"),
                "nsteps": int(row["Steps"]),
                "cumulative_average_s_per_step": numeric(row["Cumulative average"], "s/step"),
            }
            continue
        if label in CASE_LABELS.values():
            rows[label] = {
                "final_printed_s_per_step": numeric(row["Final printed timing"], "s/step"),
                "cumulative_seconds": numeric(row["Cumulative wall"], "s"),
                "nsteps": int(row["Steps"]),
                "cumulative_average_s_per_step": numeric(row["Cumulative average"], "s/step"),
            }
    return rows, total


def close_enough(left: float, right: float, tolerance: float = 5e-6) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def main() -> int:
    cases = load_cases()
    audit_rows, audit_total = parse_timing_table()
    errors: list[str] = []

    expected_totals = {
        "cumulative_seconds": 0.0,
        "nsteps": 0,
    }

    for case_id, label in CASE_LABELS.items():
        case = cases.get(case_id)
        if not case:
            errors.append(f"{EVIDENCE_PATH} is missing case `{case_id}`")
            continue
        runtime = case.get("gpu_runtime")
        if not isinstance(runtime, dict):
            errors.append(f"{EVIDENCE_PATH} case `{case_id}` has no gpu_runtime object")
            continue
        row = audit_rows.get(label)
        if row is None:
            errors.append(f"{AUDIT_PATH} timing table is missing row `{label}`")
            continue

        for field in [
            "final_printed_s_per_step",
            "cumulative_seconds",
            "cumulative_average_s_per_step",
        ]:
            expected = runtime.get(field)
            observed = row[field]
            if not isinstance(expected, (int, float)):
                errors.append(f"{EVIDENCE_PATH} case `{case_id}` field `{field}` is not numeric")
                continue
            if not close_enough(float(expected), observed):
                errors.append(
                    f"{AUDIT_PATH} row `{label}` field `{field}` is {observed}, "
                    f"but evidence has {expected}"
                )

        expected_steps = runtime.get("nsteps")
        if expected_steps != row["nsteps"]:
            errors.append(
                f"{AUDIT_PATH} row `{label}` nsteps is {row['nsteps']}, "
                f"but evidence has {expected_steps}"
            )

        expected_totals["cumulative_seconds"] += float(runtime["cumulative_seconds"])
        expected_totals["nsteps"] += int(runtime["nsteps"])

    if audit_total is None:
        errors.append(f"{AUDIT_PATH} timing table is missing `three-case total` row")
    else:
        expected_average = (
            expected_totals["cumulative_seconds"] / expected_totals["nsteps"]
            if expected_totals["nsteps"]
            else 0.0
        )
        for field, expected in [
            ("cumulative_seconds", expected_totals["cumulative_seconds"]),
            ("nsteps", expected_totals["nsteps"]),
            ("cumulative_average_s_per_step", expected_average),
        ]:
            observed = audit_total[field]
            if field == "nsteps":
                if observed != expected:
                    errors.append(
                        f"{AUDIT_PATH} three-case total `{field}` is {observed}, "
                        f"expected {expected}"
                    )
            elif not close_enough(float(observed), float(expected)):
                errors.append(
                    f"{AUDIT_PATH} three-case total `{field}` is {observed}, "
                    f"expected {expected:.8f}"
                )

    if errors:
        print("GPU timing audit consistency check failed:")
        for error in errors:
            print(f"  {error}")
        return 1

    print(
        "GPU timing audit consistency check passed "
        f"({len(CASE_LABELS)} rows plus total)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
