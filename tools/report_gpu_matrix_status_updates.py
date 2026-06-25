#!/usr/bin/env python3
"""Suggest validation-matrix status updates from the evidence ledger."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from gpu_validation_plan import CASE_VARIANTS, CORRECTNESS_ONLY_VARIANTS


MATRIX_PATH = Path("docs/gpu_validation_matrix.md")
EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
ROW_RE = re.compile(r"^\|\s*`(?P<id>[a-z0-9_]+)`\s*\|(?P<body>.*)\|$")
STATUS_RE = re.compile(r"`(?P<status>[a-z-]+)`")

EVIDENCE_TO_MATRIX = {
    "gpu_runtime_only": "validated-gpu-runtime",
    "needs_benchmark": "implemented-needs-benchmark",
    "external_record_needs_copy": "recorded-correct-and-faster",
    "paired_speedup_claimed": "recorded-correct-and-faster",
    "paired_not_faster": "recorded-paired-not-faster",
    "host_boundary_pending": "host-boundary",
    "host_boundary_verified": "host-boundary",
    "excluded": "excluded",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def matrix_statuses(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    in_matrix = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip() == "## Validation Matrix":
            in_matrix = True
            continue
        if in_matrix and line.startswith("## "):
            break
        if not in_matrix:
            continue
        match = ROW_RE.match(line)
        if not match:
            continue
        statuses = STATUS_RE.findall(match.group("body"))
        if len(statuses) == 1:
            rows[match.group("id")] = statuses[0]
    return rows


def evidence_cases(path: Path) -> dict[str, dict]:
    evidence = load_json(path)
    return {
        case["id"]: case
        for case in evidence.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }


def variant_recommendation(case_id: str, case: dict) -> str | None:
    variants = CASE_VARIANTS.get(case_id, [])
    if len(variants) <= 1:
        return None
    expected = {variant.id for variant in variants}
    results = case.get("variant_pair_results")
    if not isinstance(results, dict):
        return None
    if not expected.issubset(results):
        return None

    speedup_variants = {
        variant.id
        for variant in variants
        if (case_id, variant.id) not in CORRECTNESS_ONLY_VARIANTS
    }
    outcomes = {
        result.get("outcome")
        for variant, result in results.items()
        if variant in speedup_variants and isinstance(result, dict)
    }
    if outcomes == {"gpu_faster"}:
        return "recorded-correct-and-faster"
    if "gpu_not_faster" in outcomes:
        return "recorded-paired-not-faster"
    return None


def recommended_status(case_id: str, case: dict) -> str:
    variant_status = variant_recommendation(case_id, case)
    if variant_status is not None:
        return variant_status
    evidence_status = case.get("evidence_status")
    return EVIDENCE_TO_MATRIX.get(str(evidence_status), "unknown")


def render_markdown(matrix: dict[str, str], evidence: dict[str, dict]) -> str:
    rows = []
    for case_id in sorted(evidence):
        current = matrix.get(case_id, "<missing>")
        recommended = recommended_status(case_id, evidence[case_id])
        if current != recommended:
            rows.append((case_id, current, recommended))

    lines = [
        "# GPU Validation Matrix Status Update Report",
        "",
        "This report compares `docs/gpu_validation_matrix.md` against "
        "`docs/gpu_validation_evidence.json`.",
        "",
    ]
    if not rows:
        lines.append("No matrix status updates are currently suggested.")
        return "\n".join(lines) + "\n"

    lines.extend(
        [
            "| ID | Current matrix status | Suggested status |",
            "| --- | --- | --- |",
        ]
    )
    for case_id, current, recommended in rows:
        lines.append(f"| `{case_id}` | `{current}` | `{recommended}` |")
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Suggest validation-matrix status updates from evidence."
    )
    parser.add_argument("--matrix", type=Path, default=MATRIX_PATH)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE_PATH)
    parser.add_argument(
        "--fail-on-suggestion",
        action="store_true",
        help="return nonzero when any matrix update is suggested",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    matrix = matrix_statuses(args.matrix)
    evidence = evidence_cases(args.evidence)
    report = render_markdown(matrix, evidence)
    print(report, end="")
    has_suggestion = "| `" in report
    return 1 if args.fail_on_suggestion and has_suggestion else 0


if __name__ == "__main__":
    sys.exit(main())
