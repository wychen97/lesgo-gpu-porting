#!/usr/bin/env python3
"""Generate key-level lesgo.conf validation coverage from source and manifests."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from repo_paths import repo_path
from check_lesgo_conf_coverage_docs import source_keys


MAP_PATH = repo_path("docs", "lesgo_conf_gpu_validation_map.json")
EVIDENCE_PATH = repo_path("docs", "gpu_validation_evidence.json")
OUTPUT_PATH = repo_path("docs", "lesgo_conf_key_validation_coverage.md")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def group_rows() -> dict[str, list[str]]:
    validation_map = load_json(MAP_PATH)
    rows: dict[str, list[str]] = {}
    for entry in validation_map.get("parser_groups", []):
        if not isinstance(entry, dict):
            continue
        group = entry.get("group")
        validation_rows = entry.get("validation_rows")
        if isinstance(group, str) and isinstance(validation_rows, list):
            rows[group] = [
                row_id
                for row_id in validation_rows
                if isinstance(row_id, str) and row_id
            ]
    return rows


def evidence_statuses() -> dict[str, str]:
    evidence = load_json(EVIDENCE_PATH)
    return {
        case["id"]: case["evidence_status"]
        for case in evidence.get("cases", [])
        if (
            isinstance(case, dict)
            and isinstance(case.get("id"), str)
            and isinstance(case.get("evidence_status"), str)
        )
    }


def markdown_table() -> str:
    keys_by_group = source_keys()
    rows_by_group = group_rows()
    statuses = evidence_statuses()
    lines = [
        "# lesgo.conf Key Validation Coverage",
        "",
        "This generated table expands the non-LVLSET `input_util.f90` parser",
        "keys into their responsible GPU validation rows.  A mapping means the",
        "row is responsible for correctness and timing evidence; it is not, by",
        "itself, a speedup claim.",
        "",
        "| Parser group | Key | Validation rows | Current evidence states |",
        "| --- | --- | --- | --- |",
    ]
    for group in [
        "DOMAIN",
        "MODEL",
        "CORIOLIS",
        "FLOW_COND",
        "OUTPUT",
        "TURBINES",
        "SCALARS",
        "HIT",
    ]:
        validation_rows = rows_by_group.get(group, [])
        validation_text = ", ".join(f"`{row_id}`" for row_id in validation_rows)
        status_text = ", ".join(
            f"`{statuses.get(row_id, 'missing-evidence-row')}`"
            for row_id in validation_rows
        )
        for key in keys_by_group.get(group, []):
            lines.append(
                f"| `{group}` | `{key}` | {validation_text} | {status_text} |"
            )
    lines.extend(
        [
            "",
            "Regenerate this file with:",
            "",
            "```bash",
            "python3 tools/report_lesgo_conf_key_validation.py --write",
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report key-level lesgo.conf GPU validation coverage."
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help=f"write the generated report to {OUTPUT_PATH}",
    )
    args = parser.parse_args()

    text = markdown_table()
    if args.write:
        OUTPUT_PATH.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
