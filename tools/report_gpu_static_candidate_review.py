#!/usr/bin/env python3
"""Generate the full static unmarked-candidate GPU review report."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict

from repo_paths import repo_path
from report_gpu_static_inventory import (
    BUCKET_VALIDATION_ROWS,
    REVIEW_BUCKET_DESCRIPTIONS,
    Subprogram,
    collect,
    review_bucket,
)


OUTPUT_PATH = repo_path("docs", "gpu_static_candidate_review.md")


def candidate_items() -> list[Subprogram]:
    return [
        item
        for item in collect()
        if item.classification == "unmarked-runtime-candidate"
    ]


def markdown() -> str:
    candidates = candidate_items()
    by_bucket: dict[str, list[Subprogram]] = defaultdict(list)
    for item in candidates:
        by_bucket[review_bucket(item)].append(item)

    lines = [
        "# Static GPU Candidate Review",
        "",
        "This generated report lists every tracked non-LVLSET Fortran",
        "subroutine/function that the static scanner classifies as an",
        "`unmarked-runtime-candidate`.  These entries are not automatically",
        "missing GPU kernels: many are host-model logic, fallback compatibility,",
        "diagnostics, or setup/control code.  The validation rows show where",
        "runtime correctness and timing evidence must close each bucket.",
        "",
        "| Review bucket | Candidates | Validation rows | Meaning |",
        "| --- | ---: | --- | --- |",
    ]
    for bucket in sorted(by_bucket):
        validation_rows = ", ".join(
            f"`{row_id}`" for row_id in BUCKET_VALIDATION_ROWS.get(bucket, [])
        )
        lines.append(
            f"| `{bucket}` | {len(by_bucket[bucket])} | {validation_rows} | "
            f"{REVIEW_BUCKET_DESCRIPTIONS[bucket]} |"
        )

    lines.append("")
    lines.append("| File:line | Subprogram | Review bucket | Validation rows |")
    lines.append("| --- | --- | --- | --- |")
    for item in sorted(candidates, key=lambda entry: (str(entry.path), entry.start_line, entry.name)):
        bucket = review_bucket(item)
        validation_rows = ", ".join(
            f"`{row_id}`" for row_id in BUCKET_VALIDATION_ROWS.get(bucket, [])
        )
        lines.append(
            f"| `{item.path}:{item.start_line}` | `{item.name}` | "
            f"`{bucket}` | {validation_rows} |"
        )

    lines.extend(
        [
            "",
            "Regenerate this file with:",
            "",
            "```bash",
            "python3 tools/report_gpu_static_candidate_review.py --write",
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate static GPU candidate review documentation."
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help=f"write the generated report to {OUTPUT_PATH}",
    )
    args = parser.parse_args()

    text = markdown()
    if args.write:
        OUTPUT_PATH.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
