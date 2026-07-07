#!/usr/bin/env python3
"""Generate a full static GPU classification inventory for all subprograms."""

from __future__ import annotations

import argparse
import sys
from collections import Counter

from report_gpu_static_inventory import (
    BUCKET_VALIDATION_ROWS,
    ROOT,
    collect,
    review_bucket,
)


OUTPUT_PATH = ROOT / "docs/gpu_static_full_inventory.md"


def markdown() -> str:
    subprograms = collect()
    counts = Counter(item.classification for item in subprograms)
    lines = [
        "# Static GPU Full Inventory",
        "",
        "This generated report lists every tracked root-level non-LVLSET",
        "Fortran subroutine/function and its static GPU-audit classification.",
        "It is a source heuristic, not proof of runtime speedup.  Use the",
        "validation matrix and evidence ledger for correctness and timing",
        "claims.",
        "",
        "| Classification | Subprograms |",
        "| --- | ---: |",
    ]
    for classification in [
        "gpu-marked",
        "gpu-file-unmarked",
        "host-boundary",
        "host-or-diagnostic",
        "unmarked-runtime-candidate",
    ]:
        lines.append(f"| `{classification}` | {counts[classification]} |")

    lines.extend(
        [
            "",
            "| File:line | Subprogram | Classification | Review bucket | Validation rows |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for item in sorted(subprograms, key=lambda entry: (str(entry.path), entry.start_line, entry.name)):
        bucket = ""
        validation_rows = ""
        if item.classification == "unmarked-runtime-candidate":
            bucket = review_bucket(item)
            validation_rows = ", ".join(
                f"`{row_id}`" for row_id in BUCKET_VALIDATION_ROWS.get(bucket, [])
            )
        bucket_text = f"`{bucket}`" if bucket else ""
        lines.append(
            f"| `{item.path}:{item.start_line}` | `{item.name}` | "
            f"`{item.classification}` | {bucket_text} | {validation_rows} |"
        )

    lines.extend(
        [
            "",
            "Regenerate this file with:",
            "",
            "```bash",
            "python3 tools/report_gpu_static_full_inventory.py --write",
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate full static GPU classification documentation."
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
