#!/usr/bin/env python3
"""Verify the generated static GPU candidate review report is current."""

from __future__ import annotations

import sys

from report_gpu_static_candidate_review import OUTPUT_PATH, markdown


def main() -> int:
    expected = markdown()
    if not OUTPUT_PATH.exists():
        print(f"{OUTPUT_PATH} is missing")
        print("Regenerate it with:")
        print("  python3 tools/report_gpu_static_candidate_review.py --write")
        return 1

    actual = OUTPUT_PATH.read_text(encoding="utf-8")
    if actual != expected:
        print(f"{OUTPUT_PATH} is out of date")
        print("Regenerate it with:")
        print("  python3 tools/report_gpu_static_candidate_review.py --write")
        return 1

    row_count = sum(1 for line in actual.splitlines() if line.startswith("| `") and ".f90:" in line)
    print(
        "GPU static candidate review report check passed "
        f"({row_count} candidates)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
