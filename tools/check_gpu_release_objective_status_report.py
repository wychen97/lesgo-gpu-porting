#!/usr/bin/env python3
"""Verify the generated GPU release-objective status report is current."""

from __future__ import annotations

import sys

from report_gpu_release_objective_status import OUTPUT_PATH, markdown


def main() -> int:
    expected = markdown()
    if not OUTPUT_PATH.exists():
        print(f"{OUTPUT_PATH} is missing")
        print("Regenerate it with:")
        print("  python3 tools/report_gpu_release_objective_status.py --write")
        return 1

    actual = OUTPUT_PATH.read_text(encoding="utf-8")
    if actual != expected:
        print(f"{OUTPUT_PATH} is out of date")
        print("Regenerate it with:")
        print("  python3 tools/report_gpu_release_objective_status.py --write")
        return 1

    gap_line = next(
        line for line in actual.splitlines() if line.startswith("Strict release-gate gaps:")
    )
    print(f"GPU release-objective status report check passed ({gap_line}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
