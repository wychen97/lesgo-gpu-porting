#!/usr/bin/env python3
"""Verify CONTRIBUTING.md points to the readiness source of truth."""

from __future__ import annotations

import sys
from pathlib import Path


CONTRIBUTING_PATH = Path("CONTRIBUTING.md")

REQUIRED_SNIPPETS = (
    "python3 tools/check_branch_readiness.py",
    "python3 tools/check_branch_readiness.py --with-cmake-configure",
    "python3 tools/check_branch_readiness.py --with-hit-cmake-configure",
    'README_FINAL_OPTIMIZED_20260619.md` under "Local Readability Checks"',
    "`tools/README.md` is the index for the individual readiness tools",
)


def main() -> int:
    text = CONTRIBUTING_PATH.read_text(encoding="utf-8")
    missing = [snippet for snippet in REQUIRED_SNIPPETS if snippet not in text]

    if missing:
        print(f"{CONTRIBUTING_PATH} is missing readiness guidance:")
        for snippet in missing:
            print(f"  {snippet}")
        return 1

    print(
        f"Contributing readiness-doc check passed "
        f"({len(REQUIRED_SNIPPETS)} snippets)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
