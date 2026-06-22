#!/usr/bin/env python3
"""Verify root CMake build options have collaborator-facing table rows."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import public_knob_names


DOC_PATH = Path("docs/code_organization.md")
OPTION_ROW_RE = re.compile(r"^\|\s*`(?P<name>USE_[A-Za-z0-9_]+)`\s*\|")


def collect_cmake_options() -> set[str]:
    return public_knob_names("USE_")


def collect_documented_option_rows() -> dict[str, int]:
    rows: dict[str, int] = {}
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        match = OPTION_ROW_RE.match(line)
        if not match:
            continue
        name = match.group("name")
        rows[name] = rows.get(name, 0) + 1
    return rows


def main() -> int:
    cmake_options = collect_cmake_options()
    documented_rows = collect_documented_option_rows()
    documented = set(documented_rows)

    missing = sorted(cmake_options - documented)
    stale = sorted(documented - cmake_options)
    duplicates = sorted(name for name, count in documented_rows.items() if count > 1)

    if missing or stale or duplicates:
        print(f"{DOC_PATH} build-option table is out of sync:")
        for name in missing:
            print(f"  missing root CMake build option row: {name}")
        for name in stale:
            print(f"  stale build-option row: {name}")
        for name in duplicates:
            print(f"  duplicate build-option row: {name} ({documented_rows[name]} rows)")
        return 1

    print(
        "CMake option documentation check passed "
        f"({len(cmake_options)} root build option rows)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
