#!/usr/bin/env python3
"""Verify public CMake cache variables have collaborator-facing table rows."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import cache_variable_names


DOC_PATH = Path("docs/build_profiles.md")
CACHE_ROW_RE = re.compile(r"^\|\s*`(?P<name>[A-Za-z0-9_]+)`\s*\|")


def collect_cache_variables() -> set[str]:
    return cache_variable_names()


def collect_documented_variable_rows() -> dict[str, int]:
    rows: dict[str, int] = {}
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        match = CACHE_ROW_RE.match(line)
        if not match:
            continue
        name = match.group("name")
        rows[name] = rows.get(name, 0) + 1
    return rows


def main() -> int:
    cache_variables = collect_cache_variables()
    documented_rows = collect_documented_variable_rows()
    documented = set(documented_rows)

    missing = sorted(cache_variables - documented)
    stale = sorted(documented - cache_variables)
    duplicates = sorted(name for name, count in documented_rows.items() if count > 1)

    if missing or stale or duplicates:
        print(f"{DOC_PATH} cache-variable table is out of sync:")
        for name in missing:
            print(f"  missing public CMake cache variable row: {name}")
        for name in stale:
            print(f"  stale cache-variable row: {name}")
        for name in duplicates:
            print(f"  duplicate cache-variable row: {name} ({documented_rows[name]} rows)")
        return 1

    print(
        "CMake cache-variable documentation check passed "
        f"({len(cache_variables)} public cache-variable rows)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
