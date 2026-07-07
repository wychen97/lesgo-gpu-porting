#!/usr/bin/env python3
"""Verify large active Fortran files have a source navigation map."""

from __future__ import annotations

import sys
from pathlib import Path

from fortran_inventory import tracked_fortran_files
from repo_paths import ROOT


LINE_THRESHOLD = 900
HEADER_SCAN_LINES = 180
MAP_MARKERS = ("Navigation map", "Routine map")

# LVLSET remains outside the optimized production path for this branch.  The
# files are tracked, but this readability gate is scoped to active production
# and optional non-LVLSET modules.
EXCLUDED_PREFIXES = ("trees_",)
EXCLUDED_NAMES = {
    "level_set.f90",
    "level_set_base.f90",
}


def is_excluded(path: Path) -> bool:
    name = path.name
    if name in EXCLUDED_NAMES:
        return True
    return any(name.startswith(prefix) for prefix in EXCLUDED_PREFIXES)


def line_count(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def has_navigation_map(path: Path) -> bool:
    lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for _ in range(HEADER_SCAN_LINES):
            line = handle.readline()
            if not line:
                break
            lines.append(line)
    header = "".join(lines)
    return any(marker in header for marker in MAP_MARKERS)


def main() -> int:
    missing: list[str] = []
    checked = 0

    for path in tracked_fortran_files():
        if is_excluded(path):
            continue
        lines = line_count(path)
        if lines < LINE_THRESHOLD:
            continue
        checked += 1
        if not has_navigation_map(path):
            rel = path.relative_to(ROOT)
            missing.append(f"{rel} ({lines} lines)")

    if missing:
        print("Large Fortran navigation-map check failed:")
        for item in missing:
            print(f"  missing Navigation map/Routine map: {item}")
        return 1

    print(f"Navigation-map check passed ({checked} large active Fortran files).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
