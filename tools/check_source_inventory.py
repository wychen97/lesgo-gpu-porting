#!/usr/bin/env python3
"""Verify every tracked Fortran source has one source-inventory table row."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import public_knob_names
from fortran_inventory import normalized_repo_path, tracked_fortran_paths


ROOT = Path(__file__).resolve().parents[1]
INVENTORY_PATH = ROOT / "docs/source_file_inventory.md"
INVENTORY_ROW_RE = re.compile(
    r"^\|\s*`(?P<path>[^`]+\.(?:f|f90|F|F90))`\s*"
    r"\|\s*(?P<build>[^|]+)\|"
)
BUILD_OPTION_RE = re.compile(r"\bUSE_[A-Za-z0-9_]+\b")
LITERAL_BUILD_PATHS = {"common", "helper only"}


def tracked_fortran_sources() -> set[str]:
    return {normalized_repo_path(path) for path in tracked_fortran_paths()}


def cmake_build_options() -> set[str]:
    return public_knob_names("USE_")


def inventory_fortran_rows() -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    for line in INVENTORY_PATH.read_text(encoding="utf-8").splitlines():
        match = INVENTORY_ROW_RE.match(line)
        if not match:
            continue
        path = match.group("path").replace("\\", "/")
        build_path = " ".join(match.group("build").strip().split())
        rows.setdefault(path, []).append(build_path)
    return rows


def invalid_build_paths(
    documented_rows: dict[str, list[str]], options: set[str]
) -> list[str]:
    invalid: list[str] = []
    for path, build_paths in documented_rows.items():
        for build_path in build_paths:
            if build_path in LITERAL_BUILD_PATHS:
                continue
            referenced_options = set(BUILD_OPTION_RE.findall(build_path))
            if not referenced_options:
                invalid.append(f"{path}: unrecognized build path: {build_path}")
                continue
            unknown = sorted(referenced_options - options)
            if unknown:
                invalid.append(
                    f"{path}: build path references unknown CMake option(s): "
                    f"{', '.join(unknown)}"
                )
    return invalid


def main() -> int:
    tracked = tracked_fortran_sources()
    documented_rows = inventory_fortran_rows()
    documented = set(documented_rows)
    options = cmake_build_options()

    missing = sorted(tracked - documented)
    stale = sorted(documented - tracked)
    duplicates = sorted(path for path, rows in documented_rows.items() if len(rows) > 1)
    invalid_paths = invalid_build_paths(documented_rows, options)

    if missing or stale or duplicates or invalid_paths:
        print(f"{INVENTORY_PATH.relative_to(ROOT)} is out of sync:")
        for path in missing:
            print(f"  missing tracked source table row: {path}")
        for path in stale:
            print(f"  stale inventory entry: {path}")
        for path in duplicates:
            print(f"  duplicate inventory entry: {path} ({len(documented_rows[path])} rows)")
        for item in invalid_paths:
            print(f"  {item}")
        return 1

    print(
        "Source inventory check passed "
        f"({len(tracked)} Fortran table rows, {len(options)} CMake options)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
