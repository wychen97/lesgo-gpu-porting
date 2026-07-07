#!/usr/bin/env python3
"""Verify root Fortran sources are covered by named CMake source groups."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from cmake_metadata import CMAKE_PATH, cmake_source_groups, source_list_group_references
from repo_paths import ROOT

FORTRAN_SUFFIXES = {".f", ".f90", ".F", ".F90"}


def tracked_root_fortran_sources() -> set[str]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    paths = result.stdout.decode("utf-8").split("\0")
    return {
        path
        for path in paths
        if path
        and Path(path).parent == Path(".")
        and Path(path).suffix in FORTRAN_SUFFIXES
    }


def main() -> int:
    tracked = tracked_root_fortran_sources()
    source_groups = cmake_source_groups()

    grouped: dict[str, list[str]] = {}
    for group_name, sources in source_groups.items():
        for source in sources:
            grouped.setdefault(source, []).append(group_name)

    grouped_sources = set(grouped)
    missing = sorted(tracked - grouped_sources)
    stale = sorted(grouped_sources - tracked)
    duplicates = sorted(
        source for source, groups in grouped.items() if len(groups) > 1
    )
    empty_groups = sorted(name for name, sources in source_groups.items() if not sources)

    source_references = source_list_group_references()
    unused_groups = sorted(name for name in source_groups if name not in source_references)
    undefined_references = sorted(source_references - set(source_groups))

    if (
        missing
        or stale
        or duplicates
        or empty_groups
        or unused_groups
        or undefined_references
    ):
        print(f"{CMAKE_PATH.relative_to(ROOT)} source-group check failed:")
        for source in missing:
            print(f"  missing from named source groups: {source}")
        for source in stale:
            print(f"  source group references missing/untracked file: {source}")
        for source in duplicates:
            print(
                f"  source appears in multiple groups: {source} "
                f"({', '.join(grouped[source])})"
            )
        for group in empty_groups:
            print(f"  empty source group: {group}")
        for group in unused_groups:
            print(f"  source group is never appended to Sources: {group}")
        for group in undefined_references:
            print(f"  Sources references undefined source group: {group}")
        return 1

    print(
        "CMake source-group check passed "
        f"({len(tracked)} root Fortran files, {len(source_groups)} groups)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
