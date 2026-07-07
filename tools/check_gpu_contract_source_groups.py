#!/usr/bin/env python3
"""Verify GPU contracts document every named CMake source group."""

from __future__ import annotations

import re
import sys

from repo_paths import ROOT, repo_path

CMAKE_PATH = repo_path("CMakeLists.txt")
CONTRACT_PATH = repo_path("docs", "gpu_module_contracts.md")
SOURCE_GROUP_RE = re.compile(r"^\s*set\s*\(\s*([A-Za-z0-9_]+_SOURCES)\b")
DOCUMENTED_GROUP_RE = re.compile(r"`([A-Z0-9_]+_SOURCES)`")


def cmake_source_groups() -> list[str]:
    groups: list[str] = []
    for line in CMAKE_PATH.read_text(encoding="utf-8").splitlines():
        match = SOURCE_GROUP_RE.match(line)
        if match:
            groups.append(match.group(1))
    return groups


def documented_source_groups() -> list[str]:
    text = CONTRACT_PATH.read_text(encoding="utf-8")
    return DOCUMENTED_GROUP_RE.findall(text)


def main() -> int:
    groups = cmake_source_groups()
    documented = documented_source_groups()

    group_set = set(groups)
    documented_set = set(documented)
    missing = sorted(group_set - documented_set)
    stale = sorted(documented_set - group_set)
    duplicate_cmake = sorted(group for group in group_set if groups.count(group) > 1)
    duplicate_docs = sorted(
        group for group in documented_set if documented.count(group) > 1
    )

    if missing or stale or duplicate_cmake or duplicate_docs:
        print(f"{CONTRACT_PATH.relative_to(ROOT)} source-group map is out of sync:")
        for group in missing:
            print(f"  missing CMake source group: {group}")
        for group in stale:
            print(f"  stale documented source group: {group}")
        for group in duplicate_cmake:
            print(f"  duplicate CMake source group in {CMAKE_PATH.name}: {group}")
        for group in duplicate_docs:
            print(f"  duplicate documented source group: {group}")
        return 1

    print(
        "GPU contract source-group check passed "
        f"({len(groups)} CMake source groups)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
