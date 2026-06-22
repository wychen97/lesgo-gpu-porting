#!/usr/bin/env python3
"""Verify root Fortran sources are covered by named CMake source groups."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CMAKE_PATH = ROOT / "CMakeLists.txt"
FORTRAN_SUFFIXES = {".f", ".f90", ".F", ".F90"}
SOURCE_GROUP_RE = re.compile(r"^[A-Za-z0-9_]+_SOURCES$")
VARIABLE_REF_RE = re.compile(r"\$\{(?P<name>[A-Za-z0-9_]+)\}")
LIST_APPEND_SOURCES_RE = re.compile(
    r"^\s*list\s*\(\s*APPEND\s+Sources\b(?P<rest>.*)$",
    re.IGNORECASE,
)


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


def tokenize_cmake_items(text: str) -> list[str]:
    return [
        token.strip()
        for token in re.split(r"\s+", text.replace(")", " "))
        if token.strip()
    ]


def cmake_set_blocks() -> dict[str, list[str]]:
    blocks: dict[str, list[str]] = {}
    lines = CMAKE_PATH.read_text(encoding="utf-8").splitlines()
    index = 0

    while index < len(lines):
        stripped = lines[index].strip()
        match = re.match(r"^set\s*\(\s*(?P<name>[A-Za-z0-9_]+)\b(?P<rest>.*)$", stripped)
        if not match:
            index += 1
            continue

        name = match.group("name")
        pieces = [match.group("rest")]
        while ")" not in pieces[-1] and index + 1 < len(lines):
            index += 1
            pieces.append(lines[index].strip())

        blocks[name] = tokenize_cmake_items(" ".join(pieces))
        index += 1

    return blocks


def source_list_group_references(set_blocks: dict[str, list[str]]) -> set[str]:
    references = {
        match.group("name")
        for token in set_blocks.get("Sources", [])
        for match in VARIABLE_REF_RE.finditer(token)
    }

    lines = CMAKE_PATH.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        match = LIST_APPEND_SOURCES_RE.match(stripped)
        if not match:
            index += 1
            continue

        pieces = [match.group("rest")]
        while ")" not in pieces[-1] and index + 1 < len(lines):
            index += 1
            pieces.append(lines[index].strip())

        for piece in pieces:
            references.update(
                match.group("name") for match in VARIABLE_REF_RE.finditer(piece)
            )
        index += 1

    return references


def main() -> int:
    tracked = tracked_root_fortran_sources()
    set_blocks = cmake_set_blocks()
    source_groups = {
        name: [item for item in items if Path(item).suffix in FORTRAN_SUFFIXES]
        for name, items in set_blocks.items()
        if SOURCE_GROUP_RE.match(name)
    }

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

    source_references = source_list_group_references(set_blocks)
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
