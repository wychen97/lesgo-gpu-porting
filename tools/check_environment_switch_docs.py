#!/usr/bin/env python3
"""Verify that solver `LESGO_*` environment switches are documented.

This check intentionally scans tracked Fortran source, not benchmark launch
scripts.  Test harness variables such as `LESGO_NSTEPS` or `LESGO_CMAKE_ARGS`
belong to their runner scripts; solver runtime switches belong in
`docs/environment_switches.md`.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


FORTRAN_SUFFIXES = {".f", ".f90"}
SWITCH_RE = re.compile(r"\bLESGO_[A-Z0-9_]*[A-Z0-9]\b")
TABLE_SWITCH_RE = re.compile(r"^\|\s*`(?P<name>LESGO_[A-Z0-9_]+)`\s*\|")
DOC_PATH = Path("docs/environment_switches.md")


def tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        check=True,
        stdout=subprocess.PIPE,
    )
    names = result.stdout.decode("utf-8").split("\0")
    return [Path(name) for name in names if name]


def is_fortran_source(path: Path) -> bool:
    return path.suffix.lower() in FORTRAN_SUFFIXES


def collect_source_switches(paths: list[Path]) -> dict[str, list[str]]:
    locations: dict[str, list[str]] = defaultdict(list)
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for lineno, line in enumerate(text.splitlines(), start=1):
            for match in SWITCH_RE.finditer(line):
                locations[match.group(0)].append(f"{path}:{lineno}")
    return dict(sorted(locations.items()))


def collect_documented_switches() -> set[str]:
    text = DOC_PATH.read_text(encoding="utf-8")
    return set(SWITCH_RE.findall(text))


def collect_table_switches() -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        match = TABLE_SWITCH_RE.match(line)
        if match:
            counts[match.group("name")] += 1
    return dict(sorted(counts.items()))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check that Fortran LESGO_* switches are documented."
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="print documented source switches and their source locations",
    )
    args = parser.parse_args()

    source_paths = [path for path in tracked_files() if is_fortran_source(path)]
    source_switches = collect_source_switches(source_paths)
    documented = collect_documented_switches()
    table_switches = collect_table_switches()

    missing = sorted(set(source_switches) - documented)
    missing_table_rows = sorted(set(source_switches) - set(table_switches))
    stale_table_rows = sorted(set(table_switches) - set(source_switches))
    duplicate_table_rows = sorted(
        name for name, count in table_switches.items() if count > 1
    )

    if args.list:
        for name, locations in source_switches.items():
            print(name)
            for location in locations:
                print(f"  {location}")
        return 0

    if missing:
        print(f"{DOC_PATH} is missing solver environment switches:")
        for name in missing:
            print(f"  {name}")
            for location in source_switches[name][:5]:
                print(f"    {location}")
            if len(source_switches[name]) > 5:
                print(f"    ... {len(source_switches[name]) - 5} more")
        return 1

    if missing_table_rows:
        print(f"{DOC_PATH} is missing classified table rows:")
        for name in missing_table_rows:
            print(f"  {name}")
            for location in source_switches[name][:5]:
                print(f"    {location}")
            if len(source_switches[name]) > 5:
                print(f"    ... {len(source_switches[name]) - 5} more")
        return 1

    if stale_table_rows:
        print(f"{DOC_PATH} has stale classified table rows:")
        for name in stale_table_rows:
            print(f"  {name}")
        return 1

    if duplicate_table_rows:
        print(f"{DOC_PATH} has duplicate classified table rows:")
        for name in duplicate_table_rows:
            print(f"  {name} ({table_switches[name]} rows)")
        return 1

    print(
        "Environment switch documentation check passed "
        f"({len(source_switches)} switches, {len(source_paths)} source files, "
        f"{len(table_switches)} table rows)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
