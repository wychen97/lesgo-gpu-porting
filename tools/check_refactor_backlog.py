#!/usr/bin/env python3
"""Verify the refactor backlog lists the current source-size hotspots."""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

from fortran_inventory import tracked_fortran_files


ROOT = Path(__file__).resolve().parents[1]
BACKLOG_PATH = ROOT / "docs/refactor_backlog.md"
HOTSPOT_COUNT = 8

# LVLSET is tracked but outside the optimized production path for this branch.
EXCLUDED_PREFIXES = ("trees_",)
EXCLUDED_NAMES = {
    "level_set.f90",
    "level_set_base.f90",
}

TABLE_ROW_RE = re.compile(
    r"^\|\s*(?P<priority>\d+)\s*"
    r"\|\s*`(?P<path>[^`]+)`\s*"
    r"\|\s*(?P<size>[0-9,]+)\s+bytes\s*\|"
)


@dataclass(frozen=True)
class Hotspot:
    path: Path
    size: int


def is_excluded(path: Path) -> bool:
    name = path.name
    if name in EXCLUDED_NAMES:
        return True
    return any(name.startswith(prefix) for prefix in EXCLUDED_PREFIXES)


def expected_hotspots() -> list[Hotspot]:
    candidates = [
        Hotspot(path.relative_to(ROOT), path.stat().st_size)
        for path in tracked_fortran_files()
        if not is_excluded(path)
    ]
    candidates.sort(key=lambda item: (-item.size, str(item.path)))
    return candidates[:HOTSPOT_COUNT]


def documented_hotspots() -> list[Hotspot]:
    hotspots: list[Hotspot] = []
    for line in BACKLOG_PATH.read_text(encoding="utf-8").splitlines():
        match = TABLE_ROW_RE.match(line)
        if not match:
            continue
        size = int(match.group("size").replace(",", ""))
        hotspots.append(Hotspot(Path(match.group("path")), size))
    return hotspots


def format_hotspots(hotspots: list[Hotspot]) -> str:
    rows = []
    for priority, hotspot in enumerate(hotspots, start=1):
        rows.append(f"  {priority}. {hotspot.path} ({hotspot.size:,} bytes)")
    return "\n".join(rows)


def main() -> int:
    expected = expected_hotspots()
    documented = documented_hotspots()

    if documented != expected:
        print(f"{BACKLOG_PATH.relative_to(ROOT)} has a stale hotspot table.")
        print("Expected current top non-LVLSET Fortran hotspots:")
        print(format_hotspots(expected))
        print("Documented:")
        print(format_hotspots(documented))
        return 1

    print(f"Refactor-backlog hotspot check passed ({len(expected)} files).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
