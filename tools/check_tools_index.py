#!/usr/bin/env python3
"""Verify tools/README.md indexes all tracked top-level tool files."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


INDEX_PATH = Path("tools/README.md")
TOOL_ROW_RE = re.compile(r"^\|\s*`(?P<path>tools/[^`]+)`\s*\|")


def tracked_tool_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "tools"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [
        Path(line)
        for line in result.stdout.splitlines()
        if line and Path(line).parent == Path("tools") and Path(line) != INDEX_PATH
    ]


def documented_tool_rows() -> dict[Path, int]:
    rows: dict[Path, int] = {}
    for line in INDEX_PATH.read_text(encoding="utf-8").splitlines():
        match = TOOL_ROW_RE.match(line)
        if not match:
            continue
        path = Path(match.group("path"))
        rows[path] = rows.get(path, 0) + 1
    return rows


def main() -> int:
    files = tracked_tool_files()
    documented_rows = documented_tool_rows()
    documented = set(documented_rows)

    missing = sorted(set(files) - documented)
    stale = sorted(documented - set(files))
    duplicates = sorted(
        path for path, count in documented_rows.items() if count > 1
    )

    if missing or stale or duplicates:
        print(f"{INDEX_PATH} tooling table is out of sync:")
        for path in missing:
            print(f"  missing tracked top-level tool file: {path}")
        for path in stale:
            print(f"  stale tool row: {path}")
        for path in duplicates:
            print(f"  duplicate tool row: {path} ({documented_rows[path]} rows)")
        return 1

    print(f"Tooling index check passed ({len(files)} tool files).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
