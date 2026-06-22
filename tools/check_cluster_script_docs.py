#!/usr/bin/env python3
"""Verify root-level cluster/helper shell scripts are documented."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from script_inventory import tracked_root_shell_scripts


DOC_PATH = Path("docs/cluster_scripts.md")
SCRIPT_ROW_RE = re.compile(r"^\|\s*`(?P<path>[^`]+\.sh)`\s*\|")


def documented_script_rows() -> dict[Path, int]:
    rows: dict[Path, int] = {}
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        match = SCRIPT_ROW_RE.match(line)
        if not match:
            continue
        path = Path(match.group("path"))
        rows[path] = rows.get(path, 0) + 1
    return rows


def main() -> int:
    scripts = tracked_root_shell_scripts()
    documented_rows = documented_script_rows()
    documented = set(documented_rows)

    missing = sorted(set(scripts) - documented)
    stale = sorted(documented - set(scripts))
    duplicates = sorted(
        path for path, count in documented_rows.items() if count > 1
    )

    if missing or stale or duplicates:
        print(f"{DOC_PATH} root-level shell script table is out of sync:")
        for path in missing:
            print(f"  missing root-level shell script: {path}")
        for path in stale:
            print(f"  stale shell script row: {path}")
        for path in duplicates:
            print(f"  duplicate shell script row: {path} ({documented_rows[path]} rows)")
        return 1

    print(f"Cluster script documentation check passed ({len(scripts)} scripts).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
