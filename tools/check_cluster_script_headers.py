#!/usr/bin/env python3
"""Verify root-level shell helpers carry a visible Status header."""

from __future__ import annotations

import sys
from pathlib import Path

from script_inventory import tracked_root_shell_scripts


HEADER_LINES = 12


def has_status_header(path: Path) -> bool:
    lines = path.read_text(encoding="utf-8").splitlines()[:HEADER_LINES]
    return any(line.startswith("# Status: ") for line in lines)


def main() -> int:
    scripts = tracked_root_shell_scripts()
    missing = [path for path in scripts if not has_status_header(path)]

    if missing:
        print("Root-level shell helpers are missing a '# Status:' header:")
        for path in missing:
            print(f"  {path}")
        return 1

    print(f"Cluster script header check passed ({len(scripts)} scripts).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
