#!/usr/bin/env python3
"""Verify tracked readiness checks are included in the readiness wrapper."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from readiness_manifest import WRAPPER_RELATIVE_PATH, wrapper_script_paths

WRAPPER_ONLY = {WRAPPER_RELATIVE_PATH}


def tracked_readiness_scripts() -> set[Path]:
    result = subprocess.run(
        ["git", "ls-files", "tools/" + "check_" + "*.py"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return {Path(line) for line in result.stdout.splitlines() if line}


def wrapper_scripts() -> set[Path]:
    return set(wrapper_script_paths())


def main() -> int:
    tracked = tracked_readiness_scripts()
    expected = tracked - WRAPPER_ONLY
    configured = wrapper_scripts()

    missing = sorted(expected - configured)
    stale = sorted(configured - tracked)

    if missing or stale:
        print(f"{WRAPPER_RELATIVE_PATH} readiness coverage check failed.")
        if missing:
            print("Tracked readiness scripts missing from PYTHON_CHECKS:")
            for path in missing:
                print(f"  {path}")
        if stale:
            print("PYTHON_CHECKS entries that are not tracked readiness scripts:")
            for path in stale:
                print(f"  {path}")
        return 1

    print(f"Readiness coverage check passed ({len(configured)} scripts in wrapper).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
