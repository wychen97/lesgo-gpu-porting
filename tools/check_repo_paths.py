#!/usr/bin/env python3
"""Verify tools use the shared repository path helper."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from repo_paths import ROOT


ALLOWED_DIRECT_ROOT_FILES = {Path("tools/repo_paths.py")}
DIRECT_ROOT_EXPRESSION = "Path(__file__).resolve().parents" + "[1]"


def tracked_tool_python_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "tools/*.py"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def main() -> int:
    problems: list[str] = []
    for relative in tracked_tool_python_files():
        text = (ROOT / relative).read_text(encoding="utf-8")
        if DIRECT_ROOT_EXPRESSION not in text:
            continue
        if relative in ALLOWED_DIRECT_ROOT_FILES:
            continue
        problems.append(f"{relative}: use tools/repo_paths.py instead")

    if problems:
        print("Repository path-helper check failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(
        "Repository path-helper check passed "
        f"({len(tracked_tool_python_files())} Python tools)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
