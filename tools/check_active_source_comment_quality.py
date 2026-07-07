#!/usr/bin/env python3
"""Reject vague maintenance comments in active non-LVLSET Fortran source."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from repo_paths import ROOT

# LVLSET/tree sources are intentionally deferred from the optimized production
# path.  Their legacy comments are tracked separately by docs/refactor_backlog.md.
DEFERRED_PREFIXES = ("level_set", "trees_")

VAGUE_COMMENT_FRAGMENTS = (
    "for now",
    "temporary",
    "fix later",
    "try with",
    "blows up",
    "can remove after testing",
    "kludgy",
    "todo",
    "fixme",
    "hack",
    "some people",
    "any better ideas",
    "could maybe",
    "could arguably",
    "not sure",
    "should be ok",
    "bad results",
)


def tracked_fortran_sources() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "*.f90", "*.F90"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    paths = [Path(line) for line in result.stdout.splitlines() if line]
    return [
        path
        for path in paths
        if not path.name.startswith(DEFERRED_PREFIXES)
    ]


def comment_text(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("!"):
        return stripped
    if "!" in line:
        return line.split("!", 1)[1]
    return ""


def main() -> int:
    problems: list[str] = []
    for relative in tracked_fortran_sources():
        path = ROOT / relative
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(),
            start=1,
        ):
            comment = comment_text(line).lower()
            if not comment:
                continue
            for fragment in VAGUE_COMMENT_FRAGMENTS:
                if fragment in comment:
                    problems.append(
                        f"{relative}:{line_number}: {fragment}: {line.strip()}"
                    )

    if problems:
        print("Active source comment-quality check failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(
        f"Active source comment-quality check passed "
        f"({len(VAGUE_COMMENT_FRAGMENTS)} fragments)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
