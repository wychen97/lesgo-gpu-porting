#!/usr/bin/env python3
"""Reject vague maintenance comments in CMake build files."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

VAGUE_COMMENT_FRAGMENTS = (
    "for now",
    "temporary",
    "fix later",
    "not sure",
    "maybe",
    "todo",
    "fixme",
    "hack",
)


def tracked_cmake_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "CMakeLists.txt", "*.cmake"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def comment_text(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return stripped
    if "#" in line:
        return line.split("#", 1)[1]
    return ""


def main() -> int:
    problems: list[str] = []
    for relative in tracked_cmake_files():
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
        print("CMake comment-quality check failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(
        f"CMake comment-quality check passed "
        f"({len(VAGUE_COMMENT_FRAGMENTS)} fragments)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
