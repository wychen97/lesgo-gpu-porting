#!/usr/bin/env python3
"""Reject stale personal or ad-hoc comment fragments in source comments."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

# These are not general words to ban from documentation.  They are old source
# comment fragments that made the production code harder for collaborators to
# read because they referred to individual edits, contributor shorthand,
# misspellings, or unclear checkpoint wording instead of model behavior.
STALE_SOURCE_FRAGMENTS = (
    "Siang",
    "xiang",
    "Atharva",
    "Claude",
    "Codex",
    "collaborator's",
    "autowraped",
    "derivativex",
    "check point data",
    "check point for iwm",
    "Initialize integral wall model xiang",
)


def tracked_comment_sources() -> list[Path]:
    result = subprocess.run(
        [
            "git",
            "ls-files",
            "*.f",
            "*.F",
            "*.f90",
            "*.F90",
            "*.c",
            "*.h",
            "*.sh",
            "*.pbs",
            "*.sbatch",
            "*.cmake",
            "CMakeLists.txt",
        ],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def comment_text(path: Path, line: str) -> str:
    stripped = line.lstrip()
    if path.suffix.lower() in {".f90", ".f"}:
        if stripped.startswith("!"):
            return stripped
        if "!" in line:
            return line.split("!", 1)[1]
        return ""
    if (
        path.suffix in {".sh", ".pbs", ".sbatch", ".cmake"}
        or path.name == "CMakeLists.txt"
    ):
        if stripped.startswith("#"):
            return stripped
        return ""
    if path.suffix in {".c", ".h"}:
        if (
            stripped.startswith("//")
            or stripped.startswith("/*")
            or stripped.startswith("*")
        ):
            return stripped
    return ""


def main() -> int:
    problems: list[str] = []
    for relative in tracked_comment_sources():
        path = ROOT / relative
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(),
            start=1,
        ):
            comment = comment_text(relative, line)
            if not comment:
                continue
            for fragment in STALE_SOURCE_FRAGMENTS:
                if fragment in comment:
                    problems.append(
                        f"{relative}:{line_number}: {fragment}: {line.strip()}"
                    )

    if problems:
        print("Source comment-hygiene check failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(
        f"Source comment-hygiene check passed "
        f"({len(STALE_SOURCE_FRAGMENTS)} fragments)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
