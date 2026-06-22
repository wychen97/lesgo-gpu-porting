"""Shared tracked Fortran-source inventory helpers for readiness checks."""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORTRAN_SUFFIXES = {".f", ".f90", ".F", ".F90"}


def tracked_paths(*patterns: str) -> list[Path]:
    command = ["git", "ls-files", "-z"]
    command.extend(patterns)
    result = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    names = result.stdout.decode("utf-8").split("\0")
    return [Path(name) for name in names if name]


def tracked_fortran_paths() -> list[Path]:
    return [
        path
        for path in tracked_paths()
        if path.suffix in FORTRAN_SUFFIXES
    ]


def tracked_fortran_files() -> list[Path]:
    return [ROOT / path for path in tracked_fortran_paths()]


def normalized_repo_path(path: Path) -> str:
    return str(path).replace("\\", "/")
