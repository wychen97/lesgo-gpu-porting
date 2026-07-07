#!/usr/bin/env python3
"""Verify MPI-only imports and sync wrappers stay out of serial compile paths."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from fortran_inventory import ROOT, normalized_repo_path, tracked_fortran_files


MPI_REQUIRED_FILES = {
    "concurrent_precursor.f90",
    "cuda_mpi_debug.f90",
    "mpi_defs.f90",
    "mpi_transpose_mod.f90",
}

SYNC_SYMBOL_RE = re.compile(r"\b(mpi_sync_real_array|mpi_sync_down\w*)\b", re.I)
RAW_MPI_IMPORT_RE = re.compile(r"^\s*use\s+mpi(?:\s*,|\s*$)", re.I)
MPI_CALL_RE = re.compile(r"\bcall\s+mpi_\w+", re.I)
POSITIVE_PPMPI_GUARDS = {"PPMPI", "PPGPU_AWARE_MPI"}


def positive_guard_macros(line: str) -> set[str]:
    stripped = line.strip()
    if stripped.startswith("#ifdef"):
        parts = stripped.split()
        return {parts[1]} if len(parts) >= 2 else set()
    if stripped.startswith("#ifndef"):
        return set()
    if stripped.startswith(("#if", "#elif")):
        result: set[str] = set()
        for macro in POSITIVE_PPMPI_GUARDS:
            if re.search(rf"\bdefined\s*\(\s*{macro}\s*\)", stripped):
                result.add(macro)
            elif re.search(rf"\bdefined\s+{macro}\b", stripped):
                result.add(macro)
        return result
    return set()


def check_file(path: Path) -> list[str]:
    if path.name in MPI_REQUIRED_FILES:
        return []

    rel = normalized_repo_path(path.relative_to(ROOT))
    issues: list[str] = []
    guard_stack: list[set[str]] = []

    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith(("#if", "#ifdef", "#ifndef")):
            guard_stack.append(positive_guard_macros(stripped))
            continue
        if stripped.startswith("#elif"):
            if guard_stack:
                guard_stack.pop()
            guard_stack.append(positive_guard_macros(stripped))
            continue
        if stripped.startswith("#else"):
            if guard_stack:
                guard_stack.pop()
            guard_stack.append(set())
            continue
        if stripped.startswith("#endif"):
            if guard_stack:
                guard_stack.pop()
            continue

        code = line.split("!", 1)[0]
        guarded = any(macros & POSITIVE_PPMPI_GUARDS for macros in guard_stack)
        if RAW_MPI_IMPORT_RE.search(code) and not guarded:
            issues.append(f"{rel}:{line_no}: raw `use mpi` must be guarded by PPMPI")
            continue
        if MPI_CALL_RE.search(code) and not guarded:
            issues.append(f"{rel}:{line_no}: direct MPI call must be guarded by PPMPI")
            continue
        if not SYNC_SYMBOL_RE.search(code):
            continue
        if guarded:
            continue
        issues.append(
            f"{rel}:{line_no}: MPI sync wrapper use must be guarded by PPMPI"
        )

    return issues


def main() -> int:
    issues: list[str] = []
    for path in tracked_fortran_files():
        issues.extend(check_file(path))

    if issues:
        print("MPI sync guard hygiene check failed:")
        for issue in issues:
            print(f"  {issue}")
        return 1

    print("MPI sync guard hygiene check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
