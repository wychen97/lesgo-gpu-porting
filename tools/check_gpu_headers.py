#!/usr/bin/env python3
"""Verify GPU-specific Fortran files explain their ownership or contract."""

from __future__ import annotations

import sys
from pathlib import Path

from fortran_inventory import tracked_fortran_files


ROOT = Path(__file__).resolve().parents[1]
HEADER_SCAN_LINES = 120
GPU_FILE_MARKERS = ("_gpu", "gpu.")
REQUIRED_MARKERS = (
    "GPU implementation",
    "GPU helper module",
    "GPU FFT layer",
    "GPU port",
    "Device residency contract",
    "Data movement contract",
    "Ownership map",
    "Navigation map",
    "Routine map",
)


def is_gpu_file(path: Path) -> bool:
    lower = path.name.lower()
    return any(marker in lower for marker in GPU_FILE_MARKERS)


def header_text(path: Path) -> str:
    lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for _ in range(HEADER_SCAN_LINES):
            line = handle.readline()
            if not line:
                break
            lines.append(line)
    return "".join(lines)


def main() -> int:
    missing: list[str] = []
    checked = 0

    for path in tracked_fortran_files():
        if not is_gpu_file(path):
            continue
        checked += 1
        header = header_text(path)
        if not any(marker in header for marker in REQUIRED_MARKERS):
            missing.append(str(path.relative_to(ROOT)))

    if missing:
        print("GPU header check failed:")
        for item in missing:
            print(f"  missing ownership/contract header: {item}")
        return 1

    print(f"GPU header check passed ({checked} GPU Fortran files).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
