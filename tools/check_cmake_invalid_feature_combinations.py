#!/usr/bin/env python3
"""Verify invalid public CMake feature combinations fail at configure time."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CASES = [
    (
        "CPS without MPI",
        [
            "-DUSE_CPS=ON",
            "-DUSE_MPI=OFF",
        ],
        "USE_CPS requires USE_MPI=ON",
    ),
    (
        "CPU baseline mixed with LES GPU",
        [
            "-DUSE_CPU_BUILD=ON",
            "-DUSE_LES_GPU=ON",
        ],
        "USE_CPU_BUILD requires USE_LES_GPU=OFF",
    ),
    (
        "scalar GPU without scalar transport",
        [
            "-DUSE_SCALARS_GPU=ON",
            "-DUSE_SCALARS=OFF",
            "-DUSE_LES_GPU=ON",
        ],
        "USE_SCALARS_GPU requires USE_SCALARS=ON",
    ),
    (
        "scalar GPU without LES GPU",
        [
            "-DUSE_SCALARS_GPU=ON",
            "-DUSE_SCALARS=ON",
            "-DUSE_LES_GPU=OFF",
        ],
        "USE_SCALARS_GPU requires USE_LES_GPU=ON",
    ),
    (
        "forced GPU-aware MPI without MPI",
        [
            "-DUSE_GPU_AWARE_MPI=ON",
            "-DUSE_MPI=OFF",
            "-DUSE_LES_GPU=ON",
        ],
        "USE_GPU_AWARE_MPI=ON requires USE_MPI=ON",
    ),
    (
        "forced GPU-aware MPI without LES GPU",
        [
            "-DUSE_GPU_AWARE_MPI=ON",
            "-DUSE_MPI=ON",
            "-DUSE_LES_GPU=OFF",
        ],
        "USE_GPU_AWARE_MPI=ON requires USE_LES_GPU=ON",
    ),
]


def run_invalid_case(label: str, args: list[str], expected: str) -> bool:
    with tempfile.TemporaryDirectory(prefix="lesgo_invalid_cmake_") as build_dir:
        result = subprocess.run(
            ["cmake", "-S", str(ROOT), "-B", build_dir, *args],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

    output = result.stdout or ""
    if result.returncode == 0:
        print(f"{label}: configure unexpectedly succeeded")
        return False
    if expected not in output:
        print(f"{label}: expected diagnostic not found")
        print(f"Expected substring: {expected}")
        print("Observed output tail:")
        print(output[-4000:].rstrip())
        return False
    return True


def main() -> int:
    ok = True
    for label, args, expected in CASES:
        ok = run_invalid_case(label, args, expected) and ok

    if not ok:
        return 1

    print(
        "CMake invalid feature-combination check passed "
        f"({len(CASES)} rejected combinations)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
