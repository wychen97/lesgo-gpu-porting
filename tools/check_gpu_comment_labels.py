#!/usr/bin/env python3
"""Reject stale internal optimization labels in production GPU source comments."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

GPU_SOURCE_FILES = (
    "convec_gpu.f90",
    "derivatives_gpu.f90",
    "fft_gpu.f90",
    "lagrange_Sdep_gpu.f90",
    "press_gpu.f90",
    "sgs_gpu.f90",
    "tridag_gpu.f90",
    "turbines_gpu.f90",
)

STALE_LABELS = (
    "P0a",
    "P0b",
    "P2.3",
    "P3:",
    "post-P1",
    "Phase P1",
    "HOST BOUNCE",
    "Phase 1 of this port",
    "Phase 2 of this port",
    "left on the host for now",
    "will move it to GPU",
)


def main() -> int:
    problems: list[str] = []
    for relative in GPU_SOURCE_FILES:
        path = ROOT / relative
        text = path.read_text(encoding="utf-8")
        for line_number, line in enumerate(text.splitlines(), start=1):
            for label in STALE_LABELS:
                if label in line:
                    problems.append(f"{relative}:{line_number}: {label}: {line.strip()}")

    if problems:
        print("GPU comment-label check failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(
        f"GPU comment-label check passed "
        f"({len(GPU_SOURCE_FILES)} files, {len(STALE_LABELS)} labels)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
