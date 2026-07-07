#!/usr/bin/env python3
"""Reject stale internal optimization labels in production GPU source comments."""

from __future__ import annotations

import sys
from pathlib import Path

from report_gpu_static_inventory import GPU_MARKERS, GPU_SOURCE_FILES, tracked_fortran


ROOT = Path(__file__).resolve().parents[1]

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


def gpu_comment_source_files() -> list[Path]:
    """Return non-LVLSET root Fortran files with production GPU markers."""
    result: list[Path] = []
    for path in tracked_fortran():
        text = (ROOT / path).read_text(encoding="utf-8")
        lower_text = text.lower()
        has_gpu_marker = any(marker.lower() in lower_text for marker in GPU_MARKERS)
        if has_gpu_marker or path.name in GPU_SOURCE_FILES:
            result.append(path)
    return result


def main() -> int:
    problems: list[str] = []
    gpu_sources = gpu_comment_source_files()
    for relative in gpu_sources:
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
        f"({len(gpu_sources)} files, {len(STALE_LABELS)} labels)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
