#!/usr/bin/env python3
"""Verify IWM point-local wall-surface arrays keep `(iwm_i, iwm_j)` semantics."""

from __future__ import annotations

import re
import sys
from pathlib import Path


SOURCE_PATH = Path("iwmles.f90")

# These arrays are allocated over the wall surface.  Inside loops over
# iwm_i/iwm_j, a write or read at (iwm_i, iwm_i) is almost certainly a typo:
# it updates only the diagonal when nx == ny and can be out-of-bounds otherwise.
POINT_LOCAL_2D_ARRAYS = (
    "iwm_flt_us",
    "iwm_tR",
    "iwm_tauwx",
    "iwm_tauwy",
    "iwm_utx",
    "iwm_uty",
    "iwm_Dz",
    "iwm_z0",
    "iwm_Ax",
    "iwm_Ay",
    "iwm_u_inst",
    "iwm_v_inst",
    "iwm_w_inst",
    "iwm_p_inst",
)

DIAGONAL_INDEX_RE = re.compile(
    r"\b(?P<array>"
    + "|".join(re.escape(name) for name in POINT_LOCAL_2D_ARRAYS)
    + r")\s*\(\s*iwm_i\s*,\s*iwm_i\b"
)


def source_without_comments(line: str) -> str:
    """Drop normal Fortran comments; OpenACC/directive comments are code-like."""
    stripped = line.lstrip()
    if stripped.startswith("!$"):
        return line
    return line.split("!", 1)[0]


def continued_statement_chunks() -> list[tuple[int, str]]:
    chunks: list[tuple[int, str]] = []
    start_line = 0
    pieces: list[str] = []

    for line_number, line in enumerate(
        SOURCE_PATH.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        code = source_without_comments(line).strip()
        if not code:
            continue

        if not pieces:
            start_line = line_number

        if code.startswith("&"):
            code = code[1:].lstrip()

        continued = code.endswith("&")
        if continued:
            code = code[:-1].rstrip()
        pieces.append(code)

        if not continued:
            chunks.append((start_line, " ".join(pieces)))
            pieces = []

    if pieces:
        chunks.append((start_line, " ".join(pieces)))

    return chunks


def main() -> int:
    violations: list[tuple[int, str, str]] = []
    for line_number, statement in continued_statement_chunks():
        match = DIAGONAL_INDEX_RE.search(statement)
        if match:
            violations.append((line_number, match.group("array"), statement))

    if violations:
        print(f"{SOURCE_PATH} has diagonal IWM wall-surface indexing:")
        for line_number, array, line in violations:
            print(f"  line {line_number}: {array} should use (iwm_i, iwm_j): {line}")
        return 1

    print(
        f"IWM wall-surface indexing check passed "
        f"({len(POINT_LOCAL_2D_ARRAYS)} arrays)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
