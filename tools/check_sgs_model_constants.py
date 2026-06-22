#!/usr/bin/env python3
"""Verify active SGS dispatch uses named model constants."""

from __future__ import annotations

import re
import sys
from pathlib import Path


SOURCE_PATHS = (
    Path("initialize.f90"),
    Path("main.f90"),
    Path("sgs_gpu.f90"),
    Path("sgs_param.f90"),
    Path("sgs_stag_util.f90"),
)

RAW_SGS_MODEL_RE = re.compile(
    r"\bsgs_model\s*"
    r"(?P<op>==|/=|>=|<=|>|<|\.eq\.|\.ne\.|\.ge\.|\.le\.|\.gt\.|\.lt\.)"
    r"\s*(?P<value>[1-9])(?:\.|\b)",
    flags=re.IGNORECASE,
)


def source_without_comments(line: str) -> str:
    """Drop normal Fortran comments; OpenACC/directive comments are code-like."""
    stripped = line.lstrip()
    if stripped.startswith("!$"):
        return line
    return line.split("!", 1)[0]


def continued_statement_chunks(path: Path) -> list[tuple[int, str]]:
    chunks: list[tuple[int, str]] = []
    start_line = 0
    pieces: list[str] = []

    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
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
    violations: list[tuple[Path, int, str]] = []

    for path in SOURCE_PATHS:
        for line_number, statement in continued_statement_chunks(path):
            if RAW_SGS_MODEL_RE.search(statement):
                violations.append((path, line_number, statement))

    if violations:
        print("Active SGS dispatch has raw numeric sgs_model comparisons.")
        print("Use SGS_MODEL_* constants from sgs_param instead:")
        for path, line_number, statement in violations:
            print(f"  {path}:{line_number}: {statement}")
        return 1

    print(
        "SGS model-constant check passed "
        f"({len(SOURCE_PATHS)} active source files)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
