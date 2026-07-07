#!/usr/bin/env python3
"""Verify documented Fortran source references point to tracked files."""

from __future__ import annotations

import fnmatch
import re
import sys
from pathlib import Path

from fortran_inventory import tracked_fortran_paths, tracked_paths
from repo_paths import ROOT


SOURCE_REF_RE = re.compile(r"`([^`]+\.(?:f|f90|F|F90))`")
DOC_SUFFIXES = {".md"}


def tracked_fortran_sources() -> set[Path]:
    return set(tracked_fortran_paths())


def tracked_markdown_docs() -> list[Path]:
    return [
        path
        for path in tracked_paths()
        if path.suffix in DOC_SUFFIXES
    ]


def normalize_reference(raw: str) -> str:
    return raw.replace("\\", "/").rstrip(".,:;)]}")


def reference_matches(reference: str, sources: set[Path]) -> bool:
    if "*" in reference or "?" in reference:
        pattern = reference
        if "/" in pattern:
            return any(fnmatch.fnmatch(str(path), pattern) for path in sources)
        return any(fnmatch.fnmatch(path.name, pattern) for path in sources)

    path = Path(reference)
    if "/" in reference:
        return path in sources
    return any(source.name == reference for source in sources)


def main() -> int:
    sources = tracked_fortran_sources()
    missing: list[str] = []
    checked = 0

    for doc in tracked_markdown_docs():
        text = (ROOT / doc).read_text(encoding="utf-8")
        for match in SOURCE_REF_RE.finditer(text):
            reference = normalize_reference(match.group(1))
            checked += 1
            if not reference_matches(reference, sources):
                missing.append(f"{doc}: stale Fortran reference: {reference}")

    if missing:
        print("Documented Fortran source reference check failed:")
        for item in missing:
            print(f"  {item}")
        return 1

    print(
        "Documented Fortran source reference check passed "
        f"({checked} references, {len(sources)} source files)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
