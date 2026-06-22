#!/usr/bin/env python3
"""Verify the top-level handoff README points to the current doc entrypoints."""

from __future__ import annotations

import re
import sys
from pathlib import Path


README_PATH = Path("README_FINAL_OPTIMIZED_20260619.md")
DOCS_INDEX_PATH = Path("docs/README.md")

HANDOFF_MARKER = "Start here for collaborator handoff:"


def fenced_block_after(text: str, marker: str) -> list[str]:
    start = text.find(marker)
    if start < 0:
        raise RuntimeError(f"{README_PATH} is missing marker: {marker}")

    tail = text[start:]
    match = re.search(r"```text\n(?P<body>.*?)\n```", tail, re.DOTALL)
    if not match:
        raise RuntimeError(f"{README_PATH} is missing the handoff code block")
    return [line.strip() for line in match.group("body").splitlines() if line.strip()]


def docs_reading_order() -> list[str]:
    text = DOCS_INDEX_PATH.read_text(encoding="utf-8")
    start = text.find("## Reading Order")
    end = text.find("## Historical References")
    if start < 0 or end < start:
        raise RuntimeError(f"{DOCS_INDEX_PATH} reading-order section is malformed")

    section = text[start:end]
    docs = re.findall(r"^\d+\.\s+`([^`]+)`", section, flags=re.MULTILINE)
    return ["docs/README.md"] + [f"docs/{name}" for name in docs]


def expected_handoff_entries() -> list[str]:
    return [
        *docs_reading_order(),
        "CONTRIBUTING.md",
        "tools/README.md",
    ]


def main() -> int:
    actual = fenced_block_after(
        README_PATH.read_text(encoding="utf-8"),
        HANDOFF_MARKER,
    )
    expected = expected_handoff_entries()

    if actual != expected:
        print(f"{README_PATH} handoff list differs from {DOCS_INDEX_PATH}.")

        missing = [entry for entry in expected if entry not in actual]
        stale = [entry for entry in actual if entry not in expected]
        if missing:
            print("Missing handoff entries:")
            for entry in missing:
                print(f"  {entry}")
        if stale:
            print("Stale handoff entries:")
            for entry in stale:
                print(f"  {entry}")
        if not missing and not stale:
            print("The same entries are listed, but the order differs.")
        return 1

    print(f"Handoff README check passed ({len(expected)} entries).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
