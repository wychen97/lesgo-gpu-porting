#!/usr/bin/env python3
"""Verify docs/README.md indexes all collaborator-facing docs."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


INDEX_PATH = Path("docs/README.md")
LOCAL_DOC_RE = re.compile(r"`(?P<name>[A-Za-z0-9_-]+\.md)`")


def tracked_docs() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "docs/*.md"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def documented_local_docs() -> dict[str, int]:
    rows: dict[str, int] = {}
    text = INDEX_PATH.read_text(encoding="utf-8")
    for match in LOCAL_DOC_RE.finditer(text):
        name = match.group("name")
        rows[name] = rows.get(name, 0) + 1
    return rows


def main() -> int:
    docs = [path for path in tracked_docs() if path != INDEX_PATH]
    tracked_names = {path.name for path in docs}
    documented_rows = documented_local_docs()
    documented_names = set(documented_rows)

    missing = sorted(tracked_names - documented_names)
    stale = sorted(documented_names - tracked_names)
    duplicates = sorted(
        name for name, count in documented_rows.items() if count > 1
    )

    if missing or stale or duplicates:
        print(f"{INDEX_PATH} documentation index is out of sync:")
        for path in missing:
            print(f"  missing tracked documentation file: docs/{path}")
        for path in stale:
            print(f"  stale documentation entry: {path}")
        for path in duplicates:
            print(f"  duplicate documentation entry: {path} ({documented_rows[path]} rows)")
        return 1

    print(f"Documentation index check passed ({len(docs)} indexed docs).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
