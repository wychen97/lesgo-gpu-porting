#!/usr/bin/env python3
"""Verify contributor/readiness docs point to the readiness source of truth."""

from __future__ import annotations

import sys
from pathlib import Path


REQUIRED_SNIPPETS_BY_PATH = {
    Path("CONTRIBUTING.md"): (
        "python3 tools/check_branch_readiness.py",
        "python3 tools/check_branch_readiness.py --with-cmake-configure",
        "python3 tools/check_branch_readiness.py --with-hit-cmake-configure",
        'README_FINAL_OPTIMIZED_20260619.md` under "Local Readability Checks"',
        "`tools/README.md` is the index for the individual readiness tools",
    ),
    Path("README_FINAL_OPTIMIZED_20260619.md"): (
        "python3 tools/check_branch_readiness.py --with-hit-cmake-configure",
    ),
    Path("docs/architecture_hardening_audit.md"): (
        "python3 tools/check_branch_readiness.py --with-git-diff-check --with-cmake-configure --with-hit-cmake-configure",
    ),
    Path("docs/collaboration_readiness_status.md"): (
        "python3 tools/check_branch_readiness.py --with-hit-cmake-configure",
    ),
    Path("docs/collaborator_quickstart.md"): (
        "python3 tools/check_branch_readiness.py --with-hit-cmake-configure",
    ),
    Path("docs/gpu_development_guidelines.md"): (
        "python3 tools/check_branch_readiness.py --with-hit-cmake-configure",
    ),
}


def main() -> int:
    missing_by_path: dict[Path, list[str]] = {}
    snippet_count = 0
    for path, snippets in REQUIRED_SNIPPETS_BY_PATH.items():
        text = path.read_text(encoding="utf-8")
        snippet_count += len(snippets)
        missing = [snippet for snippet in snippets if snippet not in text]
        if missing:
            missing_by_path[path] = missing

    if missing_by_path:
        print("Readiness documentation is missing required guidance:")
        for path, missing in missing_by_path.items():
            print(f"  {path}:")
            for snippet in missing:
                print(f"    {snippet}")
        return 1

    print(
        "Contributing readiness-doc check passed "
        f"({snippet_count} snippets across {len(REQUIRED_SNIPPETS_BY_PATH)} files)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
