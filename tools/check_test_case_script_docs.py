#!/usr/bin/env python3
"""Verify script-heavy test-case trees have local README guidance."""

from __future__ import annotations

import sys
from pathlib import Path

from script_inventory import tracked_test_case_files, tracked_test_case_script_dirs


GUIDANCE_SECTION_MARKERS = (
    "## Launcher Groups",
    "## Files",
    "## Recommended sequence",
)
LAUNCH_ARTIFACT_MARKERS = (
    ".pbs",
    ".sbatch",
    ".sh",
    "qsub ",
    "bash ",
)


def readme_has_launcher_guidance(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    return (
        any(marker in text for marker in GUIDANCE_SECTION_MARKERS)
        and any(marker in text for marker in LAUNCH_ARTIFACT_MARKERS)
    )


def main() -> int:
    files = tracked_test_case_files()
    script_case_dirs = tracked_test_case_script_dirs()

    missing = sorted(
        case_dir / "README.md"
        for case_dir in script_case_dirs
        if case_dir / "README.md" not in files
    )
    weak = sorted(
        readme
        for readme in (case_dir / "README.md" for case_dir in script_case_dirs)
        if readme in files and not readme_has_launcher_guidance(readme)
    )

    if missing or weak:
        print("Test-case script documentation check failed:")
        for path in missing:
            print(f"  missing README for script-heavy case tree: {path}")
        for path in weak:
            print(f"  README lacks launcher guidance markers: {path}")
        return 1

    print(
        f"Test-case script documentation check passed "
        f"({len(script_case_dirs)} case trees)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
