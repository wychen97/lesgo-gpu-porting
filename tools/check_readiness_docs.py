#!/usr/bin/env python3
"""Verify the README lists the readiness-wrapper Python checks."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from readiness_manifest import WRAPPER_RELATIVE_PATH, wrapper_checks

README_PATH = Path("README_FINAL_OPTIMIZED_20260619.md")

CHECKLIST_MARKER = "The readiness wrapper runs these portable checks:"
SCRIPT_LINE_RE = re.compile(r"^\s*python3\s+(tools/[A-Za-z0-9_./-]+\.py)\s*$")


def collect_wrapper_checks() -> list[str]:
    return [script for _label, script in wrapper_checks()]


def collect_readme_checklist() -> list[str]:
    text = README_PATH.read_text(encoding="utf-8")
    marker_index = text.find(CHECKLIST_MARKER)
    if marker_index < 0:
        raise RuntimeError(f"{README_PATH} is missing the readiness checklist marker")

    tail = text[marker_index:]
    match = re.search(r"```bash\n(?P<body>.*?)\n```", tail, flags=re.DOTALL)
    if not match:
        raise RuntimeError(f"{README_PATH} is missing the readiness checklist code block")

    scripts: list[str] = []
    for line in match.group("body").splitlines():
        script_match = SCRIPT_LINE_RE.match(line)
        if script_match:
            scripts.append(script_match.group(1))
    return scripts


def main() -> int:
    wrapper_checks = collect_wrapper_checks()
    readme_checks = collect_readme_checklist()

    if readme_checks != wrapper_checks:
        print(f"{README_PATH} readiness checklist differs from {WRAPPER_RELATIVE_PATH}.")

        missing = [script for script in wrapper_checks if script not in readme_checks]
        stale = [script for script in readme_checks if script not in wrapper_checks]
        if missing:
            print("Missing from README checklist:")
            for script in missing:
                print(f"  {script}")
        if stale:
            print("Stale entries in README checklist:")
            for script in stale:
                print(f"  {script}")

        if not missing and not stale:
            print("The same scripts are listed, but the order differs.")
        return 1

    print(f"Readiness documentation check passed ({len(wrapper_checks)} scripts).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
