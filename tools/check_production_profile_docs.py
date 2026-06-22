#!/usr/bin/env python3
"""Verify production USE_* settings agree across collaborator docs."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import public_knob_values


ORG_PATH = Path("docs/code_organization.md")
PROFILE_PATH = Path("docs/build_profiles.md")
README_PATH = Path("README_FINAL_OPTIMIZED_20260619.md")

USE_ASSIGNMENT_RE = re.compile(r"(?:-D)?(USE_[A-Za-z0-9_]+)=([A-Za-z0-9_]+)")


def fenced_block_after(text: str, marker: str) -> str:
    start = text.find(marker)
    if start < 0:
        raise RuntimeError(f"missing marker: {marker}")

    tail = text[start:]
    match = re.search(r"```(?:[A-Za-z0-9_+-]+)?\n(?P<body>.*?)\n```", tail, re.DOTALL)
    if not match:
        raise RuntimeError(f"missing fenced block after marker: {marker}")
    return match.group("body")


def collect_use_assignments(block: str) -> dict[str, str]:
    return {name: value for name, value in USE_ASSIGNMENT_RE.findall(block)}


def cmake_use_values() -> dict[str, set[str]]:
    return public_knob_values("USE_")


def production_doc_settings() -> dict[Path, dict[str, str]]:
    org_text = ORG_PATH.read_text(encoding="utf-8")
    profile_text = PROFILE_PATH.read_text(encoding="utf-8")
    readme_text = README_PATH.read_text(encoding="utf-8")

    return {
        ORG_PATH: collect_use_assignments(
            fenced_block_after(org_text, "The validated production path is:")
        ),
        PROFILE_PATH: collect_use_assignments(
            fenced_block_after(
                profile_text,
                "Use this profile for the current optimized non-LVLSET production path:",
            )
        ),
        README_PATH: collect_use_assignments(
            fenced_block_after(readme_text, "## Main Build Configuration")
        ),
    }


def main() -> int:
    settings_by_path = production_doc_settings()
    reference_path = ORG_PATH
    reference = settings_by_path[reference_path]
    valid_cmake_values = cmake_use_values()

    mismatches: list[str] = []
    for path, settings in settings_by_path.items():
        for name in sorted(set(reference) | set(settings)):
            reference_value = reference.get(name, "<missing>")
            value = settings.get(name, "<missing>")
            if value != reference_value:
                mismatches.append(
                    f"{name}: {reference_path}={reference_value}, {path}={value}"
                )

    invalid_settings: list[str] = []
    for name, value in sorted(reference.items()):
        allowed_values = valid_cmake_values.get(name)
        if allowed_values is None:
            invalid_settings.append(f"{name}: not a root CMake USE_* knob")
        elif allowed_values and value not in allowed_values:
            invalid_settings.append(
                f"{name}: value {value!r} is not one of "
                f"{', '.join(sorted(allowed_values))}"
            )

    if mismatches or invalid_settings:
        print("Production USE_* settings differ across collaborator docs.")
        for mismatch in mismatches:
            print(f"  {mismatch}")
        for invalid in invalid_settings:
            print(f"  invalid production setting: {invalid}")
        return 1

    print(
        "Production profile documentation check passed "
        f"({len(reference)} USE_* settings, {len(settings_by_path)} docs)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
