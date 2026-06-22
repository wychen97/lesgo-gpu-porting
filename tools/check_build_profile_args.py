#!/usr/bin/env python3
"""Verify documented CMake profile arguments are valid root CMake knobs."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import public_knob_names


DOC_PATH = Path("docs/build_profiles.md")
PROFILE_ARG_RE = re.compile(r"(?<![A-Za-z0-9_])-D([A-Za-z0-9_]+)(?:=|\b)")


def collect_public_cmake_knobs() -> set[str]:
    return public_knob_names()


def collect_documented_profile_args() -> set[str]:
    text = DOC_PATH.read_text(encoding="utf-8")
    return set(PROFILE_ARG_RE.findall(text))


def main() -> int:
    public_knobs = collect_public_cmake_knobs()
    documented_args = sorted(collect_documented_profile_args())
    unknown_args = [name for name in documented_args if name not in public_knobs]

    if unknown_args:
        print(f"{DOC_PATH} contains unknown CMake -D profile arguments:")
        for name in unknown_args:
            print(f"  {name}")
        print("Known root CMake knobs:")
        for name in sorted(public_knobs):
            print(f"  {name}")
        return 1

    print(
        "Build-profile argument check passed "
        f"({len(documented_args)} documented -D arguments)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
