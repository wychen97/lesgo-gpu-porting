#!/usr/bin/env python3
"""Verify validated production paths are not described as experimental."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import cmake_descriptions


ROOT = Path(__file__).resolve().parents[1]
CODE_ORG = ROOT / "docs/code_organization.md"

# These source comments describe validated production paths.  Experimental
# labels are still allowed for real experiments such as point-owner ATM LB.
SOURCE_FORBIDDEN = {
    "scalars.f90": ("experimental scalar path", "experimental scalar GPU path"),
    "test_filtermodule.f90": ("experimental LES GPU route",),
}
DOC_FORBIDDEN = {
    "test-cases/large_windfarm_3072x384x400_60turbines/README.md": (
        "gpu-explicit-residency-wip",
        "explicit-residency WIP branch",
        "not fully managed-memory-free",
        "remaining modules are converted",
    ),
}


def production_enabled_options() -> set[str]:
    text = CODE_ORG.read_text(encoding="utf-8")
    match = re.search(
        r"The validated production path is:\n\n```text\n(?P<body>.*?)\n```",
        text,
        re.DOTALL,
    )
    if not match:
        raise RuntimeError("Could not find production configuration block")

    options: set[str] = set()
    for line in match.group("body").splitlines():
        if "=" not in line:
            continue
        name, value = [part.strip() for part in line.split("=", 1)]
        if name.startswith("USE_") and value in {"ON", "AUTO"}:
            options.add(name)
    return options


def forbidden_phrase_problems(
    paths_to_phrases: dict[str, tuple[str, ...]],
    message: str,
) -> list[str]:
    problems: list[str] = []
    for relative, forbidden_phrases in paths_to_phrases.items():
        text = (ROOT / relative).read_text(encoding="utf-8")
        lower_text = text.lower()
        for phrase in forbidden_phrases:
            if phrase.lower() in lower_text:
                problems.append(f"{relative} {message}: {phrase!r}")
    return problems


def main() -> int:
    problems: list[str] = []
    production_options = production_enabled_options()
    descriptions = cmake_descriptions()

    for name in sorted(production_options):
        description = descriptions.get(name)
        if description and "experimental" in description.lower():
            problems.append(
                f"CMake description for validated production option {name} "
                f"contains 'experimental': {description}"
            )

    problems.extend(
        forbidden_phrase_problems(
            SOURCE_FORBIDDEN,
            "labels a validated production path as experimental",
        )
    )
    problems.extend(
        forbidden_phrase_problems(
            DOC_FORBIDDEN,
            "contains stale final-branch wording",
        )
    )

    if problems:
        print("Production wording check failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(
        "Production wording check passed "
        f"({len(production_options)} production USE_* settings)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
