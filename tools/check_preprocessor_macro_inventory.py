#!/usr/bin/env python3
"""Verify active Fortran preprocessor guards use known build macros."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from cmake_metadata import cmake_source_groups, compile_definition_names
from fortran_inventory import tracked_fortran_paths
from repo_paths import ROOT

DIRECTIVE_RE = re.compile(
    r"^\s*#\s*(?P<directive>if|ifdef|ifndef|elif)\b(?P<body>.*)$"
)
DIRECT_MACRO_RE = re.compile(r"^\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)")
DEFINED_RE = re.compile(
    r"\bdefined\s*(?:\(\s*)?(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
)
BARE_MACRO_RE = re.compile(r"\b[A-Z_][A-Z0-9_]*\b")

COMPILER_PROVIDED_MACROS = {
    "__INTEL_COMPILER",
}

BANNED_LEGACY_MACROS = {
    "ENABLE_CUDA",
}

MAX_ISSUES = 100
DEFERRED_SOURCE_GROUPS = {
    "LVLSET_SOURCES",
}


def deferred_fortran_paths() -> set[Path]:
    return {
        Path(source)
        for group, sources in cmake_source_groups().items()
        if group in DEFERRED_SOURCE_GROUPS
        for source in sources
    }


def macros_in_preprocessor_line(line: str) -> set[str]:
    match = DIRECTIVE_RE.match(line)
    if not match:
        return set()

    directive = match.group("directive")
    body = match.group("body")
    if directive in {"ifdef", "ifndef"}:
        direct = DIRECT_MACRO_RE.match(body)
        return {direct.group("name")} if direct else set()

    macros = {item.group("name") for item in DEFINED_RE.finditer(body)}
    macros.update(
        item.group(0)
        for item in BARE_MACRO_RE.finditer(body)
        if item.group(0) != "DEFINED"
    )
    return macros


def used_source_macros() -> dict[str, list[str]]:
    macros: dict[str, list[str]] = {}
    deferred_paths = deferred_fortran_paths()
    for relative in tracked_fortran_paths():
        if relative in deferred_paths:
            continue
        path = ROOT / relative
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for macro in macros_in_preprocessor_line(line):
                macros.setdefault(macro, []).append(f"{relative}:{lineno}")
    return macros


def main() -> int:
    defined_by_cmake = compile_definition_names()
    allowed = defined_by_cmake | COMPILER_PROVIDED_MACROS
    source_macros = used_source_macros()

    issues: list[str] = []
    for macro in sorted(source_macros):
        locations = ", ".join(source_macros[macro][:3])
        if macro in BANNED_LEGACY_MACROS:
            issues.append(f"{macro}: banned legacy macro at {locations}")
            continue
        if macro not in allowed:
            issues.append(f"{macro}: not defined by root CMake at {locations}")

    if issues:
        print("Preprocessor macro inventory check failed:")
        for issue in issues[:MAX_ISSUES]:
            print(f"  {issue}")
        if len(issues) > MAX_ISSUES:
            print(f"  ... stopped after {MAX_ISSUES} issues")
        return 1

    print(
        "Preprocessor macro inventory check passed "
        f"({len(source_macros)} source macros, "
        f"{len(defined_by_cmake)} CMake definitions)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
