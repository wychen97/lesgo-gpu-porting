"""Shared root CMake metadata parsers for readiness checks."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CMAKE_PATH = ROOT / "CMakeLists.txt"

OPTION_NAME_RE = re.compile(
    r"^\s*option\s*\(\s*(?P<name>[A-Za-z0-9_]+)\b",
    re.MULTILINE,
)
CACHE_SET_NAME_RE = re.compile(
    r"^\s*set\s*\(\s*(?P<name>[A-Za-z0-9_]+)\b[^\n]*\bCACHE\b",
    re.MULTILINE,
)
CACHE_STRINGS_RE = re.compile(
    r"^\s*set_property\s*\(\s*CACHE\s+(?P<name>[A-Za-z0-9_]+)\s+"
    r"PROPERTY\s+STRINGS\s+(?P<values>[A-Za-z0-9_ ]+)\)",
    re.MULTILINE,
)
OPTION_DESCRIPTION_RE = re.compile(
    r"option\(\s*(?P<name>[A-Za-z0-9_]+)\s+"
    r'"(?P<description>[^"]*)"',
    re.DOTALL,
)
CACHE_DESCRIPTION_RE = re.compile(
    r"set\(\s*(?P<name>[A-Za-z0-9_]+)\s+"
    r'"[^"]*"\s+CACHE\s+STRING\s+"(?P<description>[^"]*)"',
    re.DOTALL,
)


def cmake_text() -> str:
    return CMAKE_PATH.read_text(encoding="utf-8")


def with_prefix(names: set[str], prefix: str | None) -> set[str]:
    if prefix is None:
        return names
    return {name for name in names if name.startswith(prefix)}


def option_names(prefix: str | None = None) -> set[str]:
    return with_prefix(set(OPTION_NAME_RE.findall(cmake_text())), prefix)


def cache_variable_names(prefix: str | None = None) -> set[str]:
    return with_prefix(set(CACHE_SET_NAME_RE.findall(cmake_text())), prefix)


def public_knob_names(prefix: str | None = None) -> set[str]:
    return option_names(prefix) | cache_variable_names(prefix)


def cache_string_values(prefix: str | None = None) -> dict[str, set[str]]:
    values = {
        match.group("name"): set(match.group("values").split())
        for match in CACHE_STRINGS_RE.finditer(cmake_text())
    }
    if prefix is None:
        return values
    return {name: value for name, value in values.items() if name.startswith(prefix)}


def public_knob_values(prefix: str | None = None) -> dict[str, set[str]]:
    values: dict[str, set[str]] = {
        name: {"ON", "OFF"} for name in option_names(prefix)
    }
    for name in cache_variable_names(prefix):
        values.setdefault(name, set())
    values.update(cache_string_values(prefix))
    return values


def cmake_descriptions(prefix: str | None = None) -> dict[str, str]:
    descriptions: dict[str, str] = {}
    for regex in (OPTION_DESCRIPTION_RE, CACHE_DESCRIPTION_RE):
        for match in regex.finditer(cmake_text()):
            name = match.group("name")
            if prefix is not None and not name.startswith(prefix):
                continue
            descriptions[name] = " ".join(match.group("description").split())
    return descriptions
