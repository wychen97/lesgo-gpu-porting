"""Shared root CMake metadata parsers for readiness checks."""

from __future__ import annotations

import re

from repo_paths import repo_path

CMAKE_PATH = repo_path("CMakeLists.txt")

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
COMPILE_DEFINITION_RE = re.compile(r"-D(?P<name>[A-Za-z_][A-Za-z0-9_]*)\b")
SOURCE_GROUP_RE = re.compile(r"^[A-Za-z0-9_]+_SOURCES$")
VARIABLE_REF_RE = re.compile(r"\$\{(?P<name>[A-Za-z0-9_]+)\}")
LIST_APPEND_SOURCES_RE = re.compile(
    r"^\s*list\s*\(\s*APPEND\s+Sources\b(?P<rest>.*)$",
    re.IGNORECASE,
)
FORTRAN_SOURCE_RE = re.compile(r"\.(?:f90|F90|f|F)$")


def cmake_text() -> str:
    return CMAKE_PATH.read_text(encoding="utf-8")


def strip_cmake_comment(line: str) -> str:
    return line.split("#", 1)[0]


def tokenize_cmake_items(text: str) -> list[str]:
    return [
        token.strip()
        for token in re.split(r"\s+", text.replace(")", " "))
        if token.strip()
    ]


def cmake_set_blocks() -> dict[str, list[str]]:
    blocks: dict[str, list[str]] = {}
    lines = cmake_text().splitlines()
    index = 0

    while index < len(lines):
        stripped = lines[index].strip()
        match = re.match(r"^set\s*\(\s*(?P<name>[A-Za-z0-9_]+)\b(?P<rest>.*)$", stripped)
        if not match:
            index += 1
            continue

        name = match.group("name")
        pieces = [match.group("rest")]
        while ")" not in pieces[-1] and index + 1 < len(lines):
            index += 1
            pieces.append(lines[index].strip())

        block_text = " ".join(strip_cmake_comment(piece) for piece in pieces)
        blocks[name] = tokenize_cmake_items(block_text)
        index += 1

    return blocks


def cmake_source_groups() -> dict[str, list[str]]:
    return {
        name: [item for item in items if FORTRAN_SOURCE_RE.search(item)]
        for name, items in cmake_set_blocks().items()
        if SOURCE_GROUP_RE.match(name)
    }


def source_list_group_references() -> set[str]:
    set_blocks = cmake_set_blocks()
    references = {
        match.group("name")
        for token in set_blocks.get("Sources", [])
        for match in VARIABLE_REF_RE.finditer(token)
    }

    lines = cmake_text().splitlines()
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        match = LIST_APPEND_SOURCES_RE.match(stripped)
        if not match:
            index += 1
            continue

        pieces = [match.group("rest")]
        while ")" not in pieces[-1] and index + 1 < len(lines):
            index += 1
            pieces.append(lines[index].strip())

        for piece in pieces:
            references.update(
                match.group("name") for match in VARIABLE_REF_RE.finditer(piece)
            )
        index += 1

    return references


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


def compile_definition_names() -> set[str]:
    return set(COMPILE_DEFINITION_RE.findall(cmake_text()))
