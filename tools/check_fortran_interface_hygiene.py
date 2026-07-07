#!/usr/bin/env python3
"""Check Fortran preprocessor symbols and module export boundaries.

This catches several classes of errors that are easy to miss in a large optional
Fortran codebase:

* stale ``#ifdef`` names that are no longer emitted by CMake;
* bare single-symbol ``#if``/``#elif`` conditions that should use
  ``defined(SYMBOL)`` for flag-style build macros;
* module/program scopes that do not declare ``implicit none``;
* ``use module, only: symbol`` imports from default-private modules where the
  imported symbol is not explicitly public;
* duplicate symbols in one ``use ..., only:`` import list;
* repeated imports of the same local name from the same module in one scope;
* stale or duplicate names in explicit ``public`` API lists;
* ``use module, only: symbol`` imports that are not actually exported by a
  tracked module, including symbols that were previously re-exported through a
  broad module import.
* ``use mpi, only:`` imports of MPI procedures that are not portable across
  the Derecho NVHPC/Cray MPI module interface;
* direct calls to specific ``test_filtermodule`` procedures without importing
  those specific procedure names.
"""

from __future__ import annotations

import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CMAKE_PATH = ROOT / "CMakeLists.txt"
ALLOWED_COMPILER_MACROS = {"__INTEL_COMPILER", "PPXLF"}
MAX_ISSUES = 200
MPI_PROCEDURE_NAMES = {
    "mpi_abort",
    "mpi_allreduce",
    "mpi_barrier",
    "mpi_cart_coords",
    "mpi_cart_create",
    "mpi_cart_rank",
    "mpi_cart_shift",
    "mpi_comm_rank",
    "mpi_comm_size",
    "mpi_comm_split",
    "mpi_get_processor_name",
    "mpi_init",
    "mpi_intercomm_create",
    "mpi_irecv",
    "mpi_isend",
    "mpi_recv",
    "mpi_reduce",
    "mpi_send",
    "mpi_sendrecv",
    "mpi_wait",
    "mpi_waitall",
    "mpi_wtime",
}
TEST_FILTER_SPECIFIC_PROCEDURES = {
    "test_filter_3",
    "test_filter_6",
    "test_filter_plane_gpu",
    "test_test_filter_3",
    "test_test_filter_6",
    "test_test_filter_plane_gpu",
}


@dataclass
class ModuleInfo:
    path: Path
    line: int
    private_default: bool = False
    public: set[str] = field(default_factory=set)
    public_lines: list[tuple[Path, int, list[str]]] = field(default_factory=list)
    private: set[str] = field(default_factory=set)
    local: set[str] = field(default_factory=set)
    use_only: set[str] = field(default_factory=set)
    use_all: list[str] = field(default_factory=list)


def tracked_fortran_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "*.f90", "*.F90", "*.f", "*.F"],
        check=True,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [ROOT / line for line in result.stdout.splitlines() if line]


def strip_comment(line: str) -> str:
    chars: list[str] = []
    in_single = False
    in_double = False
    for char in line:
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        if char == "!" and not in_single and not in_double:
            break
        chars.append(char)
    return "".join(chars)


def logical_lines(path: Path) -> list[tuple[int, str]]:
    rows: list[tuple[int, str]] = []
    buffer = ""
    start_line = 0

    for line_no, raw_line in enumerate(path.read_text(errors="ignore").splitlines(), 1):
        stripped = raw_line.lstrip()
        if stripped.startswith("#"):
            if buffer:
                rows.append((start_line, buffer.strip()))
                buffer = ""
                start_line = 0
            rows.append((line_no, stripped))
            continue

        line = strip_comment(raw_line).rstrip()
        if not line.strip():
            continue

        if not buffer:
            start_line = line_no

        continued = line.endswith("&")
        if continued:
            line = line[:-1]
        if buffer and line.lstrip().startswith("&"):
            line = line.lstrip()[1:]

        buffer = f"{buffer} {line.strip()}" if buffer else line.strip()
        if not continued:
            rows.append((start_line, buffer.strip()))
            buffer = ""
            start_line = 0

    if buffer:
        rows.append((start_line, buffer.strip()))
    return rows


def split_names(text: str) -> list[str]:
    names: list[str] = []
    token: list[str] = []
    depth = 0
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")" and depth:
            depth -= 1

        if char == "," and depth == 0:
            name = "".join(token).strip()
            if name:
                names.append(name)
            token = []
        else:
            token.append(char)

    name = "".join(token).strip()
    if name:
        names.append(name)
    return names


def clean_symbol(text: str) -> str:
    if "=>" in text:
        text = text.split("=>", 1)[1]
    text = text.split("=", 1)[0].strip().lower()
    return re.split(r"\s|\(", text)[0]


def local_symbol(text: str) -> str:
    if "=>" in text:
        text = text.split("=>", 1)[0]
    text = text.split("=", 1)[0].strip().lower()
    return re.split(r"\s|\(", text)[0]


def is_special_public_name(text: str) -> bool:
    text = text.strip().lower()
    return text.startswith(("operator", "assignment"))


def declaration_symbols(line: str) -> list[str]:
    if "::" not in line:
        return []
    left, right = line.split("::", 1)
    if re.match(r"\s*(use|public|private)\b", left, flags=re.I):
        return []
    declaration_type = re.search(
        r"\b(integer|real|logical|character|complex|procedure)\b",
        left,
        flags=re.I,
    ) or re.search(r"\b(?:type|class)\s*(?:\(|\b)", left, flags=re.I)
    if not declaration_type:
        return []
    return [local_symbol(name) for name in split_names(right) if local_symbol(name)]


def strip_cmake_comment(line: str) -> str:
    return line.split("#", 1)[0]


def cmake_macros() -> set[str]:
    text = CMAKE_PATH.read_text(encoding="utf-8")
    macros: set[str] = set()
    for match in re.finditer(r"add_definitions\s*\(([^)]*)\)", text, flags=re.I):
        body = "\n".join(strip_cmake_comment(line) for line in match.group(1).splitlines())
        macros.update(re.findall(r"-D([A-Z_][A-Z0-9_]*)", body))
    for match in re.finditer(r"add_compile_definitions\s*\(([^)]*)\)", text, flags=re.I):
        body = "\n".join(strip_cmake_comment(line) for line in match.group(1).splitlines())
        macros.update(re.findall(r"\b([A-Z_][A-Z0-9_]*)\b", body))
    return macros


def check_preprocessor_symbols(paths: list[Path]) -> list[str]:
    defined = cmake_macros() | ALLOWED_COMPILER_MACROS
    issues: list[str] = []
    symbol_re = re.compile(r"(?<![A-Za-z0-9_])([A-Z_][A-Z0-9_]{2,})(?![A-Za-z0-9_])")

    for path in paths:
        for line_no, line in logical_lines(path):
            if not re.match(r"^#\s*(if|ifdef|ifndef|elif)\b", line):
                continue
            for symbol in symbol_re.findall(line):
                if symbol not in defined:
                    rel = path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: preprocessor symbol `{symbol}` is "
                        "not defined by CMake or the compiler allow-list"
                    )
    return issues


def check_bare_preprocessor_condition_macros(paths: list[Path]) -> list[str]:
    issues: list[str] = []
    bare_condition_re = re.compile(r"^#\s*(if|elif)\s+([A-Z_][A-Z0-9_]*)\s*$")

    for path in paths:
        for line_no, line in logical_lines(path):
            match = bare_condition_re.match(line)
            if not match:
                continue
            rel = path.relative_to(ROOT)
            keyword, symbol = match.groups()
            issues.append(
                f"{rel}:{line_no}: use `#{keyword} defined({symbol})` for "
                f"flag-style macro `{symbol}`"
            )
    return issues


def check_unused_cmake_preprocessor_symbols(paths: list[Path]) -> list[str]:
    source_text = "\n".join(path.read_text(errors="ignore") for path in paths)
    issues: list[str] = []

    for symbol in sorted(cmake_macros()):
        if not symbol.startswith("PP"):
            continue
        pattern = re.compile(
            rf"(?<![A-Za-z0-9_]){re.escape(symbol)}(?![A-Za-z0-9_])"
        )
        if not pattern.search(source_text):
            issues.append(
                f"{CMAKE_PATH.relative_to(ROOT)}: CMake defines `{symbol}`, "
                "but no tracked Fortran source references it"
            )
    return issues


def check_module_program_implicit_none(paths: list[Path]) -> list[str]:
    issues: list[str] = []

    for path in paths:
        scope: tuple[str, str, int] | None = None
        seen_implicit_none = False

        for line_no, line in logical_lines(path):
            lower = line.lower()
            match = re.match(r"^(module|program)\s+(?!procedure\b)(\w+)", lower)
            if match:
                scope = (match.group(1), match.group(2), line_no)
                seen_implicit_none = False
                continue

            if not scope:
                continue

            if re.match(r"^implicit\s+none\b", lower):
                seen_implicit_none = True
                continue

            if lower == "contains" or re.match(r"^end\s+(module|program)\b", lower):
                if not seen_implicit_none:
                    rel = path.relative_to(ROOT)
                    scope_type, scope_name, start_line = scope
                    issues.append(
                        f"{rel}:{start_line}: {scope_type} `{scope_name}` "
                        "does not declare `implicit none`"
                    )
                scope = None
                seen_implicit_none = False

    return issues


def collect_modules(paths: list[Path]) -> dict[str, ModuleInfo]:
    modules: dict[str, ModuleInfo] = {}

    for path in paths:
        current: ModuleInfo | None = None
        in_contains = False
        type_depth = 0
        for line_no, line in logical_lines(path):
            lower = line.lower()
            match = re.match(r"^module\s+(?!procedure\b)(\w+)", lower)
            if match:
                current = ModuleInfo(path=path, line=line_no)
                modules[match.group(1)] = current
                in_contains = False
                type_depth = 0
                continue

            if current and re.match(r"^end\s+module\b", lower):
                current = None
                in_contains = False
                type_depth = 0
                continue

            if not current:
                continue

            if not in_contains and re.match(r"^end\s+type\b", lower):
                type_depth = max(0, type_depth - 1)

            if lower == "contains" and type_depth == 0:
                in_contains = True
                continue

            module_scope = not in_contains and type_depth == 0

            if module_scope:
                use_match = re.match(
                    r"^use\s*(?:,\s*[^:]+::\s*)?(\w+)(?:\s*,\s*only\s*:\s*(.+))?$",
                    line,
                    flags=re.I,
                )
                if use_match:
                    module_name = use_match.group(1).lower()
                    import_list = use_match.group(2)
                    if import_list is None:
                        current.use_all.append(module_name)
                    else:
                        for name in split_names(import_list):
                            if not is_special_public_name(name):
                                symbol = local_symbol(name)
                                if symbol:
                                    current.use_only.add(symbol)

                if re.match(r"^private\s*(::\s*)?$", lower):
                    current.private_default = True
                elif re.match(r"^private\b", lower):
                    rest = re.sub(r"^private\b\s*(::)?\s*", "", line, flags=re.I).strip()
                    for name in split_names(rest):
                        symbol = local_symbol(name)
                        if symbol and not is_special_public_name(name):
                            current.private.add(symbol)

                if re.match(r"^public\s*(::\s*)?$", lower):
                    current.private_default = False

                if re.match(r"^public\b", lower):
                    rest = re.sub(r"^public\b\s*(::)?\s*", "", line, flags=re.I).strip()
                    line_symbols: list[str] = []
                    for name in split_names(rest):
                        symbol = local_symbol(name)
                        if symbol and not is_special_public_name(name):
                            current.public.add(symbol)
                            line_symbols.append(symbol)
                    if line_symbols:
                        current.public_lines.append((path, line_no, line_symbols))

                current.local.update(declaration_symbols(line))

                type_match = re.match(r"^type\b(?:\s*,[^:]*)?\s*::\s*(\w+)", lower)
                if type_match:
                    current.local.add(type_match.group(1))
                    type_depth += 1
                legacy_type_match = re.match(r"^type\s+(\w+)", lower)
                if legacy_type_match and not lower.startswith("type("):
                    current.local.add(legacy_type_match.group(1))
                    type_depth += 1
                interface_match = re.match(r"^interface\s+(\w+)", lower)
                if interface_match:
                    current.local.add(interface_match.group(1))

            function_match = re.match(
                r"^(?:(?:recursive|pure|elemental|impure)\s+)*"
                r"(?:(?:integer|real|logical|complex|character|type|class)"
                r"\s*(?:\([^)]*\))?\s+)?(subroutine|function)\s+(\w+)",
                lower,
            )
            if function_match:
                current.local.add(function_match.group(2))

            if (
                module_scope
                and "public" in lower
                and "::" in lower
                and not re.match(r"^public\b", lower)
            ):
                left, right = line.split("::", 1)
                if re.search(r"\bpublic\b", left, flags=re.I):
                    line_symbols: list[str] = []
                    for name in split_names(right):
                        symbol = local_symbol(name)
                        if symbol:
                            current.public.add(symbol)
                            line_symbols.append(symbol)
                    if line_symbols:
                        current.public_lines.append((path, line_no, line_symbols))
                if re.search(r"\bprivate\b", left, flags=re.I):
                    current.private.update(declaration_symbols(line))

    return modules


def exported_module_symbols(
    modules: dict[str, ModuleInfo],
    module_name: str,
    memo: dict[str, set[str]] | None = None,
) -> set[str]:
    if memo is None:
        memo = {}
    if module_name in memo:
        return set(memo[module_name])
    module = modules.get(module_name)
    if not module:
        return set()

    memo[module_name] = set()
    available = available_module_symbols(modules, module_name, memo)

    if module.private_default:
        exported = module.public & available
    else:
        exported = available
    exported -= module.private
    memo[module_name] = set(exported)
    return exported


def available_module_symbols(
    modules: dict[str, ModuleInfo],
    module_name: str,
    memo: dict[str, set[str]] | None = None,
) -> set[str]:
    if memo is None:
        memo = {}
    module = modules.get(module_name)
    if not module:
        return set()

    available = set(module.local) | set(module.use_only)
    for dependency in module.use_all:
        available |= exported_module_symbols(modules, dependency, memo)
    available -= module.private
    return available


def check_module_exports(paths: list[Path]) -> list[str]:
    modules = collect_modules(paths)
    issues: list[str] = []

    for path in paths:
        for line_no, line in logical_lines(path):
            match = re.match(
                r"^use\s*(?:,\s*[^:]+::\s*)?(\w+)\s*,\s*only\s*:\s*(.+)$",
                line,
                flags=re.I,
            )
            if not match:
                continue

            module_name = match.group(1).lower()
            imported_names = match.group(2)
            module = modules.get(module_name)
            if not module or not module.private_default:
                continue

            for name in split_names(imported_names):
                lower = name.lower().strip()
                if not lower or lower.startswith(("operator", "assignment")):
                    continue
                symbol = clean_symbol(name)
                if symbol and symbol not in module.public:
                    rel = path.relative_to(ROOT)
                    module_rel = module.path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: imports `{symbol}` from default-private "
                        f"module `{module_name}`, but `{module_rel}` does not mark it public"
                    )
    return issues


def check_use_only_imports_resolve(paths: list[Path]) -> list[str]:
    modules = collect_modules(paths)
    issues: list[str] = []
    exported_cache: dict[str, set[str]] = {}

    for path in paths:
        for line_no, line in logical_lines(path):
            match = re.match(
                r"^use\s*(?:,\s*[^:]+::\s*)?(\w+)\s*,\s*only\s*:\s*(.+)$",
                line,
                flags=re.I,
            )
            if not match:
                continue

            module_name = match.group(1).lower()
            module = modules.get(module_name)
            if not module:
                continue

            exported = exported_module_symbols(modules, module_name, exported_cache)
            for name in split_names(match.group(2)):
                lower = name.lower().strip()
                if not lower or lower.startswith(("operator", "assignment")):
                    continue
                symbol = clean_symbol(name)
                if symbol and symbol not in exported:
                    rel = path.relative_to(ROOT)
                    module_rel = module.path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: imports `{symbol}` from tracked "
                        f"module `{module_name}`, but `{module_rel}` does not "
                        "export that symbol"
                    )
    return issues


def check_duplicate_use_only_symbols(paths: list[Path]) -> list[str]:
    issues: list[str] = []

    for path in paths:
        for line_no, line in logical_lines(path):
            match = re.match(
                r"^use\s*(?:,\s*[^:]+::\s*)?(\w+)\s*,\s*only\s*:\s*(.+)$",
                line,
                flags=re.I,
            )
            if not match:
                continue

            seen: set[str] = set()
            duplicates: list[str] = []
            for name in split_names(match.group(2)):
                lower = name.lower().strip()
                if not lower or lower.startswith(("operator", "assignment")):
                    continue
                symbol = local_symbol(name)
                if not symbol:
                    continue
                if symbol in seen and symbol not in duplicates:
                    duplicates.append(symbol)
                seen.add(symbol)

            if duplicates:
                rel = path.relative_to(ROOT)
                issues.append(
                    f"{rel}:{line_no}: duplicate symbol(s) in use-only list: "
                    + ", ".join(f"`{symbol}`" for symbol in duplicates)
                )
    return issues


def check_repeated_scope_imports(paths: list[Path]) -> list[str]:
    issues: list[str] = []

    for path in paths:
        scope: list[tuple[str, str, int]] = []
        imports: dict[tuple[str, str], int] = {}

        for line_no, line in logical_lines(path):
            lower = line.lower()
            module_match = re.match(r"^module\s+(?!procedure\b)(\w+)", lower)
            if module_match:
                scope = [("module", module_match.group(1), line_no)]
                imports = {}
                continue

            if re.match(r"^end\s+module\b", lower):
                scope = []
                imports = {}
                continue

            procedure_match = re.match(
                r"^(?:(?:recursive|pure|elemental|impure)\s+)*"
                r"(?:(?:integer|real|logical|complex|character|type|class)"
                r"\s*(?:\([^)]*\))?\s+)?(subroutine|function)\s+(\w+)",
                lower,
            )
            if procedure_match:
                scope.append((procedure_match.group(1), procedure_match.group(2), line_no))
                imports = {}
                continue

            if re.match(r"^end\s+(subroutine|function)\b", lower):
                if len(scope) > 1:
                    scope.pop()
                imports = {}
                continue

            use_match = re.match(
                r"^use\s*(?:,\s*[^:]+::\s*)?(\w+)\s*,\s*only\s*:\s*(.+)$",
                line,
                flags=re.I,
            )
            if not use_match or not scope:
                continue

            module_name = use_match.group(1).lower()
            for name in split_names(use_match.group(2)):
                if is_special_public_name(name):
                    continue
                symbol = local_symbol(name)
                if not symbol:
                    continue
                key = (module_name, symbol)
                first_line = imports.get(key)
                if first_line is not None and first_line != line_no:
                    rel = path.relative_to(ROOT)
                    scope_name = "::".join(item[1] for item in scope)
                    issues.append(
                        f"{rel}:{line_no}: `{scope_name}` imports `{symbol}` "
                        f"from `{module_name}` again; first import is line {first_line}"
                    )
                else:
                    imports[key] = line_no

    return issues


def check_mpi_procedure_only_imports(paths: list[Path]) -> list[str]:
    issues: list[str] = []

    for path in paths:
        for line_no, line in logical_lines(path):
            match = re.match(
                r"^use\s*(?:,\s*[^:]+::\s*)?mpi\s*,\s*only\s*:\s*(.+)$",
                line,
                flags=re.I,
            )
            if not match:
                continue

            procedures = sorted(
                {
                    clean_symbol(name)
                    for name in split_names(match.group(1))
                    if clean_symbol(name) in MPI_PROCEDURE_NAMES
                }
            )
            if procedures:
                rel = path.relative_to(ROOT)
                issues.append(
                    f"{rel}:{line_no}: imports MPI procedure(s) with `use mpi, only:` "
                    + ", ".join(f"`{procedure}`" for procedure in procedures)
                    + "; use broad `use mpi` for Derecho NVHPC/Cray MPI portability"
                )

    return issues


def check_test_filter_specific_imports(paths: list[Path]) -> list[str]:
    issues: list[str] = []
    call_re = re.compile(
        r"^\s*call\s+("
        + "|".join(re.escape(name) for name in sorted(TEST_FILTER_SPECIFIC_PROCEDURES))
        + r")\b",
        flags=re.I,
    )

    for path in paths:
        scope_name = "<file>"
        imports: set[str] = set()
        module_imports: set[str] = set()
        in_module = False
        in_module_contains = False

        for line_no, line in logical_lines(path):
            lower = line.lower()
            module_match = re.match(r"^module\s+(?!procedure\b)(\w+)", lower)
            if module_match:
                scope_name = "<file>"
                imports = set()
                module_imports = set()
                in_module = True
                in_module_contains = False
                continue

            if re.match(r"^end\s+module\b", lower):
                scope_name = "<file>"
                imports = set()
                module_imports = set()
                in_module = False
                in_module_contains = False
                continue

            if in_module and lower == "contains":
                in_module_contains = True
                continue

            procedure_match = re.match(
                r"^(?:(?:recursive|pure|elemental|impure)\s+)*"
                r"(?:(?:integer|real|logical|complex|character|type|class)"
                r"\s*(?:\([^)]*\))?\s+)?(subroutine|function)\s+(\w+)",
                lower,
            )
            if procedure_match:
                scope_name = procedure_match.group(2)
                imports = set(module_imports) if in_module else set()
                continue

            use_match = re.match(
                r"^use\s*(?:,\s*[^:]+::\s*)?test_filtermodule"
                r"(?:\s*,\s*only\s*:\s*(.+))?$",
                line,
                flags=re.I,
            )
            if use_match:
                import_list = use_match.group(1)
                if import_list is None:
                    imports.add("*")
                else:
                    imports.update(
                        local_symbol(name)
                        for name in split_names(import_list)
                        if local_symbol(name)
                    )
                if in_module and not in_module_contains and scope_name == "<file>":
                    module_imports.update(imports)
                continue

            call_match = call_re.match(line)
            if call_match:
                called = call_match.group(1).lower()
                if "*" not in imports and called not in imports:
                    rel = path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: `{scope_name}` calls `{called}` but "
                        "does not import that specific `test_filtermodule` procedure"
                    )

            if re.match(r"^end\s+(subroutine|function)\b", lower):
                scope_name = "<file>"
                imports = set()

    return issues


def check_public_api_symbols(paths: list[Path]) -> list[str]:
    modules = collect_modules(paths)
    issues: list[str] = []

    for module_name, module in modules.items():
        available = available_module_symbols(modules, module_name)

        module_seen: set[str] = set()
        for path, line_no, names in module.public_lines:
            line_seen: set[str] = set()
            for name in names:
                if name in line_seen:
                    rel = path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: duplicate symbol `{name}` in public list"
                    )
                line_seen.add(name)

                if name in module_seen:
                    rel = path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: symbol `{name}` is listed public more than once"
                    )
                module_seen.add(name)

                if name not in available:
                    rel = path.relative_to(ROOT)
                    issues.append(
                        f"{rel}:{line_no}: public symbol `{name}` is not defined "
                        f"or imported by module `{module_name}`"
                    )

    return issues


def main() -> int:
    paths = tracked_fortran_files()
    issues = (
        check_preprocessor_symbols(paths)
        + check_bare_preprocessor_condition_macros(paths)
        + check_unused_cmake_preprocessor_symbols(paths)
        + check_module_program_implicit_none(paths)
        + check_module_exports(paths)
        + check_use_only_imports_resolve(paths)
        + check_duplicate_use_only_symbols(paths)
        + check_repeated_scope_imports(paths)
        + check_mpi_procedure_only_imports(paths)
        + check_test_filter_specific_imports(paths)
        + check_public_api_symbols(paths)
    )

    if issues:
        print("Fortran interface hygiene check failed:")
        for issue in issues[:MAX_ISSUES]:
            print(f"  {issue}")
        if len(issues) > MAX_ISSUES:
            print(f"  ... {len(issues) - MAX_ISSUES} more issues omitted")
        return 1

    print(
        "Fortran interface hygiene check passed "
        f"({len(paths)} files, no stale/unused macros, missing module/program "
        "implicit-none statements, private-export mismatches, bare flag conditions, "
        "unresolved tracked-module imports, duplicate/repeated imports, unsafe MPI "
        "procedure-only imports, missing test-filter procedure bindings, or stale "
        "public API names)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
