#!/usr/bin/env python3
"""Regression tests for the Fortran interface hygiene checker."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import check_fortran_interface_hygiene as hygiene


ROOT = Path(__file__).resolve().parents[1]


def write_case(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text.strip() + "\n", encoding="utf-8")
    return path


def assert_contains(issues: list[str], needle: str) -> None:
    if not any(needle in issue for issue in issues):
        raise AssertionError(f"expected issue containing `{needle}`, got {issues!r}")


def test_derived_type_members_are_not_module_exports(tmpdir: Path) -> None:
    path = write_case(
        tmpdir,
        "type_scope.f90",
        """
        module grid_m
        implicit none
        private
        public grid
        type grid_t
            logical :: built
        contains
            procedure, public :: build
        end type grid_t
        type(grid_t) :: grid
        contains
        subroutine build(this)
            class(grid_t), intent(inout) :: this
        end subroutine build
        end module grid_m
        """,
    )

    modules = hygiene.collect_modules([path])
    module = modules["grid_m"]
    if "built" in module.local:
        raise AssertionError("derived-type component `built` leaked into module locals")
    if module.public != {"grid"}:
        raise AssertionError(f"unexpected public API for grid_m: {module.public!r}")
    issues = hygiene.check_public_api_symbols([path])
    if issues:
        raise AssertionError(f"valid derived-type case failed: {issues!r}")


def test_stale_public_symbol_is_reported(tmpdir: Path) -> None:
    path = write_case(
        tmpdir,
        "stale_public.f90",
        """
        module stale_public_m
        implicit none
        private
        public missing_symbol
        contains
        subroutine real_symbol()
        end subroutine real_symbol
        end module stale_public_m
        """,
    )

    issues = hygiene.check_public_api_symbols([path])
    assert_contains(issues, "public symbol `missing_symbol` is not defined")


def test_private_module_import_mismatch_is_reported(tmpdir: Path) -> None:
    paths = [
        write_case(
            tmpdir,
            "provider.f90",
            """
            module provider_m
            implicit none
            private
            public exported
            integer :: exported
            integer :: hidden
            end module provider_m
            """,
        ),
        write_case(
            tmpdir,
            "consumer.f90",
            """
            module consumer_m
            use provider_m, only : hidden
            implicit none
            end module consumer_m
            """,
        ),
    ]

    issues = hygiene.check_module_exports(paths)
    assert_contains(issues, "imports `hidden` from default-private module `provider_m`")


def test_use_only_reexport_mismatch_is_reported(tmpdir: Path) -> None:
    paths = [
        write_case(
            tmpdir,
            "external_like_provider.f90",
            """
            module external_like_provider_m
            implicit none
            integer :: mpi_status_size
            integer :: mpi_comm_world
            end module external_like_provider_m
            """,
        ),
        write_case(
            tmpdir,
            "wrapper.f90",
            """
            module wrapper_m
            use external_like_provider_m, only : mpi_status_size
            implicit none
            end module wrapper_m
            """,
        ),
        write_case(
            tmpdir,
            "consumer.f90",
            """
            module consumer_m
            use wrapper_m, only : mpi_comm_world
            implicit none
            end module consumer_m
            """,
        ),
    ]

    issues = hygiene.check_use_only_imports_resolve(paths)
    assert_contains(issues, "imports `mpi_comm_world` from tracked module `wrapper_m`")


def test_duplicate_use_only_local_names_are_reported(tmpdir: Path) -> None:
    path = write_case(
        tmpdir,
        "duplicate_import.f90",
        """
        module duplicate_import_m
        use provider_m, only : local_name => exported, local_name => hidden
        implicit none
        end module duplicate_import_m
        """,
    )

    issues = hygiene.check_duplicate_use_only_symbols([path])
    assert_contains(issues, "duplicate symbol(s) in use-only list: `local_name`")


def test_repeated_scope_imports_are_reported(tmpdir: Path) -> None:
    path = write_case(
        tmpdir,
        "repeated_import.f90",
        """
        module repeated_import_m
        use provider_m, only : exported
        use provider_m, only : exported, hidden
        implicit none
        end module repeated_import_m
        """,
    )

    issues = hygiene.check_repeated_scope_imports([path])
    assert_contains(issues, "imports `exported` from `provider_m` again")


def test_bare_flag_condition_is_reported(tmpdir: Path) -> None:
    path = write_case(
        tmpdir,
        "bare_condition.f90",
        """
        #if PPATM
        module bare_condition_m
        end module bare_condition_m
        #endif
        """,
    )

    issues = hygiene.check_bare_preprocessor_condition_macros([path])
    assert_contains(issues, "use `#if defined(PPATM)`")


def test_missing_module_implicit_none_is_reported(tmpdir: Path) -> None:
    path = write_case(
        tmpdir,
        "missing_implicit_none.f90",
        """
        module missing_implicit_none_m
        integer :: value
        contains
        subroutine touch()
        implicit none
        end subroutine touch
        end module missing_implicit_none_m
        """,
    )

    issues = hygiene.check_module_program_implicit_none([path])
    assert_contains(issues, "does not declare `implicit none`")


def main() -> int:
    tmpdir = Path(tempfile.mkdtemp(prefix=".fortran_hygiene_selftest_", dir=ROOT))
    try:
        test_derived_type_members_are_not_module_exports(tmpdir)
        test_stale_public_symbol_is_reported(tmpdir)
        test_private_module_import_mismatch_is_reported(tmpdir)
        test_use_only_reexport_mismatch_is_reported(tmpdir)
        test_duplicate_use_only_local_names_are_reported(tmpdir)
        test_repeated_scope_imports_are_reported(tmpdir)
        test_bare_flag_condition_is_reported(tmpdir)
        test_missing_module_implicit_none_is_reported(tmpdir)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    print("Fortran interface hygiene self-test passed (8 synthetic cases).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
