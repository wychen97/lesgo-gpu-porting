#!/usr/bin/env python3
"""Regression-test the broad Fortran import report parser."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from report_fortran_broad_imports import (
    BroadImport,
    ROOT,
    classify_broad_import,
    collect_broad_imports,
)


SYNTHETIC_SOURCE = """\
module broad_import_demo
use iso_c_binding
use, intrinsic :: iso_fortran_env, only : int64
use param, only : nx, ny
implicit none
contains
subroutine scoped_imports()
use mpi
use mpi, only : MPI_STATUS_SIZE
use, intrinsic :: iso_c_binding, only : c_int, c_ptr
implicit none
end subroutine scoped_imports
end module broad_import_demo
"""


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="lesgo_broad_import_selftest_") as tmp:
        path = Path(tmp) / "synthetic_use_imports.f90"
        path.write_text(SYNTHETIC_SOURCE, encoding="utf-8")

        observed = {
            (row.scope, row.module) for row in collect_broad_imports([path])
        }

    expected = {
        ("broad_import_demo", "iso_c_binding"),
        ("broad_import_demo::scoped_imports", "mpi"),
    }

    if observed != expected:
        print("Broad-import parser self-test failed.")
        print("Expected broad imports:")
        for row in sorted(expected):
            print(f"  {row[0]} uses {row[1]}")
        print("Observed broad imports:")
        for row in sorted(observed):
            print(f"  {row[0]} uses {row[1]}")
        return 1

    classification_cases = {
        "MPI compiler interface": BroadImport(ROOT / "main.f90", 64, "main", "mpi"),
        "optional CGNS interface": BroadImport(ROOT / "io.f90", 46, "io", "cgns"),
        "Level Set optional profile": BroadImport(
            ROOT / "level_set.f90", 32, "level_set", "level_set_base"
        ),
        "candidate": BroadImport(ROOT / "forcing.f90", 999, "forcing::demo", "param"),
    }
    for expected_category, row in classification_cases.items():
        observed_category = classify_broad_import(row).category
        if observed_category != expected_category:
            print("Broad-import classification self-test failed.")
            print(f"Expected {expected_category!r} for {row}.")
            print(f"Observed {observed_category!r}.")
            return 1

    print(
        "Fortran broad-import parser/classification self-test passed "
        "(2 synthetic broad imports, 4 classification cases)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
