#!/usr/bin/env python3
"""Verify the generated broad-import audit is current."""

from __future__ import annotations

import sys

import report_fortran_broad_imports as broad_imports


def main() -> int:
    rows = broad_imports.collect_broad_imports(
        broad_imports.hygiene.tracked_fortran_files()
    )
    expected = broad_imports.markdown_report(rows)
    current = broad_imports.REPORT_PATH.read_text(encoding="utf-8")
    if current != expected:
        print(f"{broad_imports.REPORT_PATH.relative_to(broad_imports.ROOT)} is out of date")
        print("Regenerate it with:")
        print("  python3 tools/report_fortran_broad_imports.py --write")
        return 1
    print(f"Fortran broad-import audit check passed ({len(rows)} broad imports).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
