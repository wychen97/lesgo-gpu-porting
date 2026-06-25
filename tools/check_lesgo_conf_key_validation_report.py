#!/usr/bin/env python3
"""Verify generated key-level lesgo.conf validation coverage is current."""

from __future__ import annotations

import sys
from pathlib import Path

from report_lesgo_conf_key_validation import OUTPUT_PATH, markdown_table


def main() -> int:
    expected = markdown_table()
    if not OUTPUT_PATH.exists():
        print(f"{OUTPUT_PATH} is missing")
        print("Regenerate it with:")
        print("  python3 tools/report_lesgo_conf_key_validation.py --write")
        return 1

    actual = OUTPUT_PATH.read_text(encoding="utf-8")
    if actual != expected:
        print(f"{OUTPUT_PATH} is out of date")
        print("Regenerate it with:")
        print("  python3 tools/report_lesgo_conf_key_validation.py --write")
        return 1

    row_count = sum(1 for line in actual.splitlines() if line.startswith("| `"))
    print(
        "lesgo.conf key validation report check passed "
        f"({row_count} parser keys)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
