#!/usr/bin/env python3
"""Check lightweight source hygiene rules for the shared LESGO branch.

The optimized GPU branch is maintained across several machines and terminals.
Keeping tracked source, docs, and scripts UTF-8/ASCII with LF line endings and
spaces-only indentation avoids mojibake, CRLF churn, and tab rendering drift in
HPC shell logs and code reviews.  This script intentionally checks only small,
human-edited text files that are part of the source tree.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SOURCE_SUFFIXES = {
    ".c",
    ".cmake",
    ".f",
    ".f90",
    ".h",
    ".html",
    ".md",
    ".py",
    ".sh",
}

SOURCE_NAMES = {"CMakeLists.txt"}

MAX_ISSUES = 200


def tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        check=True,
        stdout=subprocess.PIPE,
    )
    names = result.stdout.decode("utf-8").split("\0")
    return [Path(name) for name in names if name]


def is_checked_path(path: Path) -> bool:
    if path.name in SOURCE_NAMES:
        return True
    return path.suffix.lower() in SOURCE_SUFFIXES


def split_body_newline(line: str) -> tuple[str, str]:
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\r"):
        return line[:-1], "\r"
    if line.endswith("\n"):
        return line[:-1], "\n"
    return line, ""


def check_file(
    path: Path,
    fix_trailing_whitespace: bool,
    fix_tabs: bool,
    fix_line_endings: bool,
) -> tuple[list[str], int, int, int]:
    issues: list[str] = []
    data = path.read_bytes()

    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        return [f"{path}: invalid UTF-8 at byte {exc.start}"], 0, 0, 0

    trailing_fixed = 0
    tabs_fixed = 0
    line_endings_fixed = 0
    output_lines: list[str] = []

    for lineno, line in enumerate(text.splitlines(keepends=True), start=1):
        body, newline = split_body_newline(line)
        if newline in {"\r\n", "\r"}:
            line_endings_fixed += 1
            if fix_line_endings:
                newline = "\n"
            elif newline == "\r\n":
                issues.append(f"{path}:{lineno}: CRLF line ending")
            else:
                issues.append(f"{path}:{lineno}: carriage-return line ending")

        if "\t" in body:
            tabs_fixed += body.count("\t")
            if fix_tabs:
                body = body.replace("\t", "    ")
            else:
                issues.append(f"{path}:{lineno}: tab character")

        stripped = body.rstrip(" \t")
        if stripped != body:
            trailing_fixed += 1
            if fix_trailing_whitespace:
                body = stripped
            else:
                issues.append(f"{path}:{lineno}: trailing whitespace")

        for char in body:
            if ord(char) > 0x7F:
                issues.append(
                    f"{path}:{lineno}: non-ASCII character U+{ord(char):04X}"
                )
                break

        output_lines.append(body + newline)

    if (
        (trailing_fixed and fix_trailing_whitespace)
        or (tabs_fixed and fix_tabs)
        or (line_endings_fixed and fix_line_endings)
    ):
        path.write_text("".join(output_lines), encoding="utf-8", newline="")

    return issues, trailing_fixed, tabs_fixed, line_endings_fixed


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check tracked source/docs/scripts for portable text hygiene."
    )
    parser.add_argument(
        "--list-files",
        action="store_true",
        help="print the tracked files included in the hygiene scan",
    )
    parser.add_argument(
        "--fix-trailing-whitespace",
        action="store_true",
        help="rewrite scanned files to remove trailing spaces and tabs",
    )
    parser.add_argument(
        "--fix-tabs",
        action="store_true",
        help="rewrite scanned files to replace tab characters with four spaces",
    )
    parser.add_argument(
        "--fix-line-endings",
        action="store_true",
        help="rewrite scanned files to normalize CRLF/carriage returns to LF",
    )
    args = parser.parse_args()

    paths = [path for path in tracked_files() if is_checked_path(path)]

    if args.list_files:
        for path in paths:
            print(path)
        return 0

    issues: list[str] = []
    trailing_fixed = 0
    tabs_fixed = 0
    line_endings_fixed = 0
    for path in paths:
        (
            file_issues,
            file_trailing_fixed,
            file_tabs_fixed,
            file_line_endings_fixed,
        ) = check_file(
            path,
            args.fix_trailing_whitespace,
            args.fix_tabs,
            args.fix_line_endings,
        )
        issues.extend(file_issues)
        trailing_fixed += file_trailing_fixed
        tabs_fixed += file_tabs_fixed
        line_endings_fixed += file_line_endings_fixed
        if len(issues) >= MAX_ISSUES:
            break

    if issues:
        print("Source hygiene check failed:")
        for issue in issues[:MAX_ISSUES]:
            print(f"  {issue}")
        if len(issues) >= MAX_ISSUES:
            print(f"  ... stopped after {MAX_ISSUES} issues")
        return 1

    if args.fix_trailing_whitespace:
        print(f"Removed trailing whitespace from {trailing_fixed} lines.")
    if args.fix_tabs:
        print(f"Replaced {tabs_fixed} tab characters with spaces.")
    if args.fix_line_endings:
        print(f"Normalized {line_endings_fixed} line endings to LF.")
    print(f"Source hygiene check passed ({len(paths)} tracked files scanned).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
