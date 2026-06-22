"""Shared script inventory helpers for readiness checks."""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEST_CASE_ROOT = Path("test-cases")
TEST_CASE_SCRIPT_SUFFIXES = {".pbs", ".sbatch", ".sh"}


def git_ls_files(*patterns: str) -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", *patterns],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def tracked_root_shell_scripts() -> list[Path]:
    return [
        path
        for path in git_ls_files("*.sh")
        if path.parent == Path(".")
    ]


def tracked_test_case_files() -> set[Path]:
    return set(git_ls_files(str(TEST_CASE_ROOT)))


def owning_test_case_dir(path: Path) -> Path | None:
    parts = path.parts
    if len(parts) < 3 or parts[0] != TEST_CASE_ROOT.name:
        return None
    return Path(parts[0]) / parts[1]


def tracked_test_case_script_dirs() -> set[Path]:
    return {
        case_dir
        for path in tracked_test_case_files()
        if path.suffix in TEST_CASE_SCRIPT_SUFFIXES
        for case_dir in [owning_test_case_dir(path)]
        if case_dir is not None
    }
