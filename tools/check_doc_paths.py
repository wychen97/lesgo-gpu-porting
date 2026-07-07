#!/usr/bin/env python3
"""Verify documented repository-local handoff paths still exist."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

from readiness_manifest import wrapper_script_paths
from repo_paths import ROOT


PATH_RE = re.compile(
    r"(?<![A-Za-z0-9_/${}-])"
    r"((?:\.\./)*(?:(?:docs|tools)/[A-Za-z0-9_./-]+|"
    r"(?:README_FINAL_OPTIMIZED_20260619|CONTRIBUTING|[A-Za-z0-9_]+_PORT)\.md)"
    r")"
)
TEXT_SUFFIXES = {
    ".c",
    ".cmake",
    ".f",
    ".f90",
    ".F",
    ".F90",
    ".h",
    ".md",
    ".py",
    ".sh",
}
TEXT_NAMES = {"CMakeLists.txt"}
STATIC_REQUIRED_ENTRYPOINTS = {
    Path("README_FINAL_OPTIMIZED_20260619.md"),
    Path("CONTRIBUTING.md"),
    Path("docs/README.md"),
    Path("docs/build_profiles.md"),
    Path("docs/cluster_scripts.md"),
    Path("docs/code_organization.md"),
    Path("docs/gpu_module_contracts.md"),
    Path("docs/source_file_inventory.md"),
    Path("docs/gpu_development_guidelines.md"),
    Path("docs/environment_switches.md"),
    Path("docs/refactor_backlog.md"),
    Path("tools/README.md"),
    Path("tools/check_branch_readiness.py"),
}


def is_checked_path(path: Path) -> bool:
    if path.name in TEXT_NAMES:
        return True
    return path.suffix in TEXT_SUFFIXES


def tracked_text_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    names = result.stdout.decode("utf-8").split("\0")
    return [ROOT / name for name in names if name and is_checked_path(Path(name))]


def readiness_wrapper_scripts() -> set[Path]:
    return set(wrapper_script_paths())


def required_entrypoints() -> set[Path]:
    return STATIC_REQUIRED_ENTRYPOINTS | readiness_wrapper_scripts()


def trim_path(raw: str) -> Path:
    return Path(raw.rstrip(".,:;)]}"))


def resolve_documented_path(source_path: Path, documented: Path) -> Path:
    if documented.parts and documented.parts[0] == "..":
        return (source_path.parent / documented).resolve()
    return (ROOT / documented).resolve()


def format_missing_path(target: Path) -> str:
    try:
        return str(target.relative_to(ROOT.resolve()))
    except ValueError:
        return str(target)


def main() -> int:
    missing: list[str] = []

    for path in sorted(required_entrypoints()):
        if not (ROOT / path).exists():
            missing.append(f"required entry point is missing: {path}")

    for source_path in tracked_text_files():
        text = source_path.read_text(encoding="utf-8")
        for raw in PATH_RE.findall(text):
            relative = trim_path(raw)
            target = resolve_documented_path(source_path, relative)
            if not target.exists():
                source = source_path.relative_to(ROOT)
                missing.append(
                    f"{source}: referenced path does not exist: "
                    f"{format_missing_path(target)}"
                )

    if missing:
        print("Documentation path check failed:")
        for item in missing:
            print(f"  {item}")
        return 1

    print("Documentation path check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
