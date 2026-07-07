"""Shared access to the readiness wrapper's Python check manifest."""

from __future__ import annotations

import ast
from pathlib import Path

from repo_paths import ROOT

WRAPPER_RELATIVE_PATH = Path("tools/check_branch_readiness.py")
WRAPPER_PATH = ROOT / WRAPPER_RELATIVE_PATH


def wrapper_checks() -> list[tuple[str, str]]:
    """Return ``(label, script)`` pairs from ``PYTHON_CHECKS``."""
    tree = ast.parse(
        WRAPPER_PATH.read_text(encoding="utf-8"),
        filename=str(WRAPPER_RELATIVE_PATH),
    )
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        names = [target.id for target in node.targets if isinstance(target, ast.Name)]
        if "PYTHON_CHECKS" not in names:
            continue
        checks = ast.literal_eval(node.value)
        return [(str(label), str(script)) for label, script in checks]
    raise RuntimeError(f"{WRAPPER_RELATIVE_PATH} does not define PYTHON_CHECKS")


def wrapper_script_paths() -> list[Path]:
    """Return script paths from ``PYTHON_CHECKS`` in wrapper order."""
    return [Path(script) for _label, script in wrapper_checks()]
