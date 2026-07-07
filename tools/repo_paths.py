#!/usr/bin/env python3
"""Shared repository path helpers for top-level maintenance tools."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def repo_path(*parts: str) -> Path:
    """Return a path anchored at the repository root."""
    return ROOT.joinpath(*parts)
