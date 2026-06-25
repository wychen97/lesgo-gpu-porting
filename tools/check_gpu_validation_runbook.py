#!/usr/bin/env python3
"""Verify the GPU validation runbook covers the manifest and evidence rows."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


RUNBOOK_PATH = Path("docs/gpu_validation_runbook.md")
MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
ID_RE = re.compile(r"`(?P<id>[a-z][a-z0-9_]+)`")

REQUIRED_PHRASES = [
    "paired CPU/GPU evidence",
    "same source tree",
    "cumulative average",
    "final printed step",
    "docs/gpu_validation_evidence.json",
]


def load_case_ids(path: Path) -> set[str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        case["id"]
        for case in data.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }


def main() -> int:
    text = RUNBOOK_PATH.read_text(encoding="utf-8")
    ids_in_runbook = set(ID_RE.findall(text))
    manifest_ids = load_case_ids(MANIFEST_PATH)
    evidence_ids = load_case_ids(EVIDENCE_PATH)
    ok = True

    if manifest_ids != evidence_ids:
        print(f"{MANIFEST_PATH} and {EVIDENCE_PATH} case ids differ")
        ok = False

    missing = sorted(manifest_ids - ids_in_runbook)
    if missing:
        print(f"{RUNBOOK_PATH} is missing validation ids:")
        for row_id in missing:
            print(f"  {row_id}")
        ok = False

    for phrase in REQUIRED_PHRASES:
        if phrase not in text:
            print(f"{RUNBOOK_PATH} is missing required phrase: {phrase}")
            ok = False

    batch_count = len(re.findall(r"^## Batch ", text, flags=re.MULTILINE))
    if batch_count < 5:
        print(f"{RUNBOOK_PATH} must contain at least 5 benchmark batches")
        ok = False

    if not ok:
        return 1

    print(
        "GPU validation runbook check passed "
        f"({len(manifest_ids)} ids, {batch_count} batches)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
