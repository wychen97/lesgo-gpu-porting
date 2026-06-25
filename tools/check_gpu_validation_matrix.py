#!/usr/bin/env python3
"""Verify the non-LVLSET GPU validation matrix has required rows."""

from __future__ import annotations

import re
import sys
from pathlib import Path


DOC_PATH = Path("docs/gpu_validation_matrix.md")

REQUIRED_IDS = {
    "les_core_channel",
    "adm_disk",
    "atm_line",
    "large_windfarm",
    "sgs_disabled",
    "sgs_models_1_5",
    "dyn_tn",
    "iwm_wall_model",
    "scalar_passive",
    "scalar_active",
    "cps_velocity",
    "cps_scalar",
    "hit_inflow",
    "shifted_inflow",
    "sponge_coriolis",
    "adm_dynamic_controls",
    "diagnostics_output",
    "cgns_output",
    "lvlset",
}

ALLOWED_STATUSES = {
    "validated-gpu-runtime",
    "recorded-correct-and-faster",
    "recorded-paired-not-faster",
    "implemented-needs-benchmark",
    "host-boundary",
    "excluded",
}

STATUS_RE = re.compile(r"`(?P<status>[a-z-]+)`")
ROW_RE = re.compile(r"^\|\s*`(?P<id>[a-z0-9_]+)`\s*\|(?P<body>.*)\|$")


def table_rows() -> dict[str, str]:
    rows: dict[str, str] = {}
    in_matrix = False
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip() == "## Validation Matrix":
            in_matrix = True
            continue
        if in_matrix and line.startswith("## "):
            break
        if not in_matrix:
            continue
        match = ROW_RE.match(line)
        if not match:
            continue
        row_id = match.group("id")
        rows[row_id] = match.group("body")
    return rows


def main() -> int:
    rows = table_rows()
    ok = True

    missing = sorted(REQUIRED_IDS - set(rows))
    stale = sorted(set(rows) - REQUIRED_IDS)
    if missing:
        print(f"{DOC_PATH} is missing validation rows:")
        for row_id in missing:
            print(f"  {row_id}")
        ok = False
    if stale:
        print(f"{DOC_PATH} has stale validation rows:")
        for row_id in stale:
            print(f"  {row_id}")
        ok = False

    for row_id, body in sorted(rows.items()):
        statuses = STATUS_RE.findall(body)
        row_statuses = [status for status in statuses if status in ALLOWED_STATUSES]
        if len(row_statuses) != 1:
            print(
                f"{DOC_PATH} row `{row_id}` must contain exactly one allowed "
                "status token."
            )
            ok = False
        unknown = sorted(set(statuses) - ALLOWED_STATUSES)
        if unknown:
            print(f"{DOC_PATH} row `{row_id}` has unknown status tokens:")
            for status in unknown:
                print(f"  {status}")
            ok = False

    if not ok:
        return 1

    print(
        "GPU validation matrix check passed "
        f"({len(rows)} rows, {len(ALLOWED_STATUSES)} statuses)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
