#!/usr/bin/env python3
"""Verify lesgo.conf parser groups map to GPU validation rows."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from check_lesgo_conf_coverage_docs import source_keys


MAP_PATH = Path("docs/lesgo_conf_gpu_validation_map.json")
MATRIX_PATH = Path("docs/gpu_validation_matrix.md")
MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
ROW_RE = re.compile(r"^\|\s*`(?P<id>[a-z0-9_]+)`\s*\|")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def matrix_ids() -> set[str]:
    ids: set[str] = set()
    in_matrix = False
    for line in MATRIX_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip() == "## Validation Matrix":
            in_matrix = True
            continue
        if in_matrix and line.startswith("## "):
            break
        if not in_matrix:
            continue
        match = ROW_RE.match(line)
        if match:
            ids.add(match.group("id"))
    return ids


def json_case_ids(path: Path) -> set[str]:
    data = load_json(path)
    return {
        case["id"]
        for case in data.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }


def main() -> int:
    validation_map = load_json(MAP_PATH)
    ok = True
    if validation_map.get("schema_version") != 1:
        print(f"{MAP_PATH} must have schema_version 1")
        ok = False

    groups = validation_map.get("parser_groups")
    if not isinstance(groups, list):
        print(f"{MAP_PATH} must contain a parser_groups list")
        return 1

    expected_groups = set(source_keys())
    documented_groups = {
        group.get("group")
        for group in groups
        if isinstance(group, dict) and isinstance(group.get("group"), str)
    }

    missing_groups = sorted(expected_groups - documented_groups)
    stale_groups = sorted(documented_groups - expected_groups)
    if missing_groups:
        print(f"{MAP_PATH} is missing parser groups:")
        for group in missing_groups:
            print(f"  {group}")
        ok = False
    if stale_groups:
        print(f"{MAP_PATH} has stale parser groups:")
        for group in stale_groups:
            print(f"  {group}")
        ok = False

    matrix = matrix_ids()
    manifest = json_case_ids(MANIFEST_PATH)
    evidence = json_case_ids(EVIDENCE_PATH)
    known_ids = matrix & manifest & evidence
    if matrix != manifest or matrix != evidence:
        print("Validation row ids differ across matrix, manifest, and evidence")
        ok = False

    referenced: set[str] = set()
    for group in groups:
        if not isinstance(group, dict):
            print(f"{MAP_PATH} has a non-object parser group entry")
            ok = False
            continue
        group_name = group.get("group", "<missing group>")
        rows = group.get("validation_rows")
        if not isinstance(rows, list) or not rows:
            print(f"{MAP_PATH} group `{group_name}` needs validation_rows")
            ok = False
            continue
        for row_id in rows:
            if row_id not in known_ids:
                print(
                    f"{MAP_PATH} group `{group_name}` references unknown "
                    f"validation row `{row_id}`"
                )
                ok = False
            else:
                referenced.add(row_id)

    cross_cutting = validation_map.get("cross_cutting_rows", [])
    if not isinstance(cross_cutting, list):
        print(f"{MAP_PATH} cross_cutting_rows must be a list")
        ok = False
        cross_cutting = []
    for entry in cross_cutting:
        if not isinstance(entry, dict):
            print(f"{MAP_PATH} has a non-object cross_cutting_rows entry")
            ok = False
            continue
        row_id = entry.get("id")
        if row_id not in known_ids:
            print(f"{MAP_PATH} cross_cutting row references unknown id `{row_id}`")
            ok = False
        else:
            referenced.add(row_id)

    unreferenced = sorted(known_ids - referenced)
    if unreferenced:
        print(f"{MAP_PATH} does not account for validation rows:")
        for row_id in unreferenced:
            print(f"  {row_id}")
        ok = False

    if not ok:
        return 1

    total_keys = sum(len(keys) for keys in source_keys().values())
    print(
        "lesgo.conf validation map check passed "
        f"({len(expected_groups)} parser groups, {total_keys} keys, "
        f"{len(known_ids)} validation rows)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
