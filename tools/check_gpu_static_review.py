#!/usr/bin/env python3
"""Verify static GPU review buckets stay synchronized with the audit doc."""

from __future__ import annotations

import re
import sys
import json
from collections import Counter
from pathlib import Path

from report_gpu_static_inventory import BUCKET_VALIDATION_ROWS, collect, review_bucket


DOC_PATH = Path("docs/gpu_port_coverage_audit.md")
MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
BUCKET_HEADER = "| Review bucket | Candidates | Meaning |"
BUCKET_ROW_RE = re.compile(
    r"^\|\s*`(?P<bucket>[a-z0-9-]+)`\s*\|\s*(?P<count>[0-9]+)\s*\|"
)
VALIDATION_HEADER = "| Review bucket | Validation rows |"
VALIDATION_ROW_RE = re.compile(
    r"^\|\s*`(?P<bucket>[a-z0-9-]+)`\s*\|\s*(?P<rows>.*?)\s*\|"
)
VALIDATION_ID_RE = re.compile(r"`(?P<id>[a-z0-9_]+)`")

def documented_bucket_counts() -> dict[str, int]:
    rows: dict[str, int] = {}
    in_table = False
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip() == BUCKET_HEADER:
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("| ---"):
            continue
        if not line.startswith("|"):
            break
        match = BUCKET_ROW_RE.match(line)
        if match:
            rows[match.group("bucket")] = int(match.group("count"))
    return rows


def documented_bucket_validation_rows() -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    in_table = False
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip() == VALIDATION_HEADER:
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("| ---"):
            continue
        if not line.startswith("|"):
            break
        match = VALIDATION_ROW_RE.match(line)
        if match:
            rows[match.group("bucket")] = [
                id_match.group("id")
                for id_match in VALIDATION_ID_RE.finditer(match.group("rows"))
            ]
    return rows


def actual_bucket_counts() -> Counter[str]:
    counts: Counter[str] = Counter()
    for item in collect():
        if item.classification == "unmarked-runtime-candidate":
            counts[review_bucket(item)] += 1
    return counts


def manifest_ids() -> set[str]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    return {
        case["id"]
        for case in manifest.get("cases", [])
        if isinstance(case, dict) and isinstance(case.get("id"), str)
    }


def main() -> int:
    actual = actual_bucket_counts()
    documented = documented_bucket_counts()
    documented_validation_rows = documented_bucket_validation_rows()
    known_validation_ids = manifest_ids()
    ok = True

    manual_count = actual.get("needs-manual-review", 0)
    if manual_count:
        print(
            "Static GPU review has unclassified candidates in "
            "`needs-manual-review`:"
        )
        for item in collect():
            if (
                item.classification == "unmarked-runtime-candidate"
                and review_bucket(item) == "needs-manual-review"
            ):
                print(f"  {item.path}:{item.start_line} {item.name}")
        ok = False

    missing = sorted(set(actual) - set(documented))
    stale = sorted(set(documented) - set(actual))
    if missing:
        print(f"{DOC_PATH} is missing static review buckets:")
        for bucket in missing:
            print(f"  {bucket}")
        ok = False
    if stale:
        print(f"{DOC_PATH} has stale static review buckets:")
        for bucket in stale:
            print(f"  {bucket}")
        ok = False

    for bucket in sorted(set(actual) & set(documented)):
        if documented[bucket] != actual[bucket]:
            print(
                f"{DOC_PATH} bucket `{bucket}` has count "
                f"{documented[bucket]}, expected {actual[bucket]}"
            )
            ok = False

    missing_validation_buckets = sorted(set(actual) - set(documented_validation_rows))
    stale_validation_buckets = sorted(set(documented_validation_rows) - set(actual))
    if missing_validation_buckets:
        print(f"{DOC_PATH} is missing static-review validation-row mappings:")
        for bucket in missing_validation_buckets:
            print(f"  {bucket}")
        ok = False
    if stale_validation_buckets:
        print(f"{DOC_PATH} has stale static-review validation-row mappings:")
        for bucket in stale_validation_buckets:
            print(f"  {bucket}")
        ok = False

    expected_mapping_buckets = set(BUCKET_VALIDATION_ROWS)
    missing_policy = sorted(set(actual) - expected_mapping_buckets)
    stale_policy = sorted(expected_mapping_buckets - set(actual))
    if missing_policy:
        print("Static GPU review policy is missing buckets:")
        for bucket in missing_policy:
            print(f"  {bucket}")
        ok = False
    if stale_policy:
        print("Static GPU review policy has stale buckets:")
        for bucket in stale_policy:
            print(f"  {bucket}")
        ok = False

    for bucket in sorted(set(actual) & set(documented_validation_rows)):
        expected_rows = BUCKET_VALIDATION_ROWS.get(bucket)
        documented_rows = documented_validation_rows[bucket]
        if expected_rows is None:
            continue
        if documented_rows != expected_rows:
            print(
                f"{DOC_PATH} validation mapping for `{bucket}` is "
                f"{documented_rows}, expected {expected_rows}"
            )
            ok = False
        for row_id in documented_rows:
            if row_id not in known_validation_ids:
                print(
                    f"{DOC_PATH} validation mapping for `{bucket}` references "
                    f"unknown validation row `{row_id}`"
                )
                ok = False

    documented_total = sum(documented.values())
    actual_total = sum(actual.values())
    if documented_total != actual_total:
        print(
            f"{DOC_PATH} documents {documented_total} candidates, "
            f"expected {actual_total}"
        )
        ok = False

    if not ok:
        return 1

    print(
        "GPU static review check passed "
        f"({actual_total} candidates, {len(actual)} buckets)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
