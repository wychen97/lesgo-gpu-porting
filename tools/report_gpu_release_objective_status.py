#!/usr/bin/env python3
"""Generate the non-LVLSET GPU release-objective status report."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from report_gpu_static_inventory import (
    BUCKET_VALIDATION_ROWS,
    REVIEW_BUCKET_DESCRIPTIONS,
    ROOT,
    collect,
    review_bucket,
)
from report_lesgo_conf_key_validation import group_rows
from check_lesgo_conf_coverage_docs import source_keys
from require_gpu_release_objective import case_gap


OUTPUT_PATH = ROOT / "docs/gpu_release_objective_status.md"
MANIFEST_PATH = ROOT / "docs/gpu_benchmark_manifest.json"
EVIDENCE_PATH = ROOT / "docs/gpu_validation_evidence.json"


def load_cases(path: Path) -> dict[str, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        case["id"]: case
        for case in data.get("cases", [])
        if isinstance(case, dict) and isinstance(case.get("id"), str)
    }


def evidence_state(
    manifest: dict[str, dict],
    evidence: dict[str, dict],
    row_id: str,
) -> str:
    row = evidence.get(row_id)
    if not isinstance(row, dict):
        return "missing"
    if row_id == "lvlset":
        return "excluded"
    manifest_case = manifest.get(row_id)
    if isinstance(manifest_case, dict) and not case_gap(row_id, manifest_case, row):
        return "release-proven"
    if row.get("evidence_status") == "host_boundary_pending":
        return "host-boundary-pending"
    if row.get("speedup_claim"):
        return "speedup-claimed"
    return str(row.get("evidence_status", "unknown"))


def row_status_summary(
    manifest: dict[str, dict],
    evidence: dict[str, dict],
    row_ids: list[str],
) -> str:
    return ", ".join(
        f"`{row_id}`={evidence_state(manifest, evidence, row_id)}"
        for row_id in row_ids
    )


def row_release_proven(
    manifest: dict[str, dict],
    evidence: dict[str, dict],
    row_id: str,
) -> bool:
    manifest_case = manifest.get(row_id)
    if not isinstance(manifest_case, dict):
        return False
    return not case_gap(row_id, manifest_case, evidence.get(row_id))


def markdown() -> str:
    manifest = load_cases(MANIFEST_PATH)
    evidence = load_cases(EVIDENCE_PATH)
    subprograms = collect()
    by_class = Counter(item.classification for item in subprograms)
    keys_by_group = source_keys()
    rows_by_group = group_rows()

    candidates = [
        item
        for item in subprograms
        if item.classification == "unmarked-runtime-candidate"
    ]
    by_bucket: dict[str, list] = defaultdict(list)
    for item in candidates:
        by_bucket[review_bucket(item)].append(item)

    release_gaps: list[str] = []
    for case_id, manifest_case in manifest.items():
        release_gaps.extend(case_gap(case_id, manifest_case, evidence.get(case_id)))
    evidence_counts = Counter(
        evidence_state(manifest, evidence, case_id)
        for case_id, case in evidence.items()
        if case.get("id") != "lvlset"
    )

    lines = [
        "# GPU Release Objective Status",
        "",
        "Scope: non-LVLSET LESGO GPU port.  This report connects the static",
        "Fortran function inventory, `lesgo.conf` validation surface, and paired",
        "CPU/GPU evidence ledger.  It is generated and should not be edited by",
        "hand.",
        "",
        f"Release objective met: `{'yes' if not release_gaps else 'no'}`",
        f"Strict release-gate gaps: `{len(release_gaps)}`",
        "",
        "Run the strict gate directly with:",
        "",
        "```bash",
        "python3 tools/require_gpu_release_objective.py",
        "```",
        "",
        "## Static Source Inventory",
        "",
        "| Classification | Subprograms | Meaning |",
        "| --- | ---: | --- |",
        (
            f"| `gpu-marked` | {by_class['gpu-marked']} | Contains OpenACC, CUDA, "
            "GPU-aware MPI, or GPU preprocessor markers. |"
        ),
        (
            f"| `gpu-file-unmarked` | {by_class['gpu-file-unmarked']} | Helper "
            "inside a GPU source file but without a local marker. |"
        ),
        (
            f"| `host-boundary` | {by_class['host-boundary']} | Setup, I/O, "
            "configuration, MPI definitions, or other expected host boundary. |"
        ),
        (
            f"| `host-or-diagnostic` | {by_class['host-or-diagnostic']} | "
            "Diagnostics, initialization, reporting, restart, or low-frequency host work. |"
        ),
        (
            f"| `unmarked-runtime-candidate` | {by_class['unmarked-runtime-candidate']} | "
            "Needs runtime validation or profiling before broad GPU speed claims. |"
        ),
        "",
        "## Validation Evidence Summary",
        "",
        "| Evidence state | Rows |",
        "| --- | ---: |",
    ]
    for state in sorted(evidence_counts):
        lines.append(f"| `{state}` | {evidence_counts[state]} |")

    lines.extend(
        [
            "",
            "`release-proven` rows have paired CPU/GPU correctness evidence and a",
            "GPU-faster result.  Rows with only a single GPU runtime remain execution",
            "evidence only.",
            "",
            "## lesgo.conf Parser-Key Coverage",
            "",
            "| Parser group | Keys | Release-proven keys | Validation row states |",
            "| --- | ---: | ---: | --- |",
        ]
    )
    for group in [
        "DOMAIN",
        "MODEL",
        "CORIOLIS",
        "FLOW_COND",
        "OUTPUT",
        "TURBINES",
        "SCALARS",
        "HIT",
    ]:
        keys = keys_by_group.get(group, [])
        rows = rows_by_group.get(group, [])
        proven = (
            len(keys)
            if rows and all(row_release_proven(manifest, evidence, row_id) for row_id in rows)
            else 0
        )
        lines.append(
            f"| `{group}` | {len(keys)} | {proven} | "
            f"{row_status_summary(manifest, evidence, rows)} |"
        )

    total_keys = sum(len(keys) for keys in keys_by_group.values())
    proven_keys = 0
    for group, keys in keys_by_group.items():
        rows = rows_by_group.get(group, [])
        if rows and all(row_release_proven(manifest, evidence, row_id) for row_id in rows):
            proven_keys += len(keys)

    lines.extend(
        [
            "",
            f"Total parsed non-LVLSET `lesgo.conf` keys: `{total_keys}`",
            f"Release-proven parser keys: `{proven_keys}`",
            "",
            "## Candidate Buckets Linked To Evidence",
            "",
            "| Review bucket | Candidates | Validation row states | Meaning |",
            "| --- | ---: | --- | --- |",
        ]
    )
    for bucket in sorted(by_bucket):
        rows = BUCKET_VALIDATION_ROWS.get(bucket, [])
        lines.append(
            f"| `{bucket}` | {len(by_bucket[bucket])} | "
            f"{row_status_summary(manifest, evidence, rows)} | "
            f"{REVIEW_BUCKET_DESCRIPTIONS[bucket]} |"
        )

    lines.extend(
        [
            "",
            "## Open Validation Rows",
            "",
            "| Validation row | Current evidence state | Required next evidence |",
            "| --- | --- | --- |",
        ]
    )
    for case_id, manifest_case in manifest.items():
        if case_id == "lvlset":
            continue
        row = evidence.get(case_id, {})
        state = evidence_state(manifest, evidence, case_id)
        if not case_gap(case_id, manifest_case, row):
            continue
        required = ", ".join(f"`{item}`" for item in manifest_case.get("required_evidence", []))
        lines.append(f"| `{case_id}` | `{state}` | {required} |")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The source tree has broad GPU coverage, but the release objective is not",
            "  proven until the open validation rows have paired CPU/GPU correctness",
            "  and timing evidence.",
            "- Host-boundary rows are not GPU hot-path speed claims; they need",
            "  compatibility and overhead evidence.",
            "- LVLSET remains excluded from this audit.",
            "",
            "Regenerate this file with:",
            "",
            "```bash",
            "python3 tools/report_gpu_release_objective_status.py --write",
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate the GPU release-objective status report."
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help=f"write the generated report to {OUTPUT_PATH}",
    )
    args = parser.parse_args()
    text = markdown()
    if args.write:
        OUTPUT_PATH.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
