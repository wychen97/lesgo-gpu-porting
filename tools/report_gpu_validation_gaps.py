#!/usr/bin/env python3
"""Report remaining non-LVLSET GPU validation and speedup evidence gaps."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
RUNBOOK_PATH = Path("docs/gpu_validation_runbook.md")

BATCH_RE = re.compile(r"^## Batch (?P<number>\d+): (?P<title>.+)$")
ROW_ID_RE = re.compile(r"^-\s+`(?P<id>[a-z0-9_]+)`\s*$")


def load_cases(path: Path) -> dict[str, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        case["id"]: case
        for case in data.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }


def parse_runbook_batches() -> dict[str, str]:
    batches: dict[str, str] = {}
    current: str | None = None
    in_rows = False
    for line in RUNBOOK_PATH.read_text(encoding="utf-8").splitlines():
        batch_match = BATCH_RE.match(line)
        if batch_match:
            current = f"Batch {batch_match.group('number')}: {batch_match.group('title')}"
            in_rows = False
            continue
        if line.strip() == "## Excluded Scope":
            current = "Excluded Scope"
            in_rows = False
            continue
        if line.startswith("## "):
            current = None
            in_rows = False
            continue
        if line.strip() == "Rows:":
            in_rows = True
            continue
        if line.startswith("Required action:"):
            in_rows = False
            continue
        if not in_rows or current is None:
            continue
        row_match = ROW_ID_RE.match(line)
        if row_match:
            batches[row_match.group("id")] = current
    return batches


def runtime_average(runtime: object) -> str:
    if not isinstance(runtime, dict):
        return ""
    value = runtime.get("cumulative_average_s_per_step")
    if isinstance(value, (int, float)):
        return f"{value:.6g}"
    return ""


def case_state(case: dict) -> str:
    if case.get("speedup_claim"):
        return "speedup-claimed"
    if isinstance(case.get("cpu_runtime"), dict) and isinstance(
        case.get("gpu_runtime"), dict
    ):
        return "paired-no-speedup-claim"
    if isinstance(case.get("gpu_runtime"), dict):
        return "gpu-runtime-only"
    return str(case.get("evidence_status"))


def next_action(manifest: dict, evidence: dict) -> str:
    if evidence.get("speedup_claim"):
        return "verify release wording"
    status = evidence.get("evidence_status")
    manifest_status = manifest.get("status")
    if status == "gpu_runtime_only":
        return "add paired CPU log and current-source proof"
    if status == "external_record_needs_copy":
        return "copy prior record or rerun"
    if status == "host_boundary_pending":
        return "run compatibility/overhead check"
    if status == "host_boundary_verified":
        return "verify compatibility record"
    if status == "excluded":
        return "excluded"
    if manifest_status == "needs_new_case_matrix":
        return "run CPU/GPU matrix"
    if manifest_status == "needs_new_case":
        return "create and run paired CPU/GPU case"
    if manifest_status == "needs_current_cpu_gpu_reference":
        return "rerun or copy current-source CPU/GPU logs"
    return "add evidence"


def print_markdown(open_only: bool) -> None:
    manifest = load_cases(MANIFEST_PATH)
    evidence = load_cases(EVIDENCE_PATH)
    batches = parse_runbook_batches()
    ids = list(manifest)

    rows = []
    for case_id in ids:
        manifest_case = manifest[case_id]
        evidence_case = evidence[case_id]
        state = case_state(evidence_case)
        if open_only and state == "speedup-claimed":
            continue
        rows.append((case_id, manifest_case, evidence_case, state))

    by_batch: dict[str, Counter[str]] = defaultdict(Counter)
    for case_id, _, _, state in rows:
        by_batch[batches.get(case_id, "Unbatched")][state] += 1

    total = len(rows)
    speedup_count = sum(1 for _, _, evidence_case, _ in rows if evidence_case.get("speedup_claim"))
    print("# GPU Validation Gap Report")
    print()
    print(f"Cases reported: {total}")
    print(f"Speedup claims in reported set: {speedup_count}")
    print()
    print("| Batch | Cases | Main states |")
    print("| --- | ---: | --- |")
    for batch in sorted(by_batch):
        states = ", ".join(
            f"`{state}`={count}" for state, count in sorted(by_batch[batch].items())
        )
        print(f"| {batch} | {sum(by_batch[batch].values())} | {states} |")

    print()
    print(
        "| ID | Batch | Manifest status | Evidence state | CPU avg s/step | "
        "GPU avg s/step | Speedup | Next action |"
    )
    print("| --- | --- | --- | --- | ---: | ---: | ---: | --- |")
    for case_id, manifest_case, evidence_case, state in rows:
        claim = evidence_case.get("speedup_claim")
        speedup = ""
        if isinstance(claim, dict) and isinstance(claim.get("speedup"), (int, float)):
            speedup = f"{claim['speedup']:.6g}"
        print(
            f"| `{case_id}` | {batches.get(case_id, 'Unbatched')} | "
            f"`{manifest_case.get('status')}` | `{state}` | "
            f"{runtime_average(evidence_case.get('cpu_runtime'))} | "
            f"{runtime_average(evidence_case.get('gpu_runtime'))} | "
            f"{speedup} | {next_action(manifest_case, evidence_case)} |"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report remaining GPU validation and speedup evidence gaps."
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="include rows that already have speedup claims",
    )
    args = parser.parse_args()
    print_markdown(open_only=not args.all)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
