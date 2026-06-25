#!/usr/bin/env python3
"""Strict release gate for the non-LVLSET GPU-port objective.

This is intentionally stricter than the normal collaboration-readiness checks.
It answers whether the repository currently proves that every non-LVLSET
validation surface has CPU/GPU correctness evidence and an acceptable GPU
speedup claim.  Until the paired benchmark campaign is complete, this command is
expected to fail and report the remaining gaps.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from gpu_validation_plan import CASE_VARIANTS, CORRECTNESS_ONLY_VARIANTS


MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
PASSING_CORRECTNESS_STATUSES = {"passed", "accepted", "verified"}


def load_cases(path: Path) -> dict[str, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        case["id"]: case
        for case in data.get("cases", [])
        if isinstance(case, dict) and isinstance(case.get("id"), str)
    }


def speedup_claim_ok(claim: object) -> bool:
    return (
        isinstance(claim, dict)
        and isinstance(claim.get("speedup"), (int, float))
        and claim["speedup"] > 1.0
    )


def pair_result_ok(result: object) -> bool:
    return (
        isinstance(result, dict)
        and result.get("outcome") == "gpu_faster"
        and isinstance(result.get("speedup"), (int, float))
        and result["speedup"] > 1.0
    )


def pair_result_recorded(result: object) -> bool:
    return (
        isinstance(result, dict)
        and result.get("outcome") in {"gpu_faster", "gpu_not_faster"}
        and isinstance(result.get("speedup"), (int, float))
        and result["speedup"] > 0.0
    )


def correctness_ok(correctness: object) -> bool:
    if not isinstance(correctness, dict):
        return False
    if correctness.get("status") not in PASSING_CORRECTNESS_STATUSES:
        return False
    checks = correctness.get("checks")
    return (
        isinstance(checks, list)
        and bool(checks)
        and all(isinstance(check, str) and check.strip() for check in checks)
    )


def required_items_ok(correctness: object, required_items: set[str]) -> bool:
    if not required_items:
        return True
    if not isinstance(correctness, dict):
        return False
    evidence_items = correctness.get("evidence_items")
    if not isinstance(evidence_items, list):
        return False
    present = {item for item in evidence_items if isinstance(item, str)}
    return required_items <= present


def variant_gap(
    case_id: str,
    variant_id: str,
    evidence_case: dict,
    required_items: set[str],
) -> list[str]:
    gaps: list[str] = []
    claims = evidence_case.get("variant_speedup_claims")
    results = evidence_case.get("variant_pair_results")
    correctness = evidence_case.get("variant_correctness")
    claim = claims.get(variant_id) if isinstance(claims, dict) else None
    result = results.get(variant_id) if isinstance(results, dict) else None
    correct = correctness.get(variant_id) if isinstance(correctness, dict) else None

    label = f"`{case_id}` variant `{variant_id}`"
    if (case_id, variant_id) in CORRECTNESS_ONLY_VARIANTS:
        if not pair_result_recorded(result):
            gaps.append(f"{label} has no paired CPU/GPU result")
        if not correctness_ok(correct):
            gaps.append(f"{label} has no passing correctness evidence")
        elif not required_items_ok(correct, required_items):
            evidence_items = correct.get("evidence_items", [])
            present = {item for item in evidence_items if isinstance(item, str)}
            missing = sorted(required_items - present)
            gaps.append(
                f"{label} is missing required evidence items: "
                + ", ".join(f"`{item}`" for item in missing)
            )
        return gaps

    if not speedup_claim_ok(claim):
        gaps.append(f"{label} has no GPU-faster speedup claim")
    if not pair_result_ok(result):
        gaps.append(f"{label} has no gpu_faster paired result")
    if not correctness_ok(correct):
        gaps.append(f"{label} has no passing correctness evidence")
    elif not required_items_ok(correct, required_items):
        evidence_items = correct.get("evidence_items", [])
        present = {item for item in evidence_items if isinstance(item, str)}
        missing = sorted(required_items - present)
        gaps.append(
            f"{label} is missing required evidence items: "
            + ", ".join(f"`{item}`" for item in missing)
        )
    return gaps


def case_gap(case_id: str, manifest_case: dict, evidence_case: dict | None) -> list[str]:
    if evidence_case is None:
        return [f"`{case_id}` is missing from the evidence ledger"]

    status = manifest_case.get("status")
    if status == "excluded" or case_id == "lvlset":
        return []

    if status == "host_boundary_compatibility":
        if evidence_case.get("evidence_status") == "host_boundary_pending":
            return [f"`{case_id}` host-boundary compatibility evidence is still pending"]
        if not isinstance(evidence_case.get("source"), str) or not evidence_case["source"]:
            return [f"`{case_id}` host-boundary evidence has no source"]
        return []

    variants = CASE_VARIANTS.get(case_id, [])
    required_items = {
        item
        for item in manifest_case.get("required_evidence", [])
        if isinstance(item, str) and item
    }
    if len(variants) > 1:
        gaps: list[str] = []
        for variant in variants:
            gaps.extend(variant_gap(case_id, variant.id, evidence_case, required_items))
        return gaps

    gaps = []
    if evidence_case.get("evidence_status") != "paired_speedup_claimed":
        gaps.append(f"`{case_id}` evidence status is not paired_speedup_claimed")
    if not speedup_claim_ok(evidence_case.get("speedup_claim")):
        gaps.append(f"`{case_id}` has no GPU-faster speedup claim")
    if not pair_result_ok(evidence_case.get("pair_result")):
        gaps.append(f"`{case_id}` has no gpu_faster paired result")
    if not correctness_ok(evidence_case.get("correctness")):
        gaps.append(f"`{case_id}` has no passing correctness evidence")
    elif not required_items_ok(evidence_case.get("correctness"), required_items):
        correctness = evidence_case.get("correctness")
        evidence_items = (
            correctness.get("evidence_items", [])
            if isinstance(correctness, dict)
            else []
        )
        present = {item for item in evidence_items if isinstance(item, str)}
        missing = sorted(required_items - present)
        gaps.append(
            f"`{case_id}` is missing required evidence items: "
            + ", ".join(f"`{item}`" for item in missing)
        )
    return gaps


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fail unless the non-LVLSET GPU release objective is fully proven "
            "by paired CPU/GPU correctness and speedup evidence."
        )
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE_PATH)
    parser.add_argument("--json", action="store_true", help="emit machine-readable gaps")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest = load_cases(args.manifest)
    evidence = load_cases(args.evidence)

    gaps: list[str] = []
    for case_id, manifest_case in manifest.items():
        gaps.extend(case_gap(case_id, manifest_case, evidence.get(case_id)))

    if args.json:
        print(json.dumps({"release_objective_met": not gaps, "gaps": gaps}, indent=2))
    elif gaps:
        print("Non-LVLSET GPU release objective is not yet proven.")
        print()
        for gap in gaps:
            print(f"- {gap}")
    else:
        print("Non-LVLSET GPU release objective is fully proven.")

    return 0 if not gaps else 1


if __name__ == "__main__":
    sys.exit(main())
