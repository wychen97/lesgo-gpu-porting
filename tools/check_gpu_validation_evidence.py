#!/usr/bin/env python3
"""Verify the GPU validation evidence ledger is consistent and conservative."""

from __future__ import annotations

import json
import re
import sys
import argparse
from pathlib import Path

from require_gpu_release_objective import case_gap


DEFAULT_EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")
MATRIX_PATH = Path("docs/gpu_validation_matrix.md")
MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
ROW_RE = re.compile(r"^\|\s*`(?P<id>[a-z0-9_]+)`\s*\|(?P<body>.*)\|$")
STATUS_RE = re.compile(r"`(?P<status>[a-z-]+)`")

REQUIRED_FIELDS = {
    "id",
    "evidence_status",
    "public_record",
    "source",
    "gpu_runtime",
    "cpu_runtime",
    "speedup_claim",
    "notes",
}

ALLOWED_EVIDENCE_STATUSES = {
    "gpu_runtime_only",
    "needs_benchmark",
    "external_record_needs_copy",
    "paired_speedup_claimed",
    "paired_not_faster",
    "host_boundary_pending",
    "host_boundary_verified",
    "excluded",
}
PAIR_OUTCOMES = {"gpu_faster", "gpu_not_faster"}
PASSING_CORRECTNESS_STATUSES = {"passed", "accepted", "verified"}

MATRIX_TO_EVIDENCE = {
    "validated-gpu-runtime": {"gpu_runtime_only"},
    "recorded-correct-and-faster": {
        "external_record_needs_copy",
        "paired_speedup_claimed",
    },
    "recorded-paired-not-faster": {"paired_not_faster"},
    "implemented-needs-benchmark": {"needs_benchmark", "paired_not_faster"},
    "host-boundary": {"host_boundary_pending", "host_boundary_verified"},
    "excluded": {"excluded"},
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def matrix_statuses() -> dict[str, str]:
    rows: dict[str, str] = {}
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
        if not match:
            continue
        statuses = STATUS_RE.findall(match.group("body"))
        if len(statuses) == 1:
            rows[match.group("id")] = statuses[0]
    return rows


def manifest_ids() -> set[str]:
    manifest = load_json(MANIFEST_PATH)
    return {
        case["id"]
        for case in manifest.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }


def manifest_cases() -> dict[str, dict]:
    manifest = load_json(MANIFEST_PATH)
    return {
        case["id"]: case
        for case in manifest.get("cases", [])
        if isinstance(case, dict) and isinstance(case.get("id"), str)
    }


def manifest_required_evidence() -> dict[str, set[str]]:
    manifest = load_json(MANIFEST_PATH)
    rows: dict[str, set[str]] = {}
    for case in manifest.get("cases", []):
        if not isinstance(case, dict) or not isinstance(case.get("id"), str):
            continue
        required = case.get("required_evidence")
        rows[case["id"]] = {
            item
            for item in required
            if isinstance(item, str) and item
        } if isinstance(required, list) else set()
    return rows


def is_positive_number(value: object) -> bool:
    return isinstance(value, (int, float)) and value > 0


def validate_runtime(case_id: str, label: str, runtime: object) -> list[str]:
    errors: list[str] = []
    if runtime is None:
        return errors
    if not isinstance(runtime, dict):
        return [f"case `{case_id}` {label} must be null or an object"]
    nsteps = runtime.get("nsteps")
    if nsteps is not None and not is_positive_number(nsteps):
        errors.append(f"case `{case_id}` {label}.nsteps must be positive or null")
    for field in [
        "cumulative_seconds",
        "cumulative_average_s_per_step",
        "final_printed_s_per_step",
    ]:
        value = runtime.get(field)
        if value is not None and not is_positive_number(value):
            errors.append(f"case `{case_id}` {label}.{field} must be positive or null")
    return errors


def validate_correctness(
    case_id: str,
    label: str,
    correctness: object,
    *,
    required: bool,
) -> list[str]:
    if correctness is None:
        if required:
            return [f"case `{case_id}` {label} correctness evidence is required"]
        return []
    if not isinstance(correctness, dict):
        return [f"case `{case_id}` {label} correctness must be an object"]

    errors: list[str] = []
    status = correctness.get("status")
    if status not in PASSING_CORRECTNESS_STATUSES:
        errors.append(
            f"case `{case_id}` {label} correctness status must be one of "
            f"{sorted(PASSING_CORRECTNESS_STATUSES)}"
        )
    checks = correctness.get("checks")
    if not isinstance(checks, list) or not checks:
        errors.append(
            f"case `{case_id}` {label} correctness checks must be a non-empty list"
        )
    elif not all(isinstance(check, str) and check.strip() for check in checks):
        errors.append(
            f"case `{case_id}` {label} correctness checks must be non-empty strings"
        )
    source = correctness.get("source")
    if source is not None and (not isinstance(source, str) or not source.strip()):
        errors.append(f"case `{case_id}` {label} correctness source must be a string")
    evidence_items = correctness.get("evidence_items")
    if evidence_items is not None:
        if not isinstance(evidence_items, list) or not evidence_items:
            errors.append(
                f"case `{case_id}` {label} evidence_items must be a non-empty list"
            )
        elif not all(isinstance(item, str) and item.strip() for item in evidence_items):
            errors.append(
                f"case `{case_id}` {label} evidence_items must be non-empty strings"
            )
    return errors


def validate_required_evidence_items(
    case_id: str,
    label: str,
    correctness: object,
    required_items: set[str],
) -> list[str]:
    if not required_items:
        return []
    if not isinstance(correctness, dict):
        return [
            f"case `{case_id}` {label} cannot prove required evidence without "
            "correctness object"
        ]
    evidence_items = correctness.get("evidence_items")
    if not isinstance(evidence_items, list):
        return [
            f"case `{case_id}` {label} correctness has no machine-readable "
            "evidence_items list"
        ]
    present = {item for item in evidence_items if isinstance(item, str)}
    missing = sorted(required_items - present)
    if missing:
        return [
            f"case `{case_id}` {label} missing required evidence_items: "
            + ", ".join(f"`{item}`" for item in missing)
        ]
    return []


def validate_speedup(case: dict, required_items: set[str]) -> list[str]:
    case_id = case.get("id", "<missing id>")
    claim = case.get("speedup_claim")
    if claim is None:
        variant_claims = case.get("variant_speedup_claims")
        has_variant_claims = isinstance(variant_claims, dict) and bool(variant_claims)
        if case.get("evidence_status") == "paired_speedup_claimed" and not has_variant_claims:
            return [
                f"case `{case_id}` evidence_status is paired_speedup_claimed "
                "but speedup_claim and variant_speedup_claims are null"
            ]
        return []
    errors: list[str] = []
    if case.get("evidence_status") != "paired_speedup_claimed":
        errors.append(
            f"case `{case_id}` has speedup_claim but evidence_status is "
            f"`{case.get('evidence_status')}`, expected `paired_speedup_claimed`"
        )
    if not isinstance(claim, dict):
        return [f"case `{case_id}` speedup_claim must be null or an object"]
    if not is_positive_number(claim.get("speedup")):
        errors.append(f"case `{case_id}` speedup_claim.speedup must be positive")
    for label in ["cpu_runtime", "gpu_runtime"]:
        runtime = case.get(label)
        if not isinstance(runtime, dict):
            errors.append(
                f"case `{case_id}` has speedup_claim but no {label} object"
            )
            continue
        if not is_positive_number(runtime.get("cumulative_average_s_per_step")):
            errors.append(
                f"case `{case_id}` has speedup_claim but no {label} "
                "cumulative_average_s_per_step"
            )
    errors.extend(
        validate_correctness(
            case_id,
            "paired",
            case.get("correctness"),
            required=True,
        )
    )
    errors.extend(
        validate_required_evidence_items(
            case_id,
            "paired",
            case.get("correctness"),
            required_items,
        )
    )
    return errors


def validate_pair_result(case: dict) -> list[str]:
    case_id = case.get("id", "<missing id>")
    pair = case.get("pair_result")
    if pair is None:
        if case.get("evidence_status") == "paired_not_faster":
            return [
                f"case `{case_id}` evidence_status is paired_not_faster "
                "but pair_result is null"
            ]
        return []
    if not isinstance(pair, dict):
        return [f"case `{case_id}` pair_result must be null or an object"]
    errors: list[str] = []
    outcome = pair.get("outcome")
    if outcome not in PAIR_OUTCOMES:
        errors.append(f"case `{case_id}` pair_result.outcome is invalid: {outcome}")
    if not is_positive_number(pair.get("speedup")):
        errors.append(f"case `{case_id}` pair_result.speedup must be positive")
    if outcome == "gpu_faster" and case.get("speedup_claim") is None:
        errors.append(
            f"case `{case_id}` pair_result is gpu_faster but speedup_claim is null"
        )
    if outcome == "gpu_not_faster" and case.get("speedup_claim") is not None:
        errors.append(
            f"case `{case_id}` pair_result is gpu_not_faster but speedup_claim exists"
        )
    return errors


def validate_variant_runtimes(case: dict, required_items: set[str]) -> list[str]:
    case_id = case.get("id", "<missing id>")
    variants = case.get("variant_runtimes")
    if variants is None:
        return []
    if not isinstance(variants, dict):
        return [f"case `{case_id}` variant_runtimes must be an object"]

    errors: list[str] = []
    for variant, entry in variants.items():
        if not isinstance(variant, str) or not variant:
            errors.append(f"case `{case_id}` has an invalid variant key")
            continue
        if not isinstance(entry, dict):
            errors.append(
                f"case `{case_id}` variant `{variant}` entry must be an object"
            )
            continue
        if "cpu_runtime" not in entry and "gpu_runtime" not in entry:
            errors.append(
                f"case `{case_id}` variant `{variant}` has no CPU/GPU runtime"
            )
        for label in ["cpu_runtime", "gpu_runtime"]:
            errors.extend(
                validate_runtime(
                    case_id,
                    f"variant `{variant}` {label}",
                    entry.get(label),
                )
            )

    claims = case.get("variant_speedup_claims")
    if claims is not None and not isinstance(claims, dict):
        errors.append(f"case `{case_id}` variant_speedup_claims must be an object")
        claims = None
    results = case.get("variant_pair_results")
    if results is not None and not isinstance(results, dict):
        errors.append(f"case `{case_id}` variant_pair_results must be an object")
        results = None

    if results is not None:
        for variant, result in results.items():
            entry = variants.get(variant)
            if not isinstance(entry, dict):
                errors.append(
                    f"case `{case_id}` has variant pair result for unknown "
                    f"variant `{variant}`"
                )
                continue
            if not isinstance(result, dict):
                errors.append(
                    f"case `{case_id}` variant `{variant}` pair result must be an object"
                )
                continue
            outcome = result.get("outcome")
            if outcome not in PAIR_OUTCOMES:
                errors.append(
                    f"case `{case_id}` variant `{variant}` pair outcome is invalid: "
                    f"{outcome}"
                )
            if not is_positive_number(result.get("speedup")):
                errors.append(
                    f"case `{case_id}` variant `{variant}` pair speedup must be positive"
                )
            has_claim = isinstance(claims, dict) and variant in claims
            if outcome == "gpu_faster" and not has_claim:
                errors.append(
                    f"case `{case_id}` variant `{variant}` is gpu_faster "
                    "but has no speedup claim"
                )
            if outcome == "gpu_not_faster" and has_claim:
                errors.append(
                    f"case `{case_id}` variant `{variant}` is gpu_not_faster "
                    "but has a speedup claim"
                )
            for label in ["cpu_runtime", "gpu_runtime"]:
                runtime = entry.get(label)
                if not isinstance(runtime, dict):
                    errors.append(
                        f"case `{case_id}` variant `{variant}` pair result "
                        f"has no {label} object"
                    )
                    continue
                if not is_positive_number(runtime.get("cumulative_average_s_per_step")):
                    errors.append(
                        f"case `{case_id}` variant `{variant}` pair result "
                        f"has no {label} cumulative average"
                    )
    if claims is None:
        return errors
    if not isinstance(claims, dict):
        errors.append(f"case `{case_id}` variant_speedup_claims must be an object")
        return errors
    variant_correctness = case.get("variant_correctness")
    if variant_correctness is not None and not isinstance(variant_correctness, dict):
        errors.append(f"case `{case_id}` variant_correctness must be an object")
        variant_correctness = None
    for variant, claim in claims.items():
        entry = variants.get(variant)
        if not isinstance(entry, dict):
            errors.append(
                f"case `{case_id}` has variant speedup for unknown variant `{variant}`"
            )
            continue
        if not isinstance(claim, dict):
            errors.append(
                f"case `{case_id}` variant `{variant}` speedup claim must be an object"
            )
            continue
        if not is_positive_number(claim.get("speedup")):
            errors.append(
                f"case `{case_id}` variant `{variant}` speedup must be positive"
            )
        for label in ["cpu_runtime", "gpu_runtime"]:
            runtime = entry.get(label)
            if not isinstance(runtime, dict):
                errors.append(
                    f"case `{case_id}` variant `{variant}` speedup claim "
                    f"has no {label} object"
                )
                continue
            if not is_positive_number(runtime.get("cumulative_average_s_per_step")):
                errors.append(
                    f"case `{case_id}` variant `{variant}` speedup claim "
                    f"has no {label} cumulative average"
                )
        correctness = (
            variant_correctness.get(variant)
            if isinstance(variant_correctness, dict)
            else None
        )
        errors.extend(
            validate_correctness(
                case_id,
                f"variant `{variant}`",
                correctness,
                required=True,
            )
        )
        errors.extend(
            validate_required_evidence_items(
                case_id,
                f"variant `{variant}`",
                correctness,
                required_items,
            )
        )
    return errors


def matrix_status_accepts_evidence(
    case_id: str,
    matrix_status: str | None,
    evidence_status: object,
    evidence_case: dict,
    manifest_case: dict | None,
) -> bool:
    allowed_for_matrix = MATRIX_TO_EVIDENCE.get(matrix_status, set())
    if evidence_status in allowed_for_matrix:
        return True

    # Variant-backed rows keep the top-level evidence_status conservative
    # while the per-variant CPU/GPU pairs prove the matrix row.  Reuse the
    # strict release-gate logic so this checker cannot accidentally accept
    # an incomplete variant set.
    if (
        matrix_status == "recorded-correct-and-faster"
        and manifest_case is not None
        and not case_gap(case_id, manifest_case, evidence_case)
    ):
        return True

    return False


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify the GPU validation evidence ledger."
    )
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE_PATH)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    evidence_path = args.evidence
    evidence = load_json(evidence_path)
    ok = True
    if evidence.get("schema_version") != 1:
        print(f"{evidence_path} must have schema_version 1")
        ok = False

    cases = evidence.get("cases")
    if not isinstance(cases, list):
        print(f"{evidence_path} must contain a cases list")
        return 1

    matrix = matrix_statuses()
    manifest_by_id = manifest_cases()
    manifest = set(manifest_by_id)
    required_evidence = manifest_required_evidence()
    case_ids = [case.get("id") for case in cases if isinstance(case, dict)]
    evidence_ids = set(case_ids)

    if len(case_ids) != len(evidence_ids):
        print(f"{evidence_path} has duplicate case ids")
        ok = False

    for source_name, source_ids in [
        (str(MATRIX_PATH), set(matrix)),
        (str(MANIFEST_PATH), manifest),
    ]:
        missing = sorted(source_ids - evidence_ids)
        stale = sorted(evidence_ids - source_ids)
        if missing:
            print(f"{evidence_path} is missing ids from {source_name}:")
            for row_id in missing:
                print(f"  {row_id}")
            ok = False
        if stale:
            print(f"{evidence_path} has ids not present in {source_name}:")
            for row_id in stale:
                print(f"  {row_id}")
            ok = False

    for case in cases:
        if not isinstance(case, dict):
            print(f"{evidence_path} has a non-object case entry")
            ok = False
            continue
        case_id = case.get("id", "<missing id>")
        missing_fields = sorted(REQUIRED_FIELDS - set(case))
        if missing_fields:
            print(f"{evidence_path} case `{case_id}` is missing fields:")
            for field in missing_fields:
                print(f"  {field}")
            ok = False

        evidence_status = case.get("evidence_status")
        if evidence_status not in ALLOWED_EVIDENCE_STATUSES:
            print(
                f"{evidence_path} case `{case_id}` has invalid "
                f"evidence_status: {evidence_status}"
            )
            ok = False

        matrix_status = matrix.get(case_id)
        if not matrix_status_accepts_evidence(
            case_id,
            matrix_status,
            evidence_status,
            case,
            manifest_by_id.get(case_id),
        ):
            print(
                f"{evidence_path} case `{case_id}` status `{evidence_status}` "
                f"does not match matrix status `{matrix_status}`"
            )
            ok = False

        if not isinstance(case.get("public_record"), bool):
            print(f"{evidence_path} case `{case_id}` public_record must be boolean")
            ok = False

        if evidence_status == "paired_not_faster":
            if case.get("speedup_claim") is not None:
                print(
                    f"{evidence_path} case `{case_id}` is paired_not_faster "
                    "but has speedup_claim"
                )
                ok = False
            if not isinstance(case.get("cpu_runtime"), dict) or not isinstance(
                case.get("gpu_runtime"), dict
            ):
                print(
                    f"{evidence_path} case `{case_id}` is paired_not_faster "
                    "but does not have both CPU and GPU runtime objects"
                )
                ok = False

        for error in validate_runtime(case_id, "gpu_runtime", case.get("gpu_runtime")):
            print(f"{evidence_path} {error}")
            ok = False
        for error in validate_runtime(case_id, "cpu_runtime", case.get("cpu_runtime")):
            print(f"{evidence_path} {error}")
            ok = False
        case_required_evidence = required_evidence.get(case_id, set())
        for error in validate_speedup(case, case_required_evidence):
            print(f"{evidence_path} {error}")
            ok = False
        for error in validate_pair_result(case):
            print(f"{evidence_path} {error}")
            ok = False
        for error in validate_correctness(
            case_id,
            "paired",
            case.get("correctness"),
            required=False,
        ):
            print(f"{evidence_path} {error}")
            ok = False
        for error in validate_variant_runtimes(case, case_required_evidence):
            print(f"{evidence_path} {error}")
            ok = False

    if not ok:
        return 1

    speedup_claims = 0
    for case in cases:
        if not isinstance(case, dict):
            continue
        if case.get("speedup_claim"):
            speedup_claims += 1
        variant_claims = case.get("variant_speedup_claims")
        if isinstance(variant_claims, dict):
            speedup_claims += len(variant_claims)
    print(
        "GPU validation evidence check passed "
        f"({len(cases)} cases, {speedup_claims} speedup claims)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
