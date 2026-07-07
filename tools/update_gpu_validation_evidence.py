#!/usr/bin/env python3
"""Update one row in the GPU validation evidence ledger."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE_PATH = ROOT / "docs/gpu_validation_evidence.json"
PASSING_CORRECTNESS_STATUSES = {"passed", "accepted", "verified"}


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError("expected true/false")


def load_evidence(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def case_by_id(evidence: dict, case_id: str) -> dict:
    for case in evidence.get("cases", []):
        if isinstance(case, dict) and case.get("id") == case_id:
            return case
    raise KeyError(case_id)


def runtime_from_args(args: argparse.Namespace) -> dict:
    average = args.average_s_per_step
    if average is None and args.cumulative_seconds is not None:
        average = args.cumulative_seconds / args.nsteps
    return {
        "grid": args.grid,
        "nsteps": args.nsteps,
        "cumulative_seconds": args.cumulative_seconds,
        "cumulative_average_s_per_step": average,
        "final_printed_s_per_step": args.final_s_per_step,
    }


def has_passing_correctness(correctness: object) -> bool:
    if not isinstance(correctness, dict):
        return False
    if correctness.get("status") not in PASSING_CORRECTNESS_STATUSES:
        return False
    checks = correctness.get("checks")
    if not isinstance(checks, list) or not checks:
        return False
    return all(isinstance(check, str) and check.strip() for check in checks)


def require_passing_correctness(correctness: object, label: str) -> None:
    if not has_passing_correctness(correctness):
        raise ValueError(
            f"{label} speedup claim requires passing correctness evidence "
            "with a non-empty checks list"
        )


def correctness_from_args(args: argparse.Namespace) -> dict | None:
    evidence_items = [
        item.strip()
        for item in getattr(args, "evidence_item", [])
        if isinstance(item, str) and item.strip()
    ]
    has_any = (
        args.correctness_status is not None
        or args.correctness_source is not None
        or bool(args.correctness_check)
        or bool(evidence_items)
    )
    if not has_any:
        return None
    checks = [check.strip() for check in args.correctness_check if check.strip()]
    if not checks:
        raise ValueError("correctness evidence requires at least one --correctness-check")
    status = args.correctness_status or "passed"
    if status not in PASSING_CORRECTNESS_STATUSES:
        raise ValueError(
            "correctness status must be one of: "
            + ", ".join(sorted(PASSING_CORRECTNESS_STATUSES))
        )
    correctness = {
        "status": status,
        "checks": checks,
        "source": args.correctness_source or args.source,
    }
    if evidence_items:
        correctness["evidence_items"] = list(dict.fromkeys(evidence_items))
    return correctness


def set_correctness(case: dict, correctness: dict | None, variant: str | None) -> None:
    if correctness is None:
        return
    if variant:
        variants = case.setdefault("variant_correctness", {})
        if not isinstance(variants, dict):
            raise ValueError("variant_correctness exists but is not an object")
        variants[variant] = correctness
    else:
        case["correctness"] = correctness


def add_speedup_claim(case: dict) -> None:
    cpu = case.get("cpu_runtime")
    gpu = case.get("gpu_runtime")
    if not isinstance(cpu, dict) or not isinstance(gpu, dict):
        raise ValueError("speedup claim requires both cpu_runtime and gpu_runtime")
    cpu_avg = cpu.get("cumulative_average_s_per_step")
    gpu_avg = gpu.get("cumulative_average_s_per_step")
    if not isinstance(cpu_avg, (int, float)) or cpu_avg <= 0:
        raise ValueError("speedup claim requires positive CPU cumulative average")
    if not isinstance(gpu_avg, (int, float)) or gpu_avg <= 0:
        raise ValueError("speedup claim requires positive GPU cumulative average")
    require_passing_correctness(case.get("correctness"), "paired")
    speedup = cpu_avg / gpu_avg
    if speedup <= 1.0:
        raise ValueError(
            f"GPU is not faster on cumulative average: speedup={speedup:.6g}"
        )
    case["speedup_claim"] = {
        "speedup": speedup,
        "basis": "cumulative_average_s_per_step",
    }
    case["pair_result"] = {
        "speedup": speedup,
        "outcome": "gpu_faster",
        "basis": "cumulative_average_s_per_step",
    }
    case["evidence_status"] = "paired_speedup_claimed"


def add_not_faster_pair_result(case: dict) -> None:
    cpu = case.get("cpu_runtime")
    gpu = case.get("gpu_runtime")
    if not isinstance(cpu, dict) or not isinstance(gpu, dict):
        raise ValueError("paired result requires both cpu_runtime and gpu_runtime")
    cpu_avg = cpu.get("cumulative_average_s_per_step")
    gpu_avg = gpu.get("cumulative_average_s_per_step")
    if not isinstance(cpu_avg, (int, float)) or cpu_avg <= 0:
        raise ValueError("paired result requires positive CPU cumulative average")
    if not isinstance(gpu_avg, (int, float)) or gpu_avg <= 0:
        raise ValueError("paired result requires positive GPU cumulative average")
    speedup = cpu_avg / gpu_avg
    if speedup > 1.0:
        raise ValueError(
            f"GPU is faster on cumulative average; use speedup claim: {speedup:.6g}"
        )
    case["speedup_claim"] = None
    case["pair_result"] = {
        "speedup": speedup,
        "outcome": "gpu_not_faster",
        "basis": "cumulative_average_s_per_step",
    }
    case["evidence_status"] = "paired_not_faster"


def variant_entry(case: dict, variant: str) -> dict:
    variants = case.setdefault("variant_runtimes", {})
    if not isinstance(variants, dict):
        raise ValueError("variant_runtimes exists but is not an object")
    entry = variants.setdefault(variant, {})
    if not isinstance(entry, dict):
        raise ValueError(f"variant_runtimes.{variant} exists but is not an object")
    return entry


def add_variant_speedup_claim(case: dict, variant: str) -> None:
    entry = variant_entry(case, variant)
    cpu = entry.get("cpu_runtime")
    gpu = entry.get("gpu_runtime")
    if not isinstance(cpu, dict) or not isinstance(gpu, dict):
        raise ValueError("variant speedup claim requires both CPU and GPU runtime")
    cpu_avg = cpu.get("cumulative_average_s_per_step")
    gpu_avg = gpu.get("cumulative_average_s_per_step")
    if not isinstance(cpu_avg, (int, float)) or cpu_avg <= 0:
        raise ValueError("variant speedup claim requires positive CPU average")
    if not isinstance(gpu_avg, (int, float)) or gpu_avg <= 0:
        raise ValueError("variant speedup claim requires positive GPU average")
    variant_correctness = case.get("variant_correctness")
    correctness = (
        variant_correctness.get(variant)
        if isinstance(variant_correctness, dict)
        else None
    )
    require_passing_correctness(correctness, f"variant `{variant}`")
    speedup = cpu_avg / gpu_avg
    if speedup <= 1.0:
        raise ValueError(
            f"GPU is not faster for variant `{variant}`: speedup={speedup:.6g}"
        )
    claims = case.setdefault("variant_speedup_claims", {})
    if not isinstance(claims, dict):
        raise ValueError("variant_speedup_claims exists but is not an object")
    claims[variant] = {
        "speedup": speedup,
        "basis": "cumulative_average_s_per_step",
    }
    results = case.setdefault("variant_pair_results", {})
    if not isinstance(results, dict):
        raise ValueError("variant_pair_results exists but is not an object")
    results[variant] = {
        "speedup": speedup,
        "outcome": "gpu_faster",
        "basis": "cumulative_average_s_per_step",
    }


def add_variant_not_faster_pair_result(case: dict, variant: str) -> None:
    entry = variant_entry(case, variant)
    cpu = entry.get("cpu_runtime")
    gpu = entry.get("gpu_runtime")
    if not isinstance(cpu, dict) or not isinstance(gpu, dict):
        raise ValueError("variant paired result requires both CPU and GPU runtime")
    cpu_avg = cpu.get("cumulative_average_s_per_step")
    gpu_avg = gpu.get("cumulative_average_s_per_step")
    if not isinstance(cpu_avg, (int, float)) or cpu_avg <= 0:
        raise ValueError("variant paired result requires positive CPU average")
    if not isinstance(gpu_avg, (int, float)) or gpu_avg <= 0:
        raise ValueError("variant paired result requires positive GPU average")
    speedup = cpu_avg / gpu_avg
    if speedup > 1.0:
        raise ValueError(
            f"GPU is faster for variant `{variant}`; use speedup claim: {speedup:.6g}"
        )
    claims = case.setdefault("variant_speedup_claims", {})
    if not isinstance(claims, dict):
        raise ValueError("variant_speedup_claims exists but is not an object")
    claims.pop(variant, None)
    results = case.setdefault("variant_pair_results", {})
    if not isinstance(results, dict):
        raise ValueError("variant_pair_results exists but is not an object")
    results[variant] = {
        "speedup": speedup,
        "outcome": "gpu_not_faster",
        "basis": "cumulative_average_s_per_step",
    }


def write_evidence(path: Path, evidence: dict) -> None:
    path.write_text(
        json.dumps(evidence, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Update CPU or GPU timing evidence for one validation row. "
            "Run tools/check_gpu_validation_evidence.py afterwards."
        )
    )
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE_PATH)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--backend", choices=["cpu", "gpu"], required=True)
    parser.add_argument(
        "--variant",
        help=(
            "optional variant key, for rows with multiple runtime variants "
            "such as structure_off and structure_on"
        ),
    )
    parser.add_argument("--grid", required=True)
    parser.add_argument("--nsteps", type=positive_int, required=True)
    parser.add_argument("--cumulative-seconds", type=positive_float)
    parser.add_argument("--average-s-per-step", type=positive_float)
    parser.add_argument("--final-s-per-step", type=positive_float)
    parser.add_argument("--source", required=True)
    parser.add_argument("--notes")
    parser.add_argument("--public-record", type=parse_bool)
    parser.add_argument(
        "--correctness-status",
        choices=sorted(PASSING_CORRECTNESS_STATUSES),
        help="correctness status to store before claiming speedup; defaults to passed",
    )
    parser.add_argument(
        "--correctness-source",
        help="source/report for correctness evidence; defaults to --source",
    )
    parser.add_argument(
        "--correctness-check",
        action="append",
        default=[],
        help=(
            "repeatable correctness check summary, for example divergence, "
            "kinetic energy, turbine power, or field-difference comparison"
        ),
    )
    parser.add_argument(
        "--evidence-item",
        action="append",
        default=[],
        help=(
            "repeatable machine-readable evidence item from "
            "docs/gpu_benchmark_manifest.json, for example turbine_power or "
            "scalar_field_comparison"
        ),
    )
    parser.add_argument(
        "--claim-faster",
        action="store_true",
        help=(
            "add a speedup claim if CPU and GPU cumulative averages are both "
            "present and the GPU average is lower"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the updated case without writing the evidence file",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.cumulative_seconds is None and args.average_s_per_step is None:
        parser.error("provide --cumulative-seconds or --average-s-per-step")

    evidence = load_evidence(args.evidence)
    try:
        case = case_by_id(evidence, args.case_id)
    except KeyError:
        print(f"{args.evidence} has no case id `{args.case_id}`")
        return 1

    runtime = runtime_from_args(args)
    if args.variant:
        variant_entry(case, args.variant)[f"{args.backend}_runtime"] = runtime
    else:
        case[f"{args.backend}_runtime"] = runtime
    case["source"] = args.source
    if args.notes is not None:
        case["notes"] = args.notes
    if args.public_record is not None:
        case["public_record"] = args.public_record

    try:
        set_correctness(case, correctness_from_args(args), args.variant)
    except ValueError as exc:
        print(f"cannot record correctness evidence for `{args.case_id}`: {exc}")
        return 1

    if args.claim_faster:
        try:
            if args.variant:
                add_variant_speedup_claim(case, args.variant)
            else:
                add_speedup_claim(case)
        except ValueError as exc:
            print(f"cannot add speedup claim for `{args.case_id}`: {exc}")
            return 1

    if args.dry_run:
        print(json.dumps(case, indent=2, sort_keys=False))
        return 0

    write_evidence(args.evidence, evidence)
    target = (
        f"variant `{args.variant}` {args.backend}_runtime"
        if args.variant
        else f"{args.backend}_runtime"
    )
    print(f"Updated {args.evidence} case `{args.case_id}` {target}.")
    if case.get("speedup_claim"):
        speedup = case["speedup_claim"]["speedup"]
        print(f"Recorded GPU speedup claim: {speedup:.6g}x.")
    if args.variant and case.get("variant_speedup_claims", {}).get(args.variant):
        speedup = case["variant_speedup_claims"][args.variant]["speedup"]
        print(f"Recorded GPU speedup claim for variant `{args.variant}`: {speedup:.6g}x.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
