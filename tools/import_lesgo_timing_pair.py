#!/usr/bin/env python3
"""Import paired CPU/GPU LESGO timing logs into the evidence ledger."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from import_lesgo_timing_evidence import runtime_from_parsed, validate_runtime
from parse_lesgo_timing import parse_log, read_text
from update_gpu_validation_evidence import (
    DEFAULT_EVIDENCE_PATH,
    PASSING_CORRECTNESS_STATUSES,
    add_not_faster_pair_result,
    add_speedup_claim,
    add_variant_not_faster_pair_result,
    add_variant_speedup_claim,
    case_by_id,
    correctness_from_args,
    load_evidence,
    parse_bool,
    set_correctness,
    variant_entry,
    write_evidence,
)


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def diagnostic_value(parsed: dict, key: str) -> float | None:
    diagnostics = parsed.get("diagnostics")
    if not isinstance(diagnostics, dict):
        return None
    value = diagnostics.get(key)
    return value if isinstance(value, (int, float)) else None


def compare_float(
    label: str,
    cpu_value: float | None,
    gpu_value: float | None,
    *,
    rtol: float,
    atol: float,
) -> tuple[str | None, str | None]:
    if cpu_value is None or gpu_value is None:
        return None, f"missing {label} diagnostic in CPU or GPU log"
    if not math.isfinite(cpu_value) or not math.isfinite(gpu_value):
        return None, f"non-finite {label} diagnostic in CPU or GPU log"
    diff = abs(cpu_value - gpu_value)
    tolerance = max(atol, rtol * max(abs(cpu_value), abs(gpu_value)))
    if diff > tolerance:
        return (
            None,
            f"{label} mismatch: cpu={cpu_value:.8g}, gpu={gpu_value:.8g}, "
            f"abs_diff={diff:.3g}, tolerance={tolerance:.3g}",
        )
    return (
        f"{label} matched: cpu={cpu_value:.8g}, gpu={gpu_value:.8g}, "
        f"abs_diff={diff:.3g}, tolerance={tolerance:.3g}",
        None,
    )


def diagnostic_correctness_checks(args: argparse.Namespace, cpu: dict, gpu: dict) -> list[str]:
    checks: list[str] = []
    errors: list[str] = []

    cpu_status = cpu.get("mpi_exit_status")
    gpu_status = gpu.get("mpi_exit_status")
    if cpu_status != 0 or gpu_status != 0:
        errors.append(f"MPI exit status mismatch/failure: cpu={cpu_status}, gpu={gpu_status}")
    else:
        checks.append("MPI exit status matched: cpu=0, gpu=0")

    for label, key, rtol, atol in [
        (
            "divergence",
            "divergence",
            args.divergence_rtol,
            args.divergence_atol,
        ),
        (
            "kinetic energy",
            "kinetic_energy",
            args.kinetic_energy_rtol,
            args.kinetic_energy_atol,
        ),
    ]:
        check, error = compare_float(
            label,
            diagnostic_value(cpu, key),
            diagnostic_value(gpu, key),
            rtol=rtol,
            atol=atol,
        )
        if error:
            errors.append(error)
        else:
            assert check is not None
            checks.append(check)

    if errors:
        raise ValueError("; ".join(errors))
    return checks


def correctness_with_optional_diagnostics(
    args: argparse.Namespace,
    cpu: dict,
    gpu: dict,
) -> dict | None:
    manual_checks = [check.strip() for check in args.correctness_check if check.strip()]
    evidence_items = [
        item.strip()
        for item in args.evidence_item
        if isinstance(item, str) and item.strip()
    ]
    auto_items = ["cpu_log", "gpu_log"]
    if (
        isinstance(cpu.get("cumulative_average_s_per_step"), (int, float))
        and isinstance(gpu.get("cumulative_average_s_per_step"), (int, float))
    ):
        auto_items.append("cumulative_average")
    if (
        isinstance(cpu.get("last_iteration_wall_s"), (int, float))
        and isinstance(gpu.get("last_iteration_wall_s"), (int, float))
    ):
        auto_items.append("late_step_timing")
    has_manual = (
        args.correctness_status is not None
        or args.correctness_source is not None
        or bool(manual_checks)
        or bool(evidence_items)
    )
    if not args.compare_diagnostics:
        correctness = correctness_from_args(args)
        if correctness is not None:
            items = list(dict.fromkeys([*auto_items, *correctness.get("evidence_items", [])]))
            correctness["evidence_items"] = items
        return correctness

    status = args.correctness_status or "passed"
    if status not in PASSING_CORRECTNESS_STATUSES:
        raise ValueError(
            "correctness status must be one of: "
            + ", ".join(sorted(PASSING_CORRECTNESS_STATUSES))
        )
    diagnostic_checks = diagnostic_correctness_checks(args, cpu, gpu)
    checks = [*manual_checks, *diagnostic_checks]
    evidence_items = [*evidence_items, *auto_items, "divergence", "kinetic_energy"]
    if not checks and has_manual:
        raise ValueError("correctness evidence requires at least one check")
    return {
        "status": status,
        "checks": checks,
        "source": args.correctness_source or args.source,
        "evidence_items": list(dict.fromkeys(evidence_items)),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Parse matched CPU/GPU LESGO logs and record a paired result. "
            "The result is a speedup claim only when the GPU cumulative "
            "average is lower than the CPU cumulative average."
        )
    )
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE_PATH)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--variant")
    parser.add_argument("--cpu-log", required=True)
    parser.add_argument("--gpu-log", required=True)
    parser.add_argument("--grid", required=True)
    parser.add_argument("--nsteps", type=positive_int)
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
        "--compare-diagnostics",
        action="store_true",
        help=(
            "compare parsed CPU/GPU MPI exit status, divergence, and kinetic "
            "energy and store the result as structured correctness evidence"
        ),
    )
    parser.add_argument(
        "--divergence-rtol",
        type=nonnegative_float,
        default=1.0e-2,
        help="relative tolerance for parsed divergence comparison",
    )
    parser.add_argument(
        "--divergence-atol",
        type=nonnegative_float,
        default=1.0e-10,
        help="absolute tolerance for parsed divergence comparison",
    )
    parser.add_argument(
        "--kinetic-energy-rtol",
        type=nonnegative_float,
        default=1.0e-6,
        help="relative tolerance for parsed kinetic-energy comparison",
    )
    parser.add_argument(
        "--kinetic-energy-atol",
        type=nonnegative_float,
        default=1.0e-8,
        help="absolute tolerance for parsed kinetic-energy comparison",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="allow import when optional timing fields are missing",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    parsed_cpu = parse_log(read_text(args.cpu_log), args.nsteps)
    parsed_gpu = parse_log(read_text(args.gpu_log), args.nsteps)
    cpu_runtime = runtime_from_parsed(parsed_cpu, args.grid, args.nsteps)
    gpu_runtime = runtime_from_parsed(parsed_gpu, args.grid, args.nsteps)

    errors = []
    errors.extend(validate_runtime(args.case_id, cpu_runtime))
    errors.extend(validate_runtime(args.case_id, gpu_runtime))
    if errors and not args.allow_incomplete:
        for error in errors:
            print(error)
        return 1

    evidence = load_evidence(args.evidence)
    try:
        case = case_by_id(evidence, args.case_id)
    except KeyError:
        print(f"{args.evidence} has no case id `{args.case_id}`")
        return 1

    if args.variant:
        entry = variant_entry(case, args.variant)
        entry["cpu_runtime"] = cpu_runtime
        entry["gpu_runtime"] = gpu_runtime
    else:
        case["cpu_runtime"] = cpu_runtime
        case["gpu_runtime"] = gpu_runtime

    case["source"] = args.source
    case["parsed_pair_log_summary"] = {
        "cpu": parsed_cpu,
        "gpu": parsed_gpu,
    }
    if args.notes is not None:
        case["notes"] = args.notes
    if args.public_record is not None:
        case["public_record"] = args.public_record

    try:
        set_correctness(
            case,
            correctness_with_optional_diagnostics(args, parsed_cpu, parsed_gpu),
            args.variant,
        )
    except ValueError as exc:
        print(f"cannot record correctness evidence for `{args.case_id}`: {exc}")
        return 1

    cpu_avg = cpu_runtime.get("cumulative_average_s_per_step")
    gpu_avg = gpu_runtime.get("cumulative_average_s_per_step")
    if not isinstance(cpu_avg, (int, float)) or not isinstance(gpu_avg, (int, float)):
        print("cannot decide paired result without CPU and GPU cumulative averages")
        return 1

    try:
        if cpu_avg / gpu_avg > 1.0:
            if args.variant:
                add_variant_speedup_claim(case, args.variant)
            else:
                add_speedup_claim(case)
        else:
            if args.variant:
                add_variant_not_faster_pair_result(case, args.variant)
            else:
                add_not_faster_pair_result(case)
    except ValueError as exc:
        print(f"cannot record paired result for `{args.case_id}`: {exc}")
        return 1

    if args.dry_run:
        print(json.dumps(case, indent=2, sort_keys=False))
        return 0

    write_evidence(args.evidence, evidence)
    target = f" variant `{args.variant}`" if args.variant else ""
    speedup = cpu_avg / gpu_avg
    outcome = "gpu_faster" if speedup > 1.0 else "gpu_not_faster"
    print(
        f"Imported paired logs into {args.evidence} case `{args.case_id}`"
        f"{target}: {outcome}, speedup={speedup:.6g}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
