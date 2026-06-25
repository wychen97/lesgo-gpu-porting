#!/usr/bin/env python3
"""Import one LESGO timing log into the GPU validation evidence ledger."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from parse_lesgo_timing import parse_log, read_text
from update_gpu_validation_evidence import (
    DEFAULT_EVIDENCE_PATH,
    PASSING_CORRECTNESS_STATUSES,
    add_speedup_claim,
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


def runtime_from_parsed(parsed: dict, grid: str, nsteps: int | None) -> dict:
    steps = nsteps if nsteps is not None else parsed.get("last_iteration")
    return {
        "grid": grid,
        "nsteps": steps,
        "cumulative_seconds": parsed.get("cumulative_wall_s"),
        "cumulative_average_s_per_step": parsed.get(
            "cumulative_average_s_per_step"
        ),
        "final_printed_s_per_step": parsed.get("last_iteration_wall_s"),
    }


def validate_runtime(case_id: str, runtime: dict) -> list[str]:
    errors: list[str] = []
    if runtime.get("nsteps") is None:
        errors.append(f"case `{case_id}` timing log has no iteration count")
    if runtime.get("cumulative_seconds") is None:
        errors.append(f"case `{case_id}` timing log has no cumulative wall time")
    if runtime.get("cumulative_average_s_per_step") is None:
        errors.append(f"case `{case_id}` timing log has no cumulative average")
    if runtime.get("final_printed_s_per_step") is None:
        errors.append(f"case `{case_id}` timing log has no final printed step time")
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Parse a LESGO log and update one CPU or GPU timing row in "
            "docs/gpu_validation_evidence.json."
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
    parser.add_argument("--log", required=True, help="LESGO log path, or '-'")
    parser.add_argument("--grid", required=True)
    parser.add_argument(
        "--nsteps",
        type=positive_int,
        help="steps used for cumulative average; defaults to parsed iteration",
    )
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
            "docs/gpu_benchmark_manifest.json"
        ),
    )
    parser.add_argument("--claim-faster", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="allow import when optional timing fields are missing",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    parsed = parse_log(read_text(args.log), args.nsteps)
    runtime = runtime_from_parsed(parsed, args.grid, args.nsteps)

    errors = validate_runtime(args.case_id, runtime)
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
        variant_entry(case, args.variant)[f"{args.backend}_runtime"] = runtime
    else:
        case[f"{args.backend}_runtime"] = runtime
    case["source"] = args.source
    case["parsed_log_summary"] = parsed
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
    print(
        f"Imported {args.log} into {args.evidence} case "
        f"`{args.case_id}` {target}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
