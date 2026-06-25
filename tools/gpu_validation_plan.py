#!/usr/bin/env python3
"""Expand the GPU validation manifest into concrete paired run tasks."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
EVIDENCE_PATH = Path("docs/gpu_validation_evidence.json")

COMPLETE_EVIDENCE = {
    "paired_speedup_claimed",
    "paired_not_faster",
    "host_boundary_verified",
    "excluded",
}
PASSING_CORRECTNESS_STATUSES = {"passed", "accepted", "verified"}

PRIORITY_DESCRIPTIONS = {
    "p0-public-small": "paired CPU/GPU proof for the small public presentation cases",
    "p1-core-options": "core LES/SGS/scalar/HIT options with high coverage value",
    "p2-optional-coupling": "optional coupling, forcing, and larger scalar variants",
    "p3-large": "large wind-farm reference case",
    "p4-compatibility": "host-boundary output compatibility checks",
}


@dataclass(frozen=True)
class Variant:
    id: str
    note: str


CASE_VARIANTS: dict[str, list[Variant]] = {
    "les_core_channel": [Variant("baseline", "240^3 channel-flow baseline")],
    "adm_disk": [Variant("baseline", "240^3 actuator-disk baseline")],
    "atm_line": [
        Variant("structure_off", "ATM line model with structure disabled"),
        Variant("structure_on", "ATM line model with structure enabled"),
    ],
    "large_windfarm": [
        Variant("baseline", "3072x384x400, 60 turbines, matched output cadence")
    ],
    "sgs_disabled": [Variant("sgs_off", "sgs=.false. molecular branch")],
    "sgs_models_1_5": [
        Variant(f"sgs_model_{model}", f"sgs_model={model}") for model in range(1, 6)
    ],
    "dyn_tn": [
        Variant("sgs_model_4_dyn_tn", "USE_DYN_TN=ON with sgs_model=4"),
        Variant("sgs_model_5_dyn_tn", "USE_DYN_TN=ON with sgs_model=5"),
    ],
    "iwm_wall_model": [Variant("iwm_heavy", "IWM-heavy wall-model case")],
    "scalar_passive": [
        Variant("grid_128", "128^3 passive scalar"),
        Variant("grid_240", "240^3 passive scalar"),
    ],
    "scalar_active": [
        Variant("grid_128", "128^3 active scalar"),
        Variant("grid_240", "240^3 active scalar"),
    ],
    "cps_velocity": [Variant("two_color_velocity", "two-color CPS, scalar off")],
    "cps_scalar": [Variant("two_color_scalar", "two-color CPS with scalar coupling")],
    "hit_inflow": [
        Variant("hit_64_correctness", "64^3 HIT correctness check"),
        Variant("hit_128_timing", "128^3 HIT timing check"),
    ],
    "shifted_inflow": [Variant("shifted_inflow", "shifted inflow with matched input field")],
    "sponge_coriolis": [
        Variant("sponge_on", "sponge forcing enabled"),
        Variant("coriolis_on", "Coriolis forcing enabled"),
    ],
    "adm_dynamic_controls": [
        Variant("dynamic_yaw", "dynamic yaw control enabled"),
        Variant("dynamic_ct", "dynamic Ct control enabled"),
        Variant("rotation", "ADM rotation enabled"),
        Variant("adm_correction", "ADM correction enabled"),
    ],
    "diagnostics_output": [
        Variant("checkpoint_restart", "checkpoint/restart compatibility"),
        Variant("time_average", "time-average output compatibility"),
        Variant("plane_domain_output", "plane/domain output compatibility"),
    ],
    "cgns_output": [Variant("cgns_output", "CGNS output compatibility")],
    "lvlset": [],
}

CORRECTNESS_ONLY_VARIANTS = {
    ("hit_inflow", "hit_64_correctness"),
}


def task_priority(case_id: str, variant_id: str) -> str:
    if case_id in {"les_core_channel", "adm_disk", "atm_line"}:
        return "p0-public-small"
    if case_id in {"sgs_disabled", "sgs_models_1_5", "dyn_tn"}:
        return "p1-core-options"
    if case_id in {"scalar_passive", "scalar_active"} and variant_id == "grid_128":
        return "p1-core-options"
    if case_id == "hit_inflow" and variant_id == "hit_64_correctness":
        return "p1-core-options"
    if case_id == "large_windfarm":
        return "p3-large"
    if case_id in {"diagnostics_output", "cgns_output"}:
        return "p4-compatibility"
    return "p2-optional-coupling"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def evidence_by_id(evidence: dict) -> dict[str, dict]:
    return {
        case["id"]: case
        for case in evidence.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }


def case_goal(case: dict) -> str:
    status = case.get("status")
    if status == "host_boundary_compatibility":
        return "compatibility"
    if status == "needs_public_evidence_copy_or_rerun":
        return "copy-or-rerun"
    return "performance"


def run_needed(evidence_status: str | None) -> bool:
    return evidence_status not in COMPLETE_EVIDENCE


def task_status(evidence_status: str | None) -> str:
    if evidence_status == "paired_speedup_claimed":
        return "complete-speedup"
    if evidence_status == "paired_not_faster":
        return "complete-not-faster"
    if evidence_status == "external_record_needs_copy":
        return "copy-prior-record-or-rerun"
    if evidence_status == "host_boundary_pending":
        return "run-compatibility-check"
    if evidence_status == "host_boundary_verified":
        return "complete-compatibility"
    if evidence_status == "gpu_runtime_only":
        return "add-paired-cpu-gpu"
    if evidence_status == "needs_benchmark":
        return "run-paired-cpu-gpu"
    if evidence_status == "excluded":
        return "excluded"
        return "unknown"


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


def correctness_ok(correctness: object, required_items: set[str]) -> bool:
    if not isinstance(correctness, dict):
        return False
    if correctness.get("status") not in PASSING_CORRECTNESS_STATUSES:
        return False
    checks = correctness.get("checks")
    if not (
        isinstance(checks, list)
        and bool(checks)
        and all(isinstance(check, str) and check.strip() for check in checks)
    ):
        return False
    if required_items:
        evidence_items = correctness.get("evidence_items")
        if not isinstance(evidence_items, list):
            return False
        present = {item for item in evidence_items if isinstance(item, str)}
        if not required_items <= present:
            return False
    return True


def variant_task_status(
    evidence_row: dict,
    variant_id: str,
    required_items: set[str],
) -> str:
    claims = evidence_row.get("variant_speedup_claims")
    results = evidence_row.get("variant_pair_results")
    correctness = evidence_row.get("variant_correctness")
    claim = claims.get(variant_id) if isinstance(claims, dict) else None
    result = results.get(variant_id) if isinstance(results, dict) else None
    correct = correctness.get(variant_id) if isinstance(correctness, dict) else None
    if (
        (case_id := evidence_row.get("id"))
        and (case_id, variant_id) in CORRECTNESS_ONLY_VARIANTS
    ):
        if pair_result_recorded(result) and correctness_ok(correct, required_items):
            return "complete-correctness"
        return task_status(evidence_row.get("evidence_status"))
    if (
        speedup_claim_ok(claim)
        and pair_result_ok(result)
        and correctness_ok(correct, required_items)
    ):
        return "complete-speedup"
    return task_status(evidence_row.get("evidence_status"))


def task_run_needed(status: str) -> bool:
    return status not in {
        "complete-speedup",
        "complete-not-faster",
        "complete-correctness",
        "complete-compatibility",
        "excluded",
    }


def expand_tasks(manifest: dict, evidence: dict) -> list[dict]:
    evidence_rows = evidence_by_id(evidence)
    tasks: list[dict] = []
    for case in manifest.get("cases", []):
        if not isinstance(case, dict):
            continue
        case_id = case.get("id")
        if not isinstance(case_id, str):
            continue
        variants = CASE_VARIANTS.get(case_id, [])
        if not variants:
            continue
        uses_variant_runtime = len(variants) > 1
        evidence_row = evidence_rows.get(case_id, {})
        evidence_status = evidence_row.get("evidence_status")
        required_items = {
            item
            for item in case.get("required_evidence", [])
            if isinstance(item, str) and item
        }
        for variant in variants:
            status = (
                variant_task_status(evidence_row, variant.id, required_items)
                if uses_variant_runtime
                else task_status(evidence_status)
            )
            for run_kind, cmake_field in [("cpu", "cpu_cmake"), ("gpu", "gpu_cmake")]:
                cmake = case.get(cmake_field, [])
                if not cmake:
                    continue
                tasks.append(
                    {
                        "task_id": f"{case_id}.{variant.id}.{run_kind}",
                        "case_id": case_id,
                        "variant_id": variant.id,
                        "variant_runtime_key": variant.id if uses_variant_runtime else None,
                        "variant_note": variant.note,
                        "run_kind": run_kind,
                        "goal": case_goal(case),
                        "status": status,
                        "run_needed": task_run_needed(status),
                        "priority": task_priority(case_id, variant.id),
                        "base_grid": case.get("base_grid"),
                        "public_case": case.get("public_case"),
                        "cmake": cmake,
                        "runtime": case.get("runtime", []),
                        "required_evidence": case.get("required_evidence", []),
                    }
                )
    return tasks


def import_command(task: dict) -> str:
    parts = [
        "python3",
        "tools/import_lesgo_timing_evidence.py",
        "--case-id",
        task["case_id"],
    ]
    if task.get("variant_runtime_key"):
        parts.extend(["--variant", task["variant_runtime_key"]])
    parts.extend(
        [
            "--backend",
            task["run_kind"],
            "--log",
            f"LOG_PATH_{task['task_id'].replace('.', '_')}",
            "--grid",
            str(task["base_grid"]),
            "--source",
            "SOURCE_ID_OR_RUN_DIRECTORY",
            "--public-record",
            "false",
        ]
    )
    return " ".join(parts)


def render_import_commands(tasks: list[dict], *, needed_only: bool) -> str:
    visible = [task for task in tasks if (task["run_needed"] or not needed_only)]
    lines = [
        "# GPU Validation Timing Import Commands",
        "",
        "# Replace LOG_PATH_* and SOURCE_ID_OR_RUN_DIRECTORY before running.",
        "# Run CPU imports before GPU imports for each variant if adding speed claims later.",
        "",
    ]
    for task in visible:
        lines.append(f"# {task['task_id']} ({task['priority']})")
        lines.append(import_command(task))
        lines.append("")
    return "\n".join(lines)


def render_pair_import_commands(tasks: list[dict], *, needed_only: bool) -> str:
    visible = [task for task in tasks if (task["run_needed"] or not needed_only)]
    grouped: dict[tuple[str, str], dict[str, dict]] = {}
    for task in visible:
        key = (task["case_id"], task["variant_id"])
        grouped.setdefault(key, {})[task["run_kind"]] = task

    lines = [
        "# GPU Validation Paired Timing Import Commands",
        "",
        "# Replace LOG_PATH_* and SOURCE_ID_OR_RUN_DIRECTORY before running.",
        "# These commands parse CPU and GPU logs together and record either",
        "# a speedup claim or an explicit gpu_not_faster paired result.",
        "",
    ]
    for key in sorted(grouped):
        pair = grouped[key]
        if "cpu" not in pair or "gpu" not in pair:
            continue
        cpu = pair["cpu"]
        gpu = pair["gpu"]
        parts = [
            "python3",
            "tools/import_lesgo_timing_pair.py",
            "--case-id",
            cpu["case_id"],
        ]
        if cpu.get("variant_runtime_key"):
            parts.extend(["--variant", cpu["variant_runtime_key"]])
        evidence_items = [
            str(item)
            for item in cpu.get("required_evidence", [])
            if isinstance(item, str) and item
        ]
        parts.extend(
            [
                "--cpu-log",
                f"LOG_PATH_{cpu['task_id'].replace('.', '_')}",
                "--gpu-log",
                f"LOG_PATH_{gpu['task_id'].replace('.', '_')}",
                "--grid",
                str(cpu["base_grid"]),
                "--source",
                "SOURCE_ID_OR_RUN_DIRECTORY",
                "--compare-diagnostics",
                "--correctness-check",
                "MODULE_SPECIFIC_CORRECTNESS_CHECK",
            ]
        )
        for item in evidence_items:
            parts.extend(["--evidence-item", item])
        parts.extend(["--public-record", "false"])
        lines.append(f"# {cpu['case_id']}.{cpu['variant_id']} ({cpu['priority']})")
        lines.append(" ".join(parts))
        lines.append("")
    return "\n".join(lines)


P0_DERECHO_JOB_PREFIX = {
    "les_core_channel": "lesgo_channel",
    "adm_disk": "lesgo_adm",
    "atm_line": "lesgo_atm",
}


def derecho_job_name(task: dict) -> str:
    prefix = P0_DERECHO_JOB_PREFIX[task["case_id"]]
    variant = task["variant_id"]
    variant_part = "" if variant == "baseline" else f"_{variant}"
    return f"{prefix}{variant_part}_{task['run_kind']}"


def render_derecho_submit_commands(tasks: list[dict], *, needed_only: bool) -> str:
    visible = [
        task
        for task in tasks
        if (task["run_needed"] or not needed_only)
        and task["case_id"] in P0_DERECHO_JOB_PREFIX
        and task.get("public_case")
    ]
    grouped: dict[tuple[str, str], dict[str, dict]] = {}
    for task in visible:
        key = (task["case_id"], task["variant_id"])
        grouped.setdefault(key, {})[task["run_kind"]] = task

    lines = [
        "# Derecho Public-Case Submission Commands",
        "",
        "# Run from the repository root on Derecho after cloning/updating this tree.",
        "# These templates cover public p0 paired CPU/GPU cases only.",
        "# Submit jobs from the same case directory sequentially, or copy the case",
        "# tree first; the run scripts clean output directories before execution.",
        "",
    ]
    for key in sorted(grouped):
        pair = grouped[key]
        if "cpu" not in pair or "gpu" not in pair:
            continue
        cpu = pair["cpu"]
        gpu = pair["gpu"]
        public_case = cpu["public_case"]
        lines.append(f"# {cpu['case_id']}.{cpu['variant_id']}: {cpu['variant_note']}")
        if cpu["case_id"] == "atm_line" and cpu["variant_id"] in {
            "structure_off",
            "structure_on",
        }:
            lines.append(
                "# Before submitting this ATM variant, set lesgo.conf/inputATM "
                f"to the intended {cpu['variant_id']} configuration."
            )
        lines.extend(
            [
                f"cd {public_case}",
                "./compile_derecho.sh cpu",
                "./compile_derecho.sh gpu",
                (
                    f"qsub -N {derecho_job_name(gpu)} "
                    f"-v RUN_PROFILE=gpu,RUN_LABEL={gpu['variant_id']}_gpu "
                    "submit_derecho.pbs"
                ),
                "# Wait for the GPU job to finish and archive its log/output before",
                "# submitting the CPU job from the same directory.",
                (
                    f"qsub -N {derecho_job_name(cpu)} "
                    f"-v RUN_PROFILE=cpu,RUN_LABEL={cpu['variant_id']}_cpu "
                    "-l select=1:ncpus=32:mpiprocs=1:mem=120gb submit_derecho.pbs"
                ),
                "cd - >/dev/null",
                "",
            ]
        )
    return "\n".join(lines)


def summarize(tasks: Iterable[dict]) -> dict[str, int]:
    summary: dict[str, int] = {}
    for task in tasks:
        status = task["status"]
        summary[status] = summary.get(status, 0) + 1
    return summary


def filter_tasks(
    tasks: list[dict],
    *,
    priorities: set[str] | None,
    case_ids: set[str] | None,
    run_kind: str | None,
) -> list[dict]:
    visible = tasks
    if priorities:
        visible = [task for task in visible if task["priority"] in priorities]
    if case_ids:
        visible = [task for task in visible if task["case_id"] in case_ids]
    if run_kind:
        visible = [task for task in visible if task["run_kind"] == run_kind]
    return visible


def render_markdown(tasks: list[dict], *, needed_only: bool) -> str:
    visible = [task for task in tasks if (task["run_needed"] or not needed_only)]
    summary = summarize(visible)
    lines = [
        "# GPU Validation Run Plan",
        "",
        "This plan is generated from `docs/gpu_benchmark_manifest.json` and "
        "`docs/gpu_validation_evidence.json`. It is a run schedule, not evidence.",
        "",
        f"Tasks shown: {len(visible)}",
        "",
        "| Status | Tasks |",
        "| --- | ---: |",
    ]
    for status in sorted(summary):
        lines.append(f"| `{status}` | {summary[status]} |")
    lines.extend(
        [
            "",
            "| Priority | Tasks |",
            "| --- | ---: |",
        ]
    )
    by_priority = summarize({"status": task["priority"]} for task in visible)
    for priority in sorted(by_priority):
        lines.append(f"| `{priority}` | {by_priority[priority]} |")
    lines.extend(
        [
            "",
            "| Task | Priority | Goal | Grid | Build | Case/input | Required evidence |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for task in visible:
        cmake = ", ".join(f"`-D{arg}`" for arg in task["cmake"])
        runtime = "; ".join(task["runtime"])
        evidence = ", ".join(f"`{item}`" for item in task["required_evidence"])
        case_input = task["public_case"] or "new/derived case"
        lines.append(
            "| "
            f"`{task['task_id']}`<br>{task['variant_note']}<br>{runtime} | "
            f"`{task['priority']}` | "
            f"{task['goal']} / `{task['status']}` | "
            f"{task['base_grid']} | "
            f"{cmake} | "
            f"{case_input} | "
            f"{evidence} |"
        )
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Expand the GPU validation manifest into paired run tasks."
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE_PATH)
    parser.add_argument(
        "--format",
        choices=["markdown", "json"],
        default="markdown",
        help="output format",
    )
    parser.add_argument(
        "--commands",
        choices=["import", "pair-import", "derecho-submit"],
        help="emit command templates instead of the normal report",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="include tasks that already have complete paired evidence",
    )
    parser.add_argument(
        "--priority",
        action="append",
        choices=sorted(PRIORITY_DESCRIPTIONS),
        help="show only tasks in this priority; can be repeated",
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="case_ids",
        help="show only tasks for this validation case id; can be repeated",
    )
    parser.add_argument(
        "--run-kind",
        choices=["cpu", "gpu"],
        help="show only CPU or GPU tasks",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest = load_json(args.manifest)
    evidence = load_json(args.evidence)
    tasks = expand_tasks(manifest, evidence)
    tasks = filter_tasks(
        tasks,
        priorities=set(args.priority) if args.priority else None,
        case_ids=set(args.case_ids) if args.case_ids else None,
        run_kind=args.run_kind,
    )
    if args.commands == "import":
        print(render_import_commands(tasks, needed_only=not args.all), end="")
    elif args.commands == "pair-import":
        print(render_pair_import_commands(tasks, needed_only=not args.all), end="")
    elif args.commands == "derecho-submit":
        print(render_derecho_submit_commands(tasks, needed_only=not args.all), end="")
    elif args.format == "json":
        payload = {"schema_version": 1, "tasks": tasks}
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(render_markdown(tasks, needed_only=not args.all), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
