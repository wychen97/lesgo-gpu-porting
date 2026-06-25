#!/usr/bin/env python3
"""Verify the generated GPU validation run plan covers open benchmark rows."""

from __future__ import annotations

import sys
from pathlib import Path

from gpu_validation_plan import (
    CASE_VARIANTS,
    PRIORITY_DESCRIPTIONS,
    expand_tasks,
    load_json,
    render_derecho_submit_commands,
    render_import_commands,
    render_pair_import_commands,
)
from require_gpu_release_objective import case_gap


MANIFEST_PATH = "docs/gpu_benchmark_manifest.json"
EVIDENCE_PATH = "docs/gpu_validation_evidence.json"
def main() -> int:
    manifest = load_json(Path(MANIFEST_PATH))
    evidence = load_json(Path(EVIDENCE_PATH))
    cases = {
        case["id"]: case
        for case in manifest.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }
    evidence_rows = {
        case["id"]: case
        for case in evidence.get("cases", [])
        if isinstance(case, dict) and "id" in case
    }
    tasks = expand_tasks(manifest, evidence)
    ok = True

    stale = sorted(set(CASE_VARIANTS) - set(cases))
    missing = sorted(set(cases) - set(CASE_VARIANTS))
    if stale:
        print("CASE_VARIANTS contains ids not present in the benchmark manifest:")
        for case_id in stale:
            print(f"  {case_id}")
        ok = False
    if missing:
        print("CASE_VARIANTS is missing benchmark manifest ids:")
        for case_id in missing:
            print(f"  {case_id}")
        ok = False

    task_ids = [task["task_id"] for task in tasks]
    if len(task_ids) != len(set(task_ids)):
        print("Generated validation plan contains duplicate task ids")
        ok = False

    by_case: dict[str, list[dict]] = {}
    for task in tasks:
        by_case.setdefault(task["case_id"], []).append(task)
        if task["run_kind"] not in {"cpu", "gpu"}:
            print(f"task `{task['task_id']}` has invalid run_kind `{task['run_kind']}`")
            ok = False
        if task.get("priority") not in PRIORITY_DESCRIPTIONS:
            print(f"task `{task['task_id']}` has invalid priority `{task.get('priority')}`")
            ok = False
        expected_variant_key = (
            task["variant_id"]
            if len(CASE_VARIANTS.get(task["case_id"], [])) > 1
            else None
        )
        if task.get("variant_runtime_key") != expected_variant_key:
            print(
                f"task `{task['task_id']}` has variant_runtime_key "
                f"`{task.get('variant_runtime_key')}`, expected `{expected_variant_key}`"
            )
            ok = False
        if not task["cmake"]:
            print(f"task `{task['task_id']}` has no CMake settings")
            ok = False
        if not task["required_evidence"]:
            print(f"task `{task['task_id']}` has no required evidence list")
            ok = False

    for case_id, case in cases.items():
        variants = CASE_VARIANTS.get(case_id, [])
        generated = by_case.get(case_id, [])
        if case_id == "lvlset":
            if generated:
                print("LVLSET is excluded but generated run tasks")
                ok = False
            continue
        if not variants:
            print(f"case `{case_id}` has no planned variants")
            ok = False
            continue
        run_kinds = {task["run_kind"] for task in generated}
        if run_kinds != {"cpu", "gpu"}:
            print(
                f"case `{case_id}` must generate both CPU and GPU tasks; "
                f"got {sorted(run_kinds)}"
            )
            ok = False
        expected_task_count = 2 * len(variants)
        if len(generated) != expected_task_count:
            print(
                f"case `{case_id}` expected {expected_task_count} tasks "
                f"for {len(variants)} variants, got {len(generated)}"
            )
            ok = False
        gaps = case_gap(case_id, case, evidence_rows.get(case_id))
        needed_tasks = [task for task in generated if task["run_needed"]]
        if gaps and not needed_tasks:
            print(f"case `{case_id}` has release gaps but no task is run_needed")
            ok = False
        if not gaps and needed_tasks:
            print(f"case `{case_id}` is release-proven but still has run_needed tasks")
            ok = False

    if not ok:
        return 1

    p0_tasks = [task for task in tasks if task["priority"] == "p0-public-small"]
    import_text = render_import_commands(p0_tasks, needed_only=False)
    pair_import_text = render_pair_import_commands(p0_tasks, needed_only=False)
    submit_text = render_derecho_submit_commands(p0_tasks, needed_only=False)
    for label, text, markers in [
        (
            "import commands",
            import_text,
            ["tools/import_lesgo_timing_evidence.py", "LOG_PATH_les_core_channel_baseline_cpu"],
        ),
        (
            "pair-import commands",
            pair_import_text,
            [
                "tools/import_lesgo_timing_pair.py",
                "--cpu-log",
                "--gpu-log",
                "--compare-diagnostics",
                "--correctness-check",
                "MODULE_SPECIFIC_CORRECTNESS_CHECK",
                "--evidence-item",
            ],
        ),
        (
            "derecho-submit commands",
            submit_text,
            [
                "./compile_derecho.sh cpu",
                "./compile_derecho.sh gpu",
                "qsub -N lesgo_channel_gpu",
                "qsub -N lesgo_channel_cpu",
                "-v RUN_PROFILE=gpu,RUN_LABEL=baseline_gpu",
                "-v RUN_PROFILE=cpu,RUN_LABEL=baseline_cpu",
                "-l select=1:ncpus=32:mpiprocs=1:mem=120gb",
            ],
        ),
    ]:
        for marker in markers:
            if marker not in text:
                print(f"generated {label} missing marker: {marker}")
                ok = False

    if not ok:
        return 1

    needed = sum(1 for task in tasks if task["run_needed"])
    priorities = len({task["priority"] for task in tasks})
    print(
        "GPU validation run-plan check passed "
        f"({len(tasks)} tasks, {needed} still needed, {priorities} priorities)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
