#!/usr/bin/env python3
"""Run the local readiness checks for the shared optimized LESGO branch."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile

from repo_paths import ROOT

PYTHON_CHECKS = [
    # Branch-wide source, CMake, and production-profile hygiene.
    ("source hygiene", "tools/check_source_hygiene.py"),
    ("repository path helpers", "tools/check_repo_paths.py"),
    ("preprocessor macro inventory", "tools/check_preprocessor_macro_inventory.py"),
    ("CMake option documentation", "tools/check_cmake_option_docs.py"),
    ("CMake cache-variable documentation", "tools/check_cmake_cache_docs.py"),
    ("build-profile CMake arguments", "tools/check_build_profile_args.py"),
    ("CMake comment quality", "tools/check_cmake_comment_quality.py"),
    ("CMake invalid feature combinations", "tools/check_cmake_invalid_feature_combinations.py"),
    ("CMake source groups", "tools/check_cmake_source_groups.py"),
    ("production profile documentation", "tools/check_production_profile_docs.py"),
    ("production wording", "tools/check_production_wording.py"),
    ("active source comment quality", "tools/check_active_source_comment_quality.py"),
    # Shell helper and collaborator handoff documentation.
    ("cluster script documentation", "tools/check_cluster_script_docs.py"),
    ("cluster script headers", "tools/check_cluster_script_headers.py"),
    ("test-case script documentation", "tools/check_test_case_script_docs.py"),
    ("p0 paired case scripts", "tools/check_p0_paired_case_scripts.py"),
    ("contributing readiness documentation", "tools/check_contributing_readiness_docs.py"),
    ("handoff README", "tools/check_handoff_readme.py"),
    ("environment switch documentation", "tools/check_environment_switch_docs.py"),
    ("lesgo.conf GPU coverage documentation", "tools/check_lesgo_conf_coverage_docs.py"),
    ("lesgo.conf validation map", "tools/check_lesgo_conf_validation_map.py"),
    ("lesgo.conf key validation report", "tools/check_lesgo_conf_key_validation_report.py"),
    # Fortran source inventory, navigation, and GPU-module contracts.
    ("source comment hygiene", "tools/check_source_comment_hygiene.py"),
    ("source inventory", "tools/check_source_inventory.py"),
    ("Fortran interface hygiene", "tools/check_fortran_interface_hygiene.py"),
    ("MPI sync guard hygiene", "tools/check_mpi_sync_guard_hygiene.py"),
    ("Fortran interface hygiene self-test", "tools/check_fortran_interface_hygiene_selftest.py"),
    ("Fortran broad-import parser self-test", "tools/check_fortran_broad_imports_selftest.py"),
    ("Fortran broad-import audit", "tools/check_fortran_broad_import_audit.py"),
    ("SGS model constants", "tools/check_sgs_model_constants.py"),
    ("refactor backlog hotspots", "tools/check_refactor_backlog.py"),
    ("large-file navigation maps", "tools/check_navigation_maps.py"),
    ("IWM wall-surface indexing", "tools/check_iwm_surface_indexing.py"),
    ("ATM structure packed paths", "tools/check_atm_structure_packed_paths.py"),
    ("GPU comment labels", "tools/check_gpu_comment_labels.py"),
    ("GPU file headers", "tools/check_gpu_headers.py"),
    ("GPU contract source groups", "tools/check_gpu_contract_source_groups.py"),
    ("Level Set validation variants", "tools/check_level_set_validation_variants.py"),
    ("GPU validation matrix", "tools/check_gpu_validation_matrix.py"),
    ("GPU benchmark manifest", "tools/check_gpu_benchmark_manifest.py"),
    ("GPU validation run plan", "tools/check_gpu_validation_plan.py"),
    ("GPU validation evidence ledger", "tools/check_gpu_validation_evidence.py"),
    ("GPU timing audit consistency", "tools/check_gpu_timing_audit_consistency.py"),
    ("GPU validation import workflow", "tools/check_gpu_validation_import_tools.py"),
    ("p0 archived evidence importer", "tools/check_p0_archived_evidence_importer.py"),
    ("GPU validation runbook", "tools/check_gpu_validation_runbook.py"),
    ("GPU release-objective status report", "tools/check_gpu_release_objective_status_report.py"),
    ("GPU static review buckets", "tools/check_gpu_static_review.py"),
    ("GPU static full inventory report", "tools/check_gpu_static_full_inventory_report.py"),
    ("GPU static candidate review report", "tools/check_gpu_static_candidate_review_report.py"),
    ("documentation paths", "tools/check_doc_paths.py"),
    ("documented Fortran source references", "tools/check_doc_source_refs.py"),
    # Index and wrapper self-consistency checks.
    ("documentation index", "tools/check_docs_index.py"),
    ("tooling index", "tools/check_tools_index.py"),
    ("readiness coverage", "tools/check_readiness_coverage.py"),
    ("readiness documentation", "tools/check_readiness_docs.py"),
]

CMAKE_SMOKE_ARGS = [
    "-DUSE_MPI=OFF",
    "-DUSE_CPS=OFF",
    "-DUSE_ATM=OFF",
    "-DUSE_TURBINES=OFF",
    "-DUSE_LES_GPU=OFF",
    "-DUSE_SCALARS=OFF",
    "-DUSE_LVLSET=OFF",
    "-DUSE_HIT=OFF",
    "-DUSE_CGNS=OFF",
]

HIT_CMAKE_SMOKE_ARGS = [
    "-DUSE_MPI=OFF",
    "-DUSE_CPS=OFF",
    "-DUSE_ATM=OFF",
    "-DUSE_TURBINES=OFF",
    "-DUSE_LES_GPU=ON",
    "-DUSE_GPU_AWARE_MPI=OFF",
    "-DUSE_SCALARS=OFF",
    "-DUSE_LVLSET=OFF",
    "-DUSE_HIT=ON",
    "-DUSE_CGNS=OFF",
]


def run_command(label: str, command: list[str]) -> bool:
    print(f"== {label} ==")
    result = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        print(f"FAILED: {label} exited with {result.returncode}")
        return False
    print(f"OK: {label}")
    return True


def run_cmake_configure(label: str, cmake_args: list[str]) -> bool:
    with tempfile.TemporaryDirectory(prefix="lesgo_cmake_check_") as build_dir:
        return run_command(
            label,
            ["cmake", "-S", ".", "-B", build_dir, *cmake_args],
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the local collaboration-readiness checks."
    )
    parser.add_argument(
        "--with-cmake-configure",
        action="store_true",
        help="also run a CPU/offline CMake configure smoke test",
    )
    parser.add_argument(
        "--with-hit-cmake-configure",
        action="store_true",
        help="also run an offline configure smoke with HIT plus the LES GPU sources enabled",
    )
    parser.add_argument(
        "--with-git-diff-check",
        action="store_true",
        help="also run `git diff --check` with the current Git client",
    )
    args = parser.parse_args()

    ok = True
    for label, script in PYTHON_CHECKS:
        ok = run_command(label, [sys.executable, script]) and ok

    if args.with_git_diff_check:
        ok = run_command("git diff whitespace check", ["git", "diff", "--check"]) and ok

    if args.with_cmake_configure:
        ok = run_cmake_configure(
            "CMake configure smoke",
            CMAKE_SMOKE_ARGS,
        ) and ok

    if args.with_hit_cmake_configure:
        ok = run_cmake_configure(
            "CMake HIT configure smoke",
            HIT_CMAKE_SMOKE_ARGS,
        ) and ok

    if ok:
        print("All requested readiness checks passed.")
        return 0

    print("One or more readiness checks failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
