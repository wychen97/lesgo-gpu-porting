#!/usr/bin/env python3
"""Smoke-test the GPU validation timing import workflow on a temp ledger."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from repo_paths import ROOT, repo_path

SOURCE_EVIDENCE = repo_path("docs", "gpu_validation_evidence.json")
SOURCE_MATRIX = repo_path("docs", "gpu_validation_matrix.md")

CPU_LOG = """\
Time step information:
  Iteration:       200
  Time step:   0.1000000E-02
  Dimensional time:   0.2000000E+00
  CFL:   0.4100000E+00
Flow field information:
  Velocity divergence metric:   0.1200000E-11
  Kinetic energy:   0.9000000E+00
Simulation wall times (s):
  Iteration:   0.6000000E+00
  Cumulative:   0.1200000E+03
  Forcing:   0.2000000E-02
  Cumulative Forcing:   0.3000000E+00
Sub-component Cumulative Times (s):
  Derivatives:   0.3000000E+01
Simulation wall time (s) :   0.1210000E+03
MPI_EXIT_STATUS=0
"""

GPU_LOG = """\
Time step information:
  Iteration:       200
  Time step:   0.1000000E-02
  Dimensional time:   0.2000000E+00
  CFL:   0.4100000E+00
Flow field information:
  Velocity divergence metric:   0.1200000E-11
  Kinetic energy:   0.9000000E+00
Simulation wall times (s):
  Iteration:   0.6000000E-01
  Cumulative:   0.1200000E+02
  Forcing:   0.2000000E-02
  Cumulative Forcing:   0.3000000E+00
Sub-component Cumulative Times (s):
  Derivatives:   0.3000000E+01
Simulation wall time (s) :   0.1300000E+02
MPI_EXIT_STATUS=0
"""

GPU_NOT_FASTER_LOG = """\
Time step information:
  Iteration:       200
  Time step:   0.1000000E-02
  Dimensional time:   0.2000000E+00
  CFL:   0.4100000E+00
Flow field information:
  Velocity divergence metric:   0.1200000E-11
  Kinetic energy:   0.9000000E+00
Simulation wall times (s):
  Iteration:   0.7000000E+00
  Cumulative:   0.1400000E+03
  Forcing:   0.2000000E-02
  Cumulative Forcing:   0.3000000E+00
Sub-component Cumulative Times (s):
  Derivatives:   0.3000000E+01
Simulation wall time (s) :   0.1410000E+03
MPI_EXIT_STATUS=0
"""

GPU_BAD_DIAGNOSTIC_LOG = """\
Time step information:
  Iteration:       200
  Time step:   0.1000000E-02
  Dimensional time:   0.2000000E+00
  CFL:   0.4100000E+00
Flow field information:
  Velocity divergence metric:   0.1200000E-11
  Kinetic energy:   0.9900000E+00
Simulation wall times (s):
  Iteration:   0.6000000E-01
  Cumulative:   0.1200000E+02
  Forcing:   0.2000000E-02
  Cumulative Forcing:   0.3000000E+00
Sub-component Cumulative Times (s):
  Derivatives:   0.3000000E+01
Simulation wall time (s) :   0.1300000E+02
MPI_EXIT_STATUS=0
"""


def run(command: list[str]) -> str:
    result = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{result.returncode}: {' '.join(command)}\n{result.stdout}"
        )
    return result.stdout


def run_expect(command: list[str], expected_returncode: int) -> str:
    result = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != expected_returncode:
        raise RuntimeError(
            "command exited with "
            f"{result.returncode}, expected {expected_returncode}: "
            f"{' '.join(command)}\n{result.stdout}"
        )
    return result.stdout


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="lesgo_gpu_evidence_check_") as tmp:
        tmp_path = Path(tmp)
        evidence = tmp_path / "evidence.json"
        not_faster_evidence = tmp_path / "not_faster_evidence.json"
        not_faster_matrix = tmp_path / "not_faster_matrix.md"
        cpu_log = tmp_path / "cpu.log"
        gpu_log = tmp_path / "gpu.log"
        gpu_not_faster_log = tmp_path / "gpu_not_faster.log"
        gpu_bad_diagnostic_log = tmp_path / "gpu_bad_diagnostic.log"
        shutil.copyfile(SOURCE_EVIDENCE, evidence)
        shutil.copyfile(SOURCE_EVIDENCE, not_faster_evidence)
        matrix_text = SOURCE_MATRIX.read_text(encoding="utf-8")
        matrix_text = re.sub(
            r"^(\| `iwm_wall_model` \|[^|\n]*\|[^|\n]*\|\s*)`[^`]+`",
            r"\1`implemented-needs-benchmark`",
            matrix_text,
            count=1,
            flags=re.MULTILINE,
        )
        not_faster_matrix.write_text(matrix_text, encoding="utf-8")
        cpu_log.write_text(CPU_LOG, encoding="ascii")
        gpu_log.write_text(GPU_LOG, encoding="ascii")
        gpu_not_faster_log.write_text(GPU_NOT_FASTER_LOG, encoding="ascii")
        gpu_bad_diagnostic_log.write_text(GPU_BAD_DIAGNOSTIC_LOG, encoding="ascii")

        run(
            [
                sys.executable,
                "tools/import_lesgo_timing_evidence.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "hit_inflow",
                "--backend",
                "cpu",
                "--log",
                str(cpu_log),
                "--grid",
                "240^3",
                "--source",
                "temporary import workflow check",
                "--public-record",
                "false",
            ]
        )
        missing_correctness = run_expect(
            [
                sys.executable,
                "tools/import_lesgo_timing_evidence.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "hit_inflow",
                "--backend",
                "gpu",
                "--log",
                str(gpu_log),
                "--grid",
                "240^3",
                "--source",
                "temporary import workflow check",
                "--claim-faster",
            ],
            expected_returncode=1,
        )
        if "requires passing correctness evidence" not in missing_correctness:
            raise RuntimeError(
                "speedup import without correctness was not rejected explicitly"
            )
        run(
            [
                sys.executable,
                "tools/import_lesgo_timing_evidence.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "hit_inflow",
                "--backend",
                "gpu",
                "--log",
                str(gpu_log),
                "--grid",
                "240^3",
                "--source",
                "temporary import workflow check",
                "--correctness-check",
                "velocity divergence and kinetic energy matched CPU reference",
                "--evidence-item",
                "cpu_log",
                "--evidence-item",
                "gpu_log",
                "--evidence-item",
                "field_difference",
                "--evidence-item",
                "speedup_table",
                "--claim-faster",
            ]
        )
        run(
            [
                sys.executable,
                "tools/check_gpu_validation_evidence.py",
                "--evidence",
                str(evidence),
            ]
        )

        run(
            [
                sys.executable,
                "tools/import_lesgo_timing_evidence.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "atm_line",
                "--variant",
                "structure_off",
                "--backend",
                "cpu",
                "--log",
                str(cpu_log),
                "--grid",
                "240^3",
                "--source",
                "temporary variant import workflow check",
                "--public-record",
                "false",
            ]
        )
        run(
            [
                sys.executable,
                "tools/import_lesgo_timing_evidence.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "atm_line",
                "--variant",
                "structure_off",
                "--backend",
                "gpu",
                "--log",
                str(gpu_log),
                "--grid",
                "240^3",
                "--source",
                "temporary variant import workflow check",
                "--correctness-check",
                "structure-off turbine outputs matched CPU reference",
                "--evidence-item",
                "cpu_log",
                "--evidence-item",
                "gpu_log",
                "--evidence-item",
                "turbine_power",
                "--evidence-item",
                "structure_on_off_comparison",
                "--evidence-item",
                "late_step_timing",
                "--evidence-item",
                "cumulative_average",
                "--claim-faster",
            ]
        )
        run(
            [
                sys.executable,
                "tools/check_gpu_validation_evidence.py",
                "--evidence",
                str(evidence),
            ]
        )
        run(
            [
                sys.executable,
                "tools/import_lesgo_timing_pair.py",
                "--evidence",
                str(not_faster_evidence),
                "--case-id",
                "iwm_wall_model",
                "--cpu-log",
                str(cpu_log),
                "--gpu-log",
                str(gpu_not_faster_log),
                "--grid",
                "240^3",
                "--source",
                "temporary paired not-faster workflow check",
                "--compare-diagnostics",
                "--public-record",
                "false",
            ]
        )
        bad_diagnostic = run_expect(
            [
                sys.executable,
                "tools/import_lesgo_timing_pair.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "iwm_wall_model",
                "--cpu-log",
                str(cpu_log),
                "--gpu-log",
                str(gpu_bad_diagnostic_log),
                "--grid",
                "240^3",
                "--source",
                "temporary paired diagnostic mismatch check",
                "--compare-diagnostics",
                "--public-record",
                "false",
            ],
            expected_returncode=1,
        )
        if "kinetic energy mismatch" not in bad_diagnostic:
            raise RuntimeError("diagnostic mismatch was not rejected explicitly")
        run(
            [
                sys.executable,
                "tools/import_lesgo_timing_pair.py",
                "--evidence",
                str(evidence),
                "--case-id",
                "atm_line",
                "--variant",
                "structure_on",
                "--cpu-log",
                str(cpu_log),
                "--gpu-log",
                str(gpu_log),
                "--grid",
                "240^3",
                "--source",
                "temporary paired variant workflow check",
                "--compare-diagnostics",
                "--correctness-check",
                "structure-on turbine outputs matched CPU reference",
                "--evidence-item",
                "turbine_power",
                "--evidence-item",
                "structure_on_off_comparison",
                "--public-record",
                "false",
            ]
        )
        run(
            [
                sys.executable,
                "tools/check_gpu_validation_evidence.py",
                "--evidence",
                str(evidence),
            ]
        )
        matrix_report = run_expect(
            [
                sys.executable,
                "tools/report_gpu_matrix_status_updates.py",
                "--matrix",
                str(not_faster_matrix),
                "--evidence",
                str(not_faster_evidence),
                "--fail-on-suggestion",
            ],
            expected_returncode=1,
        )
        if "`iwm_wall_model`" not in matrix_report:
            raise RuntimeError("matrix update report did not include iwm_wall_model")
        if "`recorded-paired-not-faster`" not in matrix_report:
            raise RuntimeError("matrix update report did not suggest not-faster status")

        data = json.loads(evidence.read_text(encoding="utf-8"))
        not_faster_data = json.loads(not_faster_evidence.read_text(encoding="utf-8"))
        case = next(row for row in data["cases"] if row["id"] == "hit_inflow")
        cpu_avg = case["cpu_runtime"]["cumulative_average_s_per_step"]
        gpu_avg = case["gpu_runtime"]["cumulative_average_s_per_step"]
        speedup = case["speedup_claim"]["speedup"]
        if round(cpu_avg, 8) != 0.6:
            raise RuntimeError(f"unexpected CPU average: {cpu_avg}")
        if round(gpu_avg, 8) != 0.06:
            raise RuntimeError(f"unexpected GPU average: {gpu_avg}")
        if round(speedup, 8) != 10.0:
            raise RuntimeError(f"unexpected speedup: {speedup}")
        if case["evidence_status"] != "paired_speedup_claimed":
            raise RuntimeError(f"unexpected evidence status: {case['evidence_status']}")
        if case["correctness"]["status"] != "passed":
            raise RuntimeError("top-level correctness evidence was not recorded")
        if set(case["correctness"]["evidence_items"]) != {
            "cpu_log",
            "gpu_log",
            "field_difference",
            "speedup_table",
        }:
            raise RuntimeError("top-level evidence items were not recorded")
        if case["parsed_log_summary"]["mpi_exit_status"] != 0:
            raise RuntimeError("MPI exit status was not parsed")
        atm_case = next(row for row in data["cases"] if row["id"] == "atm_line")
        variant = atm_case["variant_runtimes"]["structure_off"]
        variant_speedup = atm_case["variant_speedup_claims"]["structure_off"]["speedup"]
        if round(variant["cpu_runtime"]["cumulative_average_s_per_step"], 8) != 0.6:
            raise RuntimeError("unexpected variant CPU average")
        if round(variant["gpu_runtime"]["cumulative_average_s_per_step"], 8) != 0.06:
            raise RuntimeError("unexpected variant GPU average")
        if round(variant_speedup, 8) != 10.0:
            raise RuntimeError(f"unexpected variant speedup: {variant_speedup}")
        if (
            atm_case["variant_correctness"]["structure_off"]["checks"][0]
            != "structure-off turbine outputs matched CPU reference"
        ):
            raise RuntimeError("variant correctness evidence was not recorded")
        if set(atm_case["variant_correctness"]["structure_off"]["evidence_items"]) != {
            "cpu_log",
            "gpu_log",
            "turbine_power",
            "structure_on_off_comparison",
            "late_step_timing",
            "cumulative_average",
        }:
            raise RuntimeError("variant evidence items were not recorded")
        iwm_case = next(
            row for row in not_faster_data["cases"] if row["id"] == "iwm_wall_model"
        )
        if iwm_case["evidence_status"] != "paired_not_faster":
            raise RuntimeError(
                f"unexpected not-faster status: {iwm_case['evidence_status']}"
            )
        if iwm_case["pair_result"]["outcome"] != "gpu_not_faster":
            raise RuntimeError("not-faster pair result was not recorded")
        if round(iwm_case["pair_result"]["speedup"], 8) != round(0.6 / 0.7, 8):
            raise RuntimeError("unexpected not-faster speedup ratio")
        variant_on = atm_case["variant_pair_results"]["structure_on"]
        if variant_on["outcome"] != "gpu_faster":
            raise RuntimeError("paired variant speedup was not recorded")
        variant_on_checks = atm_case["variant_correctness"]["structure_on"]["checks"]
        if not any(check.startswith("MPI exit status matched") for check in variant_on_checks):
            raise RuntimeError("paired diagnostic MPI correctness check was not recorded")
        if not any(check.startswith("divergence matched") for check in variant_on_checks):
            raise RuntimeError("paired diagnostic divergence check was not recorded")
        if not any(check.startswith("kinetic energy matched") for check in variant_on_checks):
            raise RuntimeError("paired diagnostic kinetic-energy check was not recorded")
        if not {
            "cpu_log",
            "gpu_log",
            "turbine_power",
            "structure_on_off_comparison",
            "late_step_timing",
            "cumulative_average",
            "divergence",
            "kinetic_energy",
        } <= set(atm_case["variant_correctness"]["structure_on"]["evidence_items"]):
            raise RuntimeError("paired diagnostic evidence items were not recorded")

    print("GPU validation import workflow check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
