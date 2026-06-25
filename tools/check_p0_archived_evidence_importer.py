#!/usr/bin/env python3
"""Smoke-test importing archived p0 public-case run logs."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_EVIDENCE = ROOT / "docs/gpu_validation_evidence.json"

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


def archive_path(root: Path, case_dir: str, label: str) -> Path:
    return root / case_dir / "run-archives" / label / f"lesgo_{label}.log"


def write_archive(root: Path, case_dir: str, label: str, text: str) -> None:
    path = archive_path(root, case_dir, label)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="ascii")


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


def case_by_id(data: dict, case_id: str) -> dict:
    for row in data["cases"]:
        if row["id"] == case_id:
            return row
    raise RuntimeError(f"missing case id: {case_id}")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="lesgo_p0_archives_") as tmp:
        tmp_path = Path(tmp)
        evidence = tmp_path / "evidence.json"
        archive_root = tmp_path / "repo"
        shutil.copyfile(SOURCE_EVIDENCE, evidence)

        write_archive(
            archive_root,
            "test-cases/channel_flow",
            "baseline_cpu",
            CPU_LOG,
        )
        write_archive(
            archive_root,
            "test-cases/channel_flow",
            "baseline_gpu",
            GPU_LOG,
        )
        write_archive(
            archive_root,
            "test-cases/adm_disk",
            "baseline_cpu",
            CPU_LOG,
        )
        write_archive(
            archive_root,
            "test-cases/adm_disk",
            "baseline_gpu",
            GPU_LOG,
        )

        run(
            [
                sys.executable,
                "tools/import_p0_archived_evidence.py",
                "--evidence",
                str(evidence),
                "--archive-root",
                str(archive_root),
                "--source",
                "temporary archived p0 import check",
                "--case",
                "les_core_channel",
                "--public-record",
                "false",
            ]
        )
        missing_module_check = run_expect(
            [
                sys.executable,
                "tools/import_p0_archived_evidence.py",
                "--evidence",
                str(evidence),
                "--archive-root",
                str(archive_root),
                "--source",
                "temporary archived p0 import check",
                "--case",
                "adm_disk",
            ],
            expected_returncode=1,
        )
        if "requires --module-check adm_disk=TEXT" not in missing_module_check:
            raise RuntimeError("ADM import did not require module-specific evidence")

        run(
            [
                sys.executable,
                "tools/import_p0_archived_evidence.py",
                "--evidence",
                str(evidence),
                "--archive-root",
                str(archive_root),
                "--source",
                "temporary archived p0 import check",
                "--case",
                "adm_disk",
                "--module-check",
                "adm_disk=turbine forcing and disk velocity matched CPU reference",
                "--public-record",
                "false",
            ]
        )
        data = json.loads(evidence.read_text(encoding="utf-8"))
        channel = case_by_id(data, "les_core_channel")
        adm = case_by_id(data, "adm_disk")
        for row in [channel, adm]:
            if row["evidence_status"] != "paired_speedup_claimed":
                raise RuntimeError(f"unexpected status for {row['id']}")
            if round(row["speedup_claim"]["speedup"], 8) != 10.0:
                raise RuntimeError(f"unexpected speedup for {row['id']}")
            correctness = row.get("correctness")
            if not correctness or correctness.get("status") != "passed":
                raise RuntimeError(f"missing correctness evidence for {row['id']}")
        if not {
            "cpu_log",
            "gpu_log",
            "divergence",
            "kinetic_energy",
            "late_step_timing",
            "cumulative_average",
        } <= set(channel["correctness"]["evidence_items"]):
            raise RuntimeError("channel required evidence items were not imported")
        if not {
            "cpu_log",
            "gpu_log",
            "turbine_forcing_output",
            "power_or_disk_velocity",
            "late_step_timing",
            "cumulative_average",
        } <= set(adm["correctness"]["evidence_items"]):
            raise RuntimeError("ADM required evidence items were not imported")
        if not any(
            "turbine forcing and disk velocity" in check
            for check in adm["correctness"]["checks"]
        ):
            raise RuntimeError("ADM module-specific correctness check was not recorded")

    print("P0 archived evidence importer check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
