#!/usr/bin/env python3
"""Run a prepared Level Set matrix inside an allocated CPU/GPU job."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import time
from pathlib import Path


GENERATED_DIRECTORIES = ("output", "turbineOutput", "turbine", "run-archive")
GENERATED_PATTERNS = (
    "vel.out*", "grid.out", "lvlset_validation.out*", "lvlset_beta.out*",
    "phi.out*", "norm.dat*", "total_time.dat", "*.log",
)


def clean_case(case_dir: Path, preserve_phi: bool) -> None:
    for name in GENERATED_DIRECTORIES:
        path = case_dir / name
        if path.is_dir():
            shutil.rmtree(path)
    for pattern in GENERATED_PATTERNS:
        if preserve_phi and pattern == "phi.out*":
            continue
        for path in case_dir.glob(pattern):
            if path.is_file():
                path.unlink()


def load_executables(path: Path) -> dict[str, dict[str, object]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    result: dict[str, dict[str, object]] = {}
    for profile, value in raw.items():
        if isinstance(value, str):
            result[profile] = {"path": value, "environment": {}}
        elif isinstance(value, dict) and "path" in value:
            result[profile] = {
                "path": value["path"],
                "environment": value.get("environment", {}),
            }
        else:
            raise ValueError(f"invalid executable mapping for {profile}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--executables", type=Path, required=True)
    parser.add_argument(
        "--launcher",
        default="mpiexec -n {nproc}",
        help="MPI launcher template; {nproc} is replaced per task",
    )
    parser.add_argument("--task", action="append", default=[])
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    args = parser.parse_args()

    manifest_path = args.matrix.resolve()
    matrix_root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    executables = load_executables(args.executables.resolve())
    selected = set(args.task)
    status_path = matrix_root / "run_status.json"
    status = {"schema_version": 1, "tasks": {}}
    if status_path.is_file() and not args.rerun:
        status = json.loads(status_path.read_text(encoding="utf-8"))

    failed = False
    for task in manifest["tasks"]:
        task_id = task["id"]
        if selected and task_id not in selected:
            continue
        old = status["tasks"].get(task_id, {})
        if old.get("status") == "passed" and not args.rerun:
            continue
        profile = task["profile"]
        if profile not in executables:
            raise SystemExit(f"missing executable mapping for {profile}")
        executable = Path(str(executables[profile]["path"])).expanduser().resolve()
        if not executable.is_file():
            raise SystemExit(f"missing executable for {profile}: {executable}")
        case_dir = matrix_root / task["case_dir"]
        clean_case(case_dir, preserve_phi=task["geometry"] != "trees")
        (case_dir / "output").mkdir(exist_ok=True)
        mpi_enabled = manifest["build_profiles"][profile]["USE_MPI"] == "ON"
        command = [str(executable)]
        if mpi_enabled:
            command = shlex.split(args.launcher.format(nproc=task["nproc"])) + command
        environment = os.environ.copy()
        environment.update({str(k): str(v) for k, v in task["environment"].items()})
        environment.update(
            {str(k): str(v) for k, v in executables[profile]["environment"].items()}
        )
        log_path = case_dir / "run.log"
        started = time.time()
        with log_path.open("w", encoding="utf-8") as stream:
            completed = subprocess.run(
                command,
                cwd=case_dir,
                env=environment,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
        elapsed = time.time() - started
        passed = completed.returncode == 0 and any(case_dir.glob("vel.out*"))
        status["tasks"][task_id] = {
            "status": "passed" if passed else "failed",
            "returncode": completed.returncode,
            "elapsed_seconds": elapsed,
            "command": command,
            "log": str(log_path.relative_to(matrix_root)),
        }
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        print(f"{task_id}: {status['tasks'][task_id]['status']} ({elapsed:.3f} s)")
        if not passed:
            failed = True
            if not args.keep_going:
                break
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
