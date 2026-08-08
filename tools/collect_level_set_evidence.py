#!/usr/bin/env python3
"""Collect CPU/bridge/GPU Level Set pair results into one evidence JSON file."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def conf_value(path: Path, key: str) -> float:
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.split("!", 1)[0].strip()
        if "=" not in stripped:
            continue
        left, right = stripped.split("=", 1)
        if left.strip().lower() == key.lower():
            return float(right.split()[0].replace("D", "E").replace("d", "e"))
    raise ValueError(f"{path} does not define {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--cluster", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    args = parser.parse_args()

    manifest_path = args.matrix.resolve()
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    status = json.loads((root / "run_status.json").read_text(encoding="utf-8"))
    grouped: dict[str, dict[str, dict[str, object]]] = {}
    for task in manifest["tasks"]:
        grouped.setdefault(task["variant"], {})[task["profile"]] = task

    pairs: list[dict[str, object]] = []
    overall = "passed"
    comparator = Path(__file__).with_name("compare_level_set_checkpoints.py")
    for variant, tasks in sorted(grouped.items()):
        for candidate_profile, candidate in sorted(tasks.items()):
            if candidate_profile.startswith("cpu_"):
                continue
            cpu_profile = "cpu_nompi" if candidate_profile.endswith("nompi") else "cpu_mpi"
            reference = tasks.get(cpu_profile)
            if not reference:
                continue
            reference_status = status["tasks"].get(reference["id"], {}).get("status")
            candidate_status = status["tasks"].get(candidate["id"], {}).get("status")
            row: dict[str, object] = {
                "variant": variant,
                "reference_profile": cpu_profile,
                "candidate_profile": candidate_profile,
                "reference_run_status": reference_status,
                "candidate_run_status": candidate_status,
            }
            if reference_status != "passed" or candidate_status != "passed":
                row["status"] = "not_compared"
                overall = "incomplete"
                pairs.append(row)
                continue
            reference_dir = root / reference["case_dir"]
            candidate_dir = root / candidate["case_dir"]
            conf = reference_dir / "lesgo.conf"
            nx = int(conf_value(conf, "Nx"))
            ny = int(conf_value(conf, "Ny"))
            nz = int(conf_value(conf, "Nz"))
            lx = conf_value(conf, "Lx")
            dx = lx / nx
            comparison_path = root / "comparisons" / f"{variant}__{candidate_profile}.json"
            command = [
                sys.executable, str(comparator),
                "--reference", str(reference_dir),
                "--candidate", str(candidate_dir),
                "--reference-log", str(reference_dir / "run.log"),
                "--candidate-log", str(candidate_dir / "run.log"),
                "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
                "--nproc", str(reference["nproc"]),
                "--dx", str(dx), "--dy", str(dx), "--dz", str(dx),
                "--rtol", str(args.rtol), "--atol", str(args.atol),
                "--out", str(comparison_path),
            ]
            completed = subprocess.run(command, stdout=subprocess.DEVNULL, check=False)
            comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
            row["status"] = comparison["status"]
            row["comparison"] = str(comparison_path.relative_to(root))
            if completed.returncode != 0:
                overall = "failed"
            pairs.append(row)

    evidence = {
        "schema_version": 1,
        "cluster": args.cluster,
        "commit": args.commit,
        "matrix": str(manifest_path),
        "status": overall,
        "pairs": pairs,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": overall, "pairs": len(pairs), "out": str(args.out)}, indent=2))
    return 0 if overall == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
