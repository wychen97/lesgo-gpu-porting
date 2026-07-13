#!/usr/bin/env python3
"""Compare uninterrupted and split/restarted LESGO runs.

The comparison focuses on the restart seam and final checkpoint agreement.
It intentionally uses only the Python standard library so it can run on login
nodes without a project-specific Python environment. Channel, ADM disk, and
ATM line cases are supported.
"""

from __future__ import annotations

import argparse
import array
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
from typing import Iterable


SIGNALS = ("power", "thrust", "RotSpeed")
FINAL_FILES = ("grid.out", "total_time.dat")
TURBINE_FILES = ("restart", "actuatorPoints", "structure_restart")
LES_RECORD_FIELDS = (
    "u",
    "v",
    "w",
    "RHSx",
    "RHSy",
    "RHSz",
    "Cs_opt2",
    "F_LM",
    "F_MM",
    "F_QN",
    "F_NN",
)
LES_HISTORY_FIELDS = frozenset(
    ("RHSx", "RHSy", "RHSz", "F_LM", "F_MM", "F_QN", "F_NN")
)


def read_numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            tokens = line.split()
            if not tokens:
                continue
            try:
                rows.append([float(token.replace("D", "E")) for token in tokens])
            except ValueError:
                continue
    return rows


def read_numeric_tokens(path: Path) -> list[float]:
    values: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for token in stream.read().split():
            try:
                values.append(float(token.replace("D", "E")))
            except ValueError:
                continue
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def compare_rows(
    reference: list[list[float]], candidate: list[list[float]], rtol: float, atol: float
) -> dict[str, object]:
    if len(reference) != len(candidate):
        return {
            "passed": False,
            "reason": "row_count_mismatch",
            "reference_rows": len(reference),
            "candidate_rows": len(candidate),
        }

    max_abs = 0.0
    max_rel = 0.0
    max_location = [-1, -1]
    passed = True
    for row_index, (left, right) in enumerate(zip(reference, candidate)):
        if len(left) != len(right):
            return {
                "passed": False,
                "reason": "column_count_mismatch",
                "row": row_index,
            }
        for column, (a, b) in enumerate(zip(left, right)):
            abs_error = abs(a - b)
            scale = max(abs(a), abs(b), 1.0)
            rel_error = abs_error / scale
            if abs_error > max_abs:
                max_abs = abs_error
                max_location = [row_index, column]
            max_rel = max(max_rel, rel_error)
            if not math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
                passed = False

    return {
        "passed": passed,
        "rows": len(reference),
        "max_abs": max_abs,
        "max_rel": max_rel,
        "max_location": max_location,
    }


def compare_numeric_tokens(
    reference: Path, candidate: Path, rtol: float, atol: float
) -> dict[str, object]:
    left = read_numeric_tokens(reference)
    right = read_numeric_tokens(candidate)
    result = compare_rows([left], [right], rtol, atol)
    result["byte_identical"] = sha256(reference) == sha256(candidate)
    result["values"] = len(left)
    return result


def fortran_record_layout(path: Path) -> tuple[str, int]:
    size = path.stat().st_size
    if size < 8:
        raise ValueError("file is too small to contain a Fortran record")
    with path.open("rb") as stream:
        prefix = stream.read(4)
        stream.seek(-4, 2)
        suffix = stream.read(4)
    for endian in ("<", ">"):
        record_bytes = struct.unpack(f"{endian}i", prefix)[0]
        trailer_bytes = struct.unpack(f"{endian}i", suffix)[0]
        if record_bytes == trailer_bytes == size - 8 and record_bytes % 8 == 0:
            return endian, record_bytes
    raise ValueError("unsupported or corrupt sequential Fortran record")


def compare_fortran_real64_record(
    reference: Path,
    candidate: Path,
    rtol: float,
    atol: float,
    history_rtol: float,
    history_atol: float,
    cs_rtol: float,
    cs_atol: float,
) -> dict[str, object]:
    try:
        left_endian, left_bytes = fortran_record_layout(reference)
        right_endian, right_bytes = fortran_record_layout(candidate)
    except ValueError as exc:
        return {"passed": False, "reason": str(exc)}
    if left_bytes != right_bytes:
        return {
            "passed": False,
            "reason": "record_size_mismatch",
            "reference_bytes": left_bytes,
            "candidate_bytes": right_bytes,
        }

    max_abs = 0.0
    max_rel = 0.0
    max_index = -1
    passed = True
    value_index = 0
    value_count = left_bytes // 8
    if value_count % len(LES_RECORD_FIELDS) != 0:
        return {
            "passed": False,
            "reason": "unexpected_les_record_value_count",
            "values": value_count,
        }
    field_size = value_count // len(LES_RECORD_FIELDS)
    field_metrics = {
        name: {"passed": True, "max_abs": 0.0, "max_rel": 0.0, "max_index": -1}
        for name in LES_RECORD_FIELDS
    }
    chunk_bytes = 8 * 1024 * 1024
    native_endian = "<" if sys.byteorder == "little" else ">"
    with reference.open("rb") as left_stream, candidate.open("rb") as right_stream:
        left_stream.seek(4)
        right_stream.seek(4)
        remaining = left_bytes
        while remaining:
            count = min(remaining, chunk_bytes)
            left_values = array.array("d")
            right_values = array.array("d")
            left_values.frombytes(left_stream.read(count))
            right_values.frombytes(right_stream.read(count))
            if left_endian != native_endian:
                left_values.byteswap()
            if right_endian != native_endian:
                right_values.byteswap()
            for offset, (left, right) in enumerate(zip(left_values, right_values)):
                current_index = value_index + offset
                field_name = LES_RECORD_FIELDS[current_index // field_size]
                field_result = field_metrics[field_name]
                abs_error = abs(left - right)
                scale = max(abs(left), abs(right), 1.0)
                rel_error = abs_error / scale
                if abs_error > max_abs:
                    max_abs = abs_error
                    max_index = value_index + offset
                max_rel = max(max_rel, rel_error)
                if abs_error > field_result["max_abs"]:
                    field_result["max_abs"] = abs_error
                    field_result["max_index"] = current_index % field_size
                field_result["max_rel"] = max(field_result["max_rel"], rel_error)
                if field_name == "Cs_opt2":
                    field_rtol = cs_rtol
                    field_atol = cs_atol
                elif field_name in LES_HISTORY_FIELDS:
                    field_rtol = history_rtol
                    field_atol = history_atol
                else:
                    field_rtol = rtol
                    field_atol = atol
                if not math.isclose(
                    left, right, rel_tol=field_rtol, abs_tol=field_atol
                ):
                    passed = False
                    field_result["passed"] = False
            value_index += len(left_values)
            remaining -= count

    return {
        "passed": passed,
        "byte_identical": sha256(reference) == sha256(candidate),
        "values": value_index,
        "max_abs": max_abs,
        "max_rel": max_rel,
        "max_index": max_index,
        "fields": field_metrics,
    }


def turbine_directories(root: Path) -> Iterable[Path]:
    output = root / "turbineOutput"
    if not output.is_dir():
        return []
    return sorted(path for path in output.iterdir() if path.is_dir())


def adm_output_files(root: Path) -> Iterable[Path]:
    output = root / "turbine"
    if not output.is_dir():
        return []
    return sorted(path for path in output.rglob("*") if path.is_file())


def compare_run_pair(args: argparse.Namespace) -> dict[str, object]:
    continuous = args.continuous.resolve()
    restarted = args.restarted.resolve()
    report: dict[str, object] = {
        "case_kind": args.case_kind,
        "continuous": str(continuous),
        "restarted": str(restarted),
        "seam_step": args.seam_step,
        "tolerances": {
            "solution_rtol": args.rtol,
            "solution_atol": args.atol,
            "history_rtol": args.history_rtol,
            "history_atol": args.history_atol,
            "cs_rtol": args.cs_rtol,
            "cs_atol": args.cs_atol,
            "state_rtol": args.state_rtol,
            "state_atol": args.state_atol,
        },
        "signals": {},
        "final_files": {},
        "turbines": {},
        "adm_files": {},
    }

    turbine_names = sorted(
        {path.name for path in turbine_directories(continuous)}
        | {path.name for path in turbine_directories(restarted)}
    )
    all_passed = args.case_kind != "atm" or bool(turbine_names)

    for turbine_name in turbine_names:
        turbine_report: dict[str, object] = {"signals": {}, "files": {}}
        for signal in SIGNALS:
            reference_path = continuous / "turbineOutput" / turbine_name / signal
            candidate_path = restarted / "turbineOutput" / turbine_name / signal
            if not reference_path.is_file() or not candidate_path.is_file():
                result: dict[str, object] = {"passed": False, "reason": "missing_file"}
            else:
                result = compare_rows(
                    read_numeric_rows(reference_path),
                    read_numeric_rows(candidate_path),
                    args.rtol,
                    args.atol,
                )
                rows = read_numeric_rows(candidate_path)
                if 0 < args.seam_step < len(rows):
                    result["candidate_seam_before"] = rows[args.seam_step - 1]
                    result["candidate_seam_after"] = rows[args.seam_step]
            turbine_report["signals"][signal] = result
            all_passed = all_passed and bool(result.get("passed"))

        for filename in TURBINE_FILES:
            reference_path = continuous / "turbineOutput" / turbine_name / filename
            candidate_path = restarted / "turbineOutput" / turbine_name / filename
            if not reference_path.exists() and not candidate_path.exists():
                continue
            if reference_path.is_file() and candidate_path.is_file():
                result = compare_numeric_tokens(
                    reference_path,
                    candidate_path,
                    args.state_rtol,
                    args.state_atol,
                )
            else:
                result = {"passed": False, "reason": "missing_file"}
            turbine_report["files"][filename] = result
            all_passed = all_passed and bool(result.get("passed"))
        report["turbines"][turbine_name] = turbine_report

    adm_names = sorted(
        {
            str(path.relative_to(continuous / "turbine"))
            for path in adm_output_files(continuous)
        }
        | {
            str(path.relative_to(restarted / "turbine"))
            for path in adm_output_files(restarted)
        }
    )
    if args.case_kind == "adm" and not adm_names:
        all_passed = False
    for filename in adm_names:
        reference_path = continuous / "turbine" / filename
        candidate_path = restarted / "turbine" / filename
        if not reference_path.is_file() or not candidate_path.is_file():
            result = {"passed": False, "reason": "missing_file"}
        elif Path(filename).name == "u_d_T.dat":
            result = compare_numeric_tokens(
                reference_path,
                candidate_path,
                args.state_rtol,
                args.state_atol,
            )
        else:
            result = compare_rows(
                read_numeric_rows(reference_path),
                read_numeric_rows(candidate_path),
                args.rtol,
                args.atol,
            )
        report["adm_files"][filename] = result
        all_passed = all_passed and bool(result.get("passed"))

    velocity_files = sorted(
        {path.name for path in continuous.glob("vel.out.c*")}
        | {path.name for path in restarted.glob("vel.out.c*")}
    )
    if not velocity_files:
        report["final_files"]["vel.out.c*"] = {
            "passed": False,
            "reason": "missing_file",
        }
        all_passed = False

    for filename in velocity_files + list(FINAL_FILES):
        reference_path = continuous / filename
        candidate_path = restarted / filename
        if not reference_path.is_file() or not candidate_path.is_file():
            result = {"passed": False, "reason": "missing_file"}
        elif filename.startswith("vel.out.c"):
            result = compare_fortran_real64_record(
                reference_path,
                candidate_path,
                args.rtol,
                args.atol,
                args.history_rtol,
                args.history_atol,
                args.cs_rtol,
                args.cs_atol,
            )
        elif filename == "total_time.dat":
            result = compare_numeric_tokens(
                reference_path, candidate_path, args.rtol, args.atol
            )
        else:
            identical = sha256(reference_path) == sha256(candidate_path)
            result = {"passed": identical, "byte_identical": identical}
        report["final_files"][filename] = result
        all_passed = all_passed and bool(result.get("passed"))

    report["passed"] = all_passed
    return report


def write_markdown(report: dict[str, object], path: Path) -> None:
    lines = [
        "# LESGO Restart Continuity Report",
        "",
        f"Case kind: `{report['case_kind']}`",
        "",
        f"Overall: **{'PASS' if report['passed'] else 'FAIL'}**",
        "",
        "| Turbine | Signal | Pass | Max absolute error | Max relative error |",
        "|---|---|---:|---:|---:|",
    ]
    for turbine, turbine_report in report["turbines"].items():
        for signal, result in turbine_report["signals"].items():
            lines.append(
                f"| {turbine} | {signal} | {result.get('passed', False)} | "
                f"{result.get('max_abs', 'n/a')} | {result.get('max_rel', 'n/a')} |"
            )
    if report["adm_files"]:
        lines.extend(
            [
                "",
                "## ADM Outputs",
                "",
                "| File | Pass | Max absolute error | Max relative error |",
                "|---|---:|---:|---:|",
            ]
        )
        for filename, result in report["adm_files"].items():
            lines.append(
                f"| {filename} | {result.get('passed', False)} | "
                f"{result.get('max_abs', 'n/a')} | "
                f"{result.get('max_rel', 'n/a')} |"
            )
    lines.extend(["", "## Final Checkpoints", "", "| File | Pass | Byte identical |", "|---|---:|---:|"])
    for filename, result in report["final_files"].items():
        lines.append(
            f"| {filename} | {result.get('passed', False)} | "
            f"{result.get('byte_identical', False)} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--continuous", type=Path, required=True)
    parser.add_argument("--restarted", type=Path, required=True)
    parser.add_argument(
        "--case-kind", choices=("channel", "adm", "atm"), default="atm"
    )
    parser.add_argument("--seam-step", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--history-rtol", type=float, default=1.0e-6)
    parser.add_argument("--history-atol", type=float, default=1.0e-6)
    parser.add_argument(
        "--cs-rtol",
        type=float,
        default=2.0e-4,
        help="Tolerance for derived Cs_opt2; physical/history fields stay stricter.",
    )
    parser.add_argument("--cs-atol", type=float, default=2.0e-4)
    parser.add_argument("--state-rtol", type=float, default=1.0e-4)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    parser.add_argument("--state-atol", type=float, default=1.0e-4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = compare_run_pair(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(report, args.out.with_suffix(".md"))
    print(json.dumps({"passed": report["passed"], "report": str(args.out)}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
