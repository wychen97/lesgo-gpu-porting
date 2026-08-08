#!/usr/bin/env python3
"""Compare complete Level Set CPU/GPU validation checkpoints."""

from __future__ import annotations

import argparse
import array
import json
import math
import re
import struct
import sys
from pathlib import Path


LES_FIELDS = (
    "u", "v", "w", "RHSx", "RHSy", "RHSz", "Cs_opt2",
    "F_LM", "F_MM", "F_QN", "F_NN",
)
LES_COEFFICIENT_FIELDS = frozenset({"Cs_opt2", "F_LM", "F_MM", "F_QN", "F_NN"})
LEVEL_SET_FIELDS = (
    "dpdx", "dpdy", "dpdz", "fx", "fy", "fz",
    "divtx", "divty", "divtz",
)
BOGUS_LIMIT = 0.5 * 1234567890.0
FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def rank_file(root: Path, basename: str, rank: int, nproc: int) -> Path:
    suffixed = root / f"{basename}.c{rank}"
    plain = root / basename
    if suffixed.is_file():
        return suffixed
    if rank == 0 and plain.is_file():
        return plain
    return suffixed if nproc > 1 else plain


def record_layout(path: Path) -> tuple[str, int]:
    size = path.stat().st_size
    if size < 8:
        raise ValueError(f"{path} is too small for a Fortran record")
    with path.open("rb") as stream:
        prefix = stream.read(4)
        stream.seek(-4, 2)
        suffix = stream.read(4)
    for endian in ("<", ">"):
        count = struct.unpack(f"{endian}i", prefix)[0]
        trailing = struct.unpack(f"{endian}i", suffix)[0]
        if count == trailing == size - 8 and count % 8 == 0:
            return endian, count
    raise ValueError(f"{path} has unsupported sequential-record markers")


def field_metrics() -> dict[str, object]:
    return {
        "count": 0,
        "compared_count": 0,
        "ignored_padding_count": 0,
        "sentinel_count": 0,
        "nonfinite_count": 0,
        "sentinel_mismatch_count": 0,
        "elementwise_violation_count": 0,
        "max_abs": 0.0,
        "max_rel": 0.0,
        "reference_scale": 0.0,
        "sum_sq_diff": 0.0,
        "sum_sq_reference": 0.0,
    }


def compare_record(
    reference: Path,
    candidate: Path,
    fields: tuple[str, ...],
    *,
    rtol: float,
    atol: float,
    rtol_by_field: dict[str, float] | None = None,
    storage_width: int | None = None,
    physical_width: int | None = None,
) -> dict[str, object]:
    if (storage_width is None) != (physical_width is None):
        raise ValueError("storage_width and physical_width must be provided together")
    if storage_width is not None and not 0 < physical_width <= storage_width:
        raise ValueError("physical_width must be in 1..storage_width")
    left_endian, left_bytes = record_layout(reference)
    right_endian, right_bytes = record_layout(candidate)
    if left_bytes != right_bytes:
        return {
            "passed": False,
            "reason": "record_size_mismatch",
            "reference_bytes": left_bytes,
            "candidate_bytes": right_bytes,
        }
    values = left_bytes // 8
    if values % len(fields):
        return {"passed": False, "reason": "field_count_mismatch", "values": values}
    field_size = values // len(fields)
    metrics = {name: field_metrics() for name in fields}
    native = "<" if sys.byteorder == "little" else ">"
    chunk_bytes = 8 * 1024 * 1024
    value_index = 0
    with reference.open("rb") as left, candidate.open("rb") as right:
        left.seek(4)
        right.seek(4)
        remaining = left_bytes
        while remaining:
            count = min(remaining, chunk_bytes)
            left_values = array.array("d")
            right_values = array.array("d")
            left_values.frombytes(left.read(count))
            right_values.frombytes(right.read(count))
            if left_endian != native:
                left_values.byteswap()
            if right_endian != native:
                right_values.byteswap()
            for offset, (ref, cand) in enumerate(zip(left_values, right_values)):
                absolute_index = value_index + offset
                name = fields[absolute_index // field_size]
                row = metrics[name]
                row["count"] += 1
                field_index = absolute_index % field_size
                if (
                    storage_width is not None
                    and field_index % storage_width >= physical_width
                ):
                    row["ignored_padding_count"] += 1
                    continue
                ref_bogus = abs(ref) >= BOGUS_LIMIT
                cand_bogus = abs(cand) >= BOGUS_LIMIT
                if ref_bogus or cand_bogus:
                    row["sentinel_count"] += 1
                    if ref_bogus != cand_bogus:
                        row["sentinel_mismatch_count"] += 1
                    continue
                if not math.isfinite(ref) or not math.isfinite(cand):
                    row["nonfinite_count"] += 1
                    continue
                row["compared_count"] += 1
                diff = cand - ref
                abs_diff = abs(diff)
                scale = max(abs(ref), abs(cand), 1.0)
                row["max_abs"] = max(row["max_abs"], abs_diff)
                row["max_rel"] = max(row["max_rel"], abs_diff / scale)
                row["reference_scale"] = max(row["reference_scale"], abs(ref))
                row["sum_sq_diff"] += diff * diff
                row["sum_sq_reference"] += ref * ref
                field_rtol = (rtol_by_field or {}).get(name, rtol)
                if not math.isclose(ref, cand, rel_tol=field_rtol, abs_tol=atol):
                    row["elementwise_violation_count"] += 1
            value_index += len(left_values)
            remaining -= count
    passed = True
    for name, row in metrics.items():
        count = row.pop("compared_count")
        sum_sq_diff = row.pop("sum_sq_diff")
        sum_sq_reference = row.pop("sum_sq_reference")
        row["rms_diff"] = math.sqrt(sum_sq_diff / count) if count else 0.0
        row["reference_rms"] = math.sqrt(sum_sq_reference / count) if count else 0.0
        row["relative_l2"] = (
            math.sqrt(sum_sq_diff / sum_sq_reference) if sum_sq_reference else 0.0
        )
        row["compared_count"] = count
        field_rtol = (rtol_by_field or {}).get(name, rtol)
        row["norm_tolerance"] = atol + field_rtol * row["reference_rms"]
        row["norm_close"] = row["rms_diff"] <= row["norm_tolerance"]
        row["passed"] = (
            row["nonfinite_count"] == 0
            and row["sentinel_mismatch_count"] == 0
            and (
                row["elementwise_violation_count"] == 0
                or row["norm_close"]
            )
        )
        passed = passed and row["passed"]
    return {"passed": passed, "field_size": field_size, "fields": metrics}


def read_record_fields(path: Path, fields: tuple[str, ...]) -> dict[str, array.array]:
    endian, byte_count = record_layout(path)
    raw = path.read_bytes()[4:-4]
    values = array.array("d")
    values.frombytes(raw)
    native = "<" if sys.byteorder == "little" else ">"
    if endian != native:
        values.byteswap()
    count = byte_count // 8 // len(fields)
    return {
        name: values[index * count : (index + 1) * count]
        for index, name in enumerate(fields)
    }


def derived_metrics(
    root: Path, *, nx: int, ny: int, nz: int, nproc: int, dx: float, dy: float, dz: float
) -> dict[str, object]:
    ld = 2 * (nx // 2 + 1)
    kinetic_sum = 0.0
    kinetic_count = 0
    integrated_force = [0.0, 0.0, 0.0]
    for rank in range(nproc):
        les = read_record_fields(rank_file(root, "vel.out", rank, nproc), LES_FIELDS)
        aux = read_record_fields(
            rank_file(root, "lvlset_validation.out", rank, nproc), LEVEL_SET_FIELDS
        )
        for component in ("u", "v", "w"):
            for index, value in enumerate(les[component]):
                i = index % ld
                if i >= nx or abs(value) >= BOGUS_LIMIT or not math.isfinite(value):
                    continue
                kinetic_sum += 0.5 * value * value
                kinetic_count += 1
        for direction, component in enumerate(("fx", "fy", "fz")):
            for index, value in enumerate(aux[component]):
                i = index % ld
                if i >= nx or abs(value) >= BOGUS_LIMIT or not math.isfinite(value):
                    continue
                integrated_force[direction] -= value * dx * dy * dz
    return {
        "mean_component_kinetic_energy": kinetic_sum / kinetic_count if kinetic_count else 0.0,
        "kinetic_sample_count": kinetic_count,
        "integrated_force": integrated_force,
    }


def last_log_value(path: Path, label: str) -> float | None:
    if not path.is_file():
        return None
    pattern = re.compile(rf"{label}[^\n]*?({FLOAT_PATTERN})", re.IGNORECASE)
    values = []
    for match in pattern.finditer(path.read_text(encoding="utf-8", errors="replace")):
        values.append(float(match.group(1).replace("D", "E").replace("d", "e")))
    return values[-1] if values else None


def compare_scalar(
    reference: float | None,
    candidate: float | None,
    *,
    rtol: float,
    atol: float,
) -> dict[str, object]:
    """Compare one required finite diagnostic using the field tolerances."""
    row: dict[str, object] = {
        "reference": reference,
        "candidate": candidate,
        "rtol": rtol,
        "atol": atol,
    }
    if reference is None or candidate is None:
        row.update({"passed": False, "reason": "missing_value"})
        return row
    if not math.isfinite(reference) or not math.isfinite(candidate):
        row.update({"passed": False, "reason": "nonfinite_value"})
        return row
    absolute_difference = abs(candidate - reference)
    scale = max(abs(reference), abs(candidate))
    tolerance = atol + rtol * scale
    row.update(
        {
            "absolute_difference": absolute_difference,
            "relative_difference": absolute_difference / max(scale, 1.0),
            "tolerance": tolerance,
            "passed": absolute_difference <= tolerance,
        }
    )
    return row


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--reference-log", type=Path, required=True)
    parser.add_argument("--candidate-log", type=Path, required=True)
    parser.add_argument("--nx", type=positive_int, required=True)
    parser.add_argument("--ny", type=positive_int, required=True)
    parser.add_argument("--nz", type=positive_int, required=True)
    parser.add_argument("--nproc", type=positive_int, default=1)
    parser.add_argument("--dx", type=float, required=True)
    parser.add_argument("--dy", type=float, required=True)
    parser.add_argument("--dz", type=float, required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    parser.add_argument("--coefficient-rtol", type=float, default=2.0e-5)
    parser.add_argument("--beta-rtol", type=float, default=5.0e-6)
    parser.add_argument("--out", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    storage_width = 2 * (args.nx // 2 + 1)
    result: dict[str, object] = {"status": "passed", "ranks": []}
    failures: list[str] = []
    for rank in range(args.nproc):
        rank_result: dict[str, object] = {"rank": rank}
        for key, basename, fields in (
            ("les_checkpoint", "vel.out", LES_FIELDS),
            ("level_set_snapshot", "lvlset_validation.out", LEVEL_SET_FIELDS),
        ):
            reference = rank_file(args.reference, basename, rank, args.nproc)
            candidate = rank_file(args.candidate, basename, rank, args.nproc)
            if not reference.is_file() or not candidate.is_file():
                row = {"passed": False, "reason": "missing_file"}
            else:
                row = compare_record(
                    reference,
                    candidate,
                    fields,
                    rtol=args.rtol,
                    atol=args.atol,
                    rtol_by_field={
                        name: args.coefficient_rtol
                        for name in LES_COEFFICIENT_FIELDS
                    } if key == "les_checkpoint" else None,
                    storage_width=storage_width,
                    physical_width=args.nx,
                )
            rank_result[key] = row
            if not row["passed"]:
                failures.append(f"rank {rank} {basename}")
        beta_reference = rank_file(args.reference, "lvlset_beta.out", rank, args.nproc)
        beta_candidate = rank_file(args.candidate, "lvlset_beta.out", rank, args.nproc)
        if beta_reference.exists() != beta_candidate.exists():
            beta = {"passed": False, "reason": "beta_presence_mismatch"}
            failures.append(f"rank {rank} beta presence")
        elif beta_reference.is_file():
            beta = compare_record(
                beta_reference, beta_candidate, ("Beta", "Tn_all"),
                rtol=args.beta_rtol, atol=args.atol,
                storage_width=storage_width, physical_width=args.nx,
            )
            if not beta["passed"]:
                failures.append(f"rank {rank} Beta")
        else:
            beta = {"passed": True, "status": "not_allocated"}
        rank_result["beta_snapshot"] = beta
        result["ranks"].append(rank_result)

    reference_derived = derived_metrics(
        args.reference, nx=args.nx, ny=args.ny, nz=args.nz, nproc=args.nproc,
        dx=args.dx, dy=args.dy, dz=args.dz,
    )
    candidate_derived = derived_metrics(
        args.candidate, nx=args.nx, ny=args.ny, nz=args.nz, nproc=args.nproc,
        dx=args.dx, dy=args.dy, dz=args.dz,
    )
    kinetic_comparison = compare_scalar(
        reference_derived["mean_component_kinetic_energy"],
        candidate_derived["mean_component_kinetic_energy"],
        rtol=args.rtol,
        atol=args.atol,
    )
    force_comparisons = [
        compare_scalar(reference, candidate, rtol=args.rtol, atol=args.atol)
        for reference, candidate in zip(
            reference_derived["integrated_force"],
            candidate_derived["integrated_force"],
        )
    ]
    sample_count_matches = (
        reference_derived["kinetic_sample_count"]
        == candidate_derived["kinetic_sample_count"]
    )
    result["derived"] = {
        "reference": reference_derived,
        "candidate": candidate_derived,
        "comparison": {
            "kinetic_energy": kinetic_comparison,
            "kinetic_sample_count": {
                "reference": reference_derived["kinetic_sample_count"],
                "candidate": candidate_derived["kinetic_sample_count"],
                "passed": sample_count_matches,
            },
            "integrated_force": force_comparisons,
        },
    }
    if not kinetic_comparison["passed"]:
        failures.append("derived kinetic energy")
    if not sample_count_matches:
        failures.append("derived kinetic sample count")
    for direction, compared in zip("xyz", force_comparisons):
        if not compared["passed"]:
            failures.append(f"integrated IBM force {direction}")
    if args.reference_log and args.candidate_log:
        printed_diagnostics: dict[str, object] = {}
        for label in ("divergence", "kinetic energy"):
            compared = compare_scalar(
                last_log_value(args.reference_log, label),
                last_log_value(args.candidate_log, label),
                rtol=args.rtol,
                atol=args.atol,
            )
            printed_diagnostics[label] = compared
            if not compared["passed"]:
                failures.append(f"printed {label}")
        result["printed_diagnostics"] = printed_diagnostics
    if failures:
        result["status"] = "failed"
        result["failures"] = failures
    output = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(output, encoding="utf-8")
    print(output, end="")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
