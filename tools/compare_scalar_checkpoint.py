#!/usr/bin/env python3
"""Compare scalar checkpoint files written by LESGO."""

from __future__ import annotations

import argparse
import array
import json
import math
import struct
import sys
from pathlib import Path


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def read_fortran_record(path: Path) -> tuple[bytes, str]:
    raw = path.read_bytes()
    if len(raw) < 8:
        raise ValueError(f"{path} is too small to be a Fortran unformatted file")
    first = struct.unpack("<i", raw[:4])[0]
    endian = "<"
    if first < 0 or first + 8 != len(raw):
        first = struct.unpack(">i", raw[:4])[0]
        endian = ">"
    if first < 0 or first + 8 != len(raw):
        raise ValueError(
            f"{path} has unsupported record markers: marker={first}, size={len(raw)}"
        )
    trailing = struct.unpack(f"{endian}i", raw[-4:])[0]
    if trailing != first:
        raise ValueError(f"{path} has mismatched Fortran record markers")
    return raw[4:-4], endian


def read_scalar_checkpoint(
    path: Path, *, nx: int, ny: int, nz: int
) -> dict[str, array.array]:
    ld = 2 * (nx // 2 + 1)
    payload, endian = read_fortran_record(path)
    psi_count = nx * ny
    stored_nz = nz
    theta_count = ld * ny * stored_nz
    rhs_count = ld * ny * stored_nz
    expected = (theta_count + rhs_count + psi_count) * 8
    if len(payload) != expected:
        stored_nz = nz + 1
        theta_count = ld * ny * stored_nz
        rhs_count = ld * ny * stored_nz
        expected = (theta_count + rhs_count + psi_count) * 8
    if len(payload) != expected:
        raise ValueError(f"{path} payload is {len(payload)} bytes, expected {expected}")
    values = array.array("d")
    values.frombytes(payload)
    if (endian == ">") == (sys.byteorder == "little"):
        values.byteswap()
    offset = 0
    theta = values[offset : offset + theta_count]
    offset += theta_count
    rhs_t = values[offset : offset + rhs_count]
    offset += rhs_count
    psi_m = values[offset : offset + psi_count]
    return {"theta": theta, "RHS_T": rhs_t, "psi_m": psi_m}


def compare_array(
    name: str, cpu: array.array, gpu: array.array
) -> dict[str, float | str]:
    if len(cpu) != len(gpu):
        raise ValueError(f"{name} length mismatch: cpu={len(cpu)}, gpu={len(gpu)}")
    cpu_min = math.inf
    cpu_max = -math.inf
    gpu_min = math.inf
    gpu_max = -math.inf
    max_abs_diff = 0.0
    sum_sq = 0.0
    for cpu_value, gpu_value in zip(cpu, gpu):
        cpu_min = min(cpu_min, cpu_value)
        cpu_max = max(cpu_max, cpu_value)
        gpu_min = min(gpu_min, gpu_value)
        gpu_max = max(gpu_max, gpu_value)
        diff = gpu_value - cpu_value
        abs_diff = abs(diff)
        if abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
        sum_sq += diff * diff
    return {
        "name": name,
        "count": len(cpu),
        "cpu_min": cpu_min,
        "cpu_max": cpu_max,
        "gpu_min": gpu_min,
        "gpu_max": gpu_max,
        "max_abs_diff": max_abs_diff,
        "rms_diff": math.sqrt(sum_sq / len(cpu)) if cpu else 0.0,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu", type=Path, required=True)
    parser.add_argument("--gpu", type=Path, required=True)
    parser.add_argument("--nx", type=positive_int, required=True)
    parser.add_argument("--ny", type=positive_int, required=True)
    parser.add_argument("--nz", type=positive_int, required=True)
    parser.add_argument("--max-theta-diff", type=float, default=1.0e-8)
    parser.add_argument("--max-rhs-diff", type=float, default=1.0e-8)
    parser.add_argument("--max-psi-diff", type=float, default=1.0e-8)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    cpu = read_scalar_checkpoint(args.cpu, nx=args.nx, ny=args.ny, nz=args.nz)
    gpu = read_scalar_checkpoint(args.gpu, nx=args.nx, ny=args.ny, nz=args.nz)
    comparisons = [
        compare_array("theta", cpu["theta"], gpu["theta"]),
        compare_array("RHS_T", cpu["RHS_T"], gpu["RHS_T"]),
        compare_array("psi_m", cpu["psi_m"], gpu["psi_m"]),
    ]
    result = {"status": "passed", "comparisons": comparisons}
    thresholds = {
        "theta": args.max_theta_diff,
        "RHS_T": args.max_rhs_diff,
        "psi_m": args.max_psi_diff,
    }
    failures = [
        f"{row['name']} max_abs_diff={row['max_abs_diff']:.6e} > {thresholds[row['name']]:.6e}"
        for row in comparisons
        if row["max_abs_diff"] > thresholds[row["name"]]
    ]
    if failures:
        result["status"] = "failed"
        result["failures"] = failures
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
