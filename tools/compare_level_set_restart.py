#!/usr/bin/env python3
"""Compare continuous and split Level Set model-4/5 restart runs."""

from __future__ import annotations

import argparse
import array
import json
import struct
import sys
from pathlib import Path

from compare_level_set_checkpoints import (
    LES_FIELDS,
    LEVEL_SET_FIELDS,
    compare_record,
    derived_metrics,
    rank_file,
)


def read_records(path: Path) -> tuple[str, list[bytes]]:
    raw = path.read_bytes()
    for endian in ("<", ">"):
        offset = 0
        records: list[bytes] = []
        try:
            while offset < len(raw):
                count = struct.unpack_from(f"{endian}i", raw, offset)[0]
                start = offset + 4
                stop = start + count
                trailing = struct.unpack_from(f"{endian}i", raw, stop)[0]
                if count < 0 or trailing != count:
                    raise ValueError
                records.append(raw[start:stop])
                offset = stop + 4
        except (struct.error, ValueError):
            continue
        if offset == len(raw):
            return endian, records
    raise ValueError(f"unsupported restart record structure: {path}")


def restart_values(path: Path) -> tuple[tuple[object, ...], array.array]:
    endian, records = read_records(path)
    if len(records) != 2 or len(records[0]) != 6 * 4 + 4 * 8:
        raise ValueError(f"unexpected Level Set restart schema: {path}")
    header = struct.unpack(f"{endian}6i4d", records[0])
    values = array.array("d")
    values.frombytes(records[1])
    native = "<" if sys.byteorder == "little" else ">"
    if endian != native:
        values.byteswap()
    return header, values


def compare_restart(
    reference: Path,
    candidate: Path,
    rtol: float,
    atol: float,
    storage_width: int,
    physical_width: int,
) -> dict[str, object]:
    left_header, left = restart_values(reference)
    right_header, right = restart_values(candidate)
    if left_header != right_header or len(left) != len(right):
        return {
            "passed": False,
            "reason": "header_or_size_mismatch",
            "reference_header": left_header,
            "candidate_header": right_header,
        }
    half = len(left) // 2
    passed = True
    fields: dict[str, object] = {}
    for index, name in enumerate(("Beta", "Tn_all")):
        max_abs = 0.0
        max_rel = 0.0
        ignored_padding_count = 0
        for field_index, (ref, cand) in enumerate(
            zip(
                left[index * half : (index + 1) * half],
                right[index * half : (index + 1) * half],
            )
        ):
            if field_index % storage_width >= physical_width:
                ignored_padding_count += 1
                continue
            diff = abs(cand - ref)
            max_abs = max(max_abs, diff)
            max_rel = max(max_rel, diff / max(abs(ref), abs(cand), 1.0))
            if diff > atol + rtol * max(abs(ref), abs(cand)):
                passed = False
        fields[name] = {
            "max_abs": max_abs,
            "max_rel": max_rel,
            "ignored_padding_count": ignored_padding_count,
        }
    return {"passed": passed, "header": left_header, "fields": fields}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--continuous", type=Path, required=True)
    parser.add_argument("--restarted", type=Path, required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--dx", type=float, required=True)
    parser.add_argument("--dy", type=float, required=True)
    parser.add_argument("--dz", type=float, required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    storage_width = 2 * (args.nx // 2 + 1)
    result: dict[str, object] = {"status": "passed", "ranks": []}
    failures: list[str] = []
    for rank in range(args.nproc):
        row: dict[str, object] = {"rank": rank}
        for key, basename, fields in (
            ("les_checkpoint", "vel.out", LES_FIELDS),
            ("validation_snapshot", "lvlset_validation.out", LEVEL_SET_FIELDS),
            ("beta_snapshot", "lvlset_beta.out", ("Beta", "Tn_all")),
        ):
            compared = compare_record(
                rank_file(args.continuous, basename, rank, args.nproc),
                rank_file(args.restarted, basename, rank, args.nproc),
                fields,
                rtol=args.rtol,
                atol=args.atol,
                storage_width=storage_width,
                physical_width=args.nx,
            )
            row[key] = compared
            if not compared["passed"]:
                failures.append(f"rank {rank} {basename}")
        restart = compare_restart(
            rank_file(args.continuous, "lvlset_sgs_restart.out", rank, args.nproc),
            rank_file(args.restarted, "lvlset_sgs_restart.out", rank, args.nproc),
            args.rtol,
            args.atol,
            storage_width,
            args.nx,
        )
        row["restart_sidecar"] = restart
        if not restart["passed"]:
            failures.append(f"rank {rank} restart sidecar")
        result["ranks"].append(row)
    result["derived"] = {
        "continuous": derived_metrics(
            args.continuous, nx=args.nx, ny=args.ny, nz=args.nz, nproc=args.nproc,
            dx=args.dx, dy=args.dy, dz=args.dz,
        ),
        "restarted": derived_metrics(
            args.restarted, nx=args.nx, ny=args.ny, nz=args.nz, nproc=args.nproc,
            dx=args.dx, dy=args.dy, dz=args.dz,
        ),
    }
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
