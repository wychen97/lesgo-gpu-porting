#!/usr/bin/env python3
"""Verify ATM packed/batched paths carry structural quantities.

This check protects the structure-on ATM performance path.  The packed gather
and batched Cl-correction paths must not fall back only because
LESGO_ATM_STRUCTURE is enabled; instead they must carry the extra structural
moment state needed by the torsion solve.  It also protects the per-turbine
CUDA induced-velocity correction in actuator_turbine_model.f90; that correction
writes the same du/Uinf/uy state for rigid and structural consumers.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INTERFACE_PATH = ROOT / "atm_lesgo_interface.f90"
ACTUATOR_PATH = ROOT / "actuator_turbine_model.f90"


def subroutine_block(text: str, name: str) -> str:
    start_re = re.compile(rf"^\s*subroutine\s+{re.escape(name)}\b", re.MULTILINE)
    match = start_re.search(text)
    if not match:
        raise AssertionError(f"missing subroutine {name}")

    end_re = re.compile(rf"^\s*end\s+subroutine\s+{re.escape(name)}\b", re.MULTILINE)
    end_match = end_re.search(text, match.end())
    if not end_match:
        raise AssertionError(f"missing end subroutine {name}")

    return text[match.start() : end_match.end()]


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def check_dispatch(text: str, failures: list[str]) -> None:
    block = subroutine_block(text, "atm_lesgo_mpi_gather")
    require(
        ".not. atm_structure_enabled()" not in block,
        "atm_lesgo_mpi_gather still has a structure-off-only gate",
        failures,
    )
    require(
        "call atm_lesgo_mpi_gather_packed()" in block,
        "atm_lesgo_mpi_gather does not dispatch to packed gather",
        failures,
    )
    require(
        "Structure-on runs add Cm/pitchingMoment" in block,
        "atm_lesgo_mpi_gather lacks the structure-on packed payload comment",
        failures,
    )


def check_packed_payload(text: str, failures: list[str]) -> None:
    required = {
        "atm_lesgo_mpi_gather_packed": [
            "struct_active = atm_structure_enabled()",
            "size(turbineArray(i) % Cm)",
            "size(turbineArray(i) % pitchingMoment)",
            "packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cm",
            "turbineArray(i) % Cm = reshape(packed_recv(pos:pos+nitem-1)",
            "turbineArray(i) % pitchingMoment = reshape(",
        ],
        "atm_lesgo_mpi_gather_slim_gpu": [
            "struct_active = atm_structure_enabled()",
            "size(turbineArray(i) % Cm)",
            "size(turbineArray(i) % pitchingMoment)",
            "call atm_pack_rank3(turbineArray(i) % Cm, packed_send, pos)",
            "call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cm)",
            "turbineArray(i) % pitchingMoment",
        ],
        "atm_lesgo_mpi_gather_slim_batch_gpu": [
            "struct_active = atm_structure_enabled()",
            "size(turbineArray(i) % Cm)",
            "size(turbineArray(i) % pitchingMoment)",
            "call atm_pack_rank3(turbineArray(i) % Cm, packed_send, pos)",
            "call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cm)",
            "turbineArray(i) % pitchingMoment",
        ],
        "atm_lesgo_mpi_gather_packed_gpu": [
            "struct_active = atm_structure_enabled()",
            "size(turbineArray(i) % Cm)",
            "size(turbineArray(i) % pitchingMoment)",
            "call atm_pack_rank3(turbineArray(i) % Cm, packed_send, pos)",
            "call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cm)",
            "turbineArray(i) % pitchingMoment",
        ],
    }

    for name, needles in required.items():
        block = subroutine_block(text, name)
        for needle in needles:
            require(needle in block, f"{name} is missing `{needle}`", failures)


def check_batched_cl_correction(text: str, failures: list[str]) -> None:
    block = subroutine_block(text, "atm_lesgo_forcing")
    call_index = block.find("call atm_batch_cl_correction_gpu()")
    require(call_index >= 0, "atm_lesgo_forcing does not call batched Cl correction", failures)
    if call_index >= 0:
        prefix = block[max(0, call_index - 500) : call_index]
        require(
            ".not. atm_structure_enabled()" not in prefix,
            "batched Cl correction is still guarded by structure OFF",
            failures,
        )
    require(
        "if (turbineArray(i) % sampling /= 'atPoint') then" in block,
        "per-turbine Cl-correction fallback still appears structure-dependent",
        failures,
    )


def check_reset_payload(text: str, failures: list[str]) -> None:
    block = subroutine_block(text, "atm_lesgo_reset_turbine_gpu")
    require("Cm => turbineArray(i) % Cm" in block, "GPU reset does not point to Cm", failures)
    require(
        "pitchingMoment => turbineArray(i) % pitchingMoment" in block,
        "GPU reset does not point to pitchingMoment",
        failures,
    )
    require("Cm(m,n,q) = 0._rprec" in block, "GPU reset does not clear Cm", failures)
    require(
        "pitchingMoment(m,n,q) = 0._rprec" in block,
        "GPU reset does not clear pitchingMoment",
        failures,
    )
    require(
        "atm_reset_cuda_enabled() .and. .not. atm_structure_enabled()" not in text,
        "reset GPU path is still disabled by structure ON",
        failures,
    )


def check_remaining_structure_guards(text: str, failures: list[str]) -> None:
    allowed_contexts = [
        "Point-owner LB uses a reduced at-point force path.",
        "legacy CUDA at-point blade-force kernel does not cover",
    ]
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if ".not. atm_structure_enabled()" not in line:
            continue
        if "atm_model_cuda_enabled()" not in line:
            continue
        context = "\n".join(lines[max(0, index - 5) : index + 1])
        if not any(allowed in context for allowed in allowed_contexts):
            failures.append(
                "unexpected structure-off guard in atm_lesgo_interface.f90 "
                f"near line {index + 1}: {line.strip()}"
            )


def check_actuator_model_guards(text: str, failures: list[str]) -> None:
    cl_block = subroutine_block(text, "atm_compute_cl_correction")
    call_index = cl_block.find("call atm_compute_cl_correction_gpu(i)")
    require(
        call_index >= 0,
        "atm_compute_cl_correction no longer calls the CUDA correction",
        failures,
    )
    if call_index >= 0:
        prefix = cl_block[max(0, call_index - 250) : call_index]
        require(
            ".not. atm_structure_enabled()" not in prefix,
            "per-turbine CUDA Cl correction is still guarded by structure OFF",
            failures,
        )

    allowed_contexts = [
        "CUDA rotation shortcut only updates bladePoints",
        "CUDA yaw shortcut only rotates bladePoints",
    ]
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if ".not. atm_structure_enabled()" not in line:
            continue
        if "atm_model_cuda_enabled()" not in line:
            continue
        context = "\n".join(lines[max(0, index - 5) : index + 1])
        if not any(allowed in context for allowed in allowed_contexts):
            failures.append(
                "unexpected structure-off guard in actuator_turbine_model.f90 "
                f"near line {index + 1}: {line.strip()}"
            )


def main() -> int:
    interface_text = INTERFACE_PATH.read_text(encoding="utf-8")
    actuator_text = ACTUATOR_PATH.read_text(encoding="utf-8")
    failures: list[str] = []

    check_dispatch(interface_text, failures)
    check_packed_payload(interface_text, failures)
    check_batched_cl_correction(interface_text, failures)
    check_reset_payload(interface_text, failures)
    check_remaining_structure_guards(interface_text, failures)
    check_actuator_model_guards(actuator_text, failures)

    if failures:
        print("ATM structure packed-path check failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("ATM structure packed-path check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
