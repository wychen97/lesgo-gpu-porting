#!/usr/bin/env python3
"""Verify the GPU benchmark manifest matches the validation matrix."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from cmake_metadata import public_knob_values


MANIFEST_PATH = Path("docs/gpu_benchmark_manifest.json")
MATRIX_PATH = Path("docs/gpu_validation_matrix.md")
ROW_RE = re.compile(r"^\|\s*`(?P<id>[a-z0-9_]+)`\s*\|")
CMAKE_SETTING_RE = re.compile(r"^(?P<name>[A-Za-z0-9_]+)=(?P<value>[A-Za-z0-9_./+-]+)$")

REQUIRED_CASE_FIELDS = {
    "id",
    "status",
    "public_case",
    "base_grid",
    "gpu_cmake",
    "cpu_cmake",
    "runtime",
    "required_evidence",
}

ALLOWED_STATUSES = {
    "needs_paired_cpu_gpu",
    "needs_current_cpu_gpu_reference",
    "needs_new_case",
    "needs_new_case_matrix",
    "needs_public_evidence_copy_or_rerun",
    "host_boundary_compatibility",
    "excluded",
}


def parse_cmake_settings(case_id: str, field: str, settings: object) -> tuple[dict[str, str], list[str]]:
    errors: list[str] = []
    parsed: dict[str, str] = {}
    valid_values = public_knob_values()
    if not isinstance(settings, list):
        return parsed, [f"{MANIFEST_PATH} case `{case_id}` field `{field}` must be a list"]

    for raw_setting in settings:
        if not isinstance(raw_setting, str):
            errors.append(
                f"{MANIFEST_PATH} case `{case_id}` field `{field}` has a non-string setting"
            )
            continue
        match = CMAKE_SETTING_RE.match(raw_setting)
        if not match:
            errors.append(
                f"{MANIFEST_PATH} case `{case_id}` field `{field}` has invalid "
                f"CMake setting syntax: {raw_setting}"
            )
            continue
        name = match.group("name")
        value = match.group("value")
        if name not in valid_values:
            errors.append(
                f"{MANIFEST_PATH} case `{case_id}` field `{field}` references "
                f"unknown root CMake knob `{name}`"
            )
            continue
        allowed_values = valid_values[name]
        if allowed_values and value not in allowed_values:
            errors.append(
                f"{MANIFEST_PATH} case `{case_id}` field `{field}` sets `{name}` "
                f"to `{value}`, expected one of {sorted(allowed_values)}"
            )
        if name in parsed and parsed[name] != value:
            errors.append(
                f"{MANIFEST_PATH} case `{case_id}` field `{field}` sets `{name}` "
                f"to both `{parsed[name]}` and `{value}`"
            )
        parsed[name] = value
    return parsed, errors


def setting_is_on(settings: dict[str, str], name: str) -> bool:
    return settings.get(name) == "ON"


def check_cmake_profile_invariants(
    case_id: str,
    status: str,
    gpu_settings: dict[str, str],
    cpu_settings: dict[str, str],
) -> list[str]:
    if status == "excluded":
        if setting_is_on(gpu_settings, "USE_LVLSET"):
            return []
        return [
            f"{MANIFEST_PATH} excluded case `{case_id}` should make its "
            "LVLSET scope explicit with USE_LVLSET=ON"
        ]

    errors: list[str] = []
    if not setting_is_on(gpu_settings, "USE_LES_GPU"):
        errors.append(
            f"{MANIFEST_PATH} case `{case_id}` gpu_cmake must include USE_LES_GPU=ON"
        )
    if not setting_is_on(cpu_settings, "USE_CPU_BUILD"):
        errors.append(
            f"{MANIFEST_PATH} case `{case_id}` cpu_cmake must include USE_CPU_BUILD=ON"
        )
    if cpu_settings.get("USE_LES_GPU") != "OFF":
        errors.append(
            f"{MANIFEST_PATH} case `{case_id}` cpu_cmake must include USE_LES_GPU=OFF"
        )
    if setting_is_on(gpu_settings, "USE_LVLSET") or setting_is_on(cpu_settings, "USE_LVLSET"):
        errors.append(
            f"{MANIFEST_PATH} case `{case_id}` is in the non-LVLSET validation "
            "scope but enables USE_LVLSET"
        )
    for field, settings in [("gpu_cmake", gpu_settings), ("cpu_cmake", cpu_settings)]:
        if setting_is_on(settings, "USE_SCALARS_GPU"):
            if not setting_is_on(settings, "USE_SCALARS"):
                errors.append(
                    f"{MANIFEST_PATH} case `{case_id}` {field} enables "
                    "USE_SCALARS_GPU without USE_SCALARS=ON"
                )
            if not setting_is_on(settings, "USE_LES_GPU"):
                errors.append(
                    f"{MANIFEST_PATH} case `{case_id}` {field} enables "
                    "USE_SCALARS_GPU without USE_LES_GPU=ON"
                )
        if setting_is_on(settings, "USE_ATM") and setting_is_on(settings, "USE_TURBINES"):
            errors.append(
                f"{MANIFEST_PATH} case `{case_id}` {field} enables both "
                "USE_ATM and USE_TURBINES; split ATM and ADM validation rows"
            )
    return errors


def matrix_ids() -> set[str]:
    ids: set[str] = set()
    in_matrix = False
    for line in MATRIX_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip() == "## Validation Matrix":
            in_matrix = True
            continue
        if in_matrix and line.startswith("## "):
            break
        if not in_matrix:
            continue
        match = ROW_RE.match(line)
        if match:
            ids.add(match.group("id"))
    return ids


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def main() -> int:
    ok = True
    manifest = load_manifest()
    if manifest.get("schema_version") != 1:
        print(f"{MANIFEST_PATH} must have schema_version 1")
        ok = False

    cases = manifest.get("cases")
    if not isinstance(cases, list):
        print(f"{MANIFEST_PATH} must contain a cases list")
        return 1

    matrix = matrix_ids()
    case_ids = [case.get("id") for case in cases if isinstance(case, dict)]
    manifest_ids = set(case_ids)

    if len(case_ids) != len(manifest_ids):
        print(f"{MANIFEST_PATH} has duplicate case ids")
        ok = False

    missing = sorted(matrix - manifest_ids)
    stale = sorted(manifest_ids - matrix)
    if missing:
        print(f"{MANIFEST_PATH} is missing matrix ids:")
        for row_id in missing:
            print(f"  {row_id}")
        ok = False
    if stale:
        print(f"{MANIFEST_PATH} has ids not present in {MATRIX_PATH}:")
        for row_id in stale:
            print(f"  {row_id}")
        ok = False

    for case in cases:
        if not isinstance(case, dict):
            print(f"{MANIFEST_PATH} has a non-object case entry")
            ok = False
            continue
        case_id = case.get("id", "<missing id>")
        missing_fields = sorted(REQUIRED_CASE_FIELDS - set(case))
        if missing_fields:
            print(f"{MANIFEST_PATH} case `{case_id}` is missing fields:")
            for field in missing_fields:
                print(f"  {field}")
            ok = False
        status = case.get("status")
        if status not in ALLOWED_STATUSES:
            print(f"{MANIFEST_PATH} case `{case_id}` has invalid status: {status}")
            ok = False
        for field in ["runtime", "required_evidence"]:
            value = case.get(field)
            if not isinstance(value, list):
                print(f"{MANIFEST_PATH} case `{case_id}` field `{field}` must be a list")
                ok = False
        gpu_settings, gpu_errors = parse_cmake_settings(
            case_id, "gpu_cmake", case.get("gpu_cmake")
        )
        cpu_settings, cpu_errors = parse_cmake_settings(
            case_id, "cpu_cmake", case.get("cpu_cmake")
        )
        for error in gpu_errors + cpu_errors:
            print(error)
            ok = False
        for error in check_cmake_profile_invariants(
            case_id, str(status), gpu_settings, cpu_settings
        ):
            print(error)
            ok = False
        public_case = case.get("public_case")
        if public_case is not None and not Path(public_case).exists():
            print(f"{MANIFEST_PATH} case `{case_id}` public_case does not exist: {public_case}")
            ok = False

    if not ok:
        return 1

    print(
        "GPU benchmark manifest check passed "
        f"({len(cases)} cases, {len(ALLOWED_STATUSES)} statuses)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
