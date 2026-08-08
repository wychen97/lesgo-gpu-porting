#!/usr/bin/env python3
"""Prepare the complete Level Set CPU/bridge/GPU validation matrix."""

from __future__ import annotations

import argparse
import array
import json
import math
import re
import shutil
import struct
from pathlib import Path


BUILD_PROFILES = {
    "cpu_mpi": {
        "USE_MPI": "ON", "USE_CPU_BUILD": "ON", "USE_LES_GPU": "OFF",
        "USE_LVLSET": "ON", "USE_LVLSET_GPU": "OFF", "USE_GPU_AWARE_MPI": "OFF",
    },
    "bridge_mpi": {
        "USE_MPI": "ON", "USE_CPU_BUILD": "OFF", "USE_LES_GPU": "ON",
        "USE_LVLSET": "ON", "USE_LVLSET_GPU": "OFF", "USE_GPU_AWARE_MPI": "OFF",
    },
    "gpu_mpi_staged": {
        "USE_MPI": "ON", "USE_CPU_BUILD": "OFF", "USE_LES_GPU": "ON",
        "USE_LVLSET": "ON", "USE_LVLSET_GPU": "ON", "USE_GPU_AWARE_MPI": "OFF",
    },
    "gpu_mpi_aware": {
        "USE_MPI": "ON", "USE_CPU_BUILD": "OFF", "USE_LES_GPU": "ON",
        "USE_LVLSET": "ON", "USE_LVLSET_GPU": "ON", "USE_GPU_AWARE_MPI": "ON",
    },
    "cpu_nompi": {
        "USE_MPI": "OFF", "USE_CPU_BUILD": "ON", "USE_LES_GPU": "OFF",
        "USE_LVLSET": "ON", "USE_LVLSET_GPU": "OFF", "USE_GPU_AWARE_MPI": "OFF",
    },
    "gpu_nompi": {
        "USE_MPI": "OFF", "USE_CPU_BUILD": "OFF", "USE_LES_GPU": "ON",
        "USE_LVLSET": "ON", "USE_LVLSET_GPU": "ON", "USE_GPU_AWARE_MPI": "OFF",
    },
}


def load_variants(path: Path) -> list[dict[str, object]]:
    """Load and validate the checked-in Level Set case definitions."""
    raw = json.loads(path.read_text(encoding="utf-8"))
    if raw.get("schema_version") != 1:
        raise ValueError(f"{path}: schema_version must be 1")
    base_settings = raw.get("base_settings")
    raw_variants = raw.get("variants")
    if not isinstance(base_settings, dict) or not base_settings:
        raise ValueError(f"{path}: base_settings must be a non-empty object")
    if not isinstance(raw_variants, list) or not raw_variants:
        raise ValueError(f"{path}: variants must be a non-empty array")
    if any(not isinstance(key, str) or not isinstance(value, str)
           for key, value in base_settings.items()):
        raise ValueError(f"{path}: base_settings keys and values must be strings")
    if base_settings.get("dyn_init") != "1" or base_settings.get("cs_count") != "2":
        raise ValueError(f"{path}: dyn_init=1 and cs_count=2 are required")

    variants: list[dict[str, object]] = []
    names: set[str] = set()
    for index, item in enumerate(raw_variants):
        context = f"{path}: variants[{index}]"
        if not isinstance(item, dict):
            raise ValueError(f"{context} must be an object")
        name = item.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError(f"{context}.name must be a non-empty string")
        if name in names:
            raise ValueError(f"{context}.name duplicates {name!r}")
        names.add(name)

        nproc = item.get("nproc")
        if not isinstance(nproc, int) or isinstance(nproc, bool) or nproc < 1:
            raise ValueError(f"{context}.nproc must be a positive integer")
        geometry = item.get("geometry")
        if geometry not in {"trees", "sphere", "tilted"}:
            raise ValueError(f"{context}.geometry is not supported")
        profiles = item.get("profiles")
        if not isinstance(profiles, list) or not profiles:
            raise ValueError(f"{context}.profiles must be a non-empty array")
        if any(not isinstance(profile, str) or not profile for profile in profiles):
            raise ValueError(f"{context}.profiles entries must be non-empty strings")
        if len(set(profiles)) != len(profiles):
            raise ValueError(f"{context}.profiles contains duplicates")
        unknown_profiles = [profile for profile in profiles if profile not in BUILD_PROFILES]
        if unknown_profiles:
            raise ValueError(f"{context}.profiles contains unknown entries: {unknown_profiles}")
        if nproc > 1 and any(BUILD_PROFILES[profile]["USE_MPI"] != "ON"
                             for profile in profiles):
            raise ValueError(f"{context}: multi-rank variants require MPI profiles")

        overrides = item.get("settings", {})
        if not isinstance(overrides, dict) or any(
            not isinstance(key, str) or not isinstance(value, str)
            for key, value in overrides.items()
        ):
            raise ValueError(f"{context}.settings keys and values must be strings")
        settings = dict(base_settings)
        settings.update(overrides)
        settings["nproc"] = str(nproc)
        expected_error = item.get("expected_error")
        if expected_error is not None and not isinstance(expected_error, str):
            raise ValueError(f"{context}.expected_error must be a string or null")
        variants.append({
            "name": name,
            "nproc": nproc,
            "geometry": geometry,
            "profiles": list(profiles),
            "settings": settings,
            "expected_error": expected_error,
        })
    return variants


def replace_setting(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?im)^(\s*{re.escape(key)}\s*=\s*).*$")
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise ValueError(f"lesgo.conf does not contain exactly one {key}= entry")
    return updated


def read_setting(text: str, key: str) -> float:
    match = re.search(rf"(?im)^\s*{re.escape(key)}\s*=\s*([^!#\s]+)", text)
    if not match:
        raise ValueError(f"missing {key} in lesgo.conf")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def write_record(path: Path, values: array.array) -> None:
    payload = values.tobytes()
    marker = struct.pack("<i", len(payload))
    path.write_bytes(marker + payload + marker)


def generate_geometry(
    case_dir: Path, conf: str, kind: str, nproc: int, *, use_mpi: bool,
) -> None:
    nx = int(read_setting(conf, "Nx"))
    ny = int(read_setting(conf, "Ny"))
    nz_total_cells = int(read_setting(conf, "Nz"))
    if nz_total_cells % nproc:
        raise ValueError("Nz must be divisible by nproc")
    local_nz = nz_total_cells // nproc + 1
    nz_global = (local_nz - 1) * nproc + 1
    ld = 2 * (nx // 2 + 1)
    lx = read_setting(conf, "Lx")
    dx = lx / nx
    dy = dx
    dz = dx
    ly = ny * dy
    lz = (nz_global - 1) * dz
    center = (0.5 * lx, 0.5 * ly, 0.5 * lz)
    for rank in range(nproc):
        values = array.array("d")
        for k in range(0, local_nz + 1):
            z = (rank * (local_nz - 1) + k - 0.5) * dz
            for j in range(ny):
                y = j * dy
                for i in range(ld):
                    x = i * dx
                    if kind == "sphere":
                        radius = 0.22 * min(lx, ly, lz)
                        phi = math.sqrt(
                            (x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2
                        ) - radius
                    elif kind == "tilted":
                        angle = math.radians(30.0)
                        xr = math.cos(angle) * (x - center[0]) + math.sin(angle) * (z - center[2])
                        zr = -math.sin(angle) * (x - center[0]) + math.cos(angle) * (z - center[2])
                        yr = y - center[1]
                        ax, ay, az = 0.20 * lx, 0.16 * ly, 0.28 * lz
                        phi = (math.sqrt((xr / ax) ** 2 + (yr / ay) ** 2 + (zr / az) ** 2) - 1.0) * min(ax, ay, az)
                    else:
                        raise ValueError(f"unknown geometry {kind}")
                    values.append(phi)
        name = f"phi.out.c{rank}" if use_mpi else "phi.out"
        write_record(case_dir / name, values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-case", type=Path, default=Path("test-cases/level_set_cubes"))
    parser.add_argument(
        "--variants", type=Path,
        help="checked-in variant JSON (default: BASE_CASE/validation_variants.json)",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    base = args.base_case.resolve()
    variants_path = (
        args.variants.resolve() if args.variants is not None
        else base / "validation_variants.json"
    )
    variants = load_variants(variants_path)
    output = args.out.resolve()
    if output.exists():
        if not args.force:
            raise SystemExit(f"{output} already exists; pass --force to replace it")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    base_conf = (base / "lesgo.conf").read_text(encoding="utf-8")
    tasks: list[dict[str, object]] = []
    for item in variants:
        for profile in item["profiles"]:
            task_id = f"{item['name']}__{profile}"
            case_dir = output / "cases" / task_id
            case_dir.mkdir(parents=True)
            conf = base_conf
            for key, value in dict(item["settings"]).items():
                conf = replace_setting(conf, key, value)
            if item["geometry"] == "trees":
                shutil.copy2(base / "trees.conf", case_dir / "trees.conf")
            else:
                conf = replace_setting(conf, "use_trees", ".false.")
                generate_geometry(
                    case_dir, conf, str(item["geometry"]), int(item["nproc"]),
                    use_mpi=BUILD_PROFILES[profile]["USE_MPI"] == "ON",
                )
            (case_dir / "lesgo.conf").write_text(conf, encoding="utf-8")
            tasks.append(
                {
                    "id": task_id,
                    "variant": item["name"],
                    "profile": profile,
                    "nproc": item["nproc"],
                    "geometry": item["geometry"],
                    "case_dir": str(case_dir.relative_to(output)),
                    "expected_error": item["expected_error"],
                    "environment": {
                        "LESGO_LVLSET_VALIDATION_SNAPSHOT": "ON",
                        "LESGO_LVLSET_INTERP_BOUNDS_CHECK": "ON",
                        "LESGO_RANDOM_SEED": "20260807",
                    },
                }
            )
    manifest = {
        "schema_version": 1,
        "purpose": "Level Set CPU, host-bridge, full-GPU, MPI, and non-MPI validation",
        "variant_source": str(variants_path),
        "required_runtime_settings": {"dyn_init": 1, "cs_count": 2},
        "build_profiles": BUILD_PROFILES,
        "tasks": tasks,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "tasks": len(tasks)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
