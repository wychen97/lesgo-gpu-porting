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


BASE_SETTINGS = {
    "nsteps": "20",
    "dyn_init": "1",
    "cs_count": "2",
    "vel_BC": ".false.",
    "use_log_profile": ".false.",
    "use_enforce_un": ".false.",
    "use_extrap_tau_log": ".false.",
    "use_extrap_tau_simple": ".true.",
    "use_modify_dutdn": ".false.",
    "lag_dyn_modify_beta": ".true.",
    "smooth_mode": "'xy'",
    "lbc_mom": "0",
}


def variant(name: str, *, nproc: int = 1, geometry: str = "trees", **settings: str) -> dict[str, object]:
    merged = dict(BASE_SETTINGS)
    merged.update(settings)
    merged["nproc"] = str(nproc)
    return {"name": name, "nproc": nproc, "geometry": geometry, "settings": merged}


VARIANTS = [
    variant("sgs_off", sgs=".false.", sgs_model="1"),
    *[
        variant(f"sgs_model_{model}", sgs=".true.", sgs_model=str(model))
        for model in range(1, 6)
    ],
    variant("model4_beta_off_lbc0", sgs=".true.", sgs_model="4", lag_dyn_modify_beta=".false."),
    variant("model4_beta_on_lbc1", sgs=".true.", sgs_model="4", lbc_mom="1"),
    variant("extrap_tau_log", sgs=".true.", sgs_model="2", use_extrap_tau_log=".true."),
    variant("legacy_extrapolation", sgs=".true.", sgs_model="2", use_extrap_tau_simple=".false."),
    variant("modify_dutdn", sgs=".true.", sgs_model="2", use_modify_dutdn=".true."),
    variant(
        "desired_log_velocity", sgs=".true.", sgs_model="2", vel_BC=".true.",
        use_log_profile=".true.",
    ),
    variant("smooth_3d", sgs=".true.", sgs_model="2", smooth_mode="'3d'"),
    variant("sphere_rank_crossing", nproc=2, geometry="sphere", sgs=".true.", sgs_model="4"),
    variant("tilted_rank_crossing", nproc=2, geometry="tilted", sgs=".true.", sgs_model="5", lbc_mom="1"),
]


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


def generate_geometry(case_dir: Path, conf: str, kind: str, nproc: int) -> None:
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
        name = "phi.out" if nproc == 1 else f"phi.out.c{rank}"
        write_record(case_dir / name, values)


def profiles_for_variant(name: str, nproc: int) -> list[str]:
    if name == "sgs_model_4":
        return ["cpu_mpi", "bridge_mpi", "gpu_mpi_staged"]
    if name == "sgs_model_5":
        return ["cpu_mpi", "gpu_mpi_staged", "cpu_nompi", "gpu_nompi"]
    if nproc > 1:
        return ["cpu_mpi", "gpu_mpi_staged", "gpu_mpi_aware"]
    return ["cpu_mpi", "gpu_mpi_staged"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-case", type=Path, default=Path("test-cases/level_set_cubes"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    base = args.base_case.resolve()
    output = args.out.resolve()
    if output.exists():
        if not args.force:
            raise SystemExit(f"{output} already exists; pass --force to replace it")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    base_conf = (base / "lesgo.conf").read_text(encoding="utf-8")
    tasks: list[dict[str, object]] = []
    for item in VARIANTS:
        for profile in profiles_for_variant(str(item["name"]), int(item["nproc"])):
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
                generate_geometry(case_dir, conf, str(item["geometry"]), int(item["nproc"]))
            (case_dir / "lesgo.conf").write_text(conf, encoding="utf-8", newline="\n")
            tasks.append(
                {
                    "id": task_id,
                    "variant": item["name"],
                    "profile": profile,
                    "nproc": item["nproc"],
                    "geometry": item["geometry"],
                    "case_dir": str(case_dir.relative_to(output)),
                    "environment": {"LESGO_LVLSET_VALIDATION_SNAPSHOT": "ON"},
                }
            )
    manifest = {
        "schema_version": 1,
        "purpose": "Level Set CPU, host-bridge, full-GPU, MPI, and non-MPI validation",
        "required_runtime_settings": {"dyn_init": 1, "cs_count": 2},
        "build_profiles": BUILD_PROFILES,
        "tasks": tasks,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "tasks": len(tasks)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
