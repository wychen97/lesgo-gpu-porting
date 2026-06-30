#!/usr/bin/env python3
"""Verify the p0 public case scripts match the active branch profile."""

from __future__ import annotations

import sys
from pathlib import Path


EXEC_SUFFIX_ORDER = [
    ("USE_MPI", "-mpi"),
    ("USE_CPS", "-cps"),
    ("USE_HIT", "-HIT"),
    ("USE_LVLSET", "-ls"),
    ("USE_TURBINES", "-turbines"),
    ("USE_ATM", "-ATM"),
    ("USE_LES_GPU", "-lesgpu"),
    ("OUTPUT_EXTRA", "-exout"),
    ("USE_DYN_TN", "-dyntn"),
    ("USE_CGNS", "-cgns"),
    ("USE_SCALARS_GPU", "-scalgpu"),
    ("USE_SCALARS", "-scalars"),
]

CASES = {
    "channel_flow": {
        "cmake_flags": {
            "USE_MPI": True,
            "USE_TURBINES": False,
            "USE_ATM": False,
            "USE_CPS": False,
            "USE_SCALARS": False,
            "USE_SCALARS_GPU": False,
            "USE_LVLSET": False,
            "USE_HIT": False,
            "USE_DYN_TN": False,
            "USE_CGNS": False,
        },
        "job_cpu": "lesgo_channel_cpu",
        "job_gpu": "lesgo_channel_gpu",
    },
    "adm_disk": {
        "cmake_flags": {
            "USE_MPI": True,
            "USE_TURBINES": True,
            "USE_ATM": False,
            "USE_CPS": False,
            "USE_SCALARS": False,
            "USE_SCALARS_GPU": False,
            "USE_LVLSET": False,
            "USE_HIT": False,
            "USE_DYN_TN": False,
            "USE_CGNS": False,
        },
        "job_cpu": "lesgo_adm_cpu",
        "job_gpu": "lesgo_adm_gpu",
    },
    "atm_line": {
        "cmake_flags": {
            "USE_MPI": True,
            "USE_TURBINES": False,
            "USE_ATM": True,
            "USE_CPS": False,
            "USE_SCALARS": False,
            "USE_SCALARS_GPU": False,
            "USE_LVLSET": False,
            "USE_HIT": False,
            "USE_DYN_TN": False,
            "USE_CGNS": False,
        },
        "job_cpu": "lesgo_atm_cpu",
        "job_gpu": "lesgo_atm_gpu",
    },
}


def require(text: str, path: Path, marker: str, errors: list[str]) -> None:
    if marker not in text:
        errors.append(f"{path} is missing required marker: {marker}")


def cpu_default_branch(root: Path) -> bool:
    readme = root / "README.md"
    return readme.exists() and "LESGO CPU Baseline Branch" in readme.read_text(
        encoding="utf-8"
    )


def executable_name(config: dict[str, object], use_les_gpu: bool) -> str:
    flags = dict(config["cmake_flags"])
    flags["USE_LES_GPU"] = use_les_gpu
    name = "lesgo"
    for option, suffix in EXEC_SUFFIX_ORDER:
        if flags.get(option, False):
            name += suffix
    return name


def check_cmake_target_naming(root: Path) -> list[str]:
    path = root / "CMakeLists.txt"
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    require(text, path, 'set(exec_name "lesgo")', errors)

    previous = -1
    for option, suffix in EXEC_SUFFIX_ORDER:
        marker = f'set(exec_name "${{exec_name}}{suffix}")'
        index = text.find(marker)
        if index < 0:
            errors.append(f"{path} is missing target suffix marker for {option}: {marker}")
            continue
        if index <= previous:
            errors.append(
                f"{path} has target suffix marker for {option} out of expected order: {marker}"
            )
        previous = index
    return errors


def check_compile_script(
    case_dir: Path, config: dict[str, object], cpu_default: bool
) -> list[str]:
    path = case_dir / "compile_derecho.sh"
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    cpu_exe = executable_name(config, use_les_gpu=False)
    gpu_exe = executable_name(config, use_les_gpu=True)
    default_marker = (
        "BUILD_PROFILE=\"${1:-${BUILD_PROFILE:-cpu}}\""
        if cpu_default
        else "BUILD_PROFILE=\"${1:-${BUILD_PROFILE:-gpu}}\""
    )
    for marker in [
        default_marker,
        "case \"${BUILD_PROFILE}\" in",
        "USE_CPU_BUILD=OFF",
        "USE_CPU_BUILD=ON",
        "USE_LES_GPU=ON",
        "USE_LES_GPU=OFF",
        "USE_GPU_AWARE_MPI=AUTO",
        "USE_GPU_AWARE_MPI=OFF",
        "lesgo-run-exe-${BUILD_PROFILE}",
        "-DUSE_CPU_BUILD=\"${USE_CPU_BUILD}\"",
        "-DUSE_LES_GPU=\"${USE_LES_GPU}\"",
        "-DUSE_GPU_AWARE_MPI=\"${USE_GPU_AWARE_MPI}\"",
        f"EXE_NAME={cpu_exe}",
        f"EXE_NAME={gpu_exe}",
    ]:
        require(text, path, marker, errors)
    return errors


def check_submit_script(
    case_dir: Path, config: dict[str, object], cpu_default: bool
) -> list[str]:
    path = case_dir / "submit_derecho.pbs"
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    default_marker = (
        "RUN_PROFILE=\"${RUN_PROFILE:-cpu}\""
        if cpu_default
        else "RUN_PROFILE=\"${RUN_PROFILE:-gpu}\""
    )
    markers = [
        default_marker,
        "RUN_LABEL=\"${RUN_LABEL:-${RUN_PROFILE}}\"",
        "ARCHIVE_ROOT=\"${ARCHIVE_ROOT:-run-archives}\"",
        "case \"${RUN_PROFILE}\" in",
        "RUN_EXE=./lesgo-run-exe-gpu",
        "RUN_EXE=./lesgo-run-exe-cpu",
        "if [[ \"${RUN_PROFILE}\" == \"gpu\" ]]; then",
        "module load nvhpc/25.9 cuda/12.9.0 cray-mpich/8.1.32 fftw/3.3.10",
        "export MPICH_GPU_SUPPORT_ENABLED=1",
        "export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1",
        "module load nvhpc/25.9 cray-mpich/8.1.32 fftw/3.3.10",
        "unset MPICH_GPU_SUPPORT_ENABLED",
        "unset MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED",
        "set_gpu_rank \"${RUN_EXE}\"",
        "tee \"${LOG_FILE}\"",
        "ARCHIVE_DIR=\"${ARCHIVE_ROOT}/${RUN_LABEL}\"",
        "Archived run evidence to ${ARCHIVE_DIR}",
        "mpiexec -n \"${MPI_RANKS}\" -ppn \"${MPI_PPN}\" \"${RUN_EXE}\"",
        config["job_cpu"],
    ]
    if not cpu_default:
        markers.append(config["job_gpu"])
    for marker in markers:
        require(text, path, marker, errors)
    return errors


def check_readme(
    case_dir: Path, config: dict[str, object], cpu_default: bool
) -> list[str]:
    path = case_dir / "README.md"
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    markers = [
        "./compile_derecho.sh cpu",
        "RUN_PROFILE=cpu",
        "lesgo-run-exe-cpu",
        "RUN_LABEL",
        "run-archives",
        "does not set the Cray GPU-aware MPI environment variables",
        config["job_cpu"],
    ]
    if not cpu_default:
        markers.extend(
            [
                "./compile_derecho.sh gpu",
                "RUN_PROFILE=gpu",
                "lesgo-run-exe-gpu",
                config["job_gpu"],
            ]
        )
    for marker in markers:
        require(text, path, marker, errors)
    return errors


def main() -> int:
    repo_root = Path(".")
    root = repo_root / "test-cases"
    cpu_default = cpu_default_branch(repo_root)
    errors: list[str] = []
    errors.extend(check_cmake_target_naming(repo_root))
    for case_name, config in CASES.items():
        case_dir = root / case_name
        if not case_dir.is_dir():
            errors.append(f"missing p0 case directory: {case_dir}")
            continue
        errors.extend(check_compile_script(case_dir, config, cpu_default))
        errors.extend(check_submit_script(case_dir, config, cpu_default))
        errors.extend(check_readme(case_dir, config, cpu_default))

    if errors:
        print("p0 paired case script check failed:")
        for error in errors:
            print(f"  {error}")
        return 1

    mode = "CPU-default" if cpu_default else "GPU-default"
    print(f"p0 paired case script check passed ({len(CASES)} case trees, {mode}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
