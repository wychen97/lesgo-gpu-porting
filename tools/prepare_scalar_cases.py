#!/usr/bin/env python3
"""Prepare isolated compact scalar validation cases from channel flow."""

from __future__ import annotations

import argparse
import re
import shlex
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from repo_paths import repo_path


@dataclass(frozen=True)
class ScalarCase:
    selector: str
    passive_scalar: bool


SCALAR_CASES = [
    ScalarCase("scalar_passive_grid_128", True),
    ScalarCase("scalar_active_grid_128", False),
]


GENERATED_NAMES = {
    "build-derecho-cpu",
    "build-derecho-gpu",
    "output",
    "turbineOutput",
    "run-archives",
}


def replace_key(text: str, key: str, value: str) -> str:
    pattern = re.compile(
        rf"^(?P<indent>\s*){re.escape(key)}\s*=\s*(?P<old>[^!\n]*)(?P<tail>.*)$",
        flags=re.IGNORECASE | re.MULTILINE,
    )

    def repl(match: re.Match[str]) -> str:
        return f"{match.group('indent')}{key} = {value}{match.group('tail')}"

    updated, count = pattern.subn(repl, text, count=1)
    if count != 1:
        raise ValueError(f"could not replace `{key}` in lesgo.conf")
    return updated


def ignore_generated(_directory: str, names: list[str]) -> set[str]:
    ignored = set()
    for name in names:
        if name in GENERATED_NAMES:
            ignored.add(name)
        if (
            name.startswith("lesgo_")
            and (".log" in name or ".o" in name or ".e" in name)
        ):
            ignored.add(name)
        if name in {"grid.out", "total_time.dat", "check_ke.out", "lesgo_param.out"}:
            ignored.add(name)
        if (
            name.startswith("vel.out.c")
            or name.startswith("tavg.out.c")
            or name.startswith("scal.out.c")
        ):
            ignored.add(name)
    return ignored


def clean_case(case_dir: Path) -> None:
    for name in [
        "output",
        "turbineOutput",
        "run-archives",
        "grid.out",
        "total_time.dat",
        "check_ke.out",
        "lesgo_param.out",
    ]:
        path = case_dir / name
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    for pattern in [
        "vel.out.c*",
        "tavg.out.c*",
        "scal.out.c*",
        "lesgo_*.log",
        "*.o[0-9]*",
        "*.e[0-9]*",
    ]:
        for path in case_dir.glob(pattern):
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def append_or_replace_scalar_block(text: str, scalar_case: ScalarCase) -> str:
    block = f"""

SCALARS {{

  lbc_scal = 0
  scal_bot = 300.0
  flux_bot = 0.0
  read_lbc_scal = .false.
  lapse_rate = 0.0
  ic_z = 0.0, 630.0
  ic_theta = 300.0, 301.0
  ic_no_vel_noise_z = 0.0
  g = 9.81
  zo_s = 0.00001
  T_scale = 300.0
  passive_scalar = {'.true.' if scalar_case.passive_scalar else '.false.'}
  Pr_sgs = 0.4

}}
"""
    pattern = re.compile(r"^SCALARS\s*\{.*?^\}", re.IGNORECASE | re.MULTILINE | re.DOTALL)
    if pattern.search(text):
        return pattern.sub(block.strip(), text, count=1)
    if not text.endswith("\n"):
        text += "\n"
    return text + block


def patch_generated_compile_script(case_dir: Path) -> None:
    script = case_dir / "compile_derecho.sh"
    if not script.exists():
        return
    text = script.read_text(encoding="utf-8")
    text = text.replace(
        'ROOT_DIR="$(cd "${CASE_DIR}/../.." && pwd -P)"',
        'ROOT_DIR="${LESGO_SOURCE_ROOT:-$(cd "${CASE_DIR}/../.." && pwd -P)}"',
    )
    text = text.replace(
        "    USE_GPU_AWARE_MPI=AUTO\n    EXE_NAME=lesgo-mpi-lesgpu\n",
        "    USE_GPU_AWARE_MPI=AUTO\n    USE_SCALARS_GPU=ON\n    EXE_NAME=lesgo-mpi-lesgpu-scalgpu-scalars\n",
    )
    text = text.replace(
        "    USE_GPU_AWARE_MPI=OFF\n    EXE_NAME=lesgo-mpi\n",
        "    USE_GPU_AWARE_MPI=OFF\n    USE_SCALARS_GPU=OFF\n    EXE_NAME=lesgo-mpi-scalars\n",
    )
    text = text.replace("-DUSE_SCALARS=OFF \\", "-DUSE_SCALARS=ON \\")
    text = text.replace(
        '-DUSE_SCALARS_GPU=OFF \\',
        '-DUSE_SCALARS_GPU="${USE_SCALARS_GPU}" \\',
    )
    script.write_text(text, encoding="utf-8")


def patch_generated_submit_script(case_dir: Path) -> None:
    script = case_dir / "submit_derecho.pbs"
    if not script.exists():
        return
    text = script.read_text(encoding="utf-8")
    marker = 'ARCHIVE_ROOT="${ARCHIVE_ROOT:-run-archives}"\n'
    replacement = marker + 'export LESGO_SCALAR_STAGE_TIMING="${LESGO_SCALAR_STAGE_TIMING:-1}"\n'
    if marker in text and "LESGO_SCALAR_STAGE_TIMING" not in text:
        text = text.replace(marker, replacement)
    text = text.replace(
        "rm -f total_time.dat check_ke.out grid.out vel.out.c* tavg.out.c* lesgo_param.out",
        "rm -f total_time.dat check_ke.out grid.out vel.out.c* tavg.out.c* scal.out.c* lesgo_param.out",
    )
    text = text.replace(
        "for item in vel.out.c* tavg.out.c*; do",
        "for item in vel.out.c* tavg.out.c* scal.out.c*; do",
    )
    script.write_text(text, encoding="utf-8")


def prepare_case(
    *,
    base_case: Path,
    out_root: Path,
    scalar_case: ScalarCase,
    profile: str,
    grid: int,
    nsteps: int,
    nenergy: int,
    force: bool,
) -> Path:
    selector = re.sub(r"grid_\d+", f"grid_{grid}", scalar_case.selector)
    target = out_root / f"{selector}_{profile}"
    if target.exists():
        if not force:
            raise FileExistsError(f"{target} already exists; use --force to replace it")
        shutil.rmtree(target)
    shutil.copytree(base_case, target, ignore=ignore_generated)
    clean_case(target)
    patch_generated_compile_script(target)
    patch_generated_submit_script(target)

    conf = target / "lesgo.conf"
    text = conf.read_text(encoding="utf-8")
    replacements = {
        "nproc": "1",
        "Nx": str(grid),
        "Ny": str(grid),
        "Nz": str(grid),
        "nsteps": str(nsteps),
        "nenergy": str(nenergy),
        "wbase": str(nenergy),
        "checkpoint_data": ".false.",
        "tavg_calc": ".false.",
        "point_calc": ".false.",
        "domain_calc": ".false.",
        "xplane_calc": ".false.",
        "yplane_calc": ".false.",
        "zplane_calc": ".false.",
        "sgs_hist_calc": ".false.",
    }
    for key, value in replacements.items():
        text = replace_key(text, key, value)
    text = append_or_replace_scalar_block(text, scalar_case)
    conf.write_text(text, encoding="utf-8")
    return target


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create compact scalar CPU/GPU validation case directories."
    )
    parser.add_argument(
        "--base-case",
        type=Path,
        default=repo_path("test-cases", "channel_flow"),
    )
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument(
        "--case",
        choices=[case.selector for case in SCALAR_CASES],
        action="append",
        help="case selector to create; defaults to all scalar selectors",
    )
    parser.add_argument("--grid", type=int, default=128)
    parser.add_argument("--nsteps", type=int, default=40)
    parser.add_argument("--nenergy", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    base_case = args.base_case.resolve()
    out_root = args.out_root.resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    if not (base_case / "lesgo.conf").exists():
        print(f"base case is missing lesgo.conf: {base_case}", file=sys.stderr)
        return 1

    selected = set(args.case or [case.selector for case in SCALAR_CASES])
    created: list[Path] = []
    for scalar_case in SCALAR_CASES:
        if scalar_case.selector not in selected:
            continue
        for profile in ["cpu", "gpu"]:
            created.append(
                prepare_case(
                    base_case=base_case,
                    out_root=out_root,
                    scalar_case=scalar_case,
                    profile=profile,
                    grid=args.grid,
                    nsteps=args.nsteps,
                    nenergy=args.nenergy,
                    force=args.force,
                )
            )

    print(f"Prepared {len(created)} scalar validation case directories under {out_root}")
    print()
    print("# Build examples:")
    for path in created:
        profile = "cpu" if path.name.endswith("_cpu") else "gpu"
        print(
            f"cd {shlex.quote(str(path))} && "
            f"LESGO_SOURCE_ROOT={shlex.quote(str(ROOT))} "
            f"./compile_derecho.sh {profile}"
        )
    print()
    print("# Submit examples:")
    for path in created:
        profile = "cpu" if path.name.endswith("_cpu") else "gpu"
        label = profile
        if profile == "cpu":
            print(
                f"cd {path} && qsub -q develop -N {path.name} "
                f"-v RUN_PROFILE=cpu,RUN_LABEL={label} "
                "-l select=1:ncpus=32:mpiprocs=1:mem=120gb -l place=shared "
                "submit_derecho.pbs"
            )
        else:
            print(
                f"cd {path} && qsub -N {path.name} "
                f"-v RUN_PROFILE=gpu,RUN_LABEL={label} submit_derecho.pbs"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
