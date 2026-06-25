#!/usr/bin/env python3
"""Prepare isolated compact SGS validation cases from a channel-flow case."""

from __future__ import annotations

import argparse
import re
import shlex
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SgsCase:
    selector: str
    sgs: str
    model: int
    use_dyn_tn: bool = False


SGS_BASE_CASES = [
    SgsCase("sgs_off", ".false.", 5),
    *[SgsCase(f"sgs_model_{model}", ".true.", model) for model in range(1, 6)],
]

SGS_DYN_TN_CASES = [
    SgsCase("sgs_model_4_dyn_tn", ".true.", 4, True),
    SgsCase("sgs_model_5_dyn_tn", ".true.", 5, True),
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
        if name.startswith("vel.out.c") or name.startswith("tavg.out.c"):
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
    for pattern in ["vel.out.c*", "tavg.out.c*", "lesgo_*.log", "*.o[0-9]*", "*.e[0-9]*"]:
        for path in case_dir.glob(pattern):
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def patch_generated_compile_script(case_dir: Path, sgs_case: SgsCase) -> None:
    script = case_dir / "compile_derecho.sh"
    if not script.exists():
        return
    text = script.read_text(encoding="utf-8")
    text = text.replace(
        'ROOT_DIR="$(cd "${CASE_DIR}/../.." && pwd -P)"',
        'ROOT_DIR="${LESGO_SOURCE_ROOT:-$(cd "${CASE_DIR}/../.." && pwd -P)}"',
    )
    if sgs_case.use_dyn_tn:
        text = text.replace("-DUSE_DYN_TN=OFF \\", "-DUSE_DYN_TN=ON \\")
        text = text.replace(
            "    EXE_NAME=lesgo-mpi-lesgpu\n",
            "    EXE_NAME=lesgo-mpi-lesgpu-dyntn\n",
        )
        text = text.replace(
            "    EXE_NAME=lesgo-mpi\n",
            "    EXE_NAME=lesgo-mpi-dyntn\n",
        )
    script.write_text(text, encoding="utf-8")


def prepare_case(
    *,
    base_case: Path,
    out_root: Path,
    sgs_case: SgsCase,
    profile: str,
    grid: int,
    nsteps: int,
    nenergy: int,
    cs_count: int,
    dyn_init: int,
    allow_missing_exe: bool,
    force: bool,
) -> Path:
    target = out_root / f"{sgs_case.selector}_{profile}"
    if target.exists():
        if not force:
            raise FileExistsError(f"{target} already exists; use --force to replace it")
        shutil.rmtree(target)
    shutil.copytree(base_case, target, ignore=ignore_generated)
    clean_case(target)
    patch_generated_compile_script(target, sgs_case)

    conf = target / "lesgo.conf"
    text = conf.read_text(encoding="utf-8")
    replacements = {
        "nproc": "1",
        "Nx": str(grid),
        "Ny": str(grid),
        "Nz": str(grid),
        "sgs_model": str(sgs_case.model),
        "cs_count": str(cs_count),
        "dyn_init": str(dyn_init),
        "sgs": sgs_case.sgs,
        "nsteps": str(nsteps),
        "nenergy": str(nenergy),
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
    conf.write_text(text, encoding="utf-8")

    exe = target / f"lesgo-run-exe-{profile}"
    if not exe.exists() and not (allow_missing_exe or sgs_case.use_dyn_tn):
        raise FileNotFoundError(f"{exe} is missing; compile the base case first")
    return target


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create compact SGS CPU/GPU validation case directories."
    )
    parser.add_argument("--base-case", type=Path, default=ROOT / "test-cases/channel_flow")
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--grid", type=int, default=128)
    parser.add_argument("--nsteps", type=int, default=40)
    parser.add_argument("--nenergy", type=int, default=4)
    parser.add_argument("--cs-count", type=int, default=5)
    parser.add_argument("--dyn-init", type=int, default=1)
    parser.add_argument(
        "--include-dyn-tn",
        action="store_true",
        help="also create USE_DYN_TN=ON sgs_model=4/5 cases; compile these generated cases before submitting",
    )
    parser.add_argument(
        "--allow-missing-exe",
        action="store_true",
        help="prepare cases without prebuilt executables; compile each generated case before submitting",
    )
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

    sgs_cases = list(SGS_BASE_CASES)
    if args.include_dyn_tn:
        sgs_cases.extend(SGS_DYN_TN_CASES)

    created: list[tuple[Path, SgsCase]] = []
    for sgs_case in sgs_cases:
        for profile in ["cpu", "gpu"]:
            path = prepare_case(
                    base_case=base_case,
                    out_root=out_root,
                    sgs_case=sgs_case,
                    profile=profile,
                    grid=args.grid,
                    nsteps=args.nsteps,
                    nenergy=args.nenergy,
                    cs_count=args.cs_count,
                    dyn_init=args.dyn_init,
                    allow_missing_exe=args.allow_missing_exe,
                    force=args.force,
                )
            created.append((path, sgs_case))

    print(f"Prepared {len(created)} SGS validation case directories under {out_root}")
    print()
    build_cases = [
        (path, case)
        for path, case in created
        if case.use_dyn_tn or args.allow_missing_exe
    ]
    if build_cases:
        print("# Build examples:")
        print("# Run these before submitting cases that do not already have executables.")
        print("# If the case is outside the repository tree, pass LESGO_SOURCE_ROOT.")
        for path, _case in build_cases:
            profile = "cpu" if path.name.endswith("_cpu") else "gpu"
            print(
                f"cd {shlex.quote(str(path))} && "
                f"LESGO_SOURCE_ROOT={shlex.quote(str(ROOT))} "
                f"./compile_derecho.sh {profile}"
            )
        print()
    print("# Submit examples:")
    for path, _case in created:
        profile = "cpu" if path.name.endswith("_cpu") else "gpu"
        label = profile
        queue = "-q develop " if profile == "cpu" else ""
        resources = (
            "-l select=1:ncpus=32:mpiprocs=1:mem=120gb -l place=shared "
            if profile == "cpu"
            else ""
        )
        print(
            f"cd {path} && qsub {queue}-N {path.name} "
            f"-v RUN_PROFILE={profile},RUN_LABEL={label} "
            f"{resources}submit_derecho.pbs"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
