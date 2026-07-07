#!/usr/bin/env python3
"""Report a static non-LVLSET GPU-port inventory by Fortran subprogram."""

from __future__ import annotations

import argparse
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORTRAN_RE = re.compile(r".*\.(?:f|f90)$", re.IGNORECASE)
START_RE = re.compile(
    r"^\s*(?:(?:recursive|pure|elemental|module)\s+)*"
    r"(?P<kind>subroutine|function)\s+(?P<name>[a-z_][a-z0-9_]*)\b",
    re.IGNORECASE,
)
END_RE = re.compile(r"^\s*end\s+(?:subroutine|function)\b", re.IGNORECASE)

GPU_MARKERS = [
    "!$acc",
    "!$cuf",
    "attributes(",
    "host_data use_device",
    "PPLES_GPU",
    "PPCONVEC_GPU",
    "PPDERIVS_GPU",
    "PPPRESS_GPU",
    "PPSGS_GPU",
    "PPSCALARS_GPU",
    "PPGPU_AWARE_MPI",
    "cuda",
    "cufft",
]

HOST_BOUNDARY_FILES = {
    "clocks.f90",
    "cuda_mpi_debug.f90",
    "finalize.f90",
    "grid.f90",
    "init_random_seed.f90",
    "initial.f90",
    "initialize.f90",
    "input_util.f90",
    "io.f90",
    "messages.f90",
    "mpi_defs.f90",
    "param_output.f90",
    "pid.f90",
    "stat_defs.f90",
    "string_util.f90",
    "types.f90",
}

HOST_OR_DIAGNOSTIC_NAME_RE = re.compile(
    r"(^|_)("
    r"apply_env|checkpoint|count|debug|destroy|diag|ensure|env|finalize|"
    r"init|initialize|output|print|read|report|restart|stage_report|timer|"
    r"token|validate|write"
    r")($|_)",
    re.IGNORECASE,
)

def tracked_gpu_source_names() -> set[str]:
    """Return tracked root-level GPU helper source filenames."""
    result = subprocess.run(
        ["git", "ls-files", "*_gpu.f90", "*_gpu.F90"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    names: set[str] = set()
    for name in result.stdout.splitlines():
        path = Path(name)
        lower_name = path.name.lower()
        if path.parent != Path("."):
            continue
        if lower_name.startswith("level_set") or lower_name.startswith("trees_"):
            continue
        names.add(path.name)
    return names


GPU_SOURCE_FILES = tracked_gpu_source_names()

REVIEW_BUCKET_DESCRIPTIONS = {
    "atm-host-model": (
        "ATM blade, controller, structure, and small math helpers that remain "
        "host-side in the current hybrid design."
    ),
    "atm-mirror-lb-control": (
        "ATM mirror, synchronization, cell-search, and load-balance control "
        "helpers around the GPU sampling/forcing path."
    ),
    "adm-cpu-fallback-profile": (
        "ADM/turbine CPU fallback or compatibility routines; profile before "
        "treating them as missing GPU work."
    ),
    "generic-helper-profile": (
        "Generic interpolation/math helpers that may be CPU fallback or "
        "low-cost support code."
    ),
    "inflow-fringe-profile": (
        "Inflow/fringe helpers that need targeted runtime validation for "
        "nonstandard inflow configurations."
    ),
    "scalar-init-fallback": (
        "Scalar initialization, stability helper, or CPU fallback routines; "
        "validate passive and active scalar cases separately."
    ),
    "iwm-wallmodel-profile": (
        "IWM wall-model candidate that needs an IWM-heavy correctness and "
        "timing case before broad speed claims."
    ),
    "cpu-fallback-compat": (
        "CPU fallback or host compatibility routines retained beside GPU "
        "production paths."
    ),
    "diagnostic-profiling": (
        "Profiling, timing, or audit helpers; not GPU hot-path kernels."
    ),
    "excluded-lvlset-bridge": (
        "LVLSET bridge code excluded from the current non-LVLSET scope."
    ),
    "needs-manual-review": "Candidate not matched by the review buckets.",
}

BUCKET_VALIDATION_ROWS = {
    "adm-cpu-fallback-profile": ["adm_disk", "adm_dynamic_controls"],
    "atm-host-model": ["atm_line", "large_windfarm"],
    "atm-mirror-lb-control": ["atm_line", "large_windfarm"],
    "cpu-fallback-compat": ["les_core_channel", "hit_inflow"],
    "diagnostic-profiling": ["diagnostics_output"],
    "excluded-lvlset-bridge": ["lvlset"],
    "generic-helper-profile": ["les_core_channel", "adm_disk", "atm_line"],
    "inflow-fringe-profile": ["hit_inflow", "shifted_inflow", "sponge_coriolis"],
    "iwm-wallmodel-profile": ["iwm_wall_model"],
    "scalar-init-fallback": ["scalar_passive", "scalar_active", "cps_scalar"],
}


@dataclass(frozen=True)
class Subprogram:
    path: Path
    name: str
    start_line: int
    end_line: int
    classification: str


def tracked_fortran() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    paths = []
    for name in result.stdout.splitlines():
        path = Path(name)
        if not FORTRAN_RE.match(name):
            continue
        if path.parent != Path("."):
            continue
        lower_name = path.name.lower()
        if lower_name.startswith("level_set") or lower_name.startswith("trees_"):
            continue
        paths.append(path)
    return sorted(paths)


def classify(path: Path, name: str, body: str) -> str:
    lower_body = body.lower()
    has_gpu_marker = any(marker.lower() in lower_body for marker in GPU_MARKERS)
    if has_gpu_marker:
        return "gpu-marked"
    if path.name in GPU_SOURCE_FILES:
        return "gpu-file-unmarked"
    if path.name in HOST_BOUNDARY_FILES:
        return "host-boundary"
    if HOST_OR_DIAGNOSTIC_NAME_RE.search(name):
        return "host-or-diagnostic"
    return "unmarked-runtime-candidate"


def parse_subprograms(path: Path) -> list[Subprogram]:
    lines = (ROOT / path).read_text(encoding="utf-8").splitlines()
    result: list[Subprogram] = []
    index = 0
    while index < len(lines):
        match = START_RE.match(lines[index])
        if not match:
            index += 1
            continue
        name = match.group("name")
        start = index
        index += 1
        while index < len(lines) and not END_RE.match(lines[index]):
            index += 1
        end = min(index, len(lines) - 1)
        body = "\n".join(lines[start : end + 1])
        result.append(
            Subprogram(
                path=path,
                name=name,
                start_line=start + 1,
                end_line=end + 1,
                classification=classify(path, name, body),
            )
        )
        index += 1
    return result


def collect() -> list[Subprogram]:
    subprograms: list[Subprogram] = []
    for path in tracked_fortran():
        subprograms.extend(parse_subprograms(path))
    return subprograms


def print_summary(subprograms: list[Subprogram]) -> None:
    by_class = Counter(item.classification for item in subprograms)
    by_file: dict[Path, Counter[str]] = defaultdict(Counter)
    for item in subprograms:
        by_file[item.path][item.classification] += 1

    print("# Static GPU Port Inventory")
    print()
    print("Scope: tracked non-LVLSET Fortran subroutines/functions.")
    print("This is a source heuristic, not performance proof.")
    print()
    print("| Classification | Subprograms |")
    print("| --- | ---: |")
    for label in [
        "gpu-marked",
        "gpu-file-unmarked",
        "host-boundary",
        "host-or-diagnostic",
        "unmarked-runtime-candidate",
    ]:
        print(f"| `{label}` | {by_class[label]} |")

    print()
    print(
        "| File | GPU-marked | Host-boundary | Host/diagnostic | "
        "Unmarked runtime candidates |"
    )
    print("| --- | ---: | ---: | ---: | ---: |")
    for path in sorted(by_file):
        counts = by_file[path]
        if not (
            counts["gpu-marked"]
            or counts["host-boundary"]
            or counts["host-or-diagnostic"]
            or counts["unmarked-runtime-candidate"]
        ):
            continue
        print(
            f"| `{path}` | {counts['gpu-marked']} | "
            f"{counts['host-boundary']} | {counts['host-or-diagnostic']} | "
            f"{counts['unmarked-runtime-candidate']} |"
        )


def print_candidates(subprograms: list[Subprogram], limit: int) -> None:
    candidates = [
        item
        for item in subprograms
        if item.classification == "unmarked-runtime-candidate"
    ]
    print("# Unmarked Runtime Candidates")
    print()
    if not candidates:
        print("No unmarked runtime candidates found.")
        return
    for item in candidates[:limit]:
        print(f"- `{item.path}:{item.start_line}` `{item.name}`")
    if len(candidates) > limit:
        print(f"- ... {len(candidates) - limit} more")


def review_bucket(item: Subprogram) -> str:
    filename = item.path.name
    name = item.name

    if filename == "forcing.f90" and name == "lvlset_bridge_time":
        return "excluded-lvlset-bridge"

    if filename == "sgs_stag_util.f90":
        return "diagnostic-profiling"

    if filename in {"actuator_turbine_model.f90", "atm_base.f90", "linear_simple.f90"}:
        return "atm-host-model"

    if filename == "atm_input_util.f90":
        return "atm-mirror-lb-control"

    if filename == "atm_lesgo_interface.f90":
        return "atm-mirror-lb-control"

    if filename in {"turbines.f90", "turbine_indicator.f90"}:
        return "adm-cpu-fallback-profile"

    if filename in {"functions.f90"}:
        return "generic-helper-profile"

    if filename in {"inflow.f90", "fringe.f90"}:
        return "inflow-fringe-profile"

    if filename in {"scalars.f90", "stability.f90"}:
        return "scalar-init-fallback"

    if filename == "iwmles.f90":
        return "iwm-wallmodel-profile"

    if filename in {
        "convec.f90",
        "derivatives.f90",
        "divstress_uv.f90",
        "divstress_w.f90",
        "emul_complex.f90",
        "fft.f90",
        "hit_inflow.f90",
        "interpolag_Sdep.f90",
        "interpolag_Ssim.f90",
        "lagrange_Sdep.f90",
        "lagrange_Ssim.f90",
        "mpi_transpose_mod.f90",
        "scaledep_dynamic.f90",
        "std_dynamic.f90",
        "tridag_array.f90",
    }:
        return "cpu-fallback-compat"

    if filename == "test_filtermodule.f90":
        return "diagnostic-profiling"

    return "needs-manual-review"


def print_review(subprograms: list[Subprogram], limit_per_bucket: int) -> None:
    candidates = [
        item
        for item in subprograms
        if item.classification == "unmarked-runtime-candidate"
    ]
    by_bucket: dict[str, list[Subprogram]] = defaultdict(list)
    for item in candidates:
        by_bucket[review_bucket(item)].append(item)

    print("# Static Candidate Review Buckets")
    print()
    print(
        "These buckets classify unmarked runtime candidates for human review. "
        "They are not a substitute for paired CPU/GPU benchmarks."
    )
    print()
    print("| Review bucket | Candidates | Meaning |")
    print("| --- | ---: | --- |")
    for bucket in sorted(by_bucket):
        description = REVIEW_BUCKET_DESCRIPTIONS[bucket]
        print(f"| `{bucket}` | {len(by_bucket[bucket])} | {description} |")

    print()
    for bucket in sorted(by_bucket):
        print(f"## {bucket}")
        print()
        for item in by_bucket[bucket][:limit_per_bucket]:
            print(f"- `{item.path}:{item.start_line}` `{item.name}`")
        remaining = len(by_bucket[bucket]) - limit_per_bucket
        if remaining > 0:
            print(f"- ... {remaining} more")
        print()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report static non-LVLSET GPU-port inventory."
    )
    parser.add_argument(
        "--candidates",
        action="store_true",
        help="print unmarked runtime candidate subprograms instead of summary",
    )
    parser.add_argument(
        "--review",
        action="store_true",
        help="group unmarked runtime candidates into human-review buckets",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=80,
        help="candidate rows to print with --candidates, or per bucket with --review",
    )
    args = parser.parse_args()

    subprograms = collect()
    if args.review:
        print_review(subprograms, args.limit)
    elif args.candidates:
        print_candidates(subprograms, args.limit)
    else:
        print_summary(subprograms)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
