#!/usr/bin/env bash

# Compile optional non-LVLSET source groups on Derecho from the current
# committed HEAD.  This catches module-interface drift in source sets that are
# not exercised by the four public presentation cases.
#
# Usage from the repository root:
#
#   tools/run_derecho_optional_compile_matrix.sh
#
# Optional environment variables:
#
#   DERECHO_COMPILE_ROOT   Remote output directory. Defaults to a timestamped
#                          directory under /glade/work/wchen/lesgo_versions.
#   BUILD_JOBS             Parallel build jobs passed to CMake. Defaults to 8.
#
# This script validates isolated optional features and common paired feature
# combinations.  It intentionally excludes LVLSET because the public optimized
# branch does not support or claim LVLSET validation.  A small serial/non-MPI
# group is included because `USE_MPI=OFF` changes many Fortran import guards.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "${repo_root}"

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "Working tree has uncommitted changes; commit first so Derecho tests a reproducible source snapshot." >&2
    exit 2
fi

commit="$(git rev-parse --short HEAD)"
timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
remote_root="${DERECHO_COMPILE_ROOT:-/glade/work/wchen/lesgo_versions/benchmarks/derecho_optional_compile_matrix_${commit}_${timestamp}}"
build_jobs="${BUILD_JOBS:-8}"

ssh derecho bash -s <<REMOTE
set -euo pipefail
mkdir -p "${remote_root}"
rm -rf "${remote_root}/src"
mkdir -p "${remote_root}/src"
REMOTE

git archive HEAD | ssh derecho "tar -xf - -C '${remote_root}/src'"
ssh derecho "printf '%s\n' '${commit}' > '${remote_root}/src/LOCAL_COMMIT.txt'"

printf 'Derecho optional compile root: %s\n' "${remote_root}"
printf 'Commit: %s\n' "${commit}"

ssh derecho bash -s -- "${remote_root}" "${build_jobs}" <<'REMOTE'
set -euo pipefail

remote_root="$1"
build_jobs="$2"
src="${remote_root}/src"

module reset
module load nvhpc/25.9 cuda/12.9.0 cray-mpich/8.1.32 fftw/3.3.10 cmake/3.31.8

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1

base_common=(
    -Dhostname=derecho
    -DUSE_MPI=ON
    -DUSE_LVLSET=OFF
    -DUSE_CGNS=OFF
)

gpu_common=(
    "${base_common[@]}"
    -DUSE_CPU_BUILD=OFF
    -DUSE_LES_GPU=ON
    -DUSE_GPU_AWARE_MPI=AUTO
)

cpu_common=(
    "${base_common[@]}"
    -DUSE_CPU_BUILD=ON
    -DUSE_LES_GPU=OFF
    -DUSE_GPU_AWARE_MPI=OFF
)

serial_base_common=(
    -Dhostname=derecho
    -DUSE_MPI=OFF
    -DUSE_CPS=OFF
    -DUSE_LVLSET=OFF
    -DUSE_CGNS=OFF
)

serial_gpu_common=(
    "${serial_base_common[@]}"
    -DUSE_CPU_BUILD=OFF
    -DUSE_LES_GPU=ON
    -DUSE_GPU_AWARE_MPI=OFF
)

serial_cpu_common=(
    "${serial_base_common[@]}"
    -DUSE_CPU_BUILD=ON
    -DUSE_LES_GPU=OFF
    -DUSE_GPU_AWARE_MPI=OFF
)

run_config() {
    local name="$1"
    shift
    local build_dir="${src}/build-derecho-optional-${name}"
    local log="${remote_root}/compile_optional_${name}.log"

    printf '\n== optional %s ==\n' "${name}"
    rm -rf "${build_dir}"
    if (
        cd "${src}"
        FC="${FC:-ftn}" cmake -S "${src}" -B "${build_dir}" "$@" > "${log}" 2>&1
        cmake --build "${build_dir}" -j "${build_jobs}" >> "${log}" 2>&1
    ); then
        printf 'PASS optional %s %s\n' "${name}" "${log}"
        grep -E "Built target|Linking Fortran executable" "${log}" | tail -n 5 || true
        rm -rf "${build_dir}"
    else
        local rc=$?
        printf 'FAIL optional %s rc=%s %s\n' "${name}" "${rc}" "${log}"
        tail -n 140 "${log}"
        exit "${rc}"
    fi
}

run_config hit_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=ON -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config hit_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=ON -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config cps_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=ON -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config cps_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=ON -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config scalars_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON \
    -DUSE_DYN_TN=OFF

run_config scalars_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config cps_scalars_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=ON -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON \
    -DUSE_DYN_TN=OFF

run_config cps_scalars_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=ON -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config hit_scalars_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=ON -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON \
    -DUSE_DYN_TN=OFF

run_config hit_scalars_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=ON -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config atm_scalars_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=ON \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON \
    -DUSE_DYN_TN=OFF

run_config atm_scalars_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=ON \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config turbines_scalars_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=ON -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON \
    -DUSE_DYN_TN=OFF

run_config turbines_scalars_cpu \
    "${cpu_common[@]}" \
    -DUSE_TURBINES=ON -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config dyntn_gpu \
    "${gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_CPS=OFF -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=ON

run_config serial_core_gpu \
    "${serial_gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config serial_core_cpu \
    "${serial_cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config serial_hit_gpu \
    "${serial_gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_HIT=ON -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config serial_scalars_gpu \
    "${serial_gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON \
    -DUSE_DYN_TN=OFF

run_config serial_scalars_cpu \
    "${serial_cpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=OFF \
    -DUSE_HIT=OFF -DUSE_SCALARS=ON -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config serial_atm_gpu \
    "${serial_gpu_common[@]}" \
    -DUSE_TURBINES=OFF -DUSE_ATM=ON \
    -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

run_config serial_turbines_gpu \
    "${serial_gpu_common[@]}" \
    -DUSE_TURBINES=ON -DUSE_ATM=OFF \
    -DUSE_HIT=OFF -DUSE_SCALARS=OFF -DUSE_SCALARS_GPU=OFF \
    -DUSE_DYN_TN=OFF

printf '\nAll optional Derecho compile profiles passed.\n'
REMOTE

printf '\nAll requested optional Derecho compile profiles passed for %s.\n' "${commit}"
