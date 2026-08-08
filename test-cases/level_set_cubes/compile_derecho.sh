#!/usr/bin/env bash

# Compile-only Derecho example for the Level Set CPU reference or A100 GPU
# implementation. This script does not submit or run a job.
#
# Usage:
#   ./compile_derecho.sh gpu
#   ./compile_derecho.sh cpu
#
# On another cluster, replace the module stack, compiler wrapper, hostname
# profile, and GPU-aware MPI setting with that site's supported equivalents.

set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ROOT_DIR="$(cd "${CASE_DIR}/../.." && pwd -P)"
BUILD_PROFILE="${1:-${BUILD_PROFILE:-gpu}}"
BUILD_DIR="${CASE_DIR}/build-derecho-${BUILD_PROFILE}"
BUILD_JOBS="${BUILD_JOBS:-8}"

case "${BUILD_PROFILE}" in
  gpu)
    USE_CPU_BUILD=OFF
    USE_LES_GPU=ON
    USE_LVLSET_GPU=ON
    USE_GPU_AWARE_MPI=AUTO
    EXE_NAME=lesgo-mpi-ls-lsgpu-lesgpu
    ;;
  cpu)
    USE_CPU_BUILD=ON
    USE_LES_GPU=OFF
    USE_LVLSET_GPU=OFF
    USE_GPU_AWARE_MPI=OFF
    EXE_NAME=lesgo-mpi-ls
    ;;
  *)
    echo "Usage: $0 [gpu|cpu]" >&2
    exit 2
    ;;
esac

module reset
module load nvhpc/26.1 cuda/12.9.0 cray-mpich/8.1.32 fftw/3.3.10 cmake/3.31.8

# These variables are required by Cray MPICH when the GPU profile passes
# device buffers directly to MPI. Other MPI implementations use site-specific
# controls; do not assume these names apply outside a Cray system.
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1

# FC=ftn is Derecho's MPI Fortran wrapper. Common alternatives on other
# clusters are mpifort and mpif90.
FC="${FC:-ftn}" cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -Dhostname=derecho \
  -DUSE_MPI=ON \
  -DUSE_CPU_BUILD="${USE_CPU_BUILD}" \
  -DUSE_LES_GPU="${USE_LES_GPU}" \
  -DUSE_GPU_AWARE_MPI="${USE_GPU_AWARE_MPI}" \
  -DUSE_LVLSET=ON \
  -DUSE_LVLSET_GPU="${USE_LVLSET_GPU}" \
  -DUSE_TURBINES=OFF \
  -DUSE_ATM=OFF \
  -DUSE_CPS=OFF \
  -DUSE_SCALARS=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF

cmake --build "${BUILD_DIR}" -j "${BUILD_JOBS}"

install -m 0755 "${BUILD_DIR}/${EXE_NAME}" \
  "${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}"

echo "Built ${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}"
