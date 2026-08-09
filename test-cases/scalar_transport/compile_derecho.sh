#!/usr/bin/env bash

# Compile the compact scalar-transport example on Derecho.
#
# Usage:
#   ./compile_derecho.sh gpu
#   ./compile_derecho.sh cpu
#
# Both passive and active scalar variants use the same executable. Their
# difference is entirely in lesgo.conf, selected by submit_derecho.pbs.

set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ROOT_DIR="$(cd "${CASE_DIR}/../.." && pwd -P)"
BUILD_PROFILE="${1:-gpu}"
BUILD_DIR="${CASE_DIR}/build-derecho-${BUILD_PROFILE}"
BUILD_JOBS="${BUILD_JOBS:-8}"

case "${BUILD_PROFILE}" in
  gpu)
    USE_CPU_BUILD=OFF
    USE_LES_GPU=ON
    USE_SCALARS_GPU=ON
    USE_GPU_AWARE_MPI=AUTO
    EXE_NAME=lesgo-mpi-lesgpu-scalgpu-scalars
    ;;
  cpu)
    USE_CPU_BUILD=ON
    USE_LES_GPU=OFF
    USE_SCALARS_GPU=OFF
    USE_GPU_AWARE_MPI=OFF
    EXE_NAME=lesgo-mpi-scalars
    ;;
  *)
    echo "Usage: $0 [gpu|cpu]" >&2
    exit 2
    ;;
esac

module reset
module load nvhpc/26.1 cuda/12.9.0 cray-mpich/8.1.32 fftw/3.3.10 cmake/3.31.8

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1

FC="${FC:-ftn}" cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -Dhostname=derecho \
  -DUSE_MPI=ON \
  -DUSE_CPU_BUILD="${USE_CPU_BUILD}" \
  -DUSE_LES_GPU="${USE_LES_GPU}" \
  -DUSE_GPU_AWARE_MPI="${USE_GPU_AWARE_MPI}" \
  -DUSE_TURBINES=OFF \
  -DUSE_ATM=OFF \
  -DUSE_CPS=OFF \
  -DUSE_SCALARS=ON \
  -DUSE_SCALARS_GPU="${USE_SCALARS_GPU}" \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF

cmake --build "${BUILD_DIR}" -j "${BUILD_JOBS}"
install -m 0755 "${BUILD_DIR}/${EXE_NAME}" \
  "${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}"

echo "Built ${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}"
