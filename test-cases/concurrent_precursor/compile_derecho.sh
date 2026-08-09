#!/usr/bin/env bash

# Compile the concurrent-precursor example on Derecho.
#
# Usage:
#   ./compile_derecho.sh gpu velocity
#   ./compile_derecho.sh gpu scalar
#   ./compile_derecho.sh cpu velocity
#   ./compile_derecho.sh cpu scalar

set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ROOT_DIR="$(cd "${CASE_DIR}/../.." && pwd -P)"
BUILD_PROFILE="${1:-gpu}"
CASE_VARIANT="${2:-velocity}"
BUILD_DIR="${CASE_DIR}/build-derecho-${BUILD_PROFILE}-${CASE_VARIANT}"
BUILD_JOBS="${BUILD_JOBS:-8}"

case "${CASE_VARIANT}" in
  velocity)
    USE_SCALARS=OFF
    USE_SCALARS_GPU=OFF
    SCALAR_EXE_SUFFIX=""
    ;;
  scalar)
    USE_SCALARS=ON
    SCALAR_EXE_SUFFIX="-scalars"
    ;;
  *)
    echo "Usage: $0 [gpu|cpu] [velocity|scalar]" >&2
    exit 2
    ;;
esac

case "${BUILD_PROFILE}" in
  gpu)
    USE_CPU_BUILD=OFF
    USE_LES_GPU=ON
    USE_GPU_AWARE_MPI=AUTO
    if [[ "${CASE_VARIANT}" == "scalar" ]]; then
      USE_SCALARS_GPU=ON
      EXE_NAME=lesgo-mpi-cps-lesgpu-scalgpu-scalars
    else
      EXE_NAME=lesgo-mpi-cps-lesgpu
    fi
    ;;
  cpu)
    USE_CPU_BUILD=ON
    USE_LES_GPU=OFF
    USE_GPU_AWARE_MPI=OFF
    USE_SCALARS_GPU=OFF
    EXE_NAME="lesgo-mpi-cps${SCALAR_EXE_SUFFIX}"
    ;;
  *)
    echo "Usage: $0 [gpu|cpu] [velocity|scalar]" >&2
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
  -DUSE_CPS=ON \
  -DUSE_SCALARS="${USE_SCALARS}" \
  -DUSE_SCALARS_GPU="${USE_SCALARS_GPU}" \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF

cmake --build "${BUILD_DIR}" -j "${BUILD_JOBS}"
install -m 0755 "${BUILD_DIR}/${EXE_NAME}" \
  "${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}-${CASE_VARIANT}"

echo "Built ${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}-${CASE_VARIANT}"
