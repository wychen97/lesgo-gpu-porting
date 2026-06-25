#!/usr/bin/env bash

# Compile-only script for the Derecho CPU or A100 GPU build.
# This script does not submit or run a job. It only builds the executable and
# installs it in this case directory as one of:
#
#   ./lesgo-run-exe-cpu
#   ./lesgo-run-exe-gpu
#
# If using another cluster, edit the module load line and FC compiler wrapper.
# Usage:
#   ./compile_derecho.sh gpu
#   ./compile_derecho.sh cpu
#
# The default profile is gpu for backward compatibility.

# Bash safety options:
#   -e: stop if any command fails.
#   -u: stop if an undefined variable is used.
#   -o pipefail: stop if any command inside a pipeline fails.
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
    USE_GPU_AWARE_MPI=AUTO
    EXE_NAME=lesgo-mpi-ATM-lesgpu
    ;;
  cpu)
    USE_CPU_BUILD=ON
    USE_LES_GPU=OFF
    USE_GPU_AWARE_MPI=OFF
    EXE_NAME=lesgo-mpi-ATM
    ;;
  *)
    echo "Usage: $0 [gpu|cpu]" >&2
    exit 2
    ;;
esac

cd "${CASE_DIR}"

# Derecho module stack tested with NVHPC + CUDA + Cray MPICH.
# On another cluster, replace these modules with that site's compiler, CUDA,
# MPI, FFTW, and CMake modules.
module reset
module load nvhpc/25.9 cuda/12.9.0 cray-mpich/8.1.32 fftw/3.3.10 cmake/3.31.8

# Cray MPICH needs these for GPU-aware MPI and managed-memory device buffers.
# If your MPI is not Cray MPICH, check the site documentation before keeping
# these environment variables.
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1

# FC=ftn is the Derecho/Cray MPI Fortran wrapper.
# On another cluster this may be mpifort, mpif90, or another site wrapper.
#
# The current CMake file uses "hostname=derecho" for the tested Derecho
# NVHPC/CUDA/A100 compiler and link profile. If another cluster needs different
# compiler flags, add/change the matching profile in the root CMakeLists.txt
# and update this line.
FC="${FC:-ftn}" cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" \
  -Dhostname=derecho \
  -DUSE_MPI=ON \
  -DUSE_CPU_BUILD="${USE_CPU_BUILD}" \
  -DUSE_LES_GPU="${USE_LES_GPU}" \
  -DUSE_GPU_AWARE_MPI="${USE_GPU_AWARE_MPI}" \
  -DUSE_TURBINES=OFF \
  -DUSE_ATM=ON \
  -DUSE_CPS=OFF \
  -DUSE_SCALARS=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF

cmake --build "${BUILD_DIR}" -j "${BUILD_JOBS}"

# ATM build names from the root CMakeLists.txt:
#   GPU: lesgo + mpi + ATM + lesgpu -> lesgo-mpi-ATM-lesgpu
#   CPU: lesgo + mpi + ATM          -> lesgo-mpi-ATM
install -m 0755 "${BUILD_DIR}/${EXE_NAME}" "${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}"
cp "${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}" "${CASE_DIR}/lesgo-run-exe"

echo "Built ${CASE_DIR}/lesgo-run-exe-${BUILD_PROFILE}"
echo "Updated compatibility copy ${CASE_DIR}/lesgo-run-exe"
