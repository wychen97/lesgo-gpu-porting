#!/usr/bin/env bash

# Compile-only script for the Derecho A100 GPU build.
# This script does not submit or run a job. It only builds the executable and
# installs it in this case directory as:
#
#   ./lesgo-run-exe
#
# If using another cluster, edit the module load line and FC compiler wrapper.

# Bash safety options:
#   -e: stop if any command fails.
#   -u: stop if an undefined variable is used.
#   -o pipefail: stop if any command inside a pipeline fails.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ROOT_DIR="$(cd "${CASE_DIR}/../.." && pwd -P)"
BUILD_DIR="${CASE_DIR}/build-derecho-gpu"
BUILD_JOBS="${BUILD_JOBS:-8}"

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
  -DUSE_CPU_BUILD=OFF \
  -DUSE_LES_GPU=ON \
  -DUSE_GPU_AWARE_MPI=AUTO \
  -DUSE_TURBINES=OFF \
  -DUSE_ATM=OFF \
  -DUSE_CPS=OFF \
  -DUSE_SCALARS=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF

cmake --build "${BUILD_DIR}" -j "${BUILD_JOBS}"

# Channel-flow build name from the root CMakeLists.txt:
#   lesgo + mpi + lesgpu -> lesgo-mpi-lesgpu
install -m 0755 "${BUILD_DIR}/lesgo-mpi-lesgpu" "${CASE_DIR}/lesgo-run-exe"

echo "Built ${CASE_DIR}/lesgo-run-exe"
