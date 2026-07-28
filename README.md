# LESGO CPU Baseline Branch

This branch is the CPU-default counterpart to the LESGO GPU-porting repository.
It is organized to stay close to the upstream JHU LESGO layout:

- solver source files and the main `CMakeLists.txt` remain at repository root;
- reusable documentation lives under `docs/`;
- runnable examples live under `test-cases/`.

The public test-case suite is intentionally small. It contains only the four
representative cases used for the current validation/progression discussion:

1. `test-cases/channel_flow` - LES core, no turbine model.
2. `test-cases/adm_disk` - actuator disk model.
3. `test-cases/atm_line` - actuator line / 5 MW turbine model.
4. `test-cases/large_windfarm_3072x384x400_60turbines` - 60-turbine large case.

Generated output, queue logs, local monitor files, and compiled binaries are
excluded. Each case has two Derecho example files:

```text
compile_derecho.sh
submit_derecho.pbs
```

`compile_derecho.sh` only compiles. It installs one case-local executable:

```text
lesgo-run-exe
```

`submit_derecho.pbs` only submits/runs the already-built executable. The two
steps are intentionally separate so users can inspect and modify compilation
settings independently from queue settings.

`lesgo-run-exe` is intentionally ignored by git because it is specific to the
compiler, MPI stack, and cluster.

## Minimal CPU Build Example

Run this from one of the case directories, for example
`test-cases/channel_flow`:

```bash
module reset
module load nvhpc/25.9 cray-mpich/8.1.32 fftw/3.3.10

CASE=$PWD
SRC=$(cd ../.. && pwd)
BUILD=$CASE/builds/channel_cpu240

FC=ftn cmake -S "$SRC" -B "$BUILD" \
  -Dhostname=derecho \
  -DUSE_MPI=ON \
  -DUSE_CPU_BUILD=ON \
  -DUSE_LES_GPU=OFF \
  -DUSE_GPU_AWARE_MPI=OFF \
  -DUSE_TURBINES=OFF \
  -DUSE_ATM=OFF \
  -DUSE_CPS=OFF \
  -DUSE_SCALARS=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF

cmake --build "$BUILD" -j 8
```

`BUILD=$CASE/builds/channel_cpu240` is just a shell variable naming the
out-of-source CMake build directory. It was not part of the original LESGO
CMake design; it is the standard CMake `-B` build-tree location.

## Main CPU CMake Switches

- `USE_CPU_BUILD=ON`: uses the CPU-baseline NVHPC flags on Derecho so dormant
  CUDA syntax still compiles without linking GPU runtime libraries.
- `USE_LES_GPU=OFF`: disables the OpenACC explicit-residency LES core.
- `USE_GPU_AWARE_MPI=OFF`: disables GPU-aware MPI paths.
- `USE_TURBINES=ON`: enables the actuator disk model.
- `USE_ATM=ON`: enables the actuator line / turbine structural-path model.
- `USE_SCALARS=ON`: enables scalar transport in the CPU path.
- `USE_CPS=ON`: enables concurrent precursor simulation.
- `USE_LVLSET=ON`: enables level-set physics. This path is kept in source but
  is not part of the current optimized public test suite.

## What To Change On Another Cluster

Each case has one compile file and one submission file. When porting them to a
new cluster, edit only these parts first:

1. Module loads:
   Replace the `module load nvhpc cuda cray-mpich fftw cmake` line with that
   cluster's compiler, CUDA, MPI, FFTW, and CMake modules. Do this in both
   files so the compile and runtime environments match.

2. MPI wrapper:
   Change `FC=ftn` in `compile_derecho.sh` if the cluster uses another MPI
   Fortran wrapper, for example `mpifort`, `mpif90`, or a site-specific NVHPC
   wrapper.

3. GPU-aware MPI:
   Keep `USE_GPU_AWARE_MPI=OFF` on this CPU branch.

4. Launch command:
   The CPU submit scripts use plain `mpiexec`. On Slurm systems, translate
   `submit_derecho.pbs` to an `sbatch` file and use the site's recommended
   `srun` or `mpiexec` command.

5. Case physics:
   Do not change `USE_TURBINES` or `USE_ATM` unless you intentionally switch
   between channel flow, ADM, and ATM physics.

Concrete Derecho sequence:

```bash
cd test-cases/channel_flow
./compile_derecho.sh
qsub submit_derecho.pbs
```

Edit `lesgo.conf` directly when changing grid size, timestep count, or output
controls.
