# LESGO GPU Porting

This repository is a GPU-porting branch of LESGO organized to stay close to
the upstream JHU LESGO layout:

- solver source files and the main `CMakeLists.txt` remain at repository root;
- reusable documentation lives under `docs/`;
- runnable examples live under `test-cases/`.

## Documentation

The maintained documentation website is:

**[LESGO GPU Porting Guide](https://wychen97.github.io/lesgo-gpu-docs/)**

Start with these pages:

- [Public Test Cases](https://wychen97.github.io/lesgo-gpu-docs/gpu/test-cases/)
- [Build and Runtime](https://wychen97.github.io/lesgo-gpu-docs/gpu/build-runtime/)
- [GPU Architecture](https://wychen97.github.io/lesgo-gpu-docs/gpu/architecture/)
- [Validation and Performance](https://wychen97.github.io/lesgo-gpu-docs/gpu/validation-performance/)

The website explains the current release and normal user workflow. Detailed
source audits and validation evidence remain versioned with the code under
`docs/`.

The public test-case suite keeps four progression cases and four compact
optional-physics examples:

1. `test-cases/channel_flow` - LES core, no turbine model.
2. `test-cases/adm_disk` - actuator disk model.
3. `test-cases/atm_line` - actuator line / 5 MW turbine model.
4. `test-cases/large_windfarm_3072x384x400_60turbines` - 60-turbine large case.
5. `test-cases/level_set_cubes` - immersed-surface Level Set example.
6. `test-cases/scalar_transport` - passive and active scalar transport.
7. `test-cases/concurrent_precursor` - velocity/scalar concurrent precursor.
8. `test-cases/inflow_and_forcing` - HIT, shifted inflow, Coriolis, and sponge.

Generated output, queue logs, local monitor files, and compiled binaries are
excluded. Each case has two Derecho example files:

```text
compile_derecho.sh
submit_derecho.pbs
```

`compile_derecho.sh` only compiles. It installs a case-local profile
executable, normally one of:

```text
lesgo-run-exe-cpu
lesgo-run-exe-gpu
```

`submit_derecho.pbs` only submits/runs the already-built executable. The two
steps are intentionally separate so users can inspect and modify compilation
settings independently from queue settings.

Case-local executables are intentionally ignored by git because they are specific to the
compiler, GPU architecture, MPI stack, and cluster.

## Minimal GPU Build Example

Run this from one of the case directories, for example
`test-cases/channel_flow`:

```bash
module reset
module load nvhpc/26.1 cuda/12.9.0 cray-mpich/8.1.32 fftw/3.3.10 cmake/3.31.8
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1

CASE=$PWD
SRC=$(cd ../.. && pwd)
BUILD=$CASE/builds/channel_gpu240

FC=ftn cmake -S "$SRC" -B "$BUILD" \
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

cmake --build "$BUILD" -j 8
```

`cray-mpich/8.1.32` remains the supported Derecho MPI module. The available
`cray-mpich/9.0.0` module identifies itself as functional-only pre-release
software and must not be used for benchmark or production results. CMake
recognizes the `NCAR_ROOT_FFTW` path exported by Derecho's NVHPC 26.1 FFTW
module, so the inactive legacy `ncarcompilers` wrapper is not required.

`BUILD=$CASE/builds/channel_gpu240` is just a shell variable naming the
out-of-source CMake build directory. It was not part of the original LESGO
CMake design; it is the standard CMake `-B` build-tree location.

## Main GPU CMake Switches

- `USE_LES_GPU=ON`: enables the optimized OpenACC explicit-residency LES core.
- `USE_GPU_AWARE_MPI=AUTO`: enables GPU-aware MPI when the cluster MPI stack
  exposes the required support.
- `USE_TURBINES=ON`: enables the actuator disk model.
- `USE_ATM=ON`: enables the actuator line / turbine structural-path model.
- `USE_SCALARS=ON` and `USE_SCALARS_GPU=ON`: enable scalar transport and the
  optimized scalar GPU path.
- `USE_CPS=ON`: enables concurrent precursor simulation.
- `USE_LVLSET=ON` and `USE_LVLSET_GPU=ON`: enable level-set physics and its
  optimized GPU implementation; see `test-cases/level_set_cubes`.

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
   `USE_GPU_AWARE_MPI=AUTO` works on the tested Cray/NVHPC setup. Use `ON` only
   if the MPI library is known to support device buffers. Use `OFF` if the MPI
   stack is not GPU-aware.

4. Launch command:
   Derecho provides `set_gpu_rank`. If another cluster does not, plain
   `mpiexec` may work for one GPU. On Slurm systems, translate
   `submit_derecho.pbs` to an `sbatch` file and use the site's recommended
   `srun` command.

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
