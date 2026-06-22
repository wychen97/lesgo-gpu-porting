# LESGO Build Profiles

This file gives collaborator-facing CMake entry points for the published
optimized non-LVLSET branch.  The source still supports older host blocks, but
the validated handoff path is the Derecho NVHPC/A100 GPU build with LVLSET
disabled.

## User-Facing Cache Variables

These values can be supplied with `cmake -D...`:
Every public root CMake cache variable should appear exactly once in this
table.  The readiness gate checks this table against `CMakeLists.txt`.

| Variable | Values | Default | Purpose |
| --- | --- | --- | --- |
| `hostname` | `derecho`, `default`, legacy host names, `delta` compatibility alias | `derecho` | Selects host-specific compiler and library paths. |
| `WRITE_ENDIAN` | `DEFAULT`, `LITTLE`, `BIG` | `DEFAULT` | Controls output binary endian preprocessor flags. |
| `READ_ENDIAN` | `DEFAULT`, `LITTLE`, `BIG` | `DEFAULT` | Controls input binary endian preprocessor flags. |
| `USE_GPU_AWARE_MPI` | `AUTO`, `ON`, `OFF` | `AUTO` | Controls GPU-aware MPI halo/tridiagonal paths when `USE_LES_GPU=ON`. |

The `USE_*` feature options are documented in `docs/code_organization.md`.
Their defaults in `CMakeLists.txt` are neutral global defaults, not case
definitions.  Case-specific physics choices should be made in the case compile
script or in an explicit `cmake -D...` command.

## Canonical Production GPU Build

Use this profile for the current optimized non-LVLSET production path:

```bash
cmake -S . -B build-gpu-production \
  -Dhostname=derecho \
  -DUSE_MPI=ON \
  -DUSE_CPS=ON \
  -DUSE_ATM=ON \
  -DUSE_TURBINES=ON \
  -DUSE_LES_GPU=ON \
  -DUSE_GPU_AWARE_MPI=AUTO \
  -DUSE_SCALARS=ON \
  -DUSE_SCALARS_GPU=ON \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF
```

This is the path represented by the final validation record.  It keeps the LES
core, SGS, CPS, scalar GPU path, turbines, and ATM/structural coupling enabled.

## CPU Baseline Build

Use this profile for NVHPC CPU-only comparison runs on nodes where CUDA runtime
libraries should not be linked:

```bash
cmake -S . -B build-cpu-baseline \
  -Dhostname=derecho \
  -DUSE_CPU_BUILD=ON \
  -DUSE_MPI=ON \
  -DUSE_LES_GPU=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_LVLSET=OFF
```

`USE_CPU_BUILD=ON` keeps enough NVHPC syntax support for dormant CUDA Fortran
declarations while avoiding GPU library links.

## Minimal Configure Smoke

Use this profile for a quick source-tree sanity check on a local machine.  It is
not a performance or physics validation run:

```bash
cmake -S . -B build-configure-smoke \
  -DUSE_MPI=OFF \
  -DUSE_CPS=OFF \
  -DUSE_ATM=OFF \
  -DUSE_TURBINES=OFF \
  -DUSE_LES_GPU=OFF \
  -DUSE_SCALARS=OFF \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=OFF \
  -DUSE_CGNS=OFF
```

The readiness wrapper can run this same configure smoke when CMake and a
suitable compiler wrapper are available:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

## Optional Module Profiles

Scalar transport with GPU scalar kernels requires both scalar and LES GPU
support:

```bash
-DUSE_SCALARS=ON -DUSE_SCALARS_GPU=ON -DUSE_LES_GPU=ON
```

HIT inflow is an opt-in inflow source.  It is not part of the standard four
wind-farm validation cases, but the branch now keeps the HIT CPU and GPU source
membership explicit.  A minimal GPU configure should use:

```bash
cmake -S . -B build-hit-gpu \
  -Dhostname=derecho \
  -DUSE_MPI=ON \
  -DUSE_CPS=OFF \
  -DUSE_ATM=OFF \
  -DUSE_TURBINES=OFF \
  -DUSE_LES_GPU=ON \
  -DUSE_GPU_AWARE_MPI=AUTO \
  -DUSE_SCALARS=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_LVLSET=OFF \
  -DUSE_HIT=ON \
  -DUSE_DYN_TN=OFF \
  -DUSE_CGNS=OFF
```

The corresponding `lesgo.conf` must set `inflow_type = 2` and provide the HIT
keys in the `FLOW_COND` block: `UP_IN`, `TI_OUT`, `LX_HIT`, `LY_HIT`,
`LZ_HIT`, `NX_HIT`, `NY_HIT`, `NZ_HIT`, `U_FILE`, `V_FILE`, and `W_FILE`.
The three HIT velocity files are currently read as ASCII scalar lists, and
`restartHIT.dat` stores the swept HIT plane location.

CGNS output and LVLSET are optional branches.  LVLSET remains outside the
optimized production path and should not be used as evidence that a GPU
hot-path change is validated.

## Validation Rule

Changing CMake option defaults, compiler flags, source membership, or GPU-aware
MPI behavior requires:

```bash
python3 tools/check_branch_readiness.py
```

and, on the intended compiler stack:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

For the optional HIT plus LES GPU source set, use:

```bash
python3 tools/check_branch_readiness.py --with-hit-cmake-configure
```
