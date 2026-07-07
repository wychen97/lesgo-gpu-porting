# Collaborator Quickstart

This branch is meant to be usable by people who already know the CPU LESGO
source.  The GPU work is organized as an extension of the familiar CPU layout,
not as a replacement architecture.

## What Stayed Familiar

- `main.f90` still owns timestep orchestration.
- The main physics modules keep their CPU-version names and roles:
  `sgs_stag_util.f90`, `press_stag_array.f90`, `tridag_array.f90`,
  `actuator_turbine_model.f90`, `atm_lesgo_interface.f90`, `scalars.f90`,
  `iwmles.f90`, and related helper files.
- Runtime case files still choose physics through LESGO input settings such as
  `sgs_model`, turbine inputs, inflow settings, and boundary conditions.
- CPU fallback paths remain important.  Do not remove a CPU branch just because
  the production GPU path is the current validation target.

## What The GPU Branch Adds

- CMake feature options select the production module set.  Start with
  `docs/build_profiles.md` before changing build flags.
- GPU-specific modules such as `convec_gpu.f90`, `derivatives_gpu.f90`,
  `fft_gpu.f90`, `press_gpu.f90`, `sgs_gpu.f90`, and
  `lagrange_Sdep_gpu.f90` hold the device kernels and GPU helper paths.
- `PPLES_GPU`, `PPSGS_GPU`, `PPSCALARS_GPU`, `PPATM`, and related preprocessor
  symbols route CPU/GPU behavior.
- GPU comments document data residency, ownership, and CPU/GPU synchronization
  contracts.  Treat these comments as part of the interface.
- Readiness scripts in `tools/` guard source hygiene, build-option docs,
  source inventories, GPU contracts, SGS dispatch constants, and known indexing
  mistakes.

## First Commands

Check the branch before editing:

```bash
python3 tools/check_branch_readiness.py
```

If CMake and the intended compiler wrapper are available:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

For optional inflow or GPU source-membership changes:

```bash
python3 tools/check_branch_readiness.py --with-hit-cmake-configure
```

Configure the validated non-LVLSET GPU profile:

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

Configure a CPU baseline build:

```bash
cmake -S . -B build-cpu-baseline \
  -Dhostname=derecho \
  -DUSE_CPU_BUILD=ON \
  -DUSE_MPI=ON \
  -DUSE_LES_GPU=OFF \
  -DUSE_SCALARS_GPU=OFF \
  -DUSE_LVLSET=OFF
```

## Where To Look First

- Build options and standard profiles: `docs/build_profiles.md`
- Source layout and timestep order: `docs/code_organization.md`
- Module ownership and validation expectations: `docs/gpu_module_contracts.md`
- Runtime `LESGO_*` switches: `docs/environment_switches.md`
- Presentation-oriented CPU/GPU cases: `test-cases/`
- Known cleanup order: `docs/refactor_backlog.md`
- Local checks and tooling: `tools/README.md`

## Current Validation Scope

The current production target is the optimized non-LVLSET GPU path with LES,
SGS, CPS, scalar GPU transport, turbines, ATM/structural coupling, diagnostics,
and output enabled.  The validated runtime SGS values are:

```text
sgs=false: supported
sgs_model=1: supported
sgs_model=2: supported
sgs_model=3: supported
sgs_model=4: supported
sgs_model=5: supported
sgs_model=6/7: unsupported runtime values, guard-fail as expected
```

LVLSET files remain tracked for completeness, but LVLSET is outside the current
optimized production scope.

## Editing Rules

- Keep changes small and focused.
- Separate behavior changes from readability-only cleanup.
- Keep CPU fallback behavior correct when changing a GPU path.
- Update the relevant doc in the same commit when a build option, source
  ownership boundary, runtime switch, or validation expectation changes.
- Do not commit local benchmark output, submitted-job records, or generated
  binaries.

## Known Remaining Cleanup

The branch is ready for collaborator use, but it is not a fully modular rewrite.
The largest files still need gradual cleanup in this order:

1. `atm_lesgo_interface.f90`
2. `sgs_stag_util.f90`
3. `actuator_turbine_model.f90`
4. `lagrange_Sdep_gpu.f90`
5. `scalars.f90`
6. `press_stag_array.f90`
7. `tridag_array.f90`
8. `iwmles.f90`

Use `docs/refactor_backlog.md` before starting that work.
