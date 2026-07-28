# GPU Module Contracts

This document records the working contracts for the optimized non-LVLSET GPU
branch.  Use it when deciding where a change belongs and which validation is
needed before accepting it.

The goal is not to freeze the implementation.  The goal is to make ownership
and CPU/GPU boundaries explicit so future optimizations do not accidentally
reintroduce hidden transfers, stale host data, or one-SGS-value-only behavior.

## Contract Summary

| Area | Primary files | Owns | Must preserve |
| --- | --- | --- | --- |
| Global field storage | `sim_param.f90` | Core velocity, gradient, RHS, pressure-gradient, stress, div-stress, and force arrays | Clear CPU/GPU ownership and persistent device mirrors in the PPLES GPU path. |
| Driver orchestration | `main.f90` | Timestep order, timing blocks, periodic diagnostics, and hot-path dispatch | Existing data dependency order and the ATM phase split. |
| Spectral transforms | `fft_gpu.f90`, `test_filtermodule.f90` | cuFFT plan lifecycle and batched spectral calls | Shared async queue usage and plan/scratch lifetime. |
| Velocity derivatives | `derivatives_gpu.f90` | `filt_da_gpu`, `ddz_uv_gpu`, `ddz_w_gpu`, `ddx_gpu`, `ddy_gpu`, `ddxy_gpu` | Device-resident inputs/outputs and no implicit full-field host copies. |
| Convection | `convec_gpu.f90` | Dealiased `-(u cross omega)` RHS assembly | Async overlap contract with ATM phase 1. |
| SGS and div-stress | `sgs_stag_util.f90`, `sgs_gpu.f90`, `lagrange_Sdep_gpu.f90` | Runtime SGS dispatch, stress tensors, Lagrangian averaging, div-stress kernels | Disabled SGS path plus supported runtime values `1..5`. |
| Pressure and tridiagonal solve | `press_gpu.f90`, `press_stag_array.f90`, `tridag_gpu.f90`, `tridag_array.f90` | Pressure RHS, pressure solve, GPU-aware MPI pressure pipeline | GPU-aware and host-staged fallback correctness, with timings reported separately. |
| Scalars | `scalars.f90`, `stability.f90` | Scalar transport, scalar derivatives, scalar halo handling | Scalar-off, scalar CPU, and scalar GPU build behavior. |
| CPS / inflow | `concurrent_precursor.f90`, `inflow.f90`, `shifted_inflow.f90`, `hit_inflow*.f90` | Inflow source dispatch and concurrent precursor exchange | Production CPS path and optional HIT/shifted inflow isolation. |
| ATM / structural coupling | `atm_lesgo_interface.f90`, `actuator_turbine_model.f90`, `atm_base.f90`, `atm_input_util.f90` | Turbine sampling, force application, induced-velocity correction, structural solver coupling | Structure-off and structure-on correctness, turbine power output, and phase ordering. |
| Diagnostics / output | `io.f90`, `rmsdiv.f90`, `time_average.f90`, `stat_defs.f90`, `clocks.f90`, `cuda_mpi_debug.f90` | Output, restart, statistics, divergence/energy diagnostics, timing/debug prints | Periodic diagnostic costs are reported separately from regular-step timing. |
| LVLSET | `level_set*.f90`, `trees_*_ls.f90` | Legacy LVLSET support | Not part of the optimized production path unless explicitly revalidated. |

## CMake Source Group Map

Root `CMakeLists.txt` groups source files before feature options append them to
the executable.  Keep these groups aligned with the ownership contracts above.

| CMake source group | Contract area |
| --- | --- |
| `COMMON_DRIVER_SOURCES` | Driver orchestration, global setup, shared utilities. |
| `COMMON_LES_CORE_SOURCES` | LES core, derivatives, SGS, pressure, wall stress, and IWM. |
| `COMMON_FORCING_SOURCES` | Forcing, boundary conditions, and inflow. |
| `COMMON_DIAGNOSTIC_SOURCES` | Diagnostics, output, restart, statistics, and timing. |
| `MPI_SOURCES` | MPI transpose and domain metadata. |
| `CPS_SOURCES` | CPS / inflow. |
| `HIT_SOURCES` | Optional HIT inflow. |
| `LVLSET_SOURCES` | Deferred LVLSET support. |
| `TURBINE_SOURCES` | Legacy actuator-disk turbine model. |
| `ATM_SOURCES` | ATM / structural coupling. |
| `LES_GPU_SOURCES` | Optimized LES GPU kernels and GPU pressure/SGS helpers. |
| `SCALAR_SOURCES` | Scalars and stability coupling. |

## Global Field Storage

`sim_param.f90` is the central ownership point for LES arrays.  Do not add a
new hot-path field in another module unless there is a clear reason it cannot
live with the existing field declarations or with that module's persistent
scratch.

Rules:

- PPLES GPU hot arrays should have persistent device mirrors.
- Host copies of hot arrays may be stale between explicit update boundaries.
- New full-field host/device transfers must be documented at the call site.
- New force, pressure, or stress arrays should state whether `main.f90`, SGS,
  pressure, or ATM owns the authoritative value after each timestep phase.

## Driver Orchestration

`main.f90` owns the timestep order.  Moving work across phases can be a valid
optimization, but only after checking data dependencies and periodic side
effects.

Do not casually move:

- ATM phase 1 relative to `convec_gpu`.
- RHS assembly relative to div-stress and convection.
- Pressure-gradient update relative to `press_stag_array_gpu`.
- Diagnostics relative to `output_loop`, `energy`, and `rmsdiv`.

Any change in `main.f90` should report whether measured timing is from a
regular step, an SGS-update step, or an output/statistics step.

## Spectral Transform Contract

`fft_gpu.f90` owns cuFFT plan setup, stream binding, and plan teardown.  GPU
modules should call its execution wrappers rather than creating local cuFFT
plans.

Rules:

- Keep cuFFT work on the same async queue used by the surrounding OpenACC
  kernels unless a new queue policy is documented.
- Do not destroy/recreate plans inside timestep routines.
- If a transform uses a new slab count or layout, add it to the plan inventory
  and document which phase uses it.

## Derivatives Contract

`derivatives_gpu.f90` reads device-resident fields and writes device-resident
derivatives.  It should not make host data authoritative during the regular
timestep.

Allowed host-visible boundaries:

- initialization;
- validation or debug modes;
- output/diagnostics that explicitly wait and reduce or copy data.

Changing derivative kernels requires at least a divergence/kinetic-energy
sanity check and comparison against the immediate previous accepted branch.

## Convection Contract

`convec_gpu.f90` computes the dealiased convective RHS and intentionally allows
queued GPU work to overlap with ATM phase 1 host work.

Rules:

- Preserve the five-phase algorithm and boundary handling unless a numerical
  validation record proves the new order.
- Keep padded big-array scratch persistent.
- Do not add a terminal wait unless the lost ATM overlap is measured and
  justified.

## SGS Contract

SGS ownership is split across:

- `sgs_stag_util.f90`: dispatch, runtime model selection, CPU/GPU routing, and
  shared SGS helpers.
- `sgs_gpu.f90`: GPU stress/div-stress kernels and supported runtime-model GPU
  paths.
- `lagrange_Sdep_gpu.f90`: batched Lagrangian dynamic SGS update.

Rules:

- The disabled SGS path must keep working.
- Runtime values `1..5` must be covered.  Do not optimize only `sgs=5`.
- Unsupported values should fail clearly.
- Wall/IWM fields such as `iwm_flt_us` are `(i,j)` surface fields, not diagonal
  accumulators.

Changing SGS dispatch requires at least:

- SGS disabled smoke;
- `sgs=1..5` smoke or an explicit reason a value is not compiled;
- divergence and kinetic-energy comparison;
- regular-step and SGS-update-step timing separation.

## Pressure Contract

Pressure owns one of the most expensive synchronization paths.  Keep
GPU-aware MPI and host-staged behavior distinct in code comments and benchmark
records.

Rules:

- `press_gpu.f90` owns pressure RHS assembly and pressure-gradient output for
  the GPU path.
- `tridag_gpu.f90` owns the GPU/MPI tridiagonal pipeline.
- `tridag_array.f90` remains the CPU/fallback orchestration point.
- `PPGPU_AWARE_MPI` means device pointers are passed through
  `host_data use_device`.
- Host-staged fallbacks must remain correct, but should not be compared as the
  same performance path.
- The replicated global tridiagonal solve is retained as a validation/benchmark
  variant, not as the production default.

Pressure changes require divergence checks and timing against the immediate
previous accepted branch.

## Scalar Contract

`scalars.f90` owns scalar CPU/GPU transport and scalar halo handling.  Scalar
changes must keep all three modes understandable:

- scalar module off;
- scalar module on with CPU path;
- scalar module on with GPU path.

If a scalar change touches CPS or stability coupling, validate those paths
together instead of testing scalar transport alone.

## CPS And Inflow Contract

`concurrent_precursor.f90` owns CPS exchange and timing.  HIT and shifted inflow
paths are optional and should remain isolated from the production CPS path.

Rules:

- Do not make CPS-required behavior depend on HIT-only assumptions.
- Stage-timing switches should remain diagnostic only.
- Any communication optimization should state whether it changes message count,
  payload size, or synchronization order.

HIT inflow is enabled by `USE_HIT=ON` and selected at runtime with
`inflow_type = 2`.  `hit_inflow.f90` owns file input, plane position, restart
state, and the CPU interpolation fallback.  `hit_inflow_gpu.f90` owns the
OpenACC/CUDA-resident HIT field, interpolation plane, and fringe blending.

Rules for HIT changes:

- Keep HIT velocity input and `restartHIT.dat` host-owned; this is setup and
  checkpoint I/O, not a per-step hot-path transfer.
- Keep the GPU per-step work in `hit_inflow_gpu.f90`, not duplicated inside
  `hit_inflow.f90`.
- Validate both the configure path and a small HIT runtime case before using
  HIT timing as production evidence.

## ATM And Structural Coupling Contract

ATM is intentionally mixed CPU/GPU code.  Keep boundaries explicit.

`atm_lesgo_interface.f90` owns:

- LESGO/ATM data exchange;
- turbine velocity sampling;
- force-field application;
- phase 1 / phase 2 split;
- explicit optional host/device compatibility boundaries.

`actuator_turbine_model.f90` owns:

- turbine model state;
- structural-solver controls;
- induced-velocity correction;
- turbine output writing.

Rules:

- Keep induced-velocity correction enabled in the canonical branch.
- Pair structure-off and structure-on tests for changes that touch ATM.
- Keep the standard atPoint sampler and force deposition device-resident.
- Spalart and nacelle host bridges must transfer only the state they own.
- With `updateInterval > 1`, restart checkpoints must be aligned to a multiple
  of the interval until a versioned held-force sidecar is implemented.
- Turbine power files remain the authoritative power output; stdout diagnostics
  are secondary.

## Diagnostics And Output Contract

Diagnostics are allowed to synchronize or reduce data, but their timing must be
identified.  Do not use a diagnostic/output timestep as the sole performance
representative for regular timesteps.

Rules:

- `energy` and `rmsdiv` should reduce on device where supported.
- `output_loop` and restart/checkpoint writes are allowed host-visible
  boundaries.
- Benchmark notes must identify output/statistics cadence.

## Validation Matrix By Change Area

| Change area | Minimum local checks | Minimum HPC validation |
| --- | --- | --- |
| Docs, comments, scripts | `python3 tools/check_branch_readiness.py` | None unless run scripts changed. |
| CMake options | readiness check plus CMake configure smoke | One configure on intended HPC compiler stack. |
| Main timestep order | readiness check | Structure-off smoke, structure-on smoke if ATM compiled, timing comparison. |
| Derivatives / convection | readiness check | Divergence, kinetic energy, regular-step timing comparison. |
| SGS | readiness check | SGS disabled and `sgs=1..5` smoke where compiled. |
| Pressure / tridiagonal | readiness check | GPU-aware and fallback awareness, divergence, pressure timing comparison. |
| Scalars | readiness check | scalar off/on, GPU scalar path when enabled. |
| CPS / inflow | readiness check | CPS smoke and stage-timing comparison; HIT configure/runtime smoke when `USE_HIT=ON`. |
| ATM / structural solver | readiness check | structure off/on, turbine power output, force/velocity sanity checks. |
| LVLSET | readiness check | dedicated LVLSET validation; not covered by production non-LVLSET smoke. |
