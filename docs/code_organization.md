# LESGO Code Organization

This branch contains the optimized LES GPU core and an opt-in GPU Level Set
path intended for shared development. The solver is still mostly organized as
flat Fortran modules, so this document is the source map collaborators should
read before editing code.
See `docs/gpu_module_contracts.md` for per-module ownership and validation
contracts.

## Production Configuration

The validated production path is:

```text
USE_MPI=ON
USE_CPS=ON
USE_ATM=ON
USE_TURBINES=ON
USE_LES_GPU=ON
USE_GPU_AWARE_MPI=AUTO
USE_SCALARS=ON
USE_SCALARS_GPU=ON
USE_LVLSET=OFF
USE_LVLSET_GPU=OFF
USE_HIT=OFF
USE_DYN_TN=OFF
USE_CGNS=OFF
```

The canonical turbine production profile keeps Level Set disabled because it
is different physics. Level Set GPU builds are supported and have a dedicated
validation profile; their evidence is recorded separately from turbine cases.

Use `docs/build_profiles.md` for copy-paste CMake configure commands and
user-facing cache variables such as `hostname`, `WRITE_ENDIAN`, and
`READ_ENDIAN`.

## Build Options

Every root CMake `USE_*` option should appear exactly once in this table.  The
readiness gate checks this table against `CMakeLists.txt`.

| Option | Preprocessor flag | Main files added or affected | Notes |
| --- | --- | --- | --- |
| `USE_MPI` | `PPMPI` | `mpi_defs.f90`, `mpi_transpose_mod.f90` | Required for production runs. |
| `USE_CPS` | `PPCPS` | `concurrent_precursor.f90` | Concurrent precursor path. |
| `USE_HIT` | `PPHIT` | `hit_inflow.f90`, `hit_inflow_gpu.f90` | Optional HIT inflow path. |
| `USE_LVLSET` | `PPLVLSET` | `level_set.f90`, `level_set_base.f90`, `trees_*_ls.f90` | Enables immersed-surface Level Set physics. |
| `USE_LVLSET_GPU` | `PPLVLSET_GPU` | `level_set_gpu.f90`, Level Set dispatch in SGS/forcing | GPU Level Set hot path. Requires `USE_LVLSET=ON` and `USE_LES_GPU=ON`. |
| `USE_TURBINES` | `PPTURBINES` | `turbines.f90`, `turbines_gpu.f90`, `turbine_indicator.f90` | Actuator disk style turbine support. |
| `USE_ATM` | `PPATM` | `atm_base.f90`, `atm_input_util.f90`, `actuator_turbine_model.f90`, `atm_lesgo_interface.f90` | Actuator turbine model / structural solver interface. |
| `USE_LES_GPU` | `PPLES_GPU`, `PPCONVEC_GPU`, `PPDERIVS_GPU`, `PPPRESS_GPU`, `PPSGS_GPU` | `*_gpu.f90`, pressure, SGS, derivatives, convection | Optimized OpenACC/CUDA LES core. |
| `USE_CPU_BUILD` | none | CMake compiler/link flags | NVHPC CPU baseline mode: accepts dormant CUDA Fortran syntax without linking GPU libraries. |
| `USE_SCALARS` | `PPSCALARS` | `scalars.f90`, `stability.f90` | Scalar transport, e.g. temperature/passive scalar. |
| `USE_SCALARS_GPU` | `PPSCALARS_GPU` | `scalars.f90` | GPU scalar transport path. Requires `USE_LES_GPU=ON`. |
| `USE_GPU_AWARE_MPI` | `PPGPU_AWARE_MPI` | halo exchange and tridiagonal pressure paths | `AUTO` links Cray CUDA GTL when available. |
| `USE_DYN_TN` | `PPDYN_TN` | Lagrangian SGS timescale update | Optional runtime model path. |
| `USE_SAFETYMODE` | `PPSAFETYMODE` | `main.f90`, diagnostic guards | Extra runtime safety checks. Enabled by default. |
| `USE_CGNS` | `PPCGNS` | `io.f90` and CGNS dependencies | Optional output dependency. |

## Module Map

### Driver and global state

- `main.f90`: timestep driver and module orchestration.
- `param.f90`, `sim_param.f90`, `sgs_param.f90`, `types.f90`: shared
  parameters, field storage, and SGS state.
- `initialize.f90`, `initial.f90`, `finalize.f90`: setup and shutdown.
- `input_util.f90`, `param_output.f90`, `string_util.f90`: configuration and
  user-facing parameter output.

### LES hot path

- `convec.f90`, `convec_gpu.f90`: nonlinear advection / convection.
- `derivatives.f90`, `derivatives_gpu.f90`: derivative kernels.
- `divstress_uv.f90`, `divstress_w.f90`: divergence of stress terms.
- `press_stag_array.f90`, `press_gpu.f90`, `tridag_array.f90`,
  `tridag_gpu.f90`: pressure solve and tridiagonal solver.
- `fft.f90`, `fft_gpu.f90`, `test_filtermodule.f90`: spectral transforms and
  test filters.
- `sgs_stag_util.f90`, `sgs_gpu.f90`, `std_dynamic.f90`,
  `scaledep_dynamic.f90`: SGS model dispatch and runtime SGS paths.
- `lagrange_Sdep.f90`, `lagrange_Sdep_gpu.f90`, `lagrange_Ssim.f90`:
  Lagrangian SGS averaging and batched GPU updates.
- `wallstress.f90`, `iwmles.f90`: wall stress and integral wall model.

### Forcing, inflow, and optional physics

- `forcing.f90`, `coriolis.f90`, `fringe.f90`, `sponge.f90`: external forcing
  and boundary-zone terms.
- `inflow.f90`, `shifted_inflow.f90`, `hit_inflow*.f90`,
  `concurrent_precursor.f90`: inflow sources and concurrent precursor support.
- `scalars.f90`, `stability.f90`: scalar transport and stability coupling.
- `level_set.f90`, `level_set_base.f90`, `level_set_gpu.f90`: Level Set CPU
  reference, shared state, and GPU kernels.
- `trees_*_ls.f90`: host-side Level Set geometry preprocessing at startup.

### Turbine and ATM modules

- `turbines.f90`, `turbines_gpu.f90`, `turbine_indicator.f90`: actuator disk
  support.
- `atm_base.f90`, `atm_input_util.f90`, `actuator_turbine_model.f90`,
  `atm_lesgo_interface.f90`: actuator turbine model, structural solver
  coupling, induced-velocity correction, and LESGO/ATM data exchange.

### I/O, diagnostics, and utilities

- `io.f90`, `time_average.f90`, `stat_defs.f90`: output, restart, statistics,
  and averaging.
- `rmsdiv.f90`, `clocks.f90`: diagnostics and timing.
- `cuda_mpi_debug.f90`: MPI-only CUDA/MPI diagnostics and debug wrappers.
- `grid.f90`, `functions.f90`, `linear_simple.f90`, `pid.f90`: numerical and
  utility support.

## Production Timestep Flow

The non-LVLSET GPU production path is orchestrated by `main.f90`.  When editing
the hot path, use this order as the starting point:

1. Derivatives:
   `filt_da_gpu`, `ddz_uv_gpu`, and `ddz_w_gpu` compute velocity gradients from
   resident `u`, `v`, and `w`.
2. Scalar derivatives and scalar transport:
   `scalars_deriv` and `scalars_transport` run only when scalar transport is
   compiled and enabled.
3. Wall stress and SGS:
   wall stress updates run before `sgs_stag_gpu`; SGS writes stress tensors and
   div-stress fields used by the RHS assembly.
4. Convection:
   `convec_gpu` computes the dealiased convective RHS.  It intentionally returns
   with async queue 1 still draining so ATM phase 1 can overlap host work.
5. ATM phase 1:
   in the PPLES GPU ATM path, `atm_lesgo_forcing(phase=1)` samples velocities,
   updates blade/structure state, and gathers turbine-side data while the
   convection backlog drains.
6. RHS assembly:
   `main.f90` combines convection and div-stress terms directly on the device.
7. Forcing and ATM phase 2:
   the original `forcing_applied` site remains the place where turbine forces
   are convolved into `fxa/fya/fza`, corrected, applied, and written.
8. Pressure:
   `press_stag_array_gpu` builds pressure RHS terms, runs the tridiagonal
   pressure path, and returns pressure gradients on the device.
9. Pressure-gradient RHS update:
   `main.f90` subtracts `dpdx/dpdy/dpdz` from `RHSx/RHSy/RHSz` on the device.
10. Projection:
    `forcing_induced` and `project` update `u`, `v`, and `w` for the next
    timestep.
11. Diagnostics and output:
    `energy`, `output_loop`, `rmsdiv`, and final output are the main
    host-visible consumers.  Periodic diagnostics can make a printed timestep
    slower than a regular timestep.

Do not move work across these boundaries unless the data dependencies and
periodic diagnostic side effects are checked explicitly.

## Device Ownership Boundaries

The optimized branch relies on explicit OpenACC device residency for the main
LES arrays.  The practical ownership rules are:

- `sim_param.f90` owns the core field declarations.  In the PPLES GPU path,
  `u/v/w`, gradients, RHS arrays, pressure gradients, stresses, div-stresses,
  and force arrays are expected to have persistent device mirrors.
- GPU hot-path modules should use `present(...)` data regions and should not
  allocate, copy, or free full fields inside per-step routines unless the
  boundary is documented.
- Full-field `update self` or `update device` calls are allowed at setup,
  restart, output, diagnostics, validation, and explicitly documented CPU
  fallback boundaries.  They should not appear silently inside regular
  timestep kernels.
- GPU-aware MPI paths should pass device pointers through
  `host_data use_device`.  Host-staged fallbacks must remain correct, but their
  timings should be reported separately.
- ATM is a mixed CPU/GPU boundary by design: velocity sampling and force
  application use device data, while parts of the blade/structural model still
  run on the host.  Keep phase 1 and phase 2 ordering paired with the comments
  in `main.f90` and `atm_lesgo_interface.f90`.
- The normal atPoint path stays device-resident. Spalart sampling and nacelle
  interpolation use an explicit host compatibility bridge only when requested.
- SGS model dispatch must keep the disabled path and supported runtime values
  `1..5` covered.  Do not optimize only `sgs=5` and leave another runtime value
  on a stale path.
- Level Set geometry is initialized on the host once, then `phi`, normals,
  overlap buffers, and interpolation workspaces remain device-resident when
  `USE_LVLSET_GPU=ON`.
- Level Set MPI communication uses packed velocity/stress/history messages;
  CUDA-aware MPI uses device pointers and the fallback stages only compact
  boundary buffers.

## Current Optimization Status

Use `docs/gpu_port_coverage_audit.md` for the detailed status matrix.  In
short, the validated production non-LVLSET hot path includes:

- LES GPU residency and the GPU-aware MPI pressure pipeline;
- convection, derivatives, pressure, and SGS runtime values `1..5`;
- the SGS disabled path;
- ADM and ATM turbine paths used by the published benchmark cases.

Implemented paths that still need broader production-size validation before
being used as broad speedup evidence:

- scalar transport with `USE_SCALARS_GPU=ON`;
- CPS/concurrent precursor, especially with scalar coupling;
- IWM-heavy wall-model cases;
- optional `USE_DYN_TN=ON` Lagrangian SGS timescale update;
- shifted inflow, sponge, and Coriolis forcing configurations.

Available as an opt-in non-production path:

- HIT inflow with GPU-assisted plane interpolation and fringe application.

Deferred:

- Arbitrary-seam ATM restart with `updateInterval > 1`; aligned restart points
  are enforced until a held-force checkpoint sidecar is implemented.

## Validation Expectations

Any source change that affects the timestep hot path should be checked against
the immediate previous accepted branch, not only against historical timings.
At minimum, record:

- CMake options and compiler modules.
- Grid, MPI rank count, GPU count, turbine count, and SGS model.
- Final divergence and kinetic energy.
- Final printed step timing and cumulative solver time.
- Whether the measured step is a regular step or an output/statistics/SGS
  update step.
