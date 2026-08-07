# LESGO Source File Inventory

This inventory lists every tracked Fortran source file in the published
optimized non-LVLSET branch.  It is intentionally broader than the production
GPU path: optional LVLSET, HIT, scalar, turbine, ATM, and validation helper
files are included so collaborators can see what exists and what is currently
inside or outside the optimized path.

The production ownership rules are defined in `docs/gpu_module_contracts.md`.
Use this file as a source-tree index, not as a validation claim.
Root `CMakeLists.txt` groups common and optional source files with matching
`*_SOURCES` variables before feature options append them into the build.

## Driver, Parameters, And Shared Utilities

| File | Build path | Primary role |
| --- | --- | --- |
| `main.f90` | common | Timestep driver and top-level solver orchestration. |
| `param.f90` | common | Runtime parameters and namelist-style configuration state. |
| `sim_param.f90` | common | Core LES field storage and shared simulation arrays. |
| `sgs_param.f90` | common | SGS model parameters and shared SGS state. |
| `types.f90` | common | Precision and common type definitions. |
| `initialize.f90` | common | Solver setup. |
| `initial.f90` | common | Initial condition setup. |
| `finalize.f90` | common | Solver shutdown and cleanup. |
| `input_util.f90` | common | Text input parsing utilities. |
| `param_output.f90` | common | Human-readable parameter output. |
| `string_util.f90` | common | String helpers. |
| `messages.f90` | common | Message formatting and output helpers. |
| `functions.f90` | common | Numerical helper functions. |
| `grid.f90` | common | Grid construction and grid metadata. |
| `cfl_util.f90` | common | CFL and timestep utility routines. |
| `clocks.f90` | common | Lightweight timing helper. |
| `init_random_seed.f90` | common | Random seed initialization. |
| `pid.f90` | common | Controller/PID helper routines. |
| `linear_simple.f90` | `USE_LVLSET` | Linear-system helper used by LVLSET support. |

## LES Hot Path

| File | Build path | Primary role |
| --- | --- | --- |
| `convec.f90` | common | CPU convection path plus CUDA fallback utilities. |
| `convec_gpu.f90` | `USE_LES_GPU` | Optimized GPU convection path. |
| `derivatives.f90` | common | CPU derivative path plus CUDA FFT derivative utilities. |
| `derivatives_gpu.f90` | `USE_LES_GPU` | Optimized GPU derivative kernels. |
| `divstress_uv.f90` | common | Horizontal velocity stress-divergence assembly. |
| `divstress_w.f90` | common | Vertical velocity stress-divergence assembly. |
| `fft.f90` | common | CPU FFT packing, unpacking, and wavenumber setup. |
| `fft_gpu.f90` | `USE_LES_GPU` | CUFFT setup and GPU FFT execution wrappers. |
| `press_stag_array.f90` | common | Pressure solve orchestration and CPU/GPU pressure utilities. |
| `press_gpu.f90` | `USE_LES_GPU` | Optimized GPU pressure path. |
| `tridag_array.f90` | common | Tridiagonal pressure-solver implementation and MPI pipeline. |
| `tridag_gpu.f90` | `USE_LES_GPU` | GPU tridiagonal support routines. |
| `sgs_stag_util.f90` | common | SGS runtime model dispatch and stress assembly. |
| `sgs_gpu.f90` | `USE_LES_GPU` | GPU SGS helper kernels and runtime paths. |
| `std_dynamic.f90` | common | Standard dynamic SGS path. |
| `scaledep_dynamic.f90` | common | Scale-dependent dynamic SGS path. |
| `lagrange_Sdep.f90` | common | CPU scale-dependent Lagrangian SGS path. |
| `lagrange_Sdep_gpu.f90` | `USE_LES_GPU` | Batched GPU Lagrangian SGS update path. |
| `lagrange_Ssim.f90` | common | Similarity-model Lagrangian SGS path. |
| `interpolag_Sdep.f90` | common | Interpolation helpers for scale-dependent Lagrangian SGS. |
| `interpolag_Ssim.f90` | common | Interpolation helpers for similarity Lagrangian SGS. |
| `test_filtermodule.f90` | common | Test-filter implementation used by SGS and wall-model paths. |
| `wallstress.f90` | common | Wall-stress model dispatch. |
| `iwmles.f90` | common | Integral wall-model implementation. |
| `emul_complex.f90` | common | Complex-number helper operations. |

## Forcing, Boundary Conditions, And Inflow

| File | Build path | Primary role |
| --- | --- | --- |
| `forcing.f90` | common | Forcing, projection, turbine-force application, and velocity halo sync. |
| `coriolis.f90` | common | Coriolis forcing. |
| `fringe.f90` | common | Fringe-region forcing. |
| `sponge.f90` | common | Sponge-zone forcing. |
| `inflow.f90` | common | Standard inflow handling. |
| `shifted_inflow.f90` | common | Shifted inflow handling. |
| `concurrent_precursor.f90` | `USE_CPS` | Concurrent precursor simulation coupling. |
| `hit_inflow.f90` | `USE_HIT` | Homogeneous isotropic turbulence inflow. |
| `hit_inflow_gpu.f90` | `USE_HIT` | GPU HIT inflow support. |

## Turbines, ATM, And Structural Coupling

| File | Build path | Primary role |
| --- | --- | --- |
| `turbines.f90` | `USE_TURBINES` | Actuator disk turbine model. |
| `turbines_gpu.f90` | `USE_TURBINES` | GPU actuator disk turbine support. |
| `turbine_indicator.f90` | `USE_TURBINES` | Turbine indicator functions. |
| `atm_base.f90` | `USE_ATM` | ATM math and utility routines. |
| `atm_input_util.f90` | `USE_ATM` | ATM, turbine, airfoil, and structural input parsing. |
| `actuator_turbine_model.f90` | `USE_ATM` | Actuator turbine model, structural solver, and turbine output. |
| `atm_lesgo_interface.f90` | `USE_ATM` | LESGO/ATM data exchange, GPU sampling, force coupling, and load balance. |

## Scalars And Stability

| File | Build path | Primary role |
| --- | --- | --- |
| `scalars.f90` | `USE_SCALARS` | Scalar transport, including optional GPU scalar path. |
| `stability.f90` | `USE_SCALARS` | Stability coupling for scalar-enabled cases. |

## MPI And Diagnostics

| File | Build path | Primary role |
| --- | --- | --- |
| `mpi_defs.f90` | `USE_MPI` | MPI domain metadata and shared MPI definitions. |
| `mpi_transpose_mod.f90` | `USE_MPI` | MPI transpose helpers. |
| `cuda_mpi_debug.f90` | `USE_MPI` | CUDA-aware MPI diagnostics and wrappers. |
| `rmsdiv.f90` | common | Divergence diagnostic. |
| `io.f90` | common | Checkpoint, restart, and field output. |
| `time_average.f90` | common | Time-averaging support. |
| `stat_defs.f90` | common | Statistics definitions and storage. |

## LVLSET Files

Level Set is an optional profile with a GPU timestep path. Geometry generation
and tree preprocessing remain host-side startup work.

| File | Build path | Primary role |
| --- | --- | --- |
| `level_set.f90` | `USE_LVLSET` | Level-set model implementation. |
| `level_set_base.f90` | `USE_LVLSET` | Level-set base data structures and helpers. |
| `level_set_gpu.f90` | `USE_LVLSET_GPU` | OpenACC Level Set interpolation, stress, forcing, smoothing, and SGS coupling. |
| `trees_base_ls.f90` | `USE_LVLSET` | Tree data structures for LVLSET. |
| `trees_global_fmask_ls.f90` | `USE_LVLSET` | Global fluid-mask tree support. |
| `trees_io_ls.f90` | `USE_LVLSET` | LVLSET tree I/O. |
| `trees_pre_ls.f90` | `USE_LVLSET` | LVLSET tree preprocessing. |
| `trees_setup_ls.f90` | `USE_LVLSET` | LVLSET tree setup. |

## Validation Helpers

| File | Build path | Primary role |
| --- | --- | --- |
| `tools/validate_filt_da_cufft.F90` | helper only | Standalone validation helper for CUFFT derivative behavior. |

## Maintenance Rule

Any tracked Fortran source file should appear exactly once as a table row in
this inventory.  The `Build path` column should use `common`, `helper only`, or
tracked root CMake `USE_*` options.  The readiness gate enforces those rules
with:

```bash
python3 tools/check_source_inventory.py
python3 tools/check_cmake_source_groups.py
```
