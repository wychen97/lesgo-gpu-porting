# Environment Switches

This branch still contains environment-variable controls from optimization,
diagnostics, and validation work.  This file classifies them so collaborators
know which switches are safe to use and which should eventually move to a
documented input or CMake option.

## Classification

- `production/path-control`: changes an active execution path.
- `physics/model`: changes numerical or physical model behavior.
- `diagnostic`: prints timing, validation, synchronization, or environment
  information without intentionally changing the result.
- `benchmark/experimental`: retained to reproduce a test or compare an
  alternate implementation; not part of the default production path.
- `reproducibility`: controls deterministic setup.
- `read-only system`: queried only for logs/audits.

## Production And Model Controls

In `atm_lesgo_interface.f90`, simple boolean ATM gates are parsed through the
shared `atm_env_flag()` helper and then cached by the individual gate function.
In `actuator_turbine_model.f90`, structure and diagnostic gates that previously
accepted only explicit true tokens are parsed through the private
`atm_model_env_enabled()` helper. Physics/model switches in the same file use
`atm_model_env_token()` so the token normalization is shared while each switch
keeps its own default and accepted values.
In `sgs_stag_util.f90`, diagnostic gates use `sgs_env_true_token_enabled()`;
path-selection and benchmark gates that historically treated any non-false set
value as enabled use `sgs_env_enabled_unless_false()`.
In `press_stag_array.f90`, pressure switches that use the same non-false
set-value convention are parsed through
`press_apply_env_enabled_unless_false()`.
In `tridag_array.f90`, transpose timing uses
`tridag_apply_env_enabled_unless_false()`, while benchmark variants that
require explicit true tokens use `tridag_apply_env_true_token()`.
In `tridag_gpu.f90`, `LESGO_TRIDAG_NCHUNK` keeps its positive-integer
override convention and `LESGO_TRIDAG_GPU_MPI` is disabled only by the exact
value `0`.
In `scalars.f90` and `concurrent_precursor.f90`, stage-timing diagnostics use
private true-token helpers so unset or unrecognized values remain disabled.
In `forcing.f90`, project and LVLSET bridge timing switches use the same
private true-token helper convention.
In `lagrange_Sdep.f90`, `LESGO_LAGRANGE_STRICT_SYNC` intentionally keeps its
numeric convention: unset, `0`, or invalid values disable strict sync, and any
valid nonzero integer enables it.
In `cuda_mpi_debug.f90`, CUDA/MPI diagnostic switches keep their legacy
convention: unset, empty, or exactly `0` disables the diagnostic; any other set
value enables it.
In `main.f90`, benchmark-only CPU reference timings are read through
`main_read_env_real()`; unset or invalid values fall back to `0.0`.

| Variable | Class | Source | Effect |
| --- | --- | --- | --- |
| `LESGO_ATM_STRUCTURE` | production/path-control | `actuator_turbine_model.f90` | Enables structural solver coupling. |
| `LESGO_ATM_STRUCTURE_GPU` | production/path-control | `actuator_turbine_model.f90` | Enables GPU structural-solver path when compiled. |
| `LESGO_ATM_STRUCTURE_GPU_DIRECT` | production/path-control | `actuator_turbine_model.f90` | Selects direct GPU structural data path. |
| `LESGO_ATM_FORCE_SHADOWS` | production/path-control | `atm_lesgo_interface.f90` | Enables/disables turbine force shadow buffers. |
| `LESGO_SGS_HALO_COMBINED` | production/path-control | `sgs_stag_util.f90` | Selects combined SGS tau halo path for supported rank layouts. |
| `LESGO_SGS_CALCSIJ_EXPLICIT` | production/path-control | `sgs_stag_util.f90` | Selects explicit calc_Sij kernel path. |
| `LESGO_SGS_EXPLICIT_POINTWISE` | production/path-control | `sgs_stag_util.f90` | Selects explicit pointwise SGS path. |
| `LESGO_PRESS_RHS_HALO_COMBINED` | production/path-control | `press_stag_array.f90` | Selects combined pressure RHS halo path for `nproc=2`. |
| `LESGO_TRIDAG_GPU_MPI` | production/path-control | `tridag_gpu.f90` | Allows disabling GPU/MPI tridiagonal pipeline. |
| `LESGO_TRIDAG_NCHUNK` | production/path-control | `tridag_gpu.f90` | Sets chunk count for tridiagonal GPU/MPI pipeline. |
| `LESGO_ATM_INDUCED_METHOD` | physics/model | `actuator_turbine_model.f90` | Selects induced-velocity correction method. Default is the 2024 exact panel path. |
| `LESGO_ATM_DU_INCLUDE_UX` | physics/model | `actuator_turbine_model.f90` | Includes/excludes axial induced velocity in `du`. |
| `LESGO_ATM_UX_USE_CD_EFF` | physics/model | `actuator_turbine_model.f90` | Chooses drag-effective axial induced velocity option. |
| `LESGO_ATM_STRUCTURE_VEL_FEEDBACK` | physics/model | `actuator_turbine_model.f90` | Allows disabling structural velocity feedback. |
| `LESGO_ATM_STRUCTURE_ALPHA_FEEDBACK` | physics/model | `actuator_turbine_model.f90` | Allows disabling structural angle-of-attack feedback. |
| `LESGO_RANDOM_SEED` | reproducibility | `init_random_seed.f90` | Overrides random seed initialization. |

## Diagnostic Switches

| Variable | Source | Effect |
| --- | --- | --- |
| `LESGO_PROJECT_STAGE_TIMING` | `forcing.f90` | Enables top-level stage timing prints. |
| `LESGO_CPS_STAGE_TIMING` | `concurrent_precursor.f90` | Enables concurrent precursor stage timing. |
| `LESGO_SCALAR_STAGE_TIMING` | `scalars.f90` | Enables scalar transport stage timing. |
| `LESGO_SGS_STAGE_TIMING` | `sgs_stag_util.f90` | Enables SGS stage timing. |
| `LESGO_SGS_STRICT_SYNC` | `sgs_stag_util.f90` | Adds strict synchronization for SGS kernel debugging. |
| `LESGO_LAGRANGE_STRICT_SYNC` | `lagrange_Sdep.f90` | Adds strict synchronization around Lagrangian SGS update debugging. |
| `LESGO_PRESS_STAGE_TIMING` | `press_stag_array.f90` | Enables pressure stage timing. |
| `LESGO_PRESS_TRANSPOSE_TIMING` | `tridag_array.f90` | Enables transpose timing in pressure/tridiagonal path. |
| `LESGO_LVLSET_BRIDGE_TIMING` | `forcing.f90` | Enables LVLSET bridge timing; LVLSET is not production in this branch. |
| `LESGO_ATM_DIAG_TIMING` | `atm_lesgo_interface.f90` | Enables ATM diagnostic timing. |
| `LESGO_ATM_STRUCTURE_TIMING` | `actuator_turbine_model.f90` | Enables structural solver timing report. |
| `LESGO_ATM_STRUCTURE_DIAG` | `actuator_turbine_model.f90` | Enables structural diagnostics. |
| `LESGO_ATM_POWER_STDOUT` | `actuator_turbine_model.f90` | Prints turbine power to stdout; turbine power files are still the authoritative output. |
| `LESGO_ATM_STRUCTURE_GPU_VALIDATE` | `actuator_turbine_model.f90` | Enables validation for GPU structural path. |
| `LESGO_MPI_CUDA_DEBUG` | `cuda_mpi_debug.f90` | Prints CUDA/MPI environment diagnostics. |
| `LESGO_MPI_CUDA_SYNC` | `cuda_mpi_debug.f90` | Adds CUDA synchronization for MPI/CUDA debugging. |

## Benchmark Or Experimental Switches

| Variable | Source | Effect |
| --- | --- | --- |
| `LESGO_SGS_CALCSIJ_DEVICE_BENCH` | `sgs_stag_util.f90` | Enables calc_Sij device lower-bound benchmark path. |
| `LESGO_PRESS_TRANSPOSE_GENERIC` | `tridag_array.f90` | Enables generic transpose variant for comparison. |
| `LESGO_PRESS_DIRECT_THOMAS_OUT` | `tridag_array.f90` | Enables direct Thomas-output variant for comparison. |
| `LESGO_ATM_POINT_OWNER_LB` | `atm_lesgo_interface.f90` | Enables experimental point-owner load balancing. |
| `LESGO_ATM_POINT_OWNER_TARGETED` | `atm_lesgo_interface.f90` | Enables targeted point-owner load balancing variant. |
| `LESGO_ATM_LB_AUTO_SELECT` | `atm_lesgo_interface.f90` | Alternates legacy/LB calls during a short probe window and selects the faster path. |
| `LESGO_ATM_LB_VALIDATE` | `atm_lesgo_interface.f90` | Enables ATM load-balance validation output. |
| `LESGO_CPU_REF_TIME_TOTAL` | `main.f90` | Supplies historical CPU reference time for speedup reporting only. |
| `LESGO_CPU_REF_TIME_FORCING` | `main.f90` | Supplies historical CPU forcing reference time for speedup reporting only. |

## Read-Only System Variables

| Variable | Source | Effect |
| --- | --- | --- |
| `CUDA_VISIBLE_DEVICES` | pressure, SGS, tridiagonal diagnostics | Logged during pointer/MPI audits. |
| `MPICH_GPU_SUPPORT_ENABLED` | pressure, SGS, tridiagonal diagnostics | Logged during pointer/MPI audits. |

## Cleanup Rules

1. Do not add a new environment switch without adding it to this file.
2. Production/path-control switches should eventually become CMake options or
   documented input-file parameters.
3. Diagnostic switches may remain environment variables if they do not change
   numerical behavior.
4. Benchmark/experimental switches should be removed once the historical test is
   no longer needed.
5. Physics/model switches must be reported in any benchmark or validation note,
   because they can change the simulated turbine response.

Run this before committing any source change that adds, removes, or renames a
solver `LESGO_*` switch.  The check requires every source switch to appear in
one classified table row, not only in prose:

```bash
python3 tools/check_environment_switch_docs.py
```
