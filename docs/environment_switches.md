# Environment Switches

This file lists solver `LESGO_*` environment-variable controls that are still present in tracked Fortran source. Removed legacy CUDA-Fortran branches also removed their diagnostic/benchmark switches, so this document intentionally reflects only the maintained code paths.

## Classification

- `production/path-control`: changes an active execution path.
- `physics/model`: changes numerical or physical model behavior.
- `diagnostic`: prints timing, validation, synchronization, or environment information without intentionally changing the result.
- `benchmark/experimental`: retained to reproduce a test or compare an implementation; not part of the default production path.
- `reproducibility`: controls deterministic setup.
- `read-only system`: queried only for logs/audits.

## Production And Model Controls

| Variable | Class | Source | Effect |
| --- | --- | --- | --- |
| `LESGO_ATM_STRUCTURE` | production/path-control | `actuator_turbine_model.f90` | Enables structural solver coupling. |
| `LESGO_ATM_STRUCTURE_VEL_FEEDBACK` | physics/model | `actuator_turbine_model.f90` | Allows disabling structural velocity feedback. |
| `LESGO_ATM_STRUCTURE_ALPHA_FEEDBACK` | physics/model | `actuator_turbine_model.f90` | Allows disabling structural angle-of-attack feedback. |
| `LESGO_RANDOM_SEED` | reproducibility | `init_random_seed.f90` | Overrides random seed initialization. |
| `LESGO_TRIDAG_GPU_MPI` | production/path-control | `tridag_gpu.f90` | Allows disabling GPU/MPI tridiagonal pipeline. |
| `LESGO_TRIDAG_NCHUNK` | production/path-control | `tridag_gpu.f90` | Sets chunk count for tridiagonal GPU/MPI pipeline. |

## Diagnostic Switches

| Variable | Source | Effect |
| --- | --- | --- |
| `LESGO_CPS_STAGE_TIMING` | `concurrent_precursor.f90` | Enables concurrent precursor stage timing. |
| `LESGO_SCALAR_STAGE_TIMING` | `scalars.f90` | Enables scalar transport stage timing. |
| `LESGO_LVLSET_BRIDGE_TIMING` | `forcing.f90` | Enables LVLSET bridge timing; LVLSET is not production in this branch. |
| `LESGO_LVLSET_VALIDATION_SNAPSHOT` | `io.f90` | Writes final Level Set pressure-gradient, force, stress-divergence, and optional `Beta` validation sidecars. |
| `LESGO_ATM_STRUCTURE_TIMING` | `actuator_turbine_model.f90` | Enables structural solver timing report. |
| `LESGO_ATM_STRUCTURE_DIAG` | `actuator_turbine_model.f90` | Enables structural diagnostics. |
| `LESGO_ATM_POWER_STDOUT` | `actuator_turbine_model.f90` | Prints turbine power to stdout; turbine power files are still the authoritative output. |
| `LESGO_MPI_CUDA_DEBUG` | `cuda_mpi_debug.f90` | Prints CUDA/MPI environment diagnostics. |
| `LESGO_MPI_CUDA_SYNC` | `cuda_mpi_debug.f90` | Adds CUDA synchronization for MPI/CUDA debugging. |

## Benchmark Or Experimental Switches

| Variable | Source | Effect |
| --- | --- | --- |
| `LESGO_CPU_REF_TIME_TOTAL` | `main.f90` | Supplies historical CPU reference time for speedup reporting only. |
| `LESGO_CPU_REF_TIME_FORCING` | `main.f90` | Supplies historical CPU forcing reference time for speedup reporting only. |

## Read-Only System Variables

| Variable | Source | Effect |
| --- | --- | --- |
| `CUDA_VISIBLE_DEVICES` | pressure, SGS, tridiagonal diagnostics | Logged during pointer/MPI audits. |
| `MPICH_GPU_SUPPORT_ENABLED` | pressure, SGS, tridiagonal diagnostics | Logged during pointer/MPI audits. |

## Cleanup Rules

1. Do not add a new environment switch without adding it to this file.
2. Production/path-control switches should eventually become CMake options or documented input-file parameters.
3. Diagnostic switches may remain environment variables if they do not change numerical behavior.
4. Benchmark/experimental switches should be removed once the historical test is no longer needed.
5. Physics/model switches must be reported in any benchmark or validation note, because they can change the simulated turbine response.

Run this before committing any source change that adds, removes, or renames a solver `LESGO_*` switch. The check requires every source switch to appear in one classified table row, not only in prose:

```bash
python3 tools/check_environment_switch_docs.py
```
