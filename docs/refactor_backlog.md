# Refactor Backlog

This backlog is for readability and maintainability only.  Do not combine these
items with numerical or performance changes unless the validation plan is clear.

## Current Size Hotspots

The largest files in the current source tree are the first places to improve
structure, but they should be split only when a clean ownership boundary exists.

| Priority | File | Current size | Refactor direction |
| --- | --- | ---: | --- |
| 1 | `actuator_turbine_model.f90` | 134,875 bytes | Separate structural-solver controls, induced-velocity method selection, diagnostics, restart state, and turbine output writing. |
| 2 | `atm_lesgo_interface.f90` | 101,248 bytes | Separate timing, optional sampling bridges, gather helpers, and LESGO force application after the active ownership boundary remains stable. |
| 3 | `lagrange_Sdep_gpu.f90` | 94,070 bytes | Isolate batched Lagrangian update kernels from setup/validation code. |
| 4 | `level_set_gpu.f90` | 72,757 bytes | Keep workspace ownership, stress treatment, smoothing, SGS, and packed-halo sections navigable; split only with the dedicated Level Set matrix in place. |
| 5 | `io.f90` | 64,458 bytes | Separate checkpoint, instantaneous output, and restart metadata helpers. |
| 6 | `scalars.f90` | 63,405 bytes | Separate scalar transport, halo handling, and timing diagnostics. |
| 7 | `atm_input_util.f90` | 57,365 bytes | Separate turbine/airfoil parsing from validation and defaulting logic. |
| 8 | `sgs_gpu.f90` | 56,080 bytes | Separate SGS GPU kernels by tensor assembly, wall stress, and model dispatch. |

`level_set.f90`, `level_set_gpu.f90`, and `trees_*_ls.f90` use a separate
validated optional profile. Further file splitting is lower priority than
preserving the completed CPU/GPU matrix and startup/runtime ownership boundary.
The readiness gate checks this table against the current tracked source tree so
the priority list does not silently drift as files are refactored.

## Preprocessor Cleanup

The current branch mixes several layers of build-time conditionals:

- production `PPLES_GPU`;
- module flags such as `PPSGS_GPU`, `PPSCALARS_GPU`, `PPATM`, and `PPCPS`;
- communication flags such as `PPMPI` and `PPGPU_AWARE_MPI`.

Refactor rule:

- Keep top-level build choices in `CMakeLists.txt`.
- Keep each module's CPU/GPU dispatch close to the module entry point.
- Avoid burying large algorithmic branches deep inside nested loops.
- When a fallback path is kept only for validation, label it as validation or
  fallback code in the source comment.

## Environment Switch Cleanup

Several environment variables remain because they control diagnostics,
validation, or benchmark-only timing.  Future cleanup should classify each one:

- `production`: behavior required by real runs;
- `diagnostic`: prints timing, validation, or environment details;
- `benchmark`: used only to reproduce historical timing studies;
- `remove`: no longer supported or no clear reason to keep.

After classification, production behavior should move to CMake options or
documented input parameters where possible.  Diagnostic and benchmark switches
may remain as environment variables if they do not affect numerical results.

## Next Safe Refactor Order

1. Keep the simplified `atm_lesgo_interface.f90` ownership boundary stable;
   split optional sampling bridges only with paired CPU/GPU validation.
2. Extract SGS runtime-model dispatch documentation from `sgs_stag_util.f90`.
3. Keep `docs/environment_switches.md` synchronized with source changes.
4. Add a small CMake option matrix smoke script result to the validation record.
5. Only after those steps, consider moving helper routines into new modules.

## Required Checks For Each Refactor

- `git diff --check`
- CMake configure on the intended HPC build environment.
- One structure-off and one structure-on smoke when ATM code changes.
- SGS disabled plus `sgs=1..5` smoke when SGS dispatch changes.
- Scalar off/on smoke when `scalars.f90` or CPS scalar handling changes.
- Timing comparison against the immediate previous accepted branch.
