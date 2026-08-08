# Level Set GPU Port

This document defines the implementation and validation boundary for the
optional immersed-surface Level Set GPU path. It does not change the default
turbine production profile.

## Build Contract

Enable all three required switches:

```text
USE_LES_GPU=ON
USE_LVLSET=ON
USE_LVLSET_GPU=ON
```

`USE_LVLSET_GPU` cannot be enabled without the other two options. A CPU
reference uses `USE_LES_GPU=OFF`, `USE_LVLSET_GPU=OFF`, and
`USE_CPU_BUILD=ON`. The combination `USE_LES_GPU=ON`, `USE_LVLSET=ON`, and
`USE_LVLSET_GPU=OFF` remains available as a slow host-bridge regression path.

## Data Ownership

Startup remains host-owned because it reads or generates geometry:

- tree preprocessing and `phi.out` generation;
- Level Set file input;
- normal-vector construction;
- optional `norm.dat` diagnostic output.

Before the timestep loop, `level_set_gpu_data_init()` creates persistent device
copies of geometry, normals, desired-velocity arrays, compact MPI overlap
arrays, and interpolation workspaces. During normal timesteps, velocity,
stress, forcing, SGS history, and Level Set data remain on the device.

The current backend exposes a `level_set_patch_descriptor_t` with level,
generation, spacing, origin, valid/ghost extents, and references to `phi` and
the normal field. Its backend is explicitly `uniform-grid`; it is an ownership
contract for a future AMR adapter, not an AMR implementation. A host geometry
change must call `level_set_geometry_refresh()`, which increments the
generation, invalidates halo and coefficient caches, refreshes the device copy,
and checks the bound extents. `level_set_finalize()` releases device mappings
and host allocations before MPI finalization.

The descriptor reserves cell-weight and coverage-mask references for composite
diagnostics. They remain unbound in the uniform backend. An adapter can bind
both through `level_set_bind_diagnostic_weights()` and clear them through
`level_set_clear_diagnostic_weights()`. The target arrays must outlive the
binding; every geometry refresh clears the references so an adapter cannot
silently reuse metadata from an old generation. When bound, the velocity-error
diagnostic excludes covered coarse cells and normalizes weighted raw error sums
by globally reduced volume weights while retaining raw sample-count reductions.
An AMR backend must provide equivalent composite integration before using
`global_CA`; summing each patch independently would double-count force and area.

Persistent optional storage follows active runtime physics:

- `udes`, `vdes`, and `wdes` exist only when `vel_BC` is enabled;
- the stress snapshot exists only for simple stress extrapolation;
- the `F_MM` snapshot exists only for Lagrangian SGS model 4 or 5.

Startup prints the persistent Level Set GPU allocation, split into geometry,
desired velocity, halo/pack buffers, and workspace memory.

The existing public calls remain the interface:

- `level_set_BC()` selects CPU or GPU stress treatment;
- `level_set_lag_dyn()` selects CPU or GPU model-4 preconditioning;
- `level_set_forcing()` selects CPU or GPU immersed-boundary forcing;
- `level_set_smooth_vel()` selects the matching smoother.

GPU callers therefore do not depend on implementation-specific public entry
points.

## Numerical Compatibility Notes

The Level Set hardening branch contains two explicit model-4/5 corrections
that must be included in CPU-reference comparisons:

- `lag_dyn_modify_beta` is now a real runtime option on both CPU and GPU. The
  previous CPU routine shadowed it with a local `.true.` constant, while the
  GPU routine always modified `beta`.
- With a physical lower wall, CPU `modify_beta` now measures wall distance
  from the global vertical coordinate
  `((coord * (nz - 1)) + k - 1) * dz`, matching the GPU implementation. The
  previous CPU expression restarted the height at zero on every MPI rank.

The coordinate correction leaves single-rank output unchanged. It
intentionally changes multi-rank CPU `beta` above rank zero and therefore
requires old-CPU/new-CPU and new-CPU/GPU comparisons before release.

## GPU Work

The GPU implementation covers:

- trilinear interpolation on velocity, stress, `phi`, and SGS-history fields;
- stress interpolation and simple extrapolation;
- logarithmic and legacy stress-extrapolation options;
- optional local derivative correction (`use_modify_dutdn`);
- desired normal-velocity and logarithmic-profile boundary forcing;
- immersed-boundary force assembly;
- `xy` SOR smoothing and CPU-order-equivalent `3d` wavefront smoothing;
- Smagorinsky wall-distance coefficient treatment;
- Lagrangian model-4/5 smoothing, `F_LM` masking, `F_MM` Neumann treatment,
  and beta modification;
- global force/area and velocity-error reductions.

Tree construction, file I/O, and normal construction are intentionally left on
the host because they execute once at startup and are not timestep hot paths.

Vertical interpolation never clamps a request to the nearest resident plane.
An insufficient local or halo request returns the Level Set sentinel instead
of accessing memory outside the allocation. Validation jobs set
`LESGO_LVLSET_INTERP_BOUNDS_CHECK=ON`; this enables a device error counter and
a synchronized failure at the end of each interpolation-bearing operation.
Production runs omit that diagnostic synchronization while retaining the
memory-safety guard. The startup interpolation self-test includes a deliberate
insufficient-halo request to verify this behavior.

## MPI Communication

The CPU implementation issued separate messages for individual velocity and
stress components. The GPU path packs:

- three velocity components into two directional exchanges;
- six stress components into compact face and overlap exchanges;
- model-4/5 `F_LM` and `F_MM` face state into one exchange.

With `PPGPU_AWARE_MPI`, MPI receives device pointers through
`host_data use_device`. Otherwise only the compact pack buffers are staged;
full three-dimensional fields are never copied for a halo exchange.

The following CPU restrictions are preserved:

- `smooth_mode='3d'` is single-rank only;
- `use_extrap_tau_log=.true.` is single-rank only;
- legacy extrapolation (`use_extrap_tau_simple=.false.`) is single-rank only;
- `use_modify_dutdn=.true.` is single-rank only.

## Validation Results

The release gate uses the checked-in declarative variants in
`test-cases/level_set_cubes/validation_variants.json`, expanded from the base
case onto a `64 x 64 x 32` uniform grid. Strict checkpoint comparisons run for
two timesteps with `DYN_init=1` and `cs_count=2`, so the dynamic-SGS update is
active while the comparison still isolates implementation differences from
long-horizon chaotic amplification. The checked-in base case remains a 20-step
smoke case.

The matrix contains 58 runtime tasks:

- 45 one-rank tasks covering SGS disabled and models 1-5, MPI and non-MPI,
  CPU reference, LES-GPU/CPU-Level-Set bridge, full Level Set GPU, and every
  supported stress and velocity-boundary option, plus centralized invalid-input
  rejection cases;
- 11 two-rank tasks covering CPU, bridge, host-staged GPU, GPU-aware MPI, and
  model-4/5 spherical interfaces crossing the z-rank boundary;
- two four-rank tasks covering host-staged and GPU-aware MPI with one rank per
  A100 and a spherical rank-crossing interface.

The comparator checks `u`, `v`, `w`, RHS, pressure gradients, SGS coefficient
and history fields, `Beta`, Level Set force, divergence, kinetic energy, and
integrated IBM force. It rejects NaN/Inf and sentinel mismatches. A field passes
when all physical cells satisfy the elementwise tolerance or when the
whole-field RMS difference satisfies `atol + rtol * reference_rms`. This avoids
letting a handful of near-zero values dominate a physically negligible field
comparison; a retained 20-step SGS-3 regression still fails the norm gate and
therefore confirms that the criterion remains discriminating.

Derecho passed all 58 runtime tasks and all 51 CPU/bridge/GPU evidence pairs
from the exact source archive at `084f71b605c66a0c50f54eebafed6c0b0dd7239f`.
Build job `7051629.desched1` compiled all ten required CPU, bridge, full-GPU,
MPI/non-MPI, staged/GPU-aware, DYN_TN, and non-Level-Set profiles:

| Ranks / GPUs | Job | Runtime tasks | Evidence pairs | Result |
| ---: | --- | ---: | ---: | --- |
| 1 / 1 A100 | `7051694.desched1` | 45/45 | 38/38 | Passed |
| 2 / 2 A100 | `7051695.desched1` | 11/11 | 12/12 | Passed |
| 4 / 4 A100 | `7051696.desched1` | 2/2 | 1/1 | Passed |

Derecho also passed 18 continuous-versus-split restart comparisons for models
4 and 5 at seams 1, 2, and 3 (`7051697.desched1`). Both model-4 and model-5
`USE_DYN_TN=ON` gates and the non-Level-Set isolation checkpoint passed
(`7051698.desched1`). The compiler stack was NVHPC 25.9, CUDA 12.9.0,
Cray MPICH 8.1.32, and FFTW 3.3.10. The detailed machine-readable record is
`docs/level_set_gpu_validation_evidence.json`.

A follow-up compiler-compatibility smoke used Derecho's newer NVHPC 26.1 with
CUDA 12.9.0, Cray MPICH 8.1.32, FFTW 3.3.10, and CMake 3.31.8. Build job
`7051966.desched1` compiled CPU-MPI, host-staged full-GPU, and GPU-aware
full-GPU profiles. Runtime job `7051975.desched1` passed both a one-rank
model-5 CPU/GPU comparison and a two-rank spherical interface comparison.
This limited smoke validates compiler compatibility; the 58-task release
matrix above remains the broader correctness evidence. Cray MPICH 9.0.0 was
not selected because Derecho labels that module functional-only pre-release
and unsuitable for performance measurements.

Delta RH96 independently passed the same 58 tasks and 51 evidence pairs from
the same exact source archive. Build jobs `20926782`, `20926783`, and
`20926784` compiled the same ten profiles:

| Ranks / GPUs | Job | Runtime tasks | Evidence pairs | Result |
| ---: | --- | ---: | ---: | --- |
| 1 / 1 A100 | `20926787` | 45/45 | 38/38 | Passed |
| 2 / 2 A100 | `20926788` | 11/11 | 12/12 | Passed |
| 4 / 4 A100 | `20926789` | 2/2 | 1/1 | Passed |

Delta also passed all 18 restart comparisons (`20926790`), both DYN_TN
model gates, and the non-Level-Set isolation checkpoint (`20926791`). The
RH96 stack was NVHPC 26.5, CUDA 13.2, Cray MPICH 9.1.0, and FFTW 3.3.10.11.
It compiled the unchanged batched Level Set source in under one minute per
parallel build job. The previous production stack's NVHPC 25.3 front end had
stalled while compiling the large batched OpenACC smoothing routines; no
unbatched source fallback or optimization-level reduction was retained.

Delta applies Slurm GPU cgroups for the one-rank-per-GPU launch. Cray MPI 9.1
documents that these cgroups can prevent peer IPC handles from opening, so the
multi-GPU validation retained `--gpus-per-task=1 --gpu-bind=single:1` and set
`MPICH_GPU_IPC_ENABLED=0`. The `USE_GPU_AWARE_MPI=ON` LESGO code path still
passes device buffers to MPI, while Cray MPI uses its documented two-copy
intra-node fallback instead of GPU peer IPC. This is a transport constraint of
the Delta launch environment, not a solver-source change.

These small matrices are correctness tests, not production performance or
scaling benchmarks. They intentionally include checkpoint snapshots and
diagnostics that would distort a timestep-performance claim.

## Numerical Behavior Changes

CPU-reference behavior changes are reported separately from GPU-port
roundoff:

- the bottom `RHSz` history correction and `rtnewt` admissibility correction
  were already present in the required starting commit `2705e7c`; this review
  branch does not rewrite that supplied base commit, so they are separated in
  numerical reporting rather than retroactively split in Git history;
- model-4 beta modification now honors the runtime switch instead of a local
  constant and uses global physical z distance on every MPI rank;
- stress extrapolation now uses a source-stable snapshot;
- the CPU Level Set SOR smoother uses deterministic multicolor ordering: two
  colors for even periodic extents and three for odd extents.

The old/new one-step CPU comparison measured a relative change of about
`2.7e-10` in global mean kinetic energy and about `9e-9` in integrated
streamwise IBM force. These are CPU algorithm corrections, not evidence of a
CPU/GPU mismatch.

## Acceptance Rule

Any future Level Set change must retain:

1. CPU and GPU build success with MPI and without MPI;
2. a non-Level-Set build to prove option isolation;
3. SGS-disabled and SGS-model 1-5 comparisons;
4. multi-rank model-4/5 correctness;
5. optional-path checks for every changed runtime branch;
6. finite divergence, kinetic energy, and checkpoint fields;
7. no new full-field host/device transfer in the timestep path.
