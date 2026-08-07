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

The existing public calls remain the interface:

- `level_set_BC()` selects CPU or GPU stress treatment;
- `level_set_lag_dyn()` selects CPU or GPU model-4 preconditioning;
- `level_set_forcing()` selects CPU or GPU immersed-boundary forcing;
- `level_set_smooth_vel()` selects the matching smoother.

GPU callers therefore do not depend on implementation-specific public entry
points.

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

Validation used the `test-cases/level_set_cubes` geometry on a
`64 x 64 x 32` grid with 10 timesteps. CPU and GPU builds used identical
inputs. Checkpoints include velocity, RHS, SGS coefficient, and Lagrangian
history fields.

| Coverage | Result |
| --- | --- |
| One rank, SGS disabled | Passed |
| One rank, SGS models 1-5 | Passed |
| Two ranks, SGS models 4 and 5, compact host-staged MPI | Passed |
| Non-MPI model 5 | Passed |
| LES-GPU plus CPU-Level-Set host bridge, model 5 | Passed |
| Optional logarithmic/simple/legacy stress paths | Passed |
| `use_modify_dutdn`, desired log profile, direct local treatment | Passed |
| `smooth_mode='3d'` exact CPU-order wavefront | Passed |

The worst normalized CPU/GPU field difference in the 16-run primary matrix was
`3.302591e-13` (two-rank model 4, `RHSz`). The `smooth_mode='3d'` comparison
was `2.614140e-15`. The host-bridge regression differed from the CPU reference
by at most `1.794502e-13`. Divergence and kinetic energy matched at printed
precision.

Nine clean build profiles passed: Level Set CPU/GPU with MPI and without MPI,
non-Level-Set CPU/GPU, forced GPU-aware MPI Level Set/non-Level-Set builds, and
the host-bridge fallback. A 10-step non-Level-Set channel checkpoint from the
new source was byte-for-byte identical to the untouched `origin/main` GPU
checkpoint at commit `3e83529438310bfac969f72a191b6dd864cea6fc`; the Level Set
work therefore does not alter that build path.

The `use_enforce_un` branch matched for a one-step smoke test. Longer runs of
that experimental CPU boundary condition require a separate physical-stability
study; CPU/GPU parity alone is not evidence that the model is suitable for a
production simulation.

## Performance Sanity Check

The small one-rank case includes diagnostics every timestep, so it is a
correctness-oriented performance sanity check rather than a production scaling
benchmark.

| Runtime path | CPU wall | GPU wall | CPU/GPU |
| --- | ---: | ---: | ---: |
| SGS disabled, 10 steps | 0.28837 s | 0.04796 s | 6.01x |
| SGS model 1, 10 steps | 0.29087 s | 0.04762 s | 6.11x |
| SGS model 5, 10 steps | 0.38014 s | 0.05506 s | 6.90x |

The two-rank fallback test shared one physical GPU and is correctness evidence,
not a scaling result. A forced CUDA-aware MPI build succeeds on OpenMPI+UCX,
but device-direct multi-rank runtime acceptance requires one MPI rank per GPU
on a multi-GPU node. Do not use a multi-rank/one-GPU run as GPU-aware MPI
validation.

## Acceptance Rule

Any future Level Set change must retain:

1. CPU and GPU build success with MPI and without MPI;
2. a non-Level-Set build to prove option isolation;
3. SGS-disabled and SGS-model 1-5 comparisons;
4. multi-rank model-4/5 correctness;
5. optional-path checks for every changed runtime branch;
6. finite divergence, kinetic energy, and checkpoint fields;
7. no new full-field host/device transfer in the timestep path.
