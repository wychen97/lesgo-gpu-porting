# GPU-Port Refactor History

This document is a historical record of the GPU-port architecture plan that led
to the current optimized non-LVLSET branch.  It is not the active implementation
checklist.  For current build profiles, production scope, and validation rules,
use `docs/build_profiles.md`, `docs/code_organization.md`, and
`docs/gpu_module_contracts.md`.

The original plan started from a validated GPU solver baseline and used the
external `GPU-port` branch as an architectural reference, not as a blind code
replacement.

## Baseline

- Baseline commit: `b15f539 baseline: current validated GPU solver`
- Working branch: `codex/gpu-port-architecture-refactor`
- External reference copy on Derecho:
  `/glade/work/wchen/lesgo-gpu-port-reference`

## Core Concepts Adopted Or Evaluated

1. Replace managed-memory residency for hot arrays with explicit device
   ownership.
   - Current code uses many `managed` allocatables in `sim_param`,
     `sgs_param`, pressure, derivative, and SGS/filter paths.
   - External code uses ordinary allocatables plus `!$acc declare create`
     and explicit `update device/self` at known CPU/GPU boundaries.

2. Convert active Lagrangian SGS from plane-by-plane execution to full local-z
   execution.
   - Current bottleneck: active SGS/filter workflow creates synchronization and
     residency debt before `calc_Sij`.
   - External code allocates 3D SGS scratch arrays and processes the full local
     slab in `lagrange_Sdep_gpu`.

3. Batch spectral filters over the full local z-slab.
   - Current production path batches fields but still behaves structurally like
     repeated horizontal-plane filtering.
   - External code uses `test_filter_b_gpu(f, nz)` and
     `test_test_filter_b_gpu(f, nz)` to run slab-batched filter work.

4. Use one async GPU queue for kernels and cuFFT.
   - External `fft_gpu` binds cuFFT plans to OpenACC async queue 1.
   - Explicit waits are placed only before host-visible operations or MPI.

5. Keep CPU/GPU boundaries explicit and sparse.
   - Wall stress, output, checkpointing, and rare diagnostics may read host
     fields.
   - The timestep hot path should not repeatedly expose full fields to the CPU.

## Historical Integration Order

### Phase 1: Principle Extraction, No Direct Module Copy

- Keep the external GPU-port branch outside this repository as a reference only.
- Do not add copied external modules to the LESGO source tree.
- Do not route production calls through external-module wrappers.
- Extract the design pattern instead: explicit device ownership, persistent
  scratch buffers, stream-ordered cuFFT, sparse synchronization, and staged
  CPU/GPU boundary updates.
- Goal: rewrite each existing LESGO module in-place using these patterns while
  preserving local interfaces and numerical behavior.

### Phase 2: SGS/Lagrangian Path First

- Introduce explicit device residency for SGS arrays:
  `S11..S33`, `Nu_t`, `Cs_opt2`, `F_LM/F_MM/F_QN/F_NN`, `Beta`, `Tn_all`.
- Port slab-batched active Lagrangian SGS update.
- Validate against current managed path on:
  - 480x240x240, 2 MPI / 2 GPU
  - 3072x384x50 slice, 1 MPI / 1 GPU

### Phase 3: Extend Device Residency Across Main Hot Path

- Move `u/v/w`, derivatives, RHS, pressure gradients, divstress arrays to
  explicit device ownership.
- Reduce full-array `update self/device` to known CPU boundary slices only.
- Validate each module independently.

### Phase 4: Pressure/Convection/Derivative Alignment

- Compare the algorithms and memory-layout choices in the reference branch
  against current CUDA Fortran kernels.
- Reimplement useful ideas inside the existing LESGO modules rather than
  calling copied reference routines.
- Keep the numerically stronger path in each module.

## Non-Copy Rule

The external branch is a design reference, not a source of drop-in modules.
Direct copies are allowed only for temporary private analysis outside the LESGO
source tree. Any production change must be an in-place LESGO implementation
with a before/after benchmark against the immediate previous accepted version.

## Rejection Rule

Any imported concept must be reverted unless it satisfies:

- Build succeeds.
- Numerical diagnostics match the immediate previous accepted version.
- Timing improves or the code enables a clearly necessary future improvement.
- Fallback path remains available until the new path is validated.
