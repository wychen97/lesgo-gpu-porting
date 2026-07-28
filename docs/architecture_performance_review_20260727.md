# Architecture And Performance Review: 2026-07-27

## Scope

This review compares the current working candidate with public `main` commit
`372964cdef67ec5019b1140e7f7646f9eb302ce5`. It focuses on the validated
non-LVLSET production path. No commit from this review has been pushed.

The complete source and CMake surface was checked, but source changes were kept
to the ATM boundary where the audit found redundant synchronization, dead
policy layers, optional-path correctness gaps, and repeated hot-path queries.

## Architecture Findings

### Production path

The normal `sampling = "atPoint"` GPU path now has one clear ownership model:

- `u`, `v`, `w`, `fxa`, `fya`, and `fza` remain device-resident;
- the batched sampler reads LES fields directly on the device;
- the host turbine model remains authoritative for controller, aerodynamic,
  output, and structural state;
- flattened blade state is transferred only at the explicit host/device
  boundary;
- the batched convolution deposits directly into resident LES force arrays.

The old interface contained disabled CUDA-era policy functions, shadow arrays,
blade mirrors, experimental point-owner load-balancing branches, an unreachable
legacy gather implementation, and duplicate per-turbine GPU routines. Those
paths were not selectable by the current configuration and have been removed.
The file decreased from about 142 kB to about 101 kB.

### Optional sampling

Spalart sampling and nacelle interpolation are established host algorithms.
They now use an explicit compatibility bridge only when requested. The bridge:

- waits for the required LES queue at a named boundary;
- refreshes host velocity fields;
- runs the host sampling/convolution algorithm;
- transfers only the resulting optional force entries back to resident LES
  arrays.

The standard atPoint case no longer allocates or refreshes the full host
`w_uv` volume.

### Input and restart invariants

ATM input now rejects invalid intervals, dimensions, sampling names, blade
epsilon, and nacelle parameters before entering the timestep loop. In
particular, enabling nacelle drag without a positive `nacelleEpsilon` no longer
reaches a divide-by-zero/NaN path.

`outputInterval = 0` remains a valid way to disable periodic ATM output. Its
cadence test now uses explicit control flow because Fortran does not guarantee
short-circuit evaluation of `.and.` expressions.

`updateInterval > 1` holds a previously computed grid-force field between
aerodynamic updates. The historical restart format does not serialize that
held field. Therefore:

- checkpoints at a multiple of `updateInterval` are supported because the
  first resumed step performs a normal force update;
- a non-aligned restart is rejected explicitly instead of silently injecting
  a force and power transient;
- arbitrary-seam support would require a versioned, per-rank ATM force-cache
  sidecar and is not hidden inside this cleanup.

Standard `updateInterval = 1` restarts remain unrestricted.

## Readability Changes

- Removed unreachable implementations rather than preserving false runtime
  switches.
- Reduced the ATM interface navigation map to active ownership boundaries.
- Replaced repeated structure environment lookups in blade, yaw, output, and
  rotation loops with one cached decision per routine.
- Named and validated the nacelle velocity-correction denominator.
- Made nacelle defaults deterministic when the feature is disabled.
- Kept CPU fallback and optional compatibility paths explicit.

The CPU-oriented turbine physics remains recognizable. The cleanup changes the
GPU coupling layer, not the mathematical organization of the CPU model.

## Performance Result

The paired benchmark used one Derecho A100, one MPI rank, a `240 x 240 x 240`
ATM case, SGS model 5, and 200 steps. Two repetitions were run in interleaved
baseline/candidate order with identical input and seed.

| Version | Mean cumulative solver time | Mean solver time per step |
| --- | ---: | ---: |
| Accepted `main` baseline | 12.38198 s | 0.0619099 s/step |
| Reviewed candidate | 12.25045 s | 0.0612523 s/step |

The candidate reduces cumulative solver time by about `1.06%`. Mean full
simulation wall time decreased from `13.38464 s` to `13.25396 s`, about
`0.98%`.

The ATM subtotal printed by the internal component timer decreased much more
than the end-to-end result because removing a global wait changes which timer
is charged for already queued GPU work. Cumulative solver time is the
authoritative performance measure.

## Correctness And Build Gates

| Gate | Result |
| --- | --- |
| Standard structure-OFF A/B, `240^3`, 200 steps | Velocity checkpoint and ATM restart byte-identical |
| Structure-ON A/B, `96^3`, 20 steps | Velocity checkpoint, ATM restart, and power byte-identical |
| CPU/GPU atPoint+nacelle, `96^3`, 12 steps | Finite; velocity relative L2 `1.13e-15` |
| CPU/GPU Spalart, `96^3`, 12 steps | Finite; velocity relative L2 `1.14e-15` |
| CPU/GPU `updateInterval=4`, `96^3`, 12 steps | Finite; velocity relative L2 `1.18e-15` |
| Maximum optional-case power absolute difference | `1.86e-9` |
| Standard continuous/split restart | Byte-identical |
| Aligned `updateInterval=4` continuous/split restart | Byte-identical |
| Unaligned `updateInterval=4` restart | Rejected with explicit configuration error |
| `outputInterval=0` CPU/GPU smoke | Passed without modulo-by-zero |
| Derecho ATM GPU and CPU builds | Passed |
| Delta ATM GPU and CPU builds | Passed |
| ATM+scalars GPU/CPU and serial ATM GPU builds | Passed |
| Source, interface, macro, comment, and packed-path checks | Passed |

## Remaining Boundaries

- Spalart and nacelle compatibility are correct but still host-assisted. Port
  them only after a representative case shows that this bridge is material.
- Arbitrary-seam `updateInterval > 1` restart requires a new force-cache
  checkpoint contract.
- LVLSET remains outside the validated optimized branch.
- Large wind-farm scaling should be measured separately; a periodic output or
  SGS-update timestep must not be used as the sole performance sample.
