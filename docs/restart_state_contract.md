# Restart State Contract

LESGO restart output is a synchronized state snapshot, not only a velocity
field. A resumed run is valid only when every enabled module reads state from
the same `jt_total` and `total_time`.

## Core LES State

| File | State | Continuation role |
| --- | --- | --- |
| `vel.out.c*` | `u`, `v`, `w` | Current velocity field. |
| `vel.out.c*` | `RHSx`, `RHSy`, `RHSz` | Previous Adams-Bashforth RHS history. |
| `vel.out.c*` | `Cs_opt2`, `F_LM`, `F_MM`, `F_QN`, `F_NN` | Dynamic and Lagrangian SGS history. |
| `total_time.dat` | step, physical/nondimensional time, `dt`, CFL | Cumulative timestep state. |
| `total_time.dat` | `lagran_dt` | Time accumulated since the previous Lagrangian SGS update. |

`total_time.dat` now has one integer followed by five real values. Historical
five-value files remain readable; their missing `lagran_dt` is initialized to
zero with a warning.

Dynamic SGS decisions use `jt_total`, not the process-local `jt`. This applies
to warm-up, first initialization at `DYN_init`, update cadence, and the
Lagrangian time accumulator. A process restart therefore does not restart the
physical SGS schedule.

In an explicit-residency GPU build, `initial()` reads velocity, AB2, and SGS
state on the host. `initialize()` transfers all of that persistent state to the
device before the first resumed timestep. Omitting the RHS transfer changes the
flow during the first resumed integration step even when turbine output at that
step still appears continuous.

## Level Set Lagrangian SGS State

Level Set runs with SGS model 4 or 5 write
`lvlset_sgs_restart.out.c*`. Version 1 contains:

- a format marker and local grid/decomposition dimensions;
- per-rank `phi` sum, absolute sum, minimum, and maximum signatures;
- `Beta` and `Tn_all` on physical planes `1:nz`.

The main LES checkpoint already owns `Cs_opt2`, `F_LM`, `F_MM`, `F_QN`, and
`F_NN`. The Level Set sidecar completes the model-4/5 state needed at a restart
seam between dynamic updates. A requested model-4/5 restart without this
sidecar is rejected; silently replacing the missing arrays with initialization
values can create a force/coefficient transient. Grid, decomposition, or
geometry-signature mismatches are also rejected.

The current contract is for the uniform-grid backend. An AMR restart format
must additionally serialize the level hierarchy, per-level refinement ratios,
patch valid/ghost boxes and geometry generations, coverage masks,
level-specific `phi`/normal ownership, coarse/fine correction state, and any
held defect-projection state. Restart must restore those records before device
geometry is rebound and must reject a hierarchy or generation mismatch. A
uniform sidecar must not be reused for an AMR hierarchy.

## ATM State

Each `turbineOutput/TURBINE_*/restart` file retains the historical controller
prefix and appends `ATM_RESTART_STATE_V2`. The V2 block records:

- restart version, step, time, turbine dimensions, and enabled-feature flags;
- rotor azimuth and nacelle-yaw increment;
- blade points;
- induced-velocity histories `du` and `uy_opt_vec`;
- flexible-structure geometry, displacement, velocity, and acceleration;
- previous aerodynamic blade forces, blade-aligned bases, and pitching moment.

The last three structural load fields are consumed by the next structure solve
before new aerodynamic loads are available. A fresh start or a legacy restart
without those fields skips that first structural advance instead of applying a
zero-load impulse.

Controller rate limiting is enabled only after previous torque/controller state
is valid. The legacy controller prefix is still readable by older code and is
written with explicit full double precision.

ATM files are read only when `jt_total > 0`. Stale turbine files in a fresh case
directory are ignored. A requested LES restart with a missing ATM restart file
is an error. V2 metadata must match the LES step/time, turbine dimensions,
tip-correction setting, and structure setting.

Intermediate LES checkpoints call the ATM checkpoint writer directly. Finalize
uses the same path and suppresses duplicate writes at the same step.

When `updateInterval > 1`, the ATM grid-force field is held between aerodynamic
updates. The historical restart format does not serialize this per-rank held
field. A restart is therefore accepted only when the checkpoint step is a
multiple of `updateInterval`, so the first resumed timestep performs a normal
force update. Non-aligned restarts fail during input initialization rather than
silently creating a force/power transient. Supporting arbitrary seams requires
a versioned ATM force-cache sidecar.

## Optional Sidecars

The following existing modules retain separate checkpoint ownership:

| Module | Restart state |
| --- | --- |
| ADM disk turbines | Disk-averaged turbine velocity and averaging time. |
| Scalars | `scal.out.c*`. |
| Dynamic Taylor timescale | `dyn_tn.out.c*`. |
| Integral wall model | `iwm_checkPoint.dat`. |
| HIT inflow | `restartHIT.dat`. |
| Time averaging | `tavg.out.c*`. |

Changes to any stateful optional module must add a continuous-versus-split test
for its own sidecar. File existence alone is not sufficient evidence that a
restart is physically continuous.

For Level Set, run continuous and split cases with models 4 and 5, using
`DYN_init=1` and `cs_count=2`. Include one seam immediately before and one seam
immediately after a dynamic update. Compare `vel.out`, the validation snapshot,
`lvlset_sgs_restart.out`, divergence, kinetic energy, and integrated immersed
force for CPU, host-bridge, and full-GPU paths.

## Validation Protocol

Use `tools/compare_atm_restart_runs.py` on two runs with identical input and random
seed. Select `--case-kind channel`, `adm`, or `atm` so the comparator requires
the sidecars owned by the enabled turbine model:

1. uninterrupted run from step 0 to `N`;
2. split run from step 0 to `S`, then restart from `S` to `N`;
3. choose at least one `S` before `DYN_init` and one `S` between `cs_count`
   update boundaries;
4. test CPU and GPU builds with structure OFF and ON;
5. compare power, thrust, rotor speed, final LES checkpoint, time metadata, and
   turbine/structure state.

Example:

```bash
python3 tools/compare_atm_restart_runs.py \
  --case-kind atm \
  --continuous run/continuous \
  --restarted run/restarted \
  --seam-step 97 \
  --out run/comparison.json
```

Velocity and turbine signals use `1e-6` relative and `1e-8` absolute tolerances.
AB2 and the four Lagrangian SGS accumulators use `1e-6` relative and absolute
tolerances. The derived `Cs_opt2` ratio uses a separate `2e-4` relative and
absolute tolerance: at the stress-free top boundary, near-zero dynamic-model
denominators can amplify machine-precision MPI-overlap differences in this
coefficient while velocity, RHS, and the saved accumulators remain converged to
much tighter tolerances. The formatted turbine/structure sidecars use `1e-4`
relative and absolute tolerances. Byte identity is reported but is not the
physical acceptance criterion.

## Validation Snapshot: 2026-07-12

The compact validation used a `120 x 120 x 120` ATM domain with two turbines,
`sgs_model = 5`, `cs_count = 5`, and `DYN_init = 100`. CPU/GPU structure OFF/ON
tests passed at a 20-step seam. The 130-step runs split at step 97 and also
passed, covering a checkpoint between SGS updates and the first dynamic update
at step 100.

The representative validation used the same physics on a `240 x 240 x 240`
grid. GPU cases used one A100 rank. CPU cases used eight MPI ranks and compared
all eight `vel.out.c*` files.

| Build | Structure | Steps / seam | Pass | Max signal relative error | Max velocity relative error | Max history relative error | Max sidecar relative error |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| GPU, 1 rank | OFF | 200 / 100 | yes | 0 | 0 | 0 | 0 |
| GPU, 1 rank | ON | 200 / 100 | yes | 0 | 0 | 0 | 0 |
| CPU, 8 ranks | OFF | 200 / 100 | yes | `2.58e-12` | `3.37e-14` | `2.77e-8` | `4.21e-14` |
| CPU, 8 ranks | ON | 200 / 100 | yes | `2.00e-9` | `9.38e-11` | `8.24e-8` | `3.82e-5` |

The final source snapshot compiled successfully for CPU and ATM/LES GPU targets
with Derecho NVHPC 25.9 and Delta `PrgEnv-nvidia`/NVHPC 25.3.

### Standard test-case suite

The checked-in standard channel, ADM disk, and ATM line cases were also tested
at `240 x 240 x 240`, `sgs_model = 5`, 200 final steps, and a 100-step restart
seam. CPU runs used eight MPI ranks on shared Derecho nodes; GPU runs used one
A100 and one MPI rank. The ADM case contained 24 disks and the ATM case contained
two NREL 5 MW line turbines.

| Case | Build | Pass | Max velocity abs. error | Max RHS abs. error | Max `Cs_opt2` abs. error | Max SGS accumulator abs. error | Turbine-output abs. error |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Channel | CPU, 8 ranks | yes | `1.31e-14` | `8.99e-13` | `1.15e-4` | `1.87e-10` | n/a |
| Channel | GPU, 1 A100 | yes | 0 | 0 | 0 | 0 | n/a |
| ADM, 24 disks | CPU, 8 ranks | yes | `1.47e-14` | `8.76e-13` | `3.16e-5` | `1.12e-10` | `7.11e-14` |
| ADM, 24 disks | GPU, 1 A100 | yes | 0 | 0 | 0 | 0 | 0 |
| ATM, 2 turbines | CPU, 8 ranks | yes | `2.98e-14` | `1.77e-12` | `3.66e-8` | `2.96e-10` | `2.10e-9` |
| ATM, 2 turbines | GPU, 1 A100 | yes | 0 | 0 | 0 | 0 | 0 |

The CPU `Cs_opt2` maxima occur on the top MPI rank. They are accepted under the
coefficient-specific criterion above; the actual transported flow and saved
Lagrangian histories remain several orders of magnitude tighter. No ATM power,
thrust, or rotor-speed discontinuity was detected at the restart seam.
