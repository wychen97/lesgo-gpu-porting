# Concurrent Precursor

This case demonstrates LESGO concurrent precursor simulation (CPS). An
upstream `RED` domain samples its flow and sends it to the fringe region of a
downstream `BLUE` domain during every timestep.

| Variant | Build command | Exchanged fields |
| --- | --- | --- |
| `velocity` | `./compile_derecho.sh gpu velocity` | `u`, `v`, and `w` |
| `scalar` | `./compile_derecho.sh gpu scalar` | `u`, `v`, `w`, and scalar |

Both domains are `64^3` with one local MPI rank. The job therefore launches an
MPI world of two ranks but each input correctly contains `nproc=1`. The RED
producer samples its flow while BLUE receives those planes through its fringe.
Both use stress-free initialization and no pressure-gradient forcing so this
short example isolates CPS exchange. Separate working directories prevent
output-file collisions. The submission script fixes `LESGO_RANDOM_SEED` so
CPU/GPU comparisons start from the same perturbation field.

## Files

- `lesgo.conf`, `lesgo_blue.conf`: RED and BLUE velocity-only CPS inputs.
- `lesgo_scalar.conf`, `lesgo_scalar_blue.conf`: RED and BLUE inputs with
  scalar coupling.
- `compile_derecho.sh`: builds a selected CPU/GPU CPS variant.
- `submit_derecho.pbs`: launches the red/blue MPMD pair.

## Derecho

```bash
./compile_derecho.sh gpu velocity
qsub submit_derecho.pbs

./compile_derecho.sh gpu scalar
qsub -v CASE_VARIANT=scalar submit_derecho.pbs
```

The GPU example requests two A100s because CPS runs two simulations
concurrently. A CPU run can use two CPU MPI ranks and no GPU request:

```bash
./compile_derecho.sh cpu velocity
qsub -q develop -N cps_cpu -v RUN_PROFILE=cpu,CASE_VARIANT=velocity \
  -l select=1:ncpus=4:mpiprocs=2:mem=32gb \
  -l place=pack \
  submit_derecho.pbs
```

Derecho's `develop` queue provides shared CPU nodes for compact tests. The
default two-A100 submission remains in `main`.

A successful run produces checkpoints in both `red/` and `blue/`, reports
`MPI_EXIT_STATUS=0`, and prints CPS stage timing. These short runs verify data
exchange and execution, not precursor statistical convergence.

`dyn_init=100` keeps the dynamic-SGS update outside this 30-step exchange test
so the example isolates CPS communication. Use `channel_flow` for SGS cadence
validation.

On another cluster, replace the module stack, compiler wrapper, scheduler
resources, MPMD `mpiexec` syntax, and GPU-rank binding command.
