# Scalar Transport

This compact `64^3` case demonstrates the LESGO scalar module without turbine
physics. Two configurations use the same executable:

| Variant | Input | Meaning |
| --- | --- | --- |
| `passive` | `lesgo.conf` | Transport scalar concentration without momentum feedback. |
| `active` | `lesgo_active.conf` | Treat the scalar as temperature and apply buoyancy feedback. |

The GPU build requires both `USE_SCALARS=ON` and `USE_SCALARS_GPU=ON`.
`USE_SCALARS_GPU` alone is invalid because it is only the GPU implementation of
the parent scalar module.

## Files

- `lesgo.conf`: passive-scalar configuration.
- `lesgo_active.conf`: active temperature/buoyancy configuration.
- `compile_derecho.sh`: CPU/GPU scalar build.
- `submit_derecho.pbs`: run one selected variant.

## Derecho

```bash
./compile_derecho.sh gpu
qsub submit_derecho.pbs
qsub -v CASE_VARIANT=active submit_derecho.pbs
```

For a CPU baseline:

```bash
./compile_derecho.sh cpu
qsub -q develop -N scalar_cpu -v RUN_PROFILE=cpu,CASE_VARIANT=passive \
  -l select=1:ncpus=4:mpiprocs=1:mem=16gb \
  -l place=pack \
  submit_derecho.pbs
```

Derecho's `develop` queue provides shared CPU nodes for compact tests. The
default GPU submissions remain in `main`.

Runs are isolated under `runs/<profile>-<variant>/`. A successful run writes
`scal.out.c0`, finishes with `MPI_EXIT_STATUS=0`, and reports finite kinetic
energy and divergence. The short run verifies configuration and execution; it
is not a statistically converged atmospheric scalar simulation. The submit
script also fixes `LESGO_RANDOM_SEED` for reproducible paired runs.

`dyn_init=100` intentionally keeps the dynamic-SGS update outside this 40-step
example. SGS variants and update cadence are covered by `channel_flow`; this
case isolates scalar transport and buoyancy.

On another cluster, replace the module stack, `ftn`, PBS resources, and
`mpiexec` command while preserving one MPI rank per GPU.
