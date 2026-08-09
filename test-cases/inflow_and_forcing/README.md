# Inflow And Atmospheric Forcing

This compact `64^3` case groups related inflow and atmospheric-forcing paths
behind one HIT-capable executable.

| Variant | Input behavior |
| --- | --- |
| `hit` | Blends the included `32^3` homogeneous-isotropic-turbulence field into the inlet fringe. |
| `coriolis_sponge` | Applies fixed geostrophic Coriolis forcing and damps fluctuations near the domain top. |
| `shifted` | Samples an upstream plane, shifts it laterally, and recycles it through the fringe. |

The HIT and shifted variants use stress-free initialization because their
specialized inflow types do not call the uniform-field initializer. The
submission script fixes `LESGO_RANDOM_SEED` so paired CPU/GPU runs begin from
the same perturbation field.

## Files

- `lesgo.conf`: HIT inflow configuration.
- `lesgo_coriolis_sponge.conf`: Coriolis and sponge configuration.
- `lesgo_shifted.conf`: shifted-inflow configuration.
- `HITData/`: compressed, deterministic `32^3` HIT input fields; the submit
  script expands them only in the run directory.
- `compile_derecho.sh`: common HIT-capable CPU/GPU build.
- `submit_derecho.pbs`: runs one selected configuration.

## Derecho

```bash
./compile_derecho.sh gpu
qsub submit_derecho.pbs
qsub -v CASE_VARIANT=coriolis_sponge submit_derecho.pbs
qsub -v CASE_VARIANT=shifted submit_derecho.pbs
```

The same executable supports every variant because `USE_HIT=ON` adds HIT
support without forcing `inflow_type=2`. For a CPU baseline:

```bash
./compile_derecho.sh cpu
qsub -q develop -N inflow_cpu -v RUN_PROFILE=cpu,CASE_VARIANT=hit \
  -l select=1:ncpus=4:mpiprocs=1:mem=16gb \
  -l place=pack \
  submit_derecho.pbs
```

Derecho's `develop` queue provides shared CPU nodes for compact tests. The
default GPU submissions remain in `main`.

A successful run writes `vel.out.c0`, reports `MPI_EXIT_STATUS=0`, and has
finite kinetic-energy and divergence diagnostics. HIT also writes
`restartHIT.dat`. The included HIT field is a compact functional input, not a
production atmospheric inflow database.

`dyn_init=100` keeps dynamic-SGS startup outside these 40-step examples. This
separates inflow/forcing behavior from the SGS cadence covered by
`channel_flow`.

On another cluster, replace the module stack, compiler wrapper, PBS resources,
and launch command. The physical configuration files do not depend on Derecho.
