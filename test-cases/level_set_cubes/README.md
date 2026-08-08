# Level Set Cubes

This case exercises immersed-surface Level Set treatment on a
`64 x 64 x 32` uniform grid. Four simple vertical solids are generated from
`trees.conf`; no turbine, ATM, scalar, CPS, or HIT module is enabled.

## Files

- `lesgo.conf`: 20-step Level Set input with SGS model 5.
- `trees.conf`: four deterministic geometry primitives used to generate
  `phi.out` at startup.
- `compile_derecho.sh`: compile-only CPU/GPU example.
- `submit_derecho.pbs`: run-only Derecho PBS example.

## Build And Run

```bash
./compile_derecho.sh gpu
qsub submit_derecho.pbs
```

For a CPU/GPU comparison:

```bash
./compile_derecho.sh cpu
./compile_derecho.sh gpu
qsub -N level_set_gpu -v RUN_PROFILE=gpu,RUN_LABEL=gpu submit_derecho.pbs
qsub -N level_set_cpu -v RUN_PROFILE=cpu,RUN_LABEL=cpu \
  -l select=1:ncpus=32:mpiprocs=1:mem=120gb submit_derecho.pbs
```

The GPU build enables:

```text
USE_LES_GPU=ON
USE_LVLSET=ON
USE_LVLSET_GPU=ON
```

The CPU reference keeps `USE_LVLSET=ON` but sets `USE_LES_GPU=OFF` and
`USE_LVLSET_GPU=OFF`. Both profiles use exactly the same `lesgo.conf` and
`trees.conf`.

Generated `phi.out.c*`, `norm.dat.c*`, checkpoints, output, binaries, and logs
are runtime artifacts and are not source files.

For the full CPU/host-bridge/GPU matrix, prepare isolated cases with:

```bash
python3 tools/prepare_level_set_validation.py \
  --out /path/to/lvlset-matrix
```

The generated manifest fixes `DYN_init=1` and `cs_count=2`, covers SGS off and
models 1-5, optional Level Set modes, model-4 beta/lower-wall variants,
MPI/non-MPI builds, and spherical geometries crossing a rank boundary. Strict
checkpoint variants use two timesteps so an active dynamic-SGS update is
compared before long-horizon chaotic amplification can obscure implementation
parity; the checked-in standalone case remains a 20-step smoke test.
`tools/run_level_set_validation.py` runs tasks inside an existing scheduler
allocation. `LESGO_LVLSET_VALIDATION_SNAPSHOT=ON` writes the additional fields
consumed by `tools/compare_level_set_checkpoints.py`.

## MPI And GPU Counts

The default input has `nproc=1`. For multiple GPUs, change all of the following
to the same decomposition:

1. `nproc` in `lesgo.conf`;
2. `MPI_RANKS` and `MPI_PPN` in `submit_derecho.pbs`;
3. `select`, `mpiprocs`, and `ngpus` in the PBS resource line.

Use one MPI rank per GPU. `Nz` must be divisible by `nproc`. The default
`smooth_mode='xy'` supports MPI; `smooth_mode='3d'` is single-rank only.
