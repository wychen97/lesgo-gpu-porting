# Channel Flow

This case isolates the optimized LES core on a 240 x 240 x 240 grid with no
turbine, precursor, scalar, HIT, CGNS, or level-set modules enabled.

## Files

- `lesgo.conf`: channel-flow input settings.
- `compile_derecho.sh`: compile-only Derecho/NVHPC example.
- `submit_derecho.pbs`: run-only Derecho PBS submission example.

## Recommended sequence

Compile and submit the default GPU run from this directory:

```bash
./compile_derecho.sh gpu
qsub submit_derecho.pbs
```

For paired CPU/GPU validation, build both profiles:

```bash
./compile_derecho.sh cpu
./compile_derecho.sh gpu
```

The compile script installs profile-specific executables:

```text
lesgo-run-exe-cpu
lesgo-run-exe-gpu
```

Run the paired GPU and CPU jobs with:

```bash
qsub -N lesgo_channel_gpu -v RUN_PROFILE=gpu,RUN_LABEL=gpu submit_derecho.pbs
qsub -N lesgo_channel_cpu -v RUN_PROFILE=cpu,RUN_LABEL=cpu \
  -l select=1:ncpus=32:mpiprocs=1:mem=120gb submit_derecho.pbs
```

Each run writes `lesgo_<RUN_LABEL>.log` and archives the log plus output files
under `run-archives/<RUN_LABEL>/`.  Set `RUN_LABEL` explicitly when running
multiple profiles from the same case directory so later runs do not overwrite
evidence from earlier runs.

The GPU profile configures CMake with:

```text
USE_LES_GPU=ON
USE_TURBINES=OFF
USE_ATM=OFF
USE_CPS=OFF
USE_SCALARS=OFF
USE_LVLSET=OFF
```

The CPU profile uses `USE_CPU_BUILD=ON` and `USE_LES_GPU=OFF`. Edit
`lesgo.conf` directly before compiling/running. The `TURBINES` block in
`lesgo.conf` is ignored because both turbine models are disabled at configure
time.

For `RUN_PROFILE=cpu`, `submit_derecho.pbs` loads the CPU runtime modules and
does not set the Cray GPU-aware MPI environment variables. For
`RUN_PROFILE=gpu`, it loads CUDA and enables GPU-aware MPI.
