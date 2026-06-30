# ATM Line Model

This case enables the actuator line / 5 MW turbine path with `USE_ATM=ON` and
disables the actuator disk model with `USE_TURBINES=OFF`.

The required NREL 5 MW and airfoil inputs are included in `inputATM/`.

## Files

- `lesgo.conf`: ATM input settings.
- `inputATM/`: 5 MW turbine, airfoil, and turbine-array inputs.
- `compile_derecho.sh`: compile-only Derecho/NVHPC example.
- `submit_derecho.pbs`: run-only Derecho PBS submission example.

## Recommended sequence

Compile and submit the default CPU run from this directory:

```bash
./compile_derecho.sh cpu
qsub submit_derecho.pbs
```

The compile script installs the CPU executable and updates the compatibility
copy:

```text
lesgo-run-exe-cpu
lesgo-run-exe
```

Run the CPU job explicitly with:

```bash
qsub -N lesgo_atm_cpu -v RUN_PROFILE=cpu,RUN_LABEL=cpu \
  -q develop -l select=1:ncpus=32:mpiprocs=1:mem=120gb \
  -l place=scatter:shared submit_derecho.pbs
```

Each run writes `lesgo_<RUN_LABEL>.log` and archives the log plus output files
under `run-archives/<RUN_LABEL>/`.  Set `RUN_LABEL` explicitly when running
multiple profiles from the same case directory so later runs do not overwrite
evidence from earlier runs.

The CPU profile configures CMake with:

```text
USE_CPU_BUILD=ON
USE_LES_GPU=OFF
USE_TURBINES=OFF
USE_ATM=ON
USE_CPS=OFF
USE_SCALARS=OFF
USE_LVLSET=OFF
```

Edit
`lesgo.conf` and `inputATM/*` directly before compiling/running.

For `RUN_PROFILE=cpu`, `submit_derecho.pbs` loads the CPU runtime modules and
does not set the Cray GPU-aware MPI environment variables.
