# ATM Line Model

This case enables the actuator line / 5 MW turbine path with `USE_ATM=ON` and
disables the actuator disk model with `USE_TURBINES=OFF`.

The required NREL 5 MW and airfoil inputs are included in `inputATM/`.

## Structural Solver And Restart

The executable contains both rigid and structural ATM paths. Structure is a
runtime selection, not a separate GPU build:

```bash
# Rigid turbine (default)
qsub -v RUN_PROFILE=gpu,RUN_LABEL=structure_off submit_derecho.pbs

# Flexible turbine structure
qsub -v RUN_PROFILE=gpu,RUN_LABEL=structure_on,LESGO_ATM_STRUCTURE=1 \
  submit_derecho.pbs
```

When structure is enabled, LESGO advances blade/tower structural state and
feeds displacement, velocity, and aerodynamic history back into the turbine
model. Restart validation should compare an uninterrupted run with a split run
and check both velocity checkpoints and ATM restart files. See
`docs/restart_state_contract.md` for the complete state contract.

## Files

- `lesgo.conf`: ATM input settings.
- `inputATM/`: 5 MW turbine, airfoil, and turbine-array inputs.
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
qsub -N lesgo_atm_gpu -v RUN_PROFILE=gpu,RUN_LABEL=gpu submit_derecho.pbs
qsub -N lesgo_atm_cpu -v RUN_PROFILE=cpu,RUN_LABEL=cpu \
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
USE_ATM=ON
USE_CPS=OFF
USE_SCALARS=OFF
USE_LVLSET=OFF
```

The CPU profile uses `USE_CPU_BUILD=ON` and `USE_LES_GPU=OFF`. Edit
`lesgo.conf` and `inputATM/*` directly before compiling/running.

For `RUN_PROFILE=cpu`, `submit_derecho.pbs` loads the CPU runtime modules and
does not set the Cray GPU-aware MPI environment variables. For
`RUN_PROFILE=gpu`, it loads CUDA and enables GPU-aware MPI.
