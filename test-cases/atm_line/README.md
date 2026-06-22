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

Compile and submit from this directory:

```bash
./compile_derecho.sh
qsub submit_derecho.pbs
```

`compile_derecho.sh` configures CMake with:

```text
USE_LES_GPU=ON
USE_TURBINES=OFF
USE_ATM=ON
USE_CPS=OFF
USE_SCALARS=OFF
USE_LVLSET=OFF
```

It builds one local executable, `lesgo-run-exe`. `submit_derecho.pbs` cleans
old output files and runs that executable as a Derecho PBS job. Edit
`lesgo.conf` and `inputATM/*` directly before compiling/running.
