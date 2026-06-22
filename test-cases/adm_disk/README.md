# ADM Disk Model

This case enables the actuator disk model through `USE_TURBINES=ON` and keeps
the line-model ATM path disabled with `USE_ATM=OFF`.

The ADM parameters are read from the `TURBINES` block in `lesgo.conf`; this
case does not use `inputATM/`.

## Files

- `lesgo.conf`: ADM input settings.
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
USE_TURBINES=ON
USE_ATM=OFF
USE_CPS=OFF
USE_SCALARS=OFF
USE_LVLSET=OFF
```

It builds one local executable, `lesgo-run-exe`. `submit_derecho.pbs` cleans
old output files and runs that executable as a Derecho PBS job. Edit
`lesgo.conf` directly before compiling/running.
