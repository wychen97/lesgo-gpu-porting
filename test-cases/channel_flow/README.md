# Channel Flow

This case isolates the optimized LES core on a 240 x 240 x 240 grid with no
turbine, precursor, scalar, HIT, CGNS, or level-set modules enabled.

## Files

- `lesgo.conf`: channel-flow input settings.
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
USE_ATM=OFF
USE_CPS=OFF
USE_SCALARS=OFF
USE_LVLSET=OFF
```

It builds one local executable, `lesgo-run-exe`. `submit_derecho.pbs` cleans
old output files and runs that executable as a Derecho PBS job. Edit
`lesgo.conf` directly before compiling/running. The `TURBINES` block in
`lesgo.conf` is ignored because both turbine models are disabled at configure
time.
