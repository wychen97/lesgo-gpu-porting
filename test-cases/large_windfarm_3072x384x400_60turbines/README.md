# LESGO 3072 x 384 x 400 / 60-Turbine Input Case

Purpose: provide the initial configuration files for the large wind-farm case
used during GPU scaling diagnostics. This directory is retained in the final
optimized non-LVLSET branch so collaborators can reproduce the 3072 x 384 x 400
/ 60-turbine benchmark setup from a clean checkout.

Default benchmark geometry:

- Grid: 3072 x 384 x 400
- Domain: 28224 m x 3780 m x 2000 m
- Turbines: 10 streamwise rows x 6 spanwise columns = 60 turbines
- Turbine scale: NREL 5MW-like, D = 126 m, hub height from turbine type file
- Spacing: 7D streamwise, 5D spanwise
- First row: 14D from inlet/periodic origin
- Turbine array file: `inputATM/turbineArrayProperties`

The 10 x 6 orientation is intentional even though the shorthand request was
"6 x 10": the target domain is long in x and narrow in y, so 10 streamwise rows
and 6 spanwise columns is the physically consistent interpretation.

The LESGO dimensional scale is kept at `z_i = 1000 m`; therefore the
non-dimensional domain in `lesgo.conf` is `Lx = 28.224`, `Ly = 3.78`,
`Lz = 2.0`.

## Files

- `lesgo.conf`: baseline large-case configuration.
- `inputATM/`: 60-turbine ATM input files and airfoil tables.
- `compile_derecho.sh`: compile-only Derecho/NVHPC example.
- `submit_derecho.pbs`: run-only Derecho PBS submission example.

## Recommended sequence on Derecho

```bash
cd /path/to/lesgo/test-cases/large_windfarm_3072x384x400_60turbines
./compile_derecho.sh
qsub submit_derecho.pbs
```

The default large-case submission settings are:

```text
MPI_RANKS=16
MPI_PPN=4
PBS select=4 nodes x 4 A100 GPUs/node
grid and timestep settings are read directly from lesgo.conf
```

Edit `lesgo.conf` directly when running a smaller smoke test or changing
`nsteps`, `wbase`, or output controls. If you change the GPU count, update both
the PBS resource line and the `MPI_RANKS` / `MPI_PPN` values in
`submit_derecho.pbs`.

For a physically developed turbulent wind-farm field, use longer spin-up
settings or restart from an already developed precursor. The short benchmark
job is intentionally focused on kernel/module timing and does not claim
converged wind-farm statistics.

## Public GPU-Porting Build Notes

This case is meant to exercise the optimized non-LVLSET GPU path. The GPU
benchmark executable is built through:

```bash
-DUSE_LES_GPU=ON
```

Use this README and the top-level `README.md` for the current validated CMake
profile and module scope. Keep generated run directories, binaries, and queue
logs out of git.
