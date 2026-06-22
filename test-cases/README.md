# Test Cases

This directory follows the upstream LESGO convention of keeping runnable cases
under `test-cases/`. Only the four cases used in the current GPU-porting
validation suite are included:

| Case | Purpose | Main CMake switches |
| --- | --- | --- |
| `channel_flow` | LES core baseline without turbine physics | `USE_LES_GPU=ON`, `USE_TURBINES=OFF`, `USE_ATM=OFF` |
| `adm_disk` | Actuator disk model increment | `USE_LES_GPU=ON`, `USE_TURBINES=ON`, `USE_ATM=OFF` |
| `atm_line` | Actuator line / 5 MW turbine path | `USE_LES_GPU=ON`, `USE_TURBINES=OFF`, `USE_ATM=ON` |
| `large_windfarm_3072x384x400_60turbines` | Large 60-turbine benchmark setup | `USE_LES_GPU=ON`, `USE_ATM=ON` |

Each case contains:

- `lesgo.conf` input settings;
- model input files when needed, such as `inputATM/`;
- `README.md` with the case-specific module setup;
- `compile_derecho.sh`, which only compiles the GPU executable;
- `submit_derecho.pbs`, which only submits/runs the case on Derecho.

The root `CMakeLists.txt` declares the available feature switches and their
neutral defaults. The case-specific feature choices live in
`compile_derecho.sh`; for example, `adm_disk` sets `USE_TURBINES=ON`, while
`channel_flow` sets both turbine paths `OFF`.

`compile_derecho.sh` builds one local executable file:

```text
lesgo-run-exe
```

That executable is ignored by git because it is cluster-specific. The normal
Derecho sequence is:

```bash
cd test-cases/<case-name>
./compile_derecho.sh
qsub submit_derecho.pbs
```

## Changing MPI Or GPU Counts

For the GPU examples in this repository, the intended layout is one MPI rank
per GPU. Changing the GPU count therefore requires changing the same count in
three places:

| Setting | Where it appears | Meaning |
| --- | --- | --- |
| PBS resource line | `submit_derecho.pbs` | Number of nodes, MPI ranks per node, and GPUs per node requested from the scheduler. |
| `MPI_RANKS` / `MPI_PPN` | `submit_derecho.pbs` | Total MPI ranks and MPI ranks per node passed to `mpiexec`. |
| `nproc` | `lesgo.conf` | Number of LESGO MPI partitions expected by the input file. |

Example for 16 A100 GPUs on Derecho 4-GPU nodes:

```bash
#PBS -l select=4:ncpus=32:mpiprocs=4:ngpus=4:mem=220gb:gpu_type=a100
MPI_RANKS=16
MPI_PPN=4
```

The corresponding `lesgo.conf` value is:

```text
nproc=16
```

Do not commit generated run data here. Keep `output/`, `turbineOutput/`,
`runs/`, `grid.out`, `vel.out.c*`, queue logs, and compiled binaries out of
the repository.
