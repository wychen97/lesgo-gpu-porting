# Test Cases

This directory follows the upstream LESGO convention of keeping runnable cases
under `test-cases/`. This branch is the CPU-default counterpart to the GPU
porting branch. Only the four cases used in the current validation suite are
included:

| Case | Purpose | Main CMake switches |
| --- | --- | --- |
| `channel_flow` | LES core baseline without turbine physics | `USE_CPU_BUILD=ON`, `USE_LES_GPU=OFF`, `USE_TURBINES=OFF`, `USE_ATM=OFF` |
| `adm_disk` | Actuator disk model increment | `USE_CPU_BUILD=ON`, `USE_LES_GPU=OFF`, `USE_TURBINES=ON`, `USE_ATM=OFF` |
| `atm_line` | Actuator line / 5 MW turbine path | `USE_CPU_BUILD=ON`, `USE_LES_GPU=OFF`, `USE_TURBINES=OFF`, `USE_ATM=ON` |
| `large_windfarm_3072x384x400_60turbines` | Large 60-turbine benchmark setup | `USE_CPU_BUILD=ON`, `USE_LES_GPU=OFF`, `USE_ATM=ON` |

Each case contains:

- `lesgo.conf` input settings;
- model input files when needed, such as `inputATM/`;
- `README.md` with the case-specific module setup;
- `compile_derecho.sh`, which compiles the selected CPU-default executable;
- `submit_derecho.pbs`, which only submits/runs the case on Derecho.

The root `CMakeLists.txt` declares the available feature switches and their
neutral defaults. The case-specific feature choices live in
`compile_derecho.sh`; for example, `adm_disk` sets `USE_TURBINES=ON`, while
`channel_flow` sets both turbine paths `OFF`.

`compile_derecho.sh` builds one local executable file by default:

```text
lesgo-run-exe-cpu
```

It also updates the compatibility copy `lesgo-run-exe`.

That executable is ignored by git because it is cluster-specific. The normal
Derecho sequence is:

```bash
cd test-cases/<case-name>
./compile_derecho.sh
qsub submit_derecho.pbs
```

## Changing MPI Counts

For the CPU examples in this branch, changing the MPI rank count requires
changing the same count in three places:

| Setting | Where it appears | Meaning |
| --- | --- | --- |
| PBS resource line | `submit_derecho.pbs` | Number of nodes, CPU cores, MPI ranks per node, and memory requested from the scheduler. |
| `MPI_RANKS` / `MPI_PPN` | `submit_derecho.pbs` | Total MPI ranks and MPI ranks per node passed to `mpiexec`. |
| `nproc` | `lesgo.conf` | Number of LESGO MPI partitions expected by the input file. |

Example for 16 CPU ranks split across two Derecho CPU nodes:

```bash
#PBS -l select=2:ncpus=8:mpiprocs=8:mem=180gb
MPI_RANKS=16
MPI_PPN=8
```

The corresponding `lesgo.conf` value is:

```text
nproc=16
```

Do not commit generated run data here. Keep `output/`, `turbineOutput/`,
`runs/`, `grid.out`, `vel.out.c*`, queue logs, and compiled binaries out of
the repository.
