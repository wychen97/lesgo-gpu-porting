# Test Cases

This directory follows the upstream LESGO convention of keeping runnable cases
under `test-cases/`. The maintained examples are:

| Case | Purpose | Main CMake switches | Latest validation record |
| --- | --- | --- | --- |
| `channel_flow` | LES core baseline without turbine physics | `USE_LES_GPU=ON`, `USE_TURBINES=OFF`, `USE_ATM=OFF` | Derecho paired CPU/GPU pass, `240^3`, 200 steps |
| `adm_disk` | Actuator disk model increment | `USE_LES_GPU=ON`, `USE_TURBINES=ON`, `USE_ATM=OFF` | Derecho paired CPU/GPU pass, `240^3`, 200 steps |
| `atm_line` | Actuator line / 5 MW turbine path | `USE_LES_GPU=ON`, `USE_TURBINES=OFF`, `USE_ATM=ON` | Derecho rigid and structural CPU/GPU passes, `240^3` |
| `large_windfarm_3072x384x400_60turbines` | Large 60-turbine benchmark setup | `USE_LES_GPU=ON`, `USE_ATM=ON` | Derecho paired CPU50/GPU16 pass, 50 steps |
| `level_set_cubes` | Immersed-surface Level Set CPU/GPU parity | `USE_LES_GPU=ON`, `USE_LVLSET=ON`, `USE_LVLSET_GPU=ON` | Derecho and Delta matrix/restart passes, 58 runtime tasks |
| `scalar_transport` | Passive scalar and active temperature/buoyancy coupling | `USE_SCALARS=ON`, `USE_SCALARS_GPU=ON` for GPU | Derecho compact CPU/GPU passive and active passes, `64^3` |
| `concurrent_precursor` | Red/blue concurrent precursor exchange, with optional scalar coupling | `USE_CPS=ON`, optionally `USE_SCALARS=ON` | Derecho compact CPU/GPU velocity and scalar passes, two `64^3` domains |
| `inflow_and_forcing` | HIT, shifted inflow, Coriolis, and sponge forcing | `USE_HIT=ON` | Derecho compact CPU/GPU passes for all three variants, `64^3` |

The compact `64^3` records are execution and numerical-parity checks, not
production scaling or statistically converged turbulence claims. Detailed
performance and acceptance evidence remains under `docs/`.

Each case contains:

- `lesgo.conf` input settings;
- model input files when needed, such as `inputATM/`;
- `README.md` with the case-specific module setup;
- `compile_derecho.sh`, which compiles a selected CPU or GPU profile;
- `submit_derecho.pbs`, which only submits/runs the case on Derecho.

The root `CMakeLists.txt` declares the available feature switches and their
neutral defaults. The case-specific feature choices live in
`compile_derecho.sh`; for example, `adm_disk` sets `USE_TURBINES=ON`, while
`channel_flow` sets both turbine paths `OFF`.

`compile_derecho.sh` installs a profile-specific local executable, normally:

```text
lesgo-run-exe-cpu
lesgo-run-exe-gpu
```

That executable is ignored by git because it is cluster-specific. The normal
Derecho sequence is:

```bash
cd test-cases/<case-name>
./compile_derecho.sh gpu
qsub submit_derecho.pbs
```

Cases with multiple feature variants document their second compile argument or
`CASE_VARIANT` submission value in the case README. The scripts keep runtime
outputs under `runs/` or `run-archives/` and do not overwrite tracked inputs.

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

The CPS case is the decomposition exception: two equal red/blue domains use an
even MPI world. With the default two-rank job, each domain has one local rank,
so each input contains `nproc=1` while the scheduler launches two ranks total.
