# Cluster Scripts

The repository root intentionally has no cluster helper `*.sh` scripts.

Build and submission examples live inside the individual cases under
`test-cases/`:

```text
test-cases/<case>/compile_derecho.sh
test-cases/<case>/submit_derecho.pbs
```

This keeps the public package easy to navigate:

- root: solver source, CMake, and documentation;
- `test-cases/`: case-specific input files and case-specific cluster examples.

## Derecho Workflow

From any case directory:

```bash
./compile_derecho.sh
qsub submit_derecho.pbs
```

`compile_derecho.sh` only compiles and installs `lesgo-run-exe`.
`submit_derecho.pbs` only submits/runs the already-built executable.

## Porting Rule

When moving to another cluster, copy the two case-local files and edit:

1. module loads;
2. compiler/MPI wrapper, usually `FC=ftn`;
3. GPU-aware MPI setting;
4. scheduler resource request;
5. final `mpiexec` or `srun` launch command.
