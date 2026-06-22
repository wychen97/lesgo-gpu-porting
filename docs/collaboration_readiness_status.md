# Collaboration Readiness Status

This note summarizes the current readability and maintainability state of the
optimized non-LVLSET GPU branch.  It is not a validation report for new physics
or performance claims; it is a handoff aid for collaborators editing the code.

## Protected Invariants

The default readiness wrapper now protects these collaboration invariants:

- tracked source, docs, scripts, and CMake files use stable text formatting;
- public CMake options, build profiles, runtime switches, source inventory, and
  documentation indexes stay synchronized;
- validated production paths are not described as experimental;
- root-level cluster/helper scripts stay classified in `docs/cluster_scripts.md`
  and carry visible status headers;
- script-heavy `test-cases/<case>/` trees provide local README launcher maps;
- active non-LVLSET source comments do not use vague maintenance wording such
  as "for now", "temporary", or "fix later";
- production source comments do not use stale personal tags or ad-hoc IWM
  checkpoint wording;
- production GPU comments do not use stale internal optimization-phase labels;
- GPU-specific source files state ownership or data-movement contracts near the
  file header;
- large active Fortran files expose header-level navigation maps;
- the IWM wall-surface state keeps point-local `(iwm_i, iwm_j)` semantics.

Run:

```bash
python3 tools/check_branch_readiness.py
```

When a compiler/CMake environment is available, also run:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

## Current Production Scope

The organized production path is the optimized non-LVLSET GPU configuration
documented in `README_FINAL_OPTIMIZED_20260619.md` and
`docs/code_organization.md`.

The active production modules are:

- LES GPU core;
- SGS runtime models `1..5`, plus disabled SGS;
- Lagrangian dynamic SGS GPU path;
- pressure and tridiagonal GPU/MPI paths;
- scalar GPU transport;
- CPS/concurrent precursor;
- ATM, turbine forcing, induced-velocity correction, and structural coupling;
- diagnostics, output, restart, and statistics paths needed by the validated
  benchmark cases.

## Intentionally Deferred

These areas remain in the repository but are not proof of production readiness:

- LVLSET and `trees_*` sources;
- ATM point-owner load balancing as a default path;
- replicated global pressure tridiagonal solve as a production default;
- broad source-file splitting without a validation plan.

Deferred code may still contain legacy comments.  The active-source comment
quality gate intentionally excludes LVLSET/tree files because they are outside
the optimized non-LVLSET production scope.

## Editing Rule

For readability-only changes, keep edits small and prove that the readiness
wrapper still passes.  For solver or ownership changes, use the validation
matrix in `docs/gpu_module_contracts.md`; local documentation checks alone are
not enough for GPU hot-path changes.
