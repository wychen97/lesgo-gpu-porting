# LESGO Documentation Index

This directory is the collaborator-facing map for the optimized non-LVLSET GPU
branch.  Read these files before changing solver code.

## Reading Order

1. `collaborator_quickstart.md`
   - fast orientation for CPU-version users;
   - first commands and standard build profiles;
   - current validation scope and known remaining cleanup.

2. `code_organization.md`
   - source-file map;
   - production build options;
   - production timestep order;
   - device ownership boundaries.

3. `build_profiles.md`
   - recommended CMake configure profiles;
   - user-facing cache variables;
   - CPU baseline and configure-smoke profiles.

4. `cluster_scripts.md`
   - root-level cluster helper script map;
   - current, optional, validation, historical, and legacy script labels.

5. `gpu_module_contracts.md`
   - per-module ownership contracts;
   - CPU/GPU boundary rules;
   - validation expectations by change area.

6. `source_file_inventory.md`
   - every tracked Fortran source file;
   - build-option grouping;
   - source-tree entry points for optional modules.

7. `gpu_development_guidelines.md`
   - rules for GPU changes;
   - synchronization and data-residency expectations;
   - review checklist for GPU kernels and MPI paths.

8. `environment_switches.md`
   - documented `LESGO_*` runtime switches;
   - production, diagnostic, and benchmark-only switch classification.

9. `refactor_backlog.md`
   - safe order for future readability refactors;
   - work intentionally deferred from the published source branch.

10. `collaboration_readiness_status.md`
   - current readability/readiness guard summary;
   - production-scope and deferred-work status;
   - editing rule for readability-only versus solver changes.

## Historical References

- `gpu_port_refactor_history.md`: historical architecture plan that informed
  the explicit-residency branch.  It is retained for rationale only; current
  build and validation guidance lives in the reading-order files above.
- `gpu_port_refactor_summary.html`: rendered historical summary of the same
  GPU-port refactor work.

## Contributor Entry Points

- `../README_FINAL_OPTIMIZED_20260619.md`: published branch summary and
  validation snapshot.
- `../CONTRIBUTING.md`: contribution workflow and required checks.
- `../tools/README.md`: local readiness and validation tooling map.
- `../tools/check_branch_readiness.py`: default local readiness gate.

## Production Scope

The validated production branch is the optimized non-LVLSET GPU path.  LVLSET
files remain in the repository for completeness, but LVLSET is outside the
current optimized production scope and should not be used as validation evidence
for GPU changes.

The major enabled production modules are:

- LES GPU core;
- SGS runtime models, excluding unsupported runtime values;
- CPS/concurrent precursor;
- scalar transport with the GPU scalar path;
- turbine and ATM/structural coupling;
- diagnostics and output paths needed by the validated benchmark cases.

## Keeping This Directory Current

When adding a new module-level contract, runtime switch, or required validation
step, update this index and run:

```bash
python3 tools/check_branch_readiness.py
```
