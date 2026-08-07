# LESGO Documentation Index

This directory is the collaborator-facing map for the optimized GPU branch,
including the opt-in Level Set GPU path. Read these files before changing
solver code.

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

6. `gpu_port_coverage_audit.md`
   - current GPU-port coverage by build/runtime path;
   - validation strength for each non-LVLSET module;
   - remaining optional-path benchmark targets.

7. `gpu_static_full_inventory.md`
   - generated full static classification of every tracked non-LVLSET
     Fortran subprogram;
   - separates GPU-marked, host-boundary, host/diagnostic, and candidate code.

8. `gpu_static_candidate_review.md`
   - generated full list of static unmarked runtime candidates;
   - maps each candidate to a review bucket and responsible validation rows.

9. `gpu_validation_matrix.md`
   - release-facing non-LVLSET validation checklist;
   - separates source coverage, correctness, and CPU/GPU speedup claims;
   - tracks which paths still need paired benchmark evidence.

10. `gpu_validation_runbook.md`
   - benchmark batches for closing missing evidence rows;
   - minimum evidence fields for paired CPU/GPU speedup claims;
   - completion rule for updating validation status.

11. `gpu_release_objective_status.md`
   - generated status tying the static function inventory to the validation
     evidence ledger;
   - explicit release-objective pass/fail state for the non-LVLSET GPU port.

12. `lesgo_conf_gpu_validation_map.json`
   - machine-readable map from parsed `lesgo.conf` groups to validation rows;
   - bridge between runtime configuration surface and benchmark evidence.

13. `lesgo_conf_key_validation_coverage.md`
   - generated key-level expansion of the `lesgo.conf` validation map;
   - shows responsible validation rows and current evidence state for each
     parsed non-LVLSET key.

14. `source_file_inventory.md`
   - every tracked Fortran source file;
   - build-option grouping;
   - source-tree entry points for optional modules.

15. `gpu_development_guidelines.md`
   - rules for GPU changes;
   - synchronization and data-residency expectations;
   - review checklist for GPU kernels and MPI paths.

16. `environment_switches.md`
   - documented `LESGO_*` runtime switches;
   - production, diagnostic, and benchmark-only switch classification.

17. `refactor_backlog.md`
   - safe order for future readability refactors;
   - work intentionally deferred from the published source branch.

18. `collaboration_readiness_status.md`
   - current readability/readiness guard summary;
   - production-scope and deferred-work status;
   - editing rule for readability-only versus solver changes.

19. `architecture_hardening_audit.md`
   - source-maintenance guardrails added for this hardening pass;
   - classes of similar mistakes now checked automatically;
   - limits of the source-level checks versus numerical validation.

20. `architecture_performance_review_20260727.md`
   - current ATM ownership and optional-path review;
   - measured baseline/candidate performance;
   - correctness, build, and restart acceptance results.

21. `fortran_broad_import_audit.md`
   - generated map of broad `use module` imports without `only:` lists;
   - readiness-gated guide for import-boundary cleanup: unclassified broad
     imports must be narrowed or explicitly classified;
   - helps keep readability refactors small and reviewable.

22. `restart_state_contract.md`
   - synchronized LES, SGS, ATM, controller, and structural restart state;
   - CPU/GPU device-ownership requirements at restart;
   - continuous-versus-split validation protocol.

23. `level_set_gpu_port.md`
   - Level Set CPU/GPU ownership and MPI communication contract;
   - SGS and optional-path validation matrix;
   - build flags, restrictions, and measured correctness/performance evidence.

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

The canonical turbine production profile remains the optimized non-Level-Set
GPU path. An additional opt-in Level Set GPU profile is now implemented and
validated separately. Evidence from one profile must not be presented as
coverage for the other.

The major enabled production modules are:

- LES GPU core;
- SGS runtime models, excluding unsupported runtime values;
- CPS/concurrent precursor;
- scalar transport with the GPU scalar path;
- turbine and ATM/structural coupling;
- immersed-surface Level Set physics when `USE_LVLSET_GPU=ON`;
- diagnostics and output paths needed by the validated benchmark cases.

## Keeping This Directory Current

When adding a new module-level contract, runtime switch, or required validation
step, update this index and run:

```bash
python3 tools/check_branch_readiness.py
```
