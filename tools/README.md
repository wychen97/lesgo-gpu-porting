# LESGO Tooling Index

This directory contains lightweight source-maintenance tools for the optimized
non-LVLSET GPU branch.  They are intentionally plain Python or small validation
sources so they can run on Windows/WSL checkouts and HPC login nodes.

Use the wrapper for normal collaboration checks:

```bash
python3 tools/check_branch_readiness.py
```

The wrapper groups checks by concern: source/build hygiene, helper-script and
handoff documentation, Fortran/GPU source contracts, then index self-checks.
Place new readiness checks in the closest existing group so the default gate
stays easy to audit.

When CMake and a suitable compiler wrapper are available, use:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

For the optional HIT plus LES GPU source set, use:

```bash
python3 tools/check_branch_readiness.py --with-hit-cmake-configure
```

## Readiness Checks

| Script | Purpose |
| --- | --- |
| `tools/check_active_source_comment_quality.py` | Verifies active non-LVLSET source comments avoid vague maintenance wording. |
| `tools/check_atm_structure_packed_paths.py` | Verifies ATM structure-on runs keep packed/batched paths and carry `Cm`/`pitchingMoment`. |
| `tools/check_branch_readiness.py` | Runs the default collaboration-readiness gate. |
| `tools/check_build_profile_args.py` | Verifies documented CMake profile `-D...` arguments are valid public CMake knobs. |
| `tools/check_cluster_script_headers.py` | Verifies root-level shell helpers carry a visible status header. |
| `tools/check_cluster_script_docs.py` | Verifies root-level cluster/helper shell scripts are documented and classified. |
| `tools/check_cmake_cache_docs.py` | Verifies public CMake cache variables have table rows in `docs/build_profiles.md`. |
| `tools/check_cmake_comment_quality.py` | Rejects vague maintenance comments in tracked CMake build files. |
| `tools/check_cmake_option_docs.py` | Verifies root `USE_*` CMake options have build-option table rows in `docs/code_organization.md`. |
| `tools/check_cmake_source_groups.py` | Verifies tracked root Fortran files appear in named CMake source groups. |
| `tools/check_contributing_readiness_docs.py` | Verifies `CONTRIBUTING.md` points to the readiness source of truth. |
| `tools/check_doc_paths.py` | Verifies documented local `docs/...` and `tools/...` references still exist. |
| `tools/check_doc_source_refs.py` | Verifies documented Fortran source references point to tracked files. |
| `tools/check_docs_index.py` | Verifies `docs/README.md` indexes tracked collaborator docs. |
| `tools/check_environment_switch_docs.py` | Verifies `LESGO_*` runtime switches used in source are documented. |
| `tools/check_fortran_broad_import_audit.py` | Verifies the generated broad Fortran import audit is current. |
| `tools/check_fortran_broad_imports_selftest.py` | Regression-tests the broad Fortran import audit parser on synthetic `use` statements. |
| `tools/check_fortran_interface_hygiene.py` | Verifies Fortran preprocessor symbols, module/program `implicit none`, and `use ..., only:` import/export boundaries. |
| `tools/check_fortran_interface_hygiene_selftest.py` | Regression-tests the Fortran interface hygiene checker on synthetic module snippets. |
| `tools/check_gpu_benchmark_manifest.py` | Verifies the benchmark manifest covers every GPU validation matrix row and uses valid, internally consistent CMake build settings. |
| `tools/check_gpu_comment_labels.py` | Verifies production GPU source comments do not use stale internal optimization labels. |
| `tools/check_gpu_contract_source_groups.py` | Verifies GPU contracts document every named CMake source group. |
| `tools/check_gpu_headers.py` | Verifies GPU-specific Fortran files explain ownership or data movement near the file header. |
| `tools/check_gpu_static_candidate_review_report.py` | Verifies the generated full static GPU candidate review report is current. |
| `tools/check_gpu_static_full_inventory_report.py` | Verifies the generated full static GPU subprogram inventory report is current. |
| `tools/check_gpu_release_objective_status_report.py` | Verifies the generated release-objective status report is current. |
| `tools/check_gpu_static_review.py` | Verifies static GPU review bucket counts and bucket-to-validation-row mappings in the audit match the scanner output. |
| `tools/check_gpu_validation_evidence.py` | Verifies the GPU validation evidence ledger is consistent and does not imply unsupported speedup claims. |
| `tools/check_gpu_validation_import_tools.py` | Smoke-tests parsing and importing LESGO timing logs into a temporary evidence ledger. |
| `tools/check_gpu_validation_matrix.py` | Verifies the non-LVLSET GPU validation matrix has required rows and status labels. |
| `tools/check_gpu_validation_plan.py` | Verifies the generated CPU/GPU validation run plan covers open benchmark rows. |
| `tools/check_gpu_validation_runbook.py` | Verifies the GPU validation runbook covers every benchmark-manifest row. |
| `tools/check_gpu_timing_audit_consistency.py` | Verifies timing numbers in the GPU coverage audit match the validation evidence ledger. |
| `tools/check_handoff_readme.py` | Verifies the top-level handoff README points to the current doc entrypoints. |
| `tools/check_iwm_surface_indexing.py` | Verifies IWM wall-surface arrays keep point-local `(iwm_i, iwm_j)` indexing. |
| `tools/check_lesgo_conf_coverage_docs.py` | Verifies the GPU coverage audit tracks non-LVLSET `lesgo.conf` parser keys. |
| `tools/check_lesgo_conf_key_validation_report.py` | Verifies the generated key-level `lesgo.conf` validation coverage report is current. |
| `tools/check_lesgo_conf_validation_map.py` | Verifies parsed `lesgo.conf` groups are mapped to GPU validation rows. |
| `tools/check_navigation_maps.py` | Verifies large active Fortran modules have a header-level navigation map. |
| `tools/check_p0_archived_evidence_importer.py` | Smoke-tests importing archived public p0 CPU/GPU logs into a temporary validation evidence ledger. |
| `tools/check_p0_paired_case_scripts.py` | Verifies the public p0 cases preserve clean paired CPU/GPU build and submit profiles, including executable names derived from root CMake target suffixes. |
| `tools/check_production_profile_docs.py` | Verifies canonical production `USE_*` settings match across build and organization docs. |
| `tools/check_production_wording.py` | Verifies validated production paths are not described as experimental. |
| `tools/check_refactor_backlog.py` | Verifies `docs/refactor_backlog.md` lists the current non-LVLSET source-size hotspots. |
| `tools/check_readiness_coverage.py` | Verifies tracked readiness checks are included in the wrapper. |
| `tools/check_readiness_docs.py` | Verifies the README readiness checklist matches the wrapper. |
| `tools/check_sgs_model_constants.py` | Verifies active SGS dispatch uses named `SGS_MODEL_*` constants. |
| `tools/check_source_comment_hygiene.py` | Verifies source comments do not use stale personal tags or ad-hoc wording. |
| `tools/check_source_hygiene.py` | Verifies tracked human-edited source is UTF-8, ASCII-only, LF-only, spaces-indented, and free of trailing whitespace. |
| `tools/check_test_case_script_docs.py` | Verifies script-heavy test-case trees have local README guidance. |
| `tools/check_source_inventory.py` | Verifies `docs/source_file_inventory.md` matches tracked Fortran sources. |
| `tools/check_tools_index.py` | Verifies this file indexes tracked top-level tool files. |

## Other Validation Utilities

| File | Purpose |
| --- | --- |
| `tools/cmake_metadata.py` | Shared root CMake option/cache-variable parsers for readiness checks. |
| `tools/fortran_inventory.py` | Shared tracked Fortran-source discovery helpers for readiness checks. |
| `tools/gpu_validation_plan.py` | Expands the benchmark manifest and evidence ledger into concrete paired CPU/GPU run tasks, import commands, and public-case Derecho submit templates. |
| `tools/import_lesgo_timing_evidence.py` | Imports a standard LESGO timing log into one validation evidence row or a named runtime variant; speedup claims require explicit correctness checks. |
| `tools/import_lesgo_timing_pair.py` | Imports matched CPU/GPU LESGO logs and records a faster or not-faster paired result; `--compare-diagnostics` stores parsed MPI/divergence/kinetic-energy checks, and faster claims require explicit correctness checks plus `--evidence-item` coverage for manifest requirements. |
| `tools/import_p0_archived_evidence.py` | Imports archived public p0 CPU/GPU runs from `run-archives/<label>/lesgo_<label>.log` into the paired validation evidence ledger. |
| `tools/parse_lesgo_timing.py` | Parses LESGO stdout/stderr timing blocks into JSON for validation evidence import. |
| `tools/compare_scalar_checkpoint.py` | Compares CPU/GPU `scal.out.c*` scalar checkpoints and reports scalar-field differences. |
| `tools/prepare_sgs_matrix_cases.py` | Creates isolated compact CPU/GPU channel-flow case directories for `sgs=false`, `sgs_model=1..5`, and optional `USE_DYN_TN=ON` validation. |
| `tools/prepare_scalar_cases.py` | Creates isolated compact CPU/GPU channel-flow case directories for passive and active scalar validation. |
| `tools/readiness_manifest.py` | Shared parser for the readiness wrapper's `PYTHON_CHECKS` manifest. |
| `tools/report_lesgo_conf_key_validation.py` | Generates the key-level `lesgo.conf` validation coverage report from source parser keys and evidence state. |
| `tools/report_gpu_matrix_status_updates.py` | Suggests validation-matrix status updates implied by the evidence ledger. |
| `tools/report_gpu_release_objective_status.py` | Generates the release-objective status report tying static source coverage to validation evidence. |
| `tools/report_gpu_static_candidate_review.py` | Generates the full static unmarked-candidate GPU review report. |
| `tools/report_gpu_static_full_inventory.py` | Generates the full static GPU classification inventory for every tracked non-LVLSET subprogram. |
| `tools/report_fortran_broad_imports.py` | Generates the broad Fortran import audit for future import-boundary cleanup. |
| `tools/script_inventory.py` | Shared tracked-script discovery helpers for readiness checks. |
| `tools/report_gpu_static_inventory.py` | Reports static non-LVLSET Fortran subprogram GPU-marker coverage and review buckets for audits. |
| `tools/report_gpu_validation_gaps.py` | Reports open GPU validation and paired-speedup evidence gaps by runbook batch. |
| `tools/update_gpu_validation_evidence.py` | Updates one CPU or GPU timing row in the validation evidence ledger. |
| `tools/require_gpu_release_objective.py` | Strict release gate that fails until every non-LVLSET validation surface has paired correctness evidence and an acceptable GPU-faster result. |
| `tools/validate_filt_da_cufft.F90` | Small CUDA Fortran validation source for `filt_da`/cuFFT behavior. |

## Maintenance Rule

When adding a new top-level tool file:

1. Add it to the table above.
2. If it is a readiness check, add it to `PYTHON_CHECKS` in
   `tools/check_branch_readiness.py`.
3. If it is a readiness check, add it to the README readiness checklist.
4. Run `python3 tools/check_branch_readiness.py`.
