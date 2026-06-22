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
| `tools/check_gpu_comment_labels.py` | Verifies production GPU source comments do not use stale internal optimization labels. |
| `tools/check_gpu_contract_source_groups.py` | Verifies GPU contracts document every named CMake source group. |
| `tools/check_gpu_headers.py` | Verifies GPU-specific Fortran files explain ownership or data movement near the file header. |
| `tools/check_handoff_readme.py` | Verifies the top-level handoff README points to the current doc entrypoints. |
| `tools/check_iwm_surface_indexing.py` | Verifies IWM wall-surface arrays keep point-local `(iwm_i, iwm_j)` indexing. |
| `tools/check_navigation_maps.py` | Verifies large active Fortran modules have a header-level navigation map. |
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
| `tools/readiness_manifest.py` | Shared parser for the readiness wrapper's `PYTHON_CHECKS` manifest. |
| `tools/script_inventory.py` | Shared tracked-script discovery helpers for readiness checks. |
| `tools/validate_filt_da_cufft.F90` | Small CUDA Fortran validation source for `filt_da`/cuFFT behavior. |

## Maintenance Rule

When adding a new top-level tool file:

1. Add it to the table above.
2. If it is a readiness check, add it to `PYTHON_CHECKS` in
   `tools/check_branch_readiness.py`.
3. If it is a readiness check, add it to the README readiness checklist.
4. Run `python3 tools/check_branch_readiness.py`.
