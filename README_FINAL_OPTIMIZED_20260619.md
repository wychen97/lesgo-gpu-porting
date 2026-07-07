# Final Optimized Non-LVLSET LESGO Version - 2026-06-19

This branch publishes the current final optimized non-LVLSET GPU source used on
Derecho after the all-module GPU-port work and SGS runtime-model completion.

Start here for collaborator handoff:

```text
docs/README.md
docs/collaborator_quickstart.md
docs/code_organization.md
docs/build_profiles.md
docs/cluster_scripts.md
docs/gpu_module_contracts.md
docs/gpu_port_coverage_audit.md
docs/gpu_static_full_inventory.md
docs/gpu_static_candidate_review.md
docs/gpu_validation_matrix.md
docs/gpu_validation_runbook.md
docs/gpu_release_objective_status.md
docs/lesgo_conf_gpu_validation_map.json
docs/lesgo_conf_key_validation_coverage.md
docs/source_file_inventory.md
docs/gpu_development_guidelines.md
docs/environment_switches.md
docs/refactor_backlog.md
docs/collaboration_readiness_status.md
docs/architecture_hardening_audit.md
docs/fortran_broad_import_audit.md
CONTRIBUTING.md
tools/README.md
```

The detailed reading order is maintained in `docs/README.md`.  Keep this
handoff block synchronized with that index.

Base branch:

```text
gpu-explicit-residency-wip
```

Derecho validation source:

```text
/glade/work/wchen/lesgo_versions/canonical/04_final_optimized_nonlvlset_20260619
```

Derecho final build:

```text
/glade/work/wchen/lesgo_versions/builds/04_final_optimized_nonlvlset_allmodules_20260619
```

Validated executable SHA256:

```text
fc88107031404604578bab5c9e149ba9ca3dfd52f187d5eb592fd74613a74e37
```

## Main Build Configuration

```text
USE_MPI=ON
USE_CPS=ON
USE_ATM=ON
USE_TURBINES=ON
USE_LES_GPU=ON
USE_GPU_AWARE_MPI=AUTO
USE_SCALARS=ON
USE_SCALARS_GPU=ON
USE_LVLSET=OFF
USE_HIT=OFF
USE_DYN_TN=OFF
USE_CGNS=OFF
LESGO_CUDA_CC=cc80
```

## Validation Summary

Final all-module smoke:

```text
/glade/work/wchen/lesgo_versions/benchmarks/04_final_optimized_allmodules_smoke_32_sgs5_n10_a100_20260619
```

Result:

```text
MPI_EXIT_STATUS=0
red/blue domains: Simulation complete
```

SGS runtime matrix:

```text
sgs=false: supported
sgs_model=1: supported
sgs_model=2: supported
sgs_model=3: supported
sgs_model=4: supported
sgs_model=5: supported
sgs_model=6/7: unsupported runtime values, guard-fail as expected
```

## Local Readability Checks

Run this before committing source or documentation edits:

```bash
python3 tools/check_branch_readiness.py
```

On a machine with CMake and a suitable compiler wrapper available, also run:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

The readiness wrapper runs these portable checks:

```bash
python3 tools/check_source_hygiene.py
python3 tools/check_cmake_option_docs.py
python3 tools/check_cmake_cache_docs.py
python3 tools/check_build_profile_args.py
python3 tools/check_cmake_comment_quality.py
python3 tools/check_cmake_source_groups.py
python3 tools/check_production_profile_docs.py
python3 tools/check_production_wording.py
python3 tools/check_active_source_comment_quality.py
python3 tools/check_cluster_script_docs.py
python3 tools/check_cluster_script_headers.py
python3 tools/check_test_case_script_docs.py
python3 tools/check_p0_paired_case_scripts.py
python3 tools/check_contributing_readiness_docs.py
python3 tools/check_handoff_readme.py
python3 tools/check_environment_switch_docs.py
python3 tools/check_lesgo_conf_coverage_docs.py
python3 tools/check_lesgo_conf_validation_map.py
python3 tools/check_lesgo_conf_key_validation_report.py
python3 tools/check_source_comment_hygiene.py
python3 tools/check_source_inventory.py
python3 tools/check_fortran_interface_hygiene.py
python3 tools/check_fortran_interface_hygiene_selftest.py
python3 tools/check_fortran_broad_import_audit.py
python3 tools/check_sgs_model_constants.py
python3 tools/check_refactor_backlog.py
python3 tools/check_navigation_maps.py
python3 tools/check_iwm_surface_indexing.py
python3 tools/check_atm_structure_packed_paths.py
python3 tools/check_gpu_comment_labels.py
python3 tools/check_gpu_headers.py
python3 tools/check_gpu_contract_source_groups.py
python3 tools/check_gpu_validation_matrix.py
python3 tools/check_gpu_benchmark_manifest.py
python3 tools/check_gpu_validation_plan.py
python3 tools/check_gpu_validation_evidence.py
python3 tools/check_gpu_timing_audit_consistency.py
python3 tools/check_gpu_validation_import_tools.py
python3 tools/check_p0_archived_evidence_importer.py
python3 tools/check_gpu_validation_runbook.py
python3 tools/check_gpu_release_objective_status_report.py
python3 tools/check_gpu_static_review.py
python3 tools/check_gpu_static_full_inventory_report.py
python3 tools/check_gpu_static_candidate_review_report.py
python3 tools/check_doc_paths.py
python3 tools/check_doc_source_refs.py
python3 tools/check_docs_index.py
python3 tools/check_tools_index.py
python3 tools/check_readiness_coverage.py
python3 tools/check_readiness_docs.py
```

The hygiene script scans tracked Fortran/CMake/docs/scripts for invalid UTF-8,
non-ASCII characters, non-LF line endings, tab characters, and trailing
whitespace.  It exists because this branch is edited through different
terminals and clusters where Unicode comments, CRLF churn, and tab indentation
have previously rendered inconsistently.
The CMake and environment-switch checks keep
`docs/build_profiles.md`, `docs/code_organization.md`, and
`docs/environment_switches.md` synchronized with root cache variables, build
cache-variable rows, build-option table rows, copy-paste profile `-D`
arguments, and solver source literals.
The CMake comment-quality check keeps build-file comments concrete enough for
collaborators to act on.
The CMake source-group check keeps root Fortran source files covered by exactly
one named `*_SOURCES` group in `CMakeLists.txt`.
The environment-switch check also requires every source `LESGO_*` switch to
appear in exactly one classified table row.
The lesgo.conf coverage check keeps `docs/gpu_port_coverage_audit.md`
synchronized with the non-LVLSET parser keys in `input_util.f90`.
The source-comment hygiene check keeps production source comments free of
stale personal tags, misspellings, and ad-hoc checkpoint wording.
The production-profile documentation check keeps the canonical production
`USE_*` settings synchronized across the README, build, and organization docs.
The production-wording check keeps validated production paths from being
described as experimental.
The active-source comment-quality check rejects vague maintenance comments in
non-LVLSET production source while leaving deferred LVLSET/tree files alone.
The cluster-script documentation check keeps root-level shell helpers classified
as current, optional, validation-only, historical, or legacy.
The cluster-script header check makes that status visible inside each root
shell helper.
The test-case script documentation check requires each script-heavy
`test-cases/<case>/` tree to have a local README that explains its launchers.
The contributing-readiness documentation check keeps `CONTRIBUTING.md` pointed
at the canonical readiness command and script-list source of truth.
The handoff README check keeps the top-level handoff block synchronized with
`docs/README.md`.
The source-inventory check keeps `docs/source_file_inventory.md` synchronized
with tracked Fortran files, with one table row per source file and valid
`Build path` labels.
The SGS model-constant check keeps active SGS dispatch code on named
`SGS_MODEL_*` constants instead of bare integer comparisons.
The refactor-backlog check keeps `docs/refactor_backlog.md` synchronized with
the current largest non-LVLSET source hotspots.
The navigation-map check requires large active Fortran files to expose a
header-level map for collaborators.
The IWM surface-indexing check protects the point-local `(i,j)` wall-model
state update convention that avoids diagonal-only surface writes.
The GPU comment-label check keeps production GPU source comments free of stale
internal optimization-phase labels.
The GPU-header check requires GPU-specific Fortran files to explain ownership
or their data-movement/device-residency contract near the top of the file.
The GPU contract source-group check keeps `docs/gpu_module_contracts.md`
aligned with the named CMake `*_SOURCES` groups.
The GPU validation matrix check keeps `docs/gpu_validation_matrix.md`
explicit about which non-LVLSET paths have paired CPU/GPU speed evidence and
which still need benchmarks.
The GPU benchmark manifest check keeps `docs/gpu_benchmark_manifest.json`
aligned with the validation matrix and public test-case paths.
The documentation path check catches stale `docs/...` and `tools/...`
references in Markdown handoff files.
The documented-source check catches stale Fortran source filenames in Markdown
handoff files, including wildcard source groups.
The documentation-index check keeps `docs/README.md` synchronized with tracked
collaborator docs.
The tooling-index check keeps `tools/README.md` synchronized with tracked
top-level tool files.
The readiness-coverage check keeps tracked readiness scripts in the default
wrapper.
The readiness-documentation check keeps this checklist synchronized with the
wrapper itself.

`git diff --check` is still useful with a single native Git client.  It remains
optional because shared Windows/WSL checkouts can report CRLF-only differences
inconsistently; the default source-hygiene check now directly enforces LF text
for tracked human-edited files.

## Notes

LVLSET remains intentionally deferred from the optimized production path.

The ATM PPLES reduce-to-operating-ranks load-balance branch is not included as
the default final path.  It remains experimental until a larger real
multi-rank/multi-turbine case shows a robust benefit.

The IWM filtered friction velocity update was corrected after publication:
`iwm_flt_us(i,j)` is now updated point-locally instead of repeatedly writing
the diagonal entry `iwm_flt_us(i,i)`.

## Repository Hygiene

The source branch should stay clean after a clone.  Local benchmark launch
files, submitted-job records, monitoring directories, porting patches, and
large generated binaries are ignored by `.gitignore`.  Keep durable benchmark
results in a documented validation note or an external run directory; do not
commit exploratory cluster artifacts into the source tree.

Line endings are controlled by `.gitattributes`: source, docs, scripts, CMake,
and case input files are LF text.  Generated numerical or visual artifacts are
marked binary.  This avoids false diffs when the same branch is touched from
Windows, WSL, Derecho, Delta, and other HPC login nodes.
