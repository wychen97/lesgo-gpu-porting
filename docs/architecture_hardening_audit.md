# Architecture Hardening Audit

This note records the readability and correctness guardrails added around the
optimized non-LVLSET GPU branch.  It is not a numerical validation report; use
the validation matrix and evidence ledger for runtime correctness and speedup
claims.

## Current Guardrails

The default readiness gate now checks these architecture risks:

- stale Fortran preprocessor symbols that are not emitted by CMake;
- unused `PP*` CMake definitions;
- bare single-symbol `#if` / `#elif` flag conditions, which should use
  `defined(...)`;
- module/program scopes missing top-level `implicit none`;
- `use ..., only:` imports from default-private modules that do not mark the
  imported symbol public;
- `use ..., only:` imports from tracked modules when the requested symbol is
  not actually exported, including lost re-exports after replacing broad
  module imports with explicit `only:` lists;
- duplicate names inside one `use ..., only:` list;
- repeated imports of the same local name from the same module within one
  module/subroutine/function scope;
- stale or duplicate explicit module `public` API names;
- accidental treatment of derived-type components or type-bound `public`
  procedures as module-level public symbols;
- unclassified broad `use module` imports that are neither narrowed nor
  documented as intentional compiler-interface/optional-module cases;
- IWM wall-surface indexing mistakes that collapse `(iwm_i, iwm_j)` into
  diagonal `(iwm_i, iwm_i)` indexing;
- ATM structure-on packed/gather paths accidentally falling back to
  structure-off-only behavior.

The synthetic self-test
`tools/check_fortran_interface_hygiene_selftest.py` covers the most subtle
Fortran interface rules without depending on the current production source
layout.

`docs/fortran_broad_import_audit.md` is a generated map of broad `use module`
imports.  The readiness gate allows classified broad imports that are tied to
MPI compiler interfaces, optional CGNS interfaces, or the deferred LVLSET
module, but it fails if a new broad import remains an unclassified cleanup
candidate.  Existing classified broad imports are not automatically wrong and
should not be rewritten in large mechanical patches.

## Issues Fixed In This Pass

- removed stale `ENABLE_CUDA` source paths before this hardening branch;
- removed the stale `OUTPUT_EXTRA` CMake/test-script path;
- normalized bare preprocessor branches such as `#elif PPATM` to
  `#elif defined(PPATM)`;
- removed duplicate `use ..., only:` imports in `divstress_w.f90`,
  `level_set.f90`, and `scalars.f90`;
- made `actuator_turbine_model.f90` import only the required `atm_base`
  symbols;
- added module-scope `implicit none` to modules that only had procedure-local
  `implicit none`;
- narrowed low-risk MPI imports in CFL/debug helpers and removed an unused MPI
  import from `messages.f90`;
- removed repeated scoped imports in `atm_lesgo_interface.f90`, `initial.f90`,
  and `forcing.f90`;
- narrowed ATM-side module imports so turbine metadata, actuator entry points,
  and MPI calls are imported from their owning modules explicitly;
- fixed the broad-import audit parser so `use module, only: ...` imports are not
  misclassified as broad imports;
- added a broad-import parser self-test so the generated audit keeps handling
  plain `use module`, `use module, only: ...`, and `use, intrinsic :: ...`
  forms distinctly;
- added a tracked-module import resolver after finding that narrowing `param`'s
  MPI import can otherwise drop legacy `MPI_COMM_WORLD` re-exports without a
  local compile failure;
- classified the remaining broad imports and made unclassified broad-import
  candidates a readiness-gate failure;
- regenerated static inventory and refactor-hotspot documentation after source
  changes.

## Required Evidence Before Merge

Run from the repository root:

```bash
python3 tools/check_branch_readiness.py --with-git-diff-check --with-cmake-configure --with-hit-cmake-configure
```

For solver behavior changes beyond readability/import/preprocessor cleanup,
also run the relevant benchmark smoke from `docs/gpu_validation_runbook.md`.
Before merging a cleaned branch, compile the current committed tree on Derecho
with the public case matrix and, when optional feature guards changed, the
optional non-LVLSET compile matrix.

## Limits

These checks are intentionally conservative source-maintenance checks.  They do
not prove:

- numerical equivalence of a solver change;
- GPU performance has not regressed;
- every broad `use module` import should be split immediately;
- large-file ownership boundaries are ready for file moves.

Use `docs/refactor_backlog.md` for future source splitting order and keep
behavioral changes separate from readability-only patches whenever possible.
