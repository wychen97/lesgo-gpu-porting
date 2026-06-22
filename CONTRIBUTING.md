# Contributing To This LESGO Branch

This branch is the shared optimized non-LVLSET GPU branch.  It is intended to
stay usable by collaborators who did not participate in the GPU port, so source
changes should be small, documented, and validated before they are pushed.

Start with the documentation index:

```text
docs/README.md
docs/collaborator_quickstart.md
```

## Before Editing

1. Read `docs/collaborator_quickstart.md` for the CPU-user orientation,
   standard build profiles, and current validation scope.
2. Read `docs/code_organization.md` for the source-file map and production
   timestep order.
3. Read `docs/gpu_module_contracts.md` for the ownership boundary of the module
   you are changing.
4. Check `docs/environment_switches.md` before adding or removing runtime
   environment variables.
5. Check `docs/refactor_backlog.md` before starting readability-only cleanup, so
   refactors happen in a safe order.

## Change Rules

- Keep the production non-LVLSET path as the default reference path.
- Do not use LVLSET behavior as proof that the optimized GPU path is validated.
- Keep CPU fallback behavior correct when changing a GPU module.
- Avoid adding new environment switches unless the setting is genuinely needed
  for production, diagnostics, or benchmark isolation.
- Prefer explicit module ownership comments over broad rewrites.
- Keep local benchmark launch scripts, submitted-job records, and generated
  binaries out of the source tree.

## Required Local Checks

Run this before committing source, documentation, CMake, or script edits:

```bash
python3 tools/check_branch_readiness.py
```

When CMake and a suitable compiler wrapper are available, also run:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

The exact portable script list and guard descriptions are maintained in
`README_FINAL_OPTIMIZED_20260619.md` under "Local Readability Checks".
`tools/README.md` is the index for the individual readiness tools.  Keep those
files as the source of truth instead of duplicating the full checklist here.

## Validation Expectations

Use the validation matrix in `docs/gpu_module_contracts.md` to decide which
smoke test or benchmark is required.  A local syntax or documentation check is
not enough for GPU hot-path changes.

At minimum, record:

- build options;
- compiler and cluster;
- case size and number of steps;
- whether ATM, CPS, scalars, and SGS were enabled;
- final solver status;
- relevant component timing or physics checks.

## Commit Hygiene

Use focused commits.  Separate mechanical cleanup from behavior changes, and
separate documentation updates from solver changes when possible.  If a commit
changes an ownership boundary, update `docs/gpu_module_contracts.md` in the
same commit.
