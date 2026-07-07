# GPU Development Guidelines

This file records the rules for keeping the optimized GPU branch readable and
safe for shared development.

## General Rules

- Keep production changes in the existing LESGO modules unless a new module
  removes real complexity.
- Keep CPU fallback semantics visible.  If a GPU path replaces a CPU path,
  preserve the CPU path until validation proves the GPU path is correct.
- Avoid adding environment-variable switches for production behavior.  Use
  CMake options or documented runtime parameters.
- Any remaining environment-variable control must be documented in
  `docs/environment_switches.md`.
- Remove switches after they are proven useless.  Leave a short comment only
  when the removed behavior explains a design choice future developers may
  question.
- Do not optimize LVLSET in this branch unless it becomes a production
  requirement.

## GPU Residency

- Prefer explicit device residency for hot arrays.
- Host/device updates should happen only at clear CPU/GPU boundaries:
  initialization, I/O, restart, diagnostics, and unavoidable MPI staging.
- If an array is used inside the timestep hot path, document whether the CPU,
  GPU, or both own the authoritative value at that point.
- Avoid hidden full-field `update self` or `update device` calls inside
  frequently executed routines.

## MPI and Communication

- `USE_GPU_AWARE_MPI=AUTO` is the default safe production choice.
- `PPGPU_AWARE_MPI` means MPI receives device pointers through
  `host_data use_device`.
- If GPU-aware MPI is unavailable, host-staged fallbacks must still compile and
  run, but timing should be reported separately.
- Do not compare GPU-aware and host-staged timings as if they are the same
  implementation.

## Timing

- Do not judge performance from a single timestep when output, statistics, SGS
  model updates, or turbine diagnostics run on periodic intervals.
- Report both regular-step timing and any periodic expensive-step timing.
- For benchmark tables, record the averaging window and whether startup steps
  were excluded.

## SGS and Wall Model

- SGS model dispatch must cover every supported runtime value, not only
  `sgs=5`.
- Unsupported runtime values should fail clearly rather than silently running a
  partially optimized path.
- IWM fields such as `iwm_flt_us`, `iwm_tR`, `iwm_Dz`, and wall stresses are
  wall-surface fields indexed by `(i,j)`.  Do not convert these into diagonal
  or per-row accumulators unless a separate algorithm proves that change.

## ATM and Structural Solver

- Keep induced-velocity correction enabled in source for the canonical branch.
- Keep structure-on and structure-off benchmarks paired when changing ATM or
  coupling code.
- ATM load-balancing changes should include both correctness diagnostics and a
  load-imbalance metric.  A load-balancing path should not be default unless
  projected benefit exceeds communication overhead.

## Documentation Required For New Optimizations

Every accepted optimization should leave:

- A concise source comment if the implementation is non-obvious.
- A benchmark record containing build options and run configuration.
- A before/after timing comparison.
- A numerical sanity check against the immediate previous accepted version.

## Pre-Commit Readability Checks

Run the consolidated local check before committing source, script, or
documentation edits:

```bash
python3 tools/check_branch_readiness.py
```

Run the CMake smoke variant when the local machine has CMake and a suitable
compiler wrapper:

```bash
python3 tools/check_branch_readiness.py --with-cmake-configure
```

For optional inflow or GPU source-membership changes, also run:

```bash
python3 tools/check_branch_readiness.py --with-hit-cmake-configure
```

Key checks include:

```bash
python3 tools/check_source_hygiene.py
python3 tools/check_cmake_option_docs.py
python3 tools/check_build_profile_args.py
python3 tools/check_environment_switch_docs.py
```

The complete readiness-check list is maintained in
`README_FINAL_OPTIMIZED_20260619.md`; `tools/check_readiness_docs.py` keeps that
list synchronized with `tools/check_branch_readiness.py`.

The hygiene check intentionally enforces portable ASCII comments and docs for
tracked source files.  Use plain-text forms such as `->`, `x`, `dx`, and
`omega` instead of Unicode arrows, multiplication signs, partial derivatives,
or Greek letters.  That keeps diffs readable in HPC shell logs,
Windows terminals, and GitHub reviews.

The environment-switch check scans tracked Fortran source for `LESGO_*`
literals and fails if any solver switch is missing from
`docs/environment_switches.md`.  Benchmark runner variables are intentionally
outside that check.

The CMake-option check scans the root `CMakeLists.txt` for `USE_*` build
options and fails if any option is missing from `docs/code_organization.md`.
The build-profile argument check scans `docs/build_profiles.md` for `-D...`
arguments and fails if a copy-paste profile references a non-public CMake knob.

Run `git diff --check` separately when using a single native Git client.  It is
not part of the default wrapper because mixed Windows/WSL checkouts can report
CRLF-only differences as whitespace errors.  The default source-hygiene check
still enforces LF line endings for tracked human-edited files.

## Line Endings

`.gitattributes` defines source, docs, scripts, CMake files, and case inputs as
LF text.  Keep generated numerical and visual outputs out of the source tree, or
mark them binary before tracking.  Do not hand-convert line endings to satisfy a
single platform; let Git apply the repository attributes.
