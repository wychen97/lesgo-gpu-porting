# GPU Validation Runbook

This runbook turns the non-LVLSET validation matrix into benchmark batches.
It is intentionally conservative: a case is not "faster on GPU" until paired
CPU/GPU evidence from the same source tree, grid, compiler family, runtime
settings, and output cadence is recorded in `docs/gpu_validation_evidence.json`.

Source-of-truth files:

- `docs/gpu_validation_matrix.md`: release-facing status table;
- `docs/gpu_benchmark_manifest.json`: required settings and evidence fields;
- `docs/gpu_validation_evidence.json`: current evidence and speedup claims;
- `docs/gpu_port_coverage_audit.md`: source coverage and remaining risk notes.

## Minimum Evidence Record

Every benchmark record must include:

- source identity: branch, commit, and whether the tree was dirty;
- build identity: full CMake cache arguments and compiler/MPI stack;
- runtime identity: `lesgo.conf`, turbine/scalar/HIT input files, grid, `nproc`,
  CPU/GPU count, `nsteps`, and output cadence;
- correctness evidence: divergence, kinetic energy, and module-specific output
  sanity checks;
- timing evidence: cumulative wall time, cumulative average `s/step`, late-step
  timing, and explicit notes for diagnostic or `cs_count` update steps;
- raw log paths or copied log snippets sufficient to reproduce the numbers.

Use `tools/import_lesgo_timing_evidence.py` for standard LESGO stdout/stderr
logs.  It parses the timing block and updates the evidence ledger in one step.
Use `tools/parse_lesgo_timing.py` directly only when inspecting a log without
changing the ledger.
Use `tools/gpu_validation_plan.py` before launching a batch to expand the
manifest and evidence ledger into concrete CPU/GPU run tasks.

```bash
python3 tools/gpu_validation_plan.py
python3 tools/gpu_validation_plan.py --priority p0-public-small --commands derecho-submit
python3 tools/gpu_validation_plan.py --priority p0-public-small --commands import
python3 tools/gpu_validation_plan.py --priority p0-public-small --commands pair-import
python3 tools/report_gpu_matrix_status_updates.py
```

Use priority filters to stage the remaining evidence work:

```bash
python3 tools/gpu_validation_plan.py --priority p0-public-small
python3 tools/gpu_validation_plan.py --priority p1-core-options
python3 tools/gpu_validation_plan.py --priority p2-optional-coupling
python3 tools/gpu_validation_plan.py --priority p3-large
python3 tools/gpu_validation_plan.py --priority p4-compatibility
```

The recommended order is:

1. `p0-public-small`: paired CPU/GPU baselines for the small public cases.
2. `p1-core-options`: SGS, dynamic-timescale, 128^3 scalar, and HIT checks.
3. `p2-optional-coupling`: CPS, IWM, shifted inflow, forcing, ADM controls,
   and 240^3 scalar checks.
4. `p3-large`: the expensive 60-turbine reference case.
5. `p4-compatibility`: diagnostics and CGNS output compatibility.

For the public p0 cases, `--commands derecho-submit` emits Derecho compile and
`qsub` templates for the case-local scripts.  Submit jobs from a reused case
directory sequentially, or copy the case tree first, because the run scripts
clean output directories before execution.

Evidence statuses used after paired runs:

- `paired_speedup_claimed`: CPU and GPU cumulative averages exist, the GPU
  average is lower, `speedup_claim` records the ratio, and passing
  correctness evidence is attached.
- `paired_not_faster`: CPU and GPU cumulative averages exist, but no GPU speedup
  is claimed.

Do not use a final printed step alone as an average.  Output-heavy, restart,
statistics, and SGS-update timesteps must be labeled separately from regular
timestep timing.

## Batch 1: Public Four-Case Paired Baselines

Purpose: convert the public presentation examples from GPU-runtime evidence to
paired CPU/GPU evidence.

Rows:

- `les_core_channel`
- `adm_disk`
- `atm_line`
- `large_windfarm`

Required action:

- rerun or copy immutable CPU/GPU logs for the same current source tree;
- record cumulative average and late-step timing separately;
- update `docs/gpu_validation_evidence.json` only after paired CPU/GPU evidence
  exists.

For the small public cases, the case-local Derecho compile scripts support both
profiles:

```bash
cd test-cases/channel_flow
./compile_derecho.sh cpu
./compile_derecho.sh gpu

cd ../adm_disk
./compile_derecho.sh cpu
./compile_derecho.sh gpu

cd ../atm_line
./compile_derecho.sh cpu
./compile_derecho.sh gpu
```

Each script writes `lesgo-run-exe-cpu` and `lesgo-run-exe-gpu`. The PBS
resource request still needs to match the selected executable: CPU baselines
should use CPU resources, while GPU baselines should use A100 GPU resources.
The case-local `submit_derecho.pbs` files accept `RUN_PROFILE=cpu` or
`RUN_PROFILE=gpu`; the README in each case gives the exact `qsub` commands.
For `RUN_PROFILE=cpu`, the submit script uses CPU runtime modules and unsets
the Cray GPU-aware MPI environment variables. For `RUN_PROFILE=gpu`, it loads
CUDA and enables GPU-aware MPI.

After the p0 jobs finish, the submit scripts preserve logs under
`run-archives/<RUN_LABEL>/lesgo_<RUN_LABEL>.log`.  Import the archived paired
logs with:

```bash
python3 tools/import_p0_archived_evidence.py \
  --source "Derecho p0 paired run, branch/commit/source-tree details here" \
  --case les_core_channel

python3 tools/import_p0_archived_evidence.py \
  --source "Derecho p0 paired run, branch/commit/source-tree details here" \
  --case adm_disk \
  --module-check "adm_disk=turbine forcing and disk velocity matched CPU reference"

python3 tools/import_p0_archived_evidence.py \
  --source "Derecho p0 paired run, branch/commit/source-tree details here" \
  --case atm_line.structure_off \
  --module-check "atm_line.structure_off=turbine power and structure-off outputs matched CPU reference"

python3 tools/import_p0_archived_evidence.py \
  --source "Derecho p0 paired run, branch/commit/source-tree details here" \
  --case atm_line.structure_on \
  --module-check "atm_line.structure_on=turbine power and structure-on outputs matched CPU reference"
```

The archive importer still calls `tools/import_lesgo_timing_pair.py`, so the
same diagnostic comparison, evidence-item coverage, and no-unproven-speedup
rules apply.

## Batch 2: SGS And Wall-Model Matrix

Purpose: close the highest-risk runtime switches in `lesgo.conf` model and wall
settings.

Rows:

- `sgs_disabled`
- `sgs_models_1_5`
- `dyn_tn`
- `iwm_wall_model`

Required action:

- run compact CPU/GPU cases for `sgs=false`;
- run `sgs_model=1..5`, separating regular timesteps from `cs_count` update
  timesteps;
- run `USE_DYN_TN=ON` for `sgs_model=4` and `sgs_model=5`;
- run an IWM-heavy wall-model case before making broad wall-model speed claims.

## Batch 3: Scalar And CPS Coupling

Purpose: validate optional modules that can materially change data residency and
communication cost.

Rows:

- `scalar_passive`
- `scalar_active`
- `cps_velocity`
- `cps_scalar`

Required action:

- run passive scalar CPU/GPU at 128^3 and 240^3;
- run active scalar CPU/GPU with buoyancy/stability diagnostics;
- run CPS velocity exchange without scalar coupling;
- run CPS with scalar coupling and record both velocity and scalar stage timing.

## Batch 4: Optional Inflow, Forcing, And ADM Controls

Purpose: cover less common but user-visible `lesgo.conf` paths before claiming
full runtime-option coverage.

Rows:

- `hit_inflow`
- `shifted_inflow`
- `sponge_coriolis`
- `adm_dynamic_controls`

Required action:

- copy the prior HIT CPU/GPU validation record into this repository or rerun it;
- run shifted-inflow CPU/GPU with identical input fields;
- run sponge-on and Coriolis-on variants;
- run dynamic ADM yaw/Ct/rotation/correction controls one at a time.

## Batch 5: Host-Boundary Compatibility

Purpose: validate output compatibility and overhead without treating host I/O as
a GPU hot-path speedup target.

Rows:

- `diagnostics_output`
- `cgns_output`

Required action:

- verify checkpoint/restart, time averages, and plane/domain outputs;
- verify CGNS output only when `USE_CGNS=ON` is needed;
- report overhead separately from ordinary timestep performance.

## Excluded Scope

Rows:

- `lvlset`

LVLSET is excluded from the current non-LVLSET GPU validation scope.  It needs a
separate project if users decide to optimize it later.

## Completion Rule

After each batch:

1. Copy or summarize the authoritative logs into the validation record.
2. Import matched CPU/GPU logs with `tools/import_lesgo_timing_pair.py`.
   Use `--compare-diagnostics` to compare parsed MPI exit status, divergence,
   and kinetic energy.  A speedup claim must also include any needed
   module-specific `--correctness-check`, such as turbine-output agreement,
   scalar-field comparison, or a field-difference report.  Repeat
   `--evidence-item` until every token in that row's
   `docs/gpu_benchmark_manifest.json` `required_evidence` list is covered.
3. For nonstandard records, update `docs/gpu_validation_evidence.json` with
   `tools/update_gpu_validation_evidence.py`.
4. Run `tools/report_gpu_matrix_status_updates.py` and update
   `docs/gpu_validation_matrix.md` to the matching stronger status:
   `recorded-correct-and-faster` or `recorded-paired-not-faster`.
5. Run `python3 tools/check_branch_readiness.py`.

The stricter release-objective gate is separate from normal readiness:

```bash
python3 tools/require_gpu_release_objective.py
```

That command should fail until all non-LVLSET validation rows have complete
paired CPU/GPU correctness evidence and acceptable GPU-faster results.  Use it
before release wording says the GPU port is complete.
