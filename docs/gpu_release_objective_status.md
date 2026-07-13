# GPU Release Objective Status

Scope: non-LVLSET LESGO GPU port.  This report connects the static
Fortran function inventory, `lesgo.conf` validation surface, and paired
CPU/GPU evidence ledger.  It is generated and should not be edited by
hand.

Release objective met: `yes`
Strict release-gate gaps: `0`

Run the strict gate directly with:

```bash
python3 tools/require_gpu_release_objective.py
```

## Static Source Inventory

| Classification | Subprograms | Meaning |
| --- | ---: | --- |
| `gpu-marked` | 147 | Contains OpenACC, CUDA, GPU-aware MPI, or GPU preprocessor markers. |
| `gpu-file-unmarked` | 2 | Helper inside a GPU source file but without a local marker. |
| `host-boundary` | 102 | Setup, I/O, configuration, MPI definitions, or other expected host boundary. |
| `host-or-diagnostic` | 52 | Diagnostics, initialization, reporting, restart, or low-frequency host work. |
| `unmarked-runtime-candidate` | 107 | Needs runtime validation or profiling before broad GPU speed claims. |

## Validation Evidence Summary

| Evidence state | Rows |
| --- | ---: |
| `release-proven` | 18 |

`release-proven` rows have paired CPU/GPU correctness evidence and a
GPU-faster result.  Rows with only a single GPU runtime remain execution
evidence only.

## lesgo.conf Parser-Key Coverage

| Parser group | Keys | Release-proven keys | Validation row states |
| --- | ---: | ---: | --- |
| `DOMAIN` | 9 | 9 | `les_core_channel`=release-proven, `large_windfarm`=release-proven |
| `MODEL` | 11 | 11 | `sgs_disabled`=release-proven, `sgs_models_1_5`=release-proven, `dyn_tn`=release-proven, `iwm_wall_model`=release-proven |
| `CORIOLIS` | 11 | 11 | `sponge_coriolis`=release-proven |
| `FLOW_COND` | 23 | 23 | `les_core_channel`=release-proven, `iwm_wall_model`=release-proven, `shifted_inflow`=release-proven, `sponge_coriolis`=release-proven |
| `OUTPUT` | 33 | 33 | `diagnostics_output`=release-proven |
| `TURBINES` | 22 | 22 | `adm_disk`=release-proven, `atm_line`=release-proven, `large_windfarm`=release-proven, `adm_dynamic_controls`=release-proven |
| `SCALARS` | 13 | 13 | `scalar_passive`=release-proven, `scalar_active`=release-proven, `cps_scalar`=release-proven |
| `HIT` | 11 | 11 | `hit_inflow`=release-proven |

Total parsed non-LVLSET `lesgo.conf` keys: `133`
Release-proven parser keys: `133`

## Candidate Buckets Linked To Evidence

| Review bucket | Candidates | Validation row states | Meaning |
| --- | ---: | --- | --- |
| `adm-cpu-fallback-profile` | 5 | `adm_disk`=release-proven, `adm_dynamic_controls`=release-proven | ADM/turbine CPU fallback or compatibility routines; profile before treating them as missing GPU work. |
| `atm-host-model` | 35 | `atm_line`=release-proven, `large_windfarm`=release-proven | ATM blade, controller, structure, and small math helpers that remain host-side in the current hybrid design. |
| `atm-mirror-lb-control` | 11 | `atm_line`=release-proven, `large_windfarm`=release-proven | ATM mirror, synchronization, cell-search, and load-balance control helpers around the GPU sampling/forcing path. |
| `cpu-fallback-compat` | 27 | `les_core_channel`=release-proven, `hit_inflow`=release-proven | CPU fallback or host compatibility routines retained beside GPU production paths. |
| `diagnostic-profiling` | 8 | `diagnostics_output`=release-proven | Profiling, timing, or audit helpers; not GPU hot-path kernels. |
| `excluded-lvlset-bridge` | 1 | `lvlset`=excluded | LVLSET bridge code excluded from the current non-LVLSET scope. |
| `generic-helper-profile` | 10 | `les_core_channel`=release-proven, `adm_disk`=release-proven, `atm_line`=release-proven | Generic interpolation/math helpers that may be CPU fallback or low-cost support code. |
| `inflow-fringe-profile` | 2 | `hit_inflow`=release-proven, `shifted_inflow`=release-proven, `sponge_coriolis`=release-proven | Inflow/fringe helpers that need targeted runtime validation for nonstandard inflow configurations. |
| `iwm-wallmodel-profile` | 1 | `iwm_wall_model`=release-proven | IWM wall-model candidate that needs an IWM-heavy correctness and timing case before broad speed claims. |
| `scalar-init-fallback` | 7 | `scalar_passive`=release-proven, `scalar_active`=release-proven, `cps_scalar`=release-proven | Scalar initialization, stability helper, or CPU fallback routines; validate passive and active scalar cases separately. |

## Open Validation Rows

| Validation row | Current evidence state | Required next evidence |
| --- | --- | --- |

## Interpretation

- The source tree has broad GPU coverage, but the release objective is not
  proven until the open validation rows have paired CPU/GPU correctness
  and timing evidence.
- Host-boundary rows are not GPU hot-path speed claims; they need
  compatibility and overhead evidence.
- LVLSET remains excluded from this audit.

Regenerate this file with:

```bash
python3 tools/report_gpu_release_objective_status.py --write
```
