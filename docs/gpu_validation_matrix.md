# GPU Validation Matrix

This matrix is the release-facing checklist for the non-LVLSET GPU port.  It
separates three different claims:

- source coverage: the GPU path exists in the source;
- correctness: the GPU path runs and produces acceptable diagnostics;
- performance: the GPU path has a paired CPU/GPU comparison and is faster for a
  representative case.

Do not replace a missing paired CPU/GPU benchmark with a single GPU runtime.
Single GPU runtimes prove that a path executes; they do not prove speedup.
The machine-readable benchmark plan for the rows below lives in
`docs/gpu_benchmark_manifest.json`.
The current timing/evidence ledger lives in
`docs/gpu_validation_evidence.json`; it intentionally records GPU-runtime-only
evidence separately from paired CPU/GPU speedup evidence.
The map from parsed `lesgo.conf` groups to these validation rows lives in
`docs/lesgo_conf_gpu_validation_map.json`.

## Status Vocabulary

| Status | Meaning |
| --- | --- |
| `validated-gpu-runtime` | GPU run evidence exists, but a current paired CPU/GPU speedup record is not in this repository. |
| `recorded-correct-and-faster` | Prior CPU/GPU correctness and speedup evidence exists, but should be copied into the public validation record before a release claim. |
| `recorded-paired-not-faster` | Paired CPU/GPU evidence exists and the GPU cumulative average is not lower than the CPU baseline. |
| `implemented-needs-benchmark` | Source path exists, but it still needs targeted correctness and CPU/GPU performance validation. |
| `host-boundary` | Host-side work is expected by design, usually setup, restart, diagnostics, or output I/O. |
| `excluded` | Outside the current non-LVLSET optimization scope. |

## Validation Matrix

| ID | Path | Key switches or runtime settings | Current status | Current evidence | Required next evidence |
| --- | --- | --- | --- | --- | --- |
| `les_core_channel` | LES core channel-flow path | `USE_LES_GPU=ON`, no turbine module | `recorded-correct-and-faster` | Derecho 240^3 paired run completed: CPU 4.7443 s/step, GPU 0.05961 s/step, speedup 79.59x by cumulative average. | Keep this paired run as the public small LES-core baseline; rerun only if source or runtime settings change. |
| `adm_disk` | Actuator disk model | `USE_TURBINES=ON`, `USE_ATM=OFF` | `recorded-correct-and-faster` | Derecho 240^3 paired ADM run completed: CPU 4.1239 s/step, GPU 0.06114 s/step, speedup 67.45x by cumulative average. | Keep turbine forcing output and disk-velocity diagnostics attached to future reruns. |
| `atm_line` | Actuator line / 5 MW turbine path | `USE_ATM=ON`, `USE_TURBINES=OFF` | `recorded-correct-and-faster` | Derecho 240^3 paired ATM runs completed: structure off 55.13x, structure on 66.40x by cumulative average. | Structure-on was run with `LESGO_ATM_STRUCTURE=1`; keep structure on/off paired when ATM changes. |
| `large_windfarm` | 60-turbine wind-farm path | `3072 x 384 x 400`, 60 turbines, multi-GPU | `recorded-correct-and-faster` | Derecho current-source 50-step pair completed: CPU50 cumulative average 23.34388 s/step, GPU16 cumulative average 0.27913 s/step, speedup 83.63x. | Keep CPU/GPU logs and 60 turbine power-file counts attached to future large-case reruns. |
| `sgs_disabled` | SGS-off branch | `sgs = .false.` | `recorded-correct-and-faster` | Derecho compact 128^3 paired run completed with `sgs=.false.`: CPU 0.40972 s/step, GPU 0.01167 s/step, speedup 35.11x by cumulative average. | Keep this compact no-SGS branch in future SGS smoke gates; production cases can use the same importer workflow if needed. |
| `sgs_models_1_5` | Supported SGS runtime values | `sgs_model = 1..5` | `recorded-correct-and-faster` | Derecho compact 128^3 paired matrix completed for `sgs_model=1..5`; speedups were 33.76x, 16.22x, 11.27x, 40.63x, and 48.36x. | Re-run this matrix after SGS edits; logs used `nenergy=4`, `cs_count=5`, and `dyn_init=1` to cover regular and SGS-update cadence. |
| `dyn_tn` | Dynamic Lagrangian timescale | `USE_DYN_TN=ON` | `recorded-correct-and-faster` | Compact 128^3 paired runs passed for `sgs_model = 4` and `5`; GPU cumulative-average speedups were 49.00x and 56.42x. | Keep this compact check in the release evidence set when SGS/Dyn-TN kernels change. |
| `iwm_wall_model` | Integral wall-model-heavy cases | `lbc_mom`, `ubc_mom`, IWM settings | `recorded-correct-and-faster` | Current-source IWM case passes paired CPU/GPU field, wall-stress diagnostic, divergence, and timing checks. | Keep IWM wall-stress diagnostics in future regression runs. |
| `scalar_passive` | Passive scalar transport | `USE_SCALARS=ON`, `USE_SCALARS_GPU=ON`, `passive_scalar = .true.` | `recorded-correct-and-faster` | GPU scalar source and scalar readiness contracts exist. | 128^3 and 240^3 passive scalar CPU/GPU evidence recorded. |
| `scalar_active` | Active scalar / stability coupling | `USE_SCALARS=ON`, `USE_SCALARS_GPU=ON`, `passive_scalar = .false.` | `recorded-correct-and-faster` | GPU scalar source includes stability helper routines. | Active-buoyancy CPU/GPU scalar diagnostics recorded. |
| `cps_velocity` | Concurrent precursor velocity exchange | `USE_CPS=ON`, scalar off | `recorded-correct-and-faster` | CPS velocity exchange has device-resident `host_data use_device` paths. | Current-source two-color CPS velocity CPU/GPU timing and checkpoint evidence recorded. |
| `cps_scalar` | Concurrent precursor with scalar coupling | `USE_CPS=ON`, `USE_SCALARS_GPU=ON` | `recorded-correct-and-faster` | CPS scalar buffers are device-resident only when `PPSCALARS_GPU` is enabled. | Current-source two-color CPS scalar CPU/GPU timing and scalar checkpoint evidence recorded. |
| `hit_inflow` | HIT inflow path | `USE_HIT=ON`, `inflow_type = 2` | `recorded-correct-and-faster` | Prior optional HIT validation recorded CPU/GPU correctness and speedup. | Copy the HIT validation report into the public validation record or rerun it. |
| `shifted_inflow` | Shifted inflow path | shifted-inflow runtime settings | `recorded-correct-and-faster` | Current-source 128^3 CPU/GPU shifted-inflow pair completed with masked checkpoint agreement. | Keep as CPU-equivalence validation; shifted inflow is numerically sensitive to forcing/timestep choices. |
| `sponge_coriolis` | Sponge and Coriolis forcing | `use_sponge`, `coriolis_forcing` | `recorded-correct-and-faster` | Current-source 128^3 sponge-on and Coriolis-on CPU/GPU variants completed with matching diagnostics. | Coriolis checkpoint comparison has small physical-field drift; use diagnostics and tolerance report for release evidence. |
| `adm_dynamic_controls` | Dynamic ADM controls | dynamic yaw/Ct/rotation/correction settings | `recorded-correct-and-faster` | Isolated current-source ADM dynamic-control variants pass CPU/GPU correctness and timing evidence. | Keep dynamic yaw/Ct/rotation/correction variants in future regression runs. |
| `diagnostics_output` | Diagnostics, restart, time averages, plane/domain output | output and checkpoint controls | `host-boundary` | Host-visible I/O and statistics boundaries are expected. | Validate correctness and overhead separately; do not count output steps as regular-step speed evidence. |
| `cgns_output` | CGNS output | `USE_CGNS=ON` | `host-boundary` | Optional host I/O dependency path. | Treat as compatibility validation, not GPU hot-path performance. |
| `lvlset` | Level-set module | `USE_LVLSET=ON`, `USE_LVLSET_GPU=ON` | `excluded` | Excluded from this legacy non-Level-Set release matrix; the dedicated CPU/GPU matrix passed with worst normalized field error `3.302591e-13`. | Maintain the separate acceptance record in `docs/level_set_gpu_port.md`; add a one-rank-per-GPU CUDA-aware MPI run before publishing multi-GPU scaling. |
