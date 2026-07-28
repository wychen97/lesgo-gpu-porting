# GPU Port Coverage Audit

This audit records the current non-LVLSET GPU-port coverage by build option
and runtime path.  It is intentionally stricter than the handoff summary:
source code can be GPU-ported while still needing stronger production-size
validation before it should be advertised as faster for every `lesgo.conf`
configuration.

Scope:

- included: non-LVLSET LESGO paths and optional modules;
- excluded: `USE_LVLSET=ON`, `level_set*.f90`, and `trees_*_ls.f90`;
- source baseline: current public `lesgo-gpu-porting` worktree;
- validation baseline: published Derecho validation records and local source
  checks available at the time of this audit.

## Timing Interpretation

Printed per-step timing and cumulative average timing are different quantities.
For the small published 240^3 GPU cases, the cumulative average is the
performance number to use for total runtime estimates.  The final printed step
is a late-step sample and is slower than the whole-run average because early
startup steps are cheaper and diagnostic cadence can affect individual printed
steps.

| Case | Cumulative average | Cumulative wall | Steps | Final printed timing |
| --- | ---: | ---: | ---: | ---: |
| channel flow | 0.05961 s/step | 11.92148 s | 200 | 0.1150396 s/step |
| ADM disk | 0.06114 s/step | 12.22882 s | 200 | 0.1164596 s/step |
| ATM line | 0.07957 s/step | 15.91372 s | 200 | 0.1194032 s/step |
| three-case total | 0.06677 s/step | 40.06402 s | 600 | N/A |

Use final printed timing for late-step behavior, and cumulative average for
total wall-time estimates.  Do not mix the labels.

## Coverage Summary

| Area | Build/runtime switches | GPU-port status | Validation status | Remaining action |
| --- | --- | --- | --- | --- |
| LES core | `USE_LES_GPU=ON` | GPU path implemented through explicit OpenACC residency plus cuFFT helpers. | Published channel-flow and wind-farm cases run successfully. | Continue comparing late-window timings against the immediate previous accepted branch. |
| Convection | `PPCONVEC_GPU` | GPU implementation in `convec_gpu.f90`; CPU path remains in `convec.f90`. | Covered by published channel, ADM, ATM, and large wind-farm runs. | None for ordinary production cases. |
| Derivatives/filtering | `PPDERIVS_GPU` | GPU implementation in `derivatives_gpu.f90` with batched spectral helpers. | Covered by published production runs. | Recheck if new `lesgo.conf` paths alter derivative cadence or boundary behavior. |
| Pressure/tridiagonal | `PPPRESS_GPU`, `PPGPU_AWARE_MPI` | GPU pressure and GPU-aware MPI paths implemented; host-staged fallback remains correct but is slower. | Covered by published production runs on Derecho GPU-aware MPI path. | Report GPU-aware and host-staged fallback timings separately. |
| SGS disabled | `sgs = false` | Device stress/div-stress path skips eddy viscosity and uses molecular branch. | Documented as supported by SGS matrix. | Keep in every future SGS smoke gate. |
| SGS models 1..5 | `sgs_model = 1..5` | Runtime dispatch in `sgs_gpu.f90` covers Smagorinsky, standard dynamic, scale-dependent dynamic, Lagrangian similarity, and Lagrangian scale-dependent. | Documented as supported by SGS matrix; main published cases exercise `sgs = 5`. | Re-run compact SGS matrix after future SGS edits; include regular-step and SGS-update-step timings. |
| Dynamic timescale update | `USE_DYN_TN=ON` | Source has `PPDYN_TN` branches in Lagrangian SGS code. | Not part of the canonical production validation. | Needs dedicated `USE_DYN_TN=ON` CPU/GPU correctness and timing run. |
| Wall stress and IWM | `lbc_mom`, `ubc_mom`, IWM settings | `iwmles.f90` and SGS wall branches contain GPU paths and protected `(i,j)` surface indexing. | Source indexing guard passes; small correctness evidence exists. | Needs production-size wall-model benchmark before claiming speedup for IWM-heavy cases. |
| ADM turbines | `USE_TURBINES=ON`, `USE_ATM=OFF` | `turbines_gpu.f90` and forcing integration support the GPU LES path. | Published ADM disk 240^3 case passed. | For presentations, use cumulative average and late-step timing with correct labels. |
| ATM / structural coupling | `USE_ATM=ON` | The atPoint sampler and force deposition are device-resident; induced-velocity correction and structure on/off are integrated with the host turbine model. Spalart/nacelle use an explicit compatibility bridge. | Standard structure-off/on A/B and optional atPoint+nacelle, Spalart, and held-force CPU/GPU checks passed in the 2026-07-27 architecture review. | Profile optional host bridges only on representative cases; aligned checkpoints are required when `updateInterval > 1`. |
| Scalar transport | `USE_SCALARS=ON`, `USE_SCALARS_GPU=ON` | Current `scalars.f90` includes GPU transport, scalar derivatives, halo/update, and stability helper routines. | All-module and scalar smoke evidence exists; not part of the four published presentation cases. | Re-run current public branch scalar CPU/GPU at 128^3 and 240^3, passive and active buoyancy, before broad speedup claims. |
| Scalar CPU fallback | `USE_SCALARS=ON`, `USE_SCALARS_GPU=OFF` | Intentionally CPU fallback. | Correctness path only, not a GPU speed path. | Do not compare as GPU scalar performance. |
| CPS / concurrent precursor | `USE_CPS=ON` | CPS velocity exchange and fringe application have device-resident paths; scalar CPS is GPU-resident only with `PPSCALARS_GPU`. | CPS/all-module smoke evidence exists. | Needs production-size CPS timing with and without scalar coupling. |
| HIT inflow | `USE_HIT=ON`, `inflow_type = 2` | HIT file/restart I/O remains host-owned; per-step plane interpolation and fringe blending are GPU-assisted in `hit_inflow_gpu.f90`. | Optional HIT validation showed CPU/GPU correctness and GPU speedup on the tested cases. | Keep as optional until it is included in standard release validation. |
| Shifted inflow | shifted-inflow runtime path | Source contains GPU-aware branches for active LES GPU builds. | Not covered by the published four-case validation. | Add a shifted-inflow smoke if users rely on it. |
| Sponge and Coriolis forcing | `use_sponge`, `coriolis_forcing` | Source contains GPU branches and reductions. | Not separately isolated in the published four-case validation. | Add targeted runtime checks for configurations that use these options. |
| CGNS output | `USE_CGNS=ON` | Host I/O dependency path; not a GPU hot-path optimization target. | Not part of canonical production validation. | Treat as output compatibility, not GPU performance evidence. |
| Diagnostics, restart, statistics | output cadence and restart settings | Host-visible boundaries are expected; core fields are copied explicitly before output where needed. | Covered only as required by published cases. | Do not use output/statistics timesteps alone as regular-step performance evidence. |
| LVLSET | `USE_LVLSET=ON` | Excluded from this branch's optimized production scope. | Not validated. | Leave deferred unless a separate LVLSET project starts. |

## `lesgo.conf` Runtime Coverage Map

The parser in `input_util.f90` exposes more runtime combinations than the four
published cases exercise.  The most important non-LVLSET runtime switches are:

For a generated key-by-key view of the same mapping, see
`docs/lesgo_conf_key_validation_coverage.md`.

| `lesgo.conf` area | Main keys | Current GPU coverage note |
| --- | --- | --- |
| Domain | `nproc`, `Nx`, `Ny`, `Nz`, `uniform_spacing` | Normal decompositions are covered by the GPU source path. New decompositions should still check divisibility, memory footprint, and MPI timing. |
| Model/SGS | `sgs`, `sgs_model`, `cs_count`, `DYN_init`, `ifilter`, `wall_damp_exp` | Supported runtime SGS values are `1..5`; `sgs=false` is supported. Timing must separate normal steps from `cs_count` update steps. |
| Coriolis | `coriolis_forcing`, `G`, `alpha`, PID keys | GPU branches exist, but the published presentation cases use `coriolis_forcing=.false.`. Coriolis-on needs a targeted check. |
| Flow condition | `inflow_type`, `lbc_mom`, `ubc_mom`, `use_sponge`, mean/random pressure force keys | Inflow type `0/1` and common wall settings are covered by published cases. HIT, shifted inflow, sponge, random force, and less common wall settings need targeted checks before speed claims. |
| Output | `checkpoint_data`, `tavg_calc`, point/domain/plane output controls | These are expected host-visible I/O/statistics boundaries. They need correctness and overhead checks, but are not GPU hot-path success criteria. |
| Turbines | ADM layout, yaw/tilt, dynamic Ct/yaw, rotation, correction, output cadence | A simple ADM case is validated. Dynamic ADM control options need their own run before claiming broad coverage. |
| Scalars | scalar boundary, lapse rate, initial profile, active/passive scalar, `Pr_sgs` | GPU source exists only when `USE_SCALARS_GPU=ON`. Production-size scalar validation remains a high-priority gap. |
| HIT | `inflow_type=2`, HIT dimensions/files | Optional HIT has a GPU helper path and prior validation, but is not part of the standard production four-case gate. |

## Parser Key Inventory

This table is generated from the non-LVLSET parser blocks in `input_util.f90`.
It is not a claim that every key has already been benchmarked; it is the
checklist used to keep the coverage audit aligned with the actual runtime
configuration surface.

| Parser block | Keys parsed from `input_util.f90` |
| --- | --- |
| `DOMAIN` | `LX`, `LY`, `LZ`, `NPROC`, `NX`, `NY`, `NZ`, `UNIFORM_SPACING`, `Z_I` |
| `MODEL` | `CO`, `CS_COUNT`, `DYN_INIT`, `IFILTER`, `MOLEC`, `NU_MOLEC`, `SGS`, `SGS_MODEL`, `U_STAR`, `VONK`, `WALL_DAMP_EXP` |
| `CORIOLIS` | `ALPHA`, `CORIOLIS_FORCING`, `FC`, `G`, `HEIGHT_SET`, `KD`, `KI`, `KP`, `PHI_SET`, `PID_TIME`, `REPEAT_INTERVAL` |
| `FLOW_COND` | `EVAL_MEAN_P_FORCE`, `FRINGE_REGION_END`, `FRINGE_REGION_LEN`, `INFLOW_TYPE`, `INFLOW_VELOCITY`, `INILAG`, `INITU`, `LBC_MOM`, `MEAN_P_FORCE_X`, `MEAN_P_FORCE_Y`, `RMS_RANDOM_FORCE`, `SAMPLING_REGION_END`, `SHIFT_N`, `SPONGE_FREQUENCY`, `SPONGE_HEIGHT`, `STOP_RANDOM_FORCE`, `UBC_MOM`, `UBOT`, `USE_MEAN_P_FORCE`, `USE_RANDOM_FORCE`, `USE_SPONGE`, `UTOP`, `ZO` |
| `OUTPUT` | `CHECKPOINT_DATA`, `CHECKPOINT_NSKIP`, `DOMAIN_CALC`, `DOMAIN_NEND`, `DOMAIN_NSKIP`, `DOMAIN_NSTART`, `LAG_CFL_COUNT`, `NENERGY`, `POINT_CALC`, `POINT_LOC`, `POINT_NEND`, `POINT_NSKIP`, `POINT_NSTART`, `TAVG_CALC`, `TAVG_NEND`, `TAVG_NSKIP`, `TAVG_NSTART`, `WBASE`, `XPLANE_CALC`, `XPLANE_LOC`, `XPLANE_NEND`, `XPLANE_NSKIP`, `XPLANE_NSTART`, `YPLANE_CALC`, `YPLANE_LOC`, `YPLANE_NEND`, `YPLANE_NSKIP`, `YPLANE_NSTART`, `ZPLANE_CALC`, `ZPLANE_LOC`, `ZPLANE_NEND`, `ZPLANE_NSKIP`, `ZPLANE_NSTART` |
| `TURBINES` | `ADM_CORRECTION`, `ALPHA1`, `ALPHA2`, `CT_PRIME`, `DIA_ALL`, `DYN_CT_PRIME`, `DYN_THETA1`, `DYN_THETA2`, `FILTER_CUTOFF`, `HEIGHT_ALL`, `NUM_X`, `NUM_Y`, `ORIENTATION`, `READ_PARAM`, `STAG_PERC`, `TBASE`, `THETA1_ALL`, `THETA2_ALL`, `THK_ALL`, `TIP_SPEED_RATIO`, `T_AVG_DIM`, `USE_ROTATION` |
| `SCALARS` | `FLUX_BOT`, `G`, `IC_NO_VEL_NOISE_Z`, `IC_THETA`, `IC_Z`, `LAPSE_RATE`, `LBC_SCAL`, `PASSIVE_SCALAR`, `PR_SGS`, `READ_LBC_SCAL`, `SCAL_BOT`, `T_SCALE`, `ZO_S` |
| `HIT` | `LX_HIT`, `LY_HIT`, `LZ_HIT`, `NX_HIT`, `NY_HIT`, `NZ_HIT`, `TI_OUT`, `UP_IN`, `U_FILE`, `V_FILE`, `W_FILE` |

## Static Function-Level Inventory

`tools/report_gpu_static_inventory.py` scans tracked root Fortran files,
excluding LVLSET, and classifies subroutines/functions by static GPU markers.
This is a heuristic.  It does not prove speedup, and it can mark small helper
functions as runtime candidates until a human classifies them.

For the generated full function-level classification table, see
`docs/gpu_static_full_inventory.md`.

Current result:

| Classification | Subprograms |
| --- | ---: |
| `gpu-marked` | 147 |
| `gpu-file-unmarked` | 2 |
| `host-boundary` | 102 |
| `host-or-diagnostic` | 41 |
| `unmarked-runtime-candidate` | 104 |

The 104 unmarked runtime candidates now have review buckets from
`tools/report_gpu_static_inventory.py --review`:

For the generated full function-level table, see
`docs/gpu_static_candidate_review.md`.

| Review bucket | Candidates | Meaning |
| --- | ---: | --- |
| `adm-cpu-fallback-profile` | 5 | ADM/turbine CPU fallback or compatibility routines; profile before treating them as missing GPU work. |
| `atm-host-model` | 35 | ATM blade, controller, structure, and small math helpers that remain host-side in the current hybrid design. |
| `atm-mirror-lb-control` | 8 | ATM synchronization, cell-search, and active coupling-control helpers around the GPU sampling/forcing path. Historical disabled mirror/load-balance routines were removed. |
| `cpu-fallback-compat` | 27 | CPU fallback or host compatibility routines retained beside GPU production paths. |
| `diagnostic-profiling` | 8 | Profiling, timing, or audit helpers; not GPU hot-path kernels. |
| `excluded-lvlset-bridge` | 1 | LVLSET bridge code excluded from the current non-LVLSET scope. |
| `generic-helper-profile` | 10 | Generic interpolation/math helpers that may be CPU fallback or low-cost support code. |
| `inflow-fringe-profile` | 2 | Inflow/fringe helpers that need targeted runtime validation for nonstandard inflow configurations. |
| `iwm-wallmodel-profile` | 1 | IWM wall-model candidate that needs an IWM-heavy correctness and timing case before broad speed claims. |
| `scalar-init-fallback` | 7 | Scalar initialization, stability helper, or CPU fallback routines; validate passive and active scalar cases separately. |

The important distinction is that this list is no longer interpreted as
107 missing GPU kernels.  It is a review queue: most entries are host/control
work, CPU fallback compatibility, diagnostics, or targeted benchmark gaps.

Each bucket is tied to one or more validation rows so the static review remains
connected to runtime evidence:

| Review bucket | Validation rows |
| --- | --- |
| `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `diagnostic-profiling` | `diagnostics_output` |
| `excluded-lvlset-bridge` | `lvlset` |
| `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `iwm-wallmodel-profile` | `iwm_wall_model` |
| `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |

Manual review of the highest-impact candidates:

| Candidate | Review result |
| --- | --- |
| `actuator_turbine_model.f90` blade/controller/structure routines | Host-side by current design. The PPLES path mirrors blade points/forces and applies LES sampling/forcing on the GPU, but turbine model and structural logic still execute on the host. Structure-on performance therefore needs paired CPU/GPU evidence. |
| `atm_lesgo_interface.f90` coupling helpers | Active helpers now cover batched atPoint residency, packed gather, and explicit optional Spalart/nacelle bridging. Historical disabled shadow/mirror and point-owner branches were removed. |
| `scalars.f90::to_big` | CPU FFTW fallback. The GPU scalar path uses separate `PPSCALARS_GPU` routines; scalar speedup still needs production-size validation. |
| `scalars.f90::ic_scal_*` and `stability.f90` | Scalar initialization and CPU/fallback stability helpers. They are not proof of missing GPU transport work, but active-scalar validation remains required. |
| `hit_inflow.f90::interpolate3D` | CPU HIT fallback. The per-step GPU-assisted HIT work lives in `hit_inflow_gpu.f90`. |
| `iwmles.f90::iwm_slv` | CPU/fallback nonlinear wall solve. In the `PPSGS_GPU` path, `iwm_calc_wallstress` returns before the host Newton loop that calls `iwm_slv`; IWM-heavy cases still need correctness and timing validation. |
| `tridag_array.f90::tridag_array` | CPU/fallback tridiagonal solver. The production GPU path uses `tridag_gpu.f90` and pressure GPU orchestration. |

Regenerate the static view with:

```bash
python3 tools/report_gpu_static_inventory.py
python3 tools/report_gpu_static_inventory.py --candidates
python3 tools/report_gpu_static_inventory.py --review
```

## Current Conclusion

The major non-LVLSET hot path is GPU-ported: LES core, derivatives, convection,
pressure, supported SGS runtime values, ADM, ATM, and the published wind-farm
path all have GPU implementations and validation evidence.

The codebase is not yet proven faster for every possible `lesgo.conf`
combination.  The remaining gaps are mostly optional or less common paths:
`USE_DYN_TN=ON`, production-size scalar transport, production-size CPS with
scalar coupling, IWM-heavy wall-model cases, shifted inflow, sponge/Coriolis
special cases, and CGNS output compatibility.  These are the next validation
targets, excluding LVLSET.

## Recommended Next Checks

1. Scalar CPU/GPU gate:
   `USE_SCALARS=ON`, `USE_SCALARS_GPU=ON`, `USE_LES_GPU=ON`, 128^3 and 240^3,
   passive and active-buoyancy cases.
2. CPS gate:
   `USE_CPS=ON`, 240^3, multi-rank/multi-GPU, with and without scalar coupling.
3. Dynamic-timescale SGS gate:
   `USE_DYN_TN=ON`, `sgs_model=4` and `sgs_model=5`, CPU/GPU correctness and
   timing separation for regular steps and SGS-update steps.
4. Wall-model gate:
   IWM/lower-wall and upper-wall combinations at production-enough grid size.
5. Optional forcing/inflow gate:
   shifted inflow, sponge, and Coriolis configurations used by real users.
