# Static GPU Candidate Review

This generated report lists every tracked non-LVLSET Fortran
subroutine/function that the static scanner classifies as an
`unmarked-runtime-candidate`.  These entries are not automatically
missing GPU kernels: many are host-model logic, fallback compatibility,
diagnostics, or setup/control code.  The validation rows show where
runtime correctness and timing evidence must close each bucket.

| Review bucket | Candidates | Validation rows | Meaning |
| --- | ---: | --- | --- |
| `adm-cpu-fallback-profile` | 4 | `adm_disk`, `adm_dynamic_controls` | ADM/turbine CPU fallback or compatibility routines; profile before treating them as missing GPU work. |
| `atm-host-model` | 30 | `atm_line`, `large_windfarm` | ATM blade, controller, structure, and small math helpers that remain host-side in the current hybrid design. |
| `atm-mirror-lb-control` | 11 | `atm_line`, `large_windfarm` | ATM mirror, synchronization, cell-search, and load-balance control helpers around the GPU sampling/forcing path. |
| `cpu-fallback-compat` | 10 | `les_core_channel`, `hit_inflow` | CPU fallback or host compatibility routines retained beside GPU production paths. |
| `diagnostic-profiling` | 14 | `diagnostics_output` | Profiling, timing, or audit helpers; not GPU hot-path kernels. |
| `excluded-lvlset-bridge` | 1 | `lvlset` | LVLSET bridge code excluded from the current non-LVLSET scope. |
| `generic-helper-profile` | 10 | `les_core_channel`, `adm_disk`, `atm_line` | Generic interpolation/math helpers that may be CPU fallback or low-cost support code. |
| `inflow-fringe-profile` | 2 | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` | Inflow/fringe helpers that need targeted runtime validation for nonstandard inflow configurations. |
| `iwm-wallmodel-profile` | 1 | `iwm_wall_model` | IWM wall-model candidate that needs an IWM-heavy correctness and timing case before broad speed claims. |
| `scalar-init-fallback` | 8 | `scalar_passive`, `scalar_active`, `cps_scalar` | Scalar initialization, stability helper, or CPU fallback routines; validate passive and active scalar cases separately. |

| File:line | Subprogram | Review bucket | Validation rows |
| --- | --- | --- | --- |
| `actuator_turbine_model.f90:828` | `atm_create_points` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1017` | `atm_update` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1046` | `atm_control_yaw` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1077` | `atm_computeRotorSpeed` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1959` | `s_fit` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1996` | `atm_calculate_variables` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2166` | `atm_airfoil_blend_info` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2240` | `atm_computeBladeForce` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2484` | `atm_computeNacelleForce` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2933` | `atm_compute_power` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:3540` | `solve_linear_system_banded_dp` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:3603` | `solve_linear_system_dp` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:40` | `error` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:57` | `interpolate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:122` | `vector_add` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:134` | `vector_divide` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:147` | `vector_multiply` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:161` | `vector_mag` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:171` | `rotatePoint` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:208` | `matrix_vector` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:222` | `cross_product` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:234` | `distance` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_input_util.f90:1274` | `readline` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_input_util.f90:1319` | `eat_whitespace` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:607` | `atm_lb_auto_record` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:2525` | `atm_lb_refresh_targeted_sample_slots` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:3459` | `atm_lesgo_lb_plan` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:3937` | `atm_lesgo_build_force_shadows` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4034` | `atm_lesgo_build_blade_mirrors` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4089` | `atm_sync_blade_points_to_device` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4107` | `atm_sync_blade_forces_to_device` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4126` | `atm_lesgo_findCells` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:6109` | `atm_lesgo_compute_Spalart_u` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `emul_complex.f90:56` | `mul_real_complex_imag_scalar` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:93` | `mul_real_complex_2D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:154` | `mul_real_complex_imag_2D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:213` | `mul_real_complex_real_2D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:43` | `padd` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:74` | `unpadd` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `forcing.f90:204` | `lvlset_bridge_time` | `excluded-lvlset-bridge` | `lvlset` |
| `fringe.f90:68` | `constructor` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `functions.f90:96` | `interp_to_uv_grid` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:174` | `interp_to_w_grid` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:636` | `bilinear_interp_sa_nocheck` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:672` | `bilinear_interp_sa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:699` | `bilinear_interp_aa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:757` | `linear_interp_sa_nocheck` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:783` | `linear_interp_sa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:809` | `linear_interp_aa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:844` | `cross_product` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:856` | `binary_search` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `hit_inflow.f90:153` | `extract_HIT_data` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `hit_inflow.f90:320` | `interpolate3D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `inflow.f90:137` | `apply_inflow` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `iwmles.f90:830` | `iwm_slv` | `iwm-wallmodel-profile` | `iwm_wall_model` |
| `linear_simple.f90:46` | `solve_linear` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:86` | `assert_eq2` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:106` | `assert_eq3` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:126` | `assert_eq4` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:154` | `ludcmp` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:226` | `lubksb` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:278` | `outerprod` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:288` | `swap` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `mpi_transpose_mod.f90:36` | `mpi_transpose` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `scalars.f90:1133` | `stability_acc` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1634` | `ic_scal_file` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1653` | `ic_scal_les` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1672` | `ic_scal_interp` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:2624` | `to_big` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `sgs_stag_util.f90:644` | `sgs_calc_cpu_start` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:655` | `sgs_calc_cpu_stop` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:668` | `sgs_calc_set_zrange` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:680` | `sgs_calc_set_audit` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:710` | `sgs_tau_detail_begin` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:731` | `sgs_tau_detail_start` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:742` | `sgs_tau_detail_stop` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:755` | `sgs_tau_detail_add_bytes` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:767` | `sgs_tau_detail_add_msg` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:782` | `sgs_dwdz_detail_begin` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:804` | `sgs_dwdz_detail_start` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:815` | `sgs_dwdz_detail_stop` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:828` | `sgs_dwdz_detail_add_msg` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:846` | `sgs_dwdz_path_audit` | `diagnostic-profiling` | `diagnostics_output` |
| `stability.f90:2` | `stability` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:60` | `calc_psi_m` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:96` | `calc_psi_h` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `tridag_array.f90:2080` | `tridag_array` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `turbine_indicator.f90:56` | `val` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:316` | `turbines_nodes` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:989` | `turbines_forcing_acc` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:1279` | `place_turbines` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |

Regenerate this file with:

```bash
python3 tools/report_gpu_static_candidate_review.py --write
```
