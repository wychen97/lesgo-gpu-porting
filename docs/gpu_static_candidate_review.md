# Static GPU Candidate Review

This generated report lists every tracked non-LVLSET Fortran
subroutine/function that the static scanner classifies as an
`unmarked-runtime-candidate`.  These entries are not automatically
missing GPU kernels: many are host-model logic, fallback compatibility,
diagnostics, or setup/control code.  The validation rows show where
runtime correctness and timing evidence must close each bucket.

| Review bucket | Candidates | Validation rows | Meaning |
| --- | ---: | --- | --- |
| `adm-cpu-fallback-profile` | 5 | `adm_disk`, `adm_dynamic_controls` | ADM/turbine CPU fallback or compatibility routines; profile before treating them as missing GPU work. |
| `atm-host-model` | 35 | `atm_line`, `large_windfarm` | ATM blade, controller, structure, and small math helpers that remain host-side in the current hybrid design. |
| `atm-mirror-lb-control` | 11 | `atm_line`, `large_windfarm` | ATM mirror, synchronization, cell-search, and load-balance control helpers around the GPU sampling/forcing path. |
| `cpu-fallback-compat` | 27 | `les_core_channel`, `hit_inflow` | CPU fallback or host compatibility routines retained beside GPU production paths. |
| `diagnostic-profiling` | 8 | `diagnostics_output` | Profiling, timing, or audit helpers; not GPU hot-path kernels. |
| `excluded-lvlset-bridge` | 1 | `lvlset` | LVLSET bridge code excluded from the current non-LVLSET scope. |
| `generic-helper-profile` | 10 | `les_core_channel`, `adm_disk`, `atm_line` | Generic interpolation/math helpers that may be CPU fallback or low-cost support code. |
| `inflow-fringe-profile` | 2 | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` | Inflow/fringe helpers that need targeted runtime validation for nonstandard inflow configurations. |
| `iwm-wallmodel-profile` | 1 | `iwm_wall_model` | IWM wall-model candidate that needs an IWM-heavy correctness and timing case before broad speed claims. |
| `scalar-init-fallback` | 7 | `scalar_passive`, `scalar_active`, `cps_scalar` | Scalar initialization, stability helper, or CPU fallback routines; validate passive and active scalar cases separately. |

| File:line | Subprogram | Review bucket | Validation rows |
| --- | --- | --- | --- |
| `actuator_turbine_model.f90:710` | `atm_create_points` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:899` | `atm_update` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:928` | `atm_control_yaw` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:959` | `atm_computeRotorSpeed` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1140` | `atm_rotateBlades` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1245` | `atm_compute_cl_correction` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1406` | `s_fit` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1443` | `atm_calculate_variables` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1613` | `atm_airfoil_blend_info` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1687` | `atm_computeBladeForce` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1931` | `atm_computeNacelleForce` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1980` | `atm_integrate_u` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2003` | `atm_yawNacelle` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2290` | `atm_compute_power` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2408` | `atm_solve_structure` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2701` | `solve_linear_system_banded_dp` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2764` | `solve_linear_system_dp` | `atm-host-model` | `atm_line`, `large_windfarm` |
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
| `atm_input_util.f90:1113` | `readline` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_input_util.f90:1158` | `eat_whitespace` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:590` | `atm_lesgo_force_gpu_atpoint` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:728` | `atm_lesgo_build_force_shadows` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:825` | `atm_lesgo_build_blade_mirrors` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:880` | `atm_sync_blade_points_to_device` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:898` | `atm_sync_blade_forces_to_device` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:917` | `atm_lesgo_findCells` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:1465` | `atm_lesgo_mpi_gather` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:1647` | `atm_lesgo_mpi_gather_packed` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:2060` | `atm_lesgo_compute_Spalart_u` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `convec.f90:24` | `convec` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:49` | `stress_uv_xy_derivs` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:70` | `stress_w_xy_derivs` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:88` | `ddx` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:131` | `ddy` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:174` | `ddxy` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:220` | `filt_da` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:323` | `ddz_uv` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:377` | `ddz_w` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `divstress_uv.f90:22` | `divstress_uv` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `divstress_w.f90:22` | `divstress_w` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:56` | `mul_real_complex_imag_scalar` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:93` | `mul_real_complex_2D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:154` | `mul_real_complex_imag_2D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:213` | `mul_real_complex_real_2D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:43` | `padd` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:74` | `unpadd` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `forcing.f90:100` | `lvlset_bridge_time` | `excluded-lvlset-bridge` | `lvlset` |
| `fringe.f90:59` | `constructor` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `functions.f90:62` | `interp_to_uv_grid` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:139` | `interp_to_w_grid` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:598` | `bilinear_interp_sa_nocheck` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:634` | `bilinear_interp_sa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:661` | `bilinear_interp_aa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:719` | `linear_interp_sa_nocheck` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:745` | `linear_interp_sa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:771` | `linear_interp_aa` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:806` | `cross_product` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:818` | `binary_search` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `hit_inflow.f90:153` | `extract_HIT_data` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `hit_inflow.f90:320` | `interpolate3D` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `inflow.f90:108` | `apply_inflow` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `interpolag_Sdep.f90:21` | `interpolag_Sdep` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `interpolag_Ssim.f90:21` | `interpolag_Ssim` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `iwmles.f90:633` | `iwm_slv` | `iwm-wallmodel-profile` | `iwm_wall_model` |
| `lagrange_Sdep.f90:21` | `lagrange_Sdep` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `lagrange_Ssim.f90:21` | `lagrange_Ssim` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `linear_simple.f90:46` | `solve_linear` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:86` | `assert_eq2` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:106` | `assert_eq3` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:126` | `assert_eq4` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:154` | `ludcmp` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:226` | `lubksb` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:278` | `outerprod` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:288` | `swap` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `mpi_transpose_mod.f90:36` | `mpi_transpose` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `scalars.f90:867` | `ic_scal_file` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:887` | `ic_scal_les` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:906` | `ic_scal_interp` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1726` | `to_big` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scaledep_dynamic.f90:21` | `scaledep_dynamic` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `sgs_stag_util.f90:47` | `sgs_stag` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:489` | `calc_Sij` | `diagnostic-profiling` | `diagnostics_output` |
| `stability.f90:2` | `stability` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:60` | `calc_psi_m` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:96` | `calc_psi_h` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `std_dynamic.f90:21` | `std_dynamic` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `test_filtermodule.f90:151` | `test_filter` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:175` | `test_filter_3` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:190` | `test_filter_6` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:208` | `test_test_filter` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:232` | `test_test_filter_3` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:247` | `test_test_filter_6` | `diagnostic-profiling` | `diagnostics_output` |
| `tridag_array.f90:251` | `tridag_array` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `turbine_indicator.f90:53` | `val` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:316` | `turbines_nodes` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:992` | `turbines_forcing_acc` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:1002` | `turbines_forcing` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:1267` | `place_turbines` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |

Regenerate this file with:

```bash
python3 tools/report_gpu_static_candidate_review.py --write
```
