# Static GPU Full Inventory

This generated report lists every tracked root-level non-LVLSET
Fortran subroutine/function and its static GPU-audit classification.
It is a source heuristic, not proof of runtime speedup.  Use the
validation matrix and evidence ledger for correctness and timing
claims.

| Classification | Subprograms |
| --- | ---: |
| `gpu-marked` | 147 |
| `gpu-file-unmarked` | 2 |
| `host-boundary` | 102 |
| `host-or-diagnostic` | 41 |
| `unmarked-runtime-candidate` | 107 |

| File:line | Subprogram | Classification | Review bucket | Validation rows |
| --- | --- | --- | --- | --- |
| `actuator_turbine_model.f90:97` | `atm_model_env_token` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:207` | `atm_structure_timing_report` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:226` | `atm_structure_diag_snapshot` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:302` | `atm_initialize` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:347` | `atm_read_actuator_points` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:399` | `atm_read_restart` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:480` | `atm_write_restart` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:574` | `atm_initialize_output` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:708` | `atm_create_points` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:897` | `atm_update` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:926` | `atm_control_yaw` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:957` | `atm_computeRotorSpeed` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1138` | `atm_rotateBlades` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1243` | `atm_compute_cl_correction` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1404` | `s_fit` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1441` | `atm_calculate_variables` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1611` | `atm_airfoil_blend_info` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1685` | `atm_computeBladeForce` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1929` | `atm_computeNacelleForce` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1978` | `atm_integrate_u` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2001` | `atm_yawNacelle` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2079` | `atm_output` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:2288` | `atm_compute_power` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2315` | `atm_write_blade_points` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:2348` | `atm_process_output` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:2406` | `atm_solve_structure` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2699` | `solve_linear_system_banded_dp` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2762` | `solve_linear_system_dp` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:40` | `error` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:57` | `interpolate` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:122` | `vector_add` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:134` | `vector_divide` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:147` | `vector_multiply` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:161` | `vector_mag` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:171` | `rotatePoint` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:208` | `matrix_vector` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:222` | `cross_product` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_base.f90:234` | `distance` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `atm_input_util.f90:359` | `read_input_conf` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:524` | `read_turbine_model_variables` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:1059` | `atm_print_initialize` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:1074` | `read_airfoil` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:1113` | `readline` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_input_util.f90:1158` | `eat_whitespace` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:245` | `atm_interp_w_to_uv` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:306` | `atm_prepare_direct_w` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:390` | `atm_lesgo_apply_force_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:435` | `atm_lesgo_convolute_force_gpu_atpoint` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:590` | `atm_lesgo_force_gpu_atpoint` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:602` | `atm_lesgo_initialize` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:650` | `atm_lesgo_finalize` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:685` | `atm_lesgo_report_timing` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:728` | `atm_lesgo_build_force_shadows` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:805` | `atm_lesgo_destroy_force_shadows` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:825` | `atm_lesgo_build_blade_mirrors` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:866` | `atm_lesgo_destroy_blade_mirrors` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:880` | `atm_sync_blade_points_to_device` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:898` | `atm_sync_blade_forces_to_device` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:917` | `atm_lesgo_findCells` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:1115` | `atm_lesgo_forcing` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:1465` | `atm_lesgo_mpi_gather` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:1647` | `atm_lesgo_mpi_gather_packed` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:1903` | `atm_lesgo_force` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2060` | `atm_lesgo_compute_Spalart_u` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:2142` | `atm_lesgo_convolute_force` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2430` | `atm_convolute_atpoint_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2581` | `atm_sample_velocity_atpoint_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2697` | `atm_batch_atpoint_init` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2775` | `atm_batch_sample_velocity_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2963` | `atm_batch_convolute_force_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3116` | `atm_batch_clc_init` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3184` | `atm_batch_cl_correction_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3440` | `atm_lesgo_apply_force` | `gpu-marked` |  |  |
| `cfl_util.f90:38` | `get_max_cfl` | `gpu-marked` |  |  |
| `cfl_util.f90:98` | `get_cfl_dt` | `gpu-marked` |  |  |
| `clocks.f90:47` | `start` | `host-boundary` |  |  |
| `clocks.f90:65` | `stop` | `host-boundary` |  |  |
| `concurrent_precursor.f90:102` | `cps_timer_start` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:116` | `cps_timer_accum` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:133` | `cps_stage_report` | `host-or-diagnostic` |  |  |
| `concurrent_precursor.f90:167` | `initialize_cps` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:222` | `synchronize_cps` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:434` | `inflow_cps` | `gpu-marked` |  |  |
| `convec.f90:24` | `convec` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `convec_gpu.f90:71` | `convec_gpu_init` | `gpu-marked` |  |  |
| `convec_gpu.f90:108` | `convec_gpu_finalize` | `gpu-marked` |  |  |
| `convec_gpu.f90:120` | `convec_gpu` | `gpu-marked` |  |  |
| `convec_gpu.f90:540` | `unpadd_gpu` | `gpu-marked` |  |  |
| `coriolis.f90:84` | `coriolis_init` | `host-or-diagnostic` |  |  |
| `coriolis.f90:164` | `coriolis_finalize` | `host-or-diagnostic` |  |  |
| `coriolis.f90:177` | `coriolis_calc` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:38` | `cuda_mpi_debug_init` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:89` | `print_env` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:125` | `cuda_pre` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:139` | `cuda_post` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:154` | `probe_real` | `host-boundary` |  |  |
| `cuda_mpi_debug.f90:164` | `probe_complex` | `host-boundary` |  |  |
| `cuda_mpi_debug.f90:174` | `mpi_dbg_sendrecv_r` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:212` | `mpi_dbg_sendrecv_c` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:250` | `mpi_dbg_send_r` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:281` | `mpi_dbg_recv_r` | `gpu-marked` |  |  |
| `derivatives.f90:49` | `stress_uv_xy_derivs` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:70` | `stress_w_xy_derivs` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:88` | `ddx` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:131` | `ddy` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:174` | `ddxy` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:220` | `filt_da` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:269` | `filt_da_vel` | `gpu-marked` |  |  |
| `derivatives.f90:294` | `ddz_vel` | `gpu-marked` |  |  |
| `derivatives.f90:323` | `ddz_uv` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives.f90:377` | `ddz_w` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `derivatives_gpu.f90:64` | `filt_da_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:157` | `ddz_uv_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:219` | `ddz_w_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:272` | `ddx_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:343` | `ddy_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:406` | `ddxy_gpu` | `gpu-marked` |  |  |
| `divstress_uv.f90:22` | `divstress_uv` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `divstress_w.f90:22` | `divstress_w` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:56` | `mul_real_complex_imag_scalar` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:93` | `mul_real_complex_2D` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:154` | `mul_real_complex_imag_2D` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `emul_complex.f90:213` | `mul_real_complex_real_2D` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:43` | `padd` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:74` | `unpadd` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `fft.f90:102` | `init_fft` | `host-or-diagnostic` |  |  |
| `fft.f90:131` | `init_wavenumber` | `gpu-marked` |  |  |
| `fft_gpu.f90:79` | `init_fft_gpu` | `gpu-marked` |  |  |
| `fft_gpu.f90:136` | `make_plan` | `gpu-marked` |  |  |
| `fft_gpu.f90:171` | `set_stream` | `gpu-marked` |  |  |
| `fft_gpu.f90:182` | `finalize_fft_gpu` | `gpu-file-unmarked` |  |  |
| `fft_gpu.f90:194` | `destroy_plan` | `gpu-marked` |  |  |
| `fft_gpu.f90:211` | `fft_gpu_exec_d2z` | `gpu-marked` |  |  |
| `fft_gpu.f90:222` | `fft_gpu_exec_z2d` | `gpu-marked` |  |  |
| `fft_gpu.f90:234` | `check_cufft` | `gpu-marked` |  |  |
| `finalize.f90:21` | `finalize` | `host-boundary` |  |  |
| `forcing.f90:100` | `lvlset_bridge_time` | `unmarked-runtime-candidate` | `excluded-lvlset-bridge` | `lvlset` |
| `forcing.f90:118` | `lvlset_bridge_report` | `host-or-diagnostic` |  |  |
| `forcing.f90:160` | `forcing_random` | `gpu-marked` |  |  |
| `forcing.f90:230` | `forcing_applied` | `gpu-marked` |  |  |
| `forcing.f90:323` | `forcing_induced` | `gpu-marked` |  |  |
| `forcing.f90:380` | `project` | `gpu-marked` |  |  |
| `fringe.f90:60` | `constructor` | `unmarked-runtime-candidate` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `functions.f90:62` | `interp_to_uv_grid` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:140` | `interp_to_w_grid` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:602` | `bilinear_interp_sa_nocheck` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:638` | `bilinear_interp_sa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:665` | `bilinear_interp_aa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:723` | `linear_interp_sa_nocheck` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:749` | `linear_interp_sa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:775` | `linear_interp_aa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:810` | `cross_product` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:822` | `binary_search` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:893` | `get_tau_wall_bot` | `gpu-marked` |  |  |
| `functions.f90:934` | `get_tau_wall_top` | `gpu-marked` |  |  |
| `functions.f90:974` | `count_lines` | `host-or-diagnostic` |  |  |
| `grid.f90:44` | `build` | `host-boundary` |  |  |
| `hit_inflow.f90:82` | `initialize_HIT` | `gpu-marked` |  |  |
| `hit_inflow.f90:153` | `extract_HIT_data` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `hit_inflow.f90:196` | `compute_HIT_plane_data` | `gpu-marked` |  |  |
| `hit_inflow.f90:242` | `inflow_HIT` | `gpu-marked` |  |  |
| `hit_inflow.f90:273` | `hit_write_restart` | `host-or-diagnostic` |  |  |
| `hit_inflow.f90:293` | `hit_read_restart` | `host-or-diagnostic` |  |  |
| `hit_inflow.f90:320` | `interpolate3D` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `hit_inflow_gpu.f90:62` | `hit_gpu_setup` | `gpu-marked` |  |  |
| `hit_inflow_gpu.f90:110` | `hit_fringe_setup_gpu` | `gpu-marked` |  |  |
| `hit_inflow_gpu.f90:137` | `hit_compute_plane_gpu` | `gpu-marked` |  |  |
| `hit_inflow_gpu.f90:305` | `hit_apply_fringe_gpu` | `gpu-marked` |  |  |
| `inflow.f90:64` | `inflow_init` | `gpu-marked` |  |  |
| `inflow.f90:108` | `apply_inflow` | `unmarked-runtime-candidate` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `inflow.f90:134` | `inflow_uniform` | `gpu-marked` |  |  |
| `init_random_seed.f90:21` | `init_random_seed` | `host-boundary` |  |  |
| `initial.f90:21` | `initial` | `host-boundary` |  |  |
| `initial.f90:203` | `check_for_interp` | `host-boundary` |  |  |
| `initial.f90:228` | `ic_file` | `host-boundary` |  |  |
| `initial.f90:244` | `ic_interp` | `host-boundary` |  |  |
| `initial.f90:385` | `ic_dns` | `host-boundary` |  |  |
| `initial.f90:475` | `ic_les` | `host-boundary` |  |  |
| `initialize.f90:21` | `initialize` | `gpu-marked` |  |  |
| `input_util.f90:63` | `read_input_conf` | `host-boundary` |  |  |
| `input_util.f90:258` | `model_block` | `host-boundary` |  |  |
| `input_util.f90:316` | `coriolis_block` | `host-boundary` |  |  |
| `input_util.f90:375` | `time_block` | `host-boundary` |  |  |
| `input_util.f90:432` | `flow_cond_block` | `host-boundary` |  |  |
| `input_util.f90:564` | `output_block` | `host-boundary` |  |  |
| `input_util.f90:669` | `level_set_block` | `host-boundary` |  |  |
| `input_util.f90:752` | `turbines_block` | `host-boundary` |  |  |
| `input_util.f90:832` | `scalars_block` | `host-boundary` |  |  |
| `input_util.f90:894` | `checkentry` | `host-boundary` |  |  |
| `input_util.f90:909` | `readline` | `host-boundary` |  |  |
| `input_util.f90:950` | `parse_vector_real` | `host-boundary` |  |  |
| `input_util.f90:985` | `parse_vector_point3D` | `host-boundary` |  |  |
| `interpolag_Sdep.f90:21` | `interpolag_Sdep` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `interpolag_Ssim.f90:21` | `interpolag_Ssim` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `io.f90:70` | `openfiles` | `host-boundary` |  |  |
| `io.f90:109` | `energy` | `gpu-marked` |  |  |
| `io.f90:173` | `write_tau_wall_bot` | `host-boundary` |  |  |
| `io.f90:200` | `write_tau_wall_top` | `host-boundary` |  |  |
| `io.f90:229` | `write_parallel_cgns` | `host-boundary` |  |  |
| `io.f90:398` | `write_null_cgns` | `host-boundary` |  |  |
| `io.f90:563` | `output_loop` | `host-boundary` |  |  |
| `io.f90:710` | `inst_write` | `gpu-marked` |  |  |
| `io.f90:1534` | `checkpoint` | `gpu-marked` |  |  |
| `io.f90:1640` | `output_final` | `host-boundary` |  |  |
| `io.f90:1654` | `output_init` | `host-boundary` |  |  |
| `iwmles.f90:115` | `iwm_wallstress` | `gpu-marked` |  |  |
| `iwmles.f90:174` | `iwm_init` | `gpu-marked` |  |  |
| `iwmles.f90:297` | `iwm_finalize` | `gpu-marked` |  |  |
| `iwmles.f90:352` | `iwm_calc_lhs` | `gpu-marked` |  |  |
| `iwmles.f90:630` | `iwm_slv` | `unmarked-runtime-candidate` | `iwm-wallmodel-profile` | `iwm_wall_model` |
| `iwmles.f90:655` | `iwm_calc_wallstress` | `gpu-marked` |  |  |
| `iwmles.f90:1164` | `iwm_monitor` | `gpu-marked` |  |  |
| `iwmles.f90:1206` | `iwm_checkPoint` | `gpu-marked` |  |  |
| `iwmles.f90:1240` | `iwm_read_checkPoint` | `gpu-marked` |  |  |
| `lagrange_Sdep.f90:21` | `lagrange_Sdep` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `lagrange_Sdep_gpu.f90:135` | `lagrange_Sdep_gpu_init` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:215` | `lagrange_Ssim_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:488` | `lagrange_Sdep_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:1084` | `interpolag_Ssim_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:1251` | `interpolag_Sdep_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:1898` | `sync_downup_F` | `gpu-marked` |  |  |
| `lagrange_Ssim.f90:21` | `lagrange_Ssim` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `linear_simple.f90:46` | `solve_linear` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:86` | `assert_eq2` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:106` | `assert_eq3` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:126` | `assert_eq4` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:154` | `ludcmp` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:226` | `lubksb` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:278` | `outerprod` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:288` | `swap` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `main.f90:1053` | `main_read_env_real` | `host-or-diagnostic` |  |  |
| `messages.f90:63` | `message_a` | `host-boundary` |  |  |
| `messages.f90:73` | `message_ai` | `host-boundary` |  |  |
| `messages.f90:84` | `message_aiai` | `host-boundary` |  |  |
| `messages.f90:95` | `message_aiar` | `host-boundary` |  |  |
| `messages.f90:108` | `message_al` | `host-boundary` |  |  |
| `messages.f90:119` | `message_aii` | `host-boundary` |  |  |
| `messages.f90:131` | `message_air` | `host-boundary` |  |  |
| `messages.f90:144` | `message_ai_array` | `host-boundary` |  |  |
| `messages.f90:158` | `message_aiai_array` | `host-boundary` |  |  |
| `messages.f90:173` | `message_ar` | `host-boundary` |  |  |
| `messages.f90:184` | `message_ar_array` | `host-boundary` |  |  |
| `messages.f90:197` | `message_aiar_array` | `host-boundary` |  |  |
| `messages.f90:212` | `warn` | `host-boundary` |  |  |
| `messages.f90:225` | `error_a` | `host-boundary` |  |  |
| `messages.f90:241` | `error_ai` | `host-boundary` |  |  |
| `messages.f90:258` | `error_ai_array` | `host-boundary` |  |  |
| `messages.f90:279` | `error_aia` | `host-boundary` |  |  |
| `messages.f90:296` | `error_aiai` | `host-boundary` |  |  |
| `messages.f90:313` | `error_aiar` | `host-boundary` |  |  |
| `messages.f90:331` | `error_arar` | `host-boundary` |  |  |
| `messages.f90:348` | `error_al` | `host-boundary` |  |  |
| `messages.f90:365` | `error_ar` | `host-boundary` |  |  |
| `messages.f90:382` | `error_ar_array` | `host-boundary` |  |  |
| `mpi_defs.f90:52` | `initialize_mpi` | `gpu-marked` |  |  |
| `mpi_defs.f90:196` | `create_mpi_comms_cps` | `host-boundary` |  |  |
| `mpi_defs.f90:246` | `mpi_sync_real_array` | `gpu-marked` |  |  |
| `mpi_defs.f90:351` | `sync_up` | `host-boundary` |  |  |
| `mpi_defs.f90:362` | `sync_downup_nb` | `host-boundary` |  |  |
| `mpi_transpose_mod.f90:36` | `mpi_transpose` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `param_output.f90:21` | `param_output` | `host-boundary` |  |  |
| `pid.f90:49` | `constructor` | `host-boundary` |  |  |
| `pid.f90:77` | `advance_noset` | `host-boundary` |  |  |
| `pid.f90:98` | `advance_set` | `host-boundary` |  |  |
| `press_gpu.f90:82` | `press_gpu_init` | `gpu-marked` |  |  |
| `press_gpu.f90:105` | `press_gpu_finalize` | `gpu-marked` |  |  |
| `press_gpu.f90:118` | `press_stag_array_gpu` | `gpu-marked` |  |  |
| `press_stag_array.f90:44` | `press_stag_array` | `gpu-marked` |  |  |
| `press_stag_array.f90:535` | `press_apply_env_enabled_unless_false` | `host-or-diagnostic` |  |  |
| `rmsdiv.f90:21` | `rmsdiv` | `gpu-marked` |  |  |
| `scalars.f90:165` | `scalars_acc_sync` | `gpu-marked` |  |  |
| `scalars.f90:191` | `scalars_timer_start` | `host-or-diagnostic` |  |  |
| `scalars.f90:203` | `scalars_timer_accum` | `host-or-diagnostic` |  |  |
| `scalars.f90:220` | `scalars_stage_report` | `host-or-diagnostic` |  |  |
| `scalars.f90:248` | `scalars_deriv_xy_big_acc` | `gpu-marked` |  |  |
| `scalars.f90:379` | `scalars_to_big_acc` | `gpu-marked` |  |  |
| `scalars.f90:428` | `scalars_return_rhs_acc` | `gpu-marked` |  |  |
| `scalars.f90:468` | `scalars_divergence_acc` | `gpu-marked` |  |  |
| `scalars.f90:528` | `scalars_update_device_state` | `gpu-marked` |  |  |
| `scalars.f90:540` | `scalars_copy_rhs_acc` | `gpu-marked` |  |  |
| `scalars.f90:559` | `scalars_advective_acc` | `gpu-marked` |  |  |
| `scalars.f90:609` | `scalars_flux_acc` | `gpu-marked` |  |  |
| `scalars.f90:693` | `scalars_rhs_theta_acc` | `gpu-marked` |  |  |
| `scalars.f90:745` | `scalars_init` | `gpu-marked` |  |  |
| `scalars.f90:829` | `ic_scal` | `gpu-marked` |  |  |
| `scalars.f90:868` | `ic_scal_file` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:887` | `ic_scal_les` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:906` | `ic_scal_interp` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1027` | `scalars_checkpoint` | `gpu-marked` |  |  |
| `scalars.f90:1047` | `scalars_deriv` | `gpu-marked` |  |  |
| `scalars.f90:1126` | `obukhov` | `gpu-marked` |  |  |
| `scalars.f90:1368` | `scalars_transport` | `gpu-marked` |  |  |
| `scalars.f90:1720` | `to_big` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1744` | `buoyancy_force` | `gpu-marked` |  |  |
| `scaledep_dynamic.f90:21` | `scaledep_dynamic` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `sgs_gpu.f90:109` | `divstress_gpu_init` | `gpu-marked` |  |  |
| `sgs_gpu.f90:126` | `std_dynamic_pples_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:278` | `scaledep_dynamic_pples_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:554` | `sgs_stag_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:915` | `calc_Sij_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:1082` | `divstress_uv_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:1185` | `divstress_w_gpu` | `gpu-marked` |  |  |
| `sgs_param.f90:111` | `sgs_param_init` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:47` | `sgs_stag` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:489` | `calc_Sij` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:665` | `sgs_stage_report` | `host-or-diagnostic` |  |  |
| `shifted_inflow.f90:68` | `shifted_inflow_init` | `gpu-marked` |  |  |
| `shifted_inflow.f90:111` | `inflow_shifted` | `gpu-marked` |  |  |
| `sim_param.f90:69` | `sim_param_init` | `gpu-marked` |  |  |
| `sponge.f90:43` | `sponge_init` | `gpu-marked` |  |  |
| `sponge.f90:67` | `sponge_force` | `gpu-marked` |  |  |
| `stability.f90:2` | `stability` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:60` | `calc_psi_m` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:96` | `calc_psi_h` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `std_dynamic.f90:21` | `std_dynamic` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `string_util.f90:68` | `numtostr_r` | `host-boundary` |  |  |
| `string_util.f90:101` | `numtostr_i` | `host-boundary` |  |  |
| `string_util.f90:125` | `eat_whitespace` | `host-boundary` |  |  |
| `string_util.f90:157` | `uppercase` | `host-boundary` |  |  |
| `string_util.f90:196` | `split_string` | `host-boundary` |  |  |
| `string_util.f90:276` | `count_string_occur` | `host-boundary` |  |  |
| `string_util.f90:313` | `string_concat_a` | `host-boundary` |  |  |
| `string_util.f90:326` | `string_concat_r` | `host-boundary` |  |  |
| `string_util.f90:341` | `string_concat_i` | `host-boundary` |  |  |
| `string_util.f90:356` | `string_concat_ai` | `host-boundary` |  |  |
| `string_util.f90:371` | `string_concat_ar` | `host-boundary` |  |  |
| `string_util.f90:386` | `string_concat_aia` | `host-boundary` |  |  |
| `string_util.f90:402` | `string_concat_ara` | `host-boundary` |  |  |
| `string_util.f90:418` | `string_concat_aiaia` | `host-boundary` |  |  |
| `string_util.f90:436` | `string_concat_arara` | `host-boundary` |  |  |
| `string_util.f90:454` | `string_concat_aiai` | `host-boundary` |  |  |
| `string_util.f90:471` | `string_concat_arar` | `host-boundary` |  |  |
| `string_util.f90:488` | `string_concat_araia` | `host-boundary` |  |  |
| `string_util.f90:507` | `string_concat_arai` | `host-boundary` |  |  |
| `string_util.f90:525` | `string_concat_aiaiai` | `host-boundary` |  |  |
| `string_util.f90:544` | `string_concat_ararar` | `host-boundary` |  |  |
| `string_util.f90:563` | `string_concat_aiaiaia` | `host-boundary` |  |  |
| `string_util.f90:583` | `string_concat_ararara` | `host-boundary` |  |  |
| `string_util.f90:607` | `string_splice_aa` | `host-boundary` |  |  |
| `string_util.f90:619` | `string_splice_ar` | `host-boundary` |  |  |
| `string_util.f90:636` | `string_splice_ai` | `host-boundary` |  |  |
| `string_util.f90:652` | `string_splice_aia` | `host-boundary` |  |  |
| `string_util.f90:668` | `string_splice_ara` | `host-boundary` |  |  |
| `string_util.f90:685` | `string_splice_aiai` | `host-boundary` |  |  |
| `string_util.f90:702` | `string_splice_arar` | `host-boundary` |  |  |
| `string_util.f90:720` | `string_splice_aiaia` | `host-boundary` |  |  |
| `string_util.f90:737` | `string_splice_arara` | `host-boundary` |  |  |
| `string_util.f90:755` | `string_splice_araia` | `host-boundary` |  |  |
| `string_util.f90:774` | `string_splice_arai` | `host-boundary` |  |  |
| `string_util.f90:793` | `string_splice_aiaiai` | `host-boundary` |  |  |
| `string_util.f90:811` | `string_splice_ararar` | `host-boundary` |  |  |
| `string_util.f90:831` | `string_splice_aiaiaia` | `host-boundary` |  |  |
| `string_util.f90:850` | `string_splice_ararara` | `host-boundary` |  |  |
| `test_filtermodule.f90:57` | `test_filter_init` | `gpu-marked` |  |  |
| `test_filtermodule.f90:151` | `test_filter` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:175` | `test_filter_3` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:190` | `test_filter_6` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:208` | `test_test_filter` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:232` | `test_test_filter_3` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:247` | `test_test_filter_6` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `test_filtermodule.f90:266` | `test_filter_b_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:301` | `test_test_filter_b_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:335` | `test_filter_plane_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:362` | `test_test_filter_plane_gpu` | `gpu-marked` |  |  |
| `time_average.f90:75` | `init` | `gpu-marked` |  |  |
| `time_average.f90:202` | `compute` | `gpu-marked` |  |  |
| `time_average.f90:463` | `finalize` | `gpu-marked` |  |  |
| `time_average.f90:816` | `checkpoint` | `gpu-marked` |  |  |
| `time_average.f90:882` | `write_parallel_cgns` | `host-or-diagnostic` |  |  |
| `time_average.f90:1052` | `write_null_cgns` | `host-or-diagnostic` |  |  |
| `tridag_array.f90:43` | `tridag_array` | `gpu-marked` |  |  |
| `tridag_array.f90:194` | `tridag_apply_env_enabled_unless_false` | `host-or-diagnostic` |  |  |
| `tridag_array.f90:222` | `tridag_apply_env_true_token` | `host-or-diagnostic` |  |  |
| `tridag_array.f90:250` | `tridag_array` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `tridag_gpu.f90:107` | `tridag_gpu_init` | `gpu-marked` |  |  |
| `tridag_gpu.f90:126` | `tridag_gpu_finalize` | `gpu-marked` |  |  |
| `tridag_gpu.f90:136` | `tridag_array_gpu` | `gpu-marked` |  |  |
| `turbine_indicator.f90:53` | `val` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbine_indicator.f90:77` | `init` | `host-or-diagnostic` |  |  |
| `turbines.f90:160` | `turbines_init` | `gpu-marked` |  |  |
| `turbines.f90:316` | `turbines_nodes` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:516` | `turbines_acc_metadata_init` | `gpu-marked` |  |  |
| `turbines.f90:662` | `turbines_acc_finalize` | `host-or-diagnostic` |  |  |
| `turbines.f90:703` | `turbines_acc_sync_device_field` | `gpu-marked` |  |  |
| `turbines.f90:736` | `turbines_forcing_acc` | `gpu-marked` |  |  |
| `turbines.f90:989` | `turbines_forcing_acc` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:999` | `turbines_forcing` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:1175` | `turbines_finalize` | `gpu-marked` |  |  |
| `turbines.f90:1195` | `turbines_checkpoint` | `host-or-diagnostic` |  |  |
| `turbines.f90:1220` | `turbine_vel_init` | `host-or-diagnostic` |  |  |
| `turbines.f90:1260` | `place_turbines` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:1410` | `read_control_files` | `host-or-diagnostic` |  |  |
| `turbines_gpu.f90:35` | `turbines_interp_w_to_uv_gpu` | `gpu-file-unmarked` |  |  |
| `wallstress.f90:21` | `wallstress` | `gpu-marked` |  |  |
| `wallstress.f90:137` | `ws_free_ubc` | `gpu-marked` |  |  |
| `wallstress.f90:164` | `ws_dns_lbc` | `gpu-marked` |  |  |
| `wallstress.f90:197` | `ws_dns_ubc` | `gpu-marked` |  |  |
| `wallstress.f90:230` | `ws_equilibrium_lbc` | `gpu-marked` |  |  |
| `wallstress.f90:345` | `ws_equilibrium_ubc` | `gpu-marked` |  |  |

Regenerate this file with:

```bash
python3 tools/report_gpu_static_full_inventory.py --write
```
