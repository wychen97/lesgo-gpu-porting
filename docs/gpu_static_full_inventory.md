# Static GPU Full Inventory

This generated report lists every tracked root-level non-LVLSET
Fortran subroutine/function and its static GPU-audit classification.
It is a source heuristic, not proof of runtime speedup.  Use the
validation matrix and evidence ledger for correctness and timing
claims.

| Classification | Subprograms |
| --- | ---: |
| `gpu-marked` | 330 |
| `gpu-file-unmarked` | 1 |
| `host-boundary` | 99 |
| `host-or-diagnostic` | 52 |
| `unmarked-runtime-candidate` | 91 |

| File:line | Subprogram | Classification | Review bucket | Validation rows |
| --- | --- | --- | --- | --- |
| `actuator_turbine_model.f90:101` | `atm_model_env_token` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:211` | `atm_structure_timing_report` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:230` | `atm_structure_diag_snapshot` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:379` | `atm_model_cuda_check` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:400` | `atm_model_cuda_sync` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:422` | `atm_initialize` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:467` | `atm_read_actuator_points` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:519` | `atm_read_restart` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:600` | `atm_write_restart` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:694` | `atm_initialize_output` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:828` | `atm_create_points` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1017` | `atm_update` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1046` | `atm_control_yaw` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1077` | `atm_computeRotorSpeed` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1258` | `atm_rotateBlades` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:1417` | `atm_compute_cl_correction` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:1583` | `atm_compute_cl_correction_gpu` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:1953` | `s_fit` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:1990` | `atm_calculate_variables` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2160` | `atm_airfoil_blend_info` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2234` | `atm_computeBladeForce` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2478` | `atm_computeNacelleForce` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2527` | `atm_integrate_u` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:2583` | `atm_yawNacelle` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:2715` | `atm_output` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:2924` | `atm_compute_power` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:2951` | `atm_write_blade_points` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:2984` | `atm_process_output` | `host-or-diagnostic` |  |  |
| `actuator_turbine_model.f90:3042` | `atm_solve_structure` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:3421` | `solve_linear_system_gpu_dp` | `gpu-marked` |  |  |
| `actuator_turbine_model.f90:3531` | `solve_linear_system_banded_dp` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `actuator_turbine_model.f90:3594` | `solve_linear_system_dp` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
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
| `atm_input_util.f90:520` | `read_input_conf` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:685` | `read_turbine_model_variables` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:1220` | `atm_print_initialize` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:1235` | `read_airfoil` | `host-or-diagnostic` |  |  |
| `atm_input_util.f90:1274` | `readline` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_input_util.f90:1319` | `eat_whitespace` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:607` | `atm_lb_auto_record` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:669` | `atm_diag_event_start` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:707` | `atm_diag_event_stop` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:725` | `atm_diag_event_flush` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:750` | `atm_cuda_check` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:771` | `atm_cuda_sync` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:792` | `atm_lesgo_reset_turbine_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:864` | `atm_interp_w_to_uv` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:961` | `atm_prepare_direct_w` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:1045` | `atm_lesgo_apply_force_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:1118` | `atm_lesgo_convolute_force_gpu_atpoint` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:1405` | `atm_lesgo_force_gpu_atpoint` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:1843` | `atm_lesgo_nacelle_force_gpu_atpoint` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2054` | `atm_point_owner_sample_turbine` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2217` | `atm_point_owner_force_turbine` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2489` | `atm_lb_ensure_targeted_buffers` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:2520` | `atm_lb_refresh_targeted_sample_slots` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:2585` | `atm_lb_targeted_velocity_exchange` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2669` | `atm_lb_pack_force_turbine` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2721` | `atm_lb_unpack_force_turbine` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2773` | `atm_point_owner_lb_gather_targeted` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2878` | `atm_point_owner_lb_force` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:2942` | `atm_point_owner_lb_gather` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3021` | `atm_point_owner_lb_validate` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:3117` | `atm_lesgo_initialize` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3171` | `atm_lesgo_finalize` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3209` | `atm_lesgo_report_timing` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:3337` | `atm_lesgo_diag_load` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:3454` | `atm_lesgo_lb_plan` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:3932` | `atm_lesgo_build_force_shadows` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4009` | `atm_lesgo_destroy_force_shadows` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:4029` | `atm_lesgo_build_blade_mirrors` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4070` | `atm_lesgo_destroy_blade_mirrors` | `host-or-diagnostic` |  |  |
| `atm_lesgo_interface.f90:4084` | `atm_sync_blade_points_to_device` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4102` | `atm_sync_blade_forces_to_device` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4121` | `atm_lesgo_findCells` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:4319` | `atm_lesgo_forcing` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:4798` | `atm_lesgo_mpi_gather` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:4991` | `atm_lesgo_mpi_gather_packed` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5219` | `atm_lesgo_mpi_gather_slim_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5304` | `atm_lesgo_mpi_gather_slim_batch_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5394` | `atm_lesgo_mpi_gather_packed_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5570` | `atm_pack_blade_forces_mirror` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5600` | `atm_unpack_blade_forces_mirror` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5630` | `atm_pack_rank3` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5652` | `atm_unpack_rank3` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5674` | `atm_pack_rank4` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5701` | `atm_unpack_rank4` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5728` | `atm_pack_rank5` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5758` | `atm_unpack_rank5` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5788` | `atm_pack_gather_scalars` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5823` | `atm_lesgo_force` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:5990` | `atm_lesgo_compute_Spalart_u` | `unmarked-runtime-candidate` | `atm-mirror-lb-control` | `atm_line`, `large_windfarm` |
| `atm_lesgo_interface.f90:6072` | `atm_lesgo_convolute_force` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:6367` | `atm_convolute_atpoint_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:6518` | `atm_sample_velocity_atpoint_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:6634` | `atm_batch_atpoint_init` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:6712` | `atm_batch_sample_velocity_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:6900` | `atm_batch_convolute_force_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:7053` | `atm_batch_clc_init` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:7121` | `atm_batch_cl_correction_gpu` | `gpu-marked` |  |  |
| `atm_lesgo_interface.f90:7377` | `atm_lesgo_apply_force` | `gpu-marked` |  |  |
| `cfl_util.f90:49` | `cfl_cuda_sync` | `gpu-marked` |  |  |
| `cfl_util.f90:71` | `get_max_cfl` | `gpu-marked` |  |  |
| `cfl_util.f90:157` | `get_cfl_dt` | `gpu-marked` |  |  |
| `clocks.f90:47` | `start` | `host-boundary` |  |  |
| `clocks.f90:65` | `stop` | `host-boundary` |  |  |
| `concurrent_precursor.f90:110` | `cps_timer_start` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:124` | `cps_timer_accum` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:141` | `cps_stage_report` | `host-or-diagnostic` |  |  |
| `concurrent_precursor.f90:185` | `initialize_cps` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:240` | `synchronize_cps` | `gpu-marked` |  |  |
| `concurrent_precursor.f90:452` | `inflow_cps` | `gpu-marked` |  |  |
| `convec.f90:98` | `convec_cuda_impl` | `gpu-marked` |  |  |
| `convec.f90:331` | `ensure_convec_cuda` | `gpu-marked` |  |  |
| `convec.f90:402` | `padd_3d_dp` | `gpu-marked` |  |  |
| `convec.f90:466` | `unpadd_3d_dp` | `gpu-marked` |  |  |
| `convec.f90:539` | `check_convec_cuda` | `gpu-marked` |  |  |
| `convec.f90:562` | `convec_cuda_sync` | `gpu-marked` |  |  |
| `convec.f90:583` | `require_convec_cufft` | `gpu-marked` |  |  |
| `convec.f90:601` | `convec` | `gpu-marked` |  |  |
| `convec_gpu.f90:71` | `convec_gpu_init` | `gpu-marked` |  |  |
| `convec_gpu.f90:108` | `convec_gpu_finalize` | `gpu-marked` |  |  |
| `convec_gpu.f90:120` | `convec_gpu` | `gpu-marked` |  |  |
| `convec_gpu.f90:540` | `unpadd_gpu` | `gpu-marked` |  |  |
| `coriolis.f90:96` | `coriolis_cuda_sync` | `gpu-marked` |  |  |
| `coriolis.f90:118` | `coriolis_init` | `host-or-diagnostic` |  |  |
| `coriolis.f90:198` | `coriolis_finalize` | `host-or-diagnostic` |  |  |
| `coriolis.f90:211` | `coriolis_calc` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:40` | `cuda_mpi_debug_init` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:107` | `print_env` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:143` | `cuda_pre` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:172` | `cuda_post` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:204` | `probe_real` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:230` | `probe_complex` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:256` | `mpi_dbg_sendrecv_r` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:294` | `mpi_dbg_sendrecv_c` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:332` | `mpi_dbg_send_r` | `gpu-marked` |  |  |
| `cuda_mpi_debug.f90:363` | `mpi_dbg_recv_r` | `gpu-marked` |  |  |
| `derivatives.f90:100` | `filt_da_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:180` | `xy_derivs_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:295` | `stress_uv_xy_derivs_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:395` | `stress_w_xy_derivs_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:480` | `stress_uv_div_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:605` | `stress_w_div_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:716` | `filt_da_vel_cuda` | `gpu-marked` |  |  |
| `derivatives.f90:817` | `ensure_filt_da_cuda_plan` | `gpu-marked` |  |  |
| `derivatives.f90:885` | `ensure_filt_da_vel_cuda_plan` | `gpu-marked` |  |  |
| `derivatives.f90:953` | `ensure_stress_uv_xy_cuda_plan` | `gpu-marked` |  |  |
| `derivatives.f90:1019` | `ensure_stress_w_xy_cuda_plan` | `gpu-marked` |  |  |
| `derivatives.f90:1148` | `derivatives_cuda_sync` | `gpu-marked` |  |  |
| `derivatives.f90:1169` | `require_filt_da_cufft_success` | `gpu-marked` |  |  |
| `derivatives.f90:1185` | `require_filt_da_cuda_success` | `gpu-marked` |  |  |
| `derivatives.f90:1201` | `stress_uv_xy_derivs` | `gpu-marked` |  |  |
| `derivatives.f90:1234` | `stress_w_xy_derivs` | `gpu-marked` |  |  |
| `derivatives.f90:1262` | `ddx` | `gpu-marked` |  |  |
| `derivatives.f90:1316` | `ddy` | `gpu-marked` |  |  |
| `derivatives.f90:1370` | `ddxy` | `gpu-marked` |  |  |
| `derivatives.f90:1427` | `filt_da` | `gpu-marked` |  |  |
| `derivatives.f90:1487` | `filt_da_vel` | `gpu-marked` |  |  |
| `derivatives.f90:1529` | `ddz_vel` | `gpu-marked` |  |  |
| `derivatives.f90:1629` | `ddz_uv` | `gpu-marked` |  |  |
| `derivatives.f90:1731` | `ddz_w` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:64` | `filt_da_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:157` | `ddz_uv_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:219` | `ddz_w_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:272` | `ddx_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:343` | `ddy_gpu` | `gpu-marked` |  |  |
| `derivatives_gpu.f90:406` | `ddxy_gpu` | `gpu-marked` |  |  |
| `divstress_uv.f90:22` | `divstress_uv` | `gpu-marked` |  |  |
| `divstress_w.f90:22` | `divstress_w` | `gpu-marked` |  |  |
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
| `forcing.f90:164` | `forcing_cuda_sync` | `gpu-marked` |  |  |
| `forcing.f90:204` | `lvlset_bridge_time` | `unmarked-runtime-candidate` | `excluded-lvlset-bridge` | `lvlset` |
| `forcing.f90:222` | `lvlset_bridge_report` | `host-or-diagnostic` |  |  |
| `forcing.f90:264` | `forcing_random` | `gpu-marked` |  |  |
| `forcing.f90:368` | `forcing_applied` | `gpu-marked` |  |  |
| `forcing.f90:497` | `forcing_induced` | `gpu-marked` |  |  |
| `forcing.f90:554` | `project` | `gpu-marked` |  |  |
| `forcing.f90:943` | `project_sync_velocity_halos_cuda` | `gpu-marked` |  |  |
| `forcing.f90:984` | `project_sync_velocity_direct_halos_cuda` | `gpu-marked` |  |  |
| `forcing.f90:1065` | `project_sync_velocity_direct_halos_overlap_cuda` | `gpu-marked` |  |  |
| `forcing.f90:1178` | `project_stage_report` | `host-or-diagnostic` |  |  |
| `forcing.f90:1215` | `project_ensure_halo_buffers` | `host-or-diagnostic` |  |  |
| `forcing.f90:1236` | `project_pack_velocity_halos_cuda` | `gpu-marked` |  |  |
| `forcing.f90:1266` | `project_unpack_velocity_halos_cuda` | `gpu-marked` |  |  |
| `fringe.f90:68` | `constructor` | `unmarked-runtime-candidate` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `functions.f90:74` | `tau_wall_cuda_sync` | `gpu-marked` |  |  |
| `functions.f90:96` | `interp_to_uv_grid` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:174` | `interp_to_w_grid` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:636` | `bilinear_interp_sa_nocheck` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:672` | `bilinear_interp_sa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:699` | `bilinear_interp_aa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:757` | `linear_interp_sa_nocheck` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:783` | `linear_interp_sa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:809` | `linear_interp_aa` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:844` | `cross_product` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:856` | `binary_search` | `unmarked-runtime-candidate` | `generic-helper-profile` | `les_core_channel`, `adm_disk`, `atm_line` |
| `functions.f90:927` | `get_tau_wall_bot` | `gpu-marked` |  |  |
| `functions.f90:983` | `get_tau_wall_top` | `gpu-marked` |  |  |
| `functions.f90:1038` | `count_lines` | `host-or-diagnostic` |  |  |
| `grid.f90:44` | `build` | `host-boundary` |  |  |
| `hit_inflow.f90:82` | `initialize_HIT` | `gpu-marked` |  |  |
| `hit_inflow.f90:153` | `extract_HIT_data` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `hit_inflow.f90:196` | `compute_HIT_plane_data` | `gpu-marked` |  |  |
| `hit_inflow.f90:242` | `inflow_HIT` | `gpu-marked` |  |  |
| `hit_inflow.f90:273` | `hit_write_restart` | `host-or-diagnostic` |  |  |
| `hit_inflow.f90:293` | `hit_read_restart` | `host-or-diagnostic` |  |  |
| `hit_inflow.f90:320` | `interpolate3D` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `hit_inflow_gpu.f90:75` | `hit_gpu_setup` | `gpu-marked` |  |  |
| `hit_inflow_gpu.f90:123` | `hit_fringe_setup_gpu` | `gpu-marked` |  |  |
| `hit_inflow_gpu.f90:150` | `hit_compute_plane_gpu` | `gpu-marked` |  |  |
| `hit_inflow_gpu.f90:322` | `hit_apply_fringe_gpu` | `gpu-marked` |  |  |
| `inflow.f90:70` | `inflow_cuda_sync` | `gpu-marked` |  |  |
| `inflow.f90:93` | `inflow_init` | `gpu-marked` |  |  |
| `inflow.f90:137` | `apply_inflow` | `unmarked-runtime-candidate` | `inflow-fringe-profile` | `hit_inflow`, `shifted_inflow`, `sponge_coriolis` |
| `inflow.f90:163` | `inflow_uniform` | `gpu-marked` |  |  |
| `init_random_seed.f90:21` | `init_random_seed` | `host-boundary` |  |  |
| `initial.f90:21` | `initial` | `gpu-marked` |  |  |
| `initial.f90:206` | `initial_cuda_sync` | `gpu-marked` |  |  |
| `initial.f90:227` | `initial_cuda_touch_velocity` | `gpu-marked` |  |  |
| `initial.f90:256` | `ic_uniform` | `gpu-marked` |  |  |
| `initial.f90:292` | `check_for_interp` | `host-boundary` |  |  |
| `initial.f90:317` | `ic_file` | `host-boundary` |  |  |
| `initial.f90:333` | `ic_interp` | `host-boundary` |  |  |
| `initial.f90:474` | `ic_dns` | `host-boundary` |  |  |
| `initial.f90:564` | `ic_les` | `host-boundary` |  |  |
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
| `interpolag_Sdep.f90:21` | `interpolag_Sdep` | `gpu-marked` |  |  |
| `interpolag_Sdep.f90:475` | `interpolag_sdep_cuda_sync` | `gpu-marked` |  |  |
| `interpolag_Ssim.f90:21` | `interpolag_Ssim` | `gpu-marked` |  |  |
| `interpolag_Ssim.f90:309` | `interpolag_ssim_cuda_sync` | `gpu-marked` |  |  |
| `io.f90:73` | `openfiles` | `host-boundary` |  |  |
| `io.f90:112` | `energy` | `gpu-marked` |  |  |
| `io.f90:216` | `write_tau_wall_bot` | `host-boundary` |  |  |
| `io.f90:243` | `write_tau_wall_top` | `host-boundary` |  |  |
| `io.f90:272` | `write_parallel_cgns` | `host-boundary` |  |  |
| `io.f90:441` | `write_null_cgns` | `host-boundary` |  |  |
| `io.f90:606` | `output_loop` | `host-boundary` |  |  |
| `io.f90:753` | `inst_write` | `gpu-marked` |  |  |
| `io.f90:1577` | `checkpoint` | `gpu-marked` |  |  |
| `io.f90:1683` | `output_final` | `host-boundary` |  |  |
| `io.f90:1697` | `output_init` | `host-boundary` |  |  |
| `iwmles.f90:156` | `iwm_cuda_sync` | `gpu-marked` |  |  |
| `iwmles.f90:192` | `iwm_wallstress` | `gpu-marked` |  |  |
| `iwmles.f90:268` | `iwm_init` | `gpu-marked` |  |  |
| `iwmles.f90:391` | `iwm_finalize` | `gpu-marked` |  |  |
| `iwmles.f90:446` | `iwm_calc_lhs` | `gpu-marked` |  |  |
| `iwmles.f90:830` | `iwm_slv` | `unmarked-runtime-candidate` | `iwm-wallmodel-profile` | `iwm_wall_model` |
| `iwmles.f90:855` | `iwm_calc_wallstress` | `gpu-marked` |  |  |
| `iwmles.f90:1599` | `iwm_monitor` | `gpu-marked` |  |  |
| `iwmles.f90:1641` | `iwm_checkPoint` | `gpu-marked` |  |  |
| `iwmles.f90:1675` | `iwm_read_checkPoint` | `gpu-marked` |  |  |
| `lagrange_Sdep.f90:21` | `lagrange_Sdep` | `gpu-marked` |  |  |
| `lagrange_Sdep.f90:812` | `lagrange_sdep_cuda_sync` | `gpu-marked` |  |  |
| `lagrange_Sdep.f90:844` | `lagrange_sdep_cuda_barrier` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:135` | `lagrange_Sdep_gpu_init` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:215` | `lagrange_Ssim_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:488` | `lagrange_Sdep_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:1084` | `interpolag_Ssim_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:1251` | `interpolag_Sdep_gpu` | `gpu-marked` |  |  |
| `lagrange_Sdep_gpu.f90:1898` | `sync_downup_F` | `gpu-marked` |  |  |
| `lagrange_Ssim.f90:21` | `lagrange_Ssim` | `gpu-marked` |  |  |
| `lagrange_Ssim.f90:542` | `lagrange_ssim_cuda_sync` | `gpu-marked` |  |  |
| `linear_simple.f90:46` | `solve_linear` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:86` | `assert_eq2` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:106` | `assert_eq3` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:126` | `assert_eq4` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:154` | `ludcmp` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:226` | `lubksb` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:278` | `outerprod` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `linear_simple.f90:288` | `swap` | `unmarked-runtime-candidate` | `atm-host-model` | `atm_line`, `large_windfarm` |
| `main.f90:1257` | `main_read_env_real` | `host-or-diagnostic` |  |  |
| `main.f90:1293` | `main_cuda_sync` | `gpu-marked` |  |  |
| `messages.f90:66` | `message_a` | `host-boundary` |  |  |
| `messages.f90:76` | `message_ai` | `host-boundary` |  |  |
| `messages.f90:87` | `message_aiai` | `host-boundary` |  |  |
| `messages.f90:98` | `message_aiar` | `host-boundary` |  |  |
| `messages.f90:111` | `message_al` | `host-boundary` |  |  |
| `messages.f90:122` | `message_aii` | `host-boundary` |  |  |
| `messages.f90:134` | `message_air` | `host-boundary` |  |  |
| `messages.f90:147` | `message_ai_array` | `host-boundary` |  |  |
| `messages.f90:161` | `message_aiai_array` | `host-boundary` |  |  |
| `messages.f90:176` | `message_ar` | `host-boundary` |  |  |
| `messages.f90:187` | `message_ar_array` | `host-boundary` |  |  |
| `messages.f90:200` | `message_aiar_array` | `host-boundary` |  |  |
| `messages.f90:215` | `warn` | `host-boundary` |  |  |
| `messages.f90:228` | `error_a` | `host-boundary` |  |  |
| `messages.f90:244` | `error_ai` | `host-boundary` |  |  |
| `messages.f90:261` | `error_ai_array` | `host-boundary` |  |  |
| `messages.f90:282` | `error_aia` | `host-boundary` |  |  |
| `messages.f90:299` | `error_aiai` | `host-boundary` |  |  |
| `messages.f90:316` | `error_aiar` | `host-boundary` |  |  |
| `messages.f90:334` | `error_arar` | `host-boundary` |  |  |
| `messages.f90:351` | `error_al` | `host-boundary` |  |  |
| `messages.f90:368` | `error_ar` | `host-boundary` |  |  |
| `messages.f90:385` | `error_ar_array` | `host-boundary` |  |  |
| `mpi_defs.f90:60` | `initialize_mpi` | `gpu-marked` |  |  |
| `mpi_defs.f90:151` | `bind_cuda_device` | `gpu-marked` |  |  |
| `mpi_defs.f90:258` | `create_mpi_comms_cps` | `host-boundary` |  |  |
| `mpi_defs.f90:308` | `mpi_sync_real_array` | `gpu-marked` |  |  |
| `mpi_defs.f90:413` | `sync_up` | `host-boundary` |  |  |
| `mpi_defs.f90:424` | `sync_downup_nb` | `host-boundary` |  |  |
| `mpi_transpose_mod.f90:36` | `mpi_transpose` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `param_output.f90:21` | `param_output` | `host-boundary` |  |  |
| `pid.f90:49` | `constructor` | `host-boundary` |  |  |
| `pid.f90:77` | `advance_noset` | `host-boundary` |  |  |
| `pid.f90:98` | `advance_set` | `host-boundary` |  |  |
| `press_gpu.f90:82` | `press_gpu_init` | `gpu-marked` |  |  |
| `press_gpu.f90:105` | `press_gpu_finalize` | `gpu-marked` |  |  |
| `press_gpu.f90:118` | `press_stag_array_gpu` | `gpu-marked` |  |  |
| `press_stag_array.f90:44` | `press_stag_array` | `gpu-marked` |  |  |
| `press_stag_array.f90:1760` | `press_apply_env_enabled_unless_false` | `host-or-diagnostic` |  |  |
| `press_stag_array.f90:1789` | `press_cuda_sync` | `gpu-marked` |  |  |
| `press_stag_array.f90:1828` | `press_queue_event_start` | `gpu-marked` |  |  |
| `press_stag_array.f90:1846` | `press_queue_event_stop` | `gpu-marked` |  |  |
| `press_stag_array.f90:1885` | `press_queue_report` | `host-or-diagnostic` |  |  |
| `press_stag_array.f90:1935` | `press_rhs_halo_audit` | `gpu-marked` |  |  |
| `press_stag_array.f90:1976` | `press_rhs_halo_report` | `host-or-diagnostic` |  |  |
| `press_stag_array.f90:2029` | `press_rhs_assembly_report` | `host-or-diagnostic` |  |  |
| `press_stag_array.f90:2070` | `press_stage_report` | `host-or-diagnostic` |  |  |
| `press_stag_array.f90:2115` | `press_pack_rhs_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2145` | `press_pack_rhs_halo_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2172` | `press_unpack_rhs_halo_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2199` | `press_pack_rhs_halo_combined_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2238` | `press_unpack_rhs_halo_combined_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2275` | `press_rhs_prep_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2321` | `press_assemble_rhs_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2341` | `press_assemble_rhs_range_cuda` | `gpu-marked` |  |  |
| `press_stag_array.f90:2379` | `press_zero_mode_cuda` | `gpu-marked` |  |  |
| `rmsdiv.f90:21` | `rmsdiv` | `gpu-marked` |  |  |
| `scalars.f90:199` | `scalars_acc_sync` | `gpu-marked` |  |  |
| `scalars.f90:226` | `scalars_timer_start` | `host-or-diagnostic` |  |  |
| `scalars.f90:238` | `scalars_timer_accum` | `host-or-diagnostic` |  |  |
| `scalars.f90:255` | `scalars_stage_report` | `host-or-diagnostic` |  |  |
| `scalars.f90:284` | `scalars_deriv_xy_big_acc` | `gpu-marked` |  |  |
| `scalars.f90:415` | `scalars_to_big_acc` | `gpu-marked` |  |  |
| `scalars.f90:464` | `scalars_return_rhs_acc` | `gpu-marked` |  |  |
| `scalars.f90:504` | `scalars_divergence_acc` | `gpu-marked` |  |  |
| `scalars.f90:564` | `scalars_update_device_state` | `gpu-marked` |  |  |
| `scalars.f90:576` | `scalars_copy_rhs_acc` | `gpu-marked` |  |  |
| `scalars.f90:595` | `scalars_advective_acc` | `gpu-marked` |  |  |
| `scalars.f90:645` | `scalars_flux_acc` | `gpu-marked` |  |  |
| `scalars.f90:729` | `scalars_rhs_theta_acc` | `gpu-marked` |  |  |
| `scalars.f90:790` | `scalars_cuda_sync` | `gpu-marked` |  |  |
| `scalars.f90:811` | `require_scalars_cufft` | `gpu-marked` |  |  |
| `scalars.f90:826` | `ensure_scalars_cufft` | `gpu-marked` |  |  |
| `scalars.f90:872` | `scalars_padd_3d_gpu` | `gpu-marked` |  |  |
| `scalars.f90:906` | `scalars_unpadd_3d_gpu` | `gpu-marked` |  |  |
| `scalars.f90:942` | `scalars_to_big_gpu` | `gpu-marked` |  |  |
| `scalars.f90:976` | `scalars_return_rhs_gpu` | `gpu-marked` |  |  |
| `scalars.f90:1013` | `scalars_timer_start` | `gpu-marked` |  |  |
| `scalars.f90:1025` | `scalars_timer_accum` | `gpu-marked` |  |  |
| `scalars.f90:1042` | `scalars_stage_report` | `host-or-diagnostic` |  |  |
| `scalars.f90:1133` | `stability_acc` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1282` | `scalars_copy_rhs_gpu` | `gpu-marked` |  |  |
| `scalars.f90:1301` | `scalars_advective_gpu` | `gpu-marked` |  |  |
| `scalars.f90:1347` | `scalars_flux_gpu` | `gpu-marked` |  |  |
| `scalars.f90:1430` | `scalars_rhs_theta_gpu` | `gpu-marked` |  |  |
| `scalars.f90:1480` | `buoyancy_force_gpu` | `gpu-marked` |  |  |
| `scalars.f90:1511` | `scalars_init` | `gpu-marked` |  |  |
| `scalars.f90:1595` | `ic_scal` | `gpu-marked` |  |  |
| `scalars.f90:1634` | `ic_scal_file` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1653` | `ic_scal_les` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1672` | `ic_scal_interp` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:1793` | `scalars_checkpoint` | `gpu-marked` |  |  |
| `scalars.f90:1813` | `scalars_deriv` | `gpu-marked` |  |  |
| `scalars.f90:1892` | `obukhov` | `gpu-marked` |  |  |
| `scalars.f90:2204` | `scalars_transport` | `gpu-marked` |  |  |
| `scalars.f90:2624` | `to_big` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `scalars.f90:2648` | `buoyancy_force` | `gpu-marked` |  |  |
| `scaledep_dynamic.f90:21` | `scaledep_dynamic` | `gpu-marked` |  |  |
| `scaledep_dynamic.f90:554` | `scaledep_dynamic_cuda_sync` | `gpu-marked` |  |  |
| `sgs_gpu.f90:109` | `divstress_gpu_init` | `gpu-marked` |  |  |
| `sgs_gpu.f90:126` | `std_dynamic_pples_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:278` | `scaledep_dynamic_pples_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:554` | `sgs_stag_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:915` | `calc_Sij_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:1082` | `divstress_uv_gpu` | `gpu-marked` |  |  |
| `sgs_gpu.f90:1185` | `divstress_w_gpu` | `gpu-marked` |  |  |
| `sgs_param.f90:156` | `sgs_param_init` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:458` | `sgs_cuda_sync` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:481` | `sgs_cuda_barrier` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:502` | `sgs_event_record` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:519` | `sgs_event_elapsed_seconds` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:546` | `sgs_diag_time` | `host-or-diagnostic` |  |  |
| `sgs_stag_util.f90:564` | `sgs_calc_diag_begin` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:610` | `sgs_calc_diag_start` | `host-or-diagnostic` |  |  |
| `sgs_stag_util.f90:623` | `sgs_calc_diag_stop` | `host-or-diagnostic` |  |  |
| `sgs_stag_util.f90:644` | `sgs_calc_cpu_start` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:655` | `sgs_calc_cpu_stop` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:668` | `sgs_calc_set_zrange` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:680` | `sgs_calc_set_audit` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:710` | `sgs_tau_detail_begin` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:731` | `sgs_tau_detail_start` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:742` | `sgs_tau_detail_stop` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:755` | `sgs_tau_detail_add_bytes` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:767` | `sgs_tau_detail_add_msg` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:782` | `sgs_dwdz_detail_begin` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:804` | `sgs_dwdz_detail_start` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:815` | `sgs_dwdz_detail_stop` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:828` | `sgs_dwdz_detail_add_msg` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:846` | `sgs_dwdz_path_audit` | `unmarked-runtime-candidate` | `diagnostic-profiling` | `diagnostics_output` |
| `sgs_stag_util.f90:892` | `sgs_pointer_env_audit` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:941` | `sgs_pointer_env_audit_device` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:991` | `sgs_stag` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:2385` | `sgs_calc_sij_device_lower_bound_bench` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:2587` | `calc_Sij` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:3303` | `calc_Sij_nut_dynamic_cuda` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:3392` | `sgs_sync_dwdz_down_cuda` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:3431` | `sgs_calc_sij_detail_report` | `gpu-marked` |  |  |
| `sgs_stag_util.f90:3711` | `sgs_tau_halo_detail_report` | `host-or-diagnostic` |  |  |
| `sgs_stag_util.f90:3853` | `sgs_stage_report` | `gpu-marked` |  |  |
| `shifted_inflow.f90:80` | `shifted_inflow_cuda_sync` | `gpu-marked` |  |  |
| `shifted_inflow.f90:103` | `shifted_inflow_init` | `gpu-marked` |  |  |
| `shifted_inflow.f90:146` | `inflow_shifted` | `gpu-marked` |  |  |
| `sim_param.f90:86` | `sim_param_init` | `gpu-marked` |  |  |
| `sponge.f90:59` | `sponge_cuda_sync` | `gpu-marked` |  |  |
| `sponge.f90:81` | `sponge_init` | `gpu-marked` |  |  |
| `sponge.f90:105` | `sponge_force` | `gpu-marked` |  |  |
| `stability.f90:2` | `stability` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:60` | `calc_psi_m` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `stability.f90:96` | `calc_psi_h` | `unmarked-runtime-candidate` | `scalar-init-fallback` | `scalar_passive`, `scalar_active`, `cps_scalar` |
| `std_dynamic.f90:21` | `std_dynamic` | `gpu-marked` |  |  |
| `std_dynamic.f90:306` | `std_dynamic_cuda_sync` | `gpu-marked` |  |  |
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
| `test_filtermodule.f90:94` | `require_test_filter_cufft_success` | `gpu-marked` |  |  |
| `test_filtermodule.f90:110` | `test_filter_cuda_sync` | `gpu-marked` |  |  |
| `test_filtermodule.f90:133` | `test_filter_cuda_barrier` | `gpu-marked` |  |  |
| `test_filtermodule.f90:154` | `ensure_test_filter_cuda_plan` | `gpu-marked` |  |  |
| `test_filtermodule.f90:185` | `ensure_test_filter_cuda_many_plan` | `gpu-marked` |  |  |
| `test_filtermodule.f90:217` | `ensure_test_filter_cuda_12_plan` | `gpu-marked` |  |  |
| `test_filtermodule.f90:248` | `apply_test_filter_cuda` | `gpu-marked` |  |  |
| `test_filtermodule.f90:283` | `apply_test_filter_cuda_3` | `gpu-marked` |  |  |
| `test_filtermodule.f90:329` | `apply_test_filter_cuda_6` | `gpu-marked` |  |  |
| `test_filtermodule.f90:381` | `apply_test_filter_cuda_managed` | `gpu-marked` |  |  |
| `test_filtermodule.f90:429` | `apply_test_filter_cuda_3_managed` | `gpu-marked` |  |  |
| `test_filtermodule.f90:485` | `apply_test_filter_cuda_6_managed` | `gpu-marked` |  |  |
| `test_filtermodule.f90:547` | `apply_test_filter_cuda_3_dual_managed` | `gpu-marked` |  |  |
| `test_filtermodule.f90:618` | `apply_test_filter_cuda_6_dual_managed` | `gpu-marked` |  |  |
| `test_filtermodule.f90:704` | `test_filter_init` | `gpu-marked` |  |  |
| `test_filtermodule.f90:806` | `test_filter` | `gpu-marked` |  |  |
| `test_filtermodule.f90:836` | `test_filter_3` | `gpu-marked` |  |  |
| `test_filtermodule.f90:857` | `test_filter_6` | `gpu-marked` |  |  |
| `test_filtermodule.f90:882` | `test_test_filter` | `gpu-marked` |  |  |
| `test_filtermodule.f90:912` | `test_test_filter_3` | `gpu-marked` |  |  |
| `test_filtermodule.f90:933` | `test_test_filter_6` | `gpu-marked` |  |  |
| `test_filtermodule.f90:959` | `test_filter_b_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:994` | `test_test_filter_b_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1028` | `test_filter_plane_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1055` | `test_test_filter_plane_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1083` | `test_filter_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1094` | `test_filter_3_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1105` | `test_filter_6_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1117` | `test_test_filter_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1128` | `test_test_filter_3_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1140` | `test_test_filter_6_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1152` | `test_filter_3_dual_gpu` | `gpu-marked` |  |  |
| `test_filtermodule.f90:1165` | `test_filter_6_dual_gpu` | `gpu-marked` |  |  |
| `time_average.f90:74` | `init` | `gpu-marked` |  |  |
| `time_average.f90:201` | `compute` | `gpu-marked` |  |  |
| `time_average.f90:462` | `finalize` | `gpu-marked` |  |  |
| `time_average.f90:815` | `checkpoint` | `gpu-marked` |  |  |
| `time_average.f90:881` | `write_parallel_cgns` | `host-or-diagnostic` |  |  |
| `time_average.f90:1051` | `write_null_cgns` | `host-or-diagnostic` |  |  |
| `tridag_array.f90:252` | `tridag_array` | `gpu-marked` |  |  |
| `tridag_array.f90:744` | `tridag_apply_env_enabled_unless_false` | `host-or-diagnostic` |  |  |
| `tridag_array.f90:772` | `tridag_apply_env_true_token` | `host-or-diagnostic` |  |  |
| `tridag_array.f90:799` | `tridag_array_transpose_thomas_cuda` | `gpu-marked` |  |  |
| `tridag_array.f90:1680` | `tridag_array_spike2_cuda` | `gpu-marked` |  |  |
| `tridag_array.f90:1889` | `tridag_array_replicated_cuda` | `gpu-marked` |  |  |
| `tridag_array.f90:2080` | `tridag_array` | `unmarked-runtime-candidate` | `cpu-fallback-compat` | `les_core_channel`, `hit_inflow` |
| `tridag_gpu.f90:107` | `tridag_gpu_init` | `gpu-marked` |  |  |
| `tridag_gpu.f90:126` | `tridag_gpu_finalize` | `gpu-marked` |  |  |
| `tridag_gpu.f90:136` | `tridag_array_gpu` | `gpu-marked` |  |  |
| `turbine_indicator.f90:56` | `val` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbine_indicator.f90:80` | `init` | `gpu-marked` |  |  |
| `turbines.f90:160` | `turbines_init` | `gpu-marked` |  |  |
| `turbines.f90:316` | `turbines_nodes` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:516` | `turbines_acc_metadata_init` | `gpu-marked` |  |  |
| `turbines.f90:662` | `turbines_acc_finalize` | `host-or-diagnostic` |  |  |
| `turbines.f90:703` | `turbines_acc_sync_device_field` | `gpu-marked` |  |  |
| `turbines.f90:736` | `turbines_forcing_acc` | `gpu-marked` |  |  |
| `turbines.f90:989` | `turbines_forcing_acc` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:999` | `turbines_forcing` | `gpu-marked` |  |  |
| `turbines.f90:1194` | `turbines_finalize` | `gpu-marked` |  |  |
| `turbines.f90:1214` | `turbines_checkpoint` | `host-or-diagnostic` |  |  |
| `turbines.f90:1239` | `turbine_vel_init` | `host-or-diagnostic` |  |  |
| `turbines.f90:1279` | `place_turbines` | `unmarked-runtime-candidate` | `adm-cpu-fallback-profile` | `adm_disk`, `adm_dynamic_controls` |
| `turbines.f90:1429` | `read_control_files` | `host-or-diagnostic` |  |  |
| `turbines_gpu.f90:42` | `turbines_interp_w_to_uv_gpu` | `gpu-marked` |  |  |
| `wallstress.f90:21` | `wallstress` | `gpu-marked` |  |  |
| `wallstress.f90:120` | `wallstress_cuda_sync` | `gpu-marked` |  |  |
| `wallstress.f90:142` | `ws_free_lbc` | `gpu-marked` |  |  |
| `wallstress.f90:188` | `ws_free_ubc` | `gpu-marked` |  |  |
| `wallstress.f90:232` | `ws_dns_lbc` | `gpu-marked` |  |  |
| `wallstress.f90:282` | `ws_dns_ubc` | `gpu-marked` |  |  |
| `wallstress.f90:332` | `ws_equilibrium_lbc` | `gpu-marked` |  |  |
| `wallstress.f90:506` | `ws_equilibrium_ubc` | `gpu-marked` |  |  |

Regenerate this file with:

```bash
python3 tools/report_gpu_static_full_inventory.py --write
```
