# lesgo.conf Key Validation Coverage

This generated table expands the non-LVLSET `input_util.f90` parser
keys into their responsible GPU validation rows.  A mapping means the
row is responsible for correctness and timing evidence; it is not, by
itself, a speedup claim.

| Parser group | Key | Validation rows | Current evidence states |
| --- | --- | --- | --- |
| `DOMAIN` | `LX` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `LY` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `LZ` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `NPROC` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `NX` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `NY` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `NZ` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `UNIFORM_SPACING` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `DOMAIN` | `Z_I` | `les_core_channel`, `large_windfarm` | `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `CO` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `CS_COUNT` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `DYN_INIT` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `IFILTER` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `MOLEC` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `NU_MOLEC` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `SGS` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `SGS_MODEL` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `U_STAR` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `VONK` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `MODEL` | `WALL_DAMP_EXP` | `sgs_disabled`, `sgs_models_1_5`, `dyn_tn`, `iwm_wall_model` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed` |
| `CORIOLIS` | `ALPHA` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `CORIOLIS_FORCING` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `FC` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `G` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `HEIGHT_SET` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `KD` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `KI` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `KP` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `PHI_SET` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `PID_TIME` | `sponge_coriolis` | `needs_benchmark` |
| `CORIOLIS` | `REPEAT_INTERVAL` | `sponge_coriolis` | `needs_benchmark` |
| `FLOW_COND` | `EVAL_MEAN_P_FORCE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `FRINGE_REGION_END` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `FRINGE_REGION_LEN` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `INFLOW_TYPE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `INFLOW_VELOCITY` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `INILAG` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `INITU` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `LBC_MOM` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `MEAN_P_FORCE_X` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `MEAN_P_FORCE_Y` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `RMS_RANDOM_FORCE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `SAMPLING_REGION_END` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `SHIFT_N` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `SPONGE_FREQUENCY` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `SPONGE_HEIGHT` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `STOP_RANDOM_FORCE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `UBC_MOM` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `UBOT` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `USE_MEAN_P_FORCE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `USE_RANDOM_FORCE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `USE_SPONGE` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `UTOP` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `FLOW_COND` | `ZO` | `les_core_channel`, `iwm_wall_model`, `shifted_inflow`, `sponge_coriolis` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `OUTPUT` | `CHECKPOINT_DATA` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `CHECKPOINT_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `DOMAIN_CALC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `DOMAIN_NEND` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `DOMAIN_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `DOMAIN_NSTART` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `LAG_CFL_COUNT` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `NENERGY` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `POINT_CALC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `POINT_LOC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `POINT_NEND` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `POINT_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `POINT_NSTART` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `TAVG_CALC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `TAVG_NEND` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `TAVG_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `TAVG_NSTART` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `WBASE` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `XPLANE_CALC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `XPLANE_LOC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `XPLANE_NEND` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `XPLANE_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `XPLANE_NSTART` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `YPLANE_CALC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `YPLANE_LOC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `YPLANE_NEND` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `YPLANE_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `YPLANE_NSTART` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `ZPLANE_CALC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `ZPLANE_LOC` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `ZPLANE_NEND` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `ZPLANE_NSKIP` | `diagnostics_output` | `host_boundary_verified` |
| `OUTPUT` | `ZPLANE_NSTART` | `diagnostics_output` | `host_boundary_verified` |
| `TURBINES` | `ADM_CORRECTION` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `ALPHA1` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `ALPHA2` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `CT_PRIME` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `DIA_ALL` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `DYN_CT_PRIME` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `DYN_THETA1` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `DYN_THETA2` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `FILTER_CUTOFF` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `HEIGHT_ALL` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `NUM_X` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `NUM_Y` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `ORIENTATION` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `READ_PARAM` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `STAG_PERC` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `TBASE` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `THETA1_ALL` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `THETA2_ALL` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `THK_ALL` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `TIP_SPEED_RATIO` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `T_AVG_DIM` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `TURBINES` | `USE_ROTATION` | `adm_disk`, `atm_line`, `large_windfarm`, `adm_dynamic_controls` | `paired_speedup_claimed`, `paired_speedup_claimed`, `paired_speedup_claimed`, `needs_benchmark` |
| `SCALARS` | `FLUX_BOT` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `G` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `IC_NO_VEL_NOISE_Z` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `IC_THETA` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `IC_Z` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `LAPSE_RATE` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `LBC_SCAL` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `PASSIVE_SCALAR` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `PR_SGS` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `READ_LBC_SCAL` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `SCAL_BOT` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `T_SCALE` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `SCALARS` | `ZO_S` | `scalar_passive`, `scalar_active`, `cps_scalar` | `needs_benchmark`, `needs_benchmark`, `paired_speedup_claimed` |
| `HIT` | `LX_HIT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `LY_HIT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `LZ_HIT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `NX_HIT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `NY_HIT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `NZ_HIT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `TI_OUT` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `UP_IN` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `U_FILE` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `V_FILE` | `hit_inflow` | `external_record_needs_copy` |
| `HIT` | `W_FILE` | `hit_inflow` | `external_record_needs_copy` |

Regenerate this file with:

```bash
python3 tools/report_lesgo_conf_key_validation.py --write
```
