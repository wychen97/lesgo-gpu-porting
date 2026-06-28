!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by:
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   The Actuator Turbine Model is free software: you can redistribute it
!!   and/or modify it under the terms of the GNU General Public License as
!!   published by the Free Software Foundation, either version 3 of the
!!   License, or (at your option) any later version.
!!
!!   The Actuator Turbine Model is distributed in the hope that it will be
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************************
module atm_lesgo_interface
!*******************************************************************************
! Navigation map for this large interface module:
!   - declarations and policy flags: module header before `contains`
!   - diagnostics/timing helpers: atm_diag_*, atm_lesgo_report_timing
!   - GPU utility kernels: atm_interp_w_to_uv, atm_prepare_direct_w,
!     atm_lesgo_apply_force_gpu, and the *_atpoint GPU routines
!   - point-owner and load-balance path: atm_point_owner_* and atm_lb_*
!   - lifecycle and planning: atm_lesgo_initialize/finalize, diag_load, lb_plan
!   - turbine geometry residency: force shadows, blade mirrors, findCells
!   - timestep entry point: atm_lesgo_forcing(phase)
!   - legacy/CPU gather and force paths: atm_lesgo_mpi_gather*,
!     atm_lesgo_force, atm_lesgo_convolute_force, atm_lesgo_apply_force
!
! Keep CPU and GPU versions semantically paired.  When changing a GPU routine,
! update the matching CPU/legacy path or document why the path intentionally
! differs in docs/gpu_module_contracts.md.
#ifdef ENABLE_CUDA
use cudafor
#endif
! This module interfaces actuator turbine module with lesgo
! It is a lesgo specific module, unlike the atm module
! The MPI management is done only in this section of the code
! This is very code dependent and will have to be modified according to
! the code being used. In this case LESGO has its own MPI details
! Look into mpi_defs.f90 for the details

! Remember to always dimensionalize the variables from LESGO
! Length is non-dimensionalized by z_i

! Lesgo data used regarding the grid (LESGO)
use param, only : dt ,nx,ny,nz,nz_tot,dx,dy,dz,coord,nproc, z_i, u_star, lbz,  &
                  total_time, jt_total, L_x, L_y, L_z
! nx, ny, nz - nodes in every direction
! z_i - non-dimensionalizing length
! dt - time-step

! These are the forces, and velocities on x,y, and z respectively
use sim_param, only : fxa, fya, fza, u, v, w

! Grid definition (LESGO)
use grid_m, only : grid

! MPI implementation from LESGO
#ifdef PPMPI
  use mpi_defs
  use mpi
  use param, only : ierr, mpi_rprec, comm, coord, up, down
#ifdef ENABLE_CUDA
  use cuda_mpi_debug, only : mpi_dbg_sendrecv_r
#endif
#endif

! Interpolating function for interpolating the velocity field to each
! actuator point
use functions, only : trilinear_interp, interp_to_uv_grid

use clock_m, only : clock_t

! Actuator Turbine Model module
use atm_base
use actuator_turbine_model
use atm_input_util, only : rprec, turbineArray, turbineModel, eat_whitespace, &
                           atm_print_initialize, updateInterval, outputInterval, &
                           diagnosticOutputInterval, airfoilType_t

! Used for testing time
! use clock_m

implicit none

! Variable for interpolating the velocity in w onto the uv grid
#ifdef PPLES_GPU
! PPLES GPU ownership map:
!   - LES fields u/v/w and fxa/fya/fza are owned by sim_param and are expected
!     to be device-resident during timestep forcing.
!   - w_uv and atm_nacelle_tmp are local timestep scratch owned by this module.
!   - shadow tables and blade mirrors are built from host turbine geometry and
!     kept on device for sampling/convolution.
!   - point-owner load-balance tables are experimental and selected only by
!     the documented LESGO_ATM_* environment controls.
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), device, allocatable, save, dimension(:) :: atm_nacelle_tmp
real(rprec), allocatable, save, dimension(:,:) :: atm_wuv_send_down
real(rprec), allocatable, save, dimension(:,:) :: atm_wuv_recv_up
!$acc declare create(w_uv, atm_wuv_send_down, atm_wuv_recv_up)
integer, allocatable, save :: atm_shadow_offsetUV(:), atm_shadow_offsetW(:)
integer, allocatable, save :: atm_shadow_ijkUV(:,:), atm_shadow_ijkW(:,:)
real(rprec), allocatable, save :: atm_shadow_locUV(:,:), atm_shadow_locW(:,:)
real(rprec), allocatable, save :: atm_shadow_forceUV(:,:), atm_shadow_forceW(:,:)
integer, device, allocatable, save :: atm_shadow_offsetUV_d(:), atm_shadow_offsetW_d(:)
integer, device, allocatable, save :: atm_shadow_ijkUV_d(:,:), atm_shadow_ijkW_d(:,:)
real(rprec), device, allocatable, save :: atm_shadow_locUV_d(:,:), atm_shadow_locW_d(:,:)
real(rprec), device, allocatable, save :: atm_shadow_forceUV_d(:,:), atm_shadow_forceW_d(:,:)
logical, save :: atm_shadow_ready = .false.
real(rprec), device, allocatable, save :: atm_bladePoints_d(:,:,:,:,:)
real(rprec), device, allocatable, save :: atm_bladeForces_d(:,:,:,:,:)
integer, save :: atm_max_m = 0, atm_max_n = 0, atm_max_q = 0
logical, save :: atm_blade_mirror_ready = .false.
#elif defined(ENABLE_CUDA)
real(rprec), managed, allocatable, dimension(:,:,:) :: w_uv
real(rprec), managed, allocatable, save, dimension(:) :: atm_nacelle_tmp
real(rprec), device, allocatable, save, dimension(:,:) :: atm_wuv_send_down
real(rprec), device, allocatable, save, dimension(:,:) :: atm_wuv_recv_up
#else
real(rprec), allocatable, dimension(:,:,:) :: w_uv
#endif

private
public atm_lesgo_initialize, atm_lesgo_forcing, atm_lesgo_finalize

! This is a list that stores all the points in the domain with a body
! force due to the turbines.
type bodyForce_t
    integer :: c ! Number of cells
    ! i,j,k stores the index for the point in the domain
#ifdef ENABLE_CUDA
    integer, managed, allocatable :: ijk(:,:)
    real(rprec), managed, allocatable :: force(:,:) ! Force vector on uv grid
    real(rprec), managed, allocatable :: location(:,:) ! Position vector on uv grid
#else
    integer, allocatable :: ijk(:,:)
    real(rprec), allocatable :: force(:,:) ! Force vector on uv grid
    real(rprec), allocatable :: location(:,:) ! Position vector on uv grid
#endif
end type bodyForce_t

#ifdef PPLES_GPU
! ---- Batched atPoint GPU state (round 3) ----
! Static concatenated tables for ALL turbines (built once): force-field cell
! locations/ijk/turbine-id, grid axes + autowrap, prefix offsets. Per force
! step only the flattened blade points/forces (~260 KB) and the per-turbine
! constants are re-uploaded; sampling and convolution then run as ONE kernel
! over all turbines instead of 60 data-region/kernel/sync rounds.
integer, parameter :: ATM_NTC = 14   ! per-turbine constant slots
logical :: atm_batch_ready   = .false.
logical :: atm_batch_sampled = .false.
integer :: atm_nbp_tot = 0, atm_cUV_tot = 0, atm_cW_tot = 0
integer,     allocatable :: atm_bp_off(:)                  ! (nTurb+1) blade-point prefix
real(rprec), allocatable :: atm_bp_all(:,:), atm_bf_all(:,:)   ! (3, nbp_tot)
real(rprec), allocatable :: atm_velbp_all(:,:)             ! (3, nbp_tot) sampled velocity
integer,     allocatable :: atm_inr_all(:)                 ! (nbp_tot) in-domain flag
real(rprec), allocatable :: atm_tconst(:,:)                ! (ATM_NTC, nTurb)
real(rprec), allocatable :: atm_locUV_all(:,:), atm_locW_all(:,:)  ! (3, c_tot) static
integer,     allocatable :: atm_ijkUV_all(:,:), atm_ijkW_all(:,:)  ! (3, c_tot) static
integer,     allocatable :: atm_tidUV(:), atm_tidW(:)      ! per-cell turbine id, static
real(rprec), allocatable :: atm_gx(:), atm_gy(:), atm_gz(:)
integer,     allocatable :: atm_awi(:), atm_awj(:)

! ---- Batched Cl/tip correction (round 4) ----
! GPU port of atm_compute_cl_correction: the O(N^2) blade-to-blade induced
! velocity loop (~68 ms/step on the host) batched over all turbines. The
! model arithmetic is replicated verbatim (same mixed-precision literals,
! same accumulation order); all turbineArray outputs are copied back to the
! host each step so host state stays exactly as the host routine leaves it.
logical :: atm_clc_ready = .false.
integer,     allocatable :: atm_pt_turb(:), atm_pt_q(:), atm_pt_base(:), atm_pt_qq(:)
real(rprec), allocatable :: atm_chord_all(:), atm_brad_all(:)      ! static
real(rprec), allocatable :: atm_db_all(:)                      ! static blade-section width
real(rprec), allocatable :: atm_clc_tc(:,:)                       ! (3,nTurb): eps_s, optEpsChord, active
real(rprec), allocatable :: atm_wv_all(:,:)                       ! (3,nbp) windVectors in
real(rprec), allocatable :: atm_cl_all(:), atm_cd_all(:), atm_vmag_all(:)
real(rprec), allocatable :: atm_du_all(:,:)                       ! (3,nbp) state in/out
real(rprec), allocatable :: atm_uyopt_vec_all(:,:)                ! (3,nbp) state in/out
real(rprec), allocatable :: atm_uinf_all(:,:), atm_uxles_all(:,:) ! (3,nbp) out
real(rprec), allocatable :: atm_g_all(:), atm_dg_all(:), atm_epsopt_all(:)
real(rprec), allocatable :: atm_uyles_vec_all(:,:)                ! (3,nbp) out
real(rprec), allocatable :: atm_uyles_all(:), atm_uyopt_all(:)    ! magnitudes out
#endif

! Body force field
type(bodyForce_t), allocatable, target, dimension(:) :: forceFieldUV, forceFieldW

! The very crucial parameter pi
real(rprec), parameter :: pi=acos(-1._rprec)

type(clock_t), save :: atm_clock_interp_w, atm_clock_update, atm_clock_reset
type(clock_t), save :: atm_clock_sample, atm_clock_force, atm_clock_gather, atm_clock_convolve
type(clock_t), save :: atm_clock_clcorr, atm_clock_apply, atm_clock_output
type(clock_t), save :: atm_clock_barrier
real(rprec), save :: atm_time_interp_w = 0._rprec
real(rprec), save :: atm_time_update = 0._rprec
real(rprec), save :: atm_time_reset = 0._rprec
real(rprec), save :: atm_time_sample = 0._rprec
real(rprec), save :: atm_time_force = 0._rprec
real(rprec), save :: atm_time_gather = 0._rprec
real(rprec), save :: atm_time_convolve = 0._rprec
real(rprec), save :: atm_time_clcorr = 0._rprec
real(rprec), save :: atm_time_apply = 0._rprec
real(rprec), save :: atm_time_output = 0._rprec
real(rprec), save :: atm_time_barrier = 0._rprec
integer, save :: atm_forcing_calls = 0
logical, save :: atm_diag_load_printed = .false.
#ifdef ENABLE_CUDA
type(cudaEvent), save :: atm_diag_evt_start, atm_diag_evt_stop
logical, save :: atm_diag_evt_created = .false.
logical, save :: atm_diag_evt_pending = .false.
integer, save :: atm_diag_evt_pending_call = 0
integer, save :: atm_lb_npoints = 0
logical, save :: atm_lb_owner_ready = .false.
logical, save :: atm_lb_setup_printed = .false.
logical, save :: atm_lb_auto_decision_made = .false.
logical, save :: atm_lb_auto_use_lb = .false.
logical, save :: atm_lb_auto_printed = .false.
integer, managed, allocatable, save :: atm_lb_owner(:)
integer, managed, allocatable, save :: atm_lb_offset(:)
integer, managed, allocatable, save :: atm_lb_grid_owner(:)
integer, managed, allocatable, save :: atm_lb_sample_send_slot(:)
integer, managed, allocatable, save :: atm_lb_sample_recv_slot(:)
integer, managed, allocatable, save :: atm_lb_force_send_slot(:)
integer, managed, allocatable, save :: atm_lb_force_recv_slot(:)
real(rprec), device, allocatable, save :: atm_lb_vel_local(:,:)
real(rprec), device, allocatable, save :: atm_lb_vel_global(:,:)
real(rprec), device, allocatable, save :: atm_lb_sample_sendbuf(:)
real(rprec), device, allocatable, save :: atm_lb_sample_recvbuf(:)
real(rprec), device, allocatable, save :: atm_lb_force_sendbuf(:)
real(rprec), device, allocatable, save :: atm_lb_force_recvbuf(:)
real(rprec), save :: atm_lb_time_sample = 0._rprec
real(rprec), save :: atm_lb_time_velocity_comm = 0._rprec
real(rprec), save :: atm_lb_time_point_force = 0._rprec
real(rprec), save :: atm_lb_time_gather = 0._rprec
real(rprec), save :: atm_lb_time_validate = 0._rprec
real(rprec), save :: atm_lb_time_sample_pack = 0._rprec
real(rprec), save :: atm_lb_time_sample_mpi = 0._rprec
real(rprec), save :: atm_lb_time_sample_unpack = 0._rprec
real(rprec), save :: atm_lb_time_force_pack = 0._rprec
real(rprec), save :: atm_lb_time_force_mpi = 0._rprec
real(rprec), save :: atm_lb_time_force_unpack = 0._rprec
real(rprec), save :: atm_lb_time_scalar_reduce = 0._rprec
integer, save :: atm_lb_sample_send_count = 0
integer, save :: atm_lb_sample_recv_count = 0
integer, save :: atm_lb_force_send_count = 0
integer, save :: atm_lb_force_recv_count = 0
integer, save :: atm_lb_targeted_calls = 0
integer, save :: atm_lb_auto_probe_calls = 0
integer, save :: atm_lb_auto_legacy_samples = 0
integer, save :: atm_lb_auto_lb_samples = 0
real(rprec), save :: atm_lb_planned_points_rank = 0._rprec
real(rprec), save :: atm_lb_planned_support_rank = 0._rprec
real(rprec), save :: atm_lb_planned_scatter_rank = 0._rprec
real(rprec), save :: atm_lb_auto_legacy_time = 0._rprec
real(rprec), save :: atm_lb_auto_lb_time = 0._rprec
integer(8), save :: atm_lb_sample_bytes_sent = 0_8
integer(8), save :: atm_lb_sample_bytes_recv = 0_8
integer(8), save :: atm_lb_force_bytes_sent = 0_8
integer(8), save :: atm_lb_force_bytes_recv = 0_8
#endif

contains

#ifndef ENABLE_CUDA
! ===== Option A (full-OpenACC core + CPU ATM): CPU build has no CUDA-Fortran, so
! these cuf gate functions are absent. Provide .false. stubs so the ATM takes its
! pure-CPU paths. Physics/model preserved (the CPU paths are the existing
! actuator-line model); only the ATM's GPU offload is off. Phase C later converts
! these CPU ATM paths to OpenACC.
logical function atm_wuv_cuda_enabled();          atm_wuv_cuda_enabled=.true.;           end function
#ifdef PPLES_GPU
logical function atm_direct_w_enabled();          atm_direct_w_enabled=.true.;           end function
#else
logical function atm_direct_w_enabled();          atm_direct_w_enabled=.false.;          end function
#endif
logical function atm_apply_cuda_enabled();        atm_apply_cuda_enabled=.false.;        end function
logical function atm_convolve_cuda_enabled();     atm_convolve_cuda_enabled=.false.;     end function
logical function atm_force_shadows_enabled();     atm_force_shadows_enabled=.false.;     end function
logical function atm_bladeforce_cuda_enabled();   atm_bladeforce_cuda_enabled=.false.;   end function
logical function atm_reset_cuda_enabled();        atm_reset_cuda_enabled=.false.;        end function
logical function atm_packed_gather_enabled();     atm_packed_gather_enabled=.false.;     end function
logical function atm_gpu_packed_gather_enabled(); atm_gpu_packed_gather_enabled=.false.; end function
logical function atm_slim_gather_enabled();       atm_slim_gather_enabled=.false.;       end function
logical function atm_batch_gather_enabled();      atm_batch_gather_enabled=.false.;      end function
logical function atm_full_gather_required();      atm_full_gather_required=.false.;      end function
logical function atm_output_barrier_enabled();    atm_output_barrier_enabled=.false.;    end function
logical function atm_extra_sync_enabled();        atm_extra_sync_enabled=.false.;        end function
logical function atm_diag_timing_enabled();       atm_diag_timing_enabled=.false.;       end function
logical function atm_point_owner_lb_enabled();    atm_point_owner_lb_enabled=.false.;    end function
logical function atm_point_owner_targeted_enabled();  atm_point_owner_targeted_enabled=.false.;  end function
logical function atm_lb_auto_select_enabled();    atm_lb_auto_select_enabled=.false.;    end function
logical function atm_lb_validate_enabled();       atm_lb_validate_enabled=.false.;       end function
logical function atm_point_owner_lb_supported();  atm_point_owner_lb_supported=.false.;  end function
logical function atm_point_owner_targeted_supported(); atm_point_owner_targeted_supported=.false.; end function
logical function atm_lb_plan_only_enabled();      atm_lb_plan_only_enabled=.false.;      end function
logical function atm_lb_point_detail_enabled();   atm_lb_point_detail_enabled=.false.;   end function
logical function atm_lb_auto_use_lb_for_call(candidate)
    logical, intent(in) :: candidate
    atm_lb_auto_use_lb_for_call = .false.
end function
#endif

#ifdef ENABLE_CUDA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_wuv_cuda_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_wuv_cuda_enabled = .true.

end function atm_wuv_cuda_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_direct_w_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_direct_w_enabled = .true.

end function atm_direct_w_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_apply_cuda_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_apply_cuda_enabled = .true.

end function atm_apply_cuda_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_convolve_cuda_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_convolve_cuda_enabled = .true.

end function atm_convolve_cuda_enabled

! ATM execution-policy gates.  Keep environment-controlled gates synchronized
! with docs/environment_switches.md; prefer fixed defaults for validated
! production paths and reserve env switches for diagnostics/experiments.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_env_flag(name, default_value)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(len=*), intent(in) :: name
logical, intent(in) :: default_value
character(len=32) :: env_value
integer :: env_status

atm_env_flag = default_value
env_value = ''
call get_environment_variable(name, env_value, status=env_status)
if (env_status == 0) then
    env_value = adjustl(env_value)
    select case (env_value(1:1))
    case ('1', 't', 'T', 'y', 'Y')
        atm_env_flag = .true.
    case ('0', 'f', 'F', 'n', 'N')
        atm_env_flag = .false.
    end select
end if

end function atm_env_flag

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_force_shadows_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .true.

if (.not. initialized) then
    enabled = atm_env_flag('LESGO_ATM_FORCE_SHADOWS', .true.)
    initialized = .true.
end if

atm_force_shadows_enabled = enabled

end function atm_force_shadows_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_bladeforce_cuda_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_bladeforce_cuda_enabled = .true.

end function atm_bladeforce_cuda_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_reset_cuda_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_reset_cuda_enabled = .true.

end function atm_reset_cuda_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_packed_gather_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_packed_gather_enabled = .true.

end function atm_packed_gather_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_gpu_packed_gather_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_gpu_packed_gather_enabled = .true.

end function atm_gpu_packed_gather_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_slim_gather_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_slim_gather_enabled = .true.

end function atm_slim_gather_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_batch_gather_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Experimental batched reduction was correct but not faster on the current
! 2-rank Derecho case, so it is kept out of the production path.
!
implicit none

atm_batch_gather_enabled = .false.

end function atm_batch_gather_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_full_gather_required()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i

! Full field gather is only needed for heavy blade/velocity diagnostics, not
! for scalar turbine power/thrust output. Default diagnosticOutputInterval is
! outputInterval, preserving legacy behavior unless the case opts out.
atm_full_gather_required = .false.
if (diagnosticOutputInterval > 0) then
    atm_full_gather_required = mod(jt_total-1, diagnosticOutputInterval) == 0
endif

do i = 1, numberOfTurbines
    if (turbineArray(i) % tipALMCorrection) then
        atm_full_gather_required = .true.
    endif
enddo

end function atm_full_gather_required

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_output_barrier_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_output_barrier_enabled = .true.

end function atm_output_barrier_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_extra_sync_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_extra_sync_enabled = .false.

end function atm_extra_sync_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_diag_timing_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = atm_env_flag('LESGO_ATM_DIAG_TIMING', .false.)
    initialized = .true.
end if

atm_diag_timing_enabled = enabled

end function atm_diag_timing_enabled

! Point-owner load balancing remains experimental.  The auto-selector probes
! legacy and LB paths, then keeps LB only when its measured cost is lower.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_point_owner_lb_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = atm_env_flag('LESGO_ATM_POINT_OWNER_LB', .false.)
    initialized = .true.
end if

atm_point_owner_lb_enabled = enabled

end function atm_point_owner_lb_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_point_owner_targeted_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = atm_env_flag('LESGO_ATM_POINT_OWNER_TARGETED', .false.)
    initialized = .true.
end if

atm_point_owner_targeted_enabled = enabled

end function atm_point_owner_targeted_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_lb_auto_select_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = atm_env_flag('LESGO_ATM_LB_AUTO_SELECT', .false.)
    initialized = .true.
end if

atm_lb_auto_select_enabled = enabled

end function atm_lb_auto_select_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function atm_lb_auto_probe_steps()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_lb_auto_probe_steps = 10

end function atm_lb_auto_probe_steps

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) function atm_lb_auto_margin()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_lb_auto_margin = 0._rprec

end function atm_lb_auto_margin

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_lb_auto_use_lb_for_call(candidate)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, intent(in) :: candidate

if (.not. candidate) then
    atm_lb_auto_use_lb_for_call = .false.
elseif (.not. atm_lb_auto_select_enabled()) then
    atm_lb_auto_use_lb_for_call = .true.
elseif (atm_lb_auto_decision_made) then
    atm_lb_auto_use_lb_for_call = atm_lb_auto_use_lb
else
    atm_lb_auto_use_lb_for_call = mod(atm_lb_auto_probe_calls, 2) == 1
endif

end function atm_lb_auto_use_lb_for_call

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lb_auto_record(used_lb, elapsed)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, intent(in) :: used_lb
real(rprec), intent(in) :: elapsed
integer :: probe_steps
real(rprec) :: avg_legacy, avg_lb, margin

if (.not. atm_lb_auto_select_enabled()) return
if (atm_lb_auto_decision_made) return

atm_lb_auto_probe_calls = atm_lb_auto_probe_calls + 1
if (used_lb) then
    atm_lb_auto_lb_samples = atm_lb_auto_lb_samples + 1
    atm_lb_auto_lb_time = atm_lb_auto_lb_time + elapsed
else
    atm_lb_auto_legacy_samples = atm_lb_auto_legacy_samples + 1
    atm_lb_auto_legacy_time = atm_lb_auto_legacy_time + elapsed
endif

probe_steps = atm_lb_auto_probe_steps()
if (atm_lb_auto_probe_calls < probe_steps) return
if (atm_lb_auto_legacy_samples <= 0 .or. atm_lb_auto_lb_samples <= 0) return

avg_legacy = atm_lb_auto_legacy_time / real(atm_lb_auto_legacy_samples, rprec)
avg_lb = atm_lb_auto_lb_time / real(atm_lb_auto_lb_samples, rprec)
margin = atm_lb_auto_margin()
atm_lb_auto_use_lb = avg_lb < avg_legacy * (1._rprec - margin)
atm_lb_auto_decision_made = .true.

if (coord == 0) then
    write(*,'(a,i0,a,i0,a,i0,3(a,es16.8),a)')                               &
        'ATM_LB_AUTO_DECISION probe_calls=', atm_lb_auto_probe_calls,        &
        ' legacy_samples=', atm_lb_auto_legacy_samples,                      &
        ' lb_samples=', atm_lb_auto_lb_samples,                              &
        ' avg_legacy=', avg_legacy, ' avg_lb=', avg_lb,                      &
        ' margin=', margin, ' selected=' //                                  &
        merge('point_owner_lb', 'legacy        ', atm_lb_auto_use_lb)
    flush(6)
endif

end subroutine atm_lb_auto_record

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_lb_validate_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = atm_env_flag('LESGO_ATM_LB_VALIDATE', .false.)
    initialized = .true.
end if

atm_lb_validate_enabled = enabled

end function atm_lb_validate_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_diag_event_start(call_id)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: call_id
integer :: istat
real :: elapsed_ms
real(rprec) :: elapsed_s

if (.not. atm_diag_timing_enabled()) return

if (.not. atm_diag_evt_created) then
    istat = cudaEventCreate(atm_diag_evt_start)
    if (istat /= cudaSuccess) stop 'ATM diagnostic event create failure'
    istat = cudaEventCreate(atm_diag_evt_stop)
    if (istat /= cudaSuccess) stop 'ATM diagnostic event create failure'
    atm_diag_evt_created = .true.
end if

if (atm_diag_evt_pending) then
    istat = cudaEventSynchronize(atm_diag_evt_stop)
    if (istat /= cudaSuccess) stop 'ATM diagnostic event sync failure'
    istat = cudaEventElapsedTime(elapsed_ms, atm_diag_evt_start,              &
        atm_diag_evt_stop)
    if (istat /= cudaSuccess) stop 'ATM diagnostic event elapsed failure'
    elapsed_s = real(elapsed_ms, rprec) * 1.0e-3_rprec
    write(*,'(a,i0,a,i0,a,es14.6,a,i0,a,es14.6)')                            &
        'ATM_GPU_TIMING call=', atm_diag_evt_pending_call, ' rank=', coord,   &
        ' gpu_event=', elapsed_s, ' kernels=', -1, ' sync=', -1.0_rprec
    atm_diag_evt_pending = .false.
end if

istat = cudaEventRecord(atm_diag_evt_start, 0)
if (istat /= cudaSuccess) stop 'ATM diagnostic event start record failure'

end subroutine atm_diag_event_start

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_diag_event_stop(call_id)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: call_id
integer :: istat

if (.not. atm_diag_timing_enabled()) return
if (.not. atm_diag_evt_created) return

istat = cudaEventRecord(atm_diag_evt_stop, 0)
if (istat /= cudaSuccess) stop 'ATM diagnostic event stop record failure'
atm_diag_evt_pending = .true.
atm_diag_evt_pending_call = call_id

end subroutine atm_diag_event_stop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_diag_event_flush()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: istat
real :: elapsed_ms
real(rprec) :: elapsed_s

if (.not. atm_diag_timing_enabled()) return
if (.not. atm_diag_evt_pending) return

istat = cudaEventSynchronize(atm_diag_evt_stop)
if (istat /= cudaSuccess) stop 'ATM diagnostic event final sync failure'
istat = cudaEventElapsedTime(elapsed_ms, atm_diag_evt_start,                  &
    atm_diag_evt_stop)
if (istat /= cudaSuccess) stop 'ATM diagnostic event final elapsed failure'
elapsed_s = real(elapsed_ms, rprec) * 1.0e-3_rprec
write(*,'(a,i0,a,i0,a,es14.6,a,i0,a,es14.6)')                                &
    'ATM_GPU_TIMING call=', atm_diag_evt_pending_call, ' rank=', coord,       &
    ' gpu_event=', elapsed_s, ' kernels=', -1, ' sync=', -1.0_rprec
atm_diag_evt_pending = .false.

end subroutine atm_diag_event_flush

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_cuda_check(where)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: where
integer :: istat

if (atm_extra_sync_enabled()) then
    call atm_cuda_sync(where)
    return
end if

istat = cudaGetLastError()
if (istat /= cudaSuccess) then
    print *, 'ATM CUDA kernel failure at ', trim(where), ': ', istat
    stop 1
end if

end subroutine atm_cuda_check

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_cuda_sync(where)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= cudaSuccess) then
    print *, 'ATM CUDA sync failure at ', trim(where), ': ', istat
    stop 1
end if
istat = cudaGetLastError()
if (istat /= cudaSuccess) then
    print *, 'ATM CUDA kernel failure at ', trim(where), ': ', istat
    stop 1
end if

end subroutine atm_cuda_sync

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_reset_turbine_gpu(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, m, n, q, a, b, mmend, nnend, qqend
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladeForces
real(rprec), managed, pointer, dimension(:,:,:,:) :: integratedBladeForces
real(rprec), managed, pointer, dimension(:,:,:,:) :: windVectors
real(rprec), managed, pointer, dimension(:,:,:,:,:) :: bladeAlignedVectors
real(rprec), managed, pointer, dimension(:,:,:) :: alpha, Cd, Cm, Cl, Cl_b, G
real(rprec), managed, pointer, dimension(:,:,:) :: lift, drag, Vmag
real(rprec), managed, pointer, dimension(:,:,:) :: induction_a, u_infinity
real(rprec), managed, pointer, dimension(:,:,:) :: axialForce, tangentialForce
real(rprec), managed, pointer, dimension(:,:,:) :: pitchingMoment

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints

bladeForces => turbineArray(i) % bladeForces
integratedBladeForces => turbineArray(i) % integratedBladeForces
windVectors => turbineArray(i) % windVectors
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
alpha => turbineArray(i) % alpha
Cd => turbineArray(i) % Cd
Cm => turbineArray(i) % Cm
Cl => turbineArray(i) % Cl
Cl_b => turbineArray(i) % Cl_b
G => turbineArray(i) % G
lift => turbineArray(i) % lift
drag => turbineArray(i) % drag
Vmag => turbineArray(i) % Vmag
induction_a => turbineArray(i) % induction_a
u_infinity => turbineArray(i) % u_infinity
axialForce => turbineArray(i) % axialForce
tangentialForce => turbineArray(i) % tangentialForce
pitchingMoment => turbineArray(i) % pitchingMoment

!$cuf kernel do(3) <<<*,*>>>
do q = 1, qqend
    do n = 1, nnend
        do m = 1, mmend
            do a = 1, 3
                bladeForces(m,n,q,a) = 0._rprec
                integratedBladeForces(m,n,q,a) = 0._rprec
                windVectors(m,n,q,a) = 0._rprec
                do b = 1, 3
                    bladeAlignedVectors(m,n,q,a,b) = 0._rprec
                enddo
            enddo

            alpha(m,n,q) = 0._rprec
            Cd(m,n,q) = 0._rprec
            Cm(m,n,q) = 0._rprec
            Cl(m,n,q) = 0._rprec
            Cl_b(m,n,q) = 0._rprec
            G(m,n,q) = 0._rprec
            lift(m,n,q) = 0._rprec
            drag(m,n,q) = 0._rprec
            Vmag(m,n,q) = 0._rprec
            induction_a(m,n,q) = 0._rprec
            u_infinity(m,n,q) = 0._rprec
            axialForce(m,n,q) = 0._rprec
            tangentialForce(m,n,q) = 0._rprec
            pitchingMoment(m,n,q) = 0._rprec
        enddo
    enddo
enddo

call atm_cuda_check('ATM turbine reset')

end subroutine atm_lesgo_reset_turbine_gpu
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_interp_w_to_uv()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: jx, jy, jz
#if defined(PPLES_GPU) && defined(PPMPI)
integer :: status(MPI_STATUS_SIZE)
#endif

#ifdef PPLES_GPU
if (atm_wuv_cuda_enabled()) then
    !$acc parallel loop collapse(3) default(present)
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, nx
        w_uv(jx,jy,jz) = 0.5_rprec * (w(jx,jy,jz+1) + w(jx,jy,jz))
    end do
    end do
    end do

    if (coord == nproc - 1) then
        !$acc parallel loop collapse(2) default(present)
        do jy = 1, ny
        do jx = 1, nx
            w_uv(jx,jy,nz) = w_uv(jx,jy,nz-1)
        end do
        end do
    end if

#ifdef PPMPI
    if (nproc > 1 .and. .not. atm_direct_w_enabled()) then
        ! GPU-aware MPI halo for w_uv: exchange device pointers directly (no
        ! full-field update self/device PCIe). Replicates
        ! mpi_sync_real_array(w_uv,lbz,DOWNUP): send :,:,1 down / recv :,:,nz from
        ! up, then (lbz==0) send :,:,nz-1 up / recv :,:,0 from down.
        !
        ! NOTE (overlap): no `wait(1)` here. The interpolation kernels above
        ! are synchronous (no async clause), so w_uv is complete; queue 1 only
        ! carries the deferred convection/RHS kernels at this point in the
        ! step, none of which touch w or w_uv. Draining it would serialize the
        ! ATM phase-1 host work against the convection backlog.
        !$acc host_data use_device(w_uv)
        call mpi_sendrecv(w_uv(1,1,1),  nx*ny, mpi_rprec, down, 1,                &
                          w_uv(1,1,nz), nx*ny, mpi_rprec, up,   1, comm, status, ierr)
        if (lbz == 0) then
            call mpi_sendrecv(w_uv(1,1,nz-1), nx*ny, mpi_rprec, up,   2,          &
                              w_uv(1,1,0),    nx*ny, mpi_rprec, down, 2, comm, status, ierr)
        endif
        !$acc end host_data
    endif
#endif
    return
end if
#elif defined(ENABLE_CUDA)
if (atm_wuv_cuda_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, nx
        w_uv(jx,jy,jz) = 0.5_rprec * (w(jx,jy,jz+1) + w(jx,jy,jz))
    end do
    end do
    end do

    if (coord == nproc - 1) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, nx
            w_uv(jx,jy,nz) = w_uv(jx,jy,nz-1)
        end do
        end do
    end if

#ifdef PPMPI
    if (nproc > 1) then
        call atm_cuda_sync('w -> uv interpolation')
        if (lbz == 0) then
            call mpi_sync_real_array(w_uv, lbz, MPI_SYNC_DOWNUP)
        elseif (lbz == 1) then
            call mpi_sync_real_array(w_uv, lbz, MPI_SYNC_DOWN)
        endif
    else
        call atm_cuda_check('w -> uv interpolation')
    endif
#else
    call atm_cuda_check('w -> uv interpolation')
#endif
    return
end if
#endif

w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

end subroutine atm_interp_w_to_uv

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_prepare_direct_w()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: jx, jy

#ifdef PPMPI
integer :: status(MPI_STATUS_SIZE)
#endif

if (.not. atm_direct_w_enabled()) return

if (.not. allocated(atm_wuv_send_down)) then
    allocate(atm_wuv_send_down(nx,ny), atm_wuv_recv_up(nx,ny))
#ifdef PPLES_GPU
    !$acc update device(atm_wuv_send_down, atm_wuv_recv_up)
#endif
endif

#ifdef PPMPI
if (nproc <= 1) then
#ifdef PPLES_GPU
    !$acc wait
#else
    call atm_cuda_check('ATM direct w single rank')
#endif
    return
endif

! The direct sampler only needs the old w_uv(nz) halo from the rank above.
#ifdef PPLES_GPU
!$acc parallel loop collapse(2) default(present)
do jy = 1, ny
do jx = 1, nx
    atm_wuv_send_down(jx,jy) = 0.5_rprec * (w(jx,jy,1) + w(jx,jy,2))
end do
end do

! The boundary-pack loop has no async clause, so it is complete before the
! following GPU-aware MPI call.  Avoid a global OpenACC wait here because it can
! charge unrelated deferred LES work to the ATM direct-w timer.
#else
!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, nx
    atm_wuv_send_down(jx,jy) = 0.5_rprec * (w(jx,jy,1) + w(jx,jy,2))
end do
end do

call atm_cuda_sync('ATM direct w boundary pack')
#endif

#ifdef PPLES_GPU
!$acc host_data use_device(atm_wuv_send_down, atm_wuv_recv_up)
call mpi_sendrecv(atm_wuv_send_down(1,1), nx*ny, mpi_rprec, down, 991,        &
                  atm_wuv_recv_up(1,1), nx*ny, mpi_rprec, up, 991,            &
                  comm, status, ierr)
!$acc end host_data
#else
call mpi_sendrecv(atm_wuv_send_down(1,1), nx*ny, mpi_rprec, down, 991,        &
                  atm_wuv_recv_up(1,1), nx*ny, mpi_rprec, up, 991,            &
                  comm, status, ierr)
#endif
if (ierr /= 0) stop 'ATM direct w device halo exchange failed'

#ifdef PPLES_GPU
! Blocking GPU-aware MPI sendrecv completes the device receive buffer before the
! following ATM sampling kernels are launched. Avoid a global OpenACC wait here:
! it can drain unrelated deferred LES work into the direct-w timer.
#else
call atm_cuda_check('ATM direct w boundary exchange')
#endif
#else
#ifdef PPLES_GPU
!$acc wait
#else
call atm_cuda_check('ATM direct w no mpi')
#endif
#endif

end subroutine atm_prepare_direct_w
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_apply_force_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: c, m
integer :: i, j, k
#ifdef PPLES_GPU
integer :: totalUV, totalW
#endif

#ifdef PPLES_GPU
if (atm_force_shadows_enabled() .and. atm_apply_cuda_enabled() .and.          &
    atm_shadow_ready) then
    totalUV = atm_shadow_offsetUV(numberOfTurbines+1)
    totalW = atm_shadow_offsetW(numberOfTurbines+1)

    if (totalUV > 0) then
        !$cuf kernel do(1) <<<*,*>>>
        do c = 1, totalUV
            i = atm_shadow_ijkUV_d(1,c)
            j = atm_shadow_ijkUV_d(2,c)
            k = atm_shadow_ijkUV_d(3,c)
            fxa(i,j,k) = fxa(i,j,k) + atm_shadow_forceUV_d(1,c)
            fya(i,j,k) = fya(i,j,k) + atm_shadow_forceUV_d(2,c)
        end do
    end if

    if (totalW > 0) then
        !$cuf kernel do(1) <<<*,*>>>
        do c = 1, totalW
            i = atm_shadow_ijkW_d(1,c)
            j = atm_shadow_ijkW_d(2,c)
            k = atm_shadow_ijkW_d(3,c)
            fza(i,j,k) = fza(i,j,k) + atm_shadow_forceW_d(3,c)
        end do
    end if

    call atm_cuda_check('apply ATM force shadows')
    return
end if
#endif
#ifdef ENABLE_CUDA
if (.not. atm_apply_cuda_enabled()) return

do m = 1, numberOfTurbines
    if (turbineArray(m) % operate) then
        !$cuf kernel do(1) <<<*,*>>>
        do c = 1, forceFieldUV(m) % c
            i = forceFieldUV(m) % ijk(1,c)
            j = forceFieldUV(m) % ijk(2,c)
            k = forceFieldUV(m) % ijk(3,c)

            fxa(i,j,k) = fxa(i,j,k) + forceFieldUV(m) % force(1,c)
            fya(i,j,k) = fya(i,j,k) + forceFieldUV(m) % force(2,c)
        end do

        !$cuf kernel do(1) <<<*,*>>>
        do c = 1, forceFieldW(m) % c
            i = forceFieldW(m) % ijk(1,c)
            j = forceFieldW(m) % ijk(2,c)
            k = forceFieldW(m) % ijk(3,c)

            fza(i,j,k) = fza(i,j,k) + forceFieldW(m) % force(3,c)
        end do
    end if
end do

call atm_cuda_check('apply ATM force')
#endif

end subroutine atm_lesgo_apply_force_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_convolute_force_gpu_atpoint(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i

integer :: j, m, n, q, c, mmend, nnend, qqend, cUV, cW, has_nacelle
real(rprec) :: projectradius, projectradius_sq, epsilon, epsilon_sq
real(rprec) :: nacelle1, nacelle2, nacelle3, nacelle_epsilon
real(rprec) :: nacelle_epsilon_sq, nacelle_kernel_norm
real(rprec) :: nacelleForce1, nacelleForce2, nacelleForce3
real(rprec) :: const1, const2, const3
real(rprec) :: a1, a2, a3, b1, b2, b3, d1, d2, d3, dist_sq, kernel
real(rprec) :: f1, f2, f3
#ifdef PPLES_GPU
integer :: offUV, offW
#endif
#ifdef ENABLE_CUDA
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladeForces, bladePoints
real(rprec), managed, pointer, dimension(:,:) :: bodyForceUV, bodyForceW
real(rprec), managed, pointer, dimension(:,:) :: locationUV, locationW
#endif

#ifdef PPLES_GPU
if (atm_force_shadows_enabled() .and. atm_convolve_cuda_enabled() .and.       &
    atm_shadow_ready .and. atm_blade_mirror_ready) then
    j = turbineArray(i) % turbineTypeID
    mmend = turbineModel(j) % numBl
    nnend = turbineArray(i) % numAnnulusSections
    qqend = turbineArray(i) % numBladePoints
    cUV = forceFieldUV(i) % c
    cW = forceFieldW(i) % c
    offUV = atm_shadow_offsetUV(i)
    offW = atm_shadow_offsetW(i)

    projectradius = turbineArray(i) % projectionRadius
    projectradius_sq = projectradius * projectradius
    epsilon = turbineArray(i) % epsilon
    epsilon_sq = epsilon * epsilon
    const1 = 1._rprec / ((epsilon * epsilon_sq) * (pi * sqrt(pi)))
    const2 = z_i / (u_star * u_star)
    const3 = const1 * const2
    if (turbineArray(i) % nacelle .and. turbineArray(i) % nacelleEpsilon > 0._rprec) then
        has_nacelle = 1
    else
        has_nacelle = 0
    endif
    nacelle1 = turbineArray(i) % nacelleLocation(1)
    nacelle2 = turbineArray(i) % nacelleLocation(2)
    nacelle3 = turbineArray(i) % nacelleLocation(3)
    nacelleForce1 = turbineArray(i) % nacelleForce(1)
    nacelleForce2 = turbineArray(i) % nacelleForce(2)
    nacelleForce3 = turbineArray(i) % nacelleForce(3)
    nacelle_epsilon = turbineArray(i) % nacelleEpsilon
    nacelle_epsilon_sq = nacelle_epsilon * nacelle_epsilon
    nacelle_kernel_norm = 1._rprec /                                         &
        ((nacelle_epsilon * nacelle_epsilon_sq) * (pi * sqrt(pi)))

    if (cUV > 0) then
        !$cuf kernel do(1) <<<*,*>>>
        do c = 1, cUV
            a1 = atm_shadow_locUV_d(1,offUV+c)
            a2 = atm_shadow_locUV_d(2,offUV+c)
            a3 = atm_shadow_locUV_d(3,offUV+c)
            f1 = 0._rprec
            f2 = 0._rprec

            do m = 1, mmend
            do n = 1, nnend
            do q = 1, qqend
                b1 = atm_bladePoints_d(m,n,q,1,i)
                b2 = atm_bladePoints_d(m,n,q,2,i)
                b3 = atm_bladePoints_d(m,n,q,3,i)
                d1 = a1 - b1
                d2 = a2 - b2
                d3 = a3 - b3
                dist_sq = d1*d1 + d2*d2 + d3*d3

                if (dist_sq <= projectradius_sq) then
                    kernel = exp(-dist_sq / epsilon_sq)
                    f1 = f1 + atm_bladeForces_d(m,n,q,1,i) * kernel
                    f2 = f2 + atm_bladeForces_d(m,n,q,2,i) * kernel
                end if
            end do
            end do
            end do

            f1 = f1 * const3
            f2 = f2 * const3

            if (has_nacelle == 1) then
                d1 = a1 - nacelle1
                d2 = a2 - nacelle2
                d3 = a3 - nacelle3
                dist_sq = d1*d1 + d2*d2 + d3*d3
                kernel = exp(-dist_sq / nacelle_epsilon_sq) * nacelle_kernel_norm
                f1 = f1 + nacelleForce1 * kernel * const2
                f2 = f2 + nacelleForce2 * kernel * const2
            endif

            atm_shadow_forceUV_d(1,offUV+c) = f1
            atm_shadow_forceUV_d(2,offUV+c) = f2
            atm_shadow_forceUV_d(3,offUV+c) = 0._rprec
        end do
    end if

    if (cW > 0) then
        !$cuf kernel do(1) <<<*,*>>>
        do c = 1, cW
            a1 = atm_shadow_locW_d(1,offW+c)
            a2 = atm_shadow_locW_d(2,offW+c)
            a3 = atm_shadow_locW_d(3,offW+c)
            f3 = 0._rprec

            do m = 1, mmend
            do n = 1, nnend
            do q = 1, qqend
                b1 = atm_bladePoints_d(m,n,q,1,i)
                b2 = atm_bladePoints_d(m,n,q,2,i)
                b3 = atm_bladePoints_d(m,n,q,3,i)
                d1 = a1 - b1
                d2 = a2 - b2
                d3 = a3 - b3
                dist_sq = d1*d1 + d2*d2 + d3*d3

                if (dist_sq <= projectradius_sq) then
                    kernel = exp(-dist_sq / epsilon_sq)
                    f3 = f3 + atm_bladeForces_d(m,n,q,3,i) * kernel
                end if
            end do
            end do
            end do

            f3 = f3 * const3

            if (has_nacelle == 1) then
                d1 = a1 - nacelle1
                d2 = a2 - nacelle2
                d3 = a3 - nacelle3
                dist_sq = d1*d1 + d2*d2 + d3*d3
                if (dist_sq <= projectradius_sq) then
                    kernel = exp(-dist_sq / nacelle_epsilon_sq) * nacelle_kernel_norm
                    f3 = f3 + nacelleForce3 * kernel * const2
                endif
            endif

            atm_shadow_forceW_d(1,offW+c) = 0._rprec
            atm_shadow_forceW_d(2,offW+c) = 0._rprec
            atm_shadow_forceW_d(3,offW+c) = f3
        end do
    end if

    call atm_cuda_check('ATM force shadow convolution')
    return
end if
#endif
#ifdef ENABLE_CUDA
j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
cUV = forceFieldUV(i) % c
cW = forceFieldW(i) % c

projectradius = turbineArray(i) % projectionRadius
projectradius_sq = projectradius * projectradius
epsilon = turbineArray(i) % epsilon
epsilon_sq = epsilon * epsilon
const1 = 1._rprec / ((epsilon * epsilon_sq) * (pi * sqrt(pi)))
const2 = z_i / (u_star * u_star)
const3 = const1 * const2
if (turbineArray(i) % nacelle .and. turbineArray(i) % nacelleEpsilon > 0._rprec) then
    has_nacelle = 1
else
    has_nacelle = 0
endif
nacelle1 = turbineArray(i) % nacelleLocation(1)
nacelle2 = turbineArray(i) % nacelleLocation(2)
nacelle3 = turbineArray(i) % nacelleLocation(3)
nacelleForce1 = turbineArray(i) % nacelleForce(1)
nacelleForce2 = turbineArray(i) % nacelleForce(2)
nacelleForce3 = turbineArray(i) % nacelleForce(3)
nacelle_epsilon = turbineArray(i) % nacelleEpsilon
nacelle_epsilon_sq = nacelle_epsilon * nacelle_epsilon
nacelle_kernel_norm = 1._rprec /                                             &
    ((nacelle_epsilon * nacelle_epsilon_sq) * (pi * sqrt(pi)))

bladeForces => turbineArray(i) % bladeForces
bladePoints => turbineArray(i) % bladePoints
bodyForceUV => forceFieldUV(i) % force
bodyForceW => forceFieldW(i) % force
locationUV => forceFieldUV(i) % location
locationW => forceFieldW(i) % location

!$cuf kernel do(1) <<<*,*>>>
do c = 1, cUV
    a1 = locationUV(1,c)
    a2 = locationUV(2,c)
    a3 = locationUV(3,c)
    f1 = 0._rprec
    f2 = 0._rprec

    do m = 1, mmend
    do n = 1, nnend
    do q = 1, qqend
        b1 = bladePoints(m,n,q,1)
        b2 = bladePoints(m,n,q,2)
        b3 = bladePoints(m,n,q,3)
        d1 = a1 - b1
        d2 = a2 - b2
        d3 = a3 - b3
        dist_sq = d1*d1 + d2*d2 + d3*d3

        if (dist_sq <= projectradius_sq) then
            kernel = exp(-dist_sq / epsilon_sq)
            f1 = f1 + bladeForces(m,n,q,1) * kernel
            f2 = f2 + bladeForces(m,n,q,2) * kernel
        end if
    end do
    end do
    end do

    f1 = f1 * const3
    f2 = f2 * const3

    if (has_nacelle == 1) then
        d1 = a1 - nacelle1
        d2 = a2 - nacelle2
        d3 = a3 - nacelle3
        dist_sq = d1*d1 + d2*d2 + d3*d3
        kernel = exp(-dist_sq / nacelle_epsilon_sq) * nacelle_kernel_norm
        f1 = f1 + nacelleForce1 * kernel * const2
        f2 = f2 + nacelleForce2 * kernel * const2
    endif

    bodyForceUV(1,c) = f1
    bodyForceUV(2,c) = f2
end do

!$cuf kernel do(1) <<<*,*>>>
do c = 1, cW
    a1 = locationW(1,c)
    a2 = locationW(2,c)
    a3 = locationW(3,c)
    f3 = 0._rprec

    do m = 1, mmend
    do n = 1, nnend
    do q = 1, qqend
        b1 = bladePoints(m,n,q,1)
        b2 = bladePoints(m,n,q,2)
        b3 = bladePoints(m,n,q,3)
        d1 = a1 - b1
        d2 = a2 - b2
        d3 = a3 - b3
        dist_sq = d1*d1 + d2*d2 + d3*d3

        if (dist_sq <= projectradius_sq) then
            kernel = exp(-dist_sq / epsilon_sq)
            f3 = f3 + bladeForces(m,n,q,3) * kernel
        end if
    end do
    end do
    end do

    f3 = f3 * const3

    if (has_nacelle == 1) then
        d1 = a1 - nacelle1
        d2 = a2 - nacelle2
        d3 = a3 - nacelle3
        dist_sq = d1*d1 + d2*d2 + d3*d3
        if (dist_sq <= projectradius_sq) then
            kernel = exp(-dist_sq / nacelle_epsilon_sq) * nacelle_kernel_norm
            f3 = f3 + nacelleForce3 * kernel * const2
        endif
    endif

    bodyForceW(3,c) = f3
end do

call atm_cuda_check('ATM force convolution')
#endif

end subroutine atm_lesgo_convolute_force_gpu_atpoint

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_force_gpu_atpoint(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched GPU equivalent of the active atPoint atm_computeBladeForce path.
implicit none

integer, intent(in) :: i

#ifdef ENABLE_CUDA
integer :: j, m, n, q, kk, mmend, nnend, qqend
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
integer :: sectionType_i, kair, rot_sign, direct_w_flag, mirror_force_flag
real(rprec) :: rotorApex1, rotorApex2, rotorApex3
real(rprec) :: uvShaft1, uvShaft2, uvShaft3
real(rprec) :: rotSpeed, pitch, pitchControlAngle, preCone, cosPreCone
real(rprec) :: z_local_min, z_local_max, degRadLocal
real(rprec) :: xyz1, xyz2, xyz3, px, py
real(rprec) :: xdiff, ydiff, zdiff
real(rprec) :: u1i, u2i, u3i, u4i, u5i, u6i
real(rprec) :: v1i, v2i, v3i, v4i, v5i, v6i
real(rprec) :: w1i, w2i, w3i, w4i, w5i, w6i, w7i, w8i
real(rprec) :: uLocal1, uLocal2, uLocal3
real(rprec) :: xvec1, xvec2, xvec3
real(rprec) :: yvec1, yvec2, yvec3
real(rprec) :: zvec1, zvec2, zvec3
real(rprec) :: dvec1, dvec2, dvec3
real(rprec) :: lvec1, lvec2, lvec3
real(rprec) :: mag, wind1, wind2, wind3, windAng, alpha_i
real(rprec) :: twistAng_i, chord_i, db_i, vmag_i, cl_i, cd_i
real(rprec) :: thrust_sum, torque_sum
real(rprec) :: xa, xb, ya, yb, lift_i, drag_i, sigma
real(rprec), parameter :: thresh = 1.e-9_rprec
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladePoints
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladeForces
real(rprec), managed, pointer, dimension(:,:,:,:) :: windVectors
real(rprec), managed, pointer, dimension(:,:,:,:,:) :: bladeAlignedVectors
real(rprec), managed, pointer, dimension(:,:,:) :: bladeRadius, solidity
real(rprec), managed, pointer, dimension(:,:,:) :: twistAng, chord
integer, managed, pointer, dimension(:,:,:) :: sectionType
real(rprec), managed, pointer, dimension(:,:,:) :: cl, cd, alpha, Vmag
real(rprec), managed, pointer, dimension(:,:,:) :: lift, drag
real(rprec), managed, pointer, dimension(:,:,:) :: axialForce, tangentialForce
real(rprec), managed, pointer, dimension(:,:,:) :: induction_a, u_infinity
real(rprec), managed, pointer, dimension(:,:,:,:) :: du
real(rprec), managed, pointer, dimension(:) :: db
type(airfoilType_t), managed, pointer, dimension(:) :: airfoil
logical :: tipCorrection

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints

rotorApex1 = turbineArray(i) % rotorApex(1)
rotorApex2 = turbineArray(i) % rotorApex(2)
rotorApex3 = turbineArray(i) % rotorApex(3)
uvShaft1 = turbineArray(i) % uvShaft(1)
uvShaft2 = turbineArray(i) % uvShaft(2)
uvShaft3 = turbineArray(i) % uvShaft(3)
rotSpeed = turbineArray(i) % rotSpeed
pitch = turbineArray(i) % Pitch
pitchControlAngle = turbineArray(i) % PitchControlAngle
preCone = turbineModel(j) % PreCone
cosPreCone = cos(preCone)
degRadLocal = pi / 180._rprec
z_local_min = grid % z(1)
z_local_max = grid % z(nz)
tipCorrection = turbineArray(i) % tipALMCorrection

if (turbineArray(i) % rotationDir == 'cw') then
    rot_sign = 1
else
    rot_sign = -1
endif

bladePoints => turbineArray(i) % bladePoints
bladeForces => turbineArray(i) % bladeForces
windVectors => turbineArray(i) % windVectors
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
bladeRadius => turbineArray(i) % bladeRadius
solidity => turbineArray(i) % solidity
twistAng => turbineArray(i) % twistAng
chord => turbineArray(i) % chord
sectionType => turbineArray(i) % sectionType
cl => turbineArray(i) % cl
cd => turbineArray(i) % cd
alpha => turbineArray(i) % alpha
Vmag => turbineArray(i) % Vmag
lift => turbineArray(i) % lift
drag => turbineArray(i) % drag
axialForce => turbineArray(i) % axialForce
tangentialForce => turbineArray(i) % tangentialForce
induction_a => turbineArray(i) % induction_a
u_infinity => turbineArray(i) % u_infinity
du => turbineArray(i) % du
db => turbineArray(i) % db
airfoil => turbineModel(j) % airfoilType
thrust_sum = 0._rprec
torque_sum = 0._rprec
direct_w_flag = 0
if (atm_direct_w_enabled()) direct_w_flag = 1
mirror_force_flag = 0
#ifdef PPLES_GPU
if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready) mirror_force_flag = 1
#endif

!$cuf kernel do(3) <<<*,*>>> reduction(+:thrust_sum, torque_sum)
do q = 1, qqend
do n = 1, nnend
do m = 1, mmend
    xyz1 = bladePoints(m,n,q,1) / z_i
    xyz2 = bladePoints(m,n,q,2) / z_i
    xyz3 = bladePoints(m,n,q,3) / z_i

    if (z_local_min <= xyz3 .and. xyz3 < z_local_max) then
        px = modulo(xyz1, L_x)
        if (abs(px) / L_x < thresh) then
            istart = 1
        elseif (abs(px - L_x) / L_x < thresh) then
            istart = nx
        else
            istart = floor(px / dx) + 1
        endif

        py = modulo(xyz2, L_y)
        if (abs(py) / L_y < thresh) then
            jstart = 1
        elseif (abs(py - L_y) / L_y < thresh) then
            jstart = ny
        else
            jstart = floor(py / dy) + 1
        endif

        if (abs(xyz3 - z_local_max) / L_z < thresh) then
            kstart = nz - 1
        else
            kstart = floor((xyz3 - z_local_min) / dz) + 1
        endif

        istart1 = istart + 1
        if (istart1 == nx + 1) istart1 = 1
        jstart1 = jstart + 1
        if (jstart1 == ny + 1) jstart1 = 1
        kstart1 = kstart + 1

        xdiff = xyz1 - real(istart - 1, rprec) * dx
        ydiff = xyz2 - real(jstart - 1, rprec) * dy
        zdiff = xyz3 - (z_local_min + real(kstart - 1, rprec) * dz)

        u1i = u(istart,jstart,kstart) + xdiff *                               &
            (u(istart1,jstart,kstart) - u(istart,jstart,kstart)) / dx
        u2i = u(istart,jstart1,kstart) + xdiff *                              &
            (u(istart1,jstart1,kstart) - u(istart,jstart1,kstart)) / dx
        u3i = u(istart,jstart,kstart1) + xdiff *                              &
            (u(istart1,jstart,kstart1) - u(istart,jstart,kstart1)) / dx
        u4i = u(istart,jstart1,kstart1) + xdiff *                             &
            (u(istart1,jstart1,kstart1) - u(istart,jstart1,kstart1)) / dx
        u5i = u1i + ydiff * (u2i - u1i) / dy
        u6i = u3i + ydiff * (u4i - u3i) / dy
        uLocal1 = (u5i + zdiff * (u6i - u5i) / dz) * u_star

        v1i = v(istart,jstart,kstart) + xdiff *                               &
            (v(istart1,jstart,kstart) - v(istart,jstart,kstart)) / dx
        v2i = v(istart,jstart1,kstart) + xdiff *                              &
            (v(istart1,jstart1,kstart) - v(istart,jstart1,kstart)) / dx
        v3i = v(istart,jstart,kstart1) + xdiff *                              &
            (v(istart1,jstart,kstart1) - v(istart,jstart,kstart1)) / dx
        v4i = v(istart,jstart1,kstart1) + xdiff *                             &
            (v(istart1,jstart1,kstart1) - v(istart,jstart1,kstart1)) / dx
        v5i = v1i + ydiff * (v2i - v1i) / dy
        v6i = v3i + ydiff * (v4i - v3i) / dy
        uLocal2 = (v5i + zdiff * (v6i - v5i) / dz) * u_star

        if (direct_w_flag == 1) then
            w1i = 0.5_rprec * (w(istart,jstart,kstart) +                      &
                w(istart,jstart,kstart+1))
            w2i = 0.5_rprec * (w(istart1,jstart,kstart) +                     &
                w(istart1,jstart,kstart+1))
            w3i = 0.5_rprec * (w(istart,jstart1,kstart) +                     &
                w(istart,jstart1,kstart+1))
            w4i = 0.5_rprec * (w(istart1,jstart1,kstart) +                    &
                w(istart1,jstart1,kstart+1))
            w5i = w1i + xdiff * (w2i - w1i) / dx
            w6i = w3i + xdiff * (w4i - w3i) / dx
            w7i = w5i + ydiff * (w6i - w5i) / dy

            if (kstart1 == nz .and. coord < nproc - 1) then
                w1i = atm_wuv_recv_up(istart,jstart)
                w2i = atm_wuv_recv_up(istart1,jstart)
                w3i = atm_wuv_recv_up(istart,jstart1)
                w4i = atm_wuv_recv_up(istart1,jstart1)
            elseif (kstart1 == nz) then
                w1i = 0.5_rprec * (w(istart,jstart,nz-1) + w(istart,jstart,nz))
                w2i = 0.5_rprec * (w(istart1,jstart,nz-1) +                  &
                    w(istart1,jstart,nz))
                w3i = 0.5_rprec * (w(istart,jstart1,nz-1) +                  &
                    w(istart,jstart1,nz))
                w4i = 0.5_rprec * (w(istart1,jstart1,nz-1) +                 &
                    w(istart1,jstart1,nz))
            else
                w1i = 0.5_rprec * (w(istart,jstart,kstart1) +                 &
                    w(istart,jstart,kstart1+1))
                w2i = 0.5_rprec * (w(istart1,jstart,kstart1) +                &
                    w(istart1,jstart,kstart1+1))
                w3i = 0.5_rprec * (w(istart,jstart1,kstart1) +                &
                    w(istart,jstart1,kstart1+1))
                w4i = 0.5_rprec * (w(istart1,jstart1,kstart1) +               &
                    w(istart1,jstart1,kstart1+1))
            endif
            w5i = w1i + xdiff * (w2i - w1i) / dx
            w6i = w3i + xdiff * (w4i - w3i) / dx
            w8i = w5i + ydiff * (w6i - w5i) / dy
            uLocal3 = (w7i + zdiff * (w8i - w7i) / dz) * u_star
        else
            w1i = w_uv(istart,jstart,kstart) + xdiff *                        &
                (w_uv(istart1,jstart,kstart) - w_uv(istart,jstart,kstart)) / dx
            w2i = w_uv(istart,jstart1,kstart) + xdiff *                       &
                (w_uv(istart1,jstart1,kstart) - w_uv(istart,jstart1,kstart)) / dx
            w3i = w_uv(istart,jstart,kstart1) + xdiff *                       &
                (w_uv(istart1,jstart,kstart1) - w_uv(istart,jstart,kstart1)) / dx
            w4i = w_uv(istart,jstart1,kstart1) + xdiff *                      &
                (w_uv(istart1,jstart1,kstart1) - w_uv(istart,jstart1,kstart1)) / dx
            w5i = w1i + ydiff * (w2i - w1i) / dy
            w6i = w3i + ydiff * (w4i - w3i) / dy
            uLocal3 = (w5i + zdiff * (w6i - w5i) / dz) * u_star
        endif

        if (rot_sign == 1) then
            zvec1 = bladePoints(m,n,q,1) - rotorApex1
            zvec2 = bladePoints(m,n,q,2) - rotorApex2
            zvec3 = bladePoints(m,n,q,3) - rotorApex3
        else
            zvec1 = rotorApex1 - bladePoints(m,n,q,1)
            zvec2 = rotorApex2 - bladePoints(m,n,q,2)
            zvec3 = rotorApex3 - bladePoints(m,n,q,3)
        endif

        mag = sqrt(zvec1*zvec1 + zvec2*zvec2 + zvec3*zvec3)
        zvec1 = zvec1 / mag
        zvec2 = zvec2 / mag
        zvec3 = zvec3 / mag

        yvec1 = zvec2 * uvShaft3 - zvec3 * uvShaft2
        yvec2 = zvec3 * uvShaft1 - zvec1 * uvShaft3
        yvec3 = zvec1 * uvShaft2 - zvec2 * uvShaft1
        mag = sqrt(yvec1*yvec1 + yvec2*yvec2 + yvec3*yvec3)
        yvec1 = yvec1 / mag
        yvec2 = yvec2 / mag
        yvec3 = yvec3 / mag

        xvec1 = yvec2 * zvec3 - yvec3 * zvec2
        xvec2 = yvec3 * zvec1 - yvec1 * zvec3
        xvec3 = yvec1 * zvec2 - yvec2 * zvec1
        mag = sqrt(xvec1*xvec1 + xvec2*xvec2 + xvec3*xvec3)
        xvec1 = xvec1 / mag
        xvec2 = xvec2 / mag
        xvec3 = xvec3 / mag

        bladeAlignedVectors(m,n,q,1,1) = xvec1
        bladeAlignedVectors(m,n,q,1,2) = xvec2
        bladeAlignedVectors(m,n,q,1,3) = xvec3
        bladeAlignedVectors(m,n,q,2,1) = yvec1
        bladeAlignedVectors(m,n,q,2,2) = yvec2
        bladeAlignedVectors(m,n,q,2,3) = yvec3
        bladeAlignedVectors(m,n,q,3,1) = zvec1
        bladeAlignedVectors(m,n,q,3,2) = zvec2
        bladeAlignedVectors(m,n,q,3,3) = zvec3

        wind1 = xvec1*uLocal1 + xvec2*uLocal2 + xvec3*uLocal3
        wind2 = yvec1*uLocal1 + yvec2*uLocal2 + yvec3*uLocal3 +               &
            rotSpeed * bladeRadius(m,n,q) * cosPreCone
        wind3 = zvec1*uLocal1 + zvec2*uLocal2 + zvec3*uLocal3

        windVectors(m,n,q,1) = wind1
        windVectors(m,n,q,2) = wind2
        windVectors(m,n,q,3) = wind3

        if (tipCorrection) then
            wind1 = wind1 + du(m,n,q,1)
            wind2 = wind2 + du(m,n,q,2)
            windVectors(m,n,q,1) = wind1
            windVectors(m,n,q,2) = wind2
        endif

        twistAng_i = twistAng(m,n,q)
        chord_i = chord(m,n,q)
        db_i = db(q)
        sectionType_i = sectionType(m,n,q)
        vmag_i = sqrt(wind1*wind1 + wind2*wind2)
        Vmag(m,n,q) = vmag_i

        windAng = atan2(wind1, wind2) / degRadLocal
        alpha_i = windAng - twistAng_i - pitch - pitchControlAngle
        alpha(m,n,q) = alpha_i

        kair = airfoil(sectionType_i) % n

        cl_i = 0._rprec
        if (alpha_i <= airfoil(sectionType_i) % AOA(1)) then
            cl_i = airfoil(sectionType_i) % Cl(1)
        elseif (alpha_i >= airfoil(sectionType_i) % AOA(kair)) then
            cl_i = airfoil(sectionType_i) % Cl(kair)
        else
            do kk = 2, kair
                if (alpha_i >= airfoil(sectionType_i) % AOA(kk-1) .and.       &
                    alpha_i <= airfoil(sectionType_i) % AOA(kk)) then
                    xa = airfoil(sectionType_i) % AOA(kk-1)
                    xb = airfoil(sectionType_i) % AOA(kk)
                    ya = airfoil(sectionType_i) % Cl(kk-1)
                    yb = airfoil(sectionType_i) % Cl(kk)
                    cl_i = ya + (yb - ya) * (alpha_i - xa) / (xb - xa)
                endif
            end do
        endif

        cd_i = 0._rprec
        if (alpha_i <= airfoil(sectionType_i) % AOA(1)) then
            cd_i = airfoil(sectionType_i) % Cd(1)
        elseif (alpha_i >= airfoil(sectionType_i) % AOA(kair)) then
            cd_i = airfoil(sectionType_i) % Cd(kair)
        else
            do kk = 2, kair
                if (alpha_i >= airfoil(sectionType_i) % AOA(kk-1) .and.       &
                    alpha_i <= airfoil(sectionType_i) % AOA(kk)) then
                    xa = airfoil(sectionType_i) % AOA(kk-1)
                    xb = airfoil(sectionType_i) % AOA(kk)
                    ya = airfoil(sectionType_i) % Cd(kk-1)
                    yb = airfoil(sectionType_i) % Cd(kk)
                    cd_i = ya + (yb - ya) * (alpha_i - xa) / (xb - xa)
                endif
            end do
        endif

        cl(m,n,q) = cl_i
        cd(m,n,q) = cd_i

        lift_i = 0.5_rprec * cl_i * (vmag_i*vmag_i) * chord_i * db_i *        &
            solidity(m,n,q)
        drag_i = 0.5_rprec * cd_i * (vmag_i*vmag_i) * chord_i * db_i *        &
            solidity(m,n,q)
        lift(m,n,q) = lift_i
        drag(m,n,q) = drag_i

        dvec1 = xvec1*wind1 + yvec1*wind2
        dvec2 = xvec2*wind1 + yvec2*wind2
        dvec3 = xvec3*wind1 + yvec3*wind2
        mag = sqrt(dvec1*dvec1 + dvec2*dvec2 + dvec3*dvec3)
        dvec1 = dvec1 / mag
        dvec2 = dvec2 / mag
        dvec3 = dvec3 / mag

        lvec1 = dvec2*zvec3 - dvec3*zvec2
        lvec2 = dvec3*zvec1 - dvec1*zvec3
        lvec3 = dvec1*zvec2 - dvec2*zvec1
        mag = sqrt(lvec1*lvec1 + lvec2*lvec2 + lvec3*lvec3)
        lvec1 = lvec1 / mag
        lvec2 = lvec2 / mag
        lvec3 = lvec3 / mag

        lvec1 = -lift_i * lvec1
        lvec2 = -lift_i * lvec2
        lvec3 = -lift_i * lvec3
        dvec1 = -drag_i * dvec1
        dvec2 = -drag_i * dvec2
        dvec3 = -drag_i * dvec3

        bladeForces(m,n,q,1) = lvec1 + dvec1
        bladeForces(m,n,q,2) = lvec2 + dvec2
        bladeForces(m,n,q,3) = lvec3 + dvec3
#ifdef PPLES_GPU
        if (mirror_force_flag == 1) then
            atm_bladeForces_d(m,n,q,1,i) = bladeForces(m,n,q,1)
            atm_bladeForces_d(m,n,q,2,i) = bladeForces(m,n,q,2)
            atm_bladeForces_d(m,n,q,3,i) = bladeForces(m,n,q,3)
        endif
#endif

        axialForce(m,n,q) = -(bladeForces(m,n,q,1)*uvShaft1 +                 &
            bladeForces(m,n,q,2)*uvShaft2 + bladeForces(m,n,q,3)*uvShaft3)
        tangentialForce(m,n,q) = bladeForces(m,n,q,1)*yvec1 +                 &
            bladeForces(m,n,q,2)*yvec2 + bladeForces(m,n,q,3)*yvec3
        thrust_sum = thrust_sum + axialForce(m,n,q)
        torque_sum = torque_sum + tangentialForce(m,n,q)*bladeRadius(m,n,q)

        windAng = windAng * degRadLocal
        sigma = chord_i * real(mmend, rprec) /                                &
            (2._rprec*pi*bladeRadius(m,n,q))
        induction_a(m,n,q) = 1._rprec /                                       &
            (4._rprec * sin(windAng)**2 /                                     &
            (sigma * (cl_i*cos(windAng) + cd_i*sin(windAng))) + 1._rprec)
        u_infinity(m,n,q) = wind1
    else
        bladeForces(m,n,q,1) = 0._rprec
        bladeForces(m,n,q,2) = 0._rprec
        bladeForces(m,n,q,3) = 0._rprec
#ifdef PPLES_GPU
        if (mirror_force_flag == 1) then
            atm_bladeForces_d(m,n,q,1,i) = 0._rprec
            atm_bladeForces_d(m,n,q,2,i) = 0._rprec
            atm_bladeForces_d(m,n,q,3,i) = 0._rprec
        endif
#endif
        bladeAlignedVectors(m,n,q,1,1) = 0._rprec
        bladeAlignedVectors(m,n,q,1,2) = 0._rprec
        bladeAlignedVectors(m,n,q,1,3) = 0._rprec
        bladeAlignedVectors(m,n,q,2,1) = 0._rprec
        bladeAlignedVectors(m,n,q,2,2) = 0._rprec
        bladeAlignedVectors(m,n,q,2,3) = 0._rprec
        bladeAlignedVectors(m,n,q,3,1) = 0._rprec
        bladeAlignedVectors(m,n,q,3,2) = 0._rprec
        bladeAlignedVectors(m,n,q,3,3) = 0._rprec
        windVectors(m,n,q,1) = 0._rprec
        windVectors(m,n,q,2) = 0._rprec
        windVectors(m,n,q,3) = 0._rprec
        alpha(m,n,q) = 0._rprec
        Vmag(m,n,q) = 0._rprec
        cl(m,n,q) = 0._rprec
        cd(m,n,q) = 0._rprec
        lift(m,n,q) = 0._rprec
        drag(m,n,q) = 0._rprec
        axialForce(m,n,q) = 0._rprec
        tangentialForce(m,n,q) = 0._rprec
        induction_a(m,n,q) = 0._rprec
        u_infinity(m,n,q) = 0._rprec
    endif
end do
end do
end do

call atm_cuda_sync('ATM blade force')

turbineArray(i) % thrust = thrust_sum
turbineArray(i) % torqueRotor = torque_sum * cosPreCone
#endif

end subroutine atm_lesgo_force_gpu_atpoint

#ifdef ENABLE_CUDA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_nacelle_force_gpu_atpoint(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, c, direct_w_flag
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
real(rprec) :: xyz1, xyz2, xyz3, px, py
real(rprec) :: xdiff, ydiff, zdiff
real(rprec) :: u1i, u2i, u3i, u4i, u5i, u6i
real(rprec) :: v1i, v2i, v3i, v4i, v5i, v6i
real(rprec) :: w1i, w2i, w3i, w4i, w5i, w6i, w7i, w8i
real(rprec) :: uLocal1, uLocal2, uLocal3
real(rprec) :: uv1, uv2, uv3, area, drag, vel, vel_corr
real(rprec) :: nacelleCd, nacelleEpsilon, z_local_min, z_local_max
real(rprec), parameter :: thresh = 1.e-9_rprec
#ifdef PPLES_GPU
real(rprec) :: nacelle_tmp_h(5)
#endif

j = turbineArray(i) % turbineTypeID
xyz1 = turbineArray(i) % nacelleLocation(1) / z_i
xyz2 = turbineArray(i) % nacelleLocation(2) / z_i
xyz3 = turbineArray(i) % nacelleLocation(3) / z_i
uv1 = turbineArray(i) % uvShaft(1)
uv2 = turbineArray(i) % uvShaft(2)
uv3 = turbineArray(i) % uvShaft(3)
nacelleCd = turbineArray(i) % nacelleCd
nacelleEpsilon = turbineArray(i) % nacelleEpsilon
area = pi * turbineModel(j) % hubRad * turbineModel(j) % hubRad
z_local_min = grid % z(1)
z_local_max = grid % z(nz)
direct_w_flag = 0
if (atm_direct_w_enabled()) direct_w_flag = 1

if (.not. allocated(atm_nacelle_tmp)) allocate(atm_nacelle_tmp(5))
atm_nacelle_tmp = 0._rprec

if (nacelleEpsilon <= 0._rprec) then
    turbineArray(i) % nacelleForce = 0._rprec
    turbineArray(i) % VelNacelle_sampled = 0._rprec
    turbineArray(i) % VelNacelle_corrected = 0._rprec
    return
endif

!$cuf kernel do(1) <<<*,*>>>
do c = 1, 1
    if (z_local_min <= xyz3 .and. xyz3 < z_local_max) then
        px = modulo(xyz1, L_x)
        if (abs(px) / L_x < thresh) then
            istart = 1
        elseif (abs(px - L_x) / L_x < thresh) then
            istart = nx
        else
            istart = floor(px / dx) + 1
        endif

        py = modulo(xyz2, L_y)
        if (abs(py) / L_y < thresh) then
            jstart = 1
        elseif (abs(py - L_y) / L_y < thresh) then
            jstart = ny
        else
            jstart = floor(py / dy) + 1
        endif

        if (abs(xyz3 - z_local_max) / L_z < thresh) then
            kstart = nz - 1
        else
            kstart = floor((xyz3 - z_local_min) / dz) + 1
        endif

        istart1 = istart + 1
        if (istart1 == nx + 1) istart1 = 1
        jstart1 = jstart + 1
        if (jstart1 == ny + 1) jstart1 = 1
        kstart1 = kstart + 1

        xdiff = xyz1 - real(istart - 1, rprec) * dx
        ydiff = xyz2 - real(jstart - 1, rprec) * dy
        zdiff = xyz3 - (z_local_min + real(kstart - 1, rprec) * dz)

        u1i = u(istart,jstart,kstart) + xdiff *                               &
            (u(istart1,jstart,kstart) - u(istart,jstart,kstart)) / dx
        u2i = u(istart,jstart1,kstart) + xdiff *                              &
            (u(istart1,jstart1,kstart) - u(istart,jstart1,kstart)) / dx
        u3i = u(istart,jstart,kstart1) + xdiff *                              &
            (u(istart1,jstart,kstart1) - u(istart,jstart,kstart1)) / dx
        u4i = u(istart,jstart1,kstart1) + xdiff *                             &
            (u(istart1,jstart1,kstart1) - u(istart,jstart1,kstart1)) / dx
        u5i = u1i + ydiff * (u2i - u1i) / dy
        u6i = u3i + ydiff * (u4i - u3i) / dy
        uLocal1 = (u5i + zdiff * (u6i - u5i) / dz) * u_star

        v1i = v(istart,jstart,kstart) + xdiff *                               &
            (v(istart1,jstart,kstart) - v(istart,jstart,kstart)) / dx
        v2i = v(istart,jstart1,kstart) + xdiff *                              &
            (v(istart1,jstart1,kstart) - v(istart,jstart1,kstart)) / dx
        v3i = v(istart,jstart,kstart1) + xdiff *                              &
            (v(istart1,jstart,kstart1) - v(istart,jstart,kstart1)) / dx
        v4i = v(istart,jstart1,kstart1) + xdiff *                             &
            (v(istart1,jstart1,kstart1) - v(istart,jstart1,kstart1)) / dx
        v5i = v1i + ydiff * (v2i - v1i) / dy
        v6i = v3i + ydiff * (v4i - v3i) / dy
        uLocal2 = (v5i + zdiff * (v6i - v5i) / dz) * u_star

        if (direct_w_flag == 1) then
            w1i = 0.5_rprec * (w(istart,jstart,kstart) +                      &
                w(istart,jstart,kstart+1))
            w2i = 0.5_rprec * (w(istart1,jstart,kstart) +                     &
                w(istart1,jstart,kstart+1))
            w3i = 0.5_rprec * (w(istart,jstart1,kstart) +                     &
                w(istart,jstart1,kstart+1))
            w4i = 0.5_rprec * (w(istart1,jstart1,kstart) +                    &
                w(istart1,jstart1,kstart+1))
            w5i = w1i + xdiff * (w2i - w1i) / dx
            w6i = w3i + xdiff * (w4i - w3i) / dx
            w7i = w5i + ydiff * (w6i - w5i) / dy

            if (kstart1 == nz .and. coord < nproc - 1) then
                w1i = atm_wuv_recv_up(istart,jstart)
                w2i = atm_wuv_recv_up(istart1,jstart)
                w3i = atm_wuv_recv_up(istart,jstart1)
                w4i = atm_wuv_recv_up(istart1,jstart1)
            elseif (kstart1 == nz) then
                w1i = 0.5_rprec * (w(istart,jstart,nz-1) + w(istart,jstart,nz))
                w2i = 0.5_rprec * (w(istart1,jstart,nz-1) +                  &
                    w(istart1,jstart,nz))
                w3i = 0.5_rprec * (w(istart,jstart1,nz-1) +                  &
                    w(istart,jstart1,nz))
                w4i = 0.5_rprec * (w(istart1,jstart1,nz-1) +                 &
                    w(istart1,jstart1,nz))
            else
                w1i = 0.5_rprec * (w(istart,jstart,kstart1) +                 &
                    w(istart,jstart,kstart1+1))
                w2i = 0.5_rprec * (w(istart1,jstart,kstart1) +                &
                    w(istart1,jstart,kstart1+1))
                w3i = 0.5_rprec * (w(istart,jstart1,kstart1) +                &
                    w(istart,jstart1,kstart1+1))
                w4i = 0.5_rprec * (w(istart1,jstart1,kstart1) +               &
                    w(istart1,jstart1,kstart1+1))
            endif
            w5i = w1i + xdiff * (w2i - w1i) / dx
            w6i = w3i + xdiff * (w4i - w3i) / dx
            w8i = w5i + ydiff * (w6i - w5i) / dy
            uLocal3 = (w7i + zdiff * (w8i - w7i) / dz) * u_star
        else
            w1i = w_uv(istart,jstart,kstart) + xdiff *                        &
                (w_uv(istart1,jstart,kstart) - w_uv(istart,jstart,kstart)) / dx
            w2i = w_uv(istart,jstart1,kstart) + xdiff *                       &
                (w_uv(istart1,jstart1,kstart) - w_uv(istart,jstart1,kstart)) / dx
            w3i = w_uv(istart,jstart,kstart1) + xdiff *                       &
                (w_uv(istart1,jstart,kstart1) - w_uv(istart,jstart,kstart1)) / dx
            w4i = w_uv(istart,jstart1,kstart1) + xdiff *                      &
                (w_uv(istart1,jstart1,kstart1) - w_uv(istart,jstart1,kstart1)) / dx
            w5i = w1i + ydiff * (w2i - w1i) / dy
            w6i = w3i + ydiff * (w4i - w3i) / dy
            uLocal3 = (w5i + zdiff * (w6i - w5i) / dz) * u_star
        endif

        vel = uv1*uLocal1 + uv2*uLocal2 + uv3*uLocal3
        vel_corr = vel / (1._rprec - 0.25_rprec / pi * nacelleCd * area /     &
            (nacelleEpsilon*nacelleEpsilon))

        atm_nacelle_tmp(4) = vel
        atm_nacelle_tmp(5) = vel_corr
        if (vel_corr >= 0._rprec) then
            drag = 0.5_rprec * nacelleCd * vel_corr * vel_corr * area
            atm_nacelle_tmp(1) = -drag * uv1
            atm_nacelle_tmp(2) = -drag * uv2
            atm_nacelle_tmp(3) = -drag * uv3
        endif
    endif
enddo

call atm_cuda_sync('ATM nacelle force')

#ifdef PPLES_GPU
nacelle_tmp_h = atm_nacelle_tmp
turbineArray(i) % nacelleForce(1) = nacelle_tmp_h(1)
turbineArray(i) % nacelleForce(2) = nacelle_tmp_h(2)
turbineArray(i) % nacelleForce(3) = nacelle_tmp_h(3)
turbineArray(i) % VelNacelle_sampled = nacelle_tmp_h(4)
turbineArray(i) % VelNacelle_corrected = nacelle_tmp_h(5)
#else
turbineArray(i) % nacelleForce(1) = atm_nacelle_tmp(1)
turbineArray(i) % nacelleForce(2) = atm_nacelle_tmp(2)
turbineArray(i) % nacelleForce(3) = atm_nacelle_tmp(3)
turbineArray(i) % VelNacelle_sampled = atm_nacelle_tmp(4)
turbineArray(i) % VelNacelle_corrected = atm_nacelle_tmp(5)
#endif

end subroutine atm_lesgo_nacelle_force_gpu_atpoint

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_point_owner_lb_supported()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i

atm_point_owner_lb_supported = .true.
do i = 1, numberOfTurbines
    if (trim(turbineArray(i) % sampling) /= 'atPoint') then
        atm_point_owner_lb_supported = .false.
    endif
enddo

end function atm_point_owner_lb_supported

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_point_owner_sample_turbine(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, m, n, q, pid, mmend, nnend, qqend
integer :: istart, jstart, kstart, istart1, jstart1, kstart1, direct_w_flag
real(rprec) :: z_local_min, z_local_max, xyz1, xyz2, xyz3, px, py
real(rprec) :: xdiff, ydiff, zdiff
real(rprec) :: u1i, u2i, u3i, u4i, u5i, u6i
real(rprec) :: v1i, v2i, v3i, v4i, v5i, v6i
real(rprec) :: w1i, w2i, w3i, w4i, w5i, w6i, w7i, w8i
real(rprec), parameter :: thresh = 1.e-9_rprec
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladePoints

if (.not. atm_lb_owner_ready) return

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
z_local_min = grid % z(1)
z_local_max = grid % z(nz)
bladePoints => turbineArray(i) % bladePoints
direct_w_flag = 0
if (atm_direct_w_enabled()) direct_w_flag = 1

!$cuf kernel do(3) <<<*,*>>>
do q = 1, qqend
do n = 1, nnend
do m = 1, mmend
    pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m
    xyz1 = bladePoints(m,n,q,1) / z_i
    xyz2 = bladePoints(m,n,q,2) / z_i
    xyz3 = bladePoints(m,n,q,3) / z_i

    if (z_local_min <= xyz3 .and. xyz3 < z_local_max) then
        px = modulo(xyz1, L_x)
        if (abs(px) / L_x < thresh) then
            istart = 1
        elseif (abs(px - L_x) / L_x < thresh) then
            istart = nx
        else
            istart = floor(px / dx) + 1
        endif

        py = modulo(xyz2, L_y)
        if (abs(py) / L_y < thresh) then
            jstart = 1
        elseif (abs(py - L_y) / L_y < thresh) then
            jstart = ny
        else
            jstart = floor(py / dy) + 1
        endif

        if (abs(xyz3 - z_local_max) / L_z < thresh) then
            kstart = nz - 1
        else
            kstart = floor((xyz3 - z_local_min) / dz) + 1
        endif

        istart1 = istart + 1
        if (istart1 == nx + 1) istart1 = 1
        jstart1 = jstart + 1
        if (jstart1 == ny + 1) jstart1 = 1
        kstart1 = kstart + 1

        xdiff = xyz1 - real(istart - 1, rprec) * dx
        ydiff = xyz2 - real(jstart - 1, rprec) * dy
        zdiff = xyz3 - (z_local_min + real(kstart - 1, rprec) * dz)

        u1i = u(istart,jstart,kstart) + xdiff *                              &
            (u(istart1,jstart,kstart) - u(istart,jstart,kstart)) / dx
        u2i = u(istart,jstart1,kstart) + xdiff *                             &
            (u(istart1,jstart1,kstart) - u(istart,jstart1,kstart)) / dx
        u3i = u(istart,jstart,kstart1) + xdiff *                             &
            (u(istart1,jstart,kstart1) - u(istart,jstart,kstart1)) / dx
        u4i = u(istart,jstart1,kstart1) + xdiff *                            &
            (u(istart1,jstart1,kstart1) - u(istart,jstart1,kstart1)) / dx
        u5i = u1i + ydiff * (u2i - u1i) / dy
        u6i = u3i + ydiff * (u4i - u3i) / dy
        atm_lb_vel_local(1,pid) = (u5i + zdiff * (u6i - u5i) / dz) * u_star

        v1i = v(istart,jstart,kstart) + xdiff *                              &
            (v(istart1,jstart,kstart) - v(istart,jstart,kstart)) / dx
        v2i = v(istart,jstart1,kstart) + xdiff *                             &
            (v(istart1,jstart1,kstart) - v(istart,jstart1,kstart)) / dx
        v3i = v(istart,jstart,kstart1) + xdiff *                             &
            (v(istart1,jstart,kstart1) - v(istart,jstart,kstart1)) / dx
        v4i = v(istart,jstart1,kstart1) + xdiff *                            &
            (v(istart1,jstart1,kstart1) - v(istart,jstart1,kstart1)) / dx
        v5i = v1i + ydiff * (v2i - v1i) / dy
        v6i = v3i + ydiff * (v4i - v3i) / dy
        atm_lb_vel_local(2,pid) = (v5i + zdiff * (v6i - v5i) / dz) * u_star

        if (direct_w_flag == 1) then
            w1i = 0.5_rprec * (w(istart,jstart,kstart) +                     &
                w(istart,jstart,kstart+1))
            w2i = 0.5_rprec * (w(istart1,jstart,kstart) +                    &
                w(istart1,jstart,kstart+1))
            w3i = 0.5_rprec * (w(istart,jstart1,kstart) +                    &
                w(istart,jstart1,kstart+1))
            w4i = 0.5_rprec * (w(istart1,jstart1,kstart) +                   &
                w(istart1,jstart1,kstart+1))
            w5i = w1i + xdiff * (w2i - w1i) / dx
            w6i = w3i + xdiff * (w4i - w3i) / dx
            w7i = w5i + ydiff * (w6i - w5i) / dy

            if (kstart1 == nz .and. coord < nproc - 1) then
                w1i = atm_wuv_recv_up(istart,jstart)
                w2i = atm_wuv_recv_up(istart1,jstart)
                w3i = atm_wuv_recv_up(istart,jstart1)
                w4i = atm_wuv_recv_up(istart1,jstart1)
            elseif (kstart1 == nz) then
                w1i = 0.5_rprec * (w(istart,jstart,nz-1) +                   &
                    w(istart,jstart,nz))
                w2i = 0.5_rprec * (w(istart1,jstart,nz-1) +                  &
                    w(istart1,jstart,nz))
                w3i = 0.5_rprec * (w(istart,jstart1,nz-1) +                  &
                    w(istart,jstart1,nz))
                w4i = 0.5_rprec * (w(istart1,jstart1,nz-1) +                 &
                    w(istart1,jstart1,nz))
            else
                w1i = 0.5_rprec * (w(istart,jstart,kstart1) +                &
                    w(istart,jstart,kstart1+1))
                w2i = 0.5_rprec * (w(istart1,jstart,kstart1) +               &
                    w(istart1,jstart,kstart1+1))
                w3i = 0.5_rprec * (w(istart,jstart1,kstart1) +               &
                    w(istart,jstart1,kstart1+1))
                w4i = 0.5_rprec * (w(istart1,jstart1,kstart1) +              &
                    w(istart1,jstart1,kstart1+1))
            endif
            w5i = w1i + xdiff * (w2i - w1i) / dx
            w6i = w3i + xdiff * (w4i - w3i) / dx
            w8i = w5i + ydiff * (w6i - w5i) / dy
            atm_lb_vel_local(3,pid) = (w7i + zdiff * (w8i - w7i) / dz) *      &
                u_star
        else
            w1i = w_uv(istart,jstart,kstart) + xdiff *                       &
                (w_uv(istart1,jstart,kstart) -                               &
                w_uv(istart,jstart,kstart)) / dx
            w2i = w_uv(istart,jstart1,kstart) + xdiff *                      &
                (w_uv(istart1,jstart1,kstart) -                              &
                w_uv(istart,jstart1,kstart)) / dx
            w3i = w_uv(istart,jstart,kstart1) + xdiff *                      &
                (w_uv(istart1,jstart,kstart1) -                              &
                w_uv(istart,jstart,kstart1)) / dx
            w4i = w_uv(istart,jstart1,kstart1) + xdiff *                     &
                (w_uv(istart1,jstart1,kstart1) -                             &
                w_uv(istart,jstart1,kstart1)) / dx
            w5i = w1i + ydiff * (w2i - w1i) / dy
            w6i = w3i + ydiff * (w4i - w3i) / dy
            atm_lb_vel_local(3,pid) = (w5i + zdiff * (w6i - w5i) / dz) *      &
                u_star
        endif
    endif
end do
end do
end do

end subroutine atm_point_owner_sample_turbine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_point_owner_force_turbine(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, m, n, q, kk, pid, mmend, nnend, qqend
integer :: sectionType_i, kair, rot_sign
real(rprec) :: rotorApex1, rotorApex2, rotorApex3
real(rprec) :: uvShaft1, uvShaft2, uvShaft3
real(rprec) :: rotSpeed, pitch, pitchControlAngle, preCone, cosPreCone
real(rprec) :: degRadLocal
real(rprec) :: uLocal1, uLocal2, uLocal3
real(rprec) :: xvec1, xvec2, xvec3, yvec1, yvec2, yvec3, zvec1, zvec2, zvec3
real(rprec) :: dvec1, dvec2, dvec3, lvec1, lvec2, lvec3
real(rprec) :: mag, wind1, wind2, wind3, windAng, alpha_i
real(rprec) :: twistAng_i, chord_i, db_i, vmag_i, cl_i, cd_i
real(rprec) :: thrust_sum, torque_sum
real(rprec) :: xa, xb, ya, yb, lift_i, drag_i, sigma
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladePoints
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladeForces
real(rprec), managed, pointer, dimension(:,:,:,:) :: windVectors
real(rprec), managed, pointer, dimension(:,:,:,:,:) :: bladeAlignedVectors
real(rprec), managed, pointer, dimension(:,:,:) :: bladeRadius, solidity
real(rprec), managed, pointer, dimension(:,:,:) :: twistAng, chord
integer, managed, pointer, dimension(:,:,:) :: sectionType
real(rprec), managed, pointer, dimension(:,:,:) :: cl, cd, alpha, Vmag
real(rprec), managed, pointer, dimension(:,:,:) :: lift, drag
real(rprec), managed, pointer, dimension(:,:,:) :: axialForce, tangentialForce
real(rprec), managed, pointer, dimension(:,:,:) :: induction_a, u_infinity
real(rprec), managed, pointer, dimension(:,:,:,:) :: du
real(rprec), managed, pointer, dimension(:) :: db
type(airfoilType_t), managed, pointer, dimension(:) :: airfoil
logical :: tipCorrection

if (.not. atm_lb_owner_ready) return

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
rotorApex1 = turbineArray(i) % rotorApex(1)
rotorApex2 = turbineArray(i) % rotorApex(2)
rotorApex3 = turbineArray(i) % rotorApex(3)
uvShaft1 = turbineArray(i) % uvShaft(1)
uvShaft2 = turbineArray(i) % uvShaft(2)
uvShaft3 = turbineArray(i) % uvShaft(3)
rotSpeed = turbineArray(i) % rotSpeed
pitch = turbineArray(i) % Pitch
pitchControlAngle = turbineArray(i) % PitchControlAngle
preCone = turbineModel(j) % PreCone
cosPreCone = cos(preCone)
degRadLocal = pi / 180._rprec
tipCorrection = turbineArray(i) % tipALMCorrection

if (turbineArray(i) % rotationDir == 'cw') then
    rot_sign = 1
else
    rot_sign = -1
endif

bladePoints => turbineArray(i) % bladePoints
bladeForces => turbineArray(i) % bladeForces
windVectors => turbineArray(i) % windVectors
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
bladeRadius => turbineArray(i) % bladeRadius
solidity => turbineArray(i) % solidity
twistAng => turbineArray(i) % twistAng
chord => turbineArray(i) % chord
sectionType => turbineArray(i) % sectionType
cl => turbineArray(i) % cl
cd => turbineArray(i) % cd
alpha => turbineArray(i) % alpha
Vmag => turbineArray(i) % Vmag
lift => turbineArray(i) % lift
drag => turbineArray(i) % drag
axialForce => turbineArray(i) % axialForce
tangentialForce => turbineArray(i) % tangentialForce
induction_a => turbineArray(i) % induction_a
u_infinity => turbineArray(i) % u_infinity
du => turbineArray(i) % du
db => turbineArray(i) % db
airfoil => turbineModel(j) % airfoilType
thrust_sum = 0._rprec
torque_sum = 0._rprec

!$cuf kernel do(3) <<<*,*>>> reduction(+:thrust_sum, torque_sum)
do q = 1, qqend
do n = 1, nnend
do m = 1, mmend
    pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m

    if (atm_lb_owner(pid) == coord) then
        uLocal1 = atm_lb_vel_global(1,pid)
        uLocal2 = atm_lb_vel_global(2,pid)
        uLocal3 = atm_lb_vel_global(3,pid)

        if (rot_sign == 1) then
            zvec1 = bladePoints(m,n,q,1) - rotorApex1
            zvec2 = bladePoints(m,n,q,2) - rotorApex2
            zvec3 = bladePoints(m,n,q,3) - rotorApex3
        else
            zvec1 = rotorApex1 - bladePoints(m,n,q,1)
            zvec2 = rotorApex2 - bladePoints(m,n,q,2)
            zvec3 = rotorApex3 - bladePoints(m,n,q,3)
        endif

        mag = sqrt(zvec1*zvec1 + zvec2*zvec2 + zvec3*zvec3)
        zvec1 = zvec1 / mag
        zvec2 = zvec2 / mag
        zvec3 = zvec3 / mag

        yvec1 = zvec2 * uvShaft3 - zvec3 * uvShaft2
        yvec2 = zvec3 * uvShaft1 - zvec1 * uvShaft3
        yvec3 = zvec1 * uvShaft2 - zvec2 * uvShaft1
        mag = sqrt(yvec1*yvec1 + yvec2*yvec2 + yvec3*yvec3)
        yvec1 = yvec1 / mag
        yvec2 = yvec2 / mag
        yvec3 = yvec3 / mag

        xvec1 = yvec2 * zvec3 - yvec3 * zvec2
        xvec2 = yvec3 * zvec1 - yvec1 * zvec3
        xvec3 = yvec1 * zvec2 - yvec2 * zvec1
        mag = sqrt(xvec1*xvec1 + xvec2*xvec2 + xvec3*xvec3)
        xvec1 = xvec1 / mag
        xvec2 = xvec2 / mag
        xvec3 = xvec3 / mag

        bladeAlignedVectors(m,n,q,1,1) = xvec1
        bladeAlignedVectors(m,n,q,1,2) = xvec2
        bladeAlignedVectors(m,n,q,1,3) = xvec3
        bladeAlignedVectors(m,n,q,2,1) = yvec1
        bladeAlignedVectors(m,n,q,2,2) = yvec2
        bladeAlignedVectors(m,n,q,2,3) = yvec3
        bladeAlignedVectors(m,n,q,3,1) = zvec1
        bladeAlignedVectors(m,n,q,3,2) = zvec2
        bladeAlignedVectors(m,n,q,3,3) = zvec3

        wind1 = xvec1*uLocal1 + xvec2*uLocal2 + xvec3*uLocal3
        wind2 = yvec1*uLocal1 + yvec2*uLocal2 + yvec3*uLocal3 +              &
            rotSpeed * bladeRadius(m,n,q) * cosPreCone
        wind3 = zvec1*uLocal1 + zvec2*uLocal2 + zvec3*uLocal3

        windVectors(m,n,q,1) = wind1
        windVectors(m,n,q,2) = wind2
        windVectors(m,n,q,3) = wind3

        if (tipCorrection) then
            wind1 = wind1 + du(m,n,q,1)
            wind2 = wind2 + du(m,n,q,2)
            windVectors(m,n,q,1) = wind1
            windVectors(m,n,q,2) = wind2
        endif

        twistAng_i = twistAng(m,n,q)
        chord_i = chord(m,n,q)
        db_i = db(q)
        sectionType_i = sectionType(m,n,q)
        vmag_i = sqrt(wind1*wind1 + wind2*wind2)
        Vmag(m,n,q) = vmag_i

        windAng = atan2(wind1, wind2) / degRadLocal
        alpha_i = windAng - twistAng_i - pitch - pitchControlAngle
        alpha(m,n,q) = alpha_i

        kair = airfoil(sectionType_i) % n
        cl_i = 0._rprec
        if (alpha_i <= airfoil(sectionType_i) % AOA(1)) then
            cl_i = airfoil(sectionType_i) % Cl(1)
        elseif (alpha_i >= airfoil(sectionType_i) % AOA(kair)) then
            cl_i = airfoil(sectionType_i) % Cl(kair)
        else
            do kk = 2, kair
                if (alpha_i >= airfoil(sectionType_i) % AOA(kk-1) .and.      &
                    alpha_i <= airfoil(sectionType_i) % AOA(kk)) then
                    xa = airfoil(sectionType_i) % AOA(kk-1)
                    xb = airfoil(sectionType_i) % AOA(kk)
                    ya = airfoil(sectionType_i) % Cl(kk-1)
                    yb = airfoil(sectionType_i) % Cl(kk)
                    cl_i = ya + (yb - ya) * (alpha_i - xa) / (xb - xa)
                endif
            end do
        endif

        cd_i = 0._rprec
        if (alpha_i <= airfoil(sectionType_i) % AOA(1)) then
            cd_i = airfoil(sectionType_i) % Cd(1)
        elseif (alpha_i >= airfoil(sectionType_i) % AOA(kair)) then
            cd_i = airfoil(sectionType_i) % Cd(kair)
        else
            do kk = 2, kair
                if (alpha_i >= airfoil(sectionType_i) % AOA(kk-1) .and.      &
                    alpha_i <= airfoil(sectionType_i) % AOA(kk)) then
                    xa = airfoil(sectionType_i) % AOA(kk-1)
                    xb = airfoil(sectionType_i) % AOA(kk)
                    ya = airfoil(sectionType_i) % Cd(kk-1)
                    yb = airfoil(sectionType_i) % Cd(kk)
                    cd_i = ya + (yb - ya) * (alpha_i - xa) / (xb - xa)
                endif
            end do
        endif

        cl(m,n,q) = cl_i
        cd(m,n,q) = cd_i
        lift_i = 0.5_rprec * cl_i * (vmag_i*vmag_i) * chord_i * db_i *       &
            solidity(m,n,q)
        drag_i = 0.5_rprec * cd_i * (vmag_i*vmag_i) * chord_i * db_i *       &
            solidity(m,n,q)
        lift(m,n,q) = lift_i
        drag(m,n,q) = drag_i

        dvec1 = xvec1*wind1 + yvec1*wind2
        dvec2 = xvec2*wind1 + yvec2*wind2
        dvec3 = xvec3*wind1 + yvec3*wind2
        mag = sqrt(dvec1*dvec1 + dvec2*dvec2 + dvec3*dvec3)
        dvec1 = dvec1 / mag
        dvec2 = dvec2 / mag
        dvec3 = dvec3 / mag

        lvec1 = dvec2*zvec3 - dvec3*zvec2
        lvec2 = dvec3*zvec1 - dvec1*zvec3
        lvec3 = dvec1*zvec2 - dvec2*zvec1
        mag = sqrt(lvec1*lvec1 + lvec2*lvec2 + lvec3*lvec3)
        lvec1 = lvec1 / mag
        lvec2 = lvec2 / mag
        lvec3 = lvec3 / mag

        lvec1 = -lift_i * lvec1
        lvec2 = -lift_i * lvec2
        lvec3 = -lift_i * lvec3
        dvec1 = -drag_i * dvec1
        dvec2 = -drag_i * dvec2
        dvec3 = -drag_i * dvec3

        bladeForces(m,n,q,1) = lvec1 + dvec1
        bladeForces(m,n,q,2) = lvec2 + dvec2
        bladeForces(m,n,q,3) = lvec3 + dvec3
        axialForce(m,n,q) = -(bladeForces(m,n,q,1)*uvShaft1 +                &
            bladeForces(m,n,q,2)*uvShaft2 + bladeForces(m,n,q,3)*uvShaft3)
        tangentialForce(m,n,q) = bladeForces(m,n,q,1)*yvec1 +                &
            bladeForces(m,n,q,2)*yvec2 + bladeForces(m,n,q,3)*yvec3
        thrust_sum = thrust_sum + axialForce(m,n,q)
        torque_sum = torque_sum + tangentialForce(m,n,q)*bladeRadius(m,n,q)

        windAng = windAng * degRadLocal
        sigma = chord_i * real(mmend, rprec) /                               &
            (2._rprec*pi*bladeRadius(m,n,q))
        induction_a(m,n,q) = 1._rprec /                                      &
            (4._rprec * sin(windAng)**2 /                                    &
            (sigma * (cl_i*cos(windAng) + cd_i*sin(windAng))) + 1._rprec)
        u_infinity(m,n,q) = wind1
    endif
end do
end do
end do

call atm_cuda_sync('ATM point-owner blade force')
turbineArray(i) % thrust = thrust_sum
turbineArray(i) % torqueRotor = torque_sum * cosPreCone

end subroutine atm_point_owner_force_turbine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_point_owner_targeted_supported()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_point_owner_targeted_supported = atm_point_owner_targeted_enabled() .and. &
    nproc == 2 .and. atm_lb_owner_ready

end function atm_point_owner_targeted_supported

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lb_ensure_targeted_buffers()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: nsend, nrecv

nsend = max(1, 3 * atm_lb_sample_send_count)
nrecv = max(1, 3 * atm_lb_sample_recv_count)
if (allocated(atm_lb_sample_sendbuf)) then
    if (size(atm_lb_sample_sendbuf) /= nsend) deallocate(atm_lb_sample_sendbuf)
endif
if (allocated(atm_lb_sample_recvbuf)) then
    if (size(atm_lb_sample_recvbuf) /= nrecv) deallocate(atm_lb_sample_recvbuf)
endif
if (.not. allocated(atm_lb_sample_sendbuf)) allocate(atm_lb_sample_sendbuf(nsend))
if (.not. allocated(atm_lb_sample_recvbuf)) allocate(atm_lb_sample_recvbuf(nrecv))

nsend = max(1, 3 * atm_lb_force_send_count)
nrecv = max(1, 3 * atm_lb_force_recv_count)
if (allocated(atm_lb_force_sendbuf)) then
    if (size(atm_lb_force_sendbuf) /= nsend) deallocate(atm_lb_force_sendbuf)
endif
if (allocated(atm_lb_force_recvbuf)) then
    if (size(atm_lb_force_recvbuf) /= nrecv) deallocate(atm_lb_force_recvbuf)
endif
if (.not. allocated(atm_lb_force_sendbuf)) allocate(atm_lb_force_sendbuf(nsend))
if (.not. allocated(atm_lb_force_recvbuf)) allocate(atm_lb_force_recvbuf(nrecv))

end subroutine atm_lb_ensure_targeted_buffers

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lb_refresh_targeted_sample_slots()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i, j, m, n, q, pid, mmend, nnend, qqend, peer_rank
integer, allocatable :: local_grid_owner(:), global_grid_owner(:)
real(rprec) :: zloc
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladePoints

if (.not. atm_point_owner_targeted_supported()) return

allocate(local_grid_owner(atm_lb_npoints), global_grid_owner(atm_lb_npoints))
local_grid_owner = 0
global_grid_owner = 0

do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    mmend = turbineModel(j) % numBl
    nnend = turbineArray(i) % numAnnulusSections
    qqend = turbineArray(i) % numBladePoints
    bladePoints => turbineArray(i) % bladePoints
    do q = 1, qqend
    do n = 1, nnend
    do m = 1, mmend
        pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m
        zloc = bladePoints(m,n,q,3) / z_i
        if (zloc >= grid % z(1) .and. zloc < grid % z(nz)) then
            local_grid_owner(pid) = coord + 1
        endif
    end do
    end do
    end do
enddo

#ifdef PPMPI
call mpi_allreduce(local_grid_owner, global_grid_owner, atm_lb_npoints,       &
                   MPI_INTEGER, MPI_MAX, comm, ierr)
#else
global_grid_owner = local_grid_owner
#endif

atm_lb_sample_send_slot = 0
atm_lb_sample_recv_slot = 0
atm_lb_sample_send_count = 0
atm_lb_sample_recv_count = 0
peer_rank = 1 - coord

do pid = 1, atm_lb_npoints
    atm_lb_grid_owner(pid) = global_grid_owner(pid) - 1
    if (atm_lb_grid_owner(pid) == coord .and. atm_lb_owner(pid) ==           &
        peer_rank) then
        atm_lb_sample_send_count = atm_lb_sample_send_count + 1
        atm_lb_sample_send_slot(pid) = atm_lb_sample_send_count
    elseif (atm_lb_grid_owner(pid) == peer_rank .and. atm_lb_owner(pid) ==   &
            coord) then
        atm_lb_sample_recv_count = atm_lb_sample_recv_count + 1
        atm_lb_sample_recv_slot(pid) = atm_lb_sample_recv_count
    endif
enddo

deallocate(local_grid_owner, global_grid_owner)

end subroutine atm_lb_refresh_targeted_sample_slots

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lb_targeted_velocity_exchange()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: pid, slot, peer, nsend, nrecv
real(rprec) :: t0, t1, t2, t3
#ifdef PPMPI
integer :: mpi_status(MPI_STATUS_SIZE)
#endif

call atm_lb_ensure_targeted_buffers()

#ifdef PPMPI
t0 = mpi_wtime()
#else
call cpu_time(t0)
#endif

atm_lb_vel_global = 0._rprec
!$cuf kernel do(1) <<<*,*>>>
do pid = 1, atm_lb_npoints
    if (atm_lb_grid_owner(pid) == coord .and. atm_lb_owner(pid) == coord) then
        atm_lb_vel_global(1,pid) = atm_lb_vel_local(1,pid)
        atm_lb_vel_global(2,pid) = atm_lb_vel_local(2,pid)
        atm_lb_vel_global(3,pid) = atm_lb_vel_local(3,pid)
    endif

    slot = atm_lb_sample_send_slot(pid)
    if (slot > 0) then
        atm_lb_sample_sendbuf(3*slot-2) = atm_lb_vel_local(1,pid)
        atm_lb_sample_sendbuf(3*slot-1) = atm_lb_vel_local(2,pid)
        atm_lb_sample_sendbuf(3*slot  ) = atm_lb_vel_local(3,pid)
    endif
end do

call atm_cuda_sync('ATM targeted sample pack')
#ifdef PPMPI
t1 = mpi_wtime()
#else
call cpu_time(t1)
#endif

#ifdef PPMPI
if (nproc == 2) then
    peer = 1 - coord
    nsend = 3 * atm_lb_sample_send_count
    nrecv = 3 * atm_lb_sample_recv_count
    call mpi_sendrecv(atm_lb_sample_sendbuf, nsend, mpi_rprec, peer, 7101,   &
                      atm_lb_sample_recvbuf, nrecv, mpi_rprec, peer, 7101,   &
                      comm, mpi_status, ierr)
endif
t2 = mpi_wtime()
#else
call cpu_time(t2)
#endif

!$cuf kernel do(1) <<<*,*>>>
do pid = 1, atm_lb_npoints
    slot = atm_lb_sample_recv_slot(pid)
    if (slot > 0) then
        atm_lb_vel_global(1,pid) = atm_lb_sample_recvbuf(3*slot-2)
        atm_lb_vel_global(2,pid) = atm_lb_sample_recvbuf(3*slot-1)
        atm_lb_vel_global(3,pid) = atm_lb_sample_recvbuf(3*slot  )
    endif
end do

call atm_cuda_sync('ATM targeted sample unpack')
#ifdef PPMPI
t3 = mpi_wtime()
#else
call cpu_time(t3)
#endif

atm_lb_time_sample_pack = atm_lb_time_sample_pack + (t1 - t0)
atm_lb_time_sample_mpi = atm_lb_time_sample_mpi + (t2 - t1)
atm_lb_time_sample_unpack = atm_lb_time_sample_unpack + (t3 - t2)
atm_lb_sample_bytes_sent = atm_lb_sample_bytes_sent +                        &
    int(3 * atm_lb_sample_send_count, 8) * int(storage_size(1.0_rprec) / 8, 8)
atm_lb_sample_bytes_recv = atm_lb_sample_bytes_recv +                        &
    int(3 * atm_lb_sample_recv_count, 8) * int(storage_size(1.0_rprec) / 8, 8)

end subroutine atm_lb_targeted_velocity_exchange

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lb_pack_force_turbine(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, m, n, q, pid, slot, mmend, nnend, qqend
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladeForces

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
bladeForces => turbineArray(i) % bladeForces

#ifdef PPLES_GPU
if (atm_blade_mirror_ready) then
    !$cuf kernel do(3) <<<*,*>>>
    do q = 1, qqend
    do n = 1, nnend
    do m = 1, mmend
        pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m
        slot = atm_lb_force_send_slot(pid)
        if (slot > 0) then
            atm_lb_force_sendbuf(3*slot-2) = atm_bladeForces_d(m,n,q,1,i)
            atm_lb_force_sendbuf(3*slot-1) = atm_bladeForces_d(m,n,q,2,i)
            atm_lb_force_sendbuf(3*slot  ) = atm_bladeForces_d(m,n,q,3,i)
        endif
    end do
    end do
    end do
    return
endif
#endif

!$cuf kernel do(3) <<<*,*>>>
do q = 1, qqend
do n = 1, nnend
do m = 1, mmend
    pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m
    slot = atm_lb_force_send_slot(pid)
    if (slot > 0) then
        atm_lb_force_sendbuf(3*slot-2) = bladeForces(m,n,q,1)
        atm_lb_force_sendbuf(3*slot-1) = bladeForces(m,n,q,2)
        atm_lb_force_sendbuf(3*slot  ) = bladeForces(m,n,q,3)
    endif
end do
end do
end do

end subroutine atm_lb_pack_force_turbine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lb_unpack_force_turbine(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, m, n, q, pid, slot, mmend, nnend, qqend
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladeForces

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
bladeForces => turbineArray(i) % bladeForces

#ifdef PPLES_GPU
if (atm_blade_mirror_ready) then
    !$cuf kernel do(3) <<<*,*>>>
    do q = 1, qqend
    do n = 1, nnend
    do m = 1, mmend
        pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m
        slot = atm_lb_force_recv_slot(pid)
        if (slot > 0) then
            atm_bladeForces_d(m,n,q,1,i) = atm_lb_force_recvbuf(3*slot-2)
            atm_bladeForces_d(m,n,q,2,i) = atm_lb_force_recvbuf(3*slot-1)
            atm_bladeForces_d(m,n,q,3,i) = atm_lb_force_recvbuf(3*slot  )
        endif
    end do
    end do
    end do
    return
endif
#endif

!$cuf kernel do(3) <<<*,*>>>
do q = 1, qqend
do n = 1, nnend
do m = 1, mmend
    pid = atm_lb_offset(i) + (q-1) * nnend * mmend + (n-1) * mmend + m
    slot = atm_lb_force_recv_slot(pid)
    if (slot > 0) then
        bladeForces(m,n,q,1) = atm_lb_force_recvbuf(3*slot-2)
        bladeForces(m,n,q,2) = atm_lb_force_recvbuf(3*slot-1)
        bladeForces(m,n,q,3) = atm_lb_force_recvbuf(3*slot  )
    endif
end do
end do
end do

end subroutine atm_lb_unpack_force_turbine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_point_owner_lb_gather_targeted()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i, pos, peer, nsend, nrecv, scalar_count
real(rprec), allocatable, save :: scalar_send(:), scalar_recv(:)
real(rprec) :: t0, t1, t2, t3, t4
#ifdef PPMPI
integer :: mpi_status(MPI_STATUS_SIZE)
#endif

if (.not. atm_point_owner_targeted_supported()) then
    call atm_point_owner_lb_gather()
    return
endif

call atm_lb_ensure_targeted_buffers()

#ifdef PPMPI
t0 = mpi_wtime()
#else
call cpu_time(t0)
#endif

do i = 1, numberOfTurbines
    call atm_lb_pack_force_turbine(i)
enddo
call atm_cuda_sync('ATM targeted force pack')
#ifdef PPMPI
t1 = mpi_wtime()
#else
call cpu_time(t1)
#endif

#ifdef PPMPI
peer = 1 - coord
nsend = 3 * atm_lb_force_send_count
nrecv = 3 * atm_lb_force_recv_count
call mpi_sendrecv(atm_lb_force_sendbuf, nsend, mpi_rprec, peer, 7102,        &
                  atm_lb_force_recvbuf, nrecv, mpi_rprec, peer, 7102,        &
                  comm, mpi_status, ierr)
t2 = mpi_wtime()
#else
call cpu_time(t2)
#endif

do i = 1, numberOfTurbines
    call atm_lb_unpack_force_turbine(i)
enddo
call atm_cuda_sync('ATM targeted force unpack')
#ifdef PPMPI
t3 = mpi_wtime()
#else
call cpu_time(t3)
#endif

scalar_count = 7 * numberOfTurbines
if (allocated(scalar_send)) then
    if (size(scalar_send) /= scalar_count) deallocate(scalar_send, scalar_recv)
endif
if (.not. allocated(scalar_send)) allocate(scalar_send(scalar_count),         &
    scalar_recv(scalar_count))

pos = 1
do i = 1, numberOfTurbines
    scalar_send(pos) = turbineArray(i) % torqueRotor
    scalar_send(pos+1) = turbineArray(i) % thrust
    scalar_send(pos+2:pos+4) = turbineArray(i) % nacelleForce
    scalar_send(pos+5) = turbineArray(i) % VelNacelle_sampled
    scalar_send(pos+6) = turbineArray(i) % VelNacelle_corrected
    pos = pos + 7
enddo

#ifdef PPMPI
call mpi_allreduce(scalar_send, scalar_recv, scalar_count, mpi_rprec,        &
                   mpi_sum, comm, ierr)
t4 = mpi_wtime()
#else
scalar_recv = scalar_send
call cpu_time(t4)
#endif

pos = 1
do i = 1, numberOfTurbines
    turbineArray(i) % torqueRotor = scalar_recv(pos)
    turbineArray(i) % thrust = scalar_recv(pos+1)
    turbineArray(i) % nacelleForce = scalar_recv(pos+2:pos+4)
    turbineArray(i) % VelNacelle_sampled = scalar_recv(pos+5)
    turbineArray(i) % VelNacelle_corrected = scalar_recv(pos+6)
    pos = pos + 7
enddo

atm_lb_time_force_pack = atm_lb_time_force_pack + (t1 - t0)
atm_lb_time_force_mpi = atm_lb_time_force_mpi + (t2 - t1)
atm_lb_time_force_unpack = atm_lb_time_force_unpack + (t3 - t2)
atm_lb_time_scalar_reduce = atm_lb_time_scalar_reduce + (t4 - t3)
atm_lb_time_gather = atm_lb_time_gather + (t4 - t0)
atm_lb_force_bytes_sent = atm_lb_force_bytes_sent +                          &
    int(3 * atm_lb_force_send_count, 8) * int(storage_size(1.0_rprec) / 8, 8)
atm_lb_force_bytes_recv = atm_lb_force_bytes_recv +                          &
    int(3 * atm_lb_force_recv_count, 8) * int(storage_size(1.0_rprec) / 8, 8)

end subroutine atm_point_owner_lb_gather_targeted

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_point_owner_lb_force()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i
real(rprec) :: t0, t1, t2, t3

if (.not. atm_lb_owner_ready) return

if (.not. allocated(atm_lb_vel_local)) then
    allocate(atm_lb_vel_local(3,atm_lb_npoints),                              &
             atm_lb_vel_global(3,atm_lb_npoints))
endif

if (atm_point_owner_targeted_supported()) call atm_lb_refresh_targeted_sample_slots()

#ifdef PPMPI
t0 = mpi_wtime()
#else
call cpu_time(t0)
#endif
atm_lb_vel_local = 0._rprec
do i = 1, numberOfTurbines
    call atm_point_owner_sample_turbine(i)
enddo
call atm_cuda_sync('ATM point-owner velocity sampling')
#ifdef PPMPI
t1 = mpi_wtime()
#else
call cpu_time(t1)
#endif
atm_lb_time_sample = atm_lb_time_sample + (t1 - t0)

#ifdef PPMPI
if (atm_point_owner_targeted_supported()) then
    atm_lb_targeted_calls = atm_lb_targeted_calls + 1
    call atm_lb_targeted_velocity_exchange()
elseif (nproc > 1) then
    call mpi_allreduce(atm_lb_vel_local, atm_lb_vel_global,                  &
                       size(atm_lb_vel_local), mpi_rprec, mpi_sum, comm, ierr)
else
    atm_lb_vel_global = atm_lb_vel_local
endif
t2 = mpi_wtime()
#else
atm_lb_vel_global = atm_lb_vel_local
call cpu_time(t2)
#endif
atm_lb_time_velocity_comm = atm_lb_time_velocity_comm + (t2 - t1)

do i = 1, numberOfTurbines
    call atm_point_owner_force_turbine(i)
    if (turbineArray(i) % nacelle) call atm_lesgo_nacelle_force_gpu_atpoint(i)
enddo
#ifdef PPMPI
t3 = mpi_wtime()
#else
call cpu_time(t3)
#endif
atm_lb_time_point_force = atm_lb_time_point_force + (t3 - t2)

end subroutine atm_point_owner_lb_force

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_point_owner_lb_gather()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i, nitem, pos, total_pack
real(rprec), device, allocatable, save :: packed_send(:), packed_recv(:)
real(rprec) :: t0, t1

total_pack = 0
do i = 1, numberOfTurbines
    total_pack = total_pack + size(turbineArray(i) % bladeForces) + 7
enddo
if (total_pack <= 0) return

if (allocated(packed_send)) then
    if (size(packed_send) /= total_pack) then
        deallocate(packed_send, packed_recv)
    endif
endif
if (.not. allocated(packed_send)) then
    allocate(packed_send(total_pack), packed_recv(total_pack))
endif

#ifdef PPMPI
t0 = mpi_wtime()
#else
call cpu_time(t0)
#endif

pos = 1
do i = 1, numberOfTurbines
    nitem = size(turbineArray(i) % bladeForces)
    call atm_pack_rank4(turbineArray(i) % bladeForces, packed_send, pos)
    pos = pos + nitem
    call atm_pack_gather_scalars(packed_send, pos,                            &
        turbineArray(i) % torqueRotor, turbineArray(i) % thrust,              &
        turbineArray(i) % nacelleForce(1), turbineArray(i) % nacelleForce(2), &
        turbineArray(i) % nacelleForce(3),                                    &
        turbineArray(i) % VelNacelle_sampled,                                 &
        turbineArray(i) % VelNacelle_corrected)
    pos = pos + 7
enddo

call atm_cuda_sync('ATM point-owner packed gather pack')
#ifdef PPMPI
if (nproc > 1) then
    call mpi_allreduce(packed_send, packed_recv, total_pack, mpi_rprec,       &
                       mpi_sum, comm, ierr)
else
    packed_recv = packed_send
endif
#else
packed_recv = packed_send
#endif

pos = 1
do i = 1, numberOfTurbines
    nitem = size(turbineArray(i) % bladeForces)
    call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % bladeForces)
    pos = pos + nitem
    turbineArray(i) % torqueRotor = packed_recv(pos)
    turbineArray(i) % thrust = packed_recv(pos+1)
    turbineArray(i) % nacelleForce = packed_recv(pos+2:pos+4)
    turbineArray(i) % VelNacelle_sampled = packed_recv(pos+5)
    turbineArray(i) % VelNacelle_corrected = packed_recv(pos+6)
    pos = pos + 7
enddo
call atm_cuda_check('ATM point-owner packed gather unpack')

#ifdef PPMPI
t1 = mpi_wtime()
#else
call cpu_time(t1)
#endif
atm_lb_time_gather = atm_lb_time_gather + (t1 - t0)

end subroutine atm_point_owner_lb_gather

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_point_owner_lb_validate(label)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: label
integer :: i, c, p
real(rprec) :: local_force(3), global_force(3), t0, t1
real(rprec) :: power_rotor
real(rprec) :: local_vel_checksum, global_vel_checksum, blade_checksum
real(rprec), allocatable :: vel_h(:,:)

if (.not. atm_lb_validate_enabled()) return

#ifdef PPMPI
t0 = mpi_wtime()
#else
call cpu_time(t0)
#endif

local_force = 0._rprec
do i = 1, numberOfTurbines
    if (turbineArray(i) % operate) then
        do c = 1, forceFieldUV(i) % c
            local_force(1) = local_force(1) + forceFieldUV(i) % force(1,c)
            local_force(2) = local_force(2) + forceFieldUV(i) % force(2,c)
        enddo
        do c = 1, forceFieldW(i) % c
            local_force(3) = local_force(3) + forceFieldW(i) % force(3,c)
        enddo
    endif
enddo

#ifdef PPMPI
call mpi_allreduce(local_force, global_force, 3, mpi_rprec, mpi_sum, comm, ierr)
#else
global_force = local_force
#endif

local_vel_checksum = 0._rprec
global_vel_checksum = 0._rprec
if (allocated(atm_lb_vel_global)) then
    allocate(vel_h(3,atm_lb_npoints))
    vel_h = atm_lb_vel_global
    do p = 1, atm_lb_npoints
        if (atm_lb_owner(p) == coord) then
            local_vel_checksum = local_vel_checksum + vel_h(1,p) +           &
                10._rprec * vel_h(2,p) + 100._rprec * vel_h(3,p)
        endif
    enddo
    deallocate(vel_h)
#ifdef PPMPI
    call mpi_allreduce(local_vel_checksum, global_vel_checksum, 1,           &
                       mpi_rprec, mpi_sum, comm, ierr)
#else
    global_vel_checksum = local_vel_checksum
#endif
endif

blade_checksum = 0._rprec
if (coord == 0) then
    do i = 1, numberOfTurbines
        blade_checksum = blade_checksum + sum(turbineArray(i) % bladeForces)
    enddo
endif

if (coord == 0) then
    write(*,'(a,a,a,i0,3(a,es16.8))') 'ATM_LB_VALIDATE label=', trim(label), &
        ' call=', atm_forcing_calls, ' force_sum_x=', global_force(1),        &
        ' force_sum_y=', global_force(2), ' force_sum_z=', global_force(3)
    write(*,'(a,a,a,i0,2(a,es16.8))') 'ATM_LB_VALIDATE_CHECKSUM label=',     &
        trim(label), ' call=', atm_forcing_calls,                             &
        ' sampled_velocity=', global_vel_checksum,                            &
        ' blade_force=', blade_checksum
    do i = 1, numberOfTurbines
        power_rotor = turbineArray(i) % torqueRotor *                        &
            turbineArray(i) % rotSpeed * turbineArray(i) % fluidDensity
        write(*,'(a,a,a,i0,a,i0,3(a,es16.8))')                               &
            'ATM_LB_VALIDATE_TURBINE label=', trim(label), ' call=',          &
            atm_forcing_calls, ' turbine=', i, ' thrust=',                    &
            turbineArray(i) % thrust, ' torque=', turbineArray(i) %           &
            torqueRotor, ' power=', power_rotor
    enddo
    flush(6)
endif

#ifdef PPMPI
t1 = mpi_wtime()
#else
call cpu_time(t1)
#endif
atm_lb_time_validate = atm_lb_time_validate + (t1 - t0)

end subroutine atm_point_owner_lb_validate
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_initialize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

! Counter to establish number of points which are influenced by body forces
integer ::  m

! Allocate space for the w_uv variable
allocate(w_uv(nx,ny,lbz:nz))
w_uv = 0._rprec
#ifdef PPLES_GPU
!$acc update device(w_uv)
#endif

call atm_initialize () ! Initialize the atm (ATM)

! Allocate the body force variables. It is an array with one per turbine.
allocate(forceFieldUV(numberOfTurbines))
allocate(forceFieldW(numberOfTurbines))

    do m=1, numberOfTurbines
        call atm_lesgo_findCells(m)
    enddo

    #ifdef PPLES_GPU
    call atm_lesgo_build_force_shadows()
    call atm_lesgo_build_blade_mirrors()
    #endif

    #ifdef PPMPI
        call mpi_barrier( comm, ierr )
    #endif

#ifdef ENABLE_CUDA
call atm_lesgo_diag_load()
#endif

#ifdef PPMPI
    ! This will create the output files and write initialization to the screen
    if (coord == 0) then
        call atm_initialize_output()
    endif
#else
    call atm_initialize_output()
#endif

#ifdef ENABLE_CUDA
call atm_lesgo_lb_plan()
#endif

end subroutine atm_lesgo_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_finalize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

! Counter for turbines
integer ::  i

#ifdef ENABLE_CUDA
call atm_diag_event_flush()
#endif
call atm_lesgo_report_timing()
call atm_structure_timing_report()

#ifdef PPLES_GPU
call atm_lesgo_destroy_force_shadows()
call atm_lesgo_destroy_blade_mirrors()
#endif

! Write if on main node
if (coord == 0) then
    write(*,*) 'Finalizing ATM...'
endif

    ! Loop through all turbines and finalize
    do i = 1, numberOfTurbines
        if (coord == turbineArray(i) % master) then
            call atm_write_restart(i) ! Write the restart file
        endif
    end do

if (coord == 0) then
    write(*,*) 'Done finalizing ATM'
endif

end subroutine atm_lesgo_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_report_timing()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

real(rprec) :: vals(11), maxvals(11), measured_total
#ifdef ENABLE_CUDA
real(rprec) :: lb_vals(8), lb_maxvals(8), lb_bytes(4), lb_bytes_sum(4)
#endif

vals = (/ atm_time_interp_w, atm_time_update, atm_time_reset, atm_time_sample, &
          atm_time_force, atm_time_gather, atm_time_convolve, atm_time_clcorr, &
          atm_time_apply, atm_time_output, atm_time_barrier /)

#ifdef PPMPI
call mpi_allreduce(vals, maxvals, size(vals), mpi_rprec, mpi_max, comm, ierr)
#else
maxvals = vals
#endif

measured_total = sum(maxvals)

if (coord == 0) then
    write(*,*) '==================================================='
    write(*,*) 'ATM Cumulative Times (s, max over ranks):'
    write(*,'(1a,I8)')    '  ATM forcing calls: ', atm_forcing_calls
    write(*,'(1a,E15.7)') '  w -> uv interpolation: ', maxvals(1)
    write(*,'(1a,E15.7)') '  turbine update/yaw/rotation: ', maxvals(2)
    write(*,'(1a,E15.7)') '  turbine/force reset: ', maxvals(3)
    write(*,'(1a,E15.7)') '  batched velocity sampling: ', maxvals(4)
    write(*,'(1a,E15.7)') '  blade/nacelle force: ', maxvals(5)
    write(*,'(1a,E15.7)') '  MPI gather: ', maxvals(6)
    write(*,'(1a,E15.7)') '  force convolution: ', maxvals(7)
    write(*,'(1a,E15.7)') '  tip correction: ', maxvals(8)
    write(*,'(1a,E15.7)') '  apply force to grid: ', maxvals(9)
    write(*,'(1a,E15.7)') '  ATM output: ', maxvals(10)
    write(*,'(1a,E15.7)') '  ATM output barrier: ', maxvals(11)
    write(*,'(1a,E15.7)') '  ATM measured subtotal: ', measured_total
#ifdef ENABLE_CUDA
    if (atm_point_owner_lb_enabled()) then
        write(*,*) 'ATM point-owner LB cumulative times (s):'
        write(*,'(1a,E15.7)') '  sample velocity kernels: ', atm_lb_time_sample
        write(*,'(1a,E15.7)') '  velocity exchange/allreduce: ',              &
            atm_lb_time_velocity_comm
        write(*,'(1a,E15.7)') '  point-owner blade/nacelle force: ',          &
            atm_lb_time_point_force
        write(*,'(1a,E15.7)') '  global packed force gather: ',               &
            atm_lb_time_gather
        write(*,'(1a,E15.7)') '  validation reductions: ', atm_lb_time_validate
        if (atm_lb_auto_select_enabled()) then
            write(*,*) 'ATM point-owner LB auto-selection summary:'
            write(*,'(1a,I8)') '  auto probe calls: ', atm_lb_auto_probe_calls
            write(*,'(1a,I8)') '  legacy samples: ', atm_lb_auto_legacy_samples
            write(*,'(1a,I8)') '  LB samples: ', atm_lb_auto_lb_samples
            if (atm_lb_auto_legacy_samples > 0) then
                write(*,'(1a,E15.7)') '  legacy avg force+gather: ',         &
                    atm_lb_auto_legacy_time /                                &
                    real(atm_lb_auto_legacy_samples, rprec)
            endif
            if (atm_lb_auto_lb_samples > 0) then
                write(*,'(1a,E15.7)') '  LB avg force+gather: ',             &
                    atm_lb_auto_lb_time / real(atm_lb_auto_lb_samples, rprec)
            endif
            write(*,'(1a,L1)') '  selected LB: ', atm_lb_auto_use_lb
        endif
    endif
#endif
    write(*,*) '==================================================='
end if

#ifdef ENABLE_CUDA
if (atm_point_owner_lb_enabled() .and. atm_point_owner_targeted_enabled()) then
    lb_vals = (/ atm_lb_time_sample_pack, atm_lb_time_sample_mpi,             &
                 atm_lb_time_sample_unpack, atm_lb_time_force_pack,           &
                 atm_lb_time_force_mpi, atm_lb_time_force_unpack,             &
                 atm_lb_time_scalar_reduce, atm_lb_time_gather /)
    lb_bytes = (/ real(atm_lb_sample_bytes_sent, rprec),                     &
                  real(atm_lb_sample_bytes_recv, rprec),                     &
                  real(atm_lb_force_bytes_sent, rprec),                      &
                  real(atm_lb_force_bytes_recv, rprec) /)
#ifdef PPMPI
    call mpi_allreduce(lb_vals, lb_maxvals, size(lb_vals), mpi_rprec,         &
                       mpi_max, comm, ierr)
    call mpi_allreduce(lb_bytes, lb_bytes_sum, size(lb_bytes), mpi_rprec,     &
                       mpi_sum, comm, ierr)
#else
    lb_maxvals = lb_vals
    lb_bytes_sum = lb_bytes
#endif
    if (coord == 0) then
        write(*,*) 'ATM targeted point-owner communication summary:'
        write(*,'(1a,I8)') '  targeted active calls: ', atm_lb_targeted_calls
        write(*,'(a,4(a,E15.7))') 'ATM_LB_TABLE_A path=targeted',            &
            ' sample_bytes_sent=', lb_bytes_sum(1),                          &
            ' sample_bytes_recv=', lb_bytes_sum(2),                          &
            ' force_bytes_sent=', lb_bytes_sum(3),                           &
            ' force_bytes_recv=', lb_bytes_sum(4)
        write(*,'(a,8(a,E15.7))') 'ATM_LB_TABLE_B path=targeted',            &
            ' sample_pack=', lb_maxvals(1),                                  &
            ' sample_mpi=', lb_maxvals(2),                                   &
            ' sample_unpack=', lb_maxvals(3),                                &
            ' force_pack=', lb_maxvals(4),                                   &
            ' force_mpi=', lb_maxvals(5),                                    &
            ' force_unpack=', lb_maxvals(6),                                 &
            ' scalar_reduce=', lb_maxvals(7),                                &
            ' force_exchange_total=', lb_maxvals(8)
        write(*,'(a)')                                                       &
            'ATM_LB_FORCE_RETURN note=current convolution reads full ' //     &
            'bladeForces per turbine; targeted path returns owner-computed ' // &
            'blade forces to the peer rank instead of using a full allreduce.'
        write(*,*) '==================================================='
    endif
    write(*,'(a,i0,7(a,E15.7))') 'ATM_LB_TABLE_C path=targeted rank=', coord,&
        ' point_owner_points=', atm_lb_planned_points_rank,                  &
        ' support_ops=', atm_lb_planned_support_rank,                        &
        ' sample_bytes_sent=', real(atm_lb_sample_bytes_sent, rprec),        &
        ' sample_bytes_recv=', real(atm_lb_sample_bytes_recv, rprec),        &
        ' force_bytes_sent=', real(atm_lb_force_bytes_sent, rprec),          &
        ' force_bytes_recv=', real(atm_lb_force_bytes_recv, rprec),          &
        ' mpi_time=', atm_lb_time_sample_mpi + atm_lb_time_force_mpi +       &
            atm_lb_time_scalar_reduce
    flush(6)
endif
#endif

end subroutine atm_lesgo_report_timing

#ifdef ENABLE_CUDA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_diag_load()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i, j, m, n, q, r
integer :: turbines_owned, actuator_owned, min_k, max_k, kidx
integer :: points_i
integer(8) :: support_touch, scatter_ops
real(rprec) :: zmin, zmax, znd, dz_safe

if (.not. atm_diag_timing_enabled()) return
if (atm_diag_load_printed) return

turbines_owned = 0
actuator_owned = 0
support_touch = 0_8
scatter_ops = 0_8
min_k = huge(1)
max_k = -huge(1)
zmin = grid % z(1)
zmax = grid % z(nz)
dz_safe = max(abs(dz), 1.0e-12_rprec)

do i = 1, numberOfTurbines
    if (coord == turbineArray(i) % master) turbines_owned = turbines_owned + 1
    support_touch = support_touch + int(forceFieldUV(i) % c, 8) +            &
        int(forceFieldW(i) % c, 8)
    j = turbineArray(i) % turbineTypeID
    points_i = turbineModel(j) % numBl * turbineArray(i) %                  &
        numAnnulusSections * turbineArray(i) % numBladePoints
    scatter_ops = scatter_ops + int(forceFieldUV(i) % c +                   &
        forceFieldW(i) % c, 8) * int(points_i, 8)
    do q = 1, turbineArray(i) % numBladePoints
    do n = 1, turbineArray(i) % numAnnulusSections
    do m = 1, turbineModel(j) % numBl
        znd = turbineArray(i) % bladePoints(m,n,q,3) / z_i
        if (znd >= zmin .and. znd < zmax) then
            actuator_owned = actuator_owned + 1
            kidx = int((znd - zmin) / dz_safe) + 1
            kidx = max(1, min(nz, kidx))
            min_k = min(min_k, kidx)
            max_k = max(max_k, kidx)
        end if
    end do
    end do
    end do
end do

if (min_k == huge(1)) then
    min_k = 0
    max_k = 0
end if

#ifdef PPMPI
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,'(a,i0,8(a,i0))') 'ATM_LOAD rank=', coord,                   &
            ' z_lo=', 1, ' z_hi=', nz, ' turbines_owned=', turbines_owned,   &
            ' actuator_owned=', actuator_owned, ' support_touch=',           &
            support_touch, ' min_k=', min_k, ' max_k=', max_k,               &
            ' scatter_ops=', scatter_ops
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
if (coord == 0) then
    write(*,'(a)') 'ATM_ORDER timestep_order=SGS_before_ATM_FORCING ATM_before_SGS_tau_halo=0'
    write(*,'(a)') 'ATM_ORDER queued_ATM_GPU_work_can_precede_same_step_SGS_tau_halo=0'
end if
#else
write(*,'(a,i0,8(a,i0))') 'ATM_LOAD rank=', coord, ' z_lo=', 1,              &
    ' z_hi=', nz, ' turbines_owned=', turbines_owned,                        &
    ' actuator_owned=', actuator_owned, ' support_touch=', support_touch,    &
    ' min_k=', min_k, ' max_k=', max_k, ' scatter_ops=', scatter_ops
flush(6)
write(*,'(a)') 'ATM_ORDER timestep_order=SGS_before_ATM_FORCING ATM_before_SGS_tau_halo=0'
#endif

atm_diag_load_printed = .true.

end subroutine atm_lesgo_diag_load

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_lb_plan_only_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_lb_plan_only_enabled = .false.

end function atm_lb_plan_only_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_lb_point_detail_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_lb_point_detail_enabled = .false.

end function atm_lb_point_detail_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) function atm_lb_safe_ratio(num, den)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

real(rprec), intent(in) :: num, den

if (den <= tiny(1.0_rprec)) then
    atm_lb_safe_ratio = huge(1.0_rprec)
else
    atm_lb_safe_ratio = num / den
end if

end function atm_lb_safe_ratio

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_lb_plan()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Diagnostic-only planner for future ATM point ownership.  It records the
! current point/support cost matrix and prints a greedy balanced plan, but it
! never writes force fields, turbine state, or the active ownership mapping.
implicit none

integer :: i, j, m, n, q, c, p, r, rr, best_rank, selected
integer :: npoints, nrank, rank_index, current_owner, remote_ranks
integer :: support_count, touched_ranks, point_id
integer :: points_i, env_bytes
integer :: peer_rank
integer, allocatable :: point_turb(:), point_m(:), point_n(:), point_q(:)
integer, allocatable :: point_owner_local(:), point_owner(:), point_plan(:)
integer, allocatable :: sort_order(:)
logical, allocatable :: used(:)
real(rprec), allocatable :: point_x(:), point_y(:), point_z(:)
real(rprec), allocatable :: support_local(:,:), support_global(:,:)
real(rprec), allocatable :: current_turbines(:), current_points(:)
real(rprec), allocatable :: current_support(:), current_scatter(:)
real(rprec), allocatable :: current_cost(:), current_control(:)
real(rprec), allocatable :: planned_points(:), planned_support(:)
real(rprec), allocatable :: planned_scatter(:), planned_sample_bytes(:)
real(rprec), allocatable :: planned_scatter_bytes(:), planned_cost(:)
real(rprec), allocatable :: local_rank_vals(:), global_rank_vals(:)
real(rprec), allocatable :: point_sort_cost(:)
real(rprec), allocatable :: turbine_plan_points(:,:)
real(rprec), allocatable :: turbine_plan_support(:,:)
real(rprec) :: alpha, lambda, control_cost, bytes_per_real
real(rprec) :: xloc, yloc, zloc, dist2, projectradius, projectradius2
real(rprec) :: support_total, candidate_score, best_score, remote_cells
real(rprec) :: sample_bytes, scatter_bytes, comm_words, point_cost
real(rprec) :: maxv, minv, meanv
real(rprec), dimension(3) :: point, cell_loc

if (.not. atm_lb_plan_only_enabled() .and.                                 &
    .not. atm_point_owner_lb_enabled()) return

nrank = max(1, nproc)
rank_index = min(max(coord + 1, 1), nrank)
bytes_per_real = real(storage_size(1.0_rprec) / 8, rprec)
env_bytes = int(bytes_per_real)
alpha = 64.0_rprec
lambda = 0.02_rprec
control_cost = 1024.0_rprec

npoints = 0
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    npoints = npoints + turbineModel(j) % numBl *                            &
        turbineArray(i) % numAnnulusSections *                               &
        turbineArray(i) % numBladePoints
end do

if (npoints <= 0) return

if (atm_point_owner_lb_enabled()) then
    if (allocated(atm_lb_owner)) then
        if (size(atm_lb_owner) /= npoints) then
            deallocate(atm_lb_owner)
            deallocate(atm_lb_offset)
            if (allocated(atm_lb_grid_owner)) deallocate(atm_lb_grid_owner)
            if (allocated(atm_lb_sample_send_slot))                           &
                deallocate(atm_lb_sample_send_slot)
            if (allocated(atm_lb_sample_recv_slot))                           &
                deallocate(atm_lb_sample_recv_slot)
            if (allocated(atm_lb_force_send_slot))                            &
                deallocate(atm_lb_force_send_slot)
            if (allocated(atm_lb_force_recv_slot))                            &
                deallocate(atm_lb_force_recv_slot)
            if (allocated(atm_lb_sample_sendbuf))                             &
                deallocate(atm_lb_sample_sendbuf)
            if (allocated(atm_lb_sample_recvbuf))                             &
                deallocate(atm_lb_sample_recvbuf)
            if (allocated(atm_lb_force_sendbuf))                              &
                deallocate(atm_lb_force_sendbuf)
            if (allocated(atm_lb_force_recvbuf))                              &
                deallocate(atm_lb_force_recvbuf)
            atm_lb_owner_ready = .false.
        endif
    endif
    if (.not. allocated(atm_lb_owner)) then
        allocate(atm_lb_owner(npoints), atm_lb_offset(numberOfTurbines))
    endif
    if (.not. allocated(atm_lb_grid_owner)) then
        allocate(atm_lb_grid_owner(npoints), atm_lb_sample_send_slot(npoints),&
                 atm_lb_sample_recv_slot(npoints),                           &
                 atm_lb_force_send_slot(npoints), atm_lb_force_recv_slot(npoints))
    endif
    atm_lb_npoints = npoints
    p = 0
    do i = 1, numberOfTurbines
        atm_lb_offset(i) = p
        j = turbineArray(i) % turbineTypeID
        p = p + turbineModel(j) % numBl *                                    &
            turbineArray(i) % numAnnulusSections *                           &
            turbineArray(i) % numBladePoints
    enddo
endif

allocate(point_turb(npoints), point_m(npoints), point_n(npoints),             &
    point_q(npoints))
allocate(point_owner_local(npoints), point_owner(npoints), point_plan(npoints))
allocate(sort_order(npoints), used(npoints))
allocate(point_x(npoints), point_y(npoints), point_z(npoints))
allocate(support_local(npoints,nrank), support_global(npoints,nrank))
allocate(current_turbines(nrank), current_points(nrank),                     &
    current_support(nrank), current_scatter(nrank), current_cost(nrank),      &
    current_control(nrank))
allocate(planned_points(nrank), planned_support(nrank),                      &
    planned_scatter(nrank), planned_sample_bytes(nrank),                     &
    planned_scatter_bytes(nrank), planned_cost(nrank))
allocate(local_rank_vals(nrank), global_rank_vals(nrank))
allocate(point_sort_cost(npoints))
allocate(turbine_plan_points(numberOfTurbines,nrank))
allocate(turbine_plan_support(numberOfTurbines,nrank))

point_owner_local = 0
point_owner = 0
point_plan = 1
support_local = 0.0_rprec
support_global = 0.0_rprec
current_turbines = 0.0_rprec
current_points = 0.0_rprec
current_support = 0.0_rprec
current_scatter = 0.0_rprec
current_cost = 0.0_rprec
current_control = 0.0_rprec
planned_points = 0.0_rprec
planned_support = 0.0_rprec
planned_scatter = 0.0_rprec
planned_sample_bytes = 0.0_rprec
planned_scatter_bytes = 0.0_rprec
planned_cost = 0.0_rprec
turbine_plan_points = 0.0_rprec
turbine_plan_support = 0.0_rprec
point_sort_cost = 0.0_rprec

p = 0
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    projectradius = turbineArray(i) % projectionRadius
    projectradius2 = projectradius * projectradius
    do q = 1, turbineArray(i) % numBladePoints
    do n = 1, turbineArray(i) % numAnnulusSections
    do m = 1, turbineModel(j) % numBl
        p = p + 1
        point_turb(p) = i
        point_m(p) = m
        point_n(p) = n
        point_q(p) = q
        point = turbineArray(i) % bladePoints(m,n,q,1:3)
        point_x(p) = point(1)
        point_y(p) = point(2)
        point_z(p) = point(3)

        zloc = point(3) / z_i
        if (zloc >= grid % z(1) .and. zloc < grid % z(nz)) then
            point_owner_local(p) = rank_index
        endif

        support_count = 0
        do c = 1, forceFieldUV(i) % c
            cell_loc = forceFieldUV(i) % location(1:3,c)
            xloc = cell_loc(1) - point(1)
            yloc = cell_loc(2) - point(2)
            zloc = cell_loc(3) - point(3)
            dist2 = xloc*xloc + yloc*yloc + zloc*zloc
            if (dist2 <= projectradius2) support_count = support_count + 1
        end do
        do c = 1, forceFieldW(i) % c
            cell_loc = forceFieldW(i) % location(1:3,c)
            xloc = cell_loc(1) - point(1)
            yloc = cell_loc(2) - point(2)
            zloc = cell_loc(3) - point(3)
            dist2 = xloc*xloc + yloc*yloc + zloc*zloc
            if (dist2 <= projectradius2) support_count = support_count + 1
        end do
        support_local(p,rank_index) = real(support_count, rprec)
    end do
    end do
    end do
end do

#ifdef PPMPI
call mpi_allreduce(support_local, support_global, size(support_local),        &
    mpi_rprec, mpi_sum, comm, ierr)
call mpi_allreduce(point_owner_local, point_owner, npoints, MPI_INTEGER,      &
    MPI_MAX, comm, ierr)
#else
support_global = support_local
point_owner = point_owner_local
#endif

do i = 1, numberOfTurbines
    if (turbineArray(i) % master >= 0 .and.                                  &
        turbineArray(i) % master < nrank) then
        current_turbines(turbineArray(i) % master + 1) =                     &
            current_turbines(turbineArray(i) % master + 1) + 1.0_rprec
        current_control(turbineArray(i) % master + 1) =                      &
            current_control(turbineArray(i) % master + 1) + control_cost
    endif
end do

do p = 1, npoints
    if (point_owner(p) >= 1 .and. point_owner(p) <= nrank) then
        current_points(point_owner(p)) = current_points(point_owner(p)) +     &
            1.0_rprec
    endif
    support_total = 0.0_rprec
    touched_ranks = 0
    do r = 1, nrank
        current_support(r) = current_support(r) + support_global(p,r)
        current_scatter(r) = current_scatter(r) + support_global(p,r)
        support_total = support_total + support_global(p,r)
        if (support_global(p,r) > 0.0_rprec) touched_ranks = touched_ranks + 1
    end do
    point_sort_cost(p) = support_total + alpha * real(touched_ranks, rprec)
end do

current_cost = current_control + current_support + current_points *           &
    24.0_rprec

planned_cost = current_control
used = .false.
do rr = 1, npoints
    selected = 0
    best_score = -1.0_rprec
    do p = 1, npoints
        if (.not. used(p)) then
            if (selected == 0 .or. point_sort_cost(p) > best_score) then
                selected = p
                best_score = point_sort_cost(p)
            endif
        endif
    end do
    if (selected == 0) exit
    used(selected) = .true.

    support_total = sum(support_global(selected,:))
    touched_ranks = 0
    do r = 1, nrank
        if (support_global(selected,r) > 0.0_rprec) touched_ranks =           &
            touched_ranks + 1
    end do

    best_rank = 1
    best_score = huge(1.0_rprec)
    do r = 1, nrank
        remote_cells = support_total - support_global(selected,r)
        remote_ranks = touched_ranks
        if (support_global(selected,r) > 0.0_rprec) remote_ranks =            &
            remote_ranks - 1
        sample_bytes = 0.0_rprec
        if (point_owner(selected) > 0 .and. point_owner(selected) /= r) then
            sample_bytes = 24.0_rprec * bytes_per_real
        endif
        scatter_bytes = remote_cells * 3.0_rprec * bytes_per_real
        comm_words = (sample_bytes + scatter_bytes) /                        &
            max(bytes_per_real, 1.0_rprec)
        candidate_score = planned_cost(r) + lambda * comm_words
        if (candidate_score < best_score) then
            best_score = candidate_score
            best_rank = r
        endif
    end do

    remote_cells = support_total - support_global(selected,best_rank)
    remote_ranks = touched_ranks
    if (support_global(selected,best_rank) > 0.0_rprec) remote_ranks =        &
        remote_ranks - 1
    sample_bytes = 0.0_rprec
    if (point_owner(selected) > 0 .and. point_owner(selected) /= best_rank)   &
        sample_bytes = 24.0_rprec * bytes_per_real
    scatter_bytes = remote_cells * 3.0_rprec * bytes_per_real
    comm_words = (sample_bytes + scatter_bytes) / max(bytes_per_real,        &
        1.0_rprec)
    point_cost = support_total + 24.0_rprec + alpha *                        &
        real(remote_ranks, rprec)

    point_plan(selected) = best_rank
    planned_points(best_rank) = planned_points(best_rank) + 1.0_rprec
    planned_support(best_rank) = planned_support(best_rank) + support_total
    planned_scatter(best_rank) = planned_scatter(best_rank) + support_total
    planned_sample_bytes(best_rank) = planned_sample_bytes(best_rank) +       &
        sample_bytes
    planned_scatter_bytes(best_rank) = planned_scatter_bytes(best_rank) +     &
        scatter_bytes
    planned_cost(best_rank) = planned_cost(best_rank) + point_cost +          &
        lambda * comm_words
    turbine_plan_points(point_turb(selected),best_rank) =                    &
        turbine_plan_points(point_turb(selected),best_rank) + 1.0_rprec
    turbine_plan_support(point_turb(selected),best_rank) =                   &
        turbine_plan_support(point_turb(selected),best_rank) + support_total
end do

if (atm_point_owner_lb_enabled()) then
    atm_lb_sample_send_slot = 0
    atm_lb_sample_recv_slot = 0
    atm_lb_force_send_slot = 0
    atm_lb_force_recv_slot = 0
    atm_lb_sample_send_count = 0
    atm_lb_sample_recv_count = 0
    atm_lb_force_send_count = 0
    atm_lb_force_recv_count = 0
    peer_rank = -1
    if (nrank == 2) peer_rank = 1 - coord

    do p = 1, npoints
        atm_lb_owner(p) = point_plan(p) - 1
        atm_lb_grid_owner(p) = point_owner(p) - 1
        if (nrank == 2 .and. point_owner(p) >= 1) then
            if (atm_lb_grid_owner(p) == coord .and. atm_lb_owner(p) ==        &
                peer_rank) then
                atm_lb_sample_send_count = atm_lb_sample_send_count + 1
                atm_lb_sample_send_slot(p) = atm_lb_sample_send_count
            elseif (atm_lb_grid_owner(p) == peer_rank .and.                  &
                    atm_lb_owner(p) == coord) then
                atm_lb_sample_recv_count = atm_lb_sample_recv_count + 1
                atm_lb_sample_recv_slot(p) = atm_lb_sample_recv_count
            endif

            if (atm_lb_owner(p) == coord) then
                atm_lb_force_send_count = atm_lb_force_send_count + 1
                atm_lb_force_send_slot(p) = atm_lb_force_send_count
            elseif (atm_lb_owner(p) == peer_rank) then
                atm_lb_force_recv_count = atm_lb_force_recv_count + 1
                atm_lb_force_recv_slot(p) = atm_lb_force_recv_count
            endif
        endif
    enddo
    atm_lb_planned_points_rank = planned_points(rank_index)
    atm_lb_planned_support_rank = planned_support(rank_index)
    atm_lb_planned_scatter_rank = planned_scatter(rank_index)
    atm_lb_owner_ready = .true.
endif

#ifdef PPMPI
call mpi_barrier(comm, ierr)
#endif

if (coord == 0 .and. atm_lb_plan_only_enabled()) then
    write(*,'(a)') '==================================================='
    write(*,'(a)') 'ATM Load-Balance Planner (diagnostic only)'
    write(*,'(a,i0,a,i0,a,es12.4,a,es12.4,a,es12.4,a,i0)')                  &
        'ATM_LB_SETTINGS turbines=', numberOfTurbines, ' points=', npoints,  &
        ' alpha=', alpha, ' lambda=', lambda, ' control_cost=',              &
        control_cost, ' bytes_per_real=', env_bytes
    write(*,'(a)') 'ATM_LB_TABLE_A rank turbines_owned actuator_points_owned support_points scatter_ops estimated_cost'
    do r = 1, nrank
        write(*,'(a,i0,5(a,es16.8))') 'ATM_LB_TABLE_A rank=', r-1,           &
            ' turbines_owned=', current_turbines(r),                         &
            ' actuator_points_owned=', current_points(r),                    &
            ' support_points=', current_support(r),                          &
            ' scatter_ops=', current_scatter(r),                             &
            ' estimated_cost=', current_cost(r)
    end do

    write(*,'(a)') 'ATM_LB_TABLE_B rank planned_actuator_points ' //          &
        'planned_support_ops planned_scatter_ops remote_sample_bytes ' //     &
        'remote_scatter_bytes estimated_total_cost'
    do r = 1, nrank
        write(*,'(a,i0,6(a,es16.8))') 'ATM_LB_TABLE_B rank=', r-1,           &
            ' planned_actuator_points=', planned_points(r),                  &
            ' planned_support_ops=', planned_support(r),                     &
            ' planned_scatter_ops=', planned_scatter(r),                     &
            ' remote_sample_bytes=', planned_sample_bytes(r),                &
            ' remote_scatter_bytes=', planned_scatter_bytes(r),              &
            ' estimated_total_cost=', planned_cost(r)
    end do

    write(*,'(a)') 'ATM_LB_TABLE_C turbine current_master rank planned_points planned_support_ops'
    do i = 1, numberOfTurbines
        do r = 1, nrank
            write(*,'(a,i0,a,i0,a,i0,2(a,es16.8))')                         &
                'ATM_LB_TABLE_C turbine=', i, ' current_master=',            &
                turbineArray(i) % master, ' rank=', r-1,                    &
                ' planned_points=', turbine_plan_points(i,r),                &
                ' planned_support_ops=', turbine_plan_support(i,r)
        end do
    end do

    call atm_lb_print_metric('actuator_points', current_points,              &
        planned_points, nrank)
    call atm_lb_print_metric('support_ops', current_support,                 &
        planned_support, nrank)
    call atm_lb_print_metric('scatter_ops', current_scatter,                 &
        planned_scatter, nrank)
    call atm_lb_print_metric('estimated_total_cost', current_cost,           &
        planned_cost, nrank)

    if (atm_lb_point_detail_enabled()) then
        point_id = 0
        do p = 1, npoints
            point_id = point_id + 1
            touched_ranks = 0
            do r = 1, nrank
                if (support_global(p,r) > 0.0_rprec) touched_ranks =         &
                    touched_ranks + 1
            end do
            write(*,'(a,i0,5(a,i0),3(a,es16.8),5(a,i0),a,es16.8)')          &
                'ATM_LB_POINT turbine=', point_turb(p), ' point=', point_id,&
                ' m=', point_m(p), ' n=', point_n(p), ' q=', point_q(p),    &
                ' x=', point_x(p), ' y=', point_y(p), ' z=', point_z(p),    &
                ' support_ranks=', touched_ranks, ' current_owner=',        &
                point_owner(p)-1, ' master=',                                &
                turbineArray(point_turb(p)) % master, ' planned_owner=',    &
                point_plan(p)-1, ' support_total=', sum(support_global(p,:))
            do r = 1, nrank
                write(*,'(a,i0,a,i0,a,es16.8)')                             &
                    'ATM_LB_POINT_RANK point=', point_id, ' rank=', r-1,    &
                    ' support_cells=', support_global(p,r)
            end do
        end do
    endif
    write(*,'(a)') 'ATM_LB_SYNTHETIC status=skipped reason=actual input turbine list was available'
    write(*,'(a)') 'ATM_LB_PLAN_ONLY no force arrays or turbine ownership were modified'
    write(*,'(a)') '==================================================='
    flush(6)
end if

if (coord == 0 .and. atm_point_owner_lb_enabled() .and.                      &
    .not. atm_lb_setup_printed) then
    write(*,'(a,i0,a,i0,a)') 'ATM_POINT_OWNER_LB enabled points=',           &
        atm_lb_npoints, ' ranks=', nrank, ' fallback=LESGO_ATM_POINT_OWNER_LB=0'
    write(*,'(a,es16.8,a,es16.8)')                                           &
        'ATM_POINT_OWNER_LB planned_cost_max_min=',                          &
        atm_lb_safe_ratio(maxval(planned_cost), minval(planned_cost)),       &
        ' support_max_min=',                                                 &
        atm_lb_safe_ratio(maxval(planned_support), minval(planned_support))
    if (atm_point_owner_targeted_enabled()) then
        write(*,'(a,i0,a,i0,a,i0,a,i0)')                                     &
            'ATM_POINT_OWNER_TARGETED nproc2_counts sample_send=',           &
            atm_lb_sample_send_count, ' sample_recv=',                       &
            atm_lb_sample_recv_count, ' force_send=',                        &
            atm_lb_force_send_count, ' force_recv=', atm_lb_force_recv_count
        write(*,'(a)')                                                       &
            'ATM_POINT_OWNER_TARGETED fallback=LESGO_ATM_POINT_OWNER_TARGETED=0'
    endif
    flush(6)
    atm_lb_setup_printed = .true.
endif

#ifdef PPMPI
call mpi_barrier(comm, ierr)
#endif

contains

subroutine atm_lb_print_metric(metric_name, current_vals, planned_vals, nvals)
implicit none

character(*), intent(in) :: metric_name
integer, intent(in) :: nvals
real(rprec), intent(in) :: current_vals(nvals), planned_vals(nvals)
real(rprec) :: current_max, current_min, current_mean
real(rprec) :: planned_max, planned_min, planned_mean

current_max = maxval(current_vals)
current_min = minval(current_vals)
current_mean = sum(current_vals) / max(1, nvals)
planned_max = maxval(planned_vals)
planned_min = minval(planned_vals)
planned_mean = sum(planned_vals) / max(1, nvals)

write(*,'(a,a,4(a,es16.8))') 'ATM_LB_TABLE_D metric=', trim(metric_name),    &
    ' current_max_min=', atm_lb_safe_ratio(current_max, current_min),         &
    ' planned_max_min=', atm_lb_safe_ratio(planned_max, planned_min),         &
    ' current_max_mean=', atm_lb_safe_ratio(current_max, current_mean),       &
    ' planned_max_mean=', atm_lb_safe_ratio(planned_max, planned_mean)

end subroutine atm_lb_print_metric

end subroutine atm_lesgo_lb_plan
#endif

#ifdef PPLES_GPU
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_build_force_shadows()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Flatten the per-turbine force-field metadata into persistent device-resident
! arrays. The timestep convolution/apply kernels then avoid dereferencing
! managed derived-type components for the grid scatter path.
implicit none

integer :: m, n, off, totalUV, totalW

if (atm_shadow_ready) return

allocate(atm_shadow_offsetUV(numberOfTurbines+1))
allocate(atm_shadow_offsetW(numberOfTurbines+1))
atm_shadow_offsetUV(1) = 0
atm_shadow_offsetW(1) = 0

do m = 1, numberOfTurbines
    atm_shadow_offsetUV(m+1) = atm_shadow_offsetUV(m) + forceFieldUV(m) % c
    atm_shadow_offsetW(m+1) = atm_shadow_offsetW(m) + forceFieldW(m) % c
end do

totalUV = atm_shadow_offsetUV(numberOfTurbines+1)
totalW = atm_shadow_offsetW(numberOfTurbines+1)

allocate(atm_shadow_ijkUV(3,max(totalUV,1)))
allocate(atm_shadow_ijkW(3,max(totalW,1)))
allocate(atm_shadow_locUV(3,max(totalUV,1)))
allocate(atm_shadow_locW(3,max(totalW,1)))
allocate(atm_shadow_forceUV(3,max(totalUV,1)))
allocate(atm_shadow_forceW(3,max(totalW,1)))
atm_shadow_forceUV = 0._rprec
atm_shadow_forceW = 0._rprec

do m = 1, numberOfTurbines
    n = forceFieldUV(m) % c
    if (n > 0) then
        off = atm_shadow_offsetUV(m)
        atm_shadow_ijkUV(:,off+1:off+n) = forceFieldUV(m) % ijk(:,1:n)
        atm_shadow_locUV(:,off+1:off+n) = forceFieldUV(m) % location(:,1:n)
    end if

    n = forceFieldW(m) % c
    if (n > 0) then
        off = atm_shadow_offsetW(m)
        atm_shadow_ijkW(:,off+1:off+n) = forceFieldW(m) % ijk(:,1:n)
        atm_shadow_locW(:,off+1:off+n) = forceFieldW(m) % location(:,1:n)
    end if
end do

allocate(atm_shadow_offsetUV_d(size(atm_shadow_offsetUV)))
allocate(atm_shadow_offsetW_d(size(atm_shadow_offsetW)))
allocate(atm_shadow_ijkUV_d(3,max(totalUV,1)))
allocate(atm_shadow_ijkW_d(3,max(totalW,1)))
allocate(atm_shadow_locUV_d(3,max(totalUV,1)))
allocate(atm_shadow_locW_d(3,max(totalW,1)))
allocate(atm_shadow_forceUV_d(3,max(totalUV,1)))
allocate(atm_shadow_forceW_d(3,max(totalW,1)))

atm_shadow_offsetUV_d = atm_shadow_offsetUV
atm_shadow_offsetW_d = atm_shadow_offsetW
atm_shadow_ijkUV_d = atm_shadow_ijkUV
atm_shadow_ijkW_d = atm_shadow_ijkW
atm_shadow_locUV_d = atm_shadow_locUV
atm_shadow_locW_d = atm_shadow_locW
atm_shadow_forceUV_d = 0._rprec
atm_shadow_forceW_d = 0._rprec

atm_shadow_ready = .true.
if (coord == 0) then
    write(*,'(a,2(a,i0))') 'ATM force shadows active:',                       &
        ' uv_cells=', totalUV, ' w_cells=', totalW
    flush(6)
end if

end subroutine atm_lesgo_build_force_shadows

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_destroy_force_shadows()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if (.not. atm_shadow_ready) return

deallocate(atm_shadow_forceUV_d, atm_shadow_forceW_d)
deallocate(atm_shadow_locUV_d, atm_shadow_locW_d)
deallocate(atm_shadow_ijkUV_d, atm_shadow_ijkW_d)
deallocate(atm_shadow_offsetUV_d, atm_shadow_offsetW_d)

deallocate(atm_shadow_forceUV, atm_shadow_forceW)
deallocate(atm_shadow_locUV, atm_shadow_locW)
deallocate(atm_shadow_ijkUV, atm_shadow_ijkW)
deallocate(atm_shadow_offsetUV, atm_shadow_offsetW)
atm_shadow_ready = .false.

end subroutine atm_lesgo_destroy_force_shadows

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_build_blade_mirrors()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: i, j

if (atm_blade_mirror_ready) return

atm_max_m = 0
atm_max_n = 0
atm_max_q = 0
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    atm_max_m = max(atm_max_m, turbineModel(j) % numBl)
    atm_max_n = max(atm_max_n, turbineArray(i) % numAnnulusSections)
    atm_max_q = max(atm_max_q, turbineArray(i) % numBladePoints)
end do

allocate(atm_bladePoints_d(max(atm_max_m,1), max(atm_max_n,1),                &
    max(atm_max_q,1), 3, max(numberOfTurbines,1)))
allocate(atm_bladeForces_d(max(atm_max_m,1), max(atm_max_n,1),                &
    max(atm_max_q,1), 3, max(numberOfTurbines,1)))
atm_bladePoints_d = 0._rprec
atm_bladeForces_d = 0._rprec

do i = 1, numberOfTurbines
    call atm_sync_blade_points_to_device(i)
    call atm_sync_blade_forces_to_device(i)
end do

atm_blade_mirror_ready = .true.
if (coord == 0) then
    write(*,'(a,3(a,i0))') 'ATM blade mirrors active:',                       &
        ' max_blades=', atm_max_m, ' max_sections=', atm_max_n,                &
        ' max_points=', atm_max_q
    flush(6)
end if

end subroutine atm_lesgo_build_blade_mirrors

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_destroy_blade_mirrors()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if (.not. atm_blade_mirror_ready) return
deallocate(atm_bladeForces_d, atm_bladePoints_d)
atm_blade_mirror_ready = .false.
atm_max_m = 0
atm_max_n = 0
atm_max_q = 0

end subroutine atm_lesgo_destroy_blade_mirrors

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_sync_blade_points_to_device(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, mmend, nnend, qqend

if (.not. allocated(atm_bladePoints_d)) return
j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
atm_bladePoints_d(1:mmend,1:nnend,1:qqend,1:3,i) =                           &
    turbineArray(i) % bladePoints(1:mmend,1:nnend,1:qqend,1:3)

end subroutine atm_sync_blade_points_to_device

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_sync_blade_forces_to_device(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j, mmend, nnend, qqend

if (.not. allocated(atm_bladeForces_d)) return
j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
atm_bladeForces_d(1:mmend,1:nnend,1:qqend,1:3,i) =                            &
    turbineArray(i) % bladeForces(1:mmend,1:nnend,1:qqend,1:3)

end subroutine atm_sync_blade_forces_to_device
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_findCells (m)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine finds all the cells that surround the turbines

! The awkward if statements are to only consider points in front and behind
! the turbine without having to

implicit none

! The turbine number
integer, intent(in) :: m

! Counter to establish number of points which are influenced by body forces
integer :: cUV, cW  ! Counters for number of points affected on UV and W grids
integer :: i, j, k

! Vector used to store x, y, z locations
real(rprec), dimension(3) :: vector_point
! These are the pointers to the grid arrays
real(rprec), pointer, dimension(:) :: x,y,z,zw

! Variables for MPI implementation
#ifdef PPMPI
integer :: base_group ! The base group from comm --> MPI_COMM_WORLD (all processors)
integer :: local_group  ! The local group of processors
integer :: member !  (1 or 0) yes or no
integer :: num_of_members  ! total number of members
#endif

! List of all the cores that belong to this turbine
! This variable gets allocated for each turbine
integer, allocatable, dimension(:) :: ls_of_cores

nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! Initialize internal counter to zero
forceFieldUV(m) % c = 0

! This will find all the locations that are influenced by each turbine
! It depends on a sphere centered on the rotor that extends beyond the blades
cUV=0  ! Initialize conuter
cW=0  ! Initialize conuter
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i

            ! Take into account the UV grid
            vector_point(3)=z(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cUV=cUV+1
!~ endif

                end if
                ! Take into account the W grid
                vector_point(3)=zw(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cW=cW+1
!~ endif
                end if
        enddo
    enddo
enddo

! Allocate space for the force fields in UV and W grids
forceFieldUV(m) % c = cUV  ! Counter
allocate(forceFieldUV(m) % force(3,cUV))
allocate(forceFieldUV(m) % location(3,cUV))
allocate(forceFieldUV(m) % ijk(3,cUV))

forceFieldW(m) % c = cW  ! Counter
allocate(forceFieldW(m) % force(3,cW))
allocate(forceFieldW(m) % location(3,cW))
allocate(forceFieldW(m) % ijk(3,cW))

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif
write(*,*) 'Number of cells being affected by ATM in turbine', m,              &
           ' cUV, cW = ', cUV, cW
#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

cUV=0
cW=0
! Run the same loop and save all variables
! The forceField arrays include all the forces which affect the domain
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i
            vector_point(3)=z(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cUV=cUV+1
                    forceFieldUV(m) % ijk(1,cUV) = i
                    forceFieldUV(m) % ijk(2,cUV) = j
                    forceFieldUV(m) % ijk(3,cUV) = k
                    forceFieldUV(m) % location(1:3,cUV) = vector_point(1:3)
                    forceFieldUV(m) % force(1:3,cUV) = 0_rprec
                endif
!~ endif
            vector_point(3)=zw(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cW=cW+1
                    forceFieldW(m) % ijk(1,cW) = i
                    forceFieldW(m) % ijk(2,cW) = j
                    forceFieldW(m) % ijk(3,cW) = k
                    forceFieldW(m) % location(1:3,cW) = vector_point(1:3)
                    forceFieldW(m) % force(:,cW) = 0_rprec
                endif
!~ endif
        enddo
    enddo
enddo


! MPI distribution
! This will create new communicator for each turbine
#ifdef PPMPI

! Store the base group from the global communicator mpi_comm_world
call MPI_COMM_GROUP(comm, base_group, ierr)

! Assign member
member = 0
! Flag to know if this turbine is operating or not
turbineArray(m) % operate = .FALSE.

! Assign proper values if turbine affects processors in this region
if (cUV > 0 .or. cW >0) then
member = 1
turbineArray(m) % operate = .TRUE.
endif

! Find the total number of processors for each turbine
call mpi_allreduce(member, num_of_members, 1, MPI_INTEGER , MPI_SUM, comm, ierr)

if (turbineArray(m) % operate) then
! Find the master processor for each turbine
    call mpi_allreduce(coord, turbineArray(m) % master, 1, MPI_INTEGER ,       &
                       MPI_MIN, comm, ierr)
else
    ! This is bogus since nz will always be less than number of processors
    ! This is done to ensure that the master is part of the processors
    ! that hold the turbine model
    call mpi_allreduce(nz_tot, turbineArray(m) % master, 1, MPI_INTEGER ,       &
                       MPI_MIN, comm, ierr)
endif

allocate(ls_of_cores(num_of_members))
ls_of_cores(1) = turbineArray(m) % master

! Notice this list is valid only for decomposition in 1 direction
do i = 2, num_of_members
    ls_of_cores(i) = ls_of_cores(i-1) + 1
enddo

! Write if this processor is the master
if (coord == turbineArray(m) % master) then
    write(*,*) 'Master for turbine',m, 'is processor', turbineArray(m) % master
endif

! Create the new communicator and group for this turbine
call MPI_GROUP_INCL(base_group, num_of_members, ls_of_cores, local_group, ierr)
call MPI_COMM_CREATE(comm, local_group, turbineArray(m) % TURBINE_COMM_WORLD, ierr)

if (turbineArray(m) % operate) then
    write(*,*) 'Processor', coord, 'has elements in turbine', m
else
    write(*,*) 'Processor', coord, 'does NOT have elements in turbine', m
endif

    call mpi_barrier( comm, ierr )

#endif

end subroutine atm_lesgo_findCells

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_forcing (phase)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutines calls the update function from the ATM Library
! and calculates the body forces needed in the domain
!
! Optional phase split (explicit-residency overlap): phase=1 runs the w_uv
! interpolation, blade update, velocity sampling, host blade-force model and
! MPI gather, then returns; phase=2 skips those and runs the convolution,
! Cl correction, apply and output. Called without phase, everything runs in
! order exactly as before. The split lets main.f90 run phase 1's ~30 ms of
! host work while the SGS/convection kernels are still queued on the GPU.
implicit none

integer, intent(in), optional :: phase
integer :: ph
integer :: i
logical :: atm_output_step
logical :: atm_lb_active = .false.   ! always declared; LB activation set only under ENABLE_CUDA
#ifdef ENABLE_CUDA
logical :: atm_diag_active
logical :: atm_lb_candidate, atm_lb_auto_probe_active
real(rprec) :: atm_diag_t0, atm_diag_t1
real(rprec) :: atm_auto_elapsed, atm_auto_t0, atm_auto_t1
#endif

!~ real(rprec) :: integrateNacelleForce, totForce
!~ integer :: c

!~ type(clock_t) :: myClock

ph = 0
if (present(phase)) ph = phase

if (ph /= 2) then

atm_forcing_calls = atm_forcing_calls + 1
#ifdef ENABLE_CUDA
atm_diag_active = atm_diag_timing_enabled()
! Point-owner LB uses a reduced at-point force path.  Keep structure off here
! until that path carries the structural force/moment state with parity tests.
atm_lb_candidate = atm_point_owner_lb_enabled() .and. atm_lb_owner_ready .and.&
    atm_point_owner_lb_supported() .and. .not. atm_full_gather_required() .and.&
    .not. atm_structure_enabled()
atm_lb_active = atm_lb_auto_use_lb_for_call(atm_lb_candidate)
atm_lb_auto_probe_active = atm_lb_candidate .and.                            &
    atm_lb_auto_select_enabled() .and. .not. atm_lb_auto_decision_made
atm_auto_elapsed = 0._rprec
if (atm_lb_auto_probe_active .and. coord == 0 .and.                         &
    .not. atm_lb_auto_printed) then
    write(*,'(a,i0,a,es12.4,a)') 'ATM_LB_AUTO_SELECT enabled probe_steps=',  &
        atm_lb_auto_probe_steps(), ' margin=', atm_lb_auto_margin(),         &
        ' fallback=LESGO_ATM_LB_AUTO_SELECT=0'
    flush(6)
    atm_lb_auto_printed = .true.
endif
if (atm_diag_active) then
    call atm_diag_event_start(atm_forcing_calls)
#ifdef PPMPI
    atm_diag_t0 = mpi_wtime()
#else
    call cpu_time(atm_diag_t0)
#endif
end if
#endif

!~ call myClock % start()
! Get the velocity from w onto the uv grid
call atm_clock_interp_w%start()
! atm_prepare_direct_w has PPLES_GPU paths and allocates/fills the direct-w
! halo used by the sampling kernels, so the runtime gate is sufficient here.
if (atm_direct_w_enabled()) then
    call atm_prepare_direct_w()
else
    call atm_interp_w_to_uv()
end if
call atm_clock_interp_w%stop()
atm_time_interp_w = atm_time_interp_w + atm_clock_interp_w%time


! Update the blade positions based on the time-step
! Time needs to be dimensionalized
! All processors carry the blade points
!~ call myCock%start_time();
!~ call atm_update(dt*z_i/u_star)

! Loop through all turbines and rotate the blades
call atm_clock_update%start()
do i = 1, numberOfTurbines
    ! If statement is for running code only with the processors on that turbine
        if (turbineArray(i) % operate) then
            ! Time is dimensionalize using velocity and length scale
            call atm_update(i, dt*z_i/u_star, total_time*z_i/u_star)
            #ifdef PPLES_GPU
            if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready) then
                call atm_sync_blade_points_to_device(i)
            endif
            #endif
        endif
    enddo
call atm_clock_update%stop()
atm_time_update = atm_time_update + atm_clock_update%time

!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Update ', myClock % time

end if   ! ph /= 2  (head of phase 1)

! Only calculate new forces if interval is correct
if ( mod(jt_total-1, updateInterval) == 0) then

    if (ph /= 2) then

#ifdef PPLES_GPU
    ! Batched device sampling for all atPoint turbines: one kernel + one D2H
    ! per force step instead of a data region + kernel + sync per turbine.
    ! Reads only bladePoints (already rotated by atm_update above) and the
    ! device u/v/w_uv, so it is safe to run before the reset loop.
    call atm_clock_sample%start()
    call atm_batch_sample_velocity_gpu()
    call atm_clock_sample%stop()
    atm_time_sample = atm_time_sample + atm_clock_sample%time
#endif

    ! Establish all turbine properties as zero
    ! This is essential for paralelization
    do i=1,numberOfTurbines
        call atm_clock_reset%start()
        turbineArray(i) % torqueRotor = 0._rprec
        turbineArray(i) % thrust = 0._rprec
        turbineArray(i) % nacelleForce = 0._rprec
        turbineArray(i) % VelNacelle_sampled = 0._rprec
        turbineArray(i) % VelNacelle_corrected = 0._rprec
#ifdef ENABLE_CUDA
        if (atm_reset_cuda_enabled()) then
            call atm_lesgo_reset_turbine_gpu(i)
        else
#endif
        turbineArray(i) % bladeForces = 0._rprec
        turbineArray(i) % integratedBladeForces = 0._rprec
        turbineArray(i) % alpha = 0._rprec
        turbineArray(i) % Cd = 0._rprec
        turbineArray(i) % Cm = 0._rprec
        turbineArray(i) % Cl = 0._rprec
        turbineArray(i) % Cl_b = 0._rprec
        turbineArray(i) % G = 0._rprec
        turbineArray(i) % lift = 0._rprec
        turbineArray(i) % drag = 0._rprec
        turbineArray(i) % Vmag = 0._rprec
        turbineArray(i) % windVectors = 0._rprec
        turbineArray(i) % induction_a = 0._rprec
        turbineArray(i) % u_infinity = 0._rprec
        turbineArray(i) % bladeAlignedVectors = 0._rprec
        turbineArray(i) % axialForce = 0._rprec
        turbineArray(i) % tangentialForce = 0._rprec
        turbineArray(i) % pitchingMoment = 0._rprec
#ifdef ENABLE_CUDA
        endif
#endif

        ! If statement is for running code only if grid points affected are in
        ! this processor. If not, no code is executed at all.
!~         if (forceFieldUV(i) % c .gt. 0 .or. forceFieldW(i) % c .gt. 0) then
        if (turbineArray(i) % operate) then
            ! The applied components are overwritten by the convolution step.
            ! Avoid clearing the full force-field arrays here; they are large
            ! enough that this host-side touch is visible in managed-memory runs.
        endif
        call atm_clock_reset%stop()
        atm_time_reset = atm_time_reset + atm_clock_reset%time

        if (.not. atm_lb_active .and. turbineArray(i) % operate) then
            ! Calculate forces for all turbines
            call atm_clock_force%start()
#ifdef ENABLE_CUDA
            if (atm_lb_auto_probe_active) then
#ifdef PPMPI
                atm_auto_t0 = mpi_wtime()
#else
                call cpu_time(atm_auto_t0)
#endif
            endif
#endif
            call atm_lesgo_force(i)
#ifdef ENABLE_CUDA
            if (atm_lb_auto_probe_active) then
#ifdef PPMPI
                atm_auto_t1 = mpi_wtime()
#else
                call cpu_time(atm_auto_t1)
#endif
                atm_auto_elapsed = atm_auto_elapsed +                        &
                    (atm_auto_t1 - atm_auto_t0)
            endif
#endif
            call atm_clock_force%stop()
            atm_time_force = atm_time_force + atm_clock_force%time

        endif

    enddo
    if (atm_lb_active) then
        call atm_clock_force%start()
#ifdef ENABLE_CUDA
        if (atm_lb_auto_probe_active) then
#ifdef PPMPI
            atm_auto_t0 = mpi_wtime()
#else
            call cpu_time(atm_auto_t0)
#endif
        endif
#endif
        call atm_point_owner_lb_force()
#ifdef ENABLE_CUDA
        if (atm_lb_auto_probe_active) then
#ifdef PPMPI
            atm_auto_t1 = mpi_wtime()
#else
            call cpu_time(atm_auto_t1)
#endif
            atm_auto_elapsed = atm_auto_elapsed + (atm_auto_t1 - atm_auto_t0)
        endif
#endif
        call atm_clock_force%stop()
        atm_time_force = atm_time_force + atm_clock_force%time
    endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Forces ', myClock % time


!~  call myClock % start()
! This will gather all the blade forces from all processors
#ifdef PPMPI
    ! This will gather all values used in MPI
!~     call mpi_barrier( MPI_COMM_WORLD, ierr )
    if (nproc > 1) then
        call atm_clock_gather%start()
#ifdef ENABLE_CUDA
        if (atm_lb_auto_probe_active) then
#ifdef PPMPI
            atm_auto_t0 = mpi_wtime()
#else
            call cpu_time(atm_auto_t0)
#endif
        endif
#endif
        if (atm_lb_active) then
            if (atm_point_owner_targeted_supported()) then
                call atm_point_owner_lb_gather_targeted()
            else
                call atm_point_owner_lb_gather()
            endif
        else
            call atm_lesgo_mpi_gather()
        endif
#ifdef ENABLE_CUDA
        if (atm_lb_auto_probe_active) then
#ifdef PPMPI
            atm_auto_t1 = mpi_wtime()
#else
            call cpu_time(atm_auto_t1)
#endif
            atm_auto_elapsed = atm_auto_elapsed + (atm_auto_t1 - atm_auto_t0)
        endif
#endif
        call atm_clock_gather%stop()
        atm_time_gather = atm_time_gather + atm_clock_gather%time
    endif
!~     call mpi_barrier( MPI_COMM_WORLD, ierr )

#endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  MPI Gather ', myClock % time

#ifdef ENABLE_CUDA
    if (atm_lb_auto_probe_active) then
        call atm_lb_auto_record(atm_lb_active, atm_auto_elapsed)
    endif
#endif

    end if   ! ph /= 2  (sampling + blade force + gather)

    if (ph == 1) return

!~  call myClock % start()
    #ifdef PPLES_GPU
    if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready .and.         &
        (atm_lb_active .or. .not. atm_packed_gather_enabled() .or.             &
         .not. atm_gpu_packed_gather_enabled() .or.                            &
         .not. atm_slim_gather_enabled() .or. atm_full_gather_required())) then
        do i = 1, numberOfTurbines
            if (turbineArray(i) % operate) then
                call atm_sync_blade_forces_to_device(i)
            endif
        enddo
    endif
    #endif

#ifdef PPLES_GPU
    ! Batched convolution for all atPoint turbines: two kernels (UV + W grid)
    ! over the concatenated force-field cell lists, scattering straight into
    ! the device fxa/fya/fza. Replaces 60 per-turbine data regions + kernels.
    ! Non-atPoint (Spalart) turbines keep the per-turbine path below.
    call atm_clock_convolve%start()
    call atm_batch_convolute_force_gpu()
    call atm_clock_convolve%stop()
    atm_time_convolve = atm_time_convolve + atm_clock_convolve%time

    ! Batched OpenACC Cl/tip correction updates the same induced-velocity
    ! state used by both rigid and structural turbine consumers.
    call atm_clock_clcorr%start()
    call atm_batch_cl_correction_gpu()
    call atm_clock_clcorr%stop()
    atm_time_clcorr = atm_time_clcorr + atm_clock_clcorr%time
#endif

    do i=1,numberOfTurbines
!~         if ( forceFieldUV(i) % c .gt. 0 .or. forceFieldW(i) % c .gt. 0) then

        ! Only perform is turbine is active in this processor
        if (turbineArray(i) % operate) then
#ifdef PPLES_GPU
            if (turbineArray(i) % sampling /= 'atPoint') then
#endif
            ! Convolute force onto the domain
            call atm_clock_convolve%start()
            call atm_lesgo_convolute_force(i)
            call atm_clock_convolve%stop()
            atm_time_convolve = atm_time_convolve + atm_clock_convolve%time
#ifdef PPLES_GPU
            endif
#endif

            ! Only do this if the correction is active
            if (turbineArray(i) % tipALMCorrection .eqv. .true.)  then
#ifdef PPLES_GPU
                if (turbineArray(i) % sampling /= 'atPoint') then
#endif
                ! Compute the correction for the Cl coefficient.  Structure-off
                ! and structure-on atPoint turbines are handled above by
                ! atm_batch_cl_correction_gpu().
                call atm_clock_clcorr%start()
                call atm_compute_cl_correction(i)
                call atm_clock_clcorr%stop()
                atm_time_clcorr = atm_time_clcorr + atm_clock_clcorr%time
#ifdef PPLES_GPU
                endif
#endif
            endif

        endif

!~         ! Sync the nacelle force
!~         integrateNacelleForce=0.
!~
!~         do c=1,forceFieldUV(i) % c
!~             if (turbineArray(i) % nacelle) then
!~                 integrateNacelleForce = integrateNacelleForce +  &
!~                     forceFieldUV(i) % force(1,c) * dx *dy * dz * z_i**2*u_star**2
!~             endif
!~         enddo


    enddo

#ifdef ENABLE_CUDA
    if (atm_lb_validate_enabled()) then
        if (atm_lb_active) then
            call atm_point_owner_lb_validate('point_owner_after_convolution')
        else
            call atm_point_owner_lb_validate('legacy_after_convolution')
        endif
    endif
#endif

!~         totForce=0.
!~         call mpi_allreduce( integrateNacelleForce,  totForce, 1,   &
!~                              mpi_rprec, mpi_sum, comm, ierr)

       !write(*,*) 'Integrated Nacelle Force is: ', integrateNacelleForce
!~         if (coord == 0) then
!~             write(*,*) 'Integrated Total Force is: ', totForce
!~         endif
endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Convolute force ', myClock % time

! Phase-1 calls stop here on non-update steps too (apply/output belong to
! phase 2 so they run exactly once per step).
if (ph == 1) return

    ! This will apply body forces onto the flow field if there are forces within
    ! this domain
!~  call myClock % start()
    call atm_clock_apply%start()
    call atm_lesgo_apply_force()
    call atm_clock_apply%stop()
    atm_time_apply = atm_time_apply + atm_clock_apply%time
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Apply force ', myClock % time

!!! Sync the integrated forces (used for debugging)
!do i=1,numberOfTurbines
!    j=turbineArray(i) % turbineTypeID ! The turbine type ID
!    ! Sync all the integrated blade forces
!    turbineArray(i) % bladeVectorDummy=turbineArray(i) % integratedBladeForces
!    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                   &
!                       turbineArray(i) % integratedBladeForces,              &
!                       size(turbineArray(i) % bladeVectorDummy),             &
!                       mpi_rprec, mpi_sum, comm, ierr)


!    if (coord==0) then

!    do q=1, turbineArray(i) % numBladePoints
!        do n=1, turbineArray(i) % numAnnulusSections
!            do m=1, turbineModel(j) % numBl
!                write(*,*) 'blade ',m,'section ',q, 'force ratio', &
!                turbineArray(i) % integratedBladeForces(m,n,q,1) /  &
!                turbineArray(i) % bladeForces(m,n,q,1) , &
!                turbineArray(i) % integratedBladeForces(m,n,q,2) /  &
!                turbineArray(i) % bladeForces(m,n,q,2) , &
!                turbineArray(i) % integratedBladeForces(m,n,q,3) /  &
!                turbineArray(i) % bladeForces(m,n,q,3)
!            enddo
!        enddo
!    enddo
!    endif

!enddo

atm_output_step = outputInterval > 0 .and. mod(jt_total-1, outputInterval) == 0
if (atm_output_step) then
    do i=1, numberOfTurbines
        if (coord == turbineArray(i) % master) then
        !~  call myClock % start()

            call atm_clock_output%start()
            call atm_output(i, jt_total, total_time*z_i/u_star)
            call atm_clock_output%stop()
            atm_time_output = atm_time_output + atm_clock_output%time
        !~     call myClock % stop()
        !~     write(*,*) 'coord ', coord, '  Output ', myClock % time
        endif
    enddo
endif

! Make sure all processors stop wait for the output to be completed
#ifdef PPMPI
if (atm_output_step .and. atm_output_barrier_enabled()) then
    call atm_clock_barrier%start()
    call mpi_barrier( comm, ierr )
    call atm_clock_barrier%stop()
    atm_time_barrier = atm_time_barrier + atm_clock_barrier%time
endif
#endif

#ifdef ENABLE_CUDA
if (atm_diag_active) then
    call atm_diag_event_stop(atm_forcing_calls)
#ifdef PPMPI
    atm_diag_t1 = mpi_wtime()
#else
    call cpu_time(atm_diag_t1)
#endif
    write(*,'(a,i0,a,i0,3(a,es14.6),a,i0)') 'ATM_TIMING call=',             &
        atm_forcing_calls, ' rank=', coord, ' cpu_wall=',                    &
        atm_diag_t1 - atm_diag_t0, ' gpu_event=', -1.0_rprec,                &
        ' sync=', -1.0_rprec, ' kernels=', -1
end if
#endif

end subroutine atm_lesgo_forcing


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Complie this subroutines only if MPI will be used
#ifdef PPMPI

subroutine atm_lesgo_mpi_gather()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will gather the necessary outputs from the turbine models
! so all processors have acces to it. This is done by means of all reduce SUM
implicit none
integer :: i
real(rprec) :: torqueRotor, thrust, VelNacelle_sampled, VelNacelle_corrected
real(rprec), dimension(3) :: nacelleForce

! Pointer for MPI communicator
integer, pointer :: TURBINE_COMMUNICATOR

if (nproc <= 1) return

#ifdef ENABLE_CUDA
if (atm_packed_gather_enabled()) then
    call atm_lesgo_mpi_gather_packed()
    return
endif
#else
! Use the packed gather - one allreduce per turbine instead of many per-array
! reductions.  Structure-on runs add Cm/pitchingMoment to the packed payload.
call atm_lesgo_mpi_gather_packed()
return
#endif

do i=1,numberOfTurbines

    ! Only do MPI sums if processors are operating in this turbine
    if (turbineArray(i) % operate) then

        TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

        turbineArray(i) % bladeVectorDummy = turbineArray(i) % bladeForces
        ! Sync all the blade forces
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeForces,                      &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync bladeAlignedVectors
        turbineArray(i) % bladeVectorDummy =                                   &
                              turbineArray(i) % bladeAlignedVectors(:,:,:,1,:)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeAlignedVectors(:,:,:,1,:),   &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        turbineArray(i) % bladeVectorDummy =                                   &
                              turbineArray(i) % bladeAlignedVectors(:,:,:,2,:)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeAlignedVectors(:,:,:,2,:),   &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        turbineArray(i) % bladeVectorDummy =                                   &
                              turbineArray(i) % bladeAlignedVectors(:,:,:,3,:)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % bladeAlignedVectors(:,:,:,3,:),   &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)


        ! Sync alpha
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % alpha
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % alpha,                            &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        ! Sync lift
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % lift
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % lift,                             &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        ! Sync drag
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % drag
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % drag,                             &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        ! Sync Cl
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cl
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % Cl,                               &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync Cd
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cd
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % Cd,                               &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        if (atm_structure_enabled()) then
            ! Sync Cm and pitching moment for the structural torsion solve.
            turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cm
            call mpi_allreduce(turbineArray(i) % bladeScalarDummy,             &
                               turbineArray(i) % Cm,                           &
                               size(turbineArray(i) % bladeScalarDummy),       &
                               mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

            turbineArray(i) % bladeScalarDummy =                               &
                turbineArray(i) % pitchingMoment
            call mpi_allreduce(turbineArray(i) % bladeScalarDummy,             &
                               turbineArray(i) % pitchingMoment,               &
                               size(turbineArray(i) % bladeScalarDummy),       &
                               mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)
        endif

        ! Sync Vmag
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % Vmag
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % Vmag,                             &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync axialForce
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % axialForce
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % axialForce,                       &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync tangentialForce
        turbineArray(i) % bladeScalarDummy = turbineArray(i) % tangentialForce
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % tangentialForce,                                 &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync wind Vectors (Vaxial, Vtangential, Vradial)
        turbineArray(i) % bladeVectorDummy = turbineArray(i) %                 &
                                             windVectors(:,:,:,1:3)
        call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                 &
                           turbineArray(i) % windVectors(:,:,:,1:3),                                 &
                           size(turbineArray(i) % bladeVectorDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync induction factor
        turbineArray(i) % bladeScalarDummy = turbineArray(i) %                 &
                                             induction_a(:,:,:)
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % induction_a(:,:,:),                                 &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync u infinity
        turbineArray(i) % bladeScalarDummy = turbineArray(i) %                 &
                                             u_infinity(:,:,:)
        call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                 &
                           turbineArray(i) % u_infinity(:,:,:),                                 &
                           size(turbineArray(i) % bladeScalarDummy),           &
                           mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Store the torqueRotor.
        ! Needs to be a different variable in order to do MPI Sum
        torqueRotor=turbineArray(i) % torqueRotor
        thrust=turbineArray(i) % thrust
        nacelleForce=turbineArray(i) % nacelleForce
        VelNacelle_sampled=turbineArray(i) % VelNacelle_sampled
        VelNacelle_corrected=turbineArray(i) % VelNacelle_corrected

        ! Sum all the individual torqueRotor from different blade points
        call mpi_allreduce( torqueRotor, turbineArray(i) % torqueRotor,        &
                           1, mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sum all the individual thrust from different blade points
        call mpi_allreduce( thrust, turbineArray(i) % thrust,                  &
                           1, mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync the nacelle force
        call mpi_allreduce( nacelleForce, turbineArray(i) % nacelleForce,      &
                           size(turbineArray(i) % nacelleForce), mpi_rprec,    &
                                mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync the nacelle sampled velocity
        call mpi_allreduce( VelNacelle_sampled,                                &
                               turbineArray(i) % VelNacelle_sampled,  1,       &
                                mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

        ! Sync the nacelle corrected velocity
        call mpi_allreduce( VelNacelle_corrected,                              &
                               turbineArray(i) % VelNacelle_corrected, 1,      &
                                mpi_rprec, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    endif
enddo
end subroutine atm_lesgo_mpi_gather

subroutine atm_lesgo_mpi_gather_packed()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Consolidate the ATM per-turbine reductions into one allreduce to reduce
! latency on multi-GPU runs.
implicit none

integer :: i, nitem, npack, pos
real(rprec), allocatable, save :: packed_send(:), packed_recv(:)
integer, pointer :: TURBINE_COMMUNICATOR
logical :: struct_active

#ifdef ENABLE_CUDA
if (atm_gpu_packed_gather_enabled()) then
    call atm_lesgo_mpi_gather_packed_gpu()
    return
endif
#endif

struct_active = atm_structure_enabled()

do i=1,numberOfTurbines

    if (turbineArray(i) % operate) then

        TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

        npack = size(turbineArray(i) % bladeForces) +                         &
                size(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:)) +      &
                size(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:)) +      &
                size(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:)) +      &
                size(turbineArray(i) % alpha) + size(turbineArray(i) % lift) +&
                size(turbineArray(i) % drag) + size(turbineArray(i) % Cl) +   &
                size(turbineArray(i) % Cd) + size(turbineArray(i) % Vmag) +   &
                size(turbineArray(i) % axialForce) +                         &
                size(turbineArray(i) % tangentialForce) +                    &
                size(turbineArray(i) % windVectors(:,:,:,1:3)) +             &
                size(turbineArray(i) % induction_a(:,:,:)) +                 &
                size(turbineArray(i) % u_infinity(:,:,:)) + 7
        if (struct_active) then
            npack = npack + size(turbineArray(i) % Cm) +                      &
                size(turbineArray(i) % pitchingMoment)
        endif

        if (allocated(packed_send)) then
            if (size(packed_send) /= npack) then
                deallocate(packed_send, packed_recv)
            endif
        endif
        if (.not. allocated(packed_send)) then
            allocate(packed_send(npack), packed_recv(npack))
        endif

        pos = 1

        nitem = size(turbineArray(i) % bladeForces)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % bladeForces,  &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % bladeAlignedVectors(:,:,:,1,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % bladeAlignedVectors(:,:,:,2,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % bladeAlignedVectors(:,:,:,3,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % alpha)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % alpha,        &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % lift)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % lift,         &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % drag)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % drag,         &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cl)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cl,           &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cd)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cd,           &
            (/ nitem /))
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cm,       &
                (/ nitem /))
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            packed_send(pos:pos+nitem-1) = reshape(                           &
                turbineArray(i) % pitchingMoment, (/ nitem /))
            pos = pos + nitem
        endif

        nitem = size(turbineArray(i) % Vmag)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Vmag,         &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % axialForce)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % axialForce,   &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % tangentialForce)
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % tangentialForce, (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % windVectors(:,:,:,1:3))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % windVectors(:,:,:,1:3), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % induction_a(:,:,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % induction_a(:,:,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % u_infinity(:,:,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % u_infinity(:,:,:), (/ nitem /))
        pos = pos + nitem

        packed_send(pos) = turbineArray(i) % torqueRotor
        packed_send(pos+1) = turbineArray(i) % thrust
        packed_send(pos+2:pos+4) = turbineArray(i) % nacelleForce
        packed_send(pos+5) = turbineArray(i) % VelNacelle_sampled
        packed_send(pos+6) = turbineArray(i) % VelNacelle_corrected

        call mpi_allreduce(packed_send, packed_recv, npack, mpi_rprec,         &
                           mpi_sum, TURBINE_COMMUNICATOR, ierr)

        pos = 1

        nitem = size(turbineArray(i) % bladeForces)
        turbineArray(i) % bladeForces = reshape(packed_recv(pos:pos+nitem-1),  &
            shape(turbineArray(i) % bladeForces))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:))
        turbineArray(i) % bladeAlignedVectors(:,:,:,1,:) = reshape(           &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:))
        turbineArray(i) % bladeAlignedVectors(:,:,:,2,:) = reshape(           &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:))
        turbineArray(i) % bladeAlignedVectors(:,:,:,3,:) = reshape(           &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % alpha)
        turbineArray(i) % alpha = reshape(packed_recv(pos:pos+nitem-1),        &
            shape(turbineArray(i) % alpha))
        pos = pos + nitem

        nitem = size(turbineArray(i) % lift)
        turbineArray(i) % lift = reshape(packed_recv(pos:pos+nitem-1),         &
            shape(turbineArray(i) % lift))
        pos = pos + nitem

        nitem = size(turbineArray(i) % drag)
        turbineArray(i) % drag = reshape(packed_recv(pos:pos+nitem-1),         &
            shape(turbineArray(i) % drag))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cl)
        turbineArray(i) % Cl = reshape(packed_recv(pos:pos+nitem-1),           &
            shape(turbineArray(i) % Cl))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cd)
        turbineArray(i) % Cd = reshape(packed_recv(pos:pos+nitem-1),           &
            shape(turbineArray(i) % Cd))
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            turbineArray(i) % Cm = reshape(packed_recv(pos:pos+nitem-1),       &
                shape(turbineArray(i) % Cm))
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            turbineArray(i) % pitchingMoment = reshape(                       &
                packed_recv(pos:pos+nitem-1),                                 &
                shape(turbineArray(i) % pitchingMoment))
            pos = pos + nitem
        endif

        nitem = size(turbineArray(i) % Vmag)
        turbineArray(i) % Vmag = reshape(packed_recv(pos:pos+nitem-1),         &
            shape(turbineArray(i) % Vmag))
        pos = pos + nitem

        nitem = size(turbineArray(i) % axialForce)
        turbineArray(i) % axialForce = reshape(                               &
            packed_recv(pos:pos+nitem-1), shape(turbineArray(i) % axialForce))
        pos = pos + nitem

        nitem = size(turbineArray(i) % tangentialForce)
        turbineArray(i) % tangentialForce = reshape(                          &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % tangentialForce))
        pos = pos + nitem

        nitem = size(turbineArray(i) % windVectors(:,:,:,1:3))
        turbineArray(i) % windVectors(:,:,:,1:3) = reshape(                   &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % windVectors(:,:,:,1:3)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % induction_a(:,:,:))
        turbineArray(i) % induction_a(:,:,:) = reshape(                       &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % induction_a(:,:,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % u_infinity(:,:,:))
        turbineArray(i) % u_infinity(:,:,:) = reshape(                        &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % u_infinity(:,:,:)))
        pos = pos + nitem

        turbineArray(i) % torqueRotor = packed_recv(pos)
        turbineArray(i) % thrust = packed_recv(pos+1)
        turbineArray(i) % nacelleForce = packed_recv(pos+2:pos+4)
        turbineArray(i) % VelNacelle_sampled = packed_recv(pos+5)
        turbineArray(i) % VelNacelle_corrected = packed_recv(pos+6)

    endif

enddo

end subroutine atm_lesgo_mpi_gather_packed

#ifdef ENABLE_CUDA
subroutine atm_lesgo_mpi_gather_slim_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Timestep path only needs blade forces for convolution plus rotor/nacelle
! scalars for controls.  Diagnostics and CL correction keep using the full
! packed gather when needed.
implicit none

integer :: i, nitem, npack, pos
real(rprec), device, allocatable, save :: packed_send(:), packed_recv(:)
integer, pointer :: TURBINE_COMMUNICATOR
logical :: batch_done, struct_active

if (atm_batch_gather_enabled()) then
    call atm_lesgo_mpi_gather_slim_batch_gpu(batch_done)
    if (batch_done) return
end if

struct_active = atm_structure_enabled()

do i=1,numberOfTurbines

    if (turbineArray(i) % operate) then

        TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

        npack = size(turbineArray(i) % bladeForces) + 7
        if (struct_active) then
            npack = npack + size(turbineArray(i) % Cm) +                      &
                size(turbineArray(i) % pitchingMoment)
        endif

        if (allocated(packed_send)) then
            if (size(packed_send) /= npack) then
                deallocate(packed_send, packed_recv)
            endif
        endif
        if (.not. allocated(packed_send)) then
            allocate(packed_send(npack), packed_recv(npack))
        endif

        pos = 1
        nitem = size(turbineArray(i) % bladeForces)
#ifdef PPLES_GPU
        if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready) then
            call atm_pack_blade_forces_mirror(i, packed_send, pos)
        else
            call atm_pack_rank4(turbineArray(i) % bladeForces, packed_send, pos)
        endif
#else
        call atm_pack_rank4(turbineArray(i) % bladeForces, packed_send, pos)
#endif
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            call atm_pack_rank3(turbineArray(i) % Cm, packed_send, pos)
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            call atm_pack_rank3(turbineArray(i) % pitchingMoment, packed_send, &
                pos)
            pos = pos + nitem
        endif

        call atm_pack_gather_scalars(packed_send, pos,                         &
            turbineArray(i) % torqueRotor, turbineArray(i) % thrust,            &
            turbineArray(i) % nacelleForce(1), turbineArray(i) % nacelleForce(2),&
            turbineArray(i) % nacelleForce(3),                                  &
            turbineArray(i) % VelNacelle_sampled,                               &
            turbineArray(i) % VelNacelle_corrected)
        call atm_cuda_sync('ATM slim gather pack')

        call mpi_allreduce(packed_send, packed_recv, npack, mpi_rprec,          &
                           mpi_sum, TURBINE_COMMUNICATOR, ierr)

        pos = 1
        nitem = size(turbineArray(i) % bladeForces)
#ifdef PPLES_GPU
        if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready) then
            call atm_unpack_blade_forces_mirror(packed_recv, pos, i)
        else
            call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % bladeForces)
        endif
#else
        call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % bladeForces)
#endif
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cm)
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            call atm_unpack_rank3(packed_recv, pos,                            &
                turbineArray(i) % pitchingMoment)
            pos = pos + nitem
        endif

        turbineArray(i) % torqueRotor = packed_recv(pos)
        turbineArray(i) % thrust = packed_recv(pos+1)
        turbineArray(i) % nacelleForce = packed_recv(pos+2:pos+4)
        turbineArray(i) % VelNacelle_sampled = packed_recv(pos+5)
        turbineArray(i) % VelNacelle_corrected = packed_recv(pos+6)

        call atm_cuda_check('ATM slim gather unpack')

    endif

enddo

end subroutine atm_lesgo_mpi_gather_slim_gpu

subroutine atm_lesgo_mpi_gather_slim_batch_gpu(done)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batch the common case where every turbine communicator spans all MPI ranks.
! If any turbine is rank-local/subset-local, fall back to the per-turbine path.
implicit none

logical, intent(out) :: done
integer :: i, nitem, npack, pos, total_pack, comm_size
real(rprec), device, allocatable, save :: packed_send(:), packed_recv(:)
integer, pointer :: TURBINE_COMMUNICATOR
logical :: struct_active

done = .false.
if (nproc <= 1) return

struct_active = atm_structure_enabled()

total_pack = 0
do i=1,numberOfTurbines
    if (.not. turbineArray(i) % operate) return

    TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD
    call mpi_comm_size(TURBINE_COMMUNICATOR, comm_size, ierr)
    if (comm_size /= nproc) return

    total_pack = total_pack + size(turbineArray(i) % bladeForces) + 7
    if (struct_active) then
        total_pack = total_pack + size(turbineArray(i) % Cm) +                 &
            size(turbineArray(i) % pitchingMoment)
    endif
enddo

if (total_pack <= 0) return

if (allocated(packed_send)) then
    if (size(packed_send) /= total_pack) then
        deallocate(packed_send, packed_recv)
    endif
endif
if (.not. allocated(packed_send)) then
    allocate(packed_send(total_pack), packed_recv(total_pack))
endif

pos = 1
do i=1,numberOfTurbines
    nitem = size(turbineArray(i) % bladeForces)
#ifdef PPLES_GPU
    if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready) then
        call atm_pack_blade_forces_mirror(i, packed_send, pos)
    else
        call atm_pack_rank4(turbineArray(i) % bladeForces, packed_send, pos)
    endif
#else
    call atm_pack_rank4(turbineArray(i) % bladeForces, packed_send, pos)
#endif
    pos = pos + nitem

    if (struct_active) then
        nitem = size(turbineArray(i) % Cm)
        call atm_pack_rank3(turbineArray(i) % Cm, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % pitchingMoment)
        call atm_pack_rank3(turbineArray(i) % pitchingMoment, packed_send, pos)
        pos = pos + nitem
    endif

    call atm_pack_gather_scalars(packed_send, pos,                         &
        turbineArray(i) % torqueRotor, turbineArray(i) % thrust,            &
        turbineArray(i) % nacelleForce(1), turbineArray(i) % nacelleForce(2),&
        turbineArray(i) % nacelleForce(3),                                  &
        turbineArray(i) % VelNacelle_sampled,                               &
        turbineArray(i) % VelNacelle_corrected)
    pos = pos + 7
enddo

call atm_cuda_sync('ATM slim batch gather pack')
call mpi_allreduce(packed_send, packed_recv, total_pack, mpi_rprec, mpi_sum, &
                   comm, ierr)

pos = 1
do i=1,numberOfTurbines
    nitem = size(turbineArray(i) % bladeForces)
#ifdef PPLES_GPU
    if (atm_force_shadows_enabled() .and. atm_blade_mirror_ready) then
        call atm_unpack_blade_forces_mirror(packed_recv, pos, i)
    else
        call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % bladeForces)
    endif
#else
    call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % bladeForces)
#endif
    pos = pos + nitem

    if (struct_active) then
        nitem = size(turbineArray(i) % Cm)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cm)
        pos = pos + nitem

        nitem = size(turbineArray(i) % pitchingMoment)
        call atm_unpack_rank3(packed_recv, pos,                              &
            turbineArray(i) % pitchingMoment)
        pos = pos + nitem
    endif

    turbineArray(i) % torqueRotor = packed_recv(pos)
    turbineArray(i) % thrust = packed_recv(pos+1)
    turbineArray(i) % nacelleForce = packed_recv(pos+2:pos+4)
    turbineArray(i) % VelNacelle_sampled = packed_recv(pos+5)
    turbineArray(i) % VelNacelle_corrected = packed_recv(pos+6)
    pos = pos + 7
enddo

call atm_cuda_check('ATM slim batch gather unpack')
done = .true.

end subroutine atm_lesgo_mpi_gather_slim_batch_gpu

subroutine atm_lesgo_mpi_gather_packed_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! GPU pack/unpack variant for the consolidated ATM allreduce.
implicit none

integer :: i, nitem, npack, pos
real(rprec), device, allocatable, save :: packed_send(:), packed_recv(:)
integer, pointer :: TURBINE_COMMUNICATOR
logical :: struct_active

struct_active = atm_structure_enabled()

if (atm_slim_gather_enabled() .and. .not. atm_full_gather_required()) then
    call atm_lesgo_mpi_gather_slim_gpu()
    return
endif

do i=1,numberOfTurbines

    if (turbineArray(i) % operate) then

        TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

        npack = size(turbineArray(i) % bladeForces) +                         &
                size(turbineArray(i) % bladeAlignedVectors) +                 &
                size(turbineArray(i) % alpha) + size(turbineArray(i) % lift) +&
                size(turbineArray(i) % drag) + size(turbineArray(i) % Cl) +   &
                size(turbineArray(i) % Cd) + size(turbineArray(i) % Vmag) +   &
                size(turbineArray(i) % axialForce) +                         &
                size(turbineArray(i) % tangentialForce) +                    &
                size(turbineArray(i) % windVectors) +                        &
                size(turbineArray(i) % induction_a) +                        &
                size(turbineArray(i) % u_infinity) + 7
        if (struct_active) then
            npack = npack + size(turbineArray(i) % Cm) +                      &
                size(turbineArray(i) % pitchingMoment)
        endif

        if (allocated(packed_send)) then
            if (size(packed_send) /= npack) then
                deallocate(packed_send, packed_recv)
            endif
        endif
        if (.not. allocated(packed_send)) then
            allocate(packed_send(npack), packed_recv(npack))
        endif

        pos = 1

        nitem = size(turbineArray(i) % bladeForces)
        call atm_pack_rank4(turbineArray(i) % bladeForces, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors)
        call atm_pack_rank5(turbineArray(i) % bladeAlignedVectors, packed_send,  &
            pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % alpha)
        call atm_pack_rank3(turbineArray(i) % alpha, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % lift)
        call atm_pack_rank3(turbineArray(i) % lift, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % drag)
        call atm_pack_rank3(turbineArray(i) % drag, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cl)
        call atm_pack_rank3(turbineArray(i) % Cl, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cd)
        call atm_pack_rank3(turbineArray(i) % Cd, packed_send, pos)
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            call atm_pack_rank3(turbineArray(i) % Cm, packed_send, pos)
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            call atm_pack_rank3(turbineArray(i) % pitchingMoment, packed_send, &
                pos)
            pos = pos + nitem
        endif

        nitem = size(turbineArray(i) % Vmag)
        call atm_pack_rank3(turbineArray(i) % Vmag, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % axialForce)
        call atm_pack_rank3(turbineArray(i) % axialForce, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % tangentialForce)
        call atm_pack_rank3(turbineArray(i) % tangentialForce, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % windVectors)
        call atm_pack_rank4(turbineArray(i) % windVectors, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % induction_a)
        call atm_pack_rank3(turbineArray(i) % induction_a, packed_send, pos)
        pos = pos + nitem

        nitem = size(turbineArray(i) % u_infinity)
        call atm_pack_rank3(turbineArray(i) % u_infinity, packed_send, pos)
        pos = pos + nitem

        call atm_pack_gather_scalars(packed_send, pos,                         &
            turbineArray(i) % torqueRotor, turbineArray(i) % thrust,            &
            turbineArray(i) % nacelleForce(1), turbineArray(i) % nacelleForce(2),&
            turbineArray(i) % nacelleForce(3),                                  &
            turbineArray(i) % VelNacelle_sampled,                               &
            turbineArray(i) % VelNacelle_corrected)
        call atm_cuda_sync('ATM packed gather pack')

        call mpi_allreduce(packed_send, packed_recv, npack, mpi_rprec,          &
                           mpi_sum, TURBINE_COMMUNICATOR, ierr)

        pos = 1

        nitem = size(turbineArray(i) % bladeForces)
        call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % bladeForces)
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors)
        call atm_unpack_rank5(packed_recv, pos,                                 &
            turbineArray(i) % bladeAlignedVectors)
        pos = pos + nitem

        nitem = size(turbineArray(i) % alpha)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % alpha)
        pos = pos + nitem

        nitem = size(turbineArray(i) % lift)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % lift)
        pos = pos + nitem

        nitem = size(turbineArray(i) % drag)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % drag)
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cl)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cl)
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cd)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cd)
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Cm)
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            call atm_unpack_rank3(packed_recv, pos,                            &
                turbineArray(i) % pitchingMoment)
            pos = pos + nitem
        endif

        nitem = size(turbineArray(i) % Vmag)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % Vmag)
        pos = pos + nitem

        nitem = size(turbineArray(i) % axialForce)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % axialForce)
        pos = pos + nitem

        nitem = size(turbineArray(i) % tangentialForce)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % tangentialForce)
        pos = pos + nitem

        nitem = size(turbineArray(i) % windVectors)
        call atm_unpack_rank4(packed_recv, pos, turbineArray(i) % windVectors)
        pos = pos + nitem

        nitem = size(turbineArray(i) % induction_a)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % induction_a)
        pos = pos + nitem

        nitem = size(turbineArray(i) % u_infinity)
        call atm_unpack_rank3(packed_recv, pos, turbineArray(i) % u_infinity)
        pos = pos + nitem

        turbineArray(i) % torqueRotor = packed_recv(pos)
        turbineArray(i) % thrust = packed_recv(pos+1)
        turbineArray(i) % nacelleForce = packed_recv(pos+2:pos+4)
        turbineArray(i) % VelNacelle_sampled = packed_recv(pos+5)
        turbineArray(i) % VelNacelle_corrected = packed_recv(pos+6)

        call atm_cuda_check('ATM packed gather unpack')

    endif

enddo

end subroutine atm_lesgo_mpi_gather_packed_gpu

#ifdef PPLES_GPU
subroutine atm_pack_blade_forces_mirror(iturb, buf, offset)
implicit none

integer, intent(in) :: iturb
real(rprec), device, intent(inout) :: buf(:)
integer, intent(in) :: offset
integer :: idx, m, n, q, a, mmend, nnend, qqend, rem, nitem, j

if (.not. atm_blade_mirror_ready) return

j = turbineArray(iturb) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(iturb) % numAnnulusSections
qqend = turbineArray(iturb) % numBladePoints
nitem = mmend * nnend * qqend * 3

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, nitem
    rem = idx - 1
    m = mod(rem, mmend) + 1
    rem = rem / mmend
    n = mod(rem, nnend) + 1
    rem = rem / nnend
    q = mod(rem, qqend) + 1
    a = rem / qqend + 1
    buf(offset + idx - 1) = atm_bladeForces_d(m,n,q,a,iturb)
end do

end subroutine atm_pack_blade_forces_mirror

subroutine atm_unpack_blade_forces_mirror(buf, offset, iturb)
implicit none

real(rprec), device, intent(in) :: buf(:)
integer, intent(in) :: offset, iturb
integer :: idx, m, n, q, a, mmend, nnend, qqend, rem, nitem, j

if (.not. atm_blade_mirror_ready) return

j = turbineArray(iturb) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(iturb) % numAnnulusSections
qqend = turbineArray(iturb) % numBladePoints
nitem = mmend * nnend * qqend * 3

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, nitem
    rem = idx - 1
    m = mod(rem, mmend) + 1
    rem = rem / mmend
    n = mod(rem, nnend) + 1
    rem = rem / nnend
    q = mod(rem, qqend) + 1
    a = rem / qqend + 1
    atm_bladeForces_d(m,n,q,a,iturb) = buf(offset + idx - 1)
end do

end subroutine atm_unpack_blade_forces_mirror
#endif

subroutine atm_pack_rank3(src, buf, offset)
implicit none

real(rprec), managed, intent(in) :: src(:,:,:)
real(rprec), device, intent(inout) :: buf(:)
integer, intent(in) :: offset
integer :: idx, i1, i2, i3, n1, n2, n3

n1 = size(src, 1)
n2 = size(src, 2)
n3 = size(src, 3)

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, n1*n2*n3
    i1 = mod(idx - 1, n1) + 1
    i2 = mod((idx - 1) / n1, n2) + 1
    i3 = (idx - 1) / (n1*n2) + 1
    buf(offset + idx - 1) = src(i1,i2,i3)
end do

end subroutine atm_pack_rank3

subroutine atm_unpack_rank3(buf, offset, dst)
implicit none

real(rprec), device, intent(in) :: buf(:)
real(rprec), managed, intent(inout) :: dst(:,:,:)
integer, intent(in) :: offset
integer :: idx, i1, i2, i3, n1, n2, n3

n1 = size(dst, 1)
n2 = size(dst, 2)
n3 = size(dst, 3)

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, n1*n2*n3
    i1 = mod(idx - 1, n1) + 1
    i2 = mod((idx - 1) / n1, n2) + 1
    i3 = (idx - 1) / (n1*n2) + 1
    dst(i1,i2,i3) = buf(offset + idx - 1)
end do

end subroutine atm_unpack_rank3

subroutine atm_pack_rank4(src, buf, offset)
implicit none

real(rprec), managed, intent(in) :: src(:,:,:,:)
real(rprec), device, intent(inout) :: buf(:)
integer, intent(in) :: offset
integer :: idx, i1, i2, i3, i4, n1, n2, n3, n4, rem

n1 = size(src, 1)
n2 = size(src, 2)
n3 = size(src, 3)
n4 = size(src, 4)

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, n1*n2*n3*n4
    rem = idx - 1
    i1 = mod(rem, n1) + 1
    rem = rem / n1
    i2 = mod(rem, n2) + 1
    rem = rem / n2
    i3 = mod(rem, n3) + 1
    i4 = rem / n3 + 1
    buf(offset + idx - 1) = src(i1,i2,i3,i4)
end do

end subroutine atm_pack_rank4

subroutine atm_unpack_rank4(buf, offset, dst)
implicit none

real(rprec), device, intent(in) :: buf(:)
real(rprec), managed, intent(inout) :: dst(:,:,:,:)
integer, intent(in) :: offset
integer :: idx, i1, i2, i3, i4, n1, n2, n3, n4, rem

n1 = size(dst, 1)
n2 = size(dst, 2)
n3 = size(dst, 3)
n4 = size(dst, 4)

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, n1*n2*n3*n4
    rem = idx - 1
    i1 = mod(rem, n1) + 1
    rem = rem / n1
    i2 = mod(rem, n2) + 1
    rem = rem / n2
    i3 = mod(rem, n3) + 1
    i4 = rem / n3 + 1
    dst(i1,i2,i3,i4) = buf(offset + idx - 1)
end do

end subroutine atm_unpack_rank4

subroutine atm_pack_rank5(src, buf, offset)
implicit none

real(rprec), managed, intent(in) :: src(:,:,:,:,:)
real(rprec), device, intent(inout) :: buf(:)
integer, intent(in) :: offset
integer :: idx, i1, i2, i3, i4, i5, n1, n2, n3, n4, n5, rem

n1 = size(src, 1)
n2 = size(src, 2)
n3 = size(src, 3)
n4 = size(src, 4)
n5 = size(src, 5)

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, n1*n2*n3*n4*n5
    rem = idx - 1
    i1 = mod(rem, n1) + 1
    rem = rem / n1
    i2 = mod(rem, n2) + 1
    rem = rem / n2
    i3 = mod(rem, n3) + 1
    rem = rem / n3
    i4 = mod(rem, n4) + 1
    i5 = rem / n4 + 1
    buf(offset + idx - 1) = src(i1,i2,i3,i4,i5)
end do

end subroutine atm_pack_rank5

subroutine atm_unpack_rank5(buf, offset, dst)
implicit none

real(rprec), device, intent(in) :: buf(:)
real(rprec), managed, intent(inout) :: dst(:,:,:,:,:)
integer, intent(in) :: offset
integer :: idx, i1, i2, i3, i4, i5, n1, n2, n3, n4, n5, rem

n1 = size(dst, 1)
n2 = size(dst, 2)
n3 = size(dst, 3)
n4 = size(dst, 4)
n5 = size(dst, 5)

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, n1*n2*n3*n4*n5
    rem = idx - 1
    i1 = mod(rem, n1) + 1
    rem = rem / n1
    i2 = mod(rem, n2) + 1
    rem = rem / n2
    i3 = mod(rem, n3) + 1
    rem = rem / n3
    i4 = mod(rem, n4) + 1
    i5 = rem / n4 + 1
    dst(i1,i2,i3,i4,i5) = buf(offset + idx - 1)
end do

end subroutine atm_unpack_rank5

subroutine atm_pack_gather_scalars(buf, offset, torque, thrust, nf1, nf2, nf3,  &
    vel_sampled, vel_corrected)
implicit none

real(rprec), device, intent(inout) :: buf(:)
integer, intent(in) :: offset
real(rprec), intent(in) :: torque, thrust, nf1, nf2, nf3
real(rprec), intent(in) :: vel_sampled, vel_corrected
integer :: idx

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, 7
    select case (idx)
    case (1)
        buf(offset) = torque
    case (2)
        buf(offset+1) = thrust
    case (3)
        buf(offset+2) = nf1
    case (4)
        buf(offset+3) = nf2
    case (5)
        buf(offset+4) = nf3
    case (6)
        buf(offset+5) = vel_sampled
    case (7)
        buf(offset+6) = vel_corrected
    end select
end do

end subroutine atm_pack_gather_scalars
#endif
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will feed the velocity at all the actuator points into the atm
! This is done by using trilinear interpolation from lesgo
! Force will be calculated based on the velocities and stored on forceField
implicit none

integer, intent(in) :: i ! The turbine number
integer :: m,n,q,j
! mpi_velocity only used for Spalart method
real(rprec), dimension(3) :: velocity, mpi_velocity
real(rprec), dimension(3) :: xyz    ! Point onto which to interpolate velocity
real(rprec), pointer, dimension(:) :: x,y,z,zw
#ifdef PPLES_GPU
integer :: mm, nn, qq, p0
real(rprec), allocatable :: velbp(:,:,:,:)  ! device-sampled velocity (3,m,n,q)
logical,     allocatable :: inr(:,:,:)      ! in-domain flag per blade point
#endif

! The MPI turbine communcator
integer, pointer :: TURBINE_COMM

TURBINE_COMM => turbineArray(i) % TURBINE_COMM_WORLD

j=turbineArray(i) % turbineTypeID ! The turbine type ID

#ifdef ENABLE_CUDA
! The legacy CUDA at-point blade-force kernel does not cover the structural
! moment/torsion feedback path, so structure-on still uses the host route.
if (atm_bladeforce_cuda_enabled() .and. .not. atm_structure_enabled() .and.  &
    turbineArray(i) % sampling == 'atPoint') then
    call atm_lesgo_force_gpu_atpoint(i)
    if (turbineArray(i) % nacelle) then
        call atm_lesgo_nacelle_force_gpu_atpoint(i)
    endif
    return
endif
#endif

! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

if (turbineArray(i) % sampling == 'Spalart') then
    ! This loop goes through all the blade points and calculates the respective
    ! body forces then imposes it onto the force field
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl

                ! Actuator point onto which to interpolate the velocity
                xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

                velocity = 0._rprec
                mpi_velocity = 0._rprec

                call atm_lesgo_compute_spalart_u(i, xyz, velocity)

                mpi_velocity = velocity

                ! Complie this subroutines only if MPI will be used
#ifdef PPMPI
!~                     call mpi_barrier( TURBINE_COMM, ierr )
                    ! Sync all the blade forces
                    call mpi_allreduce(mpi_velocity, velocity, size(velocity), &
                           mpi_rprec, mpi_sum, TURBINE_COMM , ierr)
#endif

                ! This will compute the blade force for the specific point
                if (  z(1) <= xyz(3)/z_i .and. xyz(3)/z_i < z(nz) ) then
                    call atm_computeBladeForce(i,m,n,q,velocity)
                else
                    velocity = 0._rprec
                endif

            enddo
        enddo
    enddo


else if (turbineArray(i) % sampling == 'atPoint') then
    ! This loop goes through all the blade points and calculates the respective
    ! body forces then imposes it onto the force field
#ifdef PPLES_GPU
    ! Sample velocity at the blade points on the DEVICE (reads device u,v,w_uv),
    ! so the host no longer needs update self(u,v,w). The airfoil force model
    ! (atm_computeBladeForce) still runs on the host using the sampled velocity.
    mm = turbineModel(j) % numBl
    nn = turbineArray(i) % numAnnulusSections
    qq = turbineArray(i) % numBladePoints
    if (atm_batch_sampled) then
        ! Consume this turbine's slice of the batched sample (same values, same
        ! call order as the per-turbine path; flatten order is m-outer/q-fastest).
        do q = 1, qq
        do n = 1, nn
        do m = 1, mm
            p0 = atm_bp_off(i) + ((m-1)*nn + (n-1))*qq + q
            if (atm_inr_all(p0) == 1)                                          &
                call atm_computeBladeForce(i, m, n, q, atm_velbp_all(1:3,p0))
        end do
        end do
        end do
    else
        allocate(velbp(3,mm,nn,qq), inr(mm,nn,qq))
        call atm_sample_velocity_atpoint_gpu(i, velbp, inr)
        do q = 1, qq
        do n = 1, nn
        do m = 1, mm
            if (inr(m,n,q)) call atm_computeBladeForce(i, m, n, q, velbp(1:3,m,n,q))
        end do
        end do
        end do
        deallocate(velbp, inr)
    end if
#else
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl

                ! Actuator point onto which to interpolate the velocity
                xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

                ! Non-dimensionalizes the point location
                xyz=xyz/z_i

                ! Interpolate velocities if inside the domain
                if (  z(1) <= xyz(3) .and. xyz(3) < z(nz) ) then
                    velocity(1)=                                               &
                    trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                    velocity(2)=                                               &
                    trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                    velocity(3)=                                               &
                    trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

                    ! This will compute the blade force for the specific point
                    call atm_computeBladeForce(i,m,n,q,velocity)

                endif

            enddo
        enddo
    enddo
#endif
endif

    ! Calculate Nacelle force
    if (turbineArray(i) % nacelle) then
        xyz=turbineArray(i) % nacelleLocation
        xyz=xyz/z_i
        if (  z(1) <= xyz(3) .and. xyz(3) < z(nz) ) then

            velocity(1)=                                                   &
            trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
            velocity(2)=                                                   &
            trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
            velocity(3)=                                                   &
            trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

            call atm_computeNacelleForce(i,velocity)

        endif
    endif

end subroutine atm_lesgo_force

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_compute_Spalart_u(i, xyz, velocity)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will calculate the sampling velocity using the proposed method
! from Spalart
! n turbine number
! xyz actuator point position vector
! velocity reference velocity for computing lift and drag

implicit none

integer, intent(in) :: i
real(rprec), intent(in) :: xyz(3)
real(rprec), intent(inout) :: velocity(3)

integer :: c, m, n, q

! Pointers for mesh
real(rprec), pointer, dimension(:) :: z,zw

! Test for time optimization
real(rprec) :: dist, a(3), projectradius, epsilon
real(rprec) :: epsilon_sq, kernel_norm

nullify(z,zw)
z => grid % z
zw => grid % zw

! Value of epsilon
epsilon=turbineArray(i) % epsilon
epsilon_sq = epsilon * epsilon
kernel_norm = 1._rprec / ((epsilon * epsilon_sq) * (pi * sqrt(pi)))

! Projection radius
projectradius = turbineArray(i) % projectionRadius

! Set the velocity to zero
velocity = 0._rprec


do c=1,forceFieldUV(i) % c

    a = forceFieldUV(i) %  location(1:3, c)
    m = forceFieldUV(i) %  ijk(1, c)
    n = forceFieldUV(i) %  ijk(2, c)
    q = forceFieldUV(i) %  ijk(3, c)

    dist = sqrt((a(1)-xyz(1))**2 + (a(2)-xyz(2))**2 + (a(3)-xyz(3))**2)
    if (dist .le. projectradius * z_i) then
        if ( z(1) <= a(3)/z_i .and. a(3)/z_i < z(nz)) then

        ! The value of the kernel. This is the actual smoothing function
        velocity(1) = velocity(1) + u(m,n,q) * exp(-dist*dist/epsilon_sq)  &
                                 * kernel_norm
        velocity(2) = velocity(2) + v(m,n,q) * exp(-dist*dist/epsilon_sq)  &
                                 * kernel_norm
        endif
    endif
enddo

do c=1,forceFieldW(i) % c
    a = forceFieldW(i) %  location(1:3, c)
    m = forceFieldW(i) %  ijk(1, c)
    n = forceFieldW(i) %  ijk(2, c)
    q = forceFieldW(i) %  ijk(3, c)

    dist = sqrt((a(1)-xyz(1))**2 + (a(2)-xyz(2))**2 + (a(3)-xyz(3))**2)

    if (dist .le. projectradius) then
        if ( z(1) <= a(3)/z_i .and. a(3)/z_i < z(nz)) then

        ! The value of the kernel. This is the actual smoothing function
        velocity(3) = velocity(3) + w(m,n,q) * exp(-dist*dist/epsilon_sq)  &
                                 * kernel_norm
        endif
    endif
enddo

velocity = velocity * u_star * z_i * dx * z_i * dy *z_i * dz

end subroutine atm_lesgo_compute_Spalart_u

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_convolute_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will convolute the forces for each turbine

implicit none

!~ type(clock_t) :: myClock

integer, intent(in) :: i
integer :: j, m, n, q, c,mmend,nnend,qqend

integer :: ii, jj, kk  ! Indices for lesgo fields

! Test for time optimization
real(rprec) :: dist,a(3),b(3),projectradius,epsilon,const1,const2,const3
real(rprec) :: epsilon_sq, nacelle_epsilon_sq, nacelle_kernel_norm
real(rprec) :: nacelleEpsilon

! Variables for convolution force
real(rprec) :: kernel, force(3)

! Pointers for the turbineArray quantities
real(rprec), pointer, dimension(:,:,:,:) :: bladeForces, bladePoints

real(rprec), pointer, dimension(:,:) :: bodyForceUV, bodyForceW

#ifdef ENABLE_CUDA
if (atm_convolve_cuda_enabled() .and.                                      &
    turbineArray(i) % sampling == 'atPoint') then
    call atm_lesgo_convolute_force_gpu_atpoint(i)
    return
endif
#endif
#ifdef PPLES_GPU
! Explicit-residency: the 'atPoint' Gaussian convolution (the ATM bottleneck,
! ~72% of ATM time) runs on the GPU. Spalart still uses the host loops below.
if (turbineArray(i) % sampling == 'atPoint') then
    call atm_convolute_atpoint_gpu(i)
    return
end if
#endif

nullify(bladeForces)
nullify(bladePoints)
nullify(bodyForceUV)
nullify(bodyForceW)

bladeForces => turbineArray(i) % bladeForces
bladePoints => turbineArray(i) % bladePoints

bodyForceUV => forceFieldUV(i) % force
bodyForceW =>  forceFieldW(i) % force

!real(rprec) :: dummyForce(3)  ! Debugging

j=turbineArray(i) % turbineTypeID ! The turbine type ID

! This will convolute the blade force onto the grid points
! affected by the turbines on both grids
! Only if the distance is less than specified value
mmend=turbineModel(j) % numBl
nnend=turbineArray(i) % numAnnulusSections
qqend=turbineArray(i) % numBladePoints
projectradius=turbineArray(i) % projectionRadius
epsilon=turbineArray(i) % epsilon
nacelleEpsilon = turbineArray(i) % nacelleEpsilon
epsilon_sq = epsilon * epsilon
nacelle_epsilon_sq = nacelleEpsilon * nacelleEpsilon
const1 = 1._rprec / ((epsilon * epsilon_sq) * (pi * sqrt(pi)))
const2 = z_i / (u_star*u_star)
const3=const1*const2
nacelle_kernel_norm = 1._rprec /                                          &
    ((nacelleEpsilon * nacelle_epsilon_sq) * (pi * sqrt(pi)))

! Body Force implementation using velocity sampling at the actuator point
if (turbineArray(i) % sampling == 'atPoint') then

    !~  call myClock % start()
    do c=1,forceFieldUV(i) % c
        a= forceFieldUV(i) %  location(1:3,c)
        force=0._rprec

        ! Blade forces
        do m=1, mmend
            do n=1, nnend
               do q=1, qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2            &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                    ! The value of the kernel. This is the actual smoothing function
                     force(1:2) = force(1:2) + bladeForces(m,n,q,1:2)     &
                                  * exp(-dist*dist/epsilon_sq)
                    endif

                enddo
            enddo
        enddo
        force(1:2)=force(1:2)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
    !~         if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                &
                         * nacelle_kernel_norm
                !write(*,*) 'kernel Value= ', kernel
                force(1:2) = force(1:2)+turbineArray(i) % nacelleForce(1:2) *  &
                             kernel *const2
    !~          integrateNacelleForce=integrateNacelleForce+force(1) * dx *dy * dz * z_i**3

    !~         endif
        endif


        bodyForceUV(1:2,c) = force(1:2)
    !~     if (abs(bodyForceUV(1,c)) .gt. 0) then
    !~                 write(*,*) 'bodyForceUV is: ', bodyForceUV(1,c)
    !~     endif
    enddo


    do c=1,forceFieldW(i) % c
        a= forceFieldW(i) %  location(1:3,c)
        force=0._rprec

        ! Blade forces
        do m=1,mmend
            do n=1,nnend
               do q=1,qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2            &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                    ! The value of the kernel. This is the actual smoothing function
                    force(3) = force(3) +  bladeForces(m,n,q,3) &
                               * exp(-dist*dist/epsilon_sq)
                    endif

                enddo
            enddo
        enddo
        force(3)=force(3)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
            if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                &
                         * nacelle_kernel_norm
                force(3) = force(3)+turbineArray(i) % nacelleForce(3) *           &
                           kernel *const2
            endif
        endif

        bodyForceW(3,c) = force(3)
    enddo

! The Spalart method uses the local velocity field.
! For this reason it needs to be done explicitly in this module
! and cannot be generally coded from the actuator_turbine_model module
elseif (turbineArray(i) % sampling == 'Spalart') then

    !~  call myClock % start()
    do c=1,forceFieldUV(i) % c
        a= forceFieldUV(i) %  location(1:3,c)
        force=0._rprec
        ! Indices for velocity field
        ii = forceFieldUV(i) % ijk(1,c)
        jj = forceFieldUV(i) % ijk(2,c)
        kk = forceFieldUV(i) % ijk(3,c)

        ! Blade forces
        do m=1, mmend
            do n=1, nnend
               do q=1, qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2            &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                        ! The value of the kernel.
                        ! This is the actual smoothing function
                        ! Divide by velocity magnitude
                         force(1) = force(1) +  bladeForces(m,n,q,1) *         &
                                      exp(-dist*dist/epsilon_sq)              &
                         / (turbineArray(i) % Vmag(m,n,q)) *                     &
                         ( u(ii,jj,kk)  * u_star +                             &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))

                         force(2) = force(2) +  bladeForces(m,n,q,2) *         &
                                      exp(-dist*dist/epsilon_sq)              &
                         / turbineArray(i) % Vmag(m,n,q) *                     &
                         ( v(ii,jj,kk) * u_star +                              &
                         turbineArray(i) % bladeAlignedVectors(m,n,q,2,2) *    &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))

                    endif

                enddo
            enddo
        enddo
        force(1:2)=force(1:2)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
    !~         if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                   &
                         * nacelle_kernel_norm
                !write(*,*) 'kernel Value= ', kernel
                force(1:2) = force(1:2)+turbineArray(i) % nacelleForce(1:2) *  &
                             kernel *const2
    !~          integrateNacelleForce=integrateNacelleForce+force(1) * dx *dy * dz * z_i**3

    !~         endif
        endif


        bodyForceUV(1:2,c) = force(1:2)
    enddo


    do c=1,forceFieldW(i) % c
        a= forceFieldW(i) %  location(1:3,c)
        force=0._rprec
        ! Indices for velocity field
        ii = forceFieldW(i) % ijk(1,c)
        jj = forceFieldW(i) % ijk(2,c)
        kk = forceFieldW(i) % ijk(3,c)

        ! Blade forces
        do m=1,mmend
            do n=1,nnend
               do q=1,qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2               &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                        ! The value of the kernel.
                        ! This is the actual smoothing function
                        force(3) = force(3) +  bladeForces(m,n,q,3) *          &
                                   exp(-dist*dist/epsilon_sq)                 &
                         / turbineArray(i) % Vmag(m,n,q) *                     &
                         ( w(ii,jj,kk) * u_star +                              &
                         turbineArray(i) % bladeAlignedVectors(m,n,q,2,3) *    &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))
                    endif

                enddo
            enddo
        enddo
        force(3)=force(3)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
            if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                  &
                         * nacelle_kernel_norm
                force(3) = force(3)+turbineArray(i) % nacelleForce(3) *           &
                           kernel *const2
            endif
        endif

        bodyForceW(3,c) = force(3)
    enddo

endif

end subroutine atm_lesgo_convolute_force

#ifdef PPLES_GPU
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_convolute_atpoint_gpu(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OpenACC version of the sampling=='atPoint' force convolution. Mirrors the host
! loops in atm_lesgo_convolute_force EXACTLY (same Gaussian kernel, same const1/2/3,
! same nacelle treatment: UV adds the nacelle term unconditionally, W gates it on
! dist<=projectradius). The only numerical difference is the floating-point
! summation ORDER over blade points (associativity), so results match to round-off.
! Per-turbine flat scratch is staged to the device for the kernels. The force is
! ALSO scattered into the device fxa/fya/fza here (the apply), so the host
! atm_lesgo_apply_force is skipped for atPoint and the full-field update
! device(fxa,fya,fza) in forcing.f90 is dropped. forceField%force is still
! written back to the host (a force-integration diagnostic reads it).
use sim_param, only : fxa, fya, fza
implicit none
integer, intent(in) :: i
integer :: j, m, n, q, c, p, nbp, cUV, cW
integer :: ii, jj, kk
real(rprec) :: projectradius, eps, eps_sq, const1, const2, const3
real(rprec) :: nacEps, nac_eps_sq, nac_norm
real(rprec) :: dist, fx1, fx2, fx3, kw
real(rprec) :: nlx, nly, nlz, nf1, nf2, nf3
logical :: has_nac
real(rprec), allocatable :: bp(:,:), bf(:,:), locU(:,:), locW(:,:), fU(:,:), fW(:)
integer,     allocatable :: ijkU(:,:), ijkW(:,:)

j   = turbineArray(i) % turbineTypeID
nbp = turbineModel(j) % numBl * turbineArray(i) % numAnnulusSections                &
                              * turbineArray(i) % numBladePoints
cUV = forceFieldUV(i) % c
cW  = forceFieldW(i)  % c

projectradius = turbineArray(i) % projectionRadius
eps           = turbineArray(i) % epsilon
eps_sq        = eps * eps
const1        = 1._rprec / ((eps * eps_sq) * (pi * sqrt(pi)))
const2        = z_i / (u_star * u_star)
const3        = const1 * const2
nacEps        = turbineArray(i) % nacelleEpsilon
nac_eps_sq    = nacEps * nacEps
nac_norm      = 1._rprec / ((nacEps * nac_eps_sq) * (pi * sqrt(pi)))
has_nac       = turbineArray(i) % nacelle
nlx = 0._rprec; nly = 0._rprec; nlz = 0._rprec
nf1 = 0._rprec; nf2 = 0._rprec; nf3 = 0._rprec
if (has_nac) then
    nlx = turbineArray(i) % nacelleLocation(1)
    nly = turbineArray(i) % nacelleLocation(2)
    nlz = turbineArray(i) % nacelleLocation(3)
    nf1 = turbineArray(i) % nacelleForce(1)
    nf2 = turbineArray(i) % nacelleForce(2)
    nf3 = turbineArray(i) % nacelleForce(3)
end if

! Flatten blade points/forces to (3, nbp) on the host (nbp is small)
allocate(bp(3,nbp), bf(3,nbp))
p = 0
do m = 1, turbineModel(j) % numBl
do n = 1, turbineArray(i) % numAnnulusSections
do q = 1, turbineArray(i) % numBladePoints
    p = p + 1
    bp(1,p) = turbineArray(i) % bladePoints(m,n,q,1)
    bp(2,p) = turbineArray(i) % bladePoints(m,n,q,2)
    bp(3,p) = turbineArray(i) % bladePoints(m,n,q,3)
    bf(1,p) = turbineArray(i) % bladeForces(m,n,q,1)
    bf(2,p) = turbineArray(i) % bladeForces(m,n,q,2)
    bf(3,p) = turbineArray(i) % bladeForces(m,n,q,3)
end do
end do
end do

! ---- UV grid (force components 1,2) ----
if (cUV > 0) then
    allocate(locU(3,cUV), fU(2,cUV), ijkU(3,cUV))
    locU(1:3,1:cUV) = forceFieldUV(i) % location(1:3,1:cUV)
    ijkU(1:3,1:cUV) = forceFieldUV(i) % ijk(1:3,1:cUV)
    !$acc data copyin(bp, bf, locU, ijkU) copyout(fU) present(fxa, fya)
    !$acc parallel loop gang vector private(fx1, fx2, dist, kw, ii, jj, kk)
    do c = 1, cUV
        fx1 = 0._rprec
        fx2 = 0._rprec
        !$acc loop seq
        do p = 1, nbp
            dist = sqrt((locU(1,c)-bp(1,p))**2 + (locU(2,c)-bp(2,p))**2          &
                      + (locU(3,c)-bp(3,p))**2)
            if (dist <= projectradius) then
                kw  = exp(-dist*dist/eps_sq)
                fx1 = fx1 + bf(1,p) * kw
                fx2 = fx2 + bf(2,p) * kw
            end if
        end do
        fx1 = fx1 * const3
        fx2 = fx2 * const3
        if (has_nac) then
            dist = sqrt((locU(1,c)-nlx)**2 + (locU(2,c)-nly)**2 + (locU(3,c)-nlz)**2)
            kw   = exp(-dist*dist/nac_eps_sq) * nac_norm
            fx1  = fx1 + nf1 * kw * const2
            fx2  = fx2 + nf2 * kw * const2
        end if
        fU(1,c) = fx1
        fU(2,c) = fx2
        ! Apply on device: scatter into fxa/fya. forceFieldUV cells are distinct
        ! within a turbine, so no atomics needed; across turbines (sequential
        ! calls) the adds accumulate, matching the host apply loop.
        ii = ijkU(1,c); jj = ijkU(2,c); kk = ijkU(3,c)
        fxa(ii,jj,kk) = fxa(ii,jj,kk) + fx1
        fya(ii,jj,kk) = fya(ii,jj,kk) + fx2
    end do
    !$acc end data
    forceFieldUV(i) % force(1:2,1:cUV) = fU(1:2,1:cUV)
    deallocate(locU, fU, ijkU)
end if

! ---- W grid (force component 3) ----
if (cW > 0) then
    allocate(locW(3,cW), fW(cW), ijkW(3,cW))
    locW(1:3,1:cW) = forceFieldW(i) % location(1:3,1:cW)
    ijkW(1:3,1:cW) = forceFieldW(i) % ijk(1:3,1:cW)
    !$acc data copyin(bp, bf, locW, ijkW) copyout(fW) present(fza)
    !$acc parallel loop gang vector private(fx3, dist, kw, ii, jj, kk)
    do c = 1, cW
        fx3 = 0._rprec
        !$acc loop seq
        do p = 1, nbp
            dist = sqrt((locW(1,c)-bp(1,p))**2 + (locW(2,c)-bp(2,p))**2          &
                      + (locW(3,c)-bp(3,p))**2)
            if (dist <= projectradius) then
                fx3 = fx3 + bf(3,p) * exp(-dist*dist/eps_sq)
            end if
        end do
        fx3 = fx3 * const3
        if (has_nac) then
            dist = sqrt((locW(1,c)-nlx)**2 + (locW(2,c)-nly)**2 + (locW(3,c)-nlz)**2)
            if (dist <= projectradius) then
                kw  = exp(-dist*dist/nac_eps_sq) * nac_norm
                fx3 = fx3 + nf3 * kw * const2
            end if
        end if
        fW(c) = fx3
        ! Apply on device: scatter into fza.
        ii = ijkW(1,c); jj = ijkW(2,c); kk = ijkW(3,c)
        fza(ii,jj,kk) = fza(ii,jj,kk) + fx3
    end do
    !$acc end data
    forceFieldW(i) % force(3,1:cW) = fW(1:cW)
    deallocate(locW, fW, ijkW)
end if

deallocate(bp, bf)

end subroutine atm_convolute_atpoint_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_sample_velocity_atpoint_gpu(i, velbp, inr)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Device velocity sampling at the actuator points for sampling=='atPoint'.
! Replicates the host atm_lesgo_force path EXACTLY: xyz = bladePoints/z_i, the
! in-domain test z(1)<=xyz(3)<z(nz), cell_indx (modulo+floor with the same
! boundary thresholds + autowrap), and the 8-corner trilinear_interp of u, v,
! w_uv scaled by u_star. Reads the device-resident u,v,w_uv so the host no
! longer needs update self(u,v,w). Output velbp(1:3,m,n,q) + inr(m,n,q) go to the
! host; the airfoil model atm_computeBladeForce then runs on the host unchanged.
use sim_param, only : u, v
use param,     only : nx, ny, nz, lbz, dx, dy, dz, L_x, L_y, L_z, u_star, z_i
use grid_m,    only : grid
implicit none
integer, intent(in)  :: i
real(rprec), intent(out) :: velbp(:,:,:,:)   ! (3, numBl, numAnnulusSections, numBladePoints)
logical,     intent(out) :: inr(:,:,:)       !    (numBl, numAnnulusSections, numBladePoints)
integer :: m, n, q, mm, nn, qq
integer :: is, js, ks, is1, js1, ks1
real(rprec) :: px, py, pz, xd, yd, zd
real(rprec) :: a1, a2, a3, a4, a5, a6
real(rprec), parameter :: thr = 1.e-9_rprec
real(rprec), allocatable :: gx(:), gy(:), gz(:), bpl(:,:,:,:)
integer,     allocatable :: awi(:), awj(:)

mm = turbineModel(turbineArray(i) % turbineTypeID) % numBl
nn = turbineArray(i) % numAnnulusSections
qq = turbineArray(i) % numBladePoints

! Local copies of the (constant) grid + this turbine's blade points for the device
allocate(gx(nx), gy(ny), gz(lbz:nz), awi(0:nx+1), awj(0:ny+1))
allocate(bpl(mm,nn,qq,3))
gx(1:nx)      = grid % x(1:nx)
gy(1:ny)      = grid % y(1:ny)
gz(lbz:nz)    = grid % z(lbz:nz)
awi(0:nx+1)   = grid % autowrap_i(0:nx+1)
awj(0:ny+1)   = grid % autowrap_j(0:ny+1)
bpl(1:mm,1:nn,1:qq,1:3) = turbineArray(i) % bladePoints(1:mm,1:nn,1:qq,1:3)

!$acc data copyin(gx, gy, gz, awi, awj, bpl) copyout(velbp, inr)                 &
!$acc      present(u, v, w_uv)
!$acc parallel loop collapse(3) gang vector                                      &
!$acc     private(px,py,pz,is,js,ks,is1,js1,ks1,xd,yd,zd,a1,a2,a3,a4,a5,a6)
do q = 1, qq
do n = 1, nn
do m = 1, mm
    px = bpl(m,n,q,1) / z_i
    py = bpl(m,n,q,2) / z_i
    pz = bpl(m,n,q,3) / z_i
    if (gz(1) <= pz .and. pz < gz(nz)) then
        inr(m,n,q) = .true.
        ! ---- cell_indx 'i' (autowrapped) ----
        px = modulo(px, L_x)
        if (abs(px)/L_x < thr) then
            is = 1
        else if (abs(px-L_x)/L_x < thr) then
            is = nx
        else
            is = floor(px/dx) + 1
        end if
        ! ---- cell_indx 'j' (autowrapped) ----
        py = modulo(py, L_y)
        if (abs(py)/L_y < thr) then
            js = 1
        else if (abs(py-L_y)/L_y < thr) then
            js = ny
        else
            js = floor(py/dy) + 1
        end if
        ! ---- cell_indx 'k' (no wrap) ----
        if (abs(pz - gz(nz))/L_z < thr) then
            ks = nz - 1
        else
            ks = floor((pz - gz(1))/dz) + 1
        end if
        is1 = awi(is+1); js1 = awj(js+1); ks1 = ks + 1
        xd  = px - gx(is); yd = py - gy(js); zd = pz - gz(ks)
        ! 8-corner trilinear of u
        a1 = u(is,js,ks)   + xd*(u(is1,js,ks)   - u(is,js,ks))  /dx
        a2 = u(is,js1,ks)  + xd*(u(is1,js1,ks)  - u(is,js1,ks)) /dx
        a3 = u(is,js,ks1)  + xd*(u(is1,js,ks1)  - u(is,js,ks1)) /dx
        a4 = u(is,js1,ks1) + xd*(u(is1,js1,ks1) - u(is,js1,ks1))/dx
        a5 = a1 + yd*(a2-a1)/dy
        a6 = a3 + yd*(a4-a3)/dy
        velbp(1,m,n,q) = (a5 + zd*(a6-a5)/dz) * u_star
        ! v
        a1 = v(is,js,ks)   + xd*(v(is1,js,ks)   - v(is,js,ks))  /dx
        a2 = v(is,js1,ks)  + xd*(v(is1,js1,ks)  - v(is,js1,ks)) /dx
        a3 = v(is,js,ks1)  + xd*(v(is1,js,ks1)  - v(is,js,ks1)) /dx
        a4 = v(is,js1,ks1) + xd*(v(is1,js1,ks1) - v(is,js1,ks1))/dx
        a5 = a1 + yd*(a2-a1)/dy
        a6 = a3 + yd*(a4-a3)/dy
        velbp(2,m,n,q) = (a5 + zd*(a6-a5)/dz) * u_star
        ! w_uv
        a1 = w_uv(is,js,ks)   + xd*(w_uv(is1,js,ks)   - w_uv(is,js,ks))  /dx
        a2 = w_uv(is,js1,ks)  + xd*(w_uv(is1,js1,ks)  - w_uv(is,js1,ks)) /dx
        a3 = w_uv(is,js,ks1)  + xd*(w_uv(is1,js,ks1)  - w_uv(is,js,ks1)) /dx
        a4 = w_uv(is,js1,ks1) + xd*(w_uv(is1,js1,ks1) - w_uv(is,js1,ks1))/dx
        a5 = a1 + yd*(a2-a1)/dy
        a6 = a3 + yd*(a4-a3)/dy
        velbp(3,m,n,q) = (a5 + zd*(a6-a5)/dz) * u_star
    else
        inr(m,n,q)     = .false.
        velbp(1,m,n,q) = 0._rprec
        velbp(2,m,n,q) = 0._rprec
        velbp(3,m,n,q) = 0._rprec
    end if
end do
end do
end do
!$acc end data

deallocate(gx, gy, gz, awi, awj, bpl)

end subroutine atm_sample_velocity_atpoint_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_atpoint_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Build the static device tables for the batched atPoint path: prefix offsets,
! concatenated force-field cell lists (location/ijk/turbine-id - static for
! fixed turbines), grid axes + autowrap, and the per-step staging buffers.
! Called once, lazily, on the first force step.
use grid_m, only : grid
use param,  only : nx, ny, nz, lbz
implicit none
integer :: i, j, c, p, nbp

if (atm_batch_ready) return

allocate(atm_bp_off(numberOfTurbines+1))
atm_bp_off(1) = 0
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    nbp = turbineModel(j) % numBl * turbineArray(i) % numAnnulusSections       &
                                  * turbineArray(i) % numBladePoints
    atm_bp_off(i+1) = atm_bp_off(i) + nbp
end do
atm_nbp_tot = atm_bp_off(numberOfTurbines+1)

atm_cUV_tot = 0
atm_cW_tot  = 0
do i = 1, numberOfTurbines
    atm_cUV_tot = atm_cUV_tot + forceFieldUV(i) % c
    atm_cW_tot  = atm_cW_tot  + forceFieldW(i)  % c
end do

allocate(atm_bp_all(3, max(atm_nbp_tot,1)), atm_bf_all(3, max(atm_nbp_tot,1)))
allocate(atm_velbp_all(3, max(atm_nbp_tot,1)), atm_inr_all(max(atm_nbp_tot,1)))
allocate(atm_tconst(ATM_NTC, numberOfTurbines))
atm_bp_all = 0._rprec; atm_bf_all = 0._rprec
atm_velbp_all = 0._rprec; atm_inr_all = 0; atm_tconst = 0._rprec

allocate(atm_locUV_all(3, max(atm_cUV_tot,1)),                                 &
         atm_ijkUV_all(3, max(atm_cUV_tot,1)), atm_tidUV(max(atm_cUV_tot,1)))
allocate(atm_locW_all(3, max(atm_cW_tot,1)),                                   &
         atm_ijkW_all(3, max(atm_cW_tot,1)), atm_tidW(max(atm_cW_tot,1)))
atm_tidUV = 0; atm_tidW = 0
p = 0
do i = 1, numberOfTurbines
    do c = 1, forceFieldUV(i) % c
        p = p + 1
        atm_locUV_all(1:3,p) = forceFieldUV(i) % location(1:3,c)
        atm_ijkUV_all(1:3,p) = forceFieldUV(i) % ijk(1:3,c)
        atm_tidUV(p) = i
    end do
end do
p = 0
do i = 1, numberOfTurbines
    do c = 1, forceFieldW(i) % c
        p = p + 1
        atm_locW_all(1:3,p) = forceFieldW(i) % location(1:3,c)
        atm_ijkW_all(1:3,p) = forceFieldW(i) % ijk(1:3,c)
        atm_tidW(p) = i
    end do
end do

allocate(atm_gx(nx), atm_gy(ny), atm_gz(lbz:nz))
allocate(atm_awi(0:nx+1), atm_awj(0:ny+1))
atm_gx(1:nx)    = grid % x(1:nx)
atm_gy(1:ny)    = grid % y(1:ny)
atm_gz(lbz:nz)  = grid % z(lbz:nz)
atm_awi(0:nx+1) = grid % autowrap_i(0:nx+1)
atm_awj(0:ny+1) = grid % autowrap_j(0:ny+1)

!$acc enter data copyin(atm_bp_off, atm_bp_all, atm_bf_all,                    &
!$acc                   atm_velbp_all, atm_inr_all, atm_tconst,                &
!$acc                   atm_locUV_all, atm_ijkUV_all, atm_tidUV,               &
!$acc                   atm_locW_all, atm_ijkW_all, atm_tidW,                  &
!$acc                   atm_gx, atm_gy, atm_gz, atm_awi, atm_awj)

atm_batch_ready = .true.
end subroutine atm_batch_atpoint_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_sample_velocity_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched device velocity sampling at the actuator points of ALL turbines in
! one kernel. Identical per-point arithmetic to atm_sample_velocity_atpoint_gpu
! (cell_indx thresholds, autowrap, 8-corner trilinear of u, v, w_uv x u_star);
! only the per-turbine data regions/kernels/syncs are collapsed into one
! update device + one kernel + one update self per force step.
use sim_param, only : u, v, w
use param,     only : nx, ny, nz, lbz, dx, dy, dz, L_x, L_y, L_z, u_star, z_i
implicit none
integer :: i, j, m, n, q, p
integer :: is, js, ks, is1, js1, ks1
real(rprec) :: px, py, pz, xd, yd, zd
real(rprec) :: a1, a2, a3, a4, a5, a6, a7, a8
real(rprec), parameter :: thr = 1.e-9_rprec

call atm_batch_atpoint_init()
atm_batch_sampled = .false.
if (atm_nbp_tot == 0) return

! Flatten the (rotated) blade points: m-outer / q-fastest order -- the same
! order the convolution kernel sums over.
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, turbineArray(i) % numAnnulusSections
    do q = 1, turbineArray(i) % numBladePoints
        p = p + 1
        atm_bp_all(1,p) = turbineArray(i) % bladePoints(m,n,q,1)
        atm_bp_all(2,p) = turbineArray(i) % bladePoints(m,n,q,2)
        atm_bp_all(3,p) = turbineArray(i) % bladePoints(m,n,q,3)
    end do
    end do
    end do
end do
!$acc update device(atm_bp_all)

if (atm_direct_w_enabled()) then
    !$acc parallel loop gang vector default(present)                           &
    !$acc     private(px,py,pz,is,js,ks,is1,js1,ks1,xd,yd,zd,                 &
    !$acc             a1,a2,a3,a4,a5,a6,a7,a8)
    do p = 1, atm_nbp_tot
        px = atm_bp_all(1,p) / z_i
        py = atm_bp_all(2,p) / z_i
        pz = atm_bp_all(3,p) / z_i
        if (atm_gz(1) <= pz .and. pz < atm_gz(nz)) then
            atm_inr_all(p) = 1
            px = modulo(px, L_x)
            if (abs(px)/L_x < thr) then
                is = 1
            else if (abs(px-L_x)/L_x < thr) then
                is = nx
            else
                is = floor(px/dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py)/L_y < thr) then
                js = 1
            else if (abs(py-L_y)/L_y < thr) then
                js = ny
            else
                js = floor(py/dy) + 1
            end if
            if (abs(pz - atm_gz(nz))/L_z < thr) then
                ks = nz - 1
            else
                ks = floor((pz - atm_gz(1))/dz) + 1
            end if
            is1 = atm_awi(is+1); js1 = atm_awj(js+1); ks1 = ks + 1
            xd  = px - atm_gx(is); yd = py - atm_gy(js); zd = pz - atm_gz(ks)
            a1 = u(is,js,ks)   + xd*(u(is1,js,ks)   - u(is,js,ks))  /dx
            a2 = u(is,js1,ks)  + xd*(u(is1,js1,ks)  - u(is,js1,ks)) /dx
            a3 = u(is,js,ks1)  + xd*(u(is1,js,ks1)  - u(is,js,ks1)) /dx
            a4 = u(is,js1,ks1) + xd*(u(is1,js1,ks1) - u(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(1,p) = (a5 + zd*(a6-a5)/dz) * u_star
            a1 = v(is,js,ks)   + xd*(v(is1,js,ks)   - v(is,js,ks))  /dx
            a2 = v(is,js1,ks)  + xd*(v(is1,js1,ks)  - v(is,js1,ks)) /dx
            a3 = v(is,js,ks1)  + xd*(v(is1,js,ks1)  - v(is,js,ks1)) /dx
            a4 = v(is,js1,ks1) + xd*(v(is1,js1,ks1) - v(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(2,p) = (a5 + zd*(a6-a5)/dz) * u_star
            a1 = 0.5_rprec * (w(is,js,ks) + w(is,js,ks+1))
            a2 = 0.5_rprec * (w(is1,js,ks) + w(is1,js,ks+1))
            a3 = 0.5_rprec * (w(is,js1,ks) + w(is,js1,ks+1))
            a4 = 0.5_rprec * (w(is1,js1,ks) + w(is1,js1,ks+1))
            a5 = a1 + xd*(a2-a1)/dx
            a6 = a3 + xd*(a4-a3)/dx
            a7 = a5 + yd*(a6-a5)/dy
            if (ks1 == nz .and. coord < nproc - 1) then
                a1 = atm_wuv_recv_up(is,js)
                a2 = atm_wuv_recv_up(is1,js)
                a3 = atm_wuv_recv_up(is,js1)
                a4 = atm_wuv_recv_up(is1,js1)
            else if (ks1 == nz) then
                a1 = 0.5_rprec * (w(is,js,nz-1) + w(is,js,nz))
                a2 = 0.5_rprec * (w(is1,js,nz-1) + w(is1,js,nz))
                a3 = 0.5_rprec * (w(is,js1,nz-1) + w(is,js1,nz))
                a4 = 0.5_rprec * (w(is1,js1,nz-1) + w(is1,js1,nz))
            else
                a1 = 0.5_rprec * (w(is,js,ks1) + w(is,js,ks1+1))
                a2 = 0.5_rprec * (w(is1,js,ks1) + w(is1,js,ks1+1))
                a3 = 0.5_rprec * (w(is,js1,ks1) + w(is,js1,ks1+1))
                a4 = 0.5_rprec * (w(is1,js1,ks1) + w(is1,js1,ks1+1))
            end if
            a5 = a1 + xd*(a2-a1)/dx
            a6 = a3 + xd*(a4-a3)/dx
            a8 = a5 + yd*(a6-a5)/dy
            atm_velbp_all(3,p) = (a7 + zd*(a8-a7)/dz) * u_star
        else
            atm_inr_all(p)     = 0
            atm_velbp_all(1,p) = 0._rprec
            atm_velbp_all(2,p) = 0._rprec
            atm_velbp_all(3,p) = 0._rprec
        end if
    end do
else
    !$acc parallel loop gang vector default(present)                           &
    !$acc     private(px,py,pz,is,js,ks,is1,js1,ks1,xd,yd,zd,                 &
    !$acc             a1,a2,a3,a4,a5,a6,a7,a8)
    do p = 1, atm_nbp_tot
        px = atm_bp_all(1,p) / z_i
        py = atm_bp_all(2,p) / z_i
        pz = atm_bp_all(3,p) / z_i
        if (atm_gz(1) <= pz .and. pz < atm_gz(nz)) then
            atm_inr_all(p) = 1
            px = modulo(px, L_x)
            if (abs(px)/L_x < thr) then
                is = 1
            else if (abs(px-L_x)/L_x < thr) then
                is = nx
            else
                is = floor(px/dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py)/L_y < thr) then
                js = 1
            else if (abs(py-L_y)/L_y < thr) then
                js = ny
            else
                js = floor(py/dy) + 1
            end if
            if (abs(pz - atm_gz(nz))/L_z < thr) then
                ks = nz - 1
            else
                ks = floor((pz - atm_gz(1))/dz) + 1
            end if
            is1 = atm_awi(is+1); js1 = atm_awj(js+1); ks1 = ks + 1
            xd  = px - atm_gx(is); yd = py - atm_gy(js); zd = pz - atm_gz(ks)
            a1 = u(is,js,ks)   + xd*(u(is1,js,ks)   - u(is,js,ks))  /dx
            a2 = u(is,js1,ks)  + xd*(u(is1,js1,ks)  - u(is,js1,ks)) /dx
            a3 = u(is,js,ks1)  + xd*(u(is1,js,ks1)  - u(is,js,ks1)) /dx
            a4 = u(is,js1,ks1) + xd*(u(is1,js1,ks1) - u(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(1,p) = (a5 + zd*(a6-a5)/dz) * u_star
            a1 = v(is,js,ks)   + xd*(v(is1,js,ks)   - v(is,js,ks))  /dx
            a2 = v(is,js1,ks)  + xd*(v(is1,js1,ks)  - v(is,js1,ks)) /dx
            a3 = v(is,js,ks1)  + xd*(v(is1,js,ks1)  - v(is,js,ks1)) /dx
            a4 = v(is,js1,ks1) + xd*(v(is1,js1,ks1) - v(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(2,p) = (a5 + zd*(a6-a5)/dz) * u_star
            a1 = w_uv(is,js,ks)   + xd*(w_uv(is1,js,ks)   - w_uv(is,js,ks))  /dx
            a2 = w_uv(is,js1,ks)  + xd*(w_uv(is1,js1,ks)  - w_uv(is,js1,ks)) /dx
            a3 = w_uv(is,js,ks1)  + xd*(w_uv(is1,js,ks1)  - w_uv(is,js,ks1)) /dx
            a4 = w_uv(is,js1,ks1) + xd*(w_uv(is1,js1,ks1) - w_uv(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(3,p) = (a5 + zd*(a6-a5)/dz) * u_star
        else
            atm_inr_all(p)     = 0
            atm_velbp_all(1,p) = 0._rprec
            atm_velbp_all(2,p) = 0._rprec
            atm_velbp_all(3,p) = 0._rprec
        end if
    end do
end if

!$acc update self(atm_velbp_all, atm_inr_all)
atm_batch_sampled = .true.

end subroutine atm_batch_sample_velocity_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_convolute_force_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched atPoint Gaussian force convolution + apply for ALL turbines: one
! kernel per grid (UV, W) over the concatenated force-field cell lists,
! scattering straight into the device fxa/fya/fza. The blade-point summation
! order is unchanged, but squared-distance arithmetic avoids a redundant sqrt
! before the Gaussian radius test and exponent. The scatter uses atomics
! because cells of DIFFERENT turbines may coincide where
! projection regions overlap (within a turbine each cell is unique).
! The forceField%force host write-back of the per-turbine path is dropped:
! its only readers in this build are the (disabled) load-balancer validators.
use sim_param, only : fxa, fya, fza
use param,     only : z_i, u_star
implicit none
integer :: i, j, c, p, m, n, q
integer :: ii, jj, kk
real(rprec) :: eps, dist_sq, kw, fx1, fx2, fx3
real(rprec) :: pr, prsq, epsq, cc2, cc3

call atm_batch_atpoint_init()

! Per-turbine constants + flattened blade forces (post-gather values).
do i = 1, numberOfTurbines
    if (turbineArray(i) % operate .and.                                        &
        turbineArray(i) % sampling == 'atPoint') then
        atm_tconst(14,i) = 1._rprec
    else
        atm_tconst(14,i) = 0._rprec
    end if
    ! Keep the same const1*const2 operation order as the per-turbine path.
    ! The convolution kernels below use squared distance directly, so final
    ! fields are roundoff-equivalent rather than byte-identical.
    eps = turbineArray(i) % epsilon
    atm_tconst(1,i) = turbineArray(i) % projectionRadius
    atm_tconst(2,i) = eps * eps
    atm_tconst(4,i) = z_i / (u_star * u_star)
    atm_tconst(3,i) = (1._rprec / ((eps * atm_tconst(2,i)) * (pi * sqrt(pi))))  &
                      * atm_tconst(4,i)
    if (turbineArray(i) % nacelle) then
        atm_tconst(13,i) = 1._rprec
        eps = turbineArray(i) % nacelleEpsilon
        atm_tconst(5,i)    = eps * eps
        atm_tconst(6,i)    = 1._rprec / ((eps * eps * eps) * (pi * sqrt(pi)))
        atm_tconst(7:9,i)  = turbineArray(i) % nacelleLocation(1:3)
        atm_tconst(10:12,i) = turbineArray(i) % nacelleForce(1:3)
    else
        atm_tconst(13,i)   = 0._rprec
        atm_tconst(5,i)    = 1._rprec
        atm_tconst(6,i)    = 0._rprec
        atm_tconst(7:12,i) = 0._rprec
    end if

    j = turbineArray(i) % turbineTypeID
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, turbineArray(i) % numAnnulusSections
    do q = 1, turbineArray(i) % numBladePoints
        p = p + 1
        atm_bf_all(1,p) = turbineArray(i) % bladeForces(m,n,q,1)
        atm_bf_all(2,p) = turbineArray(i) % bladeForces(m,n,q,2)
        atm_bf_all(3,p) = turbineArray(i) % bladeForces(m,n,q,3)
    end do
    end do
    end do
end do
!$acc update device(atm_bf_all, atm_tconst)

! ---- UV grid (force components 1,2) ----
if (atm_cUV_tot > 0) then
    !$acc parallel loop gang vector default(present)                           &
    !$acc     private(i, pr, prsq, epsq, cc2, cc3, fx1, fx2, dist_sq, kw, ii, jj, kk, p)
    do c = 1, atm_cUV_tot
        i = atm_tidUV(c)
        if (atm_tconst(14,i) > 0.5_rprec) then
            pr   = atm_tconst(1,i)
            epsq = atm_tconst(2,i)
            prsq = pr * pr
            cc3  = atm_tconst(3,i)
            cc2  = atm_tconst(4,i)
            fx1 = 0._rprec
            fx2 = 0._rprec
            !$acc loop seq
            do p = atm_bp_off(i)+1, atm_bp_off(i+1)
                dist_sq = (atm_locUV_all(1,c)-atm_bp_all(1,p))**2             &
                        + (atm_locUV_all(2,c)-atm_bp_all(2,p))**2             &
                        + (atm_locUV_all(3,c)-atm_bp_all(3,p))**2
                if (dist_sq <= prsq) then
                    kw  = exp(-dist_sq/epsq)
                    fx1 = fx1 + atm_bf_all(1,p) * kw
                    fx2 = fx2 + atm_bf_all(2,p) * kw
                end if
            end do
            fx1 = fx1 * cc3
            fx2 = fx2 * cc3
            if (atm_tconst(13,i) > 0.5_rprec) then
                dist_sq = (atm_locUV_all(1,c)-atm_tconst(7,i))**2             &
                        + (atm_locUV_all(2,c)-atm_tconst(8,i))**2             &
                        + (atm_locUV_all(3,c)-atm_tconst(9,i))**2
                kw  = exp(-dist_sq/atm_tconst(5,i)) * atm_tconst(6,i)
                fx1 = fx1 + atm_tconst(10,i) * kw * cc2
                fx2 = fx2 + atm_tconst(11,i) * kw * cc2
            end if
            ii = atm_ijkUV_all(1,c); jj = atm_ijkUV_all(2,c); kk = atm_ijkUV_all(3,c)
            !$acc atomic update
            fxa(ii,jj,kk) = fxa(ii,jj,kk) + fx1
            !$acc atomic update
            fya(ii,jj,kk) = fya(ii,jj,kk) + fx2
        end if
    end do
end if

! ---- W grid (force component 3) ----
if (atm_cW_tot > 0) then
    !$acc parallel loop gang vector default(present)                           &
    !$acc     private(i, pr, prsq, epsq, cc2, cc3, fx3, dist_sq, kw, ii, jj, kk, p)
    do c = 1, atm_cW_tot
        i = atm_tidW(c)
        if (atm_tconst(14,i) > 0.5_rprec) then
            pr   = atm_tconst(1,i)
            epsq = atm_tconst(2,i)
            prsq = pr * pr
            cc3  = atm_tconst(3,i)
            cc2  = atm_tconst(4,i)
            fx3 = 0._rprec
            !$acc loop seq
            do p = atm_bp_off(i)+1, atm_bp_off(i+1)
                dist_sq = (atm_locW_all(1,c)-atm_bp_all(1,p))**2              &
                        + (atm_locW_all(2,c)-atm_bp_all(2,p))**2              &
                        + (atm_locW_all(3,c)-atm_bp_all(3,p))**2
                if (dist_sq <= prsq) then
                    fx3 = fx3 + atm_bf_all(3,p) * exp(-dist_sq/epsq)
                end if
            end do
            fx3 = fx3 * cc3
            if (atm_tconst(13,i) > 0.5_rprec) then
                dist_sq = (atm_locW_all(1,c)-atm_tconst(7,i))**2              &
                        + (atm_locW_all(2,c)-atm_tconst(8,i))**2              &
                        + (atm_locW_all(3,c)-atm_tconst(9,i))**2
                if (dist_sq <= prsq) then
                    kw  = exp(-dist_sq/atm_tconst(5,i)) * atm_tconst(6,i)
                    fx3 = fx3 + atm_tconst(12,i) * kw * cc2
                end if
            end if
            ii = atm_ijkW_all(1,c); jj = atm_ijkW_all(2,c); kk = atm_ijkW_all(3,c)
            !$acc atomic update
            fza(ii,jj,kk) = fza(ii,jj,kk) + fx3
        end if
    end do
end if

end subroutine atm_batch_convolute_force_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_clc_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Static tables for the batched Cl/tip correction: per-point decode
! (turbine id, q index, (m,n)-block base, block length) and the static
! chord/bladeRadius flattened in the same m-outer/q-fastest order as
! atm_bp_all.
implicit none
integer :: i, j, m, n, q, p, nn2, qq2

if (atm_clc_ready) return
call atm_batch_atpoint_init()

allocate(atm_pt_turb(max(atm_nbp_tot,1)), atm_pt_q(max(atm_nbp_tot,1)))
allocate(atm_pt_base(max(atm_nbp_tot,1)), atm_pt_qq(max(atm_nbp_tot,1)))
allocate(atm_chord_all(max(atm_nbp_tot,1)), atm_brad_all(max(atm_nbp_tot,1)))
allocate(atm_db_all(max(atm_nbp_tot,1)))
allocate(atm_clc_tc(3, numberOfTurbines))
allocate(atm_wv_all(3, max(atm_nbp_tot,1)))
allocate(atm_cl_all(max(atm_nbp_tot,1)), atm_cd_all(max(atm_nbp_tot,1)))
allocate(atm_vmag_all(max(atm_nbp_tot,1)))
allocate(atm_du_all(3, max(atm_nbp_tot,1)), atm_uyopt_vec_all(3, max(atm_nbp_tot,1)))
allocate(atm_uinf_all(3, max(atm_nbp_tot,1)), atm_uxles_all(3, max(atm_nbp_tot,1)))
allocate(atm_g_all(max(atm_nbp_tot,1)), atm_dg_all(max(atm_nbp_tot,1)))
allocate(atm_epsopt_all(max(atm_nbp_tot,1)))
allocate(atm_uyles_vec_all(3, max(atm_nbp_tot,1)))
allocate(atm_uyles_all(max(atm_nbp_tot,1)), atm_uyopt_all(max(atm_nbp_tot,1)))
atm_pt_turb = 0; atm_pt_q = 0; atm_pt_base = 0; atm_pt_qq = 0
atm_chord_all = 0._rprec; atm_brad_all = 0._rprec; atm_db_all = 0._rprec
atm_clc_tc = 0._rprec
atm_wv_all = 0._rprec; atm_cl_all = 0._rprec; atm_cd_all = 0._rprec
atm_vmag_all = 1._rprec; atm_du_all = 0._rprec; atm_uyopt_vec_all = 0._rprec
atm_uinf_all = 0._rprec; atm_uxles_all = 0._rprec
atm_g_all = 0._rprec; atm_dg_all = 0._rprec; atm_epsopt_all = 0._rprec
atm_uyles_vec_all = 0._rprec; atm_uyles_all = 0._rprec; atm_uyopt_all = 0._rprec

do i = 1, numberOfTurbines
    j   = turbineArray(i) % turbineTypeID
    nn2 = turbineArray(i) % numAnnulusSections
    qq2 = turbineArray(i) % numBladePoints
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, nn2
    do q = 1, qq2
        p = p + 1
        atm_pt_turb(p)  = i
        atm_pt_q(p)     = q
        atm_pt_base(p)  = atm_bp_off(i) + ((m-1)*nn2 + (n-1))*qq2
        atm_pt_qq(p)    = qq2
        atm_chord_all(p) = turbineArray(i) % chord(m,n,q)
        atm_brad_all(p)  = turbineArray(i) % bladeRadius(m,n,q)
        atm_db_all(p)    = turbineArray(i) % db(q)
    end do
    end do
    end do
end do

!$acc enter data copyin(atm_pt_turb, atm_pt_q, atm_pt_base, atm_pt_qq,         &
!$acc                   atm_chord_all, atm_brad_all, atm_db_all, atm_clc_tc,   &
!$acc                   atm_wv_all, atm_cl_all, atm_cd_all, atm_vmag_all,      &
!$acc                   atm_du_all, atm_uyopt_vec_all,                         &
!$acc                   atm_uinf_all, atm_uxles_all,                           &
!$acc                   atm_g_all, atm_dg_all, atm_epsopt_all,                 &
!$acc                   atm_uyles_vec_all, atm_uyles_all, atm_uyopt_all)

atm_clc_ready = .true.
end subroutine atm_batch_clc_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_cl_correction_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched OpenACC port of the current exact-panel induced-velocity Cl/tip correction for
! all operating atPoint turbines with tipALMCorrection.  The LES panel face
! construction, optimal-epsilon panel integral, ux_LES contribution to du, and
! host-visible turbineArray outputs mirror atm_compute_cl_correction().
implicit none
integer :: i, j, m, n, q, p, k, pk, nn2, qq2
logical :: any_active
real(rprec) :: f
real(rprec) :: eps_s, opt_eps_chord, inv_eps_s, inv_2pi_eps_s
real(rprec) :: inv_4pi_sqrt, mag_uinf, inv_mag_tan
real(rprec) :: z_q, z_low, z_high, z_k, half_dzp, eps_opt_k
real(rprec) :: u1, u2, u3, uval, u2val, f_low, f_high, f_up, f_um
real(rprec) :: up, um, g_over_v_k, uy_les_k, uy_opt_k
real(rprec) :: acc_les1, acc_les2, acc_opt1, acc_opt2
real(rprec) :: tan1, tan2
real(rprec), parameter :: u_small = 1.0e-2_rprec

call atm_batch_clc_init()
if (atm_nbp_tot == 0) return

f = 0.1_rprec
inv_4pi_sqrt = 1._rprec / (4._rprec * sqrt(pi))

any_active = .false.
do i = 1, numberOfTurbines
    if (turbineArray(i) % operate .and.                                        &
        turbineArray(i) % sampling == 'atPoint' .and.                          &
        (turbineArray(i) % tipALMCorrection .eqv. .true.)) then
        atm_clc_tc(3,i) = 1._rprec
        any_active = .true.
    else
        atm_clc_tc(3,i) = 0._rprec
        cycle
    end if
    atm_clc_tc(1,i) = turbineArray(i) % epsilon
    atm_clc_tc(2,i) = turbineArray(i) % optimalEpsilonChord

    j   = turbineArray(i) % turbineTypeID
    nn2 = turbineArray(i) % numAnnulusSections
    qq2 = turbineArray(i) % numBladePoints
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, nn2
    do q = 1, qq2
        p = p + 1
        atm_wv_all(1,p) = turbineArray(i) % windVectors(m,n,q,1)
        atm_wv_all(2,p) = turbineArray(i) % windVectors(m,n,q,2)
        atm_wv_all(3,p) = turbineArray(i) % windVectors(m,n,q,3)
        atm_cl_all(p)   = turbineArray(i) % Cl(m,n,q)
        atm_cd_all(p)   = turbineArray(i) % Cd(m,n,q)
        atm_vmag_all(p) = turbineArray(i) % Vmag(m,n,q)
        atm_du_all(1,p) = turbineArray(i) % du(m,n,q,1)
        atm_du_all(2,p) = turbineArray(i) % du(m,n,q,2)
        atm_du_all(3,p) = turbineArray(i) % du(m,n,q,3)
        atm_uyopt_vec_all(1,p) = turbineArray(i) % uy_opt_vec(m,n,q,1)
        atm_uyopt_vec_all(2,p) = turbineArray(i) % uy_opt_vec(m,n,q,2)
        atm_uyopt_vec_all(3,p) = turbineArray(i) % uy_opt_vec(m,n,q,3)
    end do
    end do
    end do
end do
if (.not. any_active) return

!$acc update device(atm_wv_all, atm_cl_all, atm_cd_all, atm_vmag_all,          &
!$acc               atm_du_all, atm_uyopt_vec_all, atm_clc_tc)

! K1: point-local Uinf, ux_LES, G, epsilon_opt.
!$acc parallel loop gang vector default(present)                               &
!$acc     private(i, eps_s, opt_eps_chord, u1, u2, u3, mag_uinf)
do p = 1, atm_nbp_tot
    i = atm_pt_turb(p)
    if (atm_clc_tc(3,i) > 0.5_rprec) then
        eps_s = atm_clc_tc(1,i)
        opt_eps_chord = atm_clc_tc(2,i)
        u1 = atm_wv_all(1,p) - atm_uyopt_vec_all(1,p)
        u2 = atm_wv_all(2,p) - atm_uyopt_vec_all(2,p)
        u3 = atm_wv_all(3,p) - atm_uyopt_vec_all(3,p)
        atm_uinf_all(1,p) = u1
        atm_uinf_all(2,p) = u2
        atm_uinf_all(3,p) = u3
        mag_uinf = sqrt(u1*u1 + u2*u2 + u3*u3)
        atm_uxles_all(1,p) = atm_cd_all(p) * atm_chord_all(p) / eps_s          &
            * inv_4pi_sqrt * atm_vmag_all(p) * u1 / mag_uinf
        atm_uxles_all(2,p) = atm_cd_all(p) * atm_chord_all(p) / eps_s          &
            * inv_4pi_sqrt * atm_vmag_all(p) * u2 / mag_uinf
        atm_uxles_all(3,p) = atm_cd_all(p) * atm_chord_all(p) / eps_s          &
            * inv_4pi_sqrt * atm_vmag_all(p) * u3 / mag_uinf
        atm_g_all(p) = 0.5_rprec * atm_cl_all(p) * atm_chord_all(p)            &
            * atm_vmag_all(p) * atm_vmag_all(p)
        atm_epsopt_all(p) = atm_chord_all(p) * opt_eps_chord
    end if
end do

! K2: dG is retained for backward-compatible diagnostics/output.
!$acc parallel loop gang vector default(present) private(i, q)
do p = 1, atm_nbp_tot
    i = atm_pt_turb(p)
    if (atm_clc_tc(3,i) > 0.5_rprec) then
        q = atm_pt_q(p)
        if (q == 1) then
            atm_dg_all(p) = atm_g_all(p)
        else if (q == atm_pt_qq(p)) then
            atm_dg_all(p) = -atm_g_all(p)
        else
            atm_dg_all(p) = (atm_g_all(p+1) - atm_g_all(p-1)) * 0.5_rprec
        end if
    end if
end do

! K3: exact panel induced velocity, one thread per blade point, k loop in the
! same ascending order as the host routine.
!$acc parallel loop gang vector default(present)                               &
!$acc     private(i, q, k, pk, eps_s, inv_eps_s, inv_2pi_eps_s, z_q, z_low,    &
!$acc             z_high, z_k, half_dzp, eps_opt_k, uval, u2val, f_low,        &
!$acc             f_high, f_up, f_um, up, um, g_over_v_k, uy_les_k,            &
!$acc             uy_opt_k, inv_mag_tan, tan1, tan2, acc_les1, acc_les2,       &
!$acc             acc_opt1, acc_opt2)
do p = 1, atm_nbp_tot
    i = atm_pt_turb(p)
    if (atm_clc_tc(3,i) > 0.5_rprec) then
        q = atm_pt_q(p)
        eps_s = atm_clc_tc(1,i)
        inv_eps_s = 1._rprec / eps_s
        inv_2pi_eps_s = inv_eps_s / (2._rprec * pi)
        z_q = atm_brad_all(p)
        acc_les1 = 0._rprec
        acc_les2 = 0._rprec
        acc_opt1 = 0._rprec
        acc_opt2 = 0._rprec
        !$acc loop seq
        do k = 1, atm_pt_qq(p)
            pk = atm_pt_base(p) + k
            g_over_v_k = atm_g_all(pk) / atm_vmag_all(pk)

            if (k == 1) then
                z_low = atm_brad_all(pk) - 0.5_rprec * atm_db_all(pk)
            else
                z_low = atm_brad_all(pk-1) + 0.5_rprec * atm_db_all(pk-1)
            end if
            z_high = atm_brad_all(pk) + 0.5_rprec * atm_db_all(pk)

            uval = (z_q - z_low) * inv_eps_s
            u2val = uval * uval
            if (abs(uval) < u_small) then
                f_low = 0.5_rprec * uval *                                    &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_low = (1._rprec - exp(-u2val)) / (2._rprec * uval)
            endif
            uval = (z_q - z_high) * inv_eps_s
            u2val = uval * uval
            if (abs(uval) < u_small) then
                f_high = 0.5_rprec * uval *                                   &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_high = (1._rprec - exp(-u2val)) / (2._rprec * uval)
            endif
            uy_les_k = -g_over_v_k * inv_2pi_eps_s * (f_low - f_high)

            z_k = atm_brad_all(pk)
            half_dzp = 0.5_rprec * atm_db_all(pk)
            eps_opt_k = atm_epsopt_all(pk)
            up = (z_q - z_k + half_dzp) / eps_opt_k
            um = (z_q - z_k - half_dzp) / eps_opt_k
            u2val = up * up
            if (abs(up) < u_small) then
                f_up = 0.5_rprec * up *                                      &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_up = (1._rprec - exp(-u2val)) / (2._rprec * up)
            endif
            u2val = um * um
            if (abs(um) < u_small) then
                f_um = 0.5_rprec * um *                                      &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_um = (1._rprec - exp(-u2val)) / (2._rprec * um)
            endif
            uy_opt_k = -g_over_v_k / (2._rprec * pi * eps_opt_k) *            &
                (f_up - f_um)

            inv_mag_tan = 1._rprec / sqrt(atm_uinf_all(1,pk)**2 +             &
                                          atm_uinf_all(2,pk)**2)
            tan1 =  atm_uinf_all(2,pk) * inv_mag_tan
            tan2 = -atm_uinf_all(1,pk) * inv_mag_tan
            acc_les1 = acc_les1 + uy_les_k * tan1
            acc_les2 = acc_les2 + uy_les_k * tan2
            acc_opt1 = acc_opt1 + uy_opt_k * tan1
            acc_opt2 = acc_opt2 + uy_opt_k * tan2
        end do
        atm_uyles_vec_all(1,p) = acc_les1
        atm_uyles_vec_all(2,p) = acc_les2
        atm_uyles_vec_all(3,p) = 0._rprec
        atm_uyopt_vec_all(1,p) = acc_opt1
        atm_uyopt_vec_all(2,p) = acc_opt2
        atm_uyopt_vec_all(3,p) = 0._rprec
        atm_uyles_all(p) = sqrt(acc_les1*acc_les1 + acc_les2*acc_les2)
        atm_uyopt_all(p) = sqrt(acc_opt1*acc_opt1 + acc_opt2*acc_opt2)
        atm_du_all(1,p) = atm_du_all(1,p) * (1._rprec - f) + f *              &
            (atm_uyopt_vec_all(1,p) - atm_uyles_vec_all(1,p) +                &
             atm_uxles_all(1,p))
        atm_du_all(2,p) = atm_du_all(2,p) * (1._rprec - f) + f *              &
            (atm_uyopt_vec_all(2,p) - atm_uyles_vec_all(2,p) +                &
             atm_uxles_all(2,p))
        atm_du_all(3,p) = atm_du_all(3,p) * (1._rprec - f) + f *              &
            (atm_uyopt_vec_all(3,p) - atm_uyles_vec_all(3,p) +                &
             atm_uxles_all(3,p))
    end if
end do

!$acc update self(atm_uinf_all, atm_uxles_all, atm_g_all, atm_dg_all,          &
!$acc             atm_epsopt_all, atm_uyles_vec_all, atm_uyopt_vec_all,        &
!$acc             atm_uyles_all, atm_uyopt_all, atm_du_all)

! Unpack everything the host routine writes, preserving host-visible state for
! downstream force, structural, restart, and output consumers.
do i = 1, numberOfTurbines
    if (atm_clc_tc(3,i) < 0.5_rprec) cycle
    j   = turbineArray(i) % turbineTypeID
    nn2 = turbineArray(i) % numAnnulusSections
    qq2 = turbineArray(i) % numBladePoints
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, nn2
    do q = 1, qq2
        p = p + 1
        turbineArray(i) % Uinf_vec(m,n,q,1:3)   = atm_uinf_all(1:3,p)
        turbineArray(i) % ux_LES_vec(m,n,q,1:3) = atm_uxles_all(1:3,p)
        turbineArray(i) % G(m,n,q)              = atm_g_all(p)
        turbineArray(i) % dG(m,n,q)             = atm_dg_all(p)
        turbineArray(i) % epsilon_opt(m,n,q)    = atm_epsopt_all(p)
        turbineArray(i) % uy_LES_vec(m,n,q,1:3) = atm_uyles_vec_all(1:3,p)
        turbineArray(i) % uy_opt_vec(m,n,q,1:3) = atm_uyopt_vec_all(1:3,p)
        turbineArray(i) % uy_LES(m,n,q)         = atm_uyles_all(p)
        turbineArray(i) % uy_opt(m,n,q)         = atm_uyopt_all(p)
        turbineArray(i) % du(m,n,q,1:3)         = atm_du_all(1:3,p)
    end do
    end do
    end do
end do

end subroutine atm_batch_cl_correction_gpu
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_apply_force()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will apply the blade force onto the CFD grid by using the convolution
! function in the ATM library
implicit none

integer :: c,m
integer :: i,j,k

#ifdef ENABLE_CUDA
if (atm_apply_cuda_enabled()) then
    call atm_lesgo_apply_force_gpu()
    return
endif
#endif

do m=1, numberOfTurbines

#ifdef PPLES_GPU
    ! Explicit-residency: atPoint turbines were already applied on the device
    ! inside atm_convolute_atpoint_gpu (scatter into fxa/fya/fza). Skip them here
    ! to avoid double-application; Spalart still applies on the host below.
    if (turbineArray(m) % sampling == 'atPoint') cycle
#endif
    if (turbineArray(m) % operate) then
        ! Impose force field onto the flow field variables
        ! The forces are non-dimensionalized here as well
        do c=1,forceFieldUV(m) % c
            i=forceFieldUV(m) % ijk(1,c)
            j=forceFieldUV(m) % ijk(2,c)
            k=forceFieldUV(m) % ijk(3,c)

            fxa(i,j,k) = fxa(i,j,k) + forceFieldUV(m) % force(1,c)
            fya(i,j,k) = fya(i,j,k) + forceFieldUV(m) % force(2,c)

        enddo

        do c=1,forceFieldW(m) % c
            i=forceFieldW(m) % ijk(1,c)
            j=forceFieldW(m) % ijk(2,c)
            k=forceFieldW(m) % ijk(3,c)

            fza(i,j,k) = fza(i,j,k) + forceFieldW(m) % force(3,c)

        enddo
    endif
enddo

end subroutine atm_lesgo_apply_force


end module atm_lesgo_interface
