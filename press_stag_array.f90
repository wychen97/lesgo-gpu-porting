#if !defined(PPPRESS_GPU) || (defined(PPLVLSET) && defined(PPLES_GPU) && !defined(ENABLE_CUDA))
!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

! Navigation map:
!   - press_stag_array: pressure RHS, spectral solve, tridiagonal solve,
!     zero mode, and pressure-gradient output
!   - press_apply_env_enabled_unless_false: diagnostic/policy switch parser
!   - press_cuda_sync and event helpers: GPU timing and synchronization
!   - press_*_report routines: pressure-stage diagnostics
!   - press_pack/unpack_* routines: GPU RHS and halo packing
!   - press_assemble_rhs_* and press_zero_mode_cuda: device kernels
!
! The pressure path is performance-sensitive and depends on GPU-aware MPI.
! Keep host-staged fallbacks correct, but do not silently move full pressure
! fields through the host in the regular GPU timestep path.
!
! Ownership map:
!   - sim_param owns p and dpdx/dpdy/dpdz after this routine exits.
!   - this routine owns pressure RHS assembly, zero-mode handling, pressure
!     stage timing, and pressure-RHS halo staging buffers.
!   - tridag_array.f90 and tridag_gpu.f90 own the tridiagonal solve and
!     GPU-aware/host-staged MPI policy after the RHS is assembled.
!   - fft/fft_gpu own transform plans; pressure code should call their wrappers
!     instead of creating per-step FFT state.

!*******************************************************************************
subroutine press_stag_array()
!*******************************************************************************
!
! Calculate the pressure and its derivatives on exit. Everything is in physical
! space on exit.
!
use types, only : rprec
use param
use messages
use sim_param, only : u, v, w, divtz, p, dpdx, dpdy, dpdz
use fft
use emul_complex, only : OPERATOR(.MULI.)
use cuda_mpi_debug, only : mpi_dbg_sendrecv_r, mpi_dbg_send_r, mpi_dbg_recv_r
#ifdef ENABLE_CUDA
use cudafor
use cufft
#endif
#ifdef PPMPI
use mpi
#endif

implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: ir, ii
integer :: jz_min
#ifdef ENABLE_CUDA
real(rprec), managed, save, dimension(:,:,:), allocatable :: rH_x, rH_y, rH_z
real(rprec), managed, save, dimension(:,:), allocatable :: rtopw, rbottomw
real(rprec), managed, save, dimension(:,:,:), allocatable :: RHS_col
real(rprec), managed, save, dimension(:,:,:), allocatable :: a, b, c
real(rprec), managed, save, dimension(:,:,:), allocatable :: press_gam
real(rprec), managed, save, dimension(:,:,:), allocatable :: press_bet_inv
real(rprec), managed, save, dimension(:,:), allocatable :: press_bet
real(rprec), managed, save, dimension(:,:), allocatable :: press_kx, press_ky
real(rprec), managed, save, dimension(:), allocatable :: press_halo_send_up
real(rprec), managed, save, dimension(:), allocatable :: press_halo_recv_down
real(rprec), device, save, dimension(:), allocatable :: press_rhs_send_d
real(rprec), device, save, dimension(:), allocatable :: press_rhs_recv_d
type(cudaEvent), save :: press_evt_start, press_evt_stop
#else
! CPU build (no ENABLE_CUDA): plain host versions of the work arrays. The
! diagnostic scalars/flags below are shared by both branches.
real(rprec), save, dimension(:,:,:), allocatable :: rH_x, rH_y, rH_z
real(rprec), save, dimension(:,:), allocatable :: rtopw, rbottomw
real(rprec), save, dimension(:,:,:), allocatable :: RHS_col
real(rprec), save, dimension(:,:,:), allocatable :: a, b, c
real(rprec), save, dimension(:,:,:), allocatable :: press_gam
real(rprec), save, dimension(:,:,:), allocatable :: press_bet_inv
real(rprec), save, dimension(:,:), allocatable :: press_bet
real(rprec), save, dimension(:,:), allocatable :: press_kx, press_ky
real(rprec), save, dimension(:), allocatable :: press_halo_send_up
real(rprec), save, dimension(:), allocatable :: press_halo_recv_down
#endif
integer, save :: press_fw_plan_nz1 = 0
integer, save :: press_fw_plan_1 = 0
integer, save :: press_bk_plan_nz1 = 0
integer, save :: press_bk_plan_nzp1 = 0
! cuf-path feature toggles: default ON only in the legacy ENABLE_CUDA build.
! In a plain CPU build the cuf bodies of these branches are preprocessed
! out, so taking them would SKIP real work (e.g. the p(:,:,0) halo exchange,
! posted inside the cuf tridag) -> defaults to .false. so every branch falls
! through to the original blocking CPU path.
#ifdef ENABLE_CUDA
logical, parameter :: press_cuda_paths = .true.
#else
logical, parameter :: press_cuda_paths = .false.
#endif
logical, save :: press_cuda_plans_ready = .false.
logical, save :: press_cuda_fft_checked = .false.
logical, save :: press_cuda_fft_enabled = press_cuda_paths
logical, save :: press_cuda_solver_checked = .false.
logical, save :: press_cuda_solver_enabled = press_cuda_paths
logical, save :: press_cuda_inverse_checked = .false.
logical, save :: press_cuda_inverse_enabled = press_cuda_paths
logical, save :: press_stage_checked = .false.
logical, save :: press_stage_enabled = .false.
logical, save :: press_cuda_rhs_checked = .false.
logical, save :: press_cuda_rhs_enabled = press_cuda_paths
logical, save :: press_cuda_zero_checked = .false.
logical, save :: press_cuda_zero_enabled = press_cuda_paths
logical, save :: press_extra_sync_checked = .false.
logical, save :: press_extra_sync_enabled = .false.
logical, save :: press_tridag_gpu_checked = .false.
logical, save :: press_tridag_gpu_enabled = press_cuda_paths
logical, save :: press_nb_halo_checked = .false.
logical, save :: press_nb_halo_enabled = press_cuda_paths
logical, save :: press_packed_halo_checked = .false.
logical, save :: press_packed_halo_enabled = press_cuda_paths
logical, save :: press_overlap_rhs_halo_checked = .false.
logical, save :: press_overlap_rhs_halo_enabled = .false.
logical, save :: press_rhs_halo_combined_checked = .false.
logical, save :: press_rhs_halo_combined_enabled = press_cuda_paths
logical, save :: press_rhs_halo_barrier_checked = .false.
logical, save :: press_rhs_halo_barrier_enabled = .false.
logical, save :: press_rhs_halo_audit_printed = .false.
logical, save :: press_nb_p_halo_checked = .false.
logical, save :: press_nb_p_halo_enabled = press_cuda_paths
logical, save :: press_tridag_p_halo_checked = .false.
logical, save :: press_tridag_p_halo_enabled = press_cuda_paths
logical, save :: press_tridag_coeff_ready = .false.
logical, save :: press_queue_attrib_checked = .false.
logical, save :: press_queue_attrib_enabled = .false.
logical, save :: press_rhs_diag_checked = .false.
logical, save :: press_rhs_diag_enabled = .false.
logical, save :: press_sync_after_forward_fft_checked = .false.
logical, save :: press_sync_after_forward_fft_enabled = .false.
logical, save :: press_sync_after_rhs_checked = .false.
logical, save :: press_sync_after_rhs_enabled = .false.
logical, save :: press_sync_before_transpose_checked = .false.
logical, save :: press_sync_before_transpose_enabled = .false.
logical, save :: press_queue_header_printed = .false.
logical, save :: press_queue_events_ready = .false.
integer, save :: press_stage_count = 0
integer :: press_istat
integer :: press_n_s(2), press_inem_s(2), press_onem_s(2)
integer :: press_n, press_jz
real(rprec) :: press_a_val, press_b_val, press_c_prev, press_bet_val
real(rprec) :: press_lap
real(rprec) :: press_t0, press_t1
real(rprec) :: press_stage_forward, press_stage_rhs, press_stage_tridag
real(rprec) :: press_stage_zero, press_stage_inverse
real(rprec) :: press_q_t0, press_q_pack_cpu, press_q_pack_gpu
real(rprec) :: press_q_fft_cpu, press_q_fft_gpu
real(rprec) :: press_q_rhs_prep_cpu, press_q_rhs_prep_gpu
real(rprec) :: press_q_rhs_halo_cpu, press_q_rhs_assembly_cpu
real(rprec) :: press_q_rhs_assembly_gpu, press_q_sync_after_forward
real(rprec) :: press_q_sync_after_rhs, press_q_sync_before_transpose
real(rprec) :: press_q_rhs_interior_cpu, press_q_rhs_interior_gpu
real(rprec) :: press_q_rhs_boundary_cpu, press_q_rhs_boundary_gpu
real(rprec) :: press_rhs_halo_begin, press_rhs_halo_t0, press_rhs_halo_t1
real(rprec) :: press_rhs_halo_pack, press_rhs_halo_pre_sync
real(rprec) :: press_rhs_halo_mpi, press_rhs_halo_unpack
real(rprec) :: press_rhs_halo_post_sync, press_rhs_halo_total
real(rprec) :: press_rhs_halo_ready, press_rhs_halo_ready_min
real(rprec) :: press_rhs_halo_ready_max, press_rhs_halo_ready_rel
real(rprec) :: press_rhs_halo_arrival_skew, press_rhs_halo_barrier
real(rprec) :: press_queue_wtime
integer :: press_q_sync_count
integer :: press_rhs_halo_bytes, press_rhs_halo_send_bytes
integer :: press_rhs_halo_recv_bytes, press_rhs_halo_messages
integer :: press_rhs_halo_peer, press_rhs_halo_count, press_rhs_halo_gpu
logical :: press_p_halo_pending
logical :: press_rhs_halo_pending
logical :: press_rhs_halo_overlapped
logical :: press_p_zero_plane_fix_pending
integer :: press_rhs_jz_lo, press_rhs_jz_hi
character(len=24) :: press_rhs_halo_path
character(len=24) :: press_rhs_send_type, press_rhs_recv_type
character(len=128) :: press_rhs_cuda_visible, press_rhs_mpich_gpu

logical, save :: arrays_allocated = .false.
logical, save :: press_cpu_coeff_ready = .false.

real(rprec), dimension(2) :: aH_x, aH_y
#ifdef PPMPI
integer :: press_halo_req(8)
integer :: press_halo_status(MPI_STATUS_SIZE, 8)
integer :: press_halo_nreq
integer :: press_p_halo_nreq
#endif

! Specifiy cached constants
const = 1._rprec/(nx*ny)
const2 = const/tadv1/dt
const3 = 1._rprec/(dz**2)
const4 = 1._rprec/(dz)
#ifdef ENABLE_CUDA
press_p_halo_pending = .false.
press_rhs_halo_pending = .false.
press_rhs_halo_overlapped = .false.
press_p_zero_plane_fix_pending = .false.
#endif

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( rH_x(ld,ny,lbz:nz), rH_y(ld,ny,lbz:nz), rH_z(ld,ny,lbz:nz) )
    allocate ( rtopw(ld,ny), rbottomw(ld,ny) )
    allocate ( RHS_col(ld,ny,nz+1) )
    allocate ( a(lh,ny,nz+1), b(lh,ny,nz+1), c(lh,ny,nz+1) )
#ifdef ENABLE_CUDA
    allocate ( press_gam(lh,ny,nz+1), press_bet_inv(lh,ny,nz+1),              &
        press_bet(lh,ny) )
    allocate ( press_kx(lh,ny), press_ky(lh,ny) )
    allocate ( press_halo_send_up(3*ld*ny), press_halo_recv_down(3*ld*ny) )
    allocate ( press_rhs_send_d(3*ld*ny), press_rhs_recv_d(3*ld*ny) )
    press_kx = kx
    press_ky = ky
#endif

    arrays_allocated = .true.
endif

#ifdef ENABLE_CUDA
if (.not. press_cuda_fft_checked) then
    press_cuda_fft_enabled = .true.
    press_cuda_fft_checked = .true.
end if

if (.not. press_cuda_solver_checked) then
    press_cuda_solver_enabled = press_cuda_fft_enabled
    press_cuda_solver_checked = .true.
end if

if (.not. press_cuda_inverse_checked) then
    press_cuda_inverse_enabled = press_cuda_fft_enabled
    press_cuda_inverse_checked = .true.
end if

if (.not. press_stage_checked) then
    call press_apply_env_enabled_unless_false('LESGO_PRESS_STAGE_TIMING',      &
        press_stage_enabled)
    press_stage_checked = .true.
end if

if (.not. press_cuda_rhs_checked) then
    press_cuda_rhs_enabled = press_cuda_fft_enabled
    press_cuda_rhs_checked = .true.
end if

if (.not. press_cuda_zero_checked) then
    press_cuda_zero_enabled = press_cuda_fft_enabled
    press_cuda_zero_checked = .true.
end if

if (.not. press_extra_sync_checked) then
    press_extra_sync_enabled = .false.
    press_extra_sync_checked = .true.
end if

if (.not. press_tridag_gpu_checked) then
    press_tridag_gpu_enabled = .true.
    press_tridag_gpu_checked = .true.
end if

if (.not. press_nb_halo_checked) then
    press_nb_halo_enabled = .true.
    press_nb_halo_checked = .true.
end if

if (.not. press_packed_halo_checked) then
    press_packed_halo_enabled = .true.
    press_packed_halo_checked = .true.
end if

if (.not. press_overlap_rhs_halo_checked) then
    press_overlap_rhs_halo_enabled = .false.
    press_overlap_rhs_halo_checked = .true.
end if

if (.not. press_rhs_halo_combined_checked) then
    call press_apply_env_enabled_unless_false(                                &
        'LESGO_PRESS_RHS_HALO_COMBINED', press_rhs_halo_combined_enabled)
    press_rhs_halo_combined_checked = .true.
end if

if (.not. press_rhs_halo_barrier_checked) then
    press_rhs_halo_barrier_enabled = .false.
    press_rhs_halo_barrier_checked = .true.
end if

if (.not. press_tridag_p_halo_checked) then
    press_tridag_p_halo_enabled = .true.
    press_tridag_p_halo_checked = .true.
end if

if (.not. press_nb_p_halo_checked) then
    press_nb_p_halo_enabled = .true.
    press_nb_p_halo_checked = .true.
end if

if (.not. press_queue_attrib_checked) then
    press_queue_attrib_enabled = .false.
    press_queue_attrib_checked = .true.
end if

if (.not. press_rhs_diag_checked) then
    press_rhs_diag_enabled = .false.
    press_rhs_diag_checked = .true.
end if

if (.not. press_sync_after_forward_fft_checked) then
    press_sync_after_forward_fft_enabled = .false.
    press_sync_after_forward_fft_checked = .true.
end if

if (.not. press_sync_after_rhs_checked) then
    press_sync_after_rhs_enabled = .false.
    press_sync_after_rhs_checked = .true.
end if

if (.not. press_sync_before_transpose_checked) then
    press_sync_before_transpose_enabled = .false.
    press_sync_before_transpose_checked = .true.
end if

if ((press_queue_attrib_enabled .or. press_rhs_diag_enabled) .and.            &
    .not. press_queue_events_ready) then
    press_istat = cudaEventCreate(press_evt_start)
    if (press_istat /= cudaSuccess) stop 'press queue event start create failed'
    press_istat = cudaEventCreate(press_evt_stop)
    if (press_istat /= cudaSuccess) stop 'press queue event stop create failed'
    press_queue_events_ready = .true.
end if

if (press_cuda_fft_enabled .and. .not. press_cuda_plans_ready) then
    press_n_s(1) = ny
    press_n_s(2) = nx
    press_inem_s(1) = ny
    press_inem_s(2) = ld
    press_onem_s(1) = ny
    press_onem_s(2) = ld / 2

    press_istat = cufftPlanMany(press_fw_plan_nz1, 2, press_n_s,              &
        press_inem_s, 1, ld*ny, press_onem_s, 1, (ld/2)*ny, CUFFT_D2Z, nz-1)
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array cuFFT plan nz-1 failed: ', press_istat
        stop
    end if

    press_istat = cufftPlanMany(press_fw_plan_1, 2, press_n_s,                &
        press_inem_s, 1, ld*ny, press_onem_s, 1, (ld/2)*ny, CUFFT_D2Z, 1)
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array cuFFT plan 1 failed: ', press_istat
        stop
    end if

    press_istat = cufftPlanMany(press_bk_plan_nz1, 2, press_n_s,              &
        press_onem_s, 1, (ld/2)*ny, press_inem_s, 1, ld*ny, CUFFT_Z2D, nz-1)
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array cuFFT inverse plan nz-1 failed: ', press_istat
        stop
    end if

    press_istat = cufftPlanMany(press_bk_plan_nzp1, 2, press_n_s,             &
        press_onem_s, 1, (ld/2)*ny, press_inem_s, 1, ld*ny, CUFFT_Z2D, nz+1)
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array cuFFT inverse plan nz+1 failed: ', press_istat
        stop
    end if

    press_cuda_plans_ready = .true.
end if

if (press_cuda_solver_enabled .and. nproc == 1 .and.                         &
    .not. press_tridag_coeff_ready) then
    press_n = nz + 1
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, lh-1
        if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
            press_lap = press_kx(jx,jy)**2 + press_ky(jx,jy)**2
            press_bet_val = -1._rprec
            press_bet_inv(jx,jy,1) = 1._rprec/press_bet_val

            do press_jz = 2, press_n
                if (press_jz == 2) then
                    press_c_prev = 1._rprec
                else
                    press_c_prev = const3
                end if
                if (press_jz == press_n) then
                    press_a_val = -1._rprec
                    press_b_val = 1._rprec
                else
                    press_a_val = const3
                    press_b_val = -(press_lap + 2._rprec*const3)
                end if

                press_gam(jx,jy,press_jz) = press_c_prev *                    &
                    press_bet_inv(jx,jy,press_jz-1)
                press_bet_val = press_b_val - press_a_val *                   &
                    press_gam(jx,jy,press_jz)
                press_bet_inv(jx,jy,press_jz) = 1._rprec/press_bet_val
            end do
        end if
    end do
    end do
    call press_cuda_sync('pressure tridiagonal coefficients')
    press_tridag_coeff_ready = .true.
end if
#endif

#ifdef ENABLE_CUDA
if (.not. (press_cuda_solver_enabled .and. nproc == 1)) then
#endif
    if (coord == 0) then
        p(:,:,0) = 0._rprec
#ifdef PPSAFETYMODE
    else
        p(:,:,0) = BOGUS
#endif
    end if
#ifdef ENABLE_CUDA
end if
#endif

#ifdef ENABLE_CUDA
press_q_pack_cpu = 0._rprec
press_q_pack_gpu = 0._rprec
press_q_fft_cpu = 0._rprec
press_q_fft_gpu = 0._rprec
press_q_rhs_prep_cpu = 0._rprec
press_q_rhs_prep_gpu = 0._rprec
press_q_rhs_halo_cpu = 0._rprec
press_q_rhs_assembly_cpu = 0._rprec
press_q_rhs_assembly_gpu = 0._rprec
press_q_rhs_interior_cpu = 0._rprec
press_q_rhs_interior_gpu = 0._rprec
press_q_rhs_boundary_cpu = 0._rprec
press_q_rhs_boundary_gpu = 0._rprec
press_q_sync_after_forward = 0._rprec
press_q_sync_after_rhs = 0._rprec
press_q_sync_before_transpose = 0._rprec
press_q_sync_count = 0
if (press_stage_enabled) then
    press_stage_forward = 0._rprec
    press_stage_rhs = 0._rprec
    press_stage_tridag = 0._rprec
    press_stage_zero = 0._rprec
    press_stage_inverse = 0._rprec
    press_stage_count = press_stage_count + 1
#ifdef PPMPI
    press_t0 = mpi_wtime()
#else
    call cpu_time(press_t0)
#endif
elseif (press_queue_attrib_enabled) then
    press_stage_count = press_stage_count + 1
end if
#endif

! Get the right hand side ready
! Loop over levels
! Recall that the old timestep guys already contain the pressure
#ifdef ENABLE_CUDA
if (press_cuda_fft_enabled) then
    if (press_queue_attrib_enabled) then
        call press_queue_event_start(press_evt_start, press_q_t0)
    end if
    call press_pack_rhs_cuda(u, v, w, rH_x, rH_y, rH_z, const2)
    if (press_queue_attrib_enabled) then
        call press_queue_event_stop(press_evt_start, press_evt_stop,           &
            press_q_t0, press_q_pack_cpu, press_q_pack_gpu, press_q_sync_count,&
            'forward pack')
        call press_queue_event_start(press_evt_start, press_q_t0)
    end if

    press_istat = cufftExecD2Z(press_fw_plan_nz1, rH_x(:,:,1), rH_x(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array rH_x forward failed: ', press_istat
        stop
    end if
    press_istat = cufftExecD2Z(press_fw_plan_nz1, rH_y(:,:,1), rH_y(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array rH_y forward failed: ', press_istat
        stop
    end if
    press_istat = cufftExecD2Z(press_fw_plan_nz1, rH_z(:,:,1), rH_z(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array rH_z forward failed: ', press_istat
        stop
    end if
else
#endif
do jz = 1, nz-1
    rH_x(:,:,jz) = const2 * u(:,:,jz)
    rH_y(:,:,jz) = const2 * v(:,:,jz)
    rH_z(:,:,jz) = const2 * w(:,:,jz)

    call dfftw_execute_dft_r2c(forw, rH_x(:,:,jz), rH_x(:,:,jz))
    call dfftw_execute_dft_r2c(forw, rH_y(:,:,jz), rH_y(:,:,jz))
    call dfftw_execute_dft_r2c(forw, rH_z(:,:,jz), rH_z(:,:,jz))
end do
#ifdef ENABLE_CUDA
end if
#endif

#if defined(PPMPI) && defined(PPSAFETYMODE)
#ifdef ENABLE_CUDA
if (.not. (press_cuda_solver_enabled .and. nproc == 1)) then
#endif
  !Careful - only update real values (odd indicies)
  rH_x(1:ld:2,:,0) = BOGUS
  rH_y(1:ld:2,:,0) = BOGUS
  rH_z(1:ld:2,:,0) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef ENABLE_CUDA
if (.not. (press_cuda_solver_enabled .and. nproc == 1)) then
#endif
!Careful - only update real values (odd indicies)
rH_x(1:ld:2,:,nz) = BOGUS
rH_y(1:ld:2,:,nz) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif

#ifdef PPMPI
if (coord == nproc-1) then
#ifdef ENABLE_CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
    do jy = 1, ny
    do jx = 1, ld
        rH_z(jx,jy,nz) = const2 * w(jx,jy,nz)
    end do
    end do
#ifdef ENABLE_CUDA
    if (press_cuda_fft_enabled) then
        press_istat = cufftExecD2Z(press_fw_plan_1, rH_z(:,:,nz), rH_z(:,:,nz))
        if (press_istat /= CUFFT_SUCCESS) then
            print *, 'press_stag_array rH_z top forward failed: ', press_istat
            stop
        end if
    else
#endif
        call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,jz))
#ifdef ENABLE_CUDA
    end if
#endif
#ifdef PPSAFETYMODE
else
    rH_z(1:ld:2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
#endif
endif
#else
#ifdef ENABLE_CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
do jy = 1, ny
do jx = 1, ld
    rH_z(jx,jy,nz) = const2 * w(jx,jy,nz)
end do
end do
#ifdef ENABLE_CUDA
if (press_cuda_fft_enabled) then
    press_istat = cufftExecD2Z(press_fw_plan_1, rH_z(:,:,nz), rH_z(:,:,nz))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array rH_z top forward failed: ', press_istat
        stop
    end if
else
#endif
call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,jz))
#ifdef ENABLE_CUDA
end if
#endif
#endif

if (coord == 0) then
#ifdef ENABLE_CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
    do jy = 1, ny
    do jx = 1, ld
        rbottomw(jx,jy) = const * divtz(jx,jy,1)
    end do
    end do
#ifdef ENABLE_CUDA
    if (press_cuda_fft_enabled) then
        press_istat = cufftExecD2Z(press_fw_plan_1, rbottomw, rbottomw)
        if (press_istat /= CUFFT_SUCCESS) then
            print *, 'press_stag_array bottom forward failed: ', press_istat
            stop
        end if
    else
#endif
        call dfftw_execute_dft_r2c(forw, rbottomw, rbottomw )
#ifdef ENABLE_CUDA
    end if
#endif
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
#ifdef ENABLE_CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
    do jy = 1, ny
    do jx = 1, ld
        rtopw(jx,jy) = const * divtz(jx,jy,nz)
    end do
    end do
#ifdef ENABLE_CUDA
    if (press_cuda_fft_enabled) then
        press_istat = cufftExecD2Z(press_fw_plan_1, rtopw, rtopw)
        if (press_istat /= CUFFT_SUCCESS) then
            print *, 'press_stag_array top forward failed: ', press_istat
            stop
        end if
    else
#endif
        call dfftw_execute_dft_r2c(forw, rtopw, rtopw)
#ifdef ENABLE_CUDA
    end if
#endif
#ifdef PPMPI
endif
#endif

#ifdef ENABLE_CUDA
if (press_queue_attrib_enabled .and. press_cuda_fft_enabled) then
    call press_queue_event_stop(press_evt_start, press_evt_stop, press_q_t0,   &
        press_q_fft_cpu, press_q_fft_gpu, press_q_sync_count, 'forward FFT')
end if
#endif
if (press_sync_after_forward_fft_enabled) then
    press_q_t0 = press_queue_wtime()
    call press_cuda_sync('pressure forced sync after forward FFT')
    press_q_sync_after_forward = press_queue_wtime() - press_q_t0
    press_q_sync_count = press_q_sync_count + 1
end if


#ifdef ENABLE_CUDA
if (press_cuda_solver_enabled .and. nproc == 1) then

    ! Finish pressure solve on the GPU for the validated single-rank path.
    if (press_stage_enabled) then
        call press_cuda_sync('pressure single-rank forward timing')
#ifdef PPMPI
        press_t1 = mpi_wtime()
#else
        call cpu_time(press_t1)
#endif
        press_stage_forward = press_t1 - press_t0
        press_t0 = press_t1
    end if
!$cuf kernel do(3) <<<*,*>>>
    do jz = lbz, nz
    do jy = 1, ny
    do jx = 1, ld
        if ((jz >= 1) .and. (jz <= nz-1)) then
            if ((jx >= ld-1) .or. (jy == ny/2+1)) then
                rH_x(jx,jy,jz) = 0._rprec
                rH_y(jx,jy,jz) = 0._rprec
                rH_z(jx,jy,jz) = 0._rprec
            end if
        end if
        if ((jz == nz) .and. (jx >= ld-1 .or. jy == ny/2+1)) then
            rH_z(jx,jy,nz) = 0._rprec
        end if
    end do
    end do
    end do

!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        if ((jx >= ld-1) .or. (jy == ny/2+1)) then
            rtopw(jx,jy) = 0._rprec
            rbottomw(jx,jy) = 0._rprec
        end if
        RHS_col(jx,jy,1) = -dz*rbottomw(jx,jy)
        RHS_col(jx,jy,nz+1) = -dz*rtopw(jx,jy)
    end do
    end do

    jz_min = 2
!$cuf kernel do(3) <<<*,*>>>
    do jz = jz_min, nz
    do jy = 1, ny
    do jx = 1, lh-1
        if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
            ii = 2*jx
            ir = ii - 1

            RHS_col(ir,jy,jz) = -rH_x(ii,jy,jz-1)*press_kx(jx,jy)              &
                - rH_y(ii,jy,jz-1)*press_ky(jx,jy)                             &
                + (rH_z(ir,jy,jz) - rH_z(ir,jy,jz-1))*const4
            RHS_col(ii,jy,jz) =  rH_x(ir,jy,jz-1)*press_kx(jx,jy)              &
                + rH_y(ir,jy,jz-1)*press_ky(jx,jy)                             &
                + (rH_z(ii,jy,jz) - rH_z(ii,jy,jz-1))*const4
        end if
    end do
    end do
    end do
    if (press_stage_enabled) then
        call press_cuda_sync('pressure single-rank RHS timing')
#ifdef PPMPI
        press_t1 = mpi_wtime()
#else
        call cpu_time(press_t1)
#endif
        press_stage_rhs = press_t1 - press_t0
        press_t0 = press_t1
    end if
    press_n = nz + 1
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, lh-1
        if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
            ii = 2*jx
            ir = ii - 1
            p(ir,jy,0) = RHS_col(ir,jy,1)*press_bet_inv(jx,jy,1)
            p(ii,jy,0) = RHS_col(ii,jy,1)*press_bet_inv(jx,jy,1)

            do press_jz = 2, press_n
                if (press_jz == press_n) then
                    press_a_val = -1._rprec
                else
                    press_a_val = const3
                end if
                p(ir,jy,press_jz-1) = (RHS_col(ir,jy,press_jz) -              &
                    press_a_val*p(ir,jy,press_jz-2)) *                         &
                    press_bet_inv(jx,jy,press_jz)
                p(ii,jy,press_jz-1) = (RHS_col(ii,jy,press_jz) -              &
                    press_a_val*p(ii,jy,press_jz-2)) *                         &
                    press_bet_inv(jx,jy,press_jz)
            end do

            do press_jz = press_n-1, 1, -1
                p(ir,jy,press_jz-1) = p(ir,jy,press_jz-1) -                   &
                    press_gam(jx,jy,press_jz+1)*p(ir,jy,press_jz)
                p(ii,jy,press_jz-1) = p(ii,jy,press_jz-1) -                   &
                    press_gam(jx,jy,press_jz+1)*p(ii,jy,press_jz)
            end do
        end if
    end do
    end do
    if (press_stage_enabled) then
        call press_cuda_sync('pressure single-rank Thomas timing')
#ifdef PPMPI
        press_t1 = mpi_wtime()
#else
        call cpu_time(press_t1)
#endif
        press_stage_tridag = press_t1 - press_t0
        press_t0 = press_t1
    end if
!$cuf kernel do(1) <<<*,*>>>
    do jx = 1, 2
        p(jx,1,0) = 0._rprec
        p(jx,1,1) = -dz*rbottomw(jx,1)
        do jz = 2, nz
            p(jx,1,jz) = p(jx,1,jz-1) + rH_z(jx,1,jz)*dz
        end do
    end do

!$cuf kernel do(3) <<<*,*>>>
    do jz = 0, nz
    do jy = 1, ny
    do jx = 1, ld
        if ((jx >= ld-1) .or. (jy == ny/2+1)) p(jx,jy,jz) = 0._rprec
    end do
    end do
    end do
    if (press_stage_enabled) then
        call press_cuda_sync('pressure single-rank zero-mode timing')
#ifdef PPMPI
        press_t1 = mpi_wtime()
#else
        call cpu_time(press_t1)
#endif
        press_stage_zero = press_t1 - press_t0
        press_t0 = press_t1
    end if

!$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz-1
    do jy = 1, ny
    do jx = 1, lh
        ii = 2*jx
        ir = ii - 1
        dpdx(ir,jy,jz) = -p(ii,jy,jz)*press_kx(jx,jy)
        dpdx(ii,jy,jz) =  p(ir,jy,jz)*press_kx(jx,jy)
        dpdy(ir,jy,jz) = -p(ii,jy,jz)*press_ky(jx,jy)
        dpdy(ii,jy,jz) =  p(ir,jy,jz)*press_ky(jx,jy)
    end do
    end do
    end do
    press_istat = cufftExecZ2D(press_bk_plan_nz1, dpdx(:,:,1), dpdx(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array dpdx inverse failed: ', press_istat
        stop
    end if
    press_istat = cufftExecZ2D(press_bk_plan_nz1, dpdy(:,:,1), dpdy(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array dpdy inverse failed: ', press_istat
        stop
    end if
    press_istat = cufftExecZ2D(press_bk_plan_nzp1, p(:,:,0), p(:,:,0))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array pressure inverse failed: ', press_istat
        stop
    end if
!$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, nx
        if (jz < nz) then
            dpdz(jx,jy,jz) = (p(jx,jy,jz) - p(jx,jy,jz-1))/dz
        else
            dpdz(jx,jy,jz) = (p(jx,jy,nz) - p(jx,jy,nz-1))/dz
        end if
    end do
    end do
    end do
    call press_cuda_sync('pressure dpdz')
    if (press_stage_enabled) then
#ifdef PPMPI
        press_t1 = mpi_wtime()
#else
        call cpu_time(press_t1)
#endif
        press_stage_inverse = press_t1 - press_t0
        call press_stage_report(press_stage_count, press_stage_forward,       &
            press_stage_rhs, press_stage_tridag, press_stage_zero,            &
            press_stage_inverse)
    end if

else
#endif

#ifdef ENABLE_CUDA
if (press_cuda_fft_enabled) then
    if (.not. (press_cuda_rhs_enabled .and. nproc > 1 .and.                  &
        .not. press_extra_sync_enabled)) then
        call press_cuda_sync('pressure forward transforms before CPU solve')
    end if
    if (press_stage_enabled) then
#ifdef PPMPI
        press_t1 = mpi_wtime()
#else
        call cpu_time(press_t1)
#endif
        press_stage_forward = press_t1 - press_t0
        press_t0 = press_t1
    end if
end if
#endif

#ifdef ENABLE_CUDA
if (press_cuda_rhs_enabled .and. nproc > 1) then

    if (coord == 0) then
        jz_min = 2
    else
        jz_min = 1
    end if

    if (.not. press_cpu_coeff_ready) then
        if (coord == 0) then
#ifdef PPSAFETYMODE
            a(:,:,1) = BOGUS
#endif
            b(:,:,1) = -1._rprec
            c(:,:,1) = 1._rprec
        end if

#ifdef PPMPI
        if (coord == nproc-1) then
#endif
            a(:,:,nz+1) = -1._rprec
            b(:,:,nz+1) = 1._rprec
#ifdef PPSAFETYMODE
            c(:,:,nz+1) = BOGUS
#endif
#ifdef PPMPI
        endif
#endif

        do jz = jz_min, nz
        do jy = 1, ny
            if (jy == ny/2 + 1) cycle
            do jx = 1, lh-1
                if (jx*jy == 1) cycle
                a(jx, jy, jz) = const3
                b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 +              &
                    2._rprec*const3)
                c(jx, jy, jz) = const3
            end do
        end do
        end do

        press_cpu_coeff_ready = .true.
    end if

    if (press_queue_attrib_enabled) then
        call press_queue_event_start(press_evt_start, press_q_t0)
    end if
    call press_rhs_prep_cuda(rH_x, rH_y, rH_z, rtopw, rbottomw, RHS_col)
    if (press_queue_attrib_enabled) then
        call press_queue_event_stop(press_evt_start, press_evt_stop,           &
            press_q_t0, press_q_rhs_prep_cpu, press_q_rhs_prep_gpu,           &
            press_q_sync_count, 'RHS prep')
        press_q_t0 = press_queue_wtime()
    end if

#ifdef PPMPI
    if (press_packed_halo_enabled) then
        press_rhs_halo_begin = press_queue_wtime()
        press_rhs_halo_pack = 0._rprec
        press_rhs_halo_pre_sync = 0._rprec
        press_rhs_halo_mpi = 0._rprec
        press_rhs_halo_unpack = 0._rprec
        press_rhs_halo_post_sync = 0._rprec
        press_rhs_halo_barrier = 0._rprec
        press_rhs_halo_ready_rel = 0._rprec
        press_rhs_halo_arrival_skew = 0._rprec
        press_rhs_halo_send_bytes = 0
        press_rhs_halo_recv_bytes = 0
        press_rhs_halo_bytes = 0
        press_rhs_halo_messages = 0
        press_rhs_halo_path = 'packed-old'

        if (press_rhs_halo_combined_enabled .and. nproc == 2) then
            press_rhs_halo_path = 'combined-dev'
            press_rhs_halo_peer = 1 - coord
            press_rhs_halo_count = 3 * ld * ny
            press_rhs_halo_send_bytes = press_rhs_halo_count * 8
            press_rhs_halo_recv_bytes = press_rhs_halo_count * 8
            press_rhs_halo_bytes = press_rhs_halo_send_bytes +                 &
                press_rhs_halo_recv_bytes
            press_rhs_halo_messages = 1

            if (.not. press_rhs_halo_audit_printed) then
                if (coord == 0) write(*,*)                                    &
                    'Pressure RHS halo path: combined device nproc2=',         &
                    press_rhs_halo_combined_enabled
                if (press_rhs_diag_enabled .or. press_queue_attrib_enabled) then
                    press_rhs_cuda_visible = '-'
                    press_rhs_mpich_gpu = '-'
                    call get_environment_variable('CUDA_VISIBLE_DEVICES',      &
                        press_rhs_cuda_visible)
                    call get_environment_variable('MPICH_GPU_SUPPORT_ENABLED', &
                        press_rhs_mpich_gpu)
                    press_istat = cudaGetDevice(press_rhs_halo_gpu)
                    press_rhs_send_type = 'device'
                    press_rhs_recv_type = 'device'
                    call press_rhs_halo_audit('press_stag_array.f90',         &
                        'press_stag_array', trim(press_rhs_halo_path),         &
                        'MPI_Sendrecv', press_rhs_halo_messages,              &
                        press_rhs_halo_send_bytes, press_rhs_halo_recv_bytes,  &
                        press_rhs_halo_peer, trim(press_rhs_send_type),        &
                        trim(press_rhs_recv_type), 'yes', 'no',               &
                        press_rhs_halo_gpu, trim(press_rhs_cuda_visible),      &
                        trim(press_rhs_mpich_gpu))
                end if
                press_rhs_halo_audit_printed = .true.
            end if

            press_rhs_halo_t0 = press_queue_wtime()
            call press_pack_rhs_halo_combined_cuda(rH_x, rH_y, rH_z,           &
                press_rhs_send_d, coord)
            press_rhs_halo_pack = press_queue_wtime() - press_rhs_halo_t0

            press_rhs_halo_t0 = press_queue_wtime()
            call press_cuda_sync('pressure combined rhs halo pack')
            press_rhs_halo_pre_sync = press_queue_wtime() - press_rhs_halo_t0

            if (press_rhs_diag_enabled .or. press_queue_attrib_enabled) then
                press_rhs_halo_ready = press_queue_wtime()
                call mpi_allreduce(press_rhs_halo_ready,                      &
                    press_rhs_halo_ready_min, 1, MPI_RPREC, MPI_MIN, comm,    &
                    ierr)
                call mpi_allreduce(press_rhs_halo_ready,                      &
                    press_rhs_halo_ready_max, 1, MPI_RPREC, MPI_MAX, comm,    &
                    ierr)
                press_rhs_halo_ready_rel = press_rhs_halo_ready -              &
                    press_rhs_halo_ready_min
                press_rhs_halo_arrival_skew = press_rhs_halo_ready_max -       &
                    press_rhs_halo_ready_min
            end if

            if ((press_rhs_diag_enabled .or. press_queue_attrib_enabled) .and. &
                press_rhs_halo_barrier_enabled) then
                press_rhs_halo_t0 = press_queue_wtime()
                call mpi_barrier(comm, ierr)
                press_rhs_halo_barrier = press_queue_wtime() -                 &
                    press_rhs_halo_t0
            end if

            press_rhs_halo_t0 = press_queue_wtime()
            call mpi_sendrecv(press_rhs_send_d(1), press_rhs_halo_count,       &
                MPI_RPREC, press_rhs_halo_peer, 1101, press_rhs_recv_d(1),     &
                press_rhs_halo_count, MPI_RPREC, press_rhs_halo_peer, 1101,    &
                comm, status, ierr)
            press_rhs_halo_mpi = press_queue_wtime() - press_rhs_halo_t0

            press_rhs_halo_t0 = press_queue_wtime()
            call press_unpack_rhs_halo_combined_cuda(press_rhs_recv_d, rH_x,   &
                rH_y, rH_z, coord)
            press_rhs_halo_unpack = press_queue_wtime() - press_rhs_halo_t0

            if (press_queue_attrib_enabled .or. press_rhs_diag_enabled) then
                press_rhs_halo_t0 = press_queue_wtime()
                call press_cuda_sync('pressure combined rhs halo unpack')
                press_rhs_halo_post_sync = press_queue_wtime() -               &
                    press_rhs_halo_t0
                press_q_sync_count = press_q_sync_count + 1
            end if
        else
            if (coord < nproc - 1) then
                press_rhs_halo_t0 = press_queue_wtime()
                call press_pack_rhs_halo_cuda(rH_x, rH_y, rH_z,                &
                    press_halo_send_up)
                press_rhs_halo_pack = press_queue_wtime() - press_rhs_halo_t0
                press_rhs_halo_t0 = press_queue_wtime()
                call press_cuda_sync('pressure packed rhs halo pack')
                press_rhs_halo_pre_sync = press_queue_wtime() -                &
                    press_rhs_halo_t0
            else
                press_rhs_halo_t0 = press_queue_wtime()
                call press_cuda_sync('pressure rhs prep before halos')
                press_rhs_halo_pre_sync = press_queue_wtime() -                &
                    press_rhs_halo_t0
            end if

            if (.not. press_rhs_halo_audit_printed) then
                if (coord == 0) write(*,*)                                    &
                    'Pressure RHS halo path: legacy managed/asymmetric'
                if (coord == 0) then
                    press_rhs_halo_send_bytes = 3 * ld * ny * 8
                    press_rhs_halo_recv_bytes = ld * ny * 8
                    press_rhs_halo_peer = up
                else
                    press_rhs_halo_send_bytes = ld * ny * 8
                    press_rhs_halo_recv_bytes = 3 * ld * ny * 8
                    press_rhs_halo_peer = down
                end if
                if (press_rhs_diag_enabled .or. press_queue_attrib_enabled) then
                    press_rhs_cuda_visible = '-'
                    press_rhs_mpich_gpu = '-'
                    call get_environment_variable('CUDA_VISIBLE_DEVICES',      &
                        press_rhs_cuda_visible)
                    call get_environment_variable('MPICH_GPU_SUPPORT_ENABLED', &
                        press_rhs_mpich_gpu)
                    press_istat = cudaGetDevice(press_rhs_halo_gpu)
                    press_rhs_send_type = 'managed'
                    press_rhs_recv_type = 'managed'
                    call press_rhs_halo_audit('press_stag_array.f90',         &
                        'press_stag_array', trim(press_rhs_halo_path),         &
                        'Irecv/Isend/Waitall', 2, press_rhs_halo_send_bytes,   &
                        press_rhs_halo_recv_bytes, press_rhs_halo_peer,        &
                        trim(press_rhs_send_type), trim(press_rhs_recv_type),  &
                        'mixed', 'yes', press_rhs_halo_gpu,                   &
                        trim(press_rhs_cuda_visible), trim(press_rhs_mpich_gpu))
                end if
                press_rhs_halo_audit_printed = .true.
            end if

            if (press_rhs_diag_enabled .or. press_queue_attrib_enabled) then
                press_rhs_halo_ready = press_queue_wtime()
                call mpi_allreduce(press_rhs_halo_ready,                      &
                    press_rhs_halo_ready_min, 1, MPI_RPREC, MPI_MIN, comm,    &
                    ierr)
                call mpi_allreduce(press_rhs_halo_ready,                      &
                    press_rhs_halo_ready_max, 1, MPI_RPREC, MPI_MAX, comm,    &
                    ierr)
                press_rhs_halo_ready_rel = press_rhs_halo_ready -              &
                    press_rhs_halo_ready_min
                press_rhs_halo_arrival_skew = press_rhs_halo_ready_max -       &
                    press_rhs_halo_ready_min
            end if
            if ((press_rhs_diag_enabled .or. press_queue_attrib_enabled) .and. &
                press_rhs_halo_barrier_enabled) then
                press_rhs_halo_t0 = press_queue_wtime()
                call mpi_barrier(comm, ierr)
                press_rhs_halo_barrier = press_queue_wtime() -                 &
                    press_rhs_halo_t0
            end if

            press_halo_req = MPI_REQUEST_NULL
            press_halo_nreq = 0
            press_rhs_halo_t0 = press_queue_wtime()
            if (coord > 0) then
                press_halo_nreq = press_halo_nreq + 1
                call mpi_irecv(press_halo_recv_down(1), 3*ld*ny, MPI_RPREC,   &
                    down, 11, comm, press_halo_req(press_halo_nreq), ierr)
            end if
            if (coord < nproc - 1) then
                press_halo_nreq = press_halo_nreq + 1
                call mpi_irecv(rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 12,      &
                    comm, press_halo_req(press_halo_nreq), ierr)
                press_halo_nreq = press_halo_nreq + 1
                call mpi_isend(press_halo_send_up(1), 3*ld*ny, MPI_RPREC, up, &
                    11, comm, press_halo_req(press_halo_nreq), ierr)
            end if
            if (coord > 0) then
                press_halo_nreq = press_halo_nreq + 1
                call mpi_isend(rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 12,     &
                    comm, press_halo_req(press_halo_nreq), ierr)
            end if
            if (press_overlap_rhs_halo_enabled) then
                press_rhs_jz_lo = jz_min
                press_rhs_jz_hi = nz
                if (coord > 0) press_rhs_jz_lo = max(press_rhs_jz_lo, 2)
                if (coord < nproc - 1) press_rhs_jz_hi = min(press_rhs_jz_hi,  &
                    nz-1)
                if (press_rhs_jz_lo <= press_rhs_jz_hi) then
                    call press_assemble_rhs_range_cuda(rH_x, rH_y, rH_z,       &
                        RHS_col, press_kx, press_ky, const4, press_rhs_jz_lo, &
                        press_rhs_jz_hi)
                end if
                press_rhs_halo_pending = press_halo_nreq > 0
                press_rhs_halo_overlapped = .true.
            else
                if (press_halo_nreq > 0) then
                    call mpi_waitall(press_halo_nreq, press_halo_req,          &
                        press_halo_status, ierr)
                end if
                press_rhs_halo_mpi = press_queue_wtime() - press_rhs_halo_t0
                if (coord /= 0) then
                    press_rhs_halo_t0 = press_queue_wtime()
                    call press_unpack_rhs_halo_cuda(press_halo_recv_down,      &
                        rH_x, rH_y, rH_z)
                    press_rhs_halo_unpack = press_queue_wtime() -              &
                        press_rhs_halo_t0
                    if (press_queue_attrib_enabled .or.                       &
                        press_rhs_diag_enabled) then
                        press_rhs_halo_t0 = press_queue_wtime()
                        call press_cuda_sync('pressure packed rhs halo unpack')
                        press_rhs_halo_post_sync = press_queue_wtime() -       &
                            press_rhs_halo_t0
                        press_q_sync_count = press_q_sync_count + 1
                    end if
                endif
            end if

            if (coord == 0) then
                press_rhs_halo_send_bytes = 3 * ld * ny * 8
                press_rhs_halo_recv_bytes = ld * ny * 8
            else
                press_rhs_halo_send_bytes = ld * ny * 8
                press_rhs_halo_recv_bytes = 3 * ld * ny * 8
            end if
            press_rhs_halo_bytes = press_rhs_halo_send_bytes +                 &
                press_rhs_halo_recv_bytes
            press_rhs_halo_messages = 2
        end if

        press_rhs_halo_total = press_queue_wtime() - press_rhs_halo_begin
        if (press_queue_attrib_enabled .or. press_rhs_diag_enabled) then
            call press_rhs_halo_report(press_stage_count,                      &
                trim(press_rhs_halo_path), press_rhs_halo_pre_sync,            &
                press_rhs_halo_pack, press_rhs_halo_mpi,                       &
                press_rhs_halo_unpack, press_rhs_halo_post_sync,               &
                press_rhs_halo_total, press_rhs_halo_send_bytes,               &
                press_rhs_halo_recv_bytes, press_rhs_halo_messages,            &
                press_rhs_halo_ready_rel, press_rhs_halo_arrival_skew,         &
                press_rhs_halo_barrier)
        end if
    elseif (press_nb_halo_enabled) then
        call press_cuda_sync('pressure rhs prep before halos')
        press_halo_req = MPI_REQUEST_NULL
        press_halo_nreq = 0
        if (coord > 0) then
            press_halo_nreq = press_halo_nreq + 1
            call mpi_irecv(rH_x(1, 1, 0), ld*ny, MPI_RPREC, down, 1, comm,    &
                press_halo_req(press_halo_nreq), ierr)
            press_halo_nreq = press_halo_nreq + 1
            call mpi_irecv(rH_y(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm,    &
                press_halo_req(press_halo_nreq), ierr)
            press_halo_nreq = press_halo_nreq + 1
            call mpi_irecv(rH_z(1, 1, 0), ld*ny, MPI_RPREC, down, 3, comm,    &
                press_halo_req(press_halo_nreq), ierr)
            press_halo_nreq = press_halo_nreq + 1
            call mpi_isend(rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 6, comm,    &
                press_halo_req(press_halo_nreq), ierr)
        end if
        if (coord < nproc - 1) then
            press_halo_nreq = press_halo_nreq + 1
            call mpi_irecv(rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 6, comm,     &
                press_halo_req(press_halo_nreq), ierr)
            press_halo_nreq = press_halo_nreq + 1
            call mpi_isend(rH_x(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1, comm,   &
                press_halo_req(press_halo_nreq), ierr)
            press_halo_nreq = press_halo_nreq + 1
            call mpi_isend(rH_y(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2, comm,   &
                press_halo_req(press_halo_nreq), ierr)
            press_halo_nreq = press_halo_nreq + 1
            call mpi_isend(rH_z(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3, comm,   &
                press_halo_req(press_halo_nreq), ierr)
        end if
        if (press_halo_nreq > 0) then
            call mpi_waitall(press_halo_nreq, press_halo_req,                 &
                press_halo_status, ierr)
        end if
    else
        call press_cuda_sync('pressure rhs prep before halos')
        call mpi_dbg_sendrecv_r (rH_x(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,    &
            rH_x(1, 1, 0), ld*ny, MPI_RPREC, down, 1, comm, status, ierr,     &
            'press_rH_x_up')
        call mpi_dbg_sendrecv_r (rH_y(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,    &
            rH_y(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr,     &
            'press_rH_y_up')
        call mpi_dbg_sendrecv_r (rH_z(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,    &
            rH_z(1, 1, 0), ld*ny, MPI_RPREC, down, 3, comm, status, ierr,     &
            'press_rH_z_up')
        call mpi_dbg_sendrecv_r (rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 6,     &
            rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 6, comm, status, ierr,      &
            'press_rH_z_down')
    end if
    if (press_queue_attrib_enabled) then
        if (press_packed_halo_enabled) then
            press_q_rhs_halo_cpu = press_rhs_halo_total
        else
            press_q_rhs_halo_cpu = press_queue_wtime() - press_q_t0
        end if
    end if
#endif

    if (press_extra_sync_enabled) then
        call press_cuda_sync('pressure rhs halos before assembly')
    end if
    if (press_rhs_halo_overlapped) then
        if (press_queue_attrib_enabled .or. press_rhs_diag_enabled) then
            call press_queue_event_start(press_evt_start, press_q_t0)
        end if
#ifdef PPMPI
        if (press_rhs_halo_pending) then
            call mpi_waitall(press_halo_nreq, press_halo_req,                 &
                press_halo_status, ierr)
            press_rhs_halo_pending = .false.
        end if
#endif
        if (coord /= 0) then
            call press_unpack_rhs_halo_cuda(press_halo_recv_down,              &
                rH_x, rH_y, rH_z)
            if (jz_min <= 1) then
                call press_assemble_rhs_range_cuda(rH_x, rH_y, rH_z, RHS_col, &
                    press_kx, press_ky, const4, 1, 1)
            end if
        endif
        if (coord < nproc - 1) then
            call press_assemble_rhs_range_cuda(rH_x, rH_y, rH_z, RHS_col,     &
                press_kx, press_ky, const4, nz, nz)
        end if
        if (press_queue_attrib_enabled .or. press_rhs_diag_enabled) then
            call press_queue_event_stop(press_evt_start, press_evt_stop,       &
                press_q_t0, press_q_rhs_assembly_cpu,                         &
                press_q_rhs_assembly_gpu, press_q_sync_count,                 &
                'RHS assembly overlap remainder')
            press_q_rhs_boundary_cpu = press_q_rhs_assembly_cpu
            press_q_rhs_boundary_gpu = press_q_rhs_assembly_gpu
        end if
    else
        if (press_queue_attrib_enabled .or. press_rhs_diag_enabled) then
            press_rhs_jz_lo = jz_min
            press_rhs_jz_hi = nz
            if (coord > 0) press_rhs_jz_lo = max(press_rhs_jz_lo, 2)
            if (coord < nproc - 1) press_rhs_jz_hi = min(press_rhs_jz_hi,      &
                nz-1)
            if (press_rhs_jz_lo <= press_rhs_jz_hi) then
                call press_queue_event_start(press_evt_start, press_q_t0)
                call press_assemble_rhs_range_cuda(rH_x, rH_y, rH_z, RHS_col, &
                    press_kx, press_ky, const4, press_rhs_jz_lo,              &
                    press_rhs_jz_hi)
                call press_queue_event_stop(press_evt_start, press_evt_stop,   &
                    press_q_t0, press_q_rhs_interior_cpu,                     &
                    press_q_rhs_interior_gpu, press_q_sync_count,             &
                    'RHS assembly interior')
            end if

            if (coord /= 0 .and. jz_min <= 1) then
                call press_queue_event_start(press_evt_start, press_q_t0)
                call press_assemble_rhs_range_cuda(rH_x, rH_y, rH_z, RHS_col, &
                    press_kx, press_ky, const4, 1, 1)
                call press_queue_event_stop(press_evt_start, press_evt_stop,   &
                    press_q_t0, press_q_rhs_boundary_cpu,                     &
                    press_q_rhs_boundary_gpu, press_q_sync_count,             &
                    'RHS assembly lower boundary')
            end if
            if (coord < nproc - 1) then
                call press_queue_event_start(press_evt_start, press_q_t0)
                call press_assemble_rhs_range_cuda(rH_x, rH_y, rH_z, RHS_col, &
                    press_kx, press_ky, const4, nz, nz)
                call press_queue_event_stop(press_evt_start, press_evt_stop,   &
                    press_q_t0, press_q_rhs_boundary_cpu,                     &
                    press_q_rhs_boundary_gpu, press_q_sync_count,             &
                    'RHS assembly upper boundary')
            end if
            press_q_rhs_assembly_cpu = press_q_rhs_interior_cpu +              &
                press_q_rhs_boundary_cpu
            press_q_rhs_assembly_gpu = press_q_rhs_interior_gpu +             &
                press_q_rhs_boundary_gpu
        else
            call press_assemble_rhs_cuda(rH_x, rH_y, rH_z, RHS_col, press_kx,  &
                press_ky, const4, jz_min)
        end if
    end if
    if (press_queue_attrib_enabled .or. press_rhs_diag_enabled) then
        call press_rhs_assembly_report(press_stage_count,                      &
            press_q_rhs_interior_cpu, press_q_rhs_interior_gpu,                &
            press_q_rhs_boundary_cpu, press_q_rhs_boundary_gpu,                &
            press_q_rhs_assembly_cpu, press_q_rhs_assembly_gpu)
    end if
    if (press_sync_after_rhs_enabled) then
        press_q_t0 = press_queue_wtime()
        call press_cuda_sync('pressure forced sync after RHS')
        press_q_sync_after_rhs = press_queue_wtime() - press_q_t0
        press_q_sync_count = press_q_sync_count + 1
    end if
    if (press_extra_sync_enabled .or. .not. press_tridag_gpu_enabled) then
        call press_cuda_sync('pressure rhs assembly before CPU tridag')
    end if

else
#endif

! set oddballs to 0
rH_x(ld-1:ld,:,1:nz-1) = 0._rprec
rH_y(ld-1:ld,:,1:nz-1) = 0._rprec
rH_z(ld-1:ld,:,1:nz-1) = 0._rprec
rH_x(:,ny/2+1,1:nz-1) = 0._rprec
rH_y(:,ny/2+1,1:nz-1) = 0._rprec
rH_z(:,ny/2+1,1:nz-1) = 0._rprec
! should also set to zero for rH_z (nz) on coord == nproc-1
if (coord == nproc-1) then
    rH_z(ld-1:ld,:,nz) = 0._rprec
    rH_z(:,ny/2+1,nz) = 0._rprec
end if

! with MPI; topw and bottomw are only on top & bottom processes
rtopw(ld-1:ld, :) = 0._rprec
rtopw(:, ny/2+1) = 0._rprec
rbottomw(ld-1:ld, :) = 0._rprec
rbottomw(:, ny/2+1) = 0._rprec

if (coord == 0) then
    jz_min = 2
else
  jz_min = 1
end if

if (.not. press_cpu_coeff_ready) then
    !  a, b, and c are treated as the real part of a complex array
    if (coord == 0) then
#ifdef PPSAFETYMODE
        a(:,:,1) = BOGUS
#endif
        b(:,:,1) = -1._rprec
        c(:,:,1) = 1._rprec
    end if

#ifdef PPMPI
    if (coord == nproc-1) then
#endif
        !--top nodes
        a(:,:,nz+1) = -1._rprec
        b(:,:,nz+1) = 1._rprec
#ifdef PPSAFETYMODE
        c(:,:,nz+1) = BOGUS
#endif
#ifdef PPMPI
    endif
#endif

    do jz = jz_min, nz
    do jy = 1, ny
        if (jy == ny/2 + 1) cycle

        do jx = 1, lh-1

            if (jx*jy == 1) cycle

            ! JDA dissertation, eqn(2.85) a,b,c=coefficients
            a(jx, jy, jz) = const3
            b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 + 2._rprec*const3)
            c(jx, jy, jz) = const3

        end do
    end do
    end do

    press_cpu_coeff_ready = .true.
end if

if (coord == 0) then
    RHS_col(:,:,1) = -dz * rbottomw(:,:)
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    RHS_col(:,:,nz+1) = -dz * rtopw(:,:)
#ifdef PPMPI
endif
#endif

#ifdef PPMPI
    call mpi_dbg_sendrecv_r (rH_x(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,        &
        rH_x(1, 1, 0), ld*ny, MPI_RPREC, down, 1, comm, status, ierr,         &
        'press_rH_x_up')
    call mpi_dbg_sendrecv_r (rH_y(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,        &
        rH_y(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr,         &
        'press_rH_y_up')
    call mpi_dbg_sendrecv_r (rH_z(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,        &
        rH_z(1, 1, 0), ld*ny, MPI_RPREC, down, 3, comm, status, ierr,         &
        'press_rH_z_up')
    call mpi_dbg_sendrecv_r (rH_z(1, 1, 1), ld*ny, MPI_RPREC, down, 6,         &
        rH_z(1, 1, nz), ld*ny, MPI_RPREC, up, 6, comm, status, ierr,          &
        'press_rH_z_down')
#endif

do jz = jz_min, nz
do jy = 1, ny
    if (jy == ny/2 + 1) cycle

    do jx = 1, lh-1

        if (jx*jy == 1) cycle

        ii = 2*jx   ! imaginary index
        ir = ii - 1 ! real index

        !  Compute eye * kx * H_x
        aH_x(1) = -rH_x(ii,jy,jz-1) * kx(jx,jy)
        aH_x(2) =  rH_x(ir,jy,jz-1) * kx(jx,jy)
        aH_y(1) = -rH_y(ii,jy,jz-1) * ky(jx,jy)
        aH_y(2) =  rH_y(ir,jy,jz-1) * ky(jx,jy)

        RHS_col(ir:ii,jy,jz) =  aH_x + aH_y + (rH_z(ir:ii, jy, jz) -           &
            rH_z(ir:ii, jy, jz-1)) *const4

end do
end do
end do

#ifdef ENABLE_CUDA
end if
#endif

#ifdef ENABLE_CUDA
if (press_stage_enabled) then
#ifdef PPMPI
    press_t1 = mpi_wtime()
#else
    call cpu_time(press_t1)
#endif
    press_stage_rhs = press_t1 - press_t0
    press_t0 = press_t1
end if
#endif

! this skips zero wavenumber solution, nyquist freqs
if (press_sync_before_transpose_enabled) then
    press_q_t0 = press_queue_wtime()
    call press_cuda_sync('pressure forced sync before transpose')
    press_q_sync_before_transpose = press_queue_wtime() - press_q_t0
    press_q_sync_count = press_q_sync_count + 1
end if
#ifdef PPMPI
press_halo_req = MPI_REQUEST_NULL
press_p_halo_nreq = 0
call tridag_array (a, b, c, RHS_col, p, press_halo_req, press_p_halo_nreq)
#else
call tridag_array (a, b, c, RHS_col, p)
#endif

#ifdef ENABLE_CUDA
if (press_stage_enabled) then
#ifdef PPMPI
    press_t1 = mpi_wtime()
#else
    call cpu_time(press_t1)
#endif
    press_stage_tridag = press_t1 - press_t0
    press_t0 = press_t1
end if
#endif

! zero-wavenumber solution
#ifdef PPMPI
! wait for p(1, 1, 1) from "down"
if (coord > 0) then
    call mpi_dbg_recv_r (p(1, 1, 1), 2, MPI_RPREC, down, 8, comm, status,     &
        ierr, 'press_p_zero_recv')
end if
#endif

#ifdef ENABLE_CUDA
if (press_cuda_zero_enabled .and. nproc > 1) then
    call press_zero_mode_cuda(rH_z, rbottomw, p)
    call press_cuda_sync('pressure zero mode before send')
else
#endif

if (coord == 0) then
    p(1:2, 1, 0) = 0._rprec
    p(1:2, 1, 1) = p(1:2,1,0) - dz * rbottomw(1:2,1)
end if

do jz = 2, nz
    ! JDA dissertation, eqn(2.88)
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz
end do

#ifdef ENABLE_CUDA
end if
#endif

#ifdef PPMPI
! send p(1, 1, nz) to "up"
if (coord < nproc - 1) then
    call mpi_dbg_send_r (p(1, 1, nz), 2, MPI_RPREC, up, 8, comm, ierr,        &
        'press_p_zero_send')
end if
#endif

#ifdef PPMPI
! make sure 0 <-> nz-1 are syncronized
! 1 <-> nz should be in sync already
if (press_cuda_inverse_enabled .and. press_tridag_gpu_enabled .and.            &
    press_tridag_p_halo_enabled .and. nproc > 1) then
    press_p_halo_pending = press_p_halo_nreq > 0
    press_p_zero_plane_fix_pending = .true.
elseif (press_cuda_inverse_enabled .and. press_nb_p_halo_enabled) then
    press_halo_req = MPI_REQUEST_NULL
    press_p_halo_nreq = 0
    if (coord > 0) then
        press_p_halo_nreq = press_p_halo_nreq + 1
        call mpi_irecv(p(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm,           &
            press_halo_req(press_p_halo_nreq), ierr)
    end if
    if (coord < nproc - 1) then
        press_p_halo_nreq = press_p_halo_nreq + 1
        call mpi_isend(p(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2, comm,          &
            press_halo_req(press_p_halo_nreq), ierr)
    end if
    press_p_halo_pending = press_p_halo_nreq > 0
else
    call mpi_dbg_sendrecv_r (p(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,          &
        p(1, 1, 0), ld*ny, MPI_RPREC, down, 2, comm, status, ierr,            &
        'press_p_plane')
end if
#endif

#ifdef ENABLE_CUDA
if (press_cuda_zero_enabled .and. nproc > 1) then
    if (press_extra_sync_enabled .or. .not. press_cuda_inverse_enabled) then
        call press_cuda_sync('pressure p plane after halo')
    end if
end if
#endif

#ifdef ENABLE_CUDA
if (press_stage_enabled) then
#ifdef PPMPI
    press_t1 = mpi_wtime()
#else
    call cpu_time(press_t1)
#endif
    press_stage_zero = press_t1 - press_t0
    press_t0 = press_t1
end if
#endif

#ifdef ENABLE_CUDA
if (press_cuda_inverse_enabled) then

!$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz-1
    do jy = 1, ny
    do jx = 1, lh
        ii = 2*jx
        ir = ii - 1
        if ((jx == lh) .or. (jy == ny/2+1)) then
            dpdx(ir,jy,jz) = 0._rprec
            dpdx(ii,jy,jz) = 0._rprec
            dpdy(ir,jy,jz) = 0._rprec
            dpdy(ii,jy,jz) = 0._rprec
        else
            dpdx(ir,jy,jz) = -p(ii,jy,jz)*press_kx(jx,jy)
            dpdx(ii,jy,jz) =  p(ir,jy,jz)*press_kx(jx,jy)
            dpdy(ir,jy,jz) = -p(ii,jy,jz)*press_ky(jx,jy)
            dpdy(ii,jy,jz) =  p(ir,jy,jz)*press_ky(jx,jy)
        end if
    end do
    end do
    end do

    press_istat = cufftExecZ2D(press_bk_plan_nz1, dpdx(:,:,1), dpdx(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array dpdx inverse failed: ', press_istat
        stop
    end if
    press_istat = cufftExecZ2D(press_bk_plan_nz1, dpdy(:,:,1), dpdy(:,:,1))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array dpdy inverse failed: ', press_istat
        stop
    end if

#ifdef PPMPI
    if (press_p_halo_pending) then
        call mpi_waitall(press_p_halo_nreq, press_halo_req, press_halo_status,&
            ierr)
        press_p_halo_pending = .false.
    end if
    if (press_p_zero_plane_fix_pending) then
        call mpi_dbg_sendrecv_r(p(1, 1, nz-1), 2, MPI_RPREC, up, 13,         &
            p(1, 1, 0), 2, MPI_RPREC, down, 13, comm, status, ierr,          &
            'press_p_zero_plane')
        press_p_zero_plane_fix_pending = .false.
    end if
#endif

!$cuf kernel do(3) <<<*,*>>>
    do jz = 0, nz
    do jy = 1, ny
    do jx = 1, ld
        if ((jx >= ld-1) .or. (jy == ny/2+1)) p(jx,jy,jz) = 0._rprec
    end do
    end do
    end do

    press_istat = cufftExecZ2D(press_bk_plan_nzp1, p(:,:,0), p(:,:,0))
    if (press_istat /= CUFFT_SUCCESS) then
        print *, 'press_stag_array pressure inverse failed: ', press_istat
        stop
    end if

!$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, nx
        if (jz < nz) then
            dpdz(jx,jy,jz) = (p(jx,jy,jz) - p(jx,jy,jz-1))/dz
        else if (coord == nproc-1) then
            dpdz(jx,jy,jz) = (p(jx,jy,nz) - p(jx,jy,nz-1))/dz
        else
            dpdz(jx,jy,jz) = 0._rprec
        end if
    end do
    end do
    end do
    call press_cuda_sync('pressure hybrid inverse and dpdz')

else
#endif

    ! zero the nyquist freqs
    p(ld-1:ld,:,:) = 0._rprec
    p(:,ny/2+1,:) = 0._rprec

    ! Now need to get p(wave,level) to physical p(jx,jy,jz)
    ! Loop over height levels
    call dfftw_execute_dft_c2r(back,p(:,:,0), p(:,:,0))
    do jz = 1, nz-1
        do jy = 1, ny
        do jx = 1,lh
            ii = 2*jx
            ir = ii - 1
            dpdx(ir,jy,jz) = -p(ii,jy,jz) * kx(jx,jy)
            dpdx(ii,jy,jz) =  p(ir,jy,jz) * kx(jx,jy)
            dpdy(ir,jy,jz) = -p(ii,jy,jz) * ky(jx,jy)
            dpdy(ii,jy,jz) =  p(ir,jy,jz) * ky(jx,jy)
        end do
        end do

        ! The inactive Fourier modes of p are already zero before this inverse
        ! transform, so dpdx and dpdy inherit valid padding.
        call dfftw_execute_dft_c2r(back,dpdx(:,:,jz), dpdx(:,:,jz))
        call dfftw_execute_dft_c2r(back,dpdy(:,:,jz), dpdy(:,:,jz))
        call dfftw_execute_dft_c2r(back,p(:,:,jz), p(:,:,jz))
    end do

    if(coord==nproc-1) call dfftw_execute_dft_c2r(back,p(:,:,nz),p(:,:,nz))

    ! nz level is not needed elsewhere (although its valid)
#ifdef PPSAFETYMODE
    dpdx(:,:,nz) = BOGUS
    dpdy(:,:,nz) = BOGUS
    if(coord<nproc-1) p(:,:,nz) = BOGUS
#endif

    ! Final step compute the z-derivative of p
    ! note: p has additional level at z=-dz/2 for this derivative
    dpdz(1:nx, 1:ny, 1:nz-1) = (p(1:nx, 1:ny, 1:nz-1) - p(1:nx, 1:ny, 0:nz-2)) / dz
#ifdef PPSAFETYMODE
    if(coord<nproc-1)  dpdz(:,:,nz) = BOGUS
#endif
    if(coord==nproc-1) dpdz(1:nx,1:ny,nz) = (p(1:nx,1:ny,nz)-p(1:nx,1:ny,nz-1))/ dz

#ifdef ENABLE_CUDA
end if
#endif

#ifdef ENABLE_CUDA
if (press_stage_enabled) then
#ifdef PPMPI
    press_t1 = mpi_wtime()
#else
    call cpu_time(press_t1)
#endif
    press_stage_inverse = press_t1 - press_t0
    call press_stage_report(press_stage_count, press_stage_forward,            &
        press_stage_rhs, press_stage_tridag, press_stage_zero,                 &
        press_stage_inverse)
end if

if (press_queue_attrib_enabled) then
    call press_queue_report(press_stage_count, press_q_pack_cpu,               &
        press_q_pack_gpu, press_q_fft_cpu, press_q_fft_gpu,                   &
        press_q_rhs_prep_cpu, press_q_rhs_prep_gpu, press_q_rhs_halo_cpu,     &
        press_q_rhs_assembly_cpu, press_q_rhs_assembly_gpu,                   &
        press_q_sync_after_forward, press_q_sync_after_rhs,                   &
        press_q_sync_before_transpose, press_q_sync_count,                    &
        press_queue_header_printed)
end if

end if
#endif

end subroutine press_stag_array

!*******************************************************************************
subroutine press_apply_env_enabled_unless_false(name, enabled)
!*******************************************************************************
!
! Preserve the pressure-path switch convention: an unset variable keeps the
! current default, explicit false tokens disable the path, and any other set
! value enables it.
!
implicit none

character(len=*), intent(in) :: name
logical, intent(inout) :: enabled
character(len=32) :: setting
integer :: stat

call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('0', 'false', 'FALSE', 'False', 'off', 'OFF', 'Off', 'no', 'NO',    &
        'No')
        enabled = .false.
    case default
        enabled = .true.
    end select
end if

end subroutine press_apply_env_enabled_unless_false

#ifdef ENABLE_CUDA
!*******************************************************************************
subroutine press_cuda_sync(where)
!*******************************************************************************
use cudafor
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'press_stag_array CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'press_stag_array CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine press_cuda_sync

!*******************************************************************************
real(rprec) function press_queue_wtime()
!*******************************************************************************
use types, only : rprec
#ifdef PPMPI
use mpi
#endif
implicit none

#ifdef PPMPI
press_queue_wtime = mpi_wtime()
#else
call cpu_time(press_queue_wtime)
#endif

end function press_queue_wtime

!*******************************************************************************
subroutine press_queue_event_start(evt, wall_start)
!*******************************************************************************
use types, only : rprec
use cudafor
implicit none

type(cudaEvent), intent(inout) :: evt
real(rprec), intent(out) :: wall_start
real(rprec) :: press_queue_wtime
integer :: istat

wall_start = press_queue_wtime()
istat = cudaEventRecord(evt, 0)
if (istat /= cudaSuccess) stop 'press queue event record start failed'

end subroutine press_queue_event_start

!*******************************************************************************
subroutine press_queue_event_stop(evt_start, evt_stop, wall_start, cpu_wall,    &
    gpu_event, sync_count, where)
!*******************************************************************************
use types, only : rprec
use cudafor
implicit none

type(cudaEvent), intent(inout) :: evt_start, evt_stop
real(rprec), intent(in) :: wall_start
real(rprec), intent(out) :: cpu_wall, gpu_event
integer, intent(inout) :: sync_count
character(len=*), intent(in) :: where
real(rprec) :: press_queue_wtime
integer :: istat
real :: elapsed_ms

istat = cudaEventRecord(evt_stop, 0)
if (istat /= cudaSuccess) then
    print *, 'press queue event record stop failed at ', trim(where), ': ',    &
        istat
    stop
end if
istat = cudaEventSynchronize(evt_stop)
if (istat /= cudaSuccess) then
    print *, 'press queue event sync failed at ', trim(where), ': ', istat
    stop
end if
sync_count = sync_count + 1
istat = cudaEventElapsedTime(elapsed_ms, evt_start, evt_stop)
if (istat /= cudaSuccess) then
    print *, 'press queue event elapsed failed at ', trim(where), ': ', istat
    stop
end if
gpu_event = real(elapsed_ms, rprec) * 1.0e-3_rprec
cpu_wall = press_queue_wtime() - wall_start

end subroutine press_queue_event_stop

!*******************************************************************************
subroutine press_queue_report(call_count, pack_cpu, pack_gpu, fft_cpu, fft_gpu,&
    rhs_prep_cpu, rhs_prep_gpu, rhs_halo_cpu, rhs_assembly_cpu,                &
    rhs_assembly_gpu, sync_after_forward, sync_after_rhs,                      &
    sync_before_transpose, sync_count, header_printed)
!*******************************************************************************
use types, only : rprec
use param, only : coord
#ifdef PPMPI
use param, only : nproc, comm, ierr
use mpi
#endif
implicit none

integer, intent(in) :: call_count, sync_count
real(rprec), intent(in) :: pack_cpu, pack_gpu, fft_cpu, fft_gpu
real(rprec), intent(in) :: rhs_prep_cpu, rhs_prep_gpu, rhs_halo_cpu
real(rprec), intent(in) :: rhs_assembly_cpu, rhs_assembly_gpu
real(rprec), intent(in) :: sync_after_forward, sync_after_rhs
real(rprec), intent(in) :: sync_before_transpose
logical, intent(inout) :: header_printed
integer :: rr

if (coord == 0 .and. .not. header_printed) then
    write(*,'(A)') 'PRESS_QUEUE_ATTRIB fields: call rank forward_pack_cpu ' // &
        'forward_pack_gpu forward_fft_cpu forward_fft_gpu rhs_prep_cpu ' //    &
        'rhs_prep_gpu rhs_halo_cpu rhs_assembly_cpu rhs_assembly_gpu ' //      &
        'sync_after_forward sync_after_rhs sync_before_transpose sync_count'
    header_printed = .true.
end if

#ifdef PPMPI
do rr = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == rr) then
#endif
        write(*,'(A,I6,1X,I4,12(1X,ES12.5),1X,I4)')                           &
            'PRESS_QUEUE_ATTRIB ', call_count, coord, pack_cpu, pack_gpu,      &
            fft_cpu, fft_gpu, rhs_prep_cpu, rhs_prep_gpu, rhs_halo_cpu,        &
            rhs_assembly_cpu, rhs_assembly_gpu, sync_after_forward,            &
            sync_after_rhs, sync_before_transpose, sync_count
        flush(6)
#ifdef PPMPI
    end if
end do
call mpi_barrier(comm, ierr)
#endif

end subroutine press_queue_report

!*******************************************************************************
subroutine press_rhs_halo_audit(src_file, routine, path, mpi_pattern, messages,&
    send_bytes, recv_bytes, neighbor, send_type, recv_type, contiguous,         &
    direct_section, gpu_id, cuda_visible, mpich_gpu)
!*******************************************************************************
use param, only : coord
#ifdef PPMPI
use param, only : nproc, comm, ierr
use mpi
#endif
implicit none

character(len=*), intent(in) :: src_file, routine, path, mpi_pattern
character(len=*), intent(in) :: send_type, recv_type, contiguous
character(len=*), intent(in) :: direct_section, cuda_visible, mpich_gpu
integer, intent(in) :: messages, send_bytes, recv_bytes, neighbor, gpu_id
integer :: rr

#ifdef PPMPI
do rr = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == rr) then
#endif
        write(*,*) 'PRESS_RHS_HALO_AUDIT rank=', coord, ' file=',             &
            trim(src_file), ' routine=', trim(routine), ' path=', trim(path), &
            ' mpi=', trim(mpi_pattern), ' messages=', messages,               &
            ' send_bytes=', send_bytes, ' recv_bytes=', recv_bytes,           &
            ' neighbor=', neighbor, ' send_ptr=', trim(send_type),            &
            ' recv_ptr=', trim(recv_type), ' contiguous=', trim(contiguous),  &
            ' direct_section=', trim(direct_section), ' gpu=', gpu_id,        &
            ' CUDA_VISIBLE_DEVICES=', trim(cuda_visible),                     &
            ' MPICH_GPU_SUPPORT_ENABLED=', trim(mpich_gpu)
        flush(6)
#ifdef PPMPI
    end if
end do
call mpi_barrier(comm, ierr)
#endif

end subroutine press_rhs_halo_audit

!*******************************************************************************
subroutine press_rhs_halo_report(call_count, path, pre_sync, pack_time,         &
    mpi_time, unpack_time, post_sync, total_time, send_bytes, recv_bytes,       &
    messages, ready_rel, arrival_skew, barrier_time)
!*******************************************************************************
use types, only : rprec
use param, only : coord
#ifdef PPMPI
use param, only : nproc, comm, ierr
use mpi
#endif
implicit none

integer, intent(in) :: call_count, send_bytes, recv_bytes, messages
real(rprec), intent(in) :: pre_sync, pack_time, mpi_time, unpack_time
real(rprec), intent(in) :: post_sync, total_time, ready_rel, arrival_skew
real(rprec), intent(in) :: barrier_time
character(len=*), intent(in) :: path
logical, save :: header_printed = .false.
integer :: rr
real(rprec) :: eff_gbps

if (mpi_time > 0._rprec) then
    eff_gbps = real(send_bytes + recv_bytes, rprec) / mpi_time / 1.0e9_rprec
else
    eff_gbps = 0._rprec
end if

if (coord == 0 .and. .not. header_printed) then
    write(*,'(A)') 'PRESS_RHS_HALO_TIMING fields: call rank path pre_sync ' //&
        'pack mpi_wait unpack post_sync total send_bytes recv_bytes messages ' &
        // 'ready_rel arrival_skew barrier effective_GBps'
    header_printed = .true.
end if

#ifdef PPMPI
do rr = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == rr) then
#endif
        write(*,*) 'PRESS_RHS_HALO_TIMING ', call_count, coord, trim(path),   &
            pre_sync, pack_time, mpi_time, unpack_time, post_sync,            &
            total_time, send_bytes, recv_bytes, messages, ready_rel,          &
            arrival_skew, barrier_time, eff_gbps
        flush(6)
#ifdef PPMPI
    end if
end do
call mpi_barrier(comm, ierr)
#endif

end subroutine press_rhs_halo_report

!*******************************************************************************
subroutine press_rhs_assembly_report(call_count, interior_cpu, interior_gpu,    &
    boundary_cpu, boundary_gpu, total_cpu, total_gpu)
!*******************************************************************************
use types, only : rprec
use param, only : coord
#ifdef PPMPI
use param, only : nproc, comm, ierr
use mpi
#endif
implicit none

integer, intent(in) :: call_count
real(rprec), intent(in) :: interior_cpu, interior_gpu, boundary_cpu
real(rprec), intent(in) :: boundary_gpu, total_cpu, total_gpu
logical, save :: header_printed = .false.
integer :: rr

if (coord == 0 .and. .not. header_printed) then
    write(*,'(A)') 'PRESS_RHS_ASSEMBLY_SPLIT fields: call rank interior_cpu ' &
        // 'interior_gpu boundary_cpu boundary_gpu total_cpu total_gpu'
    header_printed = .true.
end if

#ifdef PPMPI
do rr = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == rr) then
#endif
        write(*,*) 'PRESS_RHS_ASSEMBLY_SPLIT ', call_count, coord,            &
            interior_cpu, interior_gpu, boundary_cpu, boundary_gpu,           &
            total_cpu, total_gpu
        flush(6)
#ifdef PPMPI
    end if
end do
call mpi_barrier(comm, ierr)
#endif

end subroutine press_rhs_assembly_report

!*******************************************************************************
subroutine press_stage_report(stage_count, forward, rhs, tridag, zero_mode,    &
    inverse)
!*******************************************************************************
use types, only : rprec
use param, only : coord
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC
use mpi
#endif
implicit none

integer, intent(in) :: stage_count
real(rprec), intent(in) :: forward, rhs, tridag, zero_mode, inverse
real(rprec) :: fmax, rmax, tmax, zmax, imax, total

#ifdef PPMPI
call mpi_allreduce(forward, fmax, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(rhs, rmax, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(tridag, tmax, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(zero_mode, zmax, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(inverse, imax, 1, MPI_RPREC, MPI_MAX, comm, ierr)
#else
fmax = forward
rmax = rhs
tmax = tridag
zmax = zero_mode
imax = inverse
#endif

total = fmax + rmax + tmax + zmax + imax

if (coord == 0) then
    write(*,'(a,i8)') 'Pressure solver stage timing (max rank), call ',        &
        stage_count
    write(*,'(1a,E15.7)') '  Forward pack/FFT/sync: ', fmax
    write(*,'(1a,E15.7)') '  RHS assembly + pressure halos: ', rmax
    write(*,'(1a,E15.7)') '  Tridiagonal solve: ', tmax
    write(*,'(1a,E15.7)') '  Zero mode + p halo: ', zmax
    write(*,'(1a,E15.7)') '  Inverse FFT + dpdz/sync: ', imax
    write(*,'(1a,E15.7)') '  Stage sum: ', total
end if

end subroutine press_stage_report

!*******************************************************************************
subroutine press_pack_rhs_cuda(u_in, v_in, w_in, rhx, rhy, rhz, scale)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz
use cudafor
implicit none

real(rprec), intent(in) :: scale
real(rprec), managed, intent(in) :: u_in(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: v_in(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: w_in(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhz(ld,ny,lbz:nz)
integer :: jx, jy, jz

!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz - 1
do jy = 1, ny
do jx = 1, ld
    rhx(jx,jy,jz) = scale * u_in(jx,jy,jz)
    rhy(jx,jy,jz) = scale * v_in(jx,jy,jz)
    rhz(jx,jy,jz) = scale * w_in(jx,jy,jz)
end do
end do
end do

end subroutine press_pack_rhs_cuda

!*******************************************************************************
subroutine press_pack_rhs_halo_cuda(rhx, rhy, rhz, halo_buf)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz
use cudafor
implicit none

real(rprec), managed, intent(in) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhz(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: halo_buf(3*ld*ny)
integer :: idx, jx, jy, plane_size

plane_size = ld * ny

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, plane_size
    jx = mod(idx - 1, ld) + 1
    jy = (idx - 1) / ld + 1
    halo_buf(idx) = rhx(jx,jy,nz-1)
    halo_buf(plane_size + idx) = rhy(jx,jy,nz-1)
    halo_buf(2*plane_size + idx) = rhz(jx,jy,nz-1)
end do

end subroutine press_pack_rhs_halo_cuda

!*******************************************************************************
subroutine press_unpack_rhs_halo_cuda(halo_buf, rhx, rhy, rhz)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz
use cudafor
implicit none

real(rprec), managed, intent(in) :: halo_buf(3*ld*ny)
real(rprec), managed, intent(inout) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhz(ld,ny,lbz:nz)
integer :: idx, jx, jy, plane_size

plane_size = ld * ny

!$cuf kernel do(1) <<<*,*>>>
do idx = 1, plane_size
    jx = mod(idx - 1, ld) + 1
    jy = (idx - 1) / ld + 1
    rhx(jx,jy,0) = halo_buf(idx)
    rhy(jx,jy,0) = halo_buf(plane_size + idx)
    rhz(jx,jy,0) = halo_buf(2*plane_size + idx)
end do

end subroutine press_unpack_rhs_halo_cuda

!*******************************************************************************
subroutine press_pack_rhs_halo_combined_cuda(rhx, rhy, rhz, halo_buf, coord_in)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz
use cudafor
implicit none

real(rprec), managed, intent(in) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhz(ld,ny,lbz:nz)
real(rprec), device, intent(inout) :: halo_buf(3*ld*ny)
integer, intent(in) :: coord_in
integer :: idx, jx, jy, plane_size

plane_size = ld * ny

if (coord_in == 0) then
!$cuf kernel do(1) <<<*,*>>>
    do idx = 1, plane_size
        jx = mod(idx - 1, ld) + 1
        jy = (idx - 1) / ld + 1
        halo_buf(idx) = rhx(jx,jy,nz-1)
        halo_buf(plane_size + idx) = rhy(jx,jy,nz-1)
        halo_buf(2*plane_size + idx) = rhz(jx,jy,nz-1)
    end do
else
!$cuf kernel do(1) <<<*,*>>>
    do idx = 1, plane_size
        jx = mod(idx - 1, ld) + 1
        jy = (idx - 1) / ld + 1
        halo_buf(idx) = 0._rprec
        halo_buf(plane_size + idx) = 0._rprec
        halo_buf(2*plane_size + idx) = rhz(jx,jy,1)
    end do
end if

end subroutine press_pack_rhs_halo_combined_cuda

!*******************************************************************************
subroutine press_unpack_rhs_halo_combined_cuda(halo_buf, rhx, rhy, rhz, coord_in)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz
use cudafor
implicit none

real(rprec), device, intent(in) :: halo_buf(3*ld*ny)
real(rprec), managed, intent(inout) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhz(ld,ny,lbz:nz)
integer, intent(in) :: coord_in
integer :: idx, jx, jy, plane_size

plane_size = ld * ny

if (coord_in == 0) then
!$cuf kernel do(1) <<<*,*>>>
    do idx = 1, plane_size
        jx = mod(idx - 1, ld) + 1
        jy = (idx - 1) / ld + 1
        rhz(jx,jy,nz) = halo_buf(2*plane_size + idx)
    end do
else
!$cuf kernel do(1) <<<*,*>>>
    do idx = 1, plane_size
        jx = mod(idx - 1, ld) + 1
        jy = (idx - 1) / ld + 1
        rhx(jx,jy,0) = halo_buf(idx)
        rhy(jx,jy,0) = halo_buf(plane_size + idx)
        rhz(jx,jy,0) = halo_buf(2*plane_size + idx)
    end do
end if

end subroutine press_unpack_rhs_halo_combined_cuda

!*******************************************************************************
subroutine press_rhs_prep_cuda(rhx, rhy, rhz, rtop, rbottom, rhs)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz, coord, nproc, dz
use cudafor
implicit none

real(rprec), managed, intent(inout) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhz(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rtop(ld,ny), rbottom(ld,ny)
real(rprec), managed, intent(inout) :: rhs(ld,ny,nz+1)
integer :: jx, jy, jz

!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    if ((jz <= nz-1) .and. ((jx >= ld-1) .or. (jy == ny/2+1))) then
        rhx(jx,jy,jz) = 0._rprec
        rhy(jx,jy,jz) = 0._rprec
        rhz(jx,jy,jz) = 0._rprec
    end if
    if ((coord == nproc-1) .and. (jz == nz) .and.                         &
        ((jx >= ld-1) .or. (jy == ny/2+1))) then
        rhz(jx,jy,nz) = 0._rprec
    end if
end do
end do
end do

!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, ld
    if ((jx >= ld-1) .or. (jy == ny/2+1)) then
        rtop(jx,jy) = 0._rprec
        rbottom(jx,jy) = 0._rprec
    end if
    if (coord == 0) rhs(jx,jy,1) = -dz * rbottom(jx,jy)
    if (coord == nproc-1) rhs(jx,jy,nz+1) = -dz * rtop(jx,jy)
end do
end do

end subroutine press_rhs_prep_cuda

!*******************************************************************************
subroutine press_assemble_rhs_cuda(rhx, rhy, rhz, rhs, kx_d, ky_d, dz_inv,     &
    jz_min)
!*******************************************************************************
use types, only : rprec
use param, only : ld, lh, ny, nz, lbz
implicit none

real(rprec), intent(in) :: dz_inv
integer, intent(in) :: jz_min
real(rprec), managed, intent(in) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhz(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhs(ld,ny,nz+1)
real(rprec), managed, intent(in) :: kx_d(lh,ny), ky_d(lh,ny)
call press_assemble_rhs_range_cuda(rhx, rhy, rhz, rhs, kx_d, ky_d, dz_inv,    &
    jz_min, nz)

end subroutine press_assemble_rhs_cuda

!*******************************************************************************
subroutine press_assemble_rhs_range_cuda(rhx, rhy, rhz, rhs, kx_d, ky_d,       &
    dz_inv, jz_first, jz_last)
!*******************************************************************************
use types, only : rprec
use param, only : ld, lh, ny, nz, lbz
use cudafor
implicit none

real(rprec), intent(in) :: dz_inv
integer, intent(in) :: jz_first, jz_last
real(rprec), managed, intent(in) :: rhx(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhy(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rhz(ld,ny,lbz:nz)
real(rprec), managed, intent(inout) :: rhs(ld,ny,nz+1)
real(rprec), managed, intent(in) :: kx_d(lh,ny), ky_d(lh,ny)
integer :: jx, jy, jz, ir, ii

!$cuf kernel do(3) <<<*,*>>>
do jz = jz_first, jz_last
do jy = 1, ny
do jx = 1, lh-1
    if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
        ii = 2*jx
        ir = ii - 1
        rhs(ir,jy,jz) = -rhx(ii,jy,jz-1) * kx_d(jx,jy)                       &
            - rhy(ii,jy,jz-1) * ky_d(jx,jy)                                  &
            + (rhz(ir,jy,jz) - rhz(ir,jy,jz-1)) * dz_inv
        rhs(ii,jy,jz) =  rhx(ir,jy,jz-1) * kx_d(jx,jy)                       &
            + rhy(ir,jy,jz-1) * ky_d(jx,jy)                                  &
            + (rhz(ii,jy,jz) - rhz(ii,jy,jz-1)) * dz_inv
    end if
end do
end do
end do

end subroutine press_assemble_rhs_range_cuda

!*******************************************************************************
subroutine press_zero_mode_cuda(rhz, rbottom, p_out)
!*******************************************************************************
use types, only : rprec
use param, only : ld, ny, nz, lbz, coord, dz
use cudafor
implicit none

real(rprec), managed, intent(in) :: rhz(ld,ny,lbz:nz)
real(rprec), managed, intent(in) :: rbottom(ld,ny)
real(rprec), managed, intent(inout) :: p_out(ld,ny,0:nz)
integer :: jx, jz

!$cuf kernel do(1) <<<*,*>>>
do jx = 1, 2
    if (coord == 0) then
        p_out(jx,1,0) = 0._rprec
        p_out(jx,1,1) = -dz * rbottom(jx,1)
    end if
    do jz = 2, nz
        p_out(jx,1,jz) = p_out(jx,1,jz-1) + rhz(jx,1,jz) * dz
    end do
end do

end subroutine press_zero_mode_cuda
#endif
#endif
