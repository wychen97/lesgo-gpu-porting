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

!*******************************************************************************
module sgs_stag_util
!*******************************************************************************
! Navigation map for the SGS module:
!   - module state: timing counters, diagnostics, and GPU event handles
!   - runtime SGS entry point: sgs_stag
!   - strain-rate assembly: calc_Sij and GPU helper kernels
!   - dynamic SGS coefficients: calc_Sij_nut_dynamic_cuda and lagrange modules
!   - halo paths and GPU-aware MPI: tau/dwdz sync helpers and path audits
!   - diagnostics: *_detail_report and sgs_stage_report
!
! Keep runtime SGS model branches aligned.  If a new optimized path is added
! for one `sgs` value, check the other supported runtime SGS values before
! declaring the module complete.
use types, only : rprec
#ifdef ENABLE_CUDA
use cudafor
#endif
implicit none

save
private

public sgs_stag, rtnewt
#ifdef ENABLE_CUDA
public sgs_combined_tau_halo_active
#endif

#ifdef ENABLE_CUDA
! SGS GPU ownership map:
!   - sim_param owns the resident velocity, stress, and div-stress fields used
!     by the timestep hot path.
!   - this module owns SGS dispatch policy, timing/audit counters, and
!     tau/dwdz halo staging buffers.
!   - lagrange_Sdep_gpu.f90 owns the batched Lagrangian dynamic SGS update;
!     keep runtime sgs values 1 through 5 aligned when routing into it.
!   - LESGO_SGS_* switches remain documented diagnostics/path selectors, not
!     hidden production defaults.
integer, parameter :: SGS_CALC_BOTTOM = 1, SGS_CALC_TOP = 2
integer, parameter :: SGS_CALC_INTERIOR = 3, SGS_CALC_HALO_PLANE = 4
integer, parameter :: SGS_CALC_DWDZ_PRE_SYNC = 5
integer, parameter :: SGS_CALC_DWDZ_MPI_WAIT = 6
integer, parameter :: SGS_CALC_COUNT = 6
integer, parameter :: SGS_TAU_PRE_SYNC = 1, SGS_TAU_PACK = 2
integer, parameter :: SGS_TAU_MPI = 3, SGS_TAU_UNPACK = 4
integer, parameter :: SGS_TAU_POST_SYNC = 5, SGS_TAU_D2H = 6
integer, parameter :: SGS_TAU_H2D = 7, SGS_TAU_COUNT = 7
integer, parameter :: SGS_DWDZ_PRE_SYNC = 1, SGS_DWDZ_PACK = 2
integer, parameter :: SGS_DWDZ_MPI = 3, SGS_DWDZ_UNPACK = 4
integer, parameter :: SGS_DWDZ_POST_SYNC = 5, SGS_DWDZ_D2H = 6
integer, parameter :: SGS_DWDZ_H2D = 7, SGS_DWDZ_DUMMY = 8
integer, parameter :: SGS_DWDZ_COUNT = 8
integer, parameter :: SGS_RPREC_BYTES = 8
type(cudaEvent), save :: sgs_calc_evt_start(SGS_CALC_COUNT)
type(cudaEvent), save :: sgs_calc_evt_stop(SGS_CALC_COUNT)
logical, save :: sgs_calc_evt_created = .false.
logical, save :: sgs_calc_evt_active(SGS_CALC_COUNT) = .false.
logical, save :: sgs_calc_audit_printed = .false.
real(rprec), save :: sgs_calc_cpu(SGS_CALC_COUNT) = 0._rprec
real(rprec), save :: sgs_calc_gpu(SGS_CALC_COUNT) = 0._rprec
real(rprec), save :: sgs_calc_t0(SGS_CALC_COUNT) = 0._rprec
integer, save :: sgs_calc_launches(SGS_CALC_COUNT) = 0
integer, save :: sgs_calc_total_launches = 0
integer(8), save :: sgs_calc_cells(SGS_CALC_COUNT) = 0_8
integer(8), save :: sgs_calc_total_cells = 0_8
integer, save :: sgs_calc_jz_min = 0, sgs_calc_jz_max = 0
integer, save :: sgs_calc_bulk_jz_min = 0, sgs_calc_bulk_jz_max = 0
integer, save :: sgs_calc_lbc_mode = -999, sgs_calc_ubc_mode = -999
integer, save :: sgs_calc_block_x = 0, sgs_calc_block_y = 0, sgs_calc_block_z = 0
integer, save :: sgs_calc_grid_x = 0, sgs_calc_grid_y = 0, sgs_calc_grid_z = 0
logical, save :: sgs_calc_bottom_branch = .false.
logical, save :: sgs_calc_top_branch = .false.
logical, save :: sgs_calc_explicit_used = .false.
real(rprec), save :: sgs_tau_detail(SGS_TAU_COUNT) = 0._rprec
real(rprec), save :: sgs_tau_t0 = 0._rprec
integer(8), save :: sgs_tau_send_bytes = 0_8, sgs_tau_recv_bytes = 0_8
integer, save :: sgs_tau_mpi_calls = 0, sgs_tau_neighbor = -999
logical, save :: sgs_tau_combined_msg = .false.
real(rprec), save :: sgs_tau_arr_enter = 0._rprec
real(rprec), save :: sgs_tau_arr_after_sync = 0._rprec
real(rprec), save :: sgs_tau_arr_before_mpi = 0._rprec
real(rprec), save :: sgs_tau_arr_after_mpi = 0._rprec
real(rprec), save :: sgs_tau_arr_after_post = 0._rprec
real(rprec), save :: sgs_tau_arr_barrier = 0._rprec
real(rprec), save :: sgs_dwdz_detail(SGS_DWDZ_COUNT) = 0._rprec
real(rprec), save :: sgs_dwdz_t0 = 0._rprec
integer(8), save :: sgs_dwdz_send_bytes = 0_8, sgs_dwdz_recv_bytes = 0_8
integer(8), save :: sgs_dwdz_dummy_bytes = 0_8
integer, save :: sgs_dwdz_mpi_calls = 0, sgs_dwdz_neighbor = -999
logical, save :: sgs_dwdz_combined_msg = .false.
real(rprec), save :: sgs_dwdz_arr_enter = 0._rprec
real(rprec), save :: sgs_dwdz_arr_after_sync = 0._rprec
real(rprec), save :: sgs_dwdz_arr_before_mpi = 0._rprec
real(rprec), save :: sgs_dwdz_arr_after_mpi = 0._rprec
real(rprec), save :: sgs_dwdz_arr_after_post = 0._rprec
real(rprec), save :: sgs_dwdz_arr_barrier = 0._rprec
real(rprec), device, allocatable, save :: sgs_tau2_send_d(:)
real(rprec), device, allocatable, save :: sgs_tau2_recv_d(:)
real(rprec), pinned, allocatable, save :: sgs_tau2_send_h(:)
real(rprec), pinned, allocatable, save :: sgs_tau2_recv_h(:)
real(rprec), device, allocatable, save :: sgs_dwdz_send_d(:)
real(rprec), device, allocatable, save :: sgs_dwdz_recv_d(:)
real(rprec), device, allocatable, save :: sgs_dwdz_dummy_send_d(:)
real(rprec), device, allocatable, save :: sgs_dwdz_dummy_recv_d(:)
real(rprec), pinned, allocatable, save :: sgs_dwdz_send_h(:)
real(rprec), pinned, allocatable, save :: sgs_dwdz_recv_h(:)
logical, save :: sgs_dwdz_ptr_audit_printed = .false.
logical, save :: sgs_tau_ptr_audit_printed = .false.
logical, save :: sgs_dwdz_path_audit_printed = .false.
#endif

contains

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function sgs_env_true_token_enabled(name)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
character(len=16) :: setting
integer :: stat

sgs_env_true_token_enabled = .false.
call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('1', 'true', 'TRUE', 'True', 'on', 'ON', 'On', 'yes', 'YES',        &
        'Yes')
        sgs_env_true_token_enabled = .true.
    case default
        sgs_env_true_token_enabled = .false.
    end select
end if

end function sgs_env_true_token_enabled

!*******************************************************************************
logical function sgs_env_enabled_unless_false(name, default_value)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
logical, intent(in) :: default_value
character(len=16) :: setting
integer :: stat

sgs_env_enabled_unless_false = default_value
call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('0', 'false', 'FALSE', 'False', 'off', 'OFF', 'Off', 'no', 'NO',    &
        'No')
        sgs_env_enabled_unless_false = .false.
    case default
        sgs_env_enabled_unless_false = .true.
    end select
end if

end function sgs_env_enabled_unless_false

!*******************************************************************************
logical function sgs_pointwise_cuda_enabled()
!*******************************************************************************
implicit none

sgs_pointwise_cuda_enabled = .true.

end function sgs_pointwise_cuda_enabled

!*******************************************************************************
logical function sgs_extra_sync_enabled()
!*******************************************************************************
implicit none

sgs_extra_sync_enabled = .false.

end function sgs_extra_sync_enabled

!*******************************************************************************
logical function sgs_packed_tau_halo_enabled()
!*******************************************************************************
implicit none

sgs_packed_tau_halo_enabled = .true.

end function sgs_packed_tau_halo_enabled

!*******************************************************************************
logical function sgs_combined_tau_halo_enabled()
!*******************************************************************************
implicit none

sgs_combined_tau_halo_enabled = .true.

end function sgs_combined_tau_halo_enabled

!*******************************************************************************
logical function sgs_combined_tau_halo_active()
!*******************************************************************************
implicit none

sgs_combined_tau_halo_active = sgs_pointwise_cuda_enabled() .and.              &
    sgs_combined_tau_halo_enabled()

end function sgs_combined_tau_halo_active

!*******************************************************************************
logical function sgs_halo_combined2_enabled()
!*******************************************************************************
use param, only : coord
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.
logical, save :: printed = .false.

if (.not. initialized) then
    enabled = sgs_env_enabled_unless_false('LESGO_SGS_HALO_COMBINED', .false.)
    initialized = .true.
end if
if (coord == 0 .and. .not. printed) then
    write(*,*) 'SGS tau halo nproc2 combined device path enabled=', enabled
    printed = .true.
end if

sgs_halo_combined2_enabled = enabled

end function sgs_halo_combined2_enabled

!*******************************************************************************
logical function sgs_halo_hostpinned_test_enabled()
!*******************************************************************************
implicit none

sgs_halo_hostpinned_test_enabled = .false.

end function sgs_halo_hostpinned_test_enabled

!*******************************************************************************
logical function sgs_barrier_before_tau_halo_enabled()
!*******************************************************************************
implicit none

sgs_barrier_before_tau_halo_enabled = .false.

end function sgs_barrier_before_tau_halo_enabled

!*******************************************************************************
logical function sgs_fused_dynamic_nut_enabled()
!*******************************************************************************
!
! Experimental dynamic-model Sij/Nu_t fusion was validated but slower on the
! current ATM Derecho A100 case; keep it compiled out of the production path.
!
implicit none

sgs_fused_dynamic_nut_enabled = .false.

end function sgs_fused_dynamic_nut_enabled

!*******************************************************************************
logical function sgs_stage_timing_enabled()
!*******************************************************************************
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = sgs_env_true_token_enabled('LESGO_SGS_STAGE_TIMING')
    initialized = .true.
end if

sgs_stage_timing_enabled = enabled

end function sgs_stage_timing_enabled

!*******************************************************************************
logical function sgs_strict_sync_enabled()
!*******************************************************************************
!
! Strict mode is a diagnostic fallback that synchronizes after every SGS CUDA
! kernel check. Non-strict mode keeps only required MPI/final barriers so wall
! timers do not absorb unrelated queued GPU work.
!
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = sgs_env_true_token_enabled('LESGO_SGS_STRICT_SYNC')
    initialized = .true.
end if

sgs_strict_sync_enabled = enabled

end function sgs_strict_sync_enabled

!*******************************************************************************
logical function sgs_direct_dwdz_halo_enabled()
!*******************************************************************************
implicit none

sgs_direct_dwdz_halo_enabled = .false.

end function sgs_direct_dwdz_halo_enabled

!*******************************************************************************
logical function sgs_dwdz_device_halo_enabled()
!*******************************************************************************
implicit none

sgs_dwdz_device_halo_enabled = .false.

end function sgs_dwdz_device_halo_enabled

!*******************************************************************************
logical function sgs_dwdz_hostpinned_test_enabled()
!*******************************************************************************
implicit none

sgs_dwdz_hostpinned_test_enabled = .false.

end function sgs_dwdz_hostpinned_test_enabled

!*******************************************************************************
logical function sgs_dwdz_symmetric_test_enabled()
!*******************************************************************************
implicit none

sgs_dwdz_symmetric_test_enabled = .false.

end function sgs_dwdz_symmetric_test_enabled

!*******************************************************************************
logical function sgs_barrier_before_dwdz_enabled()
!*******************************************************************************
implicit none

sgs_barrier_before_dwdz_enabled = .false.

end function sgs_barrier_before_dwdz_enabled

!*******************************************************************************
logical function sgs_tau_prebarrier_enabled()
!*******************************************************************************
implicit none

sgs_tau_prebarrier_enabled = .false.

end function sgs_tau_prebarrier_enabled

!*******************************************************************************
logical function sgs_overlap_tau_halo_enabled()
!*******************************************************************************
implicit none

sgs_overlap_tau_halo_enabled = .true.

end function sgs_overlap_tau_halo_enabled

!*******************************************************************************
logical function sgs_overlap_dwdz_halo_enabled()
!*******************************************************************************
implicit none

sgs_overlap_dwdz_halo_enabled = .true.

end function sgs_overlap_dwdz_halo_enabled

!*******************************************************************************
logical function sgs_calc_sij_explicit_enabled()
!*******************************************************************************
use param, only : coord
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .true.
logical, save :: printed = .false.

if (.not. initialized) then
    enabled = sgs_env_enabled_unless_false('LESGO_SGS_CALCSIJ_EXPLICIT',       &
        .true.)
    initialized = .true.
end if
if (coord == 0 .and. .not. printed) then
    write(*,*) 'SGS calc_Sij explicit kernel path enabled=', enabled
    printed = .true.
end if

sgs_calc_sij_explicit_enabled = enabled

end function sgs_calc_sij_explicit_enabled

!*******************************************************************************
logical function sgs_calc_sij_device_bench_enabled()
!*******************************************************************************
use param, only : coord
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.
logical, save :: printed = .false.

if (.not. initialized) then
    enabled = sgs_env_enabled_unless_false('LESGO_SGS_CALCSIJ_DEVICE_BENCH',   &
        .false.)
    initialized = .true.
end if
if (coord == 0 .and. .not. printed) then
    write(*,*) 'SGS calc_Sij device lower-bound bench enabled=', enabled
    printed = .true.
end if

sgs_calc_sij_device_bench_enabled = enabled

end function sgs_calc_sij_device_bench_enabled

!*******************************************************************************
logical function sgs_explicit_pointwise_enabled()
!*******************************************************************************
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .true.

if (.not. initialized) then
    enabled = sgs_env_enabled_unless_false('LESGO_SGS_EXPLICIT_POINTWISE',     &
        .true.)
    initialized = .true.
end if

sgs_explicit_pointwise_enabled = enabled

end function sgs_explicit_pointwise_enabled

subroutine sgs_cuda_sync(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

if (sgs_strict_sync_enabled() .or. sgs_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    if (istat /= 0) then
        print *, 'sgs_stag CUDA sync failure at ', trim(where), ': ', istat
        stop
    end if
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'sgs_stag CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine sgs_cuda_sync

!*******************************************************************************
subroutine sgs_cuda_barrier(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'sgs_stag CUDA barrier failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'sgs_stag CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine sgs_cuda_barrier

!*******************************************************************************
subroutine sgs_event_record(event, where)
!*******************************************************************************
implicit none

type(cudaEvent), intent(inout) :: event
character(len=*), intent(in) :: where
integer :: istat

istat = cudaEventRecord(event, 0)
if (istat /= 0) then
    print *, 'sgs_stag CUDA event record failure at ', trim(where), ': ', istat
    stop
end if

end subroutine sgs_event_record

!*******************************************************************************
subroutine sgs_event_elapsed_seconds(start_event, stop_event, elapsed, where)
!*******************************************************************************
implicit none

type(cudaEvent), intent(inout) :: start_event, stop_event
real(rprec), intent(out) :: elapsed
character(len=*), intent(in) :: where
integer :: istat
real :: elapsed_ms

istat = cudaEventSynchronize(stop_event)
if (istat /= 0) then
    print *, 'sgs_stag CUDA event synchronize failure at ', trim(where), ': ', &
        istat
    stop
end if
istat = cudaEventElapsedTime(elapsed_ms, start_event, stop_event)
if (istat /= 0) then
    print *, 'sgs_stag CUDA event elapsed failure at ', trim(where), ': ',     &
        istat
    stop
end if
elapsed = real(elapsed_ms, rprec) * 1.0e-3_rprec

end subroutine sgs_event_elapsed_seconds

!*******************************************************************************
subroutine sgs_diag_time(tnow)
!*******************************************************************************
#ifdef PPMPI
use mpi
#endif
implicit none

real(rprec), intent(out) :: tnow

#ifdef PPMPI
tnow = mpi_wtime()
#else
call cpu_time(tnow)
#endif

end subroutine sgs_diag_time

!*******************************************************************************
subroutine sgs_calc_diag_begin()
!*******************************************************************************
implicit none

integer :: i

sgs_calc_cpu = 0._rprec
sgs_calc_gpu = 0._rprec
sgs_calc_t0 = 0._rprec
sgs_calc_launches = 0
sgs_calc_total_launches = 0
sgs_calc_cells = 0_8
sgs_calc_total_cells = 0_8
sgs_calc_jz_min = 0
sgs_calc_jz_max = 0
sgs_calc_bulk_jz_min = 0
sgs_calc_bulk_jz_max = 0
sgs_calc_lbc_mode = -999
sgs_calc_ubc_mode = -999
sgs_calc_block_x = 0
sgs_calc_block_y = 0
sgs_calc_block_z = 0
sgs_calc_grid_x = 0
sgs_calc_grid_y = 0
sgs_calc_grid_z = 0
sgs_calc_bottom_branch = .false.
sgs_calc_top_branch = .false.
sgs_calc_explicit_used = .false.
sgs_calc_evt_active = .false.
if (.not. sgs_calc_evt_created) then
    do i = 1, SGS_CALC_COUNT
        if (cudaEventCreate(sgs_calc_evt_start(i)) /= 0) then
            print *, 'sgs_stag CUDA calc event create failure'
            stop
        end if
        if (cudaEventCreate(sgs_calc_evt_stop(i)) /= 0) then
            print *, 'sgs_stag CUDA calc event create failure'
            stop
        end if
    end do
    sgs_calc_evt_created = .true.
end if

end subroutine sgs_calc_diag_begin

!*******************************************************************************
subroutine sgs_calc_diag_start(stage, where)
!*******************************************************************************
implicit none

integer, intent(in) :: stage
character(len=*), intent(in) :: where

call sgs_diag_time(sgs_calc_t0(stage))
call sgs_event_record(sgs_calc_evt_start(stage), where)

end subroutine sgs_calc_diag_start

!*******************************************************************************
subroutine sgs_calc_diag_stop(stage, cells, where)
!*******************************************************************************
implicit none

integer, intent(in) :: stage
integer(8), intent(in) :: cells
character(len=*), intent(in) :: where
real(rprec) :: tnow

call sgs_event_record(sgs_calc_evt_stop(stage), where)
call sgs_diag_time(tnow)
sgs_calc_cpu(stage) = sgs_calc_cpu(stage) + tnow - sgs_calc_t0(stage)
sgs_calc_evt_active(stage) = .true.
sgs_calc_launches(stage) = sgs_calc_launches(stage) + 1
sgs_calc_total_launches = sgs_calc_total_launches + 1
sgs_calc_cells(stage) = sgs_calc_cells(stage) + cells
sgs_calc_total_cells = sgs_calc_total_cells + cells

end subroutine sgs_calc_diag_stop

!*******************************************************************************
subroutine sgs_calc_cpu_start(stage)
!*******************************************************************************
implicit none

integer, intent(in) :: stage

call sgs_diag_time(sgs_calc_t0(stage))

end subroutine sgs_calc_cpu_start

!*******************************************************************************
subroutine sgs_calc_cpu_stop(stage)
!*******************************************************************************
implicit none

integer, intent(in) :: stage
real(rprec) :: tnow

call sgs_diag_time(tnow)
sgs_calc_cpu(stage) = sgs_calc_cpu(stage) + tnow - sgs_calc_t0(stage)

end subroutine sgs_calc_cpu_stop

!*******************************************************************************
subroutine sgs_calc_set_zrange(jmin, jmax)
!*******************************************************************************
implicit none

integer, intent(in) :: jmin, jmax

sgs_calc_jz_min = jmin
sgs_calc_jz_max = jmax

end subroutine sgs_calc_set_zrange

!*******************************************************************************
subroutine sgs_calc_set_audit(jmin, jmax, bulk_min, bulk_max, bottom_branch,  &
    top_branch, lbc_mode, ubc_mode, explicit_used, block_x, block_y, block_z, &
    grid_x, grid_y, grid_z)
!*******************************************************************************
implicit none

integer, intent(in) :: jmin, jmax, bulk_min, bulk_max
integer, intent(in) :: lbc_mode, ubc_mode
integer, intent(in) :: block_x, block_y, block_z, grid_x, grid_y, grid_z
logical, intent(in) :: bottom_branch, top_branch, explicit_used

sgs_calc_jz_min = jmin
sgs_calc_jz_max = jmax
sgs_calc_bulk_jz_min = bulk_min
sgs_calc_bulk_jz_max = bulk_max
sgs_calc_bottom_branch = bottom_branch
sgs_calc_top_branch = top_branch
sgs_calc_lbc_mode = lbc_mode
sgs_calc_ubc_mode = ubc_mode
sgs_calc_explicit_used = explicit_used
sgs_calc_block_x = block_x
sgs_calc_block_y = block_y
sgs_calc_block_z = block_z
sgs_calc_grid_x = grid_x
sgs_calc_grid_y = grid_y
sgs_calc_grid_z = grid_z

end subroutine sgs_calc_set_audit

!*******************************************************************************
subroutine sgs_tau_detail_begin()
!*******************************************************************************
implicit none

sgs_tau_detail = 0._rprec
sgs_tau_t0 = 0._rprec
sgs_tau_send_bytes = 0_8
sgs_tau_recv_bytes = 0_8
sgs_tau_mpi_calls = 0
sgs_tau_neighbor = -999
sgs_tau_combined_msg = .false.
sgs_tau_arr_enter = 0._rprec
sgs_tau_arr_after_sync = 0._rprec
sgs_tau_arr_before_mpi = 0._rprec
sgs_tau_arr_after_mpi = 0._rprec
sgs_tau_arr_after_post = 0._rprec
sgs_tau_arr_barrier = 0._rprec

end subroutine sgs_tau_detail_begin

!*******************************************************************************
subroutine sgs_tau_detail_start(stage)
!*******************************************************************************
implicit none

integer, intent(in) :: stage

call sgs_diag_time(sgs_tau_t0)

end subroutine sgs_tau_detail_start

!*******************************************************************************
subroutine sgs_tau_detail_stop(stage)
!*******************************************************************************
implicit none

integer, intent(in) :: stage
real(rprec) :: tnow

call sgs_diag_time(tnow)
sgs_tau_detail(stage) = sgs_tau_detail(stage) + tnow - sgs_tau_t0

end subroutine sgs_tau_detail_stop

!*******************************************************************************
subroutine sgs_tau_detail_add_bytes(send_count, recv_count)
!*******************************************************************************
implicit none

integer, intent(in) :: send_count, recv_count

sgs_tau_send_bytes = sgs_tau_send_bytes + int(send_count, 8) * SGS_RPREC_BYTES
sgs_tau_recv_bytes = sgs_tau_recv_bytes + int(recv_count, 8) * SGS_RPREC_BYTES

end subroutine sgs_tau_detail_add_bytes

!*******************************************************************************
subroutine sgs_tau_detail_add_msg(send_count, recv_count, neighbor, combined)
!*******************************************************************************
implicit none

integer, intent(in) :: send_count, recv_count, neighbor
logical, intent(in) :: combined

call sgs_tau_detail_add_bytes(send_count, recv_count)
sgs_tau_mpi_calls = sgs_tau_mpi_calls + 1
sgs_tau_neighbor = neighbor
sgs_tau_combined_msg = combined

end subroutine sgs_tau_detail_add_msg

!*******************************************************************************
subroutine sgs_dwdz_detail_begin()
!*******************************************************************************
implicit none

sgs_dwdz_detail = 0._rprec
sgs_dwdz_t0 = 0._rprec
sgs_dwdz_send_bytes = 0_8
sgs_dwdz_recv_bytes = 0_8
sgs_dwdz_dummy_bytes = 0_8
sgs_dwdz_mpi_calls = 0
sgs_dwdz_neighbor = -999
sgs_dwdz_combined_msg = .false.
sgs_dwdz_arr_enter = 0._rprec
sgs_dwdz_arr_after_sync = 0._rprec
sgs_dwdz_arr_before_mpi = 0._rprec
sgs_dwdz_arr_after_mpi = 0._rprec
sgs_dwdz_arr_after_post = 0._rprec
sgs_dwdz_arr_barrier = 0._rprec

end subroutine sgs_dwdz_detail_begin

!*******************************************************************************
subroutine sgs_dwdz_detail_start(stage)
!*******************************************************************************
implicit none

integer, intent(in) :: stage

call sgs_diag_time(sgs_dwdz_t0)

end subroutine sgs_dwdz_detail_start

!*******************************************************************************
subroutine sgs_dwdz_detail_stop(stage)
!*******************************************************************************
implicit none

integer, intent(in) :: stage
real(rprec) :: tnow

call sgs_diag_time(tnow)
sgs_dwdz_detail(stage) = sgs_dwdz_detail(stage) + tnow - sgs_dwdz_t0

end subroutine sgs_dwdz_detail_stop

!*******************************************************************************
subroutine sgs_dwdz_detail_add_msg(send_count, recv_count, neighbor, combined)
!*******************************************************************************
implicit none

integer, intent(in) :: send_count, recv_count, neighbor
logical, intent(in) :: combined

sgs_dwdz_send_bytes = sgs_dwdz_send_bytes + int(send_count, 8) *              &
    SGS_RPREC_BYTES
sgs_dwdz_recv_bytes = sgs_dwdz_recv_bytes + int(recv_count, 8) *              &
    SGS_RPREC_BYTES
sgs_dwdz_mpi_calls = sgs_dwdz_mpi_calls + 1
sgs_dwdz_neighbor = neighbor
sgs_dwdz_combined_msg = combined

end subroutine sgs_dwdz_detail_add_msg

!*******************************************************************************
subroutine sgs_dwdz_path_audit(label, pattern, send_count, recv_count,         &
    neighbor, direct_section, compiler_temp_likely)
!*******************************************************************************
use param, only : coord
implicit none

character(len=*), intent(in) :: label, pattern
integer, intent(in) :: send_count, recv_count, neighbor
logical, intent(in) :: direct_section, compiler_temp_likely

if (sgs_dwdz_path_audit_printed) return
sgs_dwdz_path_audit_printed = .true.

write(*,'(a,i0,2a,2a,3(a,i0),3(a,l1))')                                     &
    'SGS_DWDZ_PATH_AUDIT rank=', coord,                                      &
    ' source=sgs_stag_util.f90 routine=calc_Sij label=', trim(label),        &
    ' pattern=', trim(pattern), ' send_bytes=',                              &
    send_count * SGS_RPREC_BYTES, ' recv_bytes=',                            &
    recv_count * SGS_RPREC_BYTES, ' neighbor=', neighbor,                    &
    ' contiguous=', .true., ' direct_array_section=', direct_section,        &
    ' compiler_temp_likely=', compiler_temp_likely
flush(6)

end subroutine sgs_dwdz_path_audit

!*******************************************************************************
character(len=12) function sgs_ptr_type_name(ptype)
!*******************************************************************************
implicit none

integer, intent(in) :: ptype

select case (ptype)
case (1)
    sgs_ptr_type_name = 'host'
case (2)
    sgs_ptr_type_name = 'device'
case (3)
    sgs_ptr_type_name = 'managed'
case default
    write(sgs_ptr_type_name,'(a,i0)') 'type', ptype
end select

end function sgs_ptr_type_name

!*******************************************************************************
subroutine sgs_pointer_env_audit(label, sendbuf, recvbuf, neighbor)
!*******************************************************************************
use param, only : coord
implicit none

character(len=*), intent(in) :: label
real(rprec), dimension(*) :: sendbuf, recvbuf
integer, intent(in) :: neighbor
type(cudaPointerAttributes) :: send_attr, recv_attr
integer :: send_stat, recv_stat, last, dev, dev_stat
character(len=64) :: mpich_gpu, visible
integer :: env_len, env_stat

if (index(label, 'dwdz') > 0) then
    if (sgs_dwdz_ptr_audit_printed) return
    sgs_dwdz_ptr_audit_printed = .true.
else
    if (sgs_tau_ptr_audit_printed) return
    sgs_tau_ptr_audit_printed = .true.
end if
call get_environment_variable('MPICH_GPU_SUPPORT_ENABLED', mpich_gpu,         &
    length=env_len, status=env_stat)
if (env_stat /= 0 .or. env_len <= 0) mpich_gpu = 'unset'
call get_environment_variable('CUDA_VISIBLE_DEVICES', visible, length=env_len,&
    status=env_stat)
if (env_stat /= 0 .or. env_len <= 0) visible = 'unset'
dev_stat = cudaGetDevice(dev)
send_attr%type = -1
send_attr%device = -1
recv_attr%type = -1
recv_attr%device = -1
send_stat = cudaPointerGetAttributes(send_attr, sendbuf)
last = cudaGetLastError()
recv_stat = cudaPointerGetAttributes(recv_attr, recvbuf)
last = cudaGetLastError()
write(*,*) 'SGS pointer/MPI audit label=', trim(label), ' rank=', coord,      &
    ' neighbor=', neighbor, ' selected_gpu=', dev, ' cudaGetDevice_status=',  &
    dev_stat, ' MPICH_GPU_SUPPORT=', trim(mpich_gpu),                         &
    ' CUDA_VISIBLE_DEVICES=', trim(visible), ' send_status=', send_stat,      &
    ' send_type=', send_attr%type, ' send_type_name=',                       &
    trim(sgs_ptr_type_name(send_attr%type)), ' send_device=',                &
    send_attr%device, ' recv_status=', recv_stat, ' recv_type=',             &
    recv_attr%type, ' recv_type_name=',                                      &
    trim(sgs_ptr_type_name(recv_attr%type)), ' recv_device=',                &
    recv_attr%device, ' MPI_GPU_SUPPORT_QUERY=not_available'

end subroutine sgs_pointer_env_audit

!*******************************************************************************
subroutine sgs_pointer_env_audit_device(label, sendbuf, recvbuf, neighbor)
!*******************************************************************************
use param, only : coord
implicit none

character(len=*), intent(in) :: label
real(rprec), device, dimension(*) :: sendbuf, recvbuf
integer, intent(in) :: neighbor
type(cudaPointerAttributes) :: send_attr, recv_attr
integer :: send_stat, recv_stat, last, dev, dev_stat
character(len=64) :: mpich_gpu, visible
integer :: env_len, env_stat

if (index(label, 'dwdz') > 0) then
    if (sgs_dwdz_ptr_audit_printed) return
    sgs_dwdz_ptr_audit_printed = .true.
else
    if (sgs_tau_ptr_audit_printed) return
    sgs_tau_ptr_audit_printed = .true.
end if
call get_environment_variable('MPICH_GPU_SUPPORT_ENABLED', mpich_gpu,         &
    length=env_len, status=env_stat)
if (env_stat /= 0 .or. env_len <= 0) mpich_gpu = 'unset'
call get_environment_variable('CUDA_VISIBLE_DEVICES', visible, length=env_len,&
    status=env_stat)
if (env_stat /= 0 .or. env_len <= 0) visible = 'unset'
dev_stat = cudaGetDevice(dev)
send_attr%type = -1
send_attr%device = -1
recv_attr%type = -1
recv_attr%device = -1
send_stat = cudaPointerGetAttributes(send_attr, sendbuf)
last = cudaGetLastError()
recv_stat = cudaPointerGetAttributes(recv_attr, recvbuf)
last = cudaGetLastError()
write(*,*) 'SGS pointer/MPI audit label=', trim(label), ' rank=', coord,      &
    ' neighbor=', neighbor, ' selected_gpu=', dev, ' cudaGetDevice_status=',  &
    dev_stat, ' MPICH_GPU_SUPPORT=', trim(mpich_gpu),                         &
    ' CUDA_VISIBLE_DEVICES=', trim(visible), ' send_status=', send_stat,      &
    ' send_type=', send_attr%type, ' send_type_name=',                       &
    trim(sgs_ptr_type_name(send_attr%type)), ' send_device=',                &
    send_attr%device, ' recv_status=', recv_stat, ' recv_type=',             &
    recv_attr%type, ' recv_type_name=',                                      &
    trim(sgs_ptr_type_name(recv_attr%type)), ' recv_device=',                &
    recv_attr%device, ' MPI_GPU_SUPPORT_QUERY=not_available'

end subroutine sgs_pointer_env_audit_device
#endif

!*******************************************************************************
subroutine sgs_stag ()
!*******************************************************************************
!
! Calculates turbulent (subgrid) stress for entire domain
!   using the model specified in param.f90 (Smag, LASD, etc)
!   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)
!   Sij values are stored on w-nodes (1:nz)
!
!   module is used to share Sij values b/w subroutines
!     (avoids memory error when arrays are very large)
!
! put everything onto w-nodes, follow original version

use types, only : rprec
use param
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
use sgs_param
use messages

#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#ifdef ENABLE_CUDA
use mpi
use cuda_mpi_debug, only : mpi_dbg_sendrecv_r
#endif
#endif

#ifdef PPLVLSET
use level_set, only : level_set_BC, level_set_Cs
#endif

implicit none

character (*), parameter :: sub_name = 'sgs_stag'

real(rprec), dimension(nz) :: l, ziko, zz
#ifdef ENABLE_CUDA
real(rprec), managed, save, allocatable, dimension(:) :: l_cuda
real(rprec), managed, save, allocatable, dimension(:,:,:) ::                 &
    tau_halo_send, tau_halo_recv
#endif
real(rprec) :: const, const2, const3, const4
real(rprec) :: l_delta2
real(rprec) :: S_mag
integer :: jx, jy, jz
integer :: jz_min, jz_max
#if defined(ENABLE_CUDA) && defined(PPMPI)
integer :: tau_halo_req(4)
integer :: tau_halo_status(MPI_STATUS_SIZE, 4)
integer :: tau_halo_nreq
integer :: tau_int_min, tau_int_max
integer :: tau2_count, tau2_peer, tau2_idx
logical :: tau_halo_posted
logical :: use_tau_overlap
#endif
#ifdef ENABLE_CUDA
logical :: dynamic_coeff_update
logical :: use_fused_dynamic_nut
logical :: use_explicit_pointwise
logical :: sgs_stage_enabled_local
integer, save :: sgs_stage_count = 0
integer, parameter :: SGS_EVT_CALC = 1, SGS_EVT_DYNAMIC = 2
integer, parameter :: SGS_EVT_NUT = 3, SGS_EVT_TAU_BOUNDARY = 4
integer, parameter :: SGS_EVT_TAU_INTERIOR = 5, SGS_EVT_COUNT = 5
type(cudaEvent) :: sgs_evt_start(SGS_EVT_COUNT), sgs_evt_stop(SGS_EVT_COUNT)
logical :: sgs_event_active
integer :: sgs_event_i
integer :: sgs_pw_block_x, sgs_pw_block_y, sgs_pw_block_z
integer :: sgs_pw_grid_x, sgs_pw_grid_y, sgs_pw_grid_z
type(dim3) :: sgs_pw_block, sgs_pw_grid
real(rprec) :: sgs_t0, sgs_t1, sgs_barrier_t0, sgs_barrier_t1
real(rprec) :: sgs_stage_calc, sgs_stage_dynamic, sgs_stage_nut
real(rprec) :: sgs_stage_tau_boundary, sgs_stage_tau_interior
real(rprec) :: sgs_stage_tau_halo, sgs_stage_final
real(rprec) :: sgs_gpu_calc, sgs_gpu_dynamic, sgs_gpu_nut
real(rprec) :: sgs_gpu_tau_boundary, sgs_gpu_tau_interior
#endif

! Cs is Smagorinsky's constant. l is a filter size (non-dim.)
#ifdef ENABLE_CUDA
sgs_stage_enabled_local = sgs_stage_timing_enabled()
sgs_event_active = sgs_stage_enabled_local .and. sgs_pointwise_cuda_enabled()
sgs_stage_calc = 0._rprec
sgs_stage_dynamic = 0._rprec
sgs_stage_nut = 0._rprec
sgs_stage_tau_boundary = 0._rprec
sgs_stage_tau_interior = 0._rprec
sgs_stage_tau_halo = 0._rprec
sgs_stage_final = 0._rprec
sgs_gpu_calc = 0._rprec
sgs_gpu_dynamic = 0._rprec
sgs_gpu_nut = 0._rprec
sgs_gpu_tau_boundary = 0._rprec
sgs_gpu_tau_interior = 0._rprec
if (sgs_stage_enabled_local) call sgs_tau_detail_begin()
if (sgs_stage_enabled_local) call sgs_dwdz_detail_begin()
if (sgs_event_active) then
    do sgs_event_i = 1, SGS_EVT_COUNT
        if (cudaEventCreate(sgs_evt_start(sgs_event_i)) /= 0) then
            print *, 'sgs_stag CUDA event create failure'
            stop
        end if
        if (cudaEventCreate(sgs_evt_stop(sgs_event_i)) /= 0) then
            print *, 'sgs_stag CUDA event create failure'
            stop
        end if
    end do
end if
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t0 = mpi_wtime()
#else
    call cpu_time(sgs_t0)
#endif
end if

dynamic_coeff_update = .false.
#if defined(ENABLE_CUDA) && defined(PPMPI)
tau_halo_posted = .false.
use_tau_overlap = .false.
#endif
if (sgs .and. (sgs_model /= SGS_MODEL_SMAGORINSKY)) then
    dynamic_coeff_update = ((jt >= DYN_init) .or. initu) .and.                &
        (mod(jt_total, cs_count) == 0)
end if
use_fused_dynamic_nut = sgs_pointwise_cuda_enabled() .and.                    &
    sgs_fused_dynamic_nut_enabled() .and. sgs .and.                           &
    (sgs_model /= SGS_MODEL_SMAGORINSKY) .and.                                &
    (.not. dynamic_coeff_update) .and. (lbc_mom > 0) .and. (ubc_mom == 0)
use_explicit_pointwise = sgs_pointwise_cuda_enabled() .and.                  &
    sgs_explicit_pointwise_enabled() .and. sgs .and.                          &
    (sgs_model /= SGS_MODEL_SMAGORINSKY)
#ifdef PPLVLSET
use_fused_dynamic_nut = .false.
use_explicit_pointwise = .false.
#endif
#if defined(ENABLE_CUDA) && defined(PPMPI)
use_tau_overlap = sgs_pointwise_cuda_enabled() .and.                          &
    sgs_overlap_tau_halo_enabled() .and. sgs_combined_tau_halo_active() .and. &
    (nproc > 1) .and. (.not. sgs_halo_combined2_enabled()) .and.             &
    (.not. sgs_halo_hostpinned_test_enabled())
#ifdef PPLVLSET
use_tau_overlap = .false.
#endif
#endif
if (sgs_event_active) call sgs_event_record(sgs_evt_start(SGS_EVT_CALC),      &
    'calc_Sij start')
if (.not. use_fused_dynamic_nut) call calc_Sij ()
if (sgs_event_active) call sgs_event_record(sgs_evt_stop(SGS_EVT_CALC),       &
    'calc_Sij stop')
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_calc = sgs_t1 - sgs_t0
    sgs_t0 = sgs_t1
end if
#else
call calc_Sij ()
#endif

#ifdef ENABLE_CUDA
if (sgs_event_active) call sgs_event_record(sgs_evt_start(SGS_EVT_DYNAMIC),   &
    'dynamic model start')
#endif
! This approximates the sum displacement during cs_count timesteps
! This is used with the lagrangian model only
if (use_cfl_dt) then
    if (sgs_model == SGS_MODEL_LAGRANGE_SIMILARITY .OR.                       &
        sgs_model == SGS_MODEL_LAGRANGE_SCALE_DEP) then
        if ( ( jt .GE. DYN_init-cs_count + 1 ) .OR.  initu ) then
            lagran_dt = lagran_dt + dt
        endif
    endif
else
lagran_dt = cs_count*dt
end if

if (sgs) then
    ! Traditional Smagorinsky model
    if (sgs_model == SGS_MODEL_SMAGORINSKY) then

#ifdef PPLVLSET
        l = delta
        call level_set_Cs (delta)
#else
        ! Parameters (Co and nn) for wallfunction defined in param.f90
        Cs_opt2 = Co**2  ! constant coefficient

        ! both Stress free
        if (lbc_mom == 0 .and. ubc_mom == 0) then
            l = delta

        ! top Stress free, bottom wall
        else if (lbc_mom > 0 .and. ubc_mom == 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == 0) then
                ! z's nondimensional, l here is on uv-nodes
                zz(1) = 0.5_rprec * dz
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)    &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_min = 2
            else
                jz_min = 1
            end if

            do jz = jz_min, nz
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                l(jz) = ( Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)  &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! both top and bottom walls, zz is distance to nearest wall
        else if (lbc_mom > 0 .and. ubc_mom > 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == 0) then
                ! z's nondimensional, l here is on uv-nodes
                zz(1) = 0.5_rprec * dz
                l(1) = ( Co**(wall_damp_exp)*(vonk*zz(1))**(-wall_damp_exp)&
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)

                jz_min = 2
            else
                jz_min = 1
            end if

            if (coord == nproc-1) then
                ! z's nondimensional, l here is on uv-nodes
                zz(nz) = 0.5_rprec * dz
                l(nz) = (Co**(wall_damp_exp)*(vonk*zz(nz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_max = nz-1
            else
                jz_max = nz
            end if

            do jz = jz_min, jz_max
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz
                zz(jz) = min( zz(jz), (nz-1)*nproc*dz - zz(jz) )
                l(jz) = (Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! top wall, bottom stress free, zz is distance to top
        else if (lbc_mom == 0 .and. ubc_mom > 0) then
            ! The variable "l" calculated below is l_sgs/Co
            ! l_sgs is from JDA eqn(2.30)
            if (coord == nproc-1) then
                ! z's nondimensional, l here is on uv-nodes
                zz(nz) = 0.5_rprec * dz
                l(nz) = (Co**(wall_damp_exp)*(vonk*zz(nz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
                jz_max = nz-1
            else
                jz_max = nz
            end if

            do jz = 1, jz_max
                ! z's nondimensional, l here is on w-nodes
                zz(jz) = ((nproc - coord)*(nz - 1) - (jz - 1)) * dz
                l(jz) = (Co**(wall_damp_exp)*(vonk*zz(jz))**(-wall_damp_exp)   &
                    + (delta)**(-wall_damp_exp) )**(-1._rprec/wall_damp_exp)
            end do

        ! Invalid combination
        else
            call error (sub_name, 'invalid b.c. combination')
        end if
#endif

    ! Dynamic procedures: modify/set Sij and Cs_opt2 (specific to sgs_model)
    else
        ! recall: l is the filter size
        l = delta

        ! Use the Smagorinsky model until DYN_init timestep
        if ((jt == 1) .and. (inilag)) then
            write(*,*) 'CS_opt2 initialiazed'
            Cs_opt2 = 0.03_rprec

        ! Update Sij, Cs every cs_count timesteps (specified in param)
        elseif ( ((jt.GE.DYN_init).OR.(initu)) .AND.                           &
            (mod(jt_total,cs_count)==0) ) then

            if (jt == DYN_init) then
                write(*,*) 'running dynamic sgs_model = ', sgs_model
            end if

            ! Standard dynamic model
            if (sgs_model == SGS_MODEL_STANDARD_DYNAMIC) then
                call std_dynamic(ziko)
                forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)

            ! Plane average dynamic model
            else if (sgs_model == SGS_MODEL_SCALE_DEP_DYNAMIC) then
                call scaledep_dynamic(ziko)
                do jz = 1, nz
                    Cs_opt2(:, :, jz) = ziko(jz)
                end do
            ! Lagrangian scale similarity model
            else if (sgs_model == SGS_MODEL_LAGRANGE_SIMILARITY) then
                call lagrange_Ssim()

            ! Lagrangian scale dependent model
            elseif (sgs_model == SGS_MODEL_LAGRANGE_SCALE_DEP) then
                call lagrange_Sdep()
            end if

        end if
    end if
end if

#ifdef ENABLE_CUDA
if (sgs_event_active) call sgs_event_record(sgs_evt_stop(SGS_EVT_DYNAMIC),    &
    'dynamic model stop')
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_dynamic = sgs_t1 - sgs_t0
    sgs_t0 = sgs_t1
end if
#endif


! Define |S| and eddy viscosity (nu_t= c_s^2 l^2 |S|) for entire domain
!   stored on w-nodes (on uvp node for jz=1 or nz for 'wall' BC only)
#ifdef ENABLE_CUDA
if (sgs_event_active) call sgs_event_record(sgs_evt_start(SGS_EVT_NUT),       &
    'Nu_t start')
if (use_fused_dynamic_nut) then
    call calc_Sij_nut_dynamic_cuda()
else if (sgs_pointwise_cuda_enabled()) then
    if (sgs_model == SGS_MODEL_SMAGORINSKY) then
        if (.not. allocated(l_cuda)) allocate(l_cuda(nz))
        l_cuda = l
    end if
    l_delta2 = delta * delta
    if (use_explicit_pointwise) then
        sgs_pw_block_x = 32
        sgs_pw_block_y = 4
        sgs_pw_block_z = 2
        sgs_pw_grid_x = (nx + sgs_pw_block_x - 1) / sgs_pw_block_x
        sgs_pw_grid_y = (ny + sgs_pw_block_y - 1) / sgs_pw_block_y
        sgs_pw_grid_z = (nz + sgs_pw_block_z - 1) / sgs_pw_block_z
        sgs_pw_block = dim3(sgs_pw_block_x, sgs_pw_block_y, sgs_pw_block_z)
        sgs_pw_grid = dim3(sgs_pw_grid_x, sgs_pw_grid_y, sgs_pw_grid_z)
        call sgs_nut_dynamic_kernel<<<sgs_pw_grid,sgs_pw_block>>>(S11, S12,  &
            S13, S22, S23, S33, Cs_opt2, Nu_t, S, l_delta2, ld, nx, ny, nz)
    else
!$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, nx
        S_mag = sqrt( 2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +        &
            S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +                 &
            S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2 )))
        if (sgs_model == SGS_MODEL_SMAGORINSKY) then
            Nu_t(jx,jy,jz) = S_mag * Cs_opt2(jx,jy,jz) *                     &
                l_cuda(jz) * l_cuda(jz)
        else
            Nu_t(jx,jy,jz) = S_mag * Cs_opt2(jx,jy,jz) * l_delta2
        end if
        if (jz == nz) S(jx,jy) = S_mag
    end do
    end do
    end do
    end if
    call sgs_cuda_sync('Nu_t')
else
#endif
do jz = 1, nz
do jy = 1, ny
do jx = 1, nx
    S(jx,jy) = sqrt( 2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +           &
        S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +                        &
        S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2 )))
    Nu_t(jx,jy,jz) = S(jx,jy)*Cs_opt2(jx,jy,jz)*l(jz)**2
end do
end do
end do
#ifdef ENABLE_CUDA
end if
if (sgs_event_active) call sgs_event_record(sgs_evt_stop(SGS_EVT_NUT),        &
    'Nu_t stop')
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_nut = sgs_t1 - sgs_t0
    sgs_t0 = sgs_t1
end if
#endif

! Calculate txx, txy, tyy, tzz for bottom level: jz=1 node (coord==0 only)
#ifdef ENABLE_CUDA
if (sgs_event_active) call sgs_event_record(                                 &
    sgs_evt_start(SGS_EVT_TAU_BOUNDARY), 'tau boundary start')
#endif
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on w-nodes
        case (0)
            if (sgs) then
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        const = 0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2)) + nu
                        txx(jx,jy,1) = -const*(S11(jx,jy,1) + S11(jx,jy,2))
                        txy(jx,jy,1) = -const*(S12(jx,jy,1) + S12(jx,jy,2))
                        tyy(jx,jy,1) = -const*(S22(jx,jy,1) + S22(jx,jy,2))
                        tzz(jx,jy,1) = -const*(S33(jx,jy,1) + S33(jx,jy,2))
                    end do
                    end do
                    call sgs_cuda_sync('tau bottom stress-free sgs')
                else
#endif
                do jy = 1, ny
                do jx = 1, nx
                    ! Total viscosity
                    const = 0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2)) + nu
                    txx(jx,jy,1) = -const*(S11(jx,jy,1) + S11(jx,jy,2))
                    txy(jx,jy,1) = -const*(S12(jx,jy,1) + S12(jx,jy,2))
                    tyy(jx,jy,1) = -const*(S22(jx,jy,1) + S22(jx,jy,2))
                    tzz(jx,jy,1) = -const*(S33(jx,jy,1) + S33(jx,jy,2))
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            else
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,1) = -(nu)*(S11(jx,jy,1) + S11(jx,jy,2))
                        txy(jx,jy,1) = -(nu)*(S12(jx,jy,1) + S12(jx,jy,2))
                        tyy(jx,jy,1) = -(nu)*(S22(jx,jy,1) + S22(jx,jy,2))
                        tzz(jx,jy,1) = -(nu)*(S33(jx,jy,1) + S33(jx,jy,2))
                    end do
                    end do
                    call sgs_cuda_sync('tau bottom stress-free molecular')
                else
#endif
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,1) = -(nu)*(S11(jx,jy,1) + S11(jx,jy,2))
                    txy(jx,jy,1) = -(nu)*(S12(jx,jy,1) + S12(jx,jy,2))
                    tyy(jx,jy,1) = -(nu)*(S22(jx,jy,1) + S22(jx,jy,2))
                    tzz(jx,jy,1) = -(nu)*(S33(jx,jy,1) + S33(jx,jy,2))
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            end if

        ! Wall
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on uvp-nodes
        case (1:)
            if (sgs) then
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        const = -2._rprec*(Nu_t(jx,jy,1)+nu)
                        txx(jx,jy,1) = const*S11(jx,jy,1)
                        txy(jx,jy,1) = const*S12(jx,jy,1)
                        tyy(jx,jy,1) = const*S22(jx,jy,1)
                        tzz(jx,jy,1) = const*S33(jx,jy,1)
                    end do
                    end do
                    call sgs_cuda_sync('tau bottom wall sgs')
                else
#endif
                do jy = 1, ny
                do jx = 1, nx
                    const = -2._rprec*(Nu_t(jx,jy,1)+nu)
                    txx(jx,jy,1) = const*S11(jx,jy,1)
                    txy(jx,jy,1) = const*S12(jx,jy,1)
                    tyy(jx,jy,1) = const*S22(jx,jy,1)
                    tzz(jx,jy,1) = const*S33(jx,jy,1)
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            else
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,1) = -2._rprec*(nu)*S11(jx,jy,1)
                        txy(jx,jy,1) = -2._rprec*(nu)*S12(jx,jy,1)
                        tyy(jx,jy,1) = -2._rprec*(nu)*S22(jx,jy,1)
                        tzz(jx,jy,1) = -2._rprec*(nu)*S33(jx,jy,1)
                    end do
                    end do
                    call sgs_cuda_sync('tau bottom wall molecular')
                else
#endif
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,1) = -2._rprec*(nu)*S11(jx,jy,1)
                    txy(jx,jy,1) = -2._rprec*(nu)*S12(jx,jy,1)
                    tyy(jx,jy,1) = -2._rprec*(nu)*S22(jx,jy,1)
                    tzz(jx,jy,1) = -2._rprec*(nu)*S33(jx,jy,1)
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            end if

    end select

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if


! Calculate txx, txy, tyy, tzz for bottom level: jz=nz node (coord==nproc-1)
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on w-nodes
        case (0)

            if (sgs) then
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        const = 0.5_rprec*(Nu_t(jx,jy,nz-1) + Nu_t(jx,jy,nz)) + nu
                        const2 = 2._rprec*(Nu_t(jx,jy,nz-1) + nu)
                        txx(jx,jy,nz-1) = -const*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                        txy(jx,jy,nz-1) = -const*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                        tyy(jx,jy,nz-1) = -const*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                        tzz(jx,jy,nz-1) = -const*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                        txz(jx,jy,nz-1) = -const2*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = -const2*S23(jx,jy,nz-1)
                    end do
                    end do
                    call sgs_cuda_sync('tau top stress-free sgs')
                else
#endif
                do jy = 1, ny
                do jx = 1, nx
                    ! Total viscosity
                    const = 0.5_rprec*(Nu_t(jx,jy,nz-1) + Nu_t(jx,jy,nz)) + nu
                    const2 = 2._rprec*(Nu_t(jx,jy,nz-1) + nu)

                    ! for top wall, it is nz-1 on the uv-grid
                    txx(jx,jy,nz-1) = -const*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                    txy(jx,jy,nz-1) = -const*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                    tyy(jx,jy,nz-1) = -const*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                    tzz(jx,jy,nz-1) = -const*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1) = -const2*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1) = -const2*S23(jx,jy,nz-1)
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            else
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,nz-1) = -(nu)*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                        txy(jx,jy,nz-1) = -(nu)*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                        tyy(jx,jy,nz-1) = -(nu)*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                        tzz(jx,jy,nz-1) = -(nu)*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                        txz(jx,jy,nz-1) = -2._rprec*(nu)*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = -2._rprec*(nu)*S23(jx,jy,nz-1)
                    end do
                    end do
                    call sgs_cuda_sync('tau top stress-free molecular')
                else
#endif
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,nz-1) = -(nu)*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                    txy(jx,jy,nz-1) = -(nu)*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                    tyy(jx,jy,nz-1) = -(nu)*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                    tzz(jx,jy,nz-1) = -(nu)*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1) = -2._rprec*(nu)*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1) = -2._rprec*(nu)*S23(jx,jy,nz-1)
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            end if

        ! Wall
        ! txx,txy,tyy,tzz stored on uvp-nodes (for this and all levels)
        !   recall: for this case, Sij are stored on uvp-nodes
        case (1:)
            if (sgs) then
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        const = -2._rprec*(Nu_t(jx,jy,nz)+nu)
                        const2 = -2._rprec*(Nu_t(jx,jy,nz-1) + nu)
                        txx(jx,jy,nz-1) = const*S11(jx,jy,nz)
                        txy(jx,jy,nz-1) = const*S12(jx,jy,nz)
                        tyy(jx,jy,nz-1) = const*S22(jx,jy,nz)
                        tzz(jx,jy,nz-1) = const*S33(jx,jy,nz)
                        txz(jx,jy,nz-1) = const2*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = const2*S23(jx,jy,nz-1)
                    end do
                    end do
                    call sgs_cuda_sync('tau top wall sgs')
                else
#endif
                do jy = 1, ny
                do jx = 1, nx
                    const = -2._rprec*(Nu_t(jx,jy,nz)+nu)
                    const2 = -2._rprec*(Nu_t(jx,jy,nz-1) + nu)

                    ! Note: Sij(nz) is on uvp-node at nz-1
                    txx(jx,jy,nz-1) = const*S11(jx,jy,nz)
                    txy(jx,jy,nz-1) = const*S12(jx,jy,nz)
                    tyy(jx,jy,nz-1) = const*S22(jx,jy,nz)
                    tzz(jx,jy,nz-1) = const*S33(jx,jy,nz)
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1)= const2*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1)= const2*S23(jx,jy,nz-1)
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            else
#ifdef ENABLE_CUDA
                if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
                    do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,nz-1) = -2._rprec*(nu)*S11(jx,jy,nz-1)
                        txy(jx,jy,nz-1) = -2._rprec*(nu)*S12(jx,jy,nz-1)
                        tyy(jx,jy,nz-1) = -2._rprec*(nu)*S22(jx,jy,nz-1)
                        tzz(jx,jy,nz-1) = -2._rprec*(nu)*S33(jx,jy,nz-1)
                        txz(jx,jy,nz-1) = -2._rprec*(nu)*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = -2._rprec*(nu)*S23(jx,jy,nz-1)
                    end do
                    end do
                    call sgs_cuda_sync('tau top wall molecular')
                else
#endif
                const = 0._rprec
                do jy = 1, ny
                do jx = 1, nx
                    txx(jx,jy,nz-1) = -2._rprec*(nu)*S11(jx,jy,nz-1)
                    txy(jx,jy,nz-1) = -2._rprec*(nu)*S12(jx,jy,nz-1)
                    tyy(jx,jy,nz-1) = -2._rprec*(nu)*S22(jx,jy,nz-1)
                    tzz(jx,jy,nz-1) = -2._rprec*(nu)*S33(jx,jy,nz-1)
                    ! for top wall, include w-grid stress since we touched nz-1
                    txz(jx,jy,nz-1)=-2._rprec*(nu)*S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1)=-2._rprec*(nu)*S23(jx,jy,nz-1)
                end do
                end do
#ifdef ENABLE_CUDA
                end if
#endif
            end if

    end select

    ! since last level already calculated
    jz_max = nz-2
else
    jz_max = nz-1
end if

#ifdef ENABLE_CUDA
if (sgs_event_active) call sgs_event_record(                                 &
    sgs_evt_stop(SGS_EVT_TAU_BOUNDARY), 'tau boundary stop')
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_tau_boundary = sgs_t1 - sgs_t0
    sgs_t0 = sgs_t1
end if
#endif

! Calculate all tau for the rest of the domain
!   txx, txy, tyy, tzz not needed at nz (so they aren't calculated)
!     txz, tyz at nz will be done later
!   txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)

if (sgs) then
    const3=-2._rprec*(nu)*0.5_rprec
    const4=-2._rprec*(nu)
#ifdef ENABLE_CUDA
#ifdef PPMPI
    if (use_tau_overlap) then
        if (.not. allocated(tau_halo_send)) then
            allocate(tau_halo_send(ld, ny, 3), tau_halo_recv(ld, ny, 3))
        else if (size(tau_halo_send, 1) /= ld .or. size(tau_halo_send, 3) /= 3) then
            deallocate(tau_halo_send, tau_halo_recv)
            allocate(tau_halo_send(ld, ny, 3), tau_halo_recv(ld, ny, 3))
        end if
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PACK)
!$cuf kernel do(2) <<<*,*>>>
        do jy=1,ny
        do jx=1,ld
            if (coord > 0) then
                if (jx <= nx) then
                    const =-0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2))
                    const2=-2._rprec*Nu_t(jx,jy,1)
                    txx(jx,jy,1)=(const+const3)*(S11(jx,jy,1) + S11(jx,jy,2))
                    txy(jx,jy,1)=(const+const3)*(S12(jx,jy,1) + S12(jx,jy,2))
                    tyy(jx,jy,1)=(const+const3)*(S22(jx,jy,1) + S22(jx,jy,2))
                    tzz(jx,jy,1)=(const+const3)*(S33(jx,jy,1) + S33(jx,jy,2))
                    txz(jx,jy,1)=(const2+const4)* S13(jx,jy,1)
                    tyz(jx,jy,1)=(const2+const4)* S23(jx,jy,1)
                end if
                tau_halo_send(jx,jy,1) = txz(jx,jy,1)
                tau_halo_send(jx,jy,2) = tyz(jx,jy,1)
            end if
            if (coord < nproc - 1) then
                if (jx <= nx) then
                    const =-0.5_rprec*(Nu_t(jx,jy,nz-1) + Nu_t(jx,jy,nz))
                    const2=-2._rprec*Nu_t(jx,jy,nz-1)
                    txx(jx,jy,nz-1)=(const+const3)*(S11(jx,jy,nz-1) +         &
                        S11(jx,jy,nz))
                    txy(jx,jy,nz-1)=(const+const3)*(S12(jx,jy,nz-1) +         &
                        S12(jx,jy,nz))
                    tyy(jx,jy,nz-1)=(const+const3)*(S22(jx,jy,nz-1) +         &
                        S22(jx,jy,nz))
                    tzz(jx,jy,nz-1)=(const+const3)*(S33(jx,jy,nz-1) +         &
                        S33(jx,jy,nz))
                    txz(jx,jy,nz-1)=(const2+const4)* S13(jx,jy,nz-1)
                    tyz(jx,jy,nz-1)=(const2+const4)* S23(jx,jy,nz-1)
                end if
                tau_halo_send(jx,jy,3) = tzz(jx,jy,nz-1)
            end if
        end do
        end do
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PACK)
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PRE_SYNC)
        call sgs_cuda_barrier('tau overlap send planes')
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PRE_SYNC)
        if (sgs_stage_enabled_local) call sgs_pointer_env_audit(              &
            'tau_overlap_managed', tau_halo_send(1,1,1), tau_halo_recv(1,1,1),&
            merge(up, down, coord < nproc - 1))

        tau_halo_req = MPI_REQUEST_NULL
        tau_halo_nreq = 0
        if (coord < nproc - 1) then
            if (sgs_stage_enabled_local) then
                call sgs_tau_detail_add_msg(0, 2*ld*ny, up, .false.)
                call sgs_tau_detail_add_msg(ld*ny, 0, up, .false.)
            end if
            tau_halo_nreq = tau_halo_nreq + 1
            call mpi_irecv(tau_halo_recv(1,1,1), 2*ld*ny, MPI_RPREC, up, 81,  &
                comm, tau_halo_req(tau_halo_nreq), ierr)
            tau_halo_nreq = tau_halo_nreq + 1
            call mpi_isend(tau_halo_send(1,1,3), ld*ny, MPI_RPREC, up, 82,    &
                comm, tau_halo_req(tau_halo_nreq), ierr)
        end if
        if (coord > 0) then
            if (sgs_stage_enabled_local) then
                call sgs_tau_detail_add_msg(0, ld*ny, down, .false.)
                call sgs_tau_detail_add_msg(2*ld*ny, 0, down, .false.)
            end if
            tau_halo_nreq = tau_halo_nreq + 1
            call mpi_irecv(tzz(1,1,0), ld*ny, MPI_RPREC, down, 82, comm,      &
                tau_halo_req(tau_halo_nreq), ierr)
            tau_halo_nreq = tau_halo_nreq + 1
            call mpi_isend(tau_halo_send(1,1,1), 2*ld*ny, MPI_RPREC, down, 81,&
                comm, tau_halo_req(tau_halo_nreq), ierr)
        end if
        tau_halo_posted = .true.

        tau_int_min = jz_min
        tau_int_max = jz_max
        if (coord > 0) tau_int_min = max(tau_int_min, 2)
        if (coord < nproc - 1) tau_int_max = min(tau_int_max, nz-2)
        if (tau_int_min <= tau_int_max) then
            if (sgs_event_active) call sgs_event_record(                     &
                sgs_evt_start(SGS_EVT_TAU_INTERIOR),                         &
                'tau overlap interior start')
            if (use_explicit_pointwise) then
                sgs_pw_block_x = 32
                sgs_pw_block_y = 4
                sgs_pw_block_z = 2
                sgs_pw_grid_x = (nx + sgs_pw_block_x - 1) / sgs_pw_block_x
                sgs_pw_grid_y = (ny + sgs_pw_block_y - 1) / sgs_pw_block_y
                sgs_pw_grid_z = (tau_int_max - tau_int_min + 1 +             &
                    sgs_pw_block_z - 1) / sgs_pw_block_z
                sgs_pw_block = dim3(sgs_pw_block_x, sgs_pw_block_y,          &
                    sgs_pw_block_z)
                sgs_pw_grid = dim3(sgs_pw_grid_x, sgs_pw_grid_y,             &
                    sgs_pw_grid_z)
                call sgs_tau_interior_kernel<<<sgs_pw_grid,sgs_pw_block>>>(  &
                    Nu_t, S11, S12, S13, S22, S23, S33, txx, txy, txz, tyy,  &
                    tyz, tzz, nu, ld, nx, ny, tau_int_min, tau_int_max)
            else
!$cuf kernel do(3) <<<*,*>>>
            do jz=tau_int_min, tau_int_max
            do jy=1,ny
            do jx=1,nx

               const =-0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))
               const2=-2._rprec*Nu_t(jx,jy,jz)

               txx(jx,jy,jz)=(const+const3)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
               txy(jx,jy,jz)=(const+const3)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
               tyy(jx,jy,jz)=(const+const3)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
               tzz(jx,jy,jz)=(const+const3)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
               txz(jx,jy,jz)=(const2+const4)* S13(jx,jy,jz)
               tyz(jx,jy,jz)=(const2+const4)* S23(jx,jy,jz)

            end do
            end do
            end do
            end if
            if (sgs_event_active) call sgs_event_record(                     &
                sgs_evt_stop(SGS_EVT_TAU_INTERIOR),                          &
                'tau overlap interior stop')
            call sgs_cuda_sync('tau overlap interior')
        end if
    else if (sgs_pointwise_cuda_enabled()) then
#else
    if (sgs_pointwise_cuda_enabled()) then
#endif
        if (sgs_event_active) call sgs_event_record(                         &
            sgs_evt_start(SGS_EVT_TAU_INTERIOR), 'tau interior start')
        if (use_explicit_pointwise) then
            sgs_pw_block_x = 32
            sgs_pw_block_y = 4
            sgs_pw_block_z = 2
            sgs_pw_grid_x = (nx + sgs_pw_block_x - 1) / sgs_pw_block_x
            sgs_pw_grid_y = (ny + sgs_pw_block_y - 1) / sgs_pw_block_y
            sgs_pw_grid_z = (jz_max - jz_min + 1 + sgs_pw_block_z - 1) /     &
                sgs_pw_block_z
            sgs_pw_block = dim3(sgs_pw_block_x, sgs_pw_block_y,              &
                sgs_pw_block_z)
            sgs_pw_grid = dim3(sgs_pw_grid_x, sgs_pw_grid_y, sgs_pw_grid_z)
            call sgs_tau_interior_kernel<<<sgs_pw_grid,sgs_pw_block>>>(Nu_t, &
                S11, S12, S13, S22, S23, S33, txx, txy, txz, tyy, tyz, tzz,  &
                nu, ld, nx, ny, jz_min, jz_max)
        else
!$cuf kernel do(3) <<<*,*>>>
        do jz=jz_min, jz_max
        do jy=1,ny
        do jx=1,nx

           const =-0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))
           const2=-2._rprec*Nu_t(jx,jy,jz)

           txx(jx,jy,jz)=(const+const3)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
           txy(jx,jy,jz)=(const+const3)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
           tyy(jx,jy,jz)=(const+const3)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
           tzz(jx,jy,jz)=(const+const3)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
           txz(jx,jy,jz)=(const2+const4)* S13(jx,jy,jz)
           tyz(jx,jy,jz)=(const2+const4)* S23(jx,jy,jz)

        end do
        end do
        end do
        end if
        if (sgs_event_active) call sgs_event_record(                         &
            sgs_evt_stop(SGS_EVT_TAU_INTERIOR), 'tau interior stop')
        call sgs_cuda_sync('tau interior')
    else
#endif
    do jz=jz_min, jz_max
    do jy=1,ny
    do jx=1,nx

       const =-0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))
       const2=-2._rprec*Nu_t(jx,jy,jz)

       txx(jx,jy,jz)=(const+const3)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
       txy(jx,jy,jz)=(const+const3)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
       tyy(jx,jy,jz)=(const+const3)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
       tzz(jx,jy,jz)=(const+const3)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
       txz(jx,jy,jz)=(const2+const4)* S13(jx,jy,jz)
       tyz(jx,jy,jz)=(const2+const4)* S23(jx,jy,jz)

    end do
    end do
    end do
#ifdef ENABLE_CUDA
    end if
#endif

else
    const=0._rprec  ! removed from tij expressions below since it's zero

    do jz = jz_min, jz_max
    do jy = 1, ny
    do jx = 1, nx
        txx(jx,jy,jz)=-(nu)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
        txy(jx,jy,jz)=-(nu)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
        tyy(jx,jy,jz)=-(nu)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
        tzz(jx,jy,jz)=-(nu)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
        txz(jx,jy,jz)=-2._rprec*(nu) * S13(jx,jy,jz)
        tyz(jx,jy,jz)=-2._rprec*(nu) * S23(jx,jy,jz)
    end do
    end do
    end do
end if

#ifdef ENABLE_CUDA
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_tau_interior = sgs_t1 - sgs_t0
    sgs_t0 = sgs_t1
end if
#endif

#ifdef PPLVLSET
!--at this point tij are only set for 1:nz-1
!--at this point u, v, w are set for 0:nz, except bottom process is 1:nz
!--some MPI synchronizing may be done in here, but this will be kept
!  separate from the rest of the code (at the risk of some redundancy)
call level_set_BC ()
#endif


#ifdef PPMPI
! txz,tyz calculated for 1:nz-1 (on w-nodes) except bottom process
! (only 2:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
#ifdef ENABLE_CUDA
if ((.not. tau_halo_posted) .and. sgs_pointwise_cuda_enabled() .and.          &
    sgs_tau_prebarrier_enabled()) then
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PRE_SYNC)
    call sgs_cuda_barrier('tau before mpi')
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PRE_SYNC)
end if
if (tau_halo_posted) then
    if (tau_halo_nreq > 0) then
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_MPI)
        call mpi_waitall(tau_halo_nreq, tau_halo_req, tau_halo_status, ierr)
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_MPI)
    end if
    if (ierr /= 0) call error(sub_name,                                       &
        'Error occurred during overlapped SGS tau halo sync:', ierr)

    if (up /= MPI_PROC_NULL) then
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_UNPACK)
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            txz(jx,jy,nz) = tau_halo_recv(jx,jy,1)
            tyz(jx,jy,nz) = tau_halo_recv(jx,jy,2)
        end do
        end do
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_UNPACK)
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_POST_SYNC)
        call sgs_cuda_sync('tau overlap halo unpack')
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_POST_SYNC)
    end if
else if (sgs_pointwise_cuda_enabled() .and. nproc == 2 .and.                 &
    (sgs_halo_combined2_enabled() .or. sgs_halo_hostpinned_test_enabled())) then
    if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_enter)
    tau2_count = 3 * ld * ny
    tau2_peer = 1 - coord
    if (.not. allocated(sgs_tau2_send_d)) then
        allocate(sgs_tau2_send_d(tau2_count), sgs_tau2_recv_d(tau2_count))
    else if (size(sgs_tau2_send_d) /= tau2_count) then
        deallocate(sgs_tau2_send_d, sgs_tau2_recv_d)
        allocate(sgs_tau2_send_d(tau2_count), sgs_tau2_recv_d(tau2_count))
    end if
    if (sgs_halo_hostpinned_test_enabled()) then
        if (.not. allocated(sgs_tau2_send_h)) then
            allocate(sgs_tau2_send_h(tau2_count), sgs_tau2_recv_h(tau2_count))
        else if (size(sgs_tau2_send_h) /= tau2_count) then
            deallocate(sgs_tau2_send_h, sgs_tau2_recv_h)
            allocate(sgs_tau2_send_h(tau2_count), sgs_tau2_recv_h(tau2_count))
        end if
    end if
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PACK)
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        tau2_idx = jx + (jy - 1) * ld
        sgs_tau2_send_d(tau2_idx) = txz(jx,jy,1)
        sgs_tau2_send_d(tau2_idx + ld*ny) = tyz(jx,jy,1)
        sgs_tau2_send_d(tau2_idx + 2*ld*ny) = tzz(jx,jy,nz-1)
    end do
    end do
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PACK)
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PRE_SYNC)
    call sgs_cuda_barrier('tau nproc2 combined pack')
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PRE_SYNC)
    if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_after_sync)
    if (sgs_stage_enabled_local) call sgs_tau_detail_add_msg(tau2_count,      &
        tau2_count, tau2_peer, .true.)

    if (sgs_halo_hostpinned_test_enabled()) then
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_D2H)
        sgs_tau2_send_h = sgs_tau2_send_d
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_D2H)
        if (sgs_stage_enabled_local) call sgs_pointer_env_audit(              &
            'tau_nproc2_hostpinned', sgs_tau2_send_h, sgs_tau2_recv_h,       &
            tau2_peer)
        if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_before_mpi)
        if (sgs_stage_enabled_local .and.                                    &
            sgs_barrier_before_tau_halo_enabled()) then
            call sgs_diag_time(sgs_barrier_t0)
            call mpi_barrier(comm, ierr)
            call sgs_diag_time(sgs_barrier_t1)
            sgs_tau_arr_barrier = sgs_tau_arr_barrier +                      &
                sgs_barrier_t1 - sgs_barrier_t0
        end if
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_MPI)
        call mpi_sendrecv(sgs_tau2_send_h, tau2_count, MPI_RPREC, tau2_peer,  &
            183, sgs_tau2_recv_h, tau2_count, MPI_RPREC, tau2_peer, 183,     &
            comm, status, ierr)
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_MPI)
        if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_after_mpi)
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_H2D)
        sgs_tau2_recv_d = sgs_tau2_recv_h
        call sgs_cuda_barrier('tau nproc2 hostpinned H2D')
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_H2D)
    else
        if (sgs_stage_enabled_local) call sgs_pointer_env_audit_device(       &
            'tau_nproc2_combined_device', sgs_tau2_send_d, sgs_tau2_recv_d,  &
            tau2_peer)
        if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_before_mpi)
        if (sgs_stage_enabled_local .and.                                    &
            sgs_barrier_before_tau_halo_enabled()) then
            call sgs_diag_time(sgs_barrier_t0)
            call mpi_barrier(comm, ierr)
            call sgs_diag_time(sgs_barrier_t1)
            sgs_tau_arr_barrier = sgs_tau_arr_barrier +                      &
                sgs_barrier_t1 - sgs_barrier_t0
        end if
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_MPI)
        call mpi_sendrecv(sgs_tau2_send_d, tau2_count, MPI_RPREC, tau2_peer,  &
            183, sgs_tau2_recv_d, tau2_count, MPI_RPREC, tau2_peer, 183,     &
            comm, status, ierr)
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_MPI)
        if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_after_mpi)
    end if
    if (ierr /= 0) call error(sub_name,                                       &
        'Error occurred during nproc2 combined SGS tau halo sync:', ierr)

    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_UNPACK)
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        tau2_idx = jx + (jy - 1) * ld
        if (coord == 0) then
            txz(jx,jy,nz) = sgs_tau2_recv_d(tau2_idx)
            tyz(jx,jy,nz) = sgs_tau2_recv_d(tau2_idx + ld*ny)
        else
            tzz(jx,jy,0) = sgs_tau2_recv_d(tau2_idx + 2*ld*ny)
        end if
    end do
    end do
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_UNPACK)
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_POST_SYNC)
    call sgs_cuda_sync('tau nproc2 combined unpack')
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_POST_SYNC)
    if (sgs_stage_enabled_local) call sgs_diag_time(sgs_tau_arr_after_post)
else if (sgs_combined_tau_halo_active()) then
    if (.not. allocated(tau_halo_send)) then
        allocate(tau_halo_send(ld, ny, 3), tau_halo_recv(ld, ny, 3))
    else if (size(tau_halo_send, 1) /= ld .or. size(tau_halo_send, 3) /= 3) then
        deallocate(tau_halo_send, tau_halo_recv)
        allocate(tau_halo_send(ld, ny, 3), tau_halo_recv(ld, ny, 3))
    end if
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PACK)
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        tau_halo_send(jx,jy,1) = txz(jx,jy,1)
        tau_halo_send(jx,jy,2) = tyz(jx,jy,1)
        tau_halo_send(jx,jy,3) = tzz(jx,jy,nz-1)
    end do
    end do
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PACK)
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PRE_SYNC)
    call sgs_cuda_barrier('tau combined halo pack')
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PRE_SYNC)
    if (sgs_stage_enabled_local) call sgs_pointer_env_audit(                  &
        'tau_old_combined_managed', tau_halo_send(1,1,1),                    &
        tau_halo_recv(1,1,1), merge(up, down, coord < nproc - 1))

    tau_halo_req = MPI_REQUEST_NULL
    tau_halo_nreq = 0
    if (coord < nproc - 1) then
        if (sgs_stage_enabled_local) then
            call sgs_tau_detail_add_msg(0, 2*ld*ny, up, .false.)
            call sgs_tau_detail_add_msg(ld*ny, 0, up, .false.)
        end if
        tau_halo_nreq = tau_halo_nreq + 1
        call mpi_irecv(tau_halo_recv(1,1,1), 2*ld*ny, MPI_RPREC, up, 81,      &
            comm, tau_halo_req(tau_halo_nreq), ierr)
        tau_halo_nreq = tau_halo_nreq + 1
        call mpi_isend(tau_halo_send(1,1,3), ld*ny, MPI_RPREC, up, 82, comm,  &
            tau_halo_req(tau_halo_nreq), ierr)
    end if
    if (coord > 0) then
        if (sgs_stage_enabled_local) then
            call sgs_tau_detail_add_msg(0, ld*ny, down, .false.)
            call sgs_tau_detail_add_msg(2*ld*ny, 0, down, .false.)
        end if
        tau_halo_nreq = tau_halo_nreq + 1
        call mpi_irecv(tzz(1,1,0), ld*ny, MPI_RPREC, down, 82, comm,          &
            tau_halo_req(tau_halo_nreq), ierr)
        tau_halo_nreq = tau_halo_nreq + 1
        call mpi_isend(tau_halo_send(1,1,1), 2*ld*ny, MPI_RPREC, down, 81,    &
            comm, tau_halo_req(tau_halo_nreq), ierr)
    end if
    if (tau_halo_nreq > 0) then
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_MPI)
        call mpi_waitall(tau_halo_nreq, tau_halo_req, tau_halo_status, ierr)
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_MPI)
    end if
    if (ierr /= 0) call error(sub_name,                                       &
        'Error occurred during combined SGS tau halo sync:', ierr)

    if (up /= MPI_PROC_NULL) then
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_UNPACK)
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            txz(jx,jy,nz) = tau_halo_recv(jx,jy,1)
            tyz(jx,jy,nz) = tau_halo_recv(jx,jy,2)
        end do
        end do
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_UNPACK)
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_POST_SYNC)
        call sgs_cuda_sync('tau combined halo unpack')
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_POST_SYNC)
    end if
else if (sgs_pointwise_cuda_enabled() .and. sgs_packed_tau_halo_enabled()) then
    if (.not. allocated(tau_halo_send)) then
        allocate(tau_halo_send(nx, ny, 2), tau_halo_recv(nx, ny, 2))
    else if (size(tau_halo_send, 3) /= 2) then
        deallocate(tau_halo_send, tau_halo_recv)
        allocate(tau_halo_send(nx, ny, 2), tau_halo_recv(nx, ny, 2))
    end if
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PACK)
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, nx
        tau_halo_send(jx,jy,1) = txz(jx,jy,1)
        tau_halo_send(jx,jy,2) = tyz(jx,jy,1)
    end do
    end do
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PACK)
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_PRE_SYNC)
    call sgs_cuda_barrier('tau packed halo pack')
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_PRE_SYNC)
    if (sgs_stage_enabled_local) call sgs_pointer_env_audit(                  &
        'tau_packed_managed', tau_halo_send(1,1,1), tau_halo_recv(1,1,1),    &
        merge(up, down, coord < nproc - 1))
    if (sgs_stage_enabled_local) then
        if (down /= MPI_PROC_NULL) call sgs_tau_detail_add_msg(2*nx*ny, 0,   &
            down, .true.)
        if (up /= MPI_PROC_NULL) call sgs_tau_detail_add_msg(0, 2*nx*ny, up, &
            .true.)
        call sgs_tau_detail_start(SGS_TAU_MPI)
    end if
    call mpi_dbg_sendrecv_r(tau_halo_send(1,1,1), 2*nx*ny, MPI_RPREC, down, 71, &
        tau_halo_recv(1,1,1), 2*nx*ny, MPI_RPREC, up, 71, comm, status, ierr, &
        'sgs_tau_packed_down')
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_MPI)
    if (ierr /= 0) call error(sub_name,                                       &
        'Error occurred during packed SGS tau halo sync:', ierr)
    if (up /= MPI_PROC_NULL) then
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_UNPACK)
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, nx
            txz(jx,jy,nz) = tau_halo_recv(jx,jy,1)
            tyz(jx,jy,nz) = tau_halo_recv(jx,jy,2)
        end do
        end do
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_UNPACK)
        if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_POST_SYNC)
        call sgs_cuda_sync('tau packed halo unpack')
        if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_POST_SYNC)
    end if
else
#endif
#ifdef ENABLE_CUDA
    if (sgs_stage_enabled_local) call sgs_tau_detail_start(SGS_TAU_MPI)
#endif
    call mpi_sync_real_array( txz, 0, MPI_SYNC_DOWN )
    call mpi_sync_real_array( tyz, 0, MPI_SYNC_DOWN )
#ifdef ENABLE_CUDA
    if (sgs_stage_enabled_local) call sgs_tau_detail_stop(SGS_TAU_MPI)
#endif
#ifdef ENABLE_CUDA
end if
#endif
#ifdef ENABLE_CUDA
if (sgs_stage_enabled_local) then
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_tau_halo = sgs_t1 - sgs_t0
    sgs_t0 = sgs_t1
end if
#endif
#ifdef PPSAFETYMODE
! Set bogus values (easier to catch if there's an error)
#ifdef ENABLE_CUDA
if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, nx
        txx(jx,jy,0) = BOGUS
        txy(jx,jy,0) = BOGUS
        txz(jx,jy,0) = BOGUS
        tyy(jx,jy,0) = BOGUS
        tyz(jx,jy,0) = BOGUS
        tzz(jx,jy,0) = BOGUS
    end do
    end do
    call sgs_cuda_sync('tau lower bogus plane')
else
#endif
txx(:, :, 0) = BOGUS
txy(:, :, 0) = BOGUS
txz(:, :, 0) = BOGUS
tyy(:, :, 0) = BOGUS
tyz(:, :, 0) = BOGUS
tzz(:, :, 0) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif
#endif

! Set bogus values (easier to catch if there's an error)
#ifdef PPSAFETYMODE
#ifdef ENABLE_CUDA
if (sgs_pointwise_cuda_enabled()) then
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, nx
        txx(jx,jy,nz) = BOGUS
        txy(jx,jy,nz) = BOGUS
        tyy(jx,jy,nz) = BOGUS
        tzz(jx,jy,nz) = BOGUS
    end do
    end do
    call sgs_cuda_sync('tau upper bogus plane')
else
#endif
txx(:, :, nz) = BOGUS
txy(:, :, nz) = BOGUS
tyy(:, :, nz) = BOGUS
tzz(:, :, nz) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif

#ifdef ENABLE_CUDA
if (sgs_pointwise_cuda_enabled()) call sgs_cuda_barrier('sgs_stag final')
if (sgs_stage_enabled_local) then
    if (sgs_event_active) then
        call sgs_event_elapsed_seconds(sgs_evt_start(SGS_EVT_CALC),           &
            sgs_evt_stop(SGS_EVT_CALC), sgs_gpu_calc, 'calc_Sij')
        call sgs_event_elapsed_seconds(sgs_evt_start(SGS_EVT_DYNAMIC),        &
            sgs_evt_stop(SGS_EVT_DYNAMIC), sgs_gpu_dynamic, 'dynamic model')
        call sgs_event_elapsed_seconds(sgs_evt_start(SGS_EVT_NUT),            &
            sgs_evt_stop(SGS_EVT_NUT), sgs_gpu_nut, 'Nu_t')
        call sgs_event_elapsed_seconds(sgs_evt_start(SGS_EVT_TAU_BOUNDARY),   &
            sgs_evt_stop(SGS_EVT_TAU_BOUNDARY), sgs_gpu_tau_boundary,         &
            'tau boundary')
        call sgs_event_elapsed_seconds(sgs_evt_start(SGS_EVT_TAU_INTERIOR),   &
            sgs_evt_stop(SGS_EVT_TAU_INTERIOR), sgs_gpu_tau_interior,         &
            'tau interior')
    end if
#ifdef PPMPI
    sgs_t1 = mpi_wtime()
#else
    call cpu_time(sgs_t1)
#endif
    sgs_stage_final = sgs_t1 - sgs_t0
    sgs_stage_count = sgs_stage_count + 1
    call sgs_stage_report(sgs_stage_count, sgs_stage_calc,                    &
        sgs_stage_dynamic, sgs_stage_nut, sgs_stage_tau_boundary,             &
        sgs_stage_tau_interior, sgs_stage_tau_halo, sgs_stage_final,          &
        sgs_gpu_calc, sgs_gpu_dynamic, sgs_gpu_nut, sgs_gpu_tau_boundary,     &
        sgs_gpu_tau_interior)
end if
if (sgs_event_active) then
    do sgs_event_i = 1, SGS_EVT_COUNT
        if (cudaEventDestroy(sgs_evt_start(sgs_event_i)) /= 0) then
            print *, 'sgs_stag CUDA event destroy failure'
            stop
        end if
        if (cudaEventDestroy(sgs_evt_stop(sgs_event_i)) /= 0) then
            print *, 'sgs_stag CUDA event destroy failure'
            stop
        end if
    end do
end if
#endif

end subroutine sgs_stag

#ifdef ENABLE_CUDA
!*******************************************************************************
attributes(global) subroutine calc_sij_interior_kernel(dudx_arr, dudy_arr,    &
    dudz_arr, dvdx_arr, dvdy_arr, dvdz_arr, dwdx_arr, dwdy_arr, dwdz_arr,     &
    s11_arr, s12_arr, s13_arr, s22_arr, s23_arr, s33_arr, ldv, nxv, nyv,      &
    lbzv, jminv, jmaxv)
!*******************************************************************************
use types, only : rprec
implicit none

integer, value :: ldv, nxv, nyv, lbzv, jminv, jmaxv
real(rprec), device :: dudx_arr(ldv,nyv,*), dudy_arr(ldv,nyv,*)
real(rprec), device :: dudz_arr(ldv,nyv,*), dvdx_arr(ldv,nyv,*)
real(rprec), device :: dvdy_arr(ldv,nyv,*), dvdz_arr(ldv,nyv,*)
real(rprec), device :: dwdx_arr(ldv,nyv,*), dwdy_arr(ldv,nyv,*)
real(rprec), device :: dwdz_arr(ldv,nyv,*)
real(rprec), device :: s11_arr(ldv,nyv,*), s12_arr(ldv,nyv,*)
real(rprec), device :: s13_arr(ldv,nyv,*), s22_arr(ldv,nyv,*)
real(rprec), device :: s23_arr(ldv,nyv,*), s33_arr(ldv,nyv,*)
integer :: jx, jy, jz, kp, km
real(rprec) :: uy, vx

jx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jy = (blockIdx%y - 1) * blockDim%y + threadIdx%y
jz = jminv + (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1
if (jx > nxv .or. jy > nyv .or. jz > jmaxv) return

kp = jz - lbzv + 1
km = jz - lbzv

s11_arr(jx,jy,jz) = 0.5_rprec*(dudx_arr(jx,jy,kp) +                         &
    dudx_arr(jx,jy,km))
uy = dudy_arr(jx,jy,kp) + dudy_arr(jx,jy,km)
vx = dvdx_arr(jx,jy,kp) + dvdx_arr(jx,jy,km)
s12_arr(jx,jy,jz) = 0.25_rprec*(uy+vx)
s13_arr(jx,jy,jz) = 0.5_rprec*(dudz_arr(jx,jy,kp) + dwdx_arr(jx,jy,kp))
s22_arr(jx,jy,jz) = 0.5_rprec*(dvdy_arr(jx,jy,kp) +                         &
    dvdy_arr(jx,jy,km))
s23_arr(jx,jy,jz) = 0.5_rprec*(dvdz_arr(jx,jy,kp) + dwdy_arr(jx,jy,kp))
s33_arr(jx,jy,jz) = 0.5_rprec*(dwdz_arr(jx,jy,kp) +                         &
    dwdz_arr(jx,jy,km))

end subroutine calc_sij_interior_kernel

!*******************************************************************************
subroutine sgs_calc_sij_device_lower_bound_bench(jmin, jmax)
!*******************************************************************************
use types, only : rprec
use param, only : ld, nx, ny, nz, lbz, coord
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
implicit none

integer, intent(in) :: jmin, jmax
integer, parameter :: bench_iters = 3
real(rprec), device, allocatable :: dudx_d(:,:,:), dudy_d(:,:,:)
real(rprec), device, allocatable :: dudz_d(:,:,:), dvdx_d(:,:,:)
real(rprec), device, allocatable :: dvdy_d(:,:,:), dvdz_d(:,:,:)
real(rprec), device, allocatable :: dwdx_d(:,:,:), dwdy_d(:,:,:)
real(rprec), device, allocatable :: dwdz_d(:,:,:)
real(rprec), device, allocatable :: s11_d(:,:,:), s12_d(:,:,:), s13_d(:,:,:)
real(rprec), device, allocatable :: s22_d(:,:,:), s23_d(:,:,:), s33_d(:,:,:)
type(cudaEvent) :: evt_start, evt_stop
type(dim3) :: block, grid
integer :: block_x, block_y, block_z, grid_x, grid_y, grid_z
integer :: istat, iter
real :: elapsed_ms
real(rprec) :: t0, t1, setup_wall, kernel_event, kernel_wall
integer(8) :: cells, input_bytes, output_bytes
logical, save :: done = .false.

if (done .or. jmin > jmax) return
done = .true.

call sgs_diag_time(t0)
allocate(dudx_d(ld,ny,lbz:nz), dudy_d(ld,ny,lbz:nz), dudz_d(ld,ny,lbz:nz))
allocate(dvdx_d(ld,ny,lbz:nz), dvdy_d(ld,ny,lbz:nz), dvdz_d(ld,ny,lbz:nz))
allocate(dwdx_d(ld,ny,lbz:nz), dwdy_d(ld,ny,lbz:nz), dwdz_d(ld,ny,lbz:nz))
allocate(s11_d(ld,ny,nz), s12_d(ld,ny,nz), s13_d(ld,ny,nz))
allocate(s22_d(ld,ny,nz), s23_d(ld,ny,nz), s33_d(ld,ny,nz))

dudx_d = dudx
dudy_d = dudy
dudz_d = dudz
dvdx_d = dvdx
dvdy_d = dvdy
dvdz_d = dvdz
dwdx_d = dwdx
dwdy_d = dwdy
dwdz_d = dwdz
istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'SGS calc_Sij device bench copy/sync failure:', istat
    stop
end if
call sgs_diag_time(t1)
setup_wall = t1 - t0

block_x = 32
block_y = 4
block_z = 2
grid_x = (nx + block_x - 1) / block_x
grid_y = (ny + block_y - 1) / block_y
grid_z = (jmax - jmin + 1 + block_z - 1) / block_z
block = dim3(block_x, block_y, block_z)
grid = dim3(grid_x, grid_y, grid_z)

istat = cudaEventCreate(evt_start)
if (istat /= 0) stop 'SGS calc_Sij device bench event create failed'
istat = cudaEventCreate(evt_stop)
if (istat /= 0) stop 'SGS calc_Sij device bench event create failed'

call sgs_diag_time(t0)
istat = cudaEventRecord(evt_start, 0)
if (istat /= 0) stop 'SGS calc_Sij device bench event record failed'
do iter = 1, bench_iters
    call calc_sij_interior_kernel<<<grid,block>>>(dudx_d, dudy_d, dudz_d,    &
        dvdx_d, dvdy_d, dvdz_d, dwdx_d, dwdy_d, dwdz_d, s11_d, s12_d,       &
        s13_d, s22_d, s23_d, s33_d, ld, nx, ny, lbz, jmin, jmax)
end do
istat = cudaEventRecord(evt_stop, 0)
if (istat /= 0) stop 'SGS calc_Sij device bench event record failed'
istat = cudaEventSynchronize(evt_stop)
if (istat /= 0) stop 'SGS calc_Sij device bench event sync failed'
call sgs_diag_time(t1)
kernel_wall = (t1 - t0) / real(bench_iters, rprec)
istat = cudaEventElapsedTime(elapsed_ms, evt_start, evt_stop)
if (istat /= 0) stop 'SGS calc_Sij device bench event elapsed failed'
kernel_event = real(elapsed_ms, rprec) * 1.0e-3_rprec /                     &
    real(bench_iters, rprec)

istat = cudaEventDestroy(evt_start)
if (istat /= 0) stop 'SGS calc_Sij device bench event destroy failed'
istat = cudaEventDestroy(evt_stop)
if (istat /= 0) stop 'SGS calc_Sij device bench event destroy failed'

cells = int(nx,8) * int(ny,8) * int(jmax - jmin + 1,8)
input_bytes = 9_8 * int(ld,8) * int(ny,8) * int(nz - lbz + 1,8) *            &
    SGS_RPREC_BYTES
output_bytes = 6_8 * int(ld,8) * int(ny,8) * int(nz,8) * SGS_RPREC_BYTES

if (coord == 0) then
    write(*,'(a)') 'SGS_CALCSIJ_DEVICE_LOWER_BOUND_BENCH'
    write(*,'(a,i8,a,i8,a,i8)') '  grid nx,ny,nz=', nx, ',', ny, ',', nz
    write(*,'(a,i8,a,i8,a,i12)') '  z range=', jmin, ':', jmax,              &
        ' cells=', cells
    write(*,'(a,E15.7)') '  setup copy+alloc wall: ', setup_wall
    write(*,'(a,E15.7)') '  kernel GPU event avg: ', kernel_event
    write(*,'(a,E15.7)') '  kernel CPU wall avg: ', kernel_wall
    write(*,'(a,i4)') '  iterations: ', bench_iters
    write(*,'(a,i12,a,i12)') '  input_bytes=', input_bytes,                  &
        ' output_bytes=', output_bytes
    write(*,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,i5)') '  block/grid=', block_x,     &
        'x', block_y, 'x', block_z, ' / ', grid_x, 'x', grid_y, 'x', grid_z
end if

deallocate(dudx_d, dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d)
deallocate(dwdx_d, dwdy_d, dwdz_d)
deallocate(s11_d, s12_d, s13_d, s22_d, s23_d, s33_d)
istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'SGS calc_Sij device bench deallocate/sync failure:', istat
    stop
end if

end subroutine sgs_calc_sij_device_lower_bound_bench

!*******************************************************************************
attributes(global) subroutine sgs_nut_dynamic_kernel(s11_arr, s12_arr,        &
    s13_arr, s22_arr, s23_arr, s33_arr, cs_arr, nut_arr, smag_top_arr,        &
    l_delta2v, ldv, nxv, nyv, nzv)
!*******************************************************************************
use types, only : rprec
implicit none

integer, value :: ldv, nxv, nyv, nzv
real(rprec), value :: l_delta2v
real(rprec), device :: s11_arr(ldv,nyv,*), s12_arr(ldv,nyv,*)
real(rprec), device :: s13_arr(ldv,nyv,*), s22_arr(ldv,nyv,*)
real(rprec), device :: s23_arr(ldv,nyv,*), s33_arr(ldv,nyv,*)
real(rprec), device :: cs_arr(ldv,nyv,*), nut_arr(ldv,nyv,*)
real(rprec), device :: smag_top_arr(ldv,nyv)
integer :: jx, jy, jz
real(rprec) :: s11v, s12v, s13v, s22v, s23v, s33v, smag

jx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jy = (blockIdx%y - 1) * blockDim%y + threadIdx%y
jz = (blockIdx%z - 1) * blockDim%z + threadIdx%z
if (jx > nxv .or. jy > nyv .or. jz > nzv) return

s11v = s11_arr(jx,jy,jz)
s12v = s12_arr(jx,jy,jz)
s13v = s13_arr(jx,jy,jz)
s22v = s22_arr(jx,jy,jz)
s23v = s23_arr(jx,jy,jz)
s33v = s33_arr(jx,jy,jz)
smag = sqrt(2._rprec * (s11v*s11v + s22v*s22v + s33v*s33v +                 &
    2._rprec * (s12v*s12v + s13v*s13v + s23v*s23v)))
nut_arr(jx,jy,jz) = smag * cs_arr(jx,jy,jz) * l_delta2v
if (jz == nzv) smag_top_arr(jx,jy) = smag

end subroutine sgs_nut_dynamic_kernel

!*******************************************************************************
attributes(global) subroutine sgs_tau_interior_kernel(nut_arr, s11_arr,       &
    s12_arr, s13_arr, s22_arr, s23_arr, s33_arr, txx_arr, txy_arr, txz_arr,   &
    tyy_arr, tyz_arr, tzz_arr, nuv, ldv, nxv, nyv, jminv, jmaxv)
!*******************************************************************************
use types, only : rprec
implicit none

integer, value :: ldv, nxv, nyv, jminv, jmaxv
real(rprec), value :: nuv
real(rprec), device :: nut_arr(ldv,nyv,*)
real(rprec), device :: s11_arr(ldv,nyv,*), s12_arr(ldv,nyv,*)
real(rprec), device :: s13_arr(ldv,nyv,*), s22_arr(ldv,nyv,*)
real(rprec), device :: s23_arr(ldv,nyv,*), s33_arr(ldv,nyv,*)
real(rprec), device :: txx_arr(ldv,nyv,0:*), txy_arr(ldv,nyv,0:*)
real(rprec), device :: txz_arr(ldv,nyv,0:*), tyy_arr(ldv,nyv,0:*)
real(rprec), device :: tyz_arr(ldv,nyv,0:*), tzz_arr(ldv,nyv,0:*)
integer :: jx, jy, jz
real(rprec) :: const, const2, const3, const4

jx = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jy = (blockIdx%y - 1) * blockDim%y + threadIdx%y
jz = jminv + (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1
if (jx > nxv .or. jy > nyv .or. jz > jmaxv) return

const3 = -nuv
const4 = -2._rprec * nuv
const = -0.5_rprec * (nut_arr(jx,jy,jz) + nut_arr(jx,jy,jz+1))
const2 = -2._rprec * nut_arr(jx,jy,jz)

txx_arr(jx,jy,jz) = (const + const3) *                                      &
    (s11_arr(jx,jy,jz) + s11_arr(jx,jy,jz+1))
txy_arr(jx,jy,jz) = (const + const3) *                                      &
    (s12_arr(jx,jy,jz) + s12_arr(jx,jy,jz+1))
tyy_arr(jx,jy,jz) = (const + const3) *                                      &
    (s22_arr(jx,jy,jz) + s22_arr(jx,jy,jz+1))
tzz_arr(jx,jy,jz) = (const + const3) *                                      &
    (s33_arr(jx,jy,jz) + s33_arr(jx,jy,jz+1))
txz_arr(jx,jy,jz) = (const2 + const4) * s13_arr(jx,jy,jz)
tyz_arr(jx,jy,jz) = (const2 + const4) * s23_arr(jx,jy,jz)

end subroutine sgs_tau_interior_kernel
#endif

!*******************************************************************************
subroutine calc_Sij()
!*******************************************************************************
! Calculate the resolved strain rate tensor, Sij = 0.5(djui - diuj)
!   values are stored on w-nodes (1:nz)

use types, only : rprec
use param
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sgs_param
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
#if defined(ENABLE_CUDA) && defined(PPMPI)
use mpi
#endif
implicit none

real(rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
integer :: jx, jy, jz
integer :: jz_min, jz_max
#ifdef ENABLE_CUDA
logical :: calc_diag_active
logical :: calc_sij_explicit
logical :: calc_sij_explicit_used
integer :: sij_bulk_min, sij_exp_bulk_min, sij_exp_bulk_max
integer :: sij_block_x, sij_block_y, sij_block_z
integer :: sij_grid_x, sij_grid_y, sij_grid_z
type(dim3) :: sij_block, sij_grid
#endif
#if defined(ENABLE_CUDA) && defined(PPMPI)
integer :: dwdz_req(2)
integer :: dwdz_status(MPI_STATUS_SIZE, 2)
integer :: dwdz_nreq
integer :: sij_bulk_max
integer :: dwdz_count, dwdz_peer, dwdz_send_count, dwdz_recv_count
integer :: dwdz_idx
logical :: dwdz_halo_posted, dwdz_device_halo_posted
logical :: dwdz_hostpinned_halo_posted
logical :: dwdz_symmetric_halo_posted
real(rprec) :: dwdz_barrier_t0, dwdz_barrier_t1
#endif

#ifdef ENABLE_CUDA
calc_diag_active = sgs_stage_timing_enabled() .and. sgs_pointwise_cuda_enabled()
calc_sij_explicit = sgs_calc_sij_explicit_enabled()
calc_sij_explicit_used = .false.
sij_bulk_min = 0
sij_exp_bulk_min = 0
sij_exp_bulk_max = 0
sij_block_x = 0
sij_block_y = 0
sij_block_z = 0
sij_grid_x = 0
sij_grid_y = 0
sij_grid_z = 0
if (calc_diag_active) call sgs_calc_diag_begin()
#endif

#if defined(ENABLE_CUDA) && defined(PPMPI)
dwdz_halo_posted = .false.
dwdz_device_halo_posted = .false.
dwdz_hostpinned_halo_posted = .false.
dwdz_symmetric_halo_posted = .false.
dwdz_nreq = 0
if (sgs_pointwise_cuda_enabled() .and. sgs_overlap_dwdz_halo_enabled() .and.  &
    (nproc > 1)) then
    if (calc_diag_active) call sgs_diag_time(sgs_dwdz_arr_enter)
    if (calc_diag_active) call sgs_calc_cpu_start(SGS_CALC_DWDZ_PRE_SYNC)
    if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_PRE_SYNC)
    call sgs_cuda_barrier('Sij dwdz overlap before mpi')
    if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_PRE_SYNC)
    if (calc_diag_active) call sgs_calc_cpu_stop(SGS_CALC_DWDZ_PRE_SYNC)
    if (calc_diag_active) then
        call sgs_diag_time(sgs_dwdz_arr_after_sync)
        sgs_dwdz_arr_before_mpi = sgs_dwdz_arr_after_sync
    end if
    dwdz_req = MPI_REQUEST_NULL
    dwdz_count = ld * ny
    if (nproc == 2 .and. (sgs_dwdz_device_halo_enabled() .or.                &
        sgs_dwdz_hostpinned_test_enabled() .or.                              &
        sgs_dwdz_symmetric_test_enabled())) then
        dwdz_peer = 1 - coord
        dwdz_send_count = merge(dwdz_count, 0, coord > 0)
        dwdz_recv_count = merge(dwdz_count, 0, coord < nproc - 1)
        if (.not. allocated(sgs_dwdz_send_d)) then
            allocate(sgs_dwdz_send_d(dwdz_count), sgs_dwdz_recv_d(dwdz_count))
        else if (size(sgs_dwdz_send_d) /= dwdz_count) then
            deallocate(sgs_dwdz_send_d, sgs_dwdz_recv_d)
            allocate(sgs_dwdz_send_d(dwdz_count), sgs_dwdz_recv_d(dwdz_count))
        end if
        if (sgs_dwdz_hostpinned_test_enabled()) then
            if (.not. allocated(sgs_dwdz_send_h)) then
                allocate(sgs_dwdz_send_h(dwdz_count),                       &
                    sgs_dwdz_recv_h(dwdz_count))
            else if (size(sgs_dwdz_send_h) /= dwdz_count) then
                deallocate(sgs_dwdz_send_h, sgs_dwdz_recv_h)
                allocate(sgs_dwdz_send_h(dwdz_count),                       &
                    sgs_dwdz_recv_h(dwdz_count))
            end if
        end if
        if (sgs_dwdz_symmetric_test_enabled()) then
            if (.not. allocated(sgs_dwdz_dummy_send_d)) then
                allocate(sgs_dwdz_dummy_send_d(dwdz_count),                  &
                    sgs_dwdz_dummy_recv_d(dwdz_count))
            else if (size(sgs_dwdz_dummy_send_d) /= dwdz_count) then
                deallocate(sgs_dwdz_dummy_send_d, sgs_dwdz_dummy_recv_d)
                allocate(sgs_dwdz_dummy_send_d(dwdz_count),                  &
                    sgs_dwdz_dummy_recv_d(dwdz_count))
            end if
        end if
        if (dwdz_send_count > 0) then
            if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_PACK)
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, ld
                dwdz_idx = jx + (jy - 1) * ld
                sgs_dwdz_send_d(dwdz_idx) = dwdz(jx,jy,1)
            end do
            end do
            call sgs_cuda_barrier('dwdz device halo pack')
            if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_PACK)
        end if
        if (sgs_dwdz_symmetric_test_enabled() .and. coord == 0) then
            if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_DUMMY)
!$cuf kernel do(1) <<<*,*>>>
            do dwdz_idx = 1, dwdz_count
                sgs_dwdz_dummy_send_d(dwdz_idx) = 0._rprec
            end do
            call sgs_cuda_barrier('dwdz symmetric dummy init')
            if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_DUMMY)
        end if
        if (sgs_dwdz_hostpinned_test_enabled()) then
            if (dwdz_send_count > 0) then
                if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_D2H)
                sgs_dwdz_send_h(1:dwdz_send_count) =                         &
                    sgs_dwdz_send_d(1:dwdz_send_count)
                if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_D2H)
            end if
            if (calc_diag_active) then
                call sgs_pointer_env_audit('dwdz_hostpinned',                &
                    sgs_dwdz_send_h, sgs_dwdz_recv_h, dwdz_peer)
                call sgs_dwdz_detail_add_msg(dwdz_send_count,                &
                    dwdz_recv_count, dwdz_peer, .true.)
                call sgs_dwdz_path_audit('hostpinned',                       &
                    'GPU_pack_D2H_Isend_Irecv_Wait_H2D_GPU_unpack',          &
                    dwdz_send_count, dwdz_recv_count, dwdz_peer, .false.,    &
                    .false.)
            end if
        else if (sgs_dwdz_symmetric_test_enabled()) then
            if (calc_diag_active) then
                sgs_dwdz_dummy_bytes = int(dwdz_count, 8) * SGS_RPREC_BYTES
                if (coord == 0) then
                    call sgs_pointer_env_audit_device(                       &
                        'dwdz_symmetric_device', sgs_dwdz_dummy_send_d,      &
                        sgs_dwdz_recv_d, dwdz_peer)
                else
                    call sgs_pointer_env_audit_device(                       &
                        'dwdz_symmetric_device', sgs_dwdz_send_d,            &
                        sgs_dwdz_dummy_recv_d, dwdz_peer)
                end if
                call sgs_dwdz_detail_add_msg(dwdz_send_count,                &
                    dwdz_recv_count, dwdz_peer, .true.)
                call sgs_dwdz_path_audit('symmetric_device',                 &
                    'GPU_pack_Isend_Irecv_symmetric_dummy_Wait_GPU_unpack',   &
                    dwdz_send_count, dwdz_recv_count, dwdz_peer, .false.,    &
                    .false.)
            end if
        else
            if (calc_diag_active) then
                call sgs_pointer_env_audit_device('dwdz_device_halo',        &
                    sgs_dwdz_send_d, sgs_dwdz_recv_d, dwdz_peer)
                call sgs_dwdz_detail_add_msg(dwdz_send_count,                &
                    dwdz_recv_count, dwdz_peer, .true.)
                call sgs_dwdz_path_audit('device_halo',                      &
                    'GPU_pack_Isend_Irecv_Wait_GPU_unpack',                  &
                    dwdz_send_count, dwdz_recv_count, dwdz_peer, .false.,    &
                    .false.)
            end if
        end if
        if (calc_diag_active .and. sgs_barrier_before_dwdz_enabled()) then
            call sgs_diag_time(dwdz_barrier_t0)
            call mpi_barrier(comm, ierr)
            call sgs_diag_time(dwdz_barrier_t1)
            sgs_dwdz_arr_barrier = sgs_dwdz_arr_barrier +                    &
                dwdz_barrier_t1 - dwdz_barrier_t0
        end if
        if (calc_diag_active) call sgs_diag_time(sgs_dwdz_arr_before_mpi)
        if (sgs_dwdz_symmetric_test_enabled()) then
            dwdz_nreq = dwdz_nreq + 1
            if (coord == 0) then
                call mpi_irecv(sgs_dwdz_recv_d, dwdz_count, MPI_RPREC,       &
                    dwdz_peer, 194, comm, dwdz_req(dwdz_nreq), ierr)
            else
                call mpi_irecv(sgs_dwdz_dummy_recv_d, dwdz_count,            &
                    MPI_RPREC, dwdz_peer, 195, comm, dwdz_req(dwdz_nreq),    &
                    ierr)
            end if
            dwdz_nreq = dwdz_nreq + 1
            if (coord == 0) then
                call mpi_isend(sgs_dwdz_dummy_send_d, dwdz_count,            &
                    MPI_RPREC, dwdz_peer, 195, comm, dwdz_req(dwdz_nreq),    &
                    ierr)
            else
                call mpi_isend(sgs_dwdz_send_d, dwdz_count, MPI_RPREC,       &
                    dwdz_peer, 194, comm, dwdz_req(dwdz_nreq), ierr)
            end if
        else
            if (dwdz_recv_count > 0) then
                dwdz_nreq = dwdz_nreq + 1
                if (sgs_dwdz_hostpinned_test_enabled()) then
                    call mpi_irecv(sgs_dwdz_recv_h, dwdz_recv_count,         &
                        MPI_RPREC, dwdz_peer, 193, comm,                    &
                        dwdz_req(dwdz_nreq), ierr)
                else
                    call mpi_irecv(sgs_dwdz_recv_d, dwdz_recv_count,         &
                        MPI_RPREC, dwdz_peer, 193, comm,                    &
                        dwdz_req(dwdz_nreq), ierr)
                end if
            end if
            if (dwdz_send_count > 0) then
                dwdz_nreq = dwdz_nreq + 1
                if (sgs_dwdz_hostpinned_test_enabled()) then
                    call mpi_isend(sgs_dwdz_send_h, dwdz_send_count,         &
                        MPI_RPREC, dwdz_peer, 193, comm,                    &
                        dwdz_req(dwdz_nreq), ierr)
                else
                    call mpi_isend(sgs_dwdz_send_d, dwdz_send_count,         &
                        MPI_RPREC, dwdz_peer, 193, comm,                    &
                        dwdz_req(dwdz_nreq), ierr)
                end if
            end if
        end if
        dwdz_symmetric_halo_posted = sgs_dwdz_symmetric_test_enabled()
        dwdz_device_halo_posted = .not. sgs_dwdz_hostpinned_test_enabled()
        dwdz_hostpinned_halo_posted = sgs_dwdz_hostpinned_test_enabled()
    else
        if (calc_diag_active) call sgs_pointer_env_audit('dwdz_direct',       &
            dwdz(1,1,1), dwdz(1,1,nz), merge(up, down, coord < nproc - 1))
        dwdz_peer = merge(up, down, coord < nproc - 1)
        if (calc_diag_active) then
            call sgs_dwdz_path_audit('old_direct',                           &
                'Isend_Irecv_Wait_direct_managed_plane',                     &
                merge(ld*ny, 0, coord > 0), merge(ld*ny, 0,                  &
                coord < nproc - 1), dwdz_peer, .true., .false.)
            if (coord < nproc - 1) call sgs_dwdz_detail_add_msg(0,           &
                ld*ny, up, .false.)
            if (coord > 0) call sgs_dwdz_detail_add_msg(ld*ny, 0, down,      &
                .false.)
        end if
        if (calc_diag_active .and. sgs_barrier_before_dwdz_enabled()) then
            call sgs_diag_time(dwdz_barrier_t0)
            call mpi_barrier(comm, ierr)
            call sgs_diag_time(dwdz_barrier_t1)
            sgs_dwdz_arr_barrier = sgs_dwdz_arr_barrier +                    &
                dwdz_barrier_t1 - dwdz_barrier_t0
        end if
        if (calc_diag_active) call sgs_diag_time(sgs_dwdz_arr_before_mpi)
        if (coord < nproc - 1) then
            dwdz_nreq = dwdz_nreq + 1
            call mpi_irecv(dwdz(1,1,nz), ld*ny, MPI_RPREC, up, 93, comm,     &
                dwdz_req(dwdz_nreq), ierr)
        end if
        if (coord > 0) then
            dwdz_nreq = dwdz_nreq + 1
            call mpi_isend(dwdz(1,1,1), ld*ny, MPI_RPREC, down, 93, comm,    &
                dwdz_req(dwdz_nreq), ierr)
        end if
    end if
    dwdz_halo_posted = .true.
end if
#endif


! Calculate Sij for jz=1 (coord==0 only)
!   stored on uvp-nodes (this level only) for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free
        case (0)
#ifdef ENABLE_CUDA
            if (sgs_pointwise_cuda_enabled()) then
                if (calc_diag_active) call sgs_calc_diag_start(              &
                    SGS_CALC_BOTTOM, 'Sij bottom stress-free start')
!$cuf kernel do(2) <<<*,*>>>
                do jy = 1, ny
                do jx = 1, nx
                    S11(jx,jy,1) = dudx(jx,jy,1)
                    S12(jx,jy,1) = 0.5_rprec*(dudy(jx,jy,1)+dvdx(jx,jy,1))
                    S13(jx,jy,1) = 0.5_rprec*(dudz(jx,jy,1)+dwdx(jx,jy,1))
                    S22(jx,jy,1) = dvdy(jx,jy,1)
                    S23(jx,jy,1) = 0.5_rprec*(dvdz(jx,jy,1)+dwdy(jx,jy,1))
                    S33(jx,jy,1) = 0.5_rprec*dwdz(jx,jy,1)
                end do
                end do
                if (calc_diag_active) call sgs_calc_diag_stop(               &
                    SGS_CALC_BOTTOM, int(nx,8)*int(ny,8),                    &
                    'Sij bottom stress-free stop')
                call sgs_cuda_sync('Sij bottom stress-free')
            else
#endif
            do jy = 1, ny
            do jx = 1, nx
                ! Sij values are supposed to be on w-nodes for this case
                !   does that mean they (Sij) should all be zero?
                ! Check ux, uy, vx, vy, and qz
                ux = dudx(jx,jy,1) ! was 0.5_rprec*(dudx(jx,jy,1) + dudx(jx,jy,1))
                uy = dudy(jx,jy,1)
                uz = dudz(jx,jy,1)
                vx = dvdx(jx,jy,1)
                vy = dvdy(jx,jy,1)
                vz = dvdz(jx,jy,1)
                wx = dwdx(jx,jy,1)
                wy = dwdy(jx,jy,1)
                wz = 0.5_rprec*(dwdz(jx,jy,1) + 0._rprec)

                ! these values are stored on w-nodes
                S11(jx,jy,1) = ux
                S12(jx,jy,1) = 0.5_rprec*(uy+vx)
                S13(jx,jy,1) = 0.5_rprec*(uz+wx)
                S22(jx,jy,1) = vy
                S23(jx,jy,1) = 0.5_rprec*(vz+wy)
                S33(jx,jy,1) = wz
            end do
            end do
#ifdef ENABLE_CUDA
            end if
#endif

        ! Wall
        ! recall dudz and dvdz are stored on uvp-nodes for first level only,
        !   'wall' only
        ! recall dwdx and dwdy are stored on w-nodes (always)
        case (1:)
#ifdef ENABLE_CUDA
            if (sgs_pointwise_cuda_enabled()) then
                if (calc_diag_active) call sgs_calc_diag_start(              &
                    SGS_CALC_BOTTOM, 'Sij bottom wall start')
!$cuf kernel do(2) <<<*,*>>>
                do jy = 1, ny
                do jx = 1, nx
                    S11(jx,jy,1) = dudx(jx,jy,1)
                    S12(jx,jy,1) = 0.5_rprec*(dudy(jx,jy,1)+dvdx(jx,jy,1))
                    S13(jx,jy,1) = 0.5_rprec*(dudz(jx,jy,1) +                 &
                        0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2)))
                    S22(jx,jy,1) = dvdy(jx,jy,1)
                    S23(jx,jy,1) = 0.5_rprec*(dvdz(jx,jy,1) +                 &
                        0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2)))
                    S33(jx,jy,1) = dwdz(jx,jy,1)
                end do
                end do
                if (calc_diag_active) call sgs_calc_diag_stop(               &
                    SGS_CALC_BOTTOM, int(nx,8)*int(ny,8),                    &
                    'Sij bottom wall stop')
                call sgs_cuda_sync('Sij bottom wall')
            else
#endif
            do jy=1,ny
            do jx=1,nx
                ! these values stored on uvp-nodes
                S11(jx,jy,1) = dudx(jx,jy,1)
                S12(jx,jy,1) = 0.5_rprec*(dudy(jx,jy,1)+dvdx(jx,jy,1))
                wx = 0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))
                S13(jx,jy,1) = 0.5_rprec*(dudz(jx,jy,1)+wx)
                S22(jx,jy,1) = dvdy(jx,jy,1)
                wy = 0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))
                S23(jx,jy,1) = 0.5_rprec*(dvdz(jx,jy,1)+wy)
                S33(jx,jy,1) = dwdz(jx,jy,1)
            end do
            end do
#ifdef ENABLE_CUDA
            end if
#endif

    end select

    ! since first level already calculated
    jz_min = 2
else
    jz_min = 1
end if

! Calculate Sij for jz=nz (coord==nproc-1 only)
!   stored on uvp-nodes (this level only nz on w-grid --> nz-1 on uvp-grid)
!       for 'wall'
!   stored on w-nodes (all) for 'stress free'
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        case (0)

#ifdef ENABLE_CUDA
            if (sgs_pointwise_cuda_enabled()) then
                if (calc_diag_active) call sgs_calc_diag_start(              &
                    SGS_CALC_TOP, 'Sij top stress-free start')
!$cuf kernel do(2) <<<*,*>>>
                do jy = 1, ny
                do jx = 1, nx
                    S11(jx,jy,nz) = dudx(jx,jy,nz-1)
                    S12(jx,jy,nz) = 0.5_rprec*(dudy(jx,jy,nz-1) +            &
                        dvdx(jx,jy,nz-1))
                    S13(jx,jy,nz) = 0.5_rprec*(dudz(jx,jy,nz) +              &
                        dwdx(jx,jy,nz))
                    S22(jx,jy,nz) = dvdy(jx,jy,nz-1)
                    S23(jx,jy,nz) = 0.5_rprec*(dvdz(jx,jy,nz) +              &
                        dwdy(jx,jy,nz))
                    S33(jx,jy,nz) = 0.5_rprec*dwdz(jx,jy,nz-1)
                end do
                end do
                if (calc_diag_active) call sgs_calc_diag_stop(               &
                    SGS_CALC_TOP, int(nx,8)*int(ny,8),                       &
                    'Sij top stress-free stop')
                call sgs_cuda_sync('Sij top stress-free')
            else
#endif
            do jy=1,ny
            do jx=1,nx
                ! Sij values are supposed to be on w-nodes for this case
                !   does that mean they (Sij) should all be zero?
                ux = dudx(jx,jy,nz-1)
                uy = dudy(jx,jy,nz-1)
                uz = dudz(jx,jy,nz)   ! this comes from wallstress() i.e. zero
                vx = dvdx(jx,jy,nz-1)
                vy = dvdy(jx,jy,nz-1)
                vz = dvdz(jx,jy,nz)   ! this comes from wallstress() i.e. zero
                wx = dwdx(jx,jy,nz)
                wy = dwdy(jx,jy,nz)
                wz = 0.5_rprec*(dwdz(jx,jy,nz-1) + 0._rprec)

                ! these values are stored on w-nodes
                S11(jx,jy,nz) = ux
                S12(jx,jy,nz) = 0.5_rprec*(uy+vx)
                S13(jx,jy,nz) = 0.5_rprec*(uz+wx)
                S22(jx,jy,nz) = vy
                S23(jx,jy,nz) = 0.5_rprec*(vz+wy)
                S33(jx,jy,nz) = wz
            end do
            end do
#ifdef ENABLE_CUDA
            end if
#endif

        ! Wall
        ! recall dudz and dvdz are stored on uvp-nodes for first level only,
        !   'wall' only
        ! recall dwdx and dwdy are stored on w-nodes (always)
        case (1:)
#ifdef ENABLE_CUDA
            if (sgs_pointwise_cuda_enabled()) then
                if (calc_diag_active) call sgs_calc_diag_start(              &
                    SGS_CALC_TOP, 'Sij top wall start')
!$cuf kernel do(2) <<<*,*>>>
                do jy = 1, ny
                do jx = 1, nx
                    S11(jx,jy,nz) = dudx(jx,jy,nz-1)
                    S12(jx,jy,nz) = 0.5_rprec*(dudy(jx,jy,nz-1) +            &
                        dvdx(jx,jy,nz-1))
                    S13(jx,jy,nz) = 0.5_rprec*(dudz(jx,jy,nz) +              &
                        0.5_rprec*(dwdx(jx,jy,nz-1)+dwdx(jx,jy,nz)))
                    S22(jx,jy,nz) = dvdy(jx,jy,nz-1)
                    S23(jx,jy,nz) = 0.5_rprec*(dvdz(jx,jy,nz) +              &
                        0.5_rprec*(dwdy(jx,jy,nz-1)+dwdy(jx,jy,nz)))
                    S33(jx,jy,nz) = dwdz(jx,jy,nz-1)
                end do
                end do
                if (calc_diag_active) call sgs_calc_diag_stop(               &
                    SGS_CALC_TOP, int(nx,8)*int(ny,8),                       &
                    'Sij top wall stop')
                call sgs_cuda_sync('Sij top wall')
            else
#endif
            do jy = 1, ny
            do jx = 1, nx
                ! these values stored on uvp-nodes
                S11(jx,jy,nz) = dudx(jx,jy,nz-1)
                S12(jx,jy,nz) = 0.5_rprec*(dudy(jx,jy,nz-1)+dvdx(jx,jy,nz-1))
                wx = 0.5_rprec*(dwdx(jx,jy,nz-1)+dwdx(jx,jy,nz))
                ! dudz from wallstress()
                S13(jx,jy,nz) = 0.5_rprec*(dudz(jx,jy,nz)+wx)
                S22(jx,jy,nz) = dvdy(jx,jy,nz-1)
                wy = 0.5_rprec*(dwdy(jx,jy,nz-1)+dwdy(jx,jy,nz))
                ! dvdz from wallstress()
                S23(jx,jy,nz) = 0.5_rprec*(dvdz(jx,jy,nz)+wy)
                S33(jx,jy,nz) = dwdz(jx,jy,nz-1)
            end do
            end do
#ifdef ENABLE_CUDA
            end if
#endif

    end select

    ! since last level already calculated
    jz_max = nz-1
else
    jz_max = nz
end if

#ifdef PPMPI
! dudz calculated for 0:nz-1 (on w-nodes) except bottom process
! (only 1:nz-1) exchange information between processors to set
! values at nz from jz=1 above to jz=nz below
#ifdef ENABLE_CUDA
if (dwdz_halo_posted) then
    ! The actual wait is deferred until after independent Sij planes are done.
else if (sgs_pointwise_cuda_enabled() .and. sgs_direct_dwdz_halo_enabled()) then
    if (calc_diag_active) call sgs_calc_cpu_start(SGS_CALC_DWDZ_MPI_WAIT)
    if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_MPI)
    call sgs_sync_dwdz_down_cuda(dwdz(:,:,1:))
    if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_MPI)
    if (calc_diag_active) call sgs_calc_cpu_stop(SGS_CALC_DWDZ_MPI_WAIT)
else
#endif
#ifdef ENABLE_CUDA
    if (calc_diag_active) call sgs_calc_cpu_start(SGS_CALC_DWDZ_MPI_WAIT)
    if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_MPI)
#endif
    call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#ifdef ENABLE_CUDA
    if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_MPI)
    if (calc_diag_active) call sgs_calc_cpu_stop(SGS_CALC_DWDZ_MPI_WAIT)
#endif
#ifdef ENABLE_CUDA
end if
#endif
#endif

! Calculate Sij for the rest of the domain
!   values are stored on w-nodes
!   dudz, dvdz, dwdx, dwdy are already stored on w-nodes
#ifdef ENABLE_CUDA
#ifdef PPMPI
if (dwdz_halo_posted) then
    sij_bulk_max = jz_max
    if (coord < nproc - 1) sij_bulk_max = min(sij_bulk_max, nz-1)
    sij_bulk_min = jz_min
    sij_exp_bulk_min = jz_min
    sij_exp_bulk_max = sij_bulk_max
    if (jz_min <= sij_bulk_max) then
        if (calc_diag_active) call sgs_calc_diag_start(SGS_CALC_INTERIOR,     &
            'Sij overlap interior start')
        if (calc_sij_explicit) then
            sij_block_x = 32
            sij_block_y = 4
            sij_block_z = 2
            sij_grid_x = (nx + sij_block_x - 1) / sij_block_x
            sij_grid_y = (ny + sij_block_y - 1) / sij_block_y
            sij_grid_z = (sij_bulk_max - jz_min + 1 + sij_block_z - 1) /      &
                sij_block_z
            sij_block = dim3(sij_block_x, sij_block_y, sij_block_z)
            sij_grid = dim3(sij_grid_x, sij_grid_y, sij_grid_z)
            call calc_sij_interior_kernel<<<sij_grid,sij_block>>>(dudx, dudy,&
                dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, S11, S12, S13,    &
                S22, S23, S33, ld, nx, ny, lbz, jz_min, sij_bulk_max)
            calc_sij_explicit_used = .true.
            sij_exp_bulk_min = jz_min
            sij_exp_bulk_max = sij_bulk_max
        else
!$cuf kernel do(3) <<<*,*>>>
            do jz = jz_min, sij_bulk_max
            do jy = 1, ny
            do jx = 1, nx
                S11(jx,jy,jz) = 0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))
                uy = (dudy(jx,jy,jz) + dudy(jx,jy,jz-1))
                vx = (dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))
                S12(jx,jy,jz) = 0.25_rprec*(uy+vx)
                S13(jx,jy,jz) = 0.5_rprec*(dudz(jx,jy,jz) + dwdx(jx,jy,jz))
                S22(jx,jy,jz) = 0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))
                S23(jx,jy,jz) = 0.5_rprec*(dvdz(jx,jy,jz) + dwdy(jx,jy,jz))
                S33(jx,jy,jz) = 0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))
            end do
            end do
            end do
        end if
        if (calc_diag_active) call sgs_calc_diag_stop(SGS_CALC_INTERIOR,      &
            int(nx,8)*int(ny,8)*int(sij_bulk_max-jz_min+1,8),                &
            'Sij overlap interior stop')
        call sgs_cuda_sync('Sij overlap interior')
    end if

    if (dwdz_nreq > 0) then
        if (calc_diag_active) call sgs_calc_cpu_start(SGS_CALC_DWDZ_MPI_WAIT)
        if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_MPI)
        call mpi_waitall(dwdz_nreq, dwdz_req, dwdz_status, ierr)
        if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_MPI)
        if (calc_diag_active) call sgs_calc_cpu_stop(SGS_CALC_DWDZ_MPI_WAIT)
        if (calc_diag_active) then
            call sgs_diag_time(sgs_dwdz_arr_after_mpi)
            if (.not. dwdz_device_halo_posted .and.                         &
                .not. dwdz_hostpinned_halo_posted)                           &
                sgs_dwdz_arr_after_post = sgs_dwdz_arr_after_mpi
        end if
    end if
    if (dwdz_hostpinned_halo_posted .and. dwdz_recv_count > 0) then
        if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_H2D)
        sgs_dwdz_recv_d(1:dwdz_recv_count) =                                 &
            sgs_dwdz_recv_h(1:dwdz_recv_count)
        call sgs_cuda_barrier('dwdz hostpinned H2D')
        if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_H2D)
    end if
    if ((dwdz_device_halo_posted .or. dwdz_hostpinned_halo_posted) .and.     &
        dwdz_recv_count > 0) then
        if (calc_diag_active) call sgs_dwdz_detail_start(SGS_DWDZ_UNPACK)
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            dwdz_idx = jx + (jy - 1) * ld
            dwdz(jx,jy,nz) = sgs_dwdz_recv_d(dwdz_idx)
        end do
        end do
        call sgs_cuda_barrier('dwdz device halo unpack')
        if (calc_diag_active) call sgs_dwdz_detail_stop(SGS_DWDZ_UNPACK)
        if (calc_diag_active) call sgs_diag_time(sgs_dwdz_arr_after_post)
    end if
    if ((coord < nproc - 1) .and. (jz_min <= nz) .and. (jz_max >= nz)) then
        if (calc_diag_active) call sgs_calc_diag_start(SGS_CALC_HALO_PLANE,   &
            'Sij overlap halo plane start')
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, nx
            S11(jx,jy,nz) = 0.5_rprec*(dudx(jx,jy,nz) + dudx(jx,jy,nz-1))
            uy = (dudy(jx,jy,nz) + dudy(jx,jy,nz-1))
            vx = (dvdx(jx,jy,nz) + dvdx(jx,jy,nz-1))
            S12(jx,jy,nz) = 0.25_rprec*(uy+vx)
            S13(jx,jy,nz) = 0.5_rprec*(dudz(jx,jy,nz) + dwdx(jx,jy,nz))
            S22(jx,jy,nz) = 0.5_rprec*(dvdy(jx,jy,nz) + dvdy(jx,jy,nz-1))
            S23(jx,jy,nz) = 0.5_rprec*(dvdz(jx,jy,nz) + dwdy(jx,jy,nz))
            S33(jx,jy,nz) = 0.5_rprec*(dwdz(jx,jy,nz) + dwdz(jx,jy,nz-1))
        end do
        end do
        if (calc_diag_active) call sgs_calc_diag_stop(SGS_CALC_HALO_PLANE,    &
            int(nx,8)*int(ny,8), 'Sij overlap halo plane stop')
        call sgs_cuda_sync('Sij overlap halo plane')
    end if
else if (sgs_pointwise_cuda_enabled()) then
#else
if (sgs_pointwise_cuda_enabled()) then
#endif
    if (calc_diag_active) call sgs_calc_diag_start(SGS_CALC_INTERIOR,         &
        'Sij interior start')
    sij_bulk_min = jz_min
    sij_exp_bulk_min = jz_min
    sij_exp_bulk_max = jz_max
    if (calc_sij_explicit) then
        sij_block_x = 32
        sij_block_y = 4
        sij_block_z = 2
        sij_grid_x = (nx + sij_block_x - 1) / sij_block_x
        sij_grid_y = (ny + sij_block_y - 1) / sij_block_y
        sij_grid_z = (jz_max - jz_min + 1 + sij_block_z - 1) / sij_block_z
        sij_block = dim3(sij_block_x, sij_block_y, sij_block_z)
        sij_grid = dim3(sij_grid_x, sij_grid_y, sij_grid_z)
        call calc_sij_interior_kernel<<<sij_grid,sij_block>>>(dudx, dudy,    &
            dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, S11, S12, S13, S22,   &
            S23, S33, ld, nx, ny, lbz, jz_min, jz_max)
        calc_sij_explicit_used = .true.
        sij_exp_bulk_min = jz_min
        sij_exp_bulk_max = jz_max
    else
!$cuf kernel do(3) <<<*,*>>>
        do jz = jz_min, jz_max
        do jy = 1, ny
        do jx = 1, nx
            S11(jx,jy,jz) = 0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))
            uy = (dudy(jx,jy,jz) + dudy(jx,jy,jz-1))
            vx = (dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))
            S12(jx,jy,jz) = 0.25_rprec*(uy+vx)
            S13(jx,jy,jz) = 0.5_rprec*(dudz(jx,jy,jz) + dwdx(jx,jy,jz))
            S22(jx,jy,jz) = 0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))
            S23(jx,jy,jz) = 0.5_rprec*(dvdz(jx,jy,jz) + dwdy(jx,jy,jz))
            S33(jx,jy,jz) = 0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))
        end do
        end do
        end do
    end if
    if (calc_diag_active) call sgs_calc_diag_stop(SGS_CALC_INTERIOR,          &
        int(nx,8)*int(ny,8)*int(jz_max-jz_min+1,8), 'Sij interior stop')
    call sgs_cuda_sync('Sij interior')
else
#endif
do jz = jz_min, jz_max
do jy = 1, ny
do jx = 1, nx
    S11(jx,jy,jz) = 0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))
    uy = (dudy(jx,jy,jz) + dudy(jx,jy,jz-1))
    vx = (dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))
    S12(jx,jy,jz) = 0.25_rprec*(uy+vx)
    S13(jx,jy,jz) = 0.5_rprec*(dudz(jx,jy,jz) + dwdx(jx,jy,jz))
    S22(jx,jy,jz) = 0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))
    S23(jx,jy,jz) = 0.5_rprec*(dvdz(jx,jy,jz) + dwdy(jx,jy,jz))
    S33(jx,jy,jz) = 0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))
end do
end do
end do
#ifdef ENABLE_CUDA
end if
#endif

#ifdef ENABLE_CUDA
if (sij_bulk_min == 0) sij_bulk_min = jz_min
if (sij_exp_bulk_min == 0) sij_exp_bulk_min = sij_bulk_min
if (sij_exp_bulk_max == 0) sij_exp_bulk_max = jz_max
if (calc_diag_active) call sgs_calc_set_audit(jz_min, jz_max, sij_bulk_min,  &
    sij_exp_bulk_max, coord == 0, coord == nproc - 1, lbc_mom, ubc_mom,      &
    calc_sij_explicit_used, sij_block_x, sij_block_y, sij_block_z,           &
    sij_grid_x, sij_grid_y, sij_grid_z)
if (sgs_calc_sij_device_bench_enabled())                                     &
    call sgs_calc_sij_device_lower_bound_bench(jz_min, jz_max)
#endif

end subroutine calc_Sij

#ifdef ENABLE_CUDA
!*******************************************************************************
subroutine calc_Sij_nut_dynamic_cuda()
!*******************************************************************************
!
! Active dynamic-SGS fast path for non-coefficient-update timesteps.  The
! dynamic models use l=delta, so Sij and Nu_t can be constructed in one pass.
! Boundary handling matches the lower-wall / upper-stress-free case used by the
! ATM validation configuration. Other boundary/model combinations use the
! original separate calc_Sij and Nu_t path.
!
use types, only : rprec
use param
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sgs_param
use cudafor
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
implicit none

real(rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
real(rprec) :: S_mag, l_delta2
integer :: jx, jy, jz

#ifdef PPMPI
call mpi_sync_real_array( dwdz(:,:,1:), 1, MPI_SYNC_DOWN )
#endif

l_delta2 = delta * delta

!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz
do jy = 1, ny
do jx = 1, nx
    if ((coord == 0) .and. (jz == 1)) then
        ux = dudx(jx,jy,1)
        uy = dudy(jx,jy,1)
        uz = dudz(jx,jy,1)
        vx = dvdx(jx,jy,1)
        vy = dvdy(jx,jy,1)
        vz = dvdz(jx,jy,1)
        wx = 0.5_rprec * (dwdx(jx,jy,1) + dwdx(jx,jy,2))
        wy = 0.5_rprec * (dwdy(jx,jy,1) + dwdy(jx,jy,2))
        wz = dwdz(jx,jy,1)
    else if ((coord == nproc - 1) .and. (jz == nz)) then
        ux = dudx(jx,jy,nz-1)
        uy = dudy(jx,jy,nz-1)
        uz = dudz(jx,jy,nz)
        vx = dvdx(jx,jy,nz-1)
        vy = dvdy(jx,jy,nz-1)
        vz = dvdz(jx,jy,nz)
        wx = dwdx(jx,jy,nz)
        wy = dwdy(jx,jy,nz)
        wz = 0.5_rprec * dwdz(jx,jy,nz-1)
    else
        ux = 0.5_rprec * (dudx(jx,jy,jz) + dudx(jx,jy,jz-1))
        uy = 0.5_rprec * (dudy(jx,jy,jz) + dudy(jx,jy,jz-1))
        uz = dudz(jx,jy,jz)
        vx = 0.5_rprec * (dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))
        vy = 0.5_rprec * (dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))
        vz = dvdz(jx,jy,jz)
        wx = dwdx(jx,jy,jz)
        wy = dwdy(jx,jy,jz)
        wz = 0.5_rprec * (dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))
    end if

    S11(jx,jy,jz) = ux
    S12(jx,jy,jz) = 0.5_rprec * (uy + vx)
    S13(jx,jy,jz) = 0.5_rprec * (uz + wx)
    S22(jx,jy,jz) = vy
    S23(jx,jy,jz) = 0.5_rprec * (vz + wy)
    S33(jx,jy,jz) = wz

    S_mag = sqrt(2._rprec * (S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +            &
        S33(jx,jy,jz)**2 + 2._rprec * (S12(jx,jy,jz)**2 +                    &
        S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
    Nu_t(jx,jy,jz) = S_mag * Cs_opt2(jx,jy,jz) * l_delta2
    if (jz == nz) S(jx,jy) = S_mag
end do
end do
end do

call sgs_cuda_sync('Sij/Nu_t fused dynamic')

end subroutine calc_Sij_nut_dynamic_cuda
#endif

#ifdef ENABLE_CUDA
#ifdef PPMPI
!*******************************************************************************
subroutine sgs_sync_dwdz_down_cuda(dwdz_arr)
!*******************************************************************************
!
! SGS only needs the downward dwdz halo: plane 1 from coord+1 into local nz.
! This direct exchange avoids the generic sync wrapper and skips null neighbors.
!
use types, only : rprec
use param, only : ld, ny, nz, coord, nproc, down, up, comm, ierr, MPI_RPREC
use mpi
implicit none

real(rprec), managed, intent(inout) :: dwdz_arr(ld,ny,1:nz)
integer :: req(2)
integer :: statuses(MPI_STATUS_SIZE, 2)
integer :: nreq

req = MPI_REQUEST_NULL
nreq = 0

if (coord < nproc - 1) then
    nreq = nreq + 1
    call mpi_irecv(dwdz_arr(1,1,nz), ld*ny, MPI_RPREC, up, 91, comm,          &
        req(nreq), ierr)
end if
if (coord > 0) then
    nreq = nreq + 1
    call mpi_isend(dwdz_arr(1,1,1), ld*ny, MPI_RPREC, down, 91, comm,         &
        req(nreq), ierr)
end if
if (nreq > 0) then
    call mpi_waitall(nreq, req, statuses, ierr)
end if

end subroutine sgs_sync_dwdz_down_cuda
#endif
#endif

#ifdef ENABLE_CUDA
!*******************************************************************************
subroutine sgs_calc_sij_detail_report(stage_count)
!*******************************************************************************
use param, only : coord, nproc
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC
use mpi
#endif
implicit none

integer, intent(in) :: stage_count
integer :: i, r
integer :: launch_max(SGS_CALC_COUNT), total_launch_max
real(rprec) :: cpu_max(SGS_CALC_COUNT), gpu_max(SGS_CALC_COUNT)
real(rprec) :: total_gpu, total_gpu_max, avg_gpu
real(rprec) :: local_avg_gpu
real(rprec) :: dwdz_total, dwdz_bytes, dwdz_gbs
real(rprec) :: dwdz_ready_min, dwdz_ready_max, dwdz_ready_rel, dwdz_skew
character(len=32) :: labels(SGS_CALC_COUNT)

labels(SGS_CALC_BOTTOM) = 'bottom boundary'
labels(SGS_CALC_TOP) = 'top boundary'
labels(SGS_CALC_INTERIOR) = 'interior/bulk'
labels(SGS_CALC_HALO_PLANE) = 'halo plane'
labels(SGS_CALC_DWDZ_PRE_SYNC) = 'dwdz pre-MPI sync'
labels(SGS_CALC_DWDZ_MPI_WAIT) = 'dwdz MPI wait'

do i = 1, SGS_CALC_COUNT
    if (sgs_calc_evt_active(i)) then
        call sgs_event_elapsed_seconds(sgs_calc_evt_start(i),                 &
            sgs_calc_evt_stop(i), sgs_calc_gpu(i), trim(labels(i)))
    else
        sgs_calc_gpu(i) = 0._rprec
    end if
end do
total_gpu = sum(sgs_calc_gpu)

#ifdef PPMPI
call mpi_allreduce(sgs_calc_cpu, cpu_max, SGS_CALC_COUNT, MPI_RPREC, MPI_MAX,&
    comm, ierr)
call mpi_allreduce(sgs_calc_gpu, gpu_max, SGS_CALC_COUNT, MPI_RPREC, MPI_MAX,&
    comm, ierr)
call mpi_allreduce(sgs_calc_launches, launch_max, SGS_CALC_COUNT,            &
    MPI_INTEGER, MPI_MAX, comm, ierr)
call mpi_allreduce(sgs_calc_total_launches, total_launch_max, 1,             &
    MPI_INTEGER, MPI_MAX, comm, ierr)
call mpi_allreduce(total_gpu, total_gpu_max, 1, MPI_RPREC, MPI_MAX, comm,    &
    ierr)
#else
cpu_max = sgs_calc_cpu
gpu_max = sgs_calc_gpu
launch_max = sgs_calc_launches
total_launch_max = sgs_calc_total_launches
total_gpu_max = total_gpu
#endif

if (total_launch_max > 0) then
    avg_gpu = total_gpu_max / real(total_launch_max, rprec)
else
    avg_gpu = 0._rprec
end if
if (sgs_calc_total_launches > 0) then
    local_avg_gpu = total_gpu / real(sgs_calc_total_launches, rprec)
else
    local_avg_gpu = 0._rprec
end if

#ifdef PPMPI
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,'(a,i8,a,i6,a,i6,a,i6,a,i12,a,i8,a,E15.7,a,E15.7)')         &
            'SGS calc_Sij rank timing call=', stage_count, ' rank=', coord,  &
            ' z_range=', sgs_calc_jz_min, ':', sgs_calc_jz_max,              &
            ' cells=', sgs_calc_total_cells, ' launches=',                   &
            sgs_calc_total_launches, ' total_gpu=', total_gpu,               &
            ' avg_gpu_per_launch=', local_avg_gpu
        write(*,*) 'SGS_CALCSIJ_RANK_SPLIT call=', stage_count,              &
            ' rank=', coord,                                                 &
            ' bottom_gpu=', sgs_calc_gpu(SGS_CALC_BOTTOM),                   &
            ' top_gpu=', sgs_calc_gpu(SGS_CALC_TOP),                         &
            ' interior_gpu=', sgs_calc_gpu(SGS_CALC_INTERIOR),               &
            ' halo_gpu=', sgs_calc_gpu(SGS_CALC_HALO_PLANE),                 &
            ' dwdz_mpi=', sgs_calc_cpu(SGS_CALC_DWDZ_MPI_WAIT),              &
            ' total_gpu=', total_gpu,                                        &
            ' bottom_cpu=', sgs_calc_cpu(SGS_CALC_BOTTOM),                   &
            ' top_cpu=', sgs_calc_cpu(SGS_CALC_TOP),                         &
            ' interior_cpu=', sgs_calc_cpu(SGS_CALC_INTERIOR),               &
            ' halo_cpu=', sgs_calc_cpu(SGS_CALC_HALO_PLANE),                 &
            ' launches_bottom=', sgs_calc_launches(SGS_CALC_BOTTOM),         &
            ' launches_top=', sgs_calc_launches(SGS_CALC_TOP),               &
            ' launches_interior=', sgs_calc_launches(SGS_CALC_INTERIOR),     &
            ' launches_halo=', sgs_calc_launches(SGS_CALC_HALO_PLANE)
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
#else
write(*,'(a,i8,a,i6,a,i6,a,i6,a,i12,a,i8,a,E15.7,a,E15.7)')                 &
    'SGS calc_Sij rank timing call=', stage_count, ' rank=', coord,          &
    ' z_range=', sgs_calc_jz_min, ':', sgs_calc_jz_max,                      &
    ' cells=', sgs_calc_total_cells, ' launches=', sgs_calc_total_launches,  &
    ' total_gpu=', total_gpu, ' avg_gpu_per_launch=', local_avg_gpu
write(*,*) 'SGS_CALCSIJ_RANK_SPLIT call=', stage_count, ' rank=', coord,     &
    ' bottom_gpu=', sgs_calc_gpu(SGS_CALC_BOTTOM),                           &
    ' top_gpu=', sgs_calc_gpu(SGS_CALC_TOP),                                 &
    ' interior_gpu=', sgs_calc_gpu(SGS_CALC_INTERIOR),                       &
    ' halo_gpu=', sgs_calc_gpu(SGS_CALC_HALO_PLANE),                         &
    ' dwdz_mpi=', sgs_calc_cpu(SGS_CALC_DWDZ_MPI_WAIT),                      &
    ' total_gpu=', total_gpu,                                                &
    ' bottom_cpu=', sgs_calc_cpu(SGS_CALC_BOTTOM),                           &
    ' top_cpu=', sgs_calc_cpu(SGS_CALC_TOP),                                 &
    ' interior_cpu=', sgs_calc_cpu(SGS_CALC_INTERIOR),                       &
    ' halo_cpu=', sgs_calc_cpu(SGS_CALC_HALO_PLANE),                         &
    ' launches_bottom=', sgs_calc_launches(SGS_CALC_BOTTOM),                 &
    ' launches_top=', sgs_calc_launches(SGS_CALC_TOP),                       &
    ' launches_interior=', sgs_calc_launches(SGS_CALC_INTERIOR),             &
    ' launches_halo=', sgs_calc_launches(SGS_CALC_HALO_PLANE)
flush(6)
#endif

dwdz_total = sum(sgs_dwdz_detail)
dwdz_bytes = real(sgs_dwdz_send_bytes + sgs_dwdz_recv_bytes, rprec)
dwdz_gbs = 0._rprec
if (sgs_dwdz_detail(SGS_DWDZ_MPI) > 0._rprec)                               &
    dwdz_gbs = dwdz_bytes / sgs_dwdz_detail(SGS_DWDZ_MPI) / 1.0e9_rprec
#ifdef PPMPI
call mpi_allreduce(sgs_dwdz_arr_after_sync, dwdz_ready_min, 1, MPI_RPREC,   &
    MPI_MIN, comm, ierr)
call mpi_allreduce(sgs_dwdz_arr_after_sync, dwdz_ready_max, 1, MPI_RPREC,   &
    MPI_MAX, comm, ierr)
dwdz_ready_rel = 0._rprec
dwdz_skew = 0._rprec
if (dwdz_ready_max > 0._rprec) then
    dwdz_ready_rel = sgs_dwdz_arr_after_sync - dwdz_ready_min
    dwdz_skew = dwdz_ready_max - dwdz_ready_min
end if
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,'(a,i0,a,i0,5(a,es14.6),3(a,i0))')                          &
            'SGS_DWDZ_ARR call=', stage_count, ' rank=', coord,             &
            ' pre=', sgs_dwdz_detail(SGS_DWDZ_PRE_SYNC),                    &
            ' ready_rel=', dwdz_ready_rel,                                   &
            ' mpi=', sgs_dwdz_detail(SGS_DWDZ_MPI),                         &
            ' skew=', dwdz_skew,                                             &
            ' barrier=', sgs_dwdz_arr_barrier, ' bytes=',                    &
            sgs_dwdz_send_bytes + sgs_dwdz_recv_bytes,                       &
            ' dummy_bytes=', sgs_dwdz_dummy_bytes, ' mpi_calls=',            &
            sgs_dwdz_mpi_calls
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,*) 'SGS dwdz halo rank comm call=', stage_count, ' rank=',   &
            coord,                                                           &
            ' pre=', sgs_dwdz_detail(SGS_DWDZ_PRE_SYNC),                    &
            ' pack=', sgs_dwdz_detail(SGS_DWDZ_PACK),                       &
            ' mpi=', sgs_dwdz_detail(SGS_DWDZ_MPI),                         &
            ' unpack=', sgs_dwdz_detail(SGS_DWDZ_UNPACK),                   &
            ' post=', sgs_dwdz_detail(SGS_DWDZ_POST_SYNC),                  &
            ' d2h=', sgs_dwdz_detail(SGS_DWDZ_D2H),                         &
            ' h2d=', sgs_dwdz_detail(SGS_DWDZ_H2D),                         &
            ' dummy=', sgs_dwdz_detail(SGS_DWDZ_DUMMY),                     &
            ' barrier=', sgs_dwdz_arr_barrier,                              &
            ' total=', dwdz_total, ' mpi_calls=', sgs_dwdz_mpi_calls,        &
            ' send_bytes=', sgs_dwdz_send_bytes, ' recv_bytes=',             &
            sgs_dwdz_recv_bytes, ' dummy_bytes=', sgs_dwdz_dummy_bytes,      &
            ' neighbor=', sgs_dwdz_neighbor,                                 &
            ' combined=', sgs_dwdz_combined_msg, ' mpi_GBps=', dwdz_gbs
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
#else
write(*,*) 'SGS dwdz halo rank comm call=', stage_count, ' rank=', coord,    &
    ' pre=', sgs_dwdz_detail(SGS_DWDZ_PRE_SYNC),                            &
    ' pack=', sgs_dwdz_detail(SGS_DWDZ_PACK),                               &
    ' mpi=', sgs_dwdz_detail(SGS_DWDZ_MPI),                                 &
    ' unpack=', sgs_dwdz_detail(SGS_DWDZ_UNPACK),                           &
    ' post=', sgs_dwdz_detail(SGS_DWDZ_POST_SYNC),                          &
    ' d2h=', sgs_dwdz_detail(SGS_DWDZ_D2H),                                  &
    ' h2d=', sgs_dwdz_detail(SGS_DWDZ_H2D),                                  &
    ' dummy=', sgs_dwdz_detail(SGS_DWDZ_DUMMY),                              &
    ' barrier=', sgs_dwdz_arr_barrier,                                      &
    ' total=', dwdz_total, ' mpi_calls=', sgs_dwdz_mpi_calls,                &
    ' send_bytes=', sgs_dwdz_send_bytes, ' recv_bytes=', sgs_dwdz_recv_bytes,&
    ' dummy_bytes=', sgs_dwdz_dummy_bytes, ' neighbor=', sgs_dwdz_neighbor,  &
    ' combined=', sgs_dwdz_combined_msg,                                     &
    ' mpi_GBps=', dwdz_gbs
flush(6)
#endif

if (.not. sgs_calc_audit_printed) then
#ifdef PPMPI
    do r = 0, nproc - 1
        call mpi_barrier(comm, ierr)
        if (coord == r) then
            write(*,'(a,i6,a,i6,a,i6,a,i6,a,i12,a,i8)')                     &
                'SGS calc_Sij work audit rank=', coord, ' nproc=', nproc,    &
                ' z_range=', sgs_calc_jz_min, ':', sgs_calc_jz_max,          &
                ' cells=', sgs_calc_total_cells, ' launches=',               &
                sgs_calc_total_launches
            write(*,*) 'SGS_CALCSIJ_WORK_AUDIT rank=', coord,                &
                ' nproc=', nproc,                                            &
                ' local_z_min=', sgs_calc_jz_min,                            &
                ' local_z_max=', sgs_calc_jz_max,                            &
                ' interior_z_min=', sgs_calc_bulk_jz_min,                    &
                ' interior_z_max=', sgs_calc_bulk_jz_max,                    &
                ' bottom_branch=', sgs_calc_bottom_branch,                   &
                ' top_branch=', sgs_calc_top_branch,                         &
                ' lbc_mom=', sgs_calc_lbc_mode,                              &
                ' ubc_mom=', sgs_calc_ubc_mode,                              &
                ' explicit=', sgs_calc_explicit_used,                        &
                ' cells_bottom=', sgs_calc_cells(SGS_CALC_BOTTOM),           &
                ' cells_top=', sgs_calc_cells(SGS_CALC_TOP),                 &
                ' cells_interior=', sgs_calc_cells(SGS_CALC_INTERIOR),       &
                ' cells_halo=', sgs_calc_cells(SGS_CALC_HALO_PLANE),         &
                ' launches_bottom=', sgs_calc_launches(SGS_CALC_BOTTOM),     &
                ' launches_top=', sgs_calc_launches(SGS_CALC_TOP),           &
                ' launches_interior=', sgs_calc_launches(SGS_CALC_INTERIOR), &
                ' launches_halo=', sgs_calc_launches(SGS_CALC_HALO_PLANE),   &
                ' block=', sgs_calc_block_x, 'x', sgs_calc_block_y, 'x',     &
                sgs_calc_block_z, ' grid=', sgs_calc_grid_x, 'x',            &
                sgs_calc_grid_y, 'x', sgs_calc_grid_z
            write(*,'(a,i6,a)') 'SGS calc_Sij bounds audit rank=', coord,    &
                ' all active CUDA kernels use local rank z bounds; no global z loop bounds were observed.'
        end if
    end do
    call mpi_barrier(comm, ierr)
#else
    write(*,'(a,i6,a,i6,a,i6,a,i6,a,i12,a,i8)')                             &
        'SGS calc_Sij work audit rank=', coord, ' nproc=', nproc,            &
        ' z_range=', sgs_calc_jz_min, ':', sgs_calc_jz_max,                  &
        ' cells=', sgs_calc_total_cells, ' launches=',                       &
        sgs_calc_total_launches
    write(*,*) 'SGS_CALCSIJ_WORK_AUDIT rank=', coord, ' nproc=', nproc,      &
        ' local_z_min=', sgs_calc_jz_min,                                    &
        ' local_z_max=', sgs_calc_jz_max,                                    &
        ' interior_z_min=', sgs_calc_bulk_jz_min,                            &
        ' interior_z_max=', sgs_calc_bulk_jz_max,                            &
        ' bottom_branch=', sgs_calc_bottom_branch,                           &
        ' top_branch=', sgs_calc_top_branch,                                 &
        ' lbc_mom=', sgs_calc_lbc_mode,                                      &
        ' ubc_mom=', sgs_calc_ubc_mode,                                      &
        ' explicit=', sgs_calc_explicit_used,                                &
        ' cells_bottom=', sgs_calc_cells(SGS_CALC_BOTTOM),                   &
        ' cells_top=', sgs_calc_cells(SGS_CALC_TOP),                         &
        ' cells_interior=', sgs_calc_cells(SGS_CALC_INTERIOR),               &
        ' cells_halo=', sgs_calc_cells(SGS_CALC_HALO_PLANE),                 &
        ' launches_bottom=', sgs_calc_launches(SGS_CALC_BOTTOM),             &
        ' launches_top=', sgs_calc_launches(SGS_CALC_TOP),                   &
        ' launches_interior=', sgs_calc_launches(SGS_CALC_INTERIOR),         &
        ' launches_halo=', sgs_calc_launches(SGS_CALC_HALO_PLANE),           &
        ' block=', sgs_calc_block_x, 'x', sgs_calc_block_y, 'x',             &
        sgs_calc_block_z, ' grid=', sgs_calc_grid_x, 'x', sgs_calc_grid_y,   &
        'x', sgs_calc_grid_z
    write(*,'(a,i6,a)') 'SGS calc_Sij bounds audit rank=', coord,            &
        ' all active CUDA kernels use local rank z bounds; no global z loop bounds were observed.'
#endif
    sgs_calc_audit_printed = .true.
end if

if (coord == 0) then
    write(*,'(a,i8)') 'SGS calc_Sij split timing (max rank), call ',          &
        stage_count
    do i = 1, SGS_CALC_COUNT
        write(*,'(3a,E15.7)') '  ', trim(labels(i)), ' CPU wall: ', cpu_max(i)
        write(*,'(3a,E15.7)') '  ', trim(labels(i)), ' GPU event: ', gpu_max(i)
    end do
    write(*,'(a,i8)') '  calc_Sij CUDA/CUF launches: ', total_launch_max
    write(*,'(1a,E15.7)') '  calc_Sij avg GPU event per launch: ', avg_gpu
    write(*,'(1a,E15.7)') '  calc_Sij total GPU event: ', total_gpu_max
end if

end subroutine sgs_calc_sij_detail_report

!*******************************************************************************
subroutine sgs_tau_halo_detail_report(stage_count)
!*******************************************************************************
use param, only : coord
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC, nproc
use mpi
#endif
implicit none

integer, intent(in) :: stage_count
integer :: r
real(rprec) :: detail_max(SGS_TAU_COUNT)
real(rprec) :: send_bytes, recv_bytes, send_bytes_max, recv_bytes_max
real(rprec) :: pack_bytes, mpi_bytes, unpack_bytes
real(rprec) :: pack_gbs, mpi_gbs, unpack_gbs
real(rprec) :: tau_total, tau_bytes, tau_gbs
real(rprec) :: tau_ready_min, tau_ready_max, tau_ready_rel, tau_skew
real(rprec) :: tau_sync_to_mpi

send_bytes = real(sgs_tau_send_bytes, rprec)
recv_bytes = real(sgs_tau_recv_bytes, rprec)
#ifdef PPMPI
call mpi_allreduce(sgs_tau_detail, detail_max, SGS_TAU_COUNT, MPI_RPREC,     &
    MPI_MAX, comm, ierr)
call mpi_allreduce(send_bytes, send_bytes_max, 1, MPI_RPREC, MPI_MAX, comm,  &
    ierr)
call mpi_allreduce(recv_bytes, recv_bytes_max, 1, MPI_RPREC, MPI_MAX, comm,  &
    ierr)
#else
detail_max = sgs_tau_detail
send_bytes_max = send_bytes
recv_bytes_max = recv_bytes
#endif

pack_bytes = send_bytes_max
mpi_bytes = send_bytes_max + recv_bytes_max
unpack_bytes = recv_bytes_max
pack_gbs = 0._rprec
mpi_gbs = 0._rprec
unpack_gbs = 0._rprec
if (detail_max(SGS_TAU_PACK) > 0._rprec)                                     &
    pack_gbs = pack_bytes / detail_max(SGS_TAU_PACK) / 1.0e9_rprec
if (detail_max(SGS_TAU_MPI) > 0._rprec)                                      &
    mpi_gbs = mpi_bytes / detail_max(SGS_TAU_MPI) / 1.0e9_rprec
if (detail_max(SGS_TAU_UNPACK) > 0._rprec)                                   &
    unpack_gbs = unpack_bytes / detail_max(SGS_TAU_UNPACK) / 1.0e9_rprec

tau_total = sum(sgs_tau_detail)
tau_bytes = real(sgs_tau_send_bytes + sgs_tau_recv_bytes, rprec)
tau_gbs = 0._rprec
if (sgs_tau_detail(SGS_TAU_MPI) > 0._rprec)                                  &
    tau_gbs = tau_bytes / sgs_tau_detail(SGS_TAU_MPI) / 1.0e9_rprec
#ifdef PPMPI
call mpi_allreduce(sgs_tau_arr_after_sync, tau_ready_min, 1, MPI_RPREC,     &
    MPI_MIN, comm, ierr)
call mpi_allreduce(sgs_tau_arr_after_sync, tau_ready_max, 1, MPI_RPREC,     &
    MPI_MAX, comm, ierr)
tau_ready_rel = 0._rprec
tau_skew = 0._rprec
if (tau_ready_max > 0._rprec) then
    tau_ready_rel = sgs_tau_arr_after_sync - tau_ready_min
    tau_skew = tau_ready_max - tau_ready_min
end if
tau_sync_to_mpi = 0._rprec
if (sgs_tau_arr_before_mpi > sgs_tau_arr_after_sync)                         &
    tau_sync_to_mpi = sgs_tau_arr_before_mpi - sgs_tau_arr_after_sync
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,'(a,i0,a,i0,6(a,es14.6),a,i0)')                              &
            'SGS_TAU_ARR call=', stage_count, ' rank=', coord,              &
            ' pre=', sgs_tau_detail(SGS_TAU_PRE_SYNC),                      &
            ' ready_rel=', tau_ready_rel,                                    &
            ' mpi=', sgs_tau_detail(SGS_TAU_MPI),                           &
            ' skew=', tau_skew,                                              &
            ' barrier=', sgs_tau_arr_barrier,                                &
            ' sync_to_mpi=', tau_sync_to_mpi,                                &
            ' bytes=', sgs_tau_send_bytes + sgs_tau_recv_bytes
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,*) 'SGS tau halo rank comm call=', stage_count, ' rank=',     &
            coord, ' pre=', sgs_tau_detail(SGS_TAU_PRE_SYNC),                &
            ' pack=', sgs_tau_detail(SGS_TAU_PACK),                          &
            ' mpi=', sgs_tau_detail(SGS_TAU_MPI),                            &
            ' unpack=', sgs_tau_detail(SGS_TAU_UNPACK),                      &
            ' post=', sgs_tau_detail(SGS_TAU_POST_SYNC),                     &
            ' d2h=', sgs_tau_detail(SGS_TAU_D2H),                            &
            ' h2d=', sgs_tau_detail(SGS_TAU_H2D),                            &
            ' total=', tau_total, ' mpi_calls=', sgs_tau_mpi_calls,          &
            ' send_bytes=', sgs_tau_send_bytes, ' recv_bytes=',              &
            sgs_tau_recv_bytes, ' neighbor=', sgs_tau_neighbor,              &
            ' combined=', sgs_tau_combined_msg, ' mpi_GBps=', tau_gbs
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
#else
write(*,*) 'SGS tau halo rank comm call=', stage_count, ' rank=', coord,      &
    ' pre=', sgs_tau_detail(SGS_TAU_PRE_SYNC),                               &
    ' pack=', sgs_tau_detail(SGS_TAU_PACK),                                  &
    ' mpi=', sgs_tau_detail(SGS_TAU_MPI),                                    &
    ' unpack=', sgs_tau_detail(SGS_TAU_UNPACK),                              &
    ' post=', sgs_tau_detail(SGS_TAU_POST_SYNC),                             &
    ' d2h=', sgs_tau_detail(SGS_TAU_D2H),                                    &
    ' h2d=', sgs_tau_detail(SGS_TAU_H2D),                                    &
    ' total=', tau_total, ' mpi_calls=', sgs_tau_mpi_calls,                  &
    ' send_bytes=', sgs_tau_send_bytes, ' recv_bytes=', sgs_tau_recv_bytes,  &
    ' neighbor=', sgs_tau_neighbor, ' combined=', sgs_tau_combined_msg,       &
    ' mpi_GBps=', tau_gbs
flush(6)
#endif

if (coord == 0) then
    write(*,'(a,i8)') 'SGS tau halo breakdown (max rank), call ', stage_count
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  pre-MPI sync: time=',            &
        detail_max(SGS_TAU_PRE_SYNC), ' bytes=', 0._rprec, ' GB/s=',         &
        0._rprec
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  pack: time=',                    &
        detail_max(SGS_TAU_PACK), ' bytes=', pack_bytes, ' GB/s=', pack_gbs
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  MPI exchange/wait: time=',       &
        detail_max(SGS_TAU_MPI), ' bytes=', mpi_bytes, ' GB/s=', mpi_gbs
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  unpack: time=',                  &
        detail_max(SGS_TAU_UNPACK), ' bytes=', unpack_bytes, ' GB/s=',       &
        unpack_gbs
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  post-MPI sync/check: time=',     &
        detail_max(SGS_TAU_POST_SYNC), ' bytes=', 0._rprec, ' GB/s=',        &
        0._rprec
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  host D2H copy: time=',           &
        detail_max(SGS_TAU_D2H), ' bytes=', pack_bytes, ' GB/s=', 0._rprec
    write(*,'(a,E15.7,a,E15.7,a,E15.7)') '  host H2D copy: time=',           &
        detail_max(SGS_TAU_H2D), ' bytes=', unpack_bytes, ' GB/s=', 0._rprec
end if

end subroutine sgs_tau_halo_detail_report
#endif

!*******************************************************************************
subroutine sgs_stage_report(stage_count, calc_sij, dynamic_model, nut,         &
    tau_boundary, tau_interior, tau_halo, final_sync, gpu_calc_sij,           &
    gpu_dynamic_model, gpu_nut, gpu_tau_boundary, gpu_tau_interior)
!*******************************************************************************
use types, only : rprec
use param, only : coord
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC, nproc
use mpi
#endif
implicit none

integer, intent(in) :: stage_count
real(rprec), intent(in) :: calc_sij, dynamic_model, nut, tau_boundary
real(rprec), intent(in) :: tau_interior, tau_halo, final_sync
real(rprec), intent(in) :: gpu_calc_sij, gpu_dynamic_model, gpu_nut
real(rprec), intent(in) :: gpu_tau_boundary, gpu_tau_interior
real(rprec) :: calc_max, dyn_max, nut_max, bnd_max, int_max, halo_max
real(rprec) :: final_max, total
real(rprec) :: gpu_calc_max, gpu_dyn_max, gpu_nut_max, gpu_bnd_max
real(rprec) :: gpu_int_max
real(rprec) :: queued_gpu
#ifdef PPMPI
integer :: r
#endif

#ifdef PPMPI
call mpi_allreduce(calc_sij, calc_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(dynamic_model, dyn_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(nut, nut_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(tau_boundary, bnd_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(tau_interior, int_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(tau_halo, halo_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(final_sync, final_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(gpu_calc_sij, gpu_calc_max, 1, MPI_RPREC, MPI_MAX, comm,   &
    ierr)
call mpi_allreduce(gpu_dynamic_model, gpu_dyn_max, 1, MPI_RPREC, MPI_MAX,    &
    comm, ierr)
call mpi_allreduce(gpu_nut, gpu_nut_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(gpu_tau_boundary, gpu_bnd_max, 1, MPI_RPREC, MPI_MAX,     &
    comm, ierr)
call mpi_allreduce(gpu_tau_interior, gpu_int_max, 1, MPI_RPREC, MPI_MAX,     &
    comm, ierr)
#else
calc_max = calc_sij
dyn_max = dynamic_model
nut_max = nut
bnd_max = tau_boundary
int_max = tau_interior
halo_max = tau_halo
final_max = final_sync
gpu_calc_max = gpu_calc_sij
gpu_dyn_max = gpu_dynamic_model
gpu_nut_max = gpu_nut
gpu_bnd_max = gpu_tau_boundary
gpu_int_max = gpu_tau_interior
#endif

total = calc_max + dyn_max + nut_max + bnd_max + int_max + halo_max +         &
    final_max
queued_gpu = gpu_calc_sij + gpu_dynamic_model + gpu_nut + gpu_tau_boundary + &
    gpu_tau_interior

#ifdef PPMPI
do r = 0, nproc - 1
    call mpi_barrier(comm, ierr)
    if (coord == r) then
        write(*,'(a,i0,a,i0,6(a,es14.6))') 'SGS_QUEUE_GPU call=',            &
            stage_count, ' rank=', coord, ' calc=', gpu_calc_sij,            &
            ' dynamic=', gpu_dynamic_model, ' nut=', gpu_nut,                &
            ' tau_bnd=', gpu_tau_boundary, ' tau_int=', gpu_tau_interior,    &
            ' total=', queued_gpu
        flush(6)
    end if
end do
call mpi_barrier(comm, ierr)
#else
write(*,'(a,i0,a,i0,6(a,es14.6))') 'SGS_QUEUE_GPU call=', stage_count,       &
    ' rank=', coord, ' calc=', gpu_calc_sij, ' dynamic=', gpu_dynamic_model, &
    ' nut=', gpu_nut, ' tau_bnd=', gpu_tau_boundary,                         &
    ' tau_int=', gpu_tau_interior, ' total=', queued_gpu
flush(6)
#endif

if (coord == 0) then
    write(*,'(a,i8)') 'SGS stage timing (max rank), call ', stage_count
    write(*,'(1a,E15.7)') '  calc_Sij: ', calc_max
    write(*,'(1a,E15.7)') '  dynamic model/Cs: ', dyn_max
    write(*,'(1a,E15.7)') '  Nu_t: ', nut_max
    write(*,'(1a,E15.7)') '  tau boundary planes: ', bnd_max
    write(*,'(1a,E15.7)') '  tau interior: ', int_max
    write(*,'(1a,E15.7)') '  tau halo pack/MPI/unpack: ', halo_max
    write(*,'(1a,E15.7)') '  safety/final sync: ', final_max
    write(*,'(1a,E15.7)') '  SGS stage sum: ', total
    write(*,'(a,i8)') 'SGS GPU event timing (max rank), call ', stage_count
    write(*,'(1a,E15.7)') '  calc_Sij GPU event: ', gpu_calc_max
    write(*,'(1a,E15.7)') '  dynamic model/Cs GPU event: ', gpu_dyn_max
    write(*,'(1a,E15.7)') '  Nu_t GPU event: ', gpu_nut_max
    write(*,'(1a,E15.7)') '  tau boundary planes GPU event: ', gpu_bnd_max
    write(*,'(1a,E15.7)') '  tau interior GPU event: ', gpu_int_max
end if

#ifdef ENABLE_CUDA
call sgs_calc_sij_detail_report(stage_count)
call sgs_tau_halo_detail_report(stage_count)
#endif

end subroutine sgs_stage_report

!*******************************************************************************
real(rprec) function rtnewt(A, jz)
!*******************************************************************************
use types, only : rprec
integer, parameter :: jmax=100
real(rprec) :: x1, x2, xacc
real(rprec) :: df,dx,f
integer :: j, jz
real(rprec), dimension(0:5) :: A
x1 = 0._rprec
x2 = 15._rprec                  ! try to find the largest root first
xacc = 0.001_rprec              ! doesn't need to be that accurate
rtnewt = 0.5_rprec*(x1+x2)
do j = 1, jmax
    f = A(0)+rtnewt*(A(1)+rtnewt*(A(2)+rtnewt*(A(3)+rtnewt*(A(4)+rtnewt*A(5)))))
    df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) +                &
        rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))
    dx = f/df
    rtnewt = rtnewt - dx
    if (abs(dx) < xacc) return
end do
rtnewt = 1._rprec  ! if don't converge fast enough
write(6,*) 'using beta=1 at jz= ', jz

end function rtnewt

end module sgs_stag_util
