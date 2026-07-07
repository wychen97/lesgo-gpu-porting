#if !defined(PPPRESS_GPU) || (defined(PPLVLSET) && defined(PPLES_GPU))
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
use sim_param, only : u, v, w, divtz, p, dpdx, dpdy, dpdz
use fft, only : back, forw, kx, ky
use emul_complex, only : OPERATOR(.MULI.)
#ifdef PPMPI
use cuda_mpi_debug, only : mpi_dbg_sendrecv_r, mpi_dbg_send_r, mpi_dbg_recv_r
use mpi
#endif

implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz
integer :: ir, ii
integer :: jz_min
! CPU build: plain host versions of the work arrays. The
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
integer, save :: press_fw_plan_nz1 = 0
integer, save :: press_fw_plan_1 = 0
integer, save :: press_bk_plan_nz1 = 0
integer, save :: press_bk_plan_nzp1 = 0
! Legacy cuf-path feature toggles: default OFF in maintained builds.
! In a plain CPU build the cuf bodies of these branches are preprocessed
! out, so taking them would SKIP real work (e.g. the p(:,:,0) halo exchange,
! posted inside the cuf tridag) -> defaults to .false. so every branch falls
! through to the original blocking CPU path.
logical, parameter :: press_cuda_paths = .false.
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

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( rH_x(ld,ny,lbz:nz), rH_y(ld,ny,lbz:nz), rH_z(ld,ny,lbz:nz) )
    allocate ( rtopw(ld,ny), rbottomw(ld,ny) )
    allocate ( RHS_col(ld,ny,nz+1) )
    allocate ( a(lh,ny,nz+1), b(lh,ny,nz+1), c(lh,ny,nz+1) )

    arrays_allocated = .true.
endif


    if (coord == 0) then
        p(:,:,0) = 0._rprec
#ifdef PPSAFETYMODE
    else
        p(:,:,0) = BOGUS
#endif
    end if


! Get the right hand side ready
! Loop over levels
! Recall that the old timestep guys already contain the pressure
do jz = 1, nz-1
    rH_x(:,:,jz) = const2 * u(:,:,jz)
    rH_y(:,:,jz) = const2 * v(:,:,jz)
    rH_z(:,:,jz) = const2 * w(:,:,jz)

    call dfftw_execute_dft_r2c(forw, rH_x(:,:,jz), rH_x(:,:,jz))
    call dfftw_execute_dft_r2c(forw, rH_y(:,:,jz), rH_y(:,:,jz))
    call dfftw_execute_dft_r2c(forw, rH_z(:,:,jz), rH_z(:,:,jz))
end do

#if defined(PPMPI) && defined(PPSAFETYMODE)
  !Careful - only update real values (odd indicies)
  rH_x(1:ld:2,:,0) = BOGUS
  rH_y(1:ld:2,:,0) = BOGUS
  rH_z(1:ld:2,:,0) = BOGUS
#endif

#ifdef PPSAFETYMODE
!Careful - only update real values (odd indicies)
rH_x(1:ld:2,:,nz) = BOGUS
rH_y(1:ld:2,:,nz) = BOGUS
#endif

#ifdef PPMPI
if (coord == nproc-1) then
    do jy = 1, ny
    do jx = 1, ld
        rH_z(jx,jy,nz) = const2 * w(jx,jy,nz)
    end do
    end do
        call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,jz))
#ifdef PPSAFETYMODE
else
    rH_z(1:ld:2,:,nz) = BOGUS !--perhaps this should be 0 on top process?
#endif
endif
#else
do jy = 1, ny
do jx = 1, ld
    rH_z(jx,jy,nz) = const2 * w(jx,jy,nz)
end do
end do
call dfftw_execute_dft_r2c(forw, rH_z(:,:,nz), rH_z(:,:,jz))
#endif

if (coord == 0) then
    do jy = 1, ny
    do jx = 1, ld
        rbottomw(jx,jy) = const * divtz(jx,jy,1)
    end do
    end do
        call dfftw_execute_dft_r2c(forw, rbottomw, rbottomw )
end if

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    do jy = 1, ny
    do jx = 1, ld
        rtopw(jx,jy) = const * divtz(jx,jy,nz)
    end do
    end do
        call dfftw_execute_dft_r2c(forw, rtopw, rtopw)
#ifdef PPMPI
endif
#endif

if (press_sync_after_forward_fft_enabled) then
    press_q_t0 = press_queue_wtime()
    call press_cuda_sync('pressure forced sync after forward FFT')
    press_q_sync_after_forward = press_queue_wtime() - press_q_t0
    press_q_sync_count = press_q_sync_count + 1
end if





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


! zero-wavenumber solution
#ifdef PPMPI
! wait for p(1, 1, 1) from "down"
if (coord > 0) then
    call mpi_dbg_recv_r (p(1, 1, 1), 2, MPI_RPREC, down, 8, comm, status,     &
        ierr, 'press_p_zero_recv')
end if
#endif


if (coord == 0) then
    p(1:2, 1, 0) = 0._rprec
    p(1:2, 1, 1) = p(1:2,1,0) - dz * rbottomw(1:2,1)
end if

do jz = 2, nz
    ! JDA dissertation, eqn(2.88)
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z(1:2, 1, jz) * dz
end do


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

#endif
