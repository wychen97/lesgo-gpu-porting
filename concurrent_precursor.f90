!!
!!  Copyright (C) 2011-2020  Johns Hopkins University
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
module concurrent_precursor
!*******************************************************************************
use types, only : rprec
use mpi_defs
use fringe
#ifdef PPLES_GPU
use openacc
#endif
implicit none

save
private

public :: vel_sample_t
public :: initialize_cps, synchronize_cps, inflow_cps

character(*), parameter :: mod_name = 'concurrent_precursor'

type vel_sample_type
    #ifdef ENABLE_CUDA
    real(rprec), managed, allocatable, dimension(:,:,:) :: u, v, w
    #else
    real(rprec), allocatable, dimension(:,:,:) :: u, v, w
    #endif
#ifdef PPSCALARS
    #ifdef ENABLE_CUDA
    real(rprec), managed, allocatable, dimension(:,:,:) :: theta
    #else
    real(rprec), allocatable, dimension(:,:,:) :: theta
    #endif
#endif
end type vel_sample_type

type(vel_sample_type), target :: vel_sample_t
type(fringe_t), target :: cps_fringe
#ifdef PPLES_GPU
integer, allocatable, dimension(:) :: cps_iwrap_acc
real(rprec), allocatable, dimension(:) :: cps_alpha_acc, cps_beta_acc
#endif
integer, save :: cps_stage_count = 0
real(rprec), save :: cps_time_sync = 0._rprec
real(rprec), save :: cps_time_inflow = 0._rprec
real(rprec), save :: cps_time_sync_red_sample = 0._rprec
real(rprec), save :: cps_time_sync_red_send = 0._rprec
real(rprec), save :: cps_time_sync_blue_recv = 0._rprec
real(rprec), save :: cps_time_sync_coriolis = 0._rprec

contains

!*******************************************************************************
logical function cps_env_true_token_enabled(name)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
character(len=32) :: setting
integer :: stat

cps_env_true_token_enabled = .false.
call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('1', 'true', 'TRUE', 'True', 'on', 'ON', 'On', 'yes', 'YES', 'Yes')
        cps_env_true_token_enabled = .true.
    case default
        cps_env_true_token_enabled = .false.
    end select
end if

end function cps_env_true_token_enabled

!*******************************************************************************
logical function cps_stage_timing_enabled()
!*******************************************************************************
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = cps_env_true_token_enabled('LESGO_CPS_STAGE_TIMING')
    initialized = .true.
end if

cps_stage_timing_enabled = enabled

end function cps_stage_timing_enabled

!*******************************************************************************
subroutine cps_timer_start(t0)
!*******************************************************************************
implicit none

real(rprec), intent(out) :: t0

#ifdef PPLES_GPU
!$acc wait
#endif
call cpu_time(t0)

end subroutine cps_timer_start

!*******************************************************************************
subroutine cps_timer_accum(t0, accum)
!*******************************************************************************
implicit none

real(rprec), intent(in) :: t0
real(rprec), intent(inout) :: accum
real(rprec) :: t1

#ifdef PPLES_GPU
!$acc wait
#endif
call cpu_time(t1)
accum = accum + max(t1 - t0, 0._rprec)

end subroutine cps_timer_accum

!*******************************************************************************
subroutine cps_stage_report(coord_in, wbase_in)
!*******************************************************************************
implicit none

integer, intent(in) :: coord_in, wbase_in
character(len=4) :: color_name

if (coord_in /= 0) return
if (wbase_in > 0) then
    if (mod(cps_stage_count, wbase_in) /= 0) return
end if

if (color == RED) then
    color_name = 'RED'
else if (color == BLUE) then
    color_name = 'BLUE'
else
    color_name = 'UNK'
end if

write(*,'(a,a,a,i8)') 'CPS stage timing (', trim(color_name),                 &
    '), call ', cps_stage_count
write(*,'(1a,E15.7)') '  synchronize_cps: ', cps_time_sync
write(*,'(1a,E15.7)') '    sync red sample: ', cps_time_sync_red_sample
write(*,'(1a,E15.7)') '    sync red send: ', cps_time_sync_red_send
write(*,'(1a,E15.7)') '    sync blue recv: ', cps_time_sync_blue_recv
write(*,'(1a,E15.7)') '    sync coriolis: ', cps_time_sync_coriolis
write(*,'(1a,E15.7)') '  inflow_cps: ', cps_time_inflow
write(*,'(1a,E15.7)') '  cps total: ', cps_time_sync + cps_time_inflow

end subroutine cps_stage_report

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function cps_cuda_enabled()
!*******************************************************************************
implicit none

cps_cuda_enabled = .true.

end function cps_cuda_enabled
#endif

    !*******************************************************************************
subroutine initialize_cps()
!*******************************************************************************
use param, only : nx, ny, nz, dx, L_x, coord, rank_of_coord, status, ierr
use param, only : fringe_region_end, fringe_region_len, sampling_region_end
use messages
use mpi
implicit none

character (*), parameter :: sub_name = mod_name // '.initialize_cps'
integer :: i

if( color == BLUE ) then
    cps_fringe = fringe_t(fringe_region_end, fringe_region_len)
    call mpi_send(cps_fringe%nx , 1, MPI_INTEGER,                              &
        rank_of_coord(coord), 1, interComm, ierr )
elseif( color == RED) then
    call mpi_recv(i , 1, MPI_INTEGER,                                          &
        rank_of_coord(coord), 1, interComm, status, ierr)
    fringe_region_len = (i - 0.5_rprec)*dx/L_x
    cps_fringe = fringe_t(sampling_region_end, fringe_region_len)
else
   call error(sub_name, 'Erroneous color specification')
endif

! Allocate the sample block
allocate(vel_sample_t%u(cps_fringe%nx, ny, nz ))
allocate(vel_sample_t%v(cps_fringe%nx, ny, nz ))
allocate(vel_sample_t%w(cps_fringe%nx, ny, nz ))
#ifdef PPSCALARS
allocate(vel_sample_t%theta(cps_fringe%nx, ny, nz))
#endif

#ifdef PPLES_GPU
! Active LES GPU builds use OpenACC explicit residency rather than the older
! ENABLE_CUDA managed-memory branch.  Keep CPS sample buffers and fringe
! metadata resident so sampling, MPI exchange, and fringe application do not
! force velocity host round-trips.
allocate(cps_iwrap_acc(cps_fringe%nx))
allocate(cps_alpha_acc(cps_fringe%nx))
allocate(cps_beta_acc(cps_fringe%nx))
cps_iwrap_acc = cps_fringe%iwrap
cps_alpha_acc = cps_fringe%alpha
cps_beta_acc = cps_fringe%beta
!$acc enter data copyin(cps_iwrap_acc, cps_alpha_acc, cps_beta_acc)
!$acc enter data create(vel_sample_t%u(1:cps_fringe%nx,1:ny,1:nz),            &
!$acc                   vel_sample_t%v(1:cps_fringe%nx,1:ny,1:nz),            &
!$acc                   vel_sample_t%w(1:cps_fringe%nx,1:ny,1:nz))
#ifdef PPSCALARS_GPU
!$acc enter data create(vel_sample_t%theta(1:cps_fringe%nx,1:ny,1:nz))
#endif
#endif

end subroutine initialize_cps

!*******************************************************************************
subroutine synchronize_cps()
!*******************************************************************************
use types, only : rprec
use messages
use param, only : ny, nz
use param, only : coord, rank_of_coord, status, ierr, MPI_RPREC
use sim_param, only : u,v,w
#ifdef PPSCALARS
use scalars, only : theta
#endif
use coriolis, only : coriolis_forcing, alpha, G
implicit none

character (*), parameter :: sub_name = mod_name // '.synchronize_cps'
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
#ifdef PPSCALARS
real(rprec), pointer, dimension(:,:,:) :: theta_p
#endif
integer :: sendsize, recvsize
integer :: i, j, k
logical :: cps_timing
real(rprec) :: cps_t0, cps_t_sub

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif

u_p => vel_sample_t%u
v_p => vel_sample_t%v
w_p => vel_sample_t%w
#ifdef PPSCALARS
theta_p => vel_sample_t%theta
#endif

sendsize = cps_fringe%nx * ny * nz
recvsize = sendsize
cps_timing = cps_stage_timing_enabled()
if (cps_timing) call cps_timer_start(cps_t0)

#ifdef PPLES_GPU
if( color == BLUE ) then
    ! Receive sampled velocities directly into device-resident CPS buffers.
    if (cps_timing) call cps_timer_start(cps_t_sub)
    !$acc host_data use_device(u_p)
    call mpi_recv( u_p(1,1,1) , recvsize, MPI_RPREC,                          &
        rank_of_coord(coord), 1, interComm, status, ierr)
    !$acc end host_data
    !$acc host_data use_device(v_p)
    call mpi_recv( v_p(1,1,1) , recvsize, MPI_RPREC,                          &
        rank_of_coord(coord), 2, interComm, status, ierr)
    !$acc end host_data
    !$acc host_data use_device(w_p)
    call mpi_recv( w_p(1,1,1) , recvsize, MPI_RPREC,                          &
        rank_of_coord(coord), 3, interComm, status, ierr)
    !$acc end host_data
#ifdef PPSCALARS_GPU
    !$acc host_data use_device(theta_p)
    call mpi_recv( theta_p(1,1,1) , recvsize, MPI_RPREC,                      &
        rank_of_coord(coord), 4, interComm, status, ierr)
    !$acc end host_data
#elif defined(PPSCALARS)
    call mpi_recv( theta_p(1,1,1) , recvsize, MPI_RPREC,                      &
        rank_of_coord(coord), 4, interComm, status, ierr)
#endif
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_blue_recv)
if (coriolis_forcing>0) then
    if (cps_timing) call cps_timer_start(cps_t_sub)
    call mpi_recv(G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm,        &
        status, ierr)
    call mpi_recv(alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm,    &
        status, ierr)
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_coriolis)
end if

elseif( color == RED ) then
    ! Sample the active device-resident velocity field into device CPS buffers.
    if (cps_timing) call cps_timer_start(cps_t_sub)
    !$acc parallel loop collapse(3) present(u_p, v_p, w_p, u, v, w, cps_iwrap_acc) async(1)
    do k = 1, nz
    do j = 1, ny
    do i = 1, cps_fringe%nx
        u_p(i,j,k) = u(cps_iwrap_acc(i),j,k)
        v_p(i,j,k) = v(cps_iwrap_acc(i),j,k)
        w_p(i,j,k) = w(cps_iwrap_acc(i),j,k)
    end do
    end do
    end do
#ifdef PPSCALARS_GPU
    !$acc parallel loop collapse(3) present(theta_p, theta, cps_iwrap_acc) async(1)
    do k = 1, nz
    do j = 1, ny
    do i = 1, cps_fringe%nx
        theta_p(i,j,k) = theta(cps_iwrap_acc(i),j,k)
    end do
    end do
    end do
#elif defined(PPSCALARS)
    theta_p(:,:,:) = theta(cps_fringe%iwrap(:),1:ny,1:nz)
#endif
    !$acc wait(1)
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_red_sample)

    if (cps_timing) call cps_timer_start(cps_t_sub)
    !$acc host_data use_device(u_p)
    call mpi_send( u_p(1,1,1), sendsize, MPI_RPREC,                           &
        rank_of_coord(coord), 1, interComm, ierr )
    !$acc end host_data
    !$acc host_data use_device(v_p)
    call mpi_send( v_p(1,1,1), sendsize, MPI_RPREC,                           &
        rank_of_coord(coord), 2, interComm, ierr )
    !$acc end host_data
    !$acc host_data use_device(w_p)
    call mpi_send( w_p(1,1,1), sendsize, MPI_RPREC,                           &
        rank_of_coord(coord), 3, interComm, ierr )
    !$acc end host_data
#ifdef PPSCALARS_GPU
    !$acc host_data use_device(theta_p)
    call mpi_send( theta_p(1,1,1), sendsize, MPI_RPREC,                       &
        rank_of_coord(coord), 4, interComm, ierr )
    !$acc end host_data
#elif defined(PPSCALARS)
    call mpi_send( theta_p(1,1,1), sendsize, MPI_RPREC,                       &
        rank_of_coord(coord), 4, interComm, ierr )
#endif
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_red_send)
if (coriolis_forcing>0) then
    if (cps_timing) call cps_timer_start(cps_t_sub)
    call mpi_send(G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm, ierr)
    call mpi_send(alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm, ierr)
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_coriolis)
end if

else
   call error(sub_name, 'Erroneous color specification')
endif

nullify(u_p, v_p, w_p)
#ifdef PPSCALARS
nullify(theta_p)
#endif
if (cps_timing) call cps_timer_accum(cps_t0, cps_time_sync)
return
#endif

if( color == BLUE ) then
    ! Recieve sampled velocities from upstream (RED)
    if (cps_timing) call cps_timer_start(cps_t_sub)
    call mpi_recv( u_p(1,1,1) , recvsize, MPI_RPREC,                           &
        rank_of_coord(coord), 1, interComm, status, ierr)
    call mpi_recv( v_p(1,1,1) , recvsize, MPI_RPREC,                           &
        rank_of_coord(coord), 2, interComm, status, ierr)
    call mpi_recv( w_p(1,1,1) , recvsize, MPI_RPREC,                           &
        rank_of_coord(coord), 3, interComm, status, ierr)
#ifdef PPSCALARS
    call mpi_recv( theta_p(1,1,1) , recvsize, MPI_RPREC,                       &
        rank_of_coord(coord), 4, interComm, status, ierr)
#endif
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_blue_recv)
if (coriolis_forcing>0) then
    if (cps_timing) call cps_timer_start(cps_t_sub)
    call mpi_recv(G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm,         &
        status, ierr)
    call mpi_recv(alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm,     &
        status, ierr)
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_coriolis)
end if

elseif( color == RED ) then
    ! Sample velocity and copy to buffers
    if (cps_timing) call cps_timer_start(cps_t_sub)
    u_p(:,:,:) = u(cps_fringe%iwrap(:),1:ny,1:nz)
    v_p(:,:,:) = v(cps_fringe%iwrap(:),1:ny,1:nz)
    w_p(:,:,:) = w(cps_fringe%iwrap(:),1:ny,1:nz)
#ifdef PPSCALARS
    theta_p(:,:,:) = theta(cps_fringe%iwrap(:),1:ny,1:nz)
#endif
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_red_sample)

    ! Send sampled velocities to downstream domain (BLUE)
    if (cps_timing) call cps_timer_start(cps_t_sub)
    call mpi_send( u_p(1,1,1), sendsize, MPI_RPREC,                            &
        rank_of_coord(coord), 1, interComm, ierr )
    call mpi_send( v_p(1,1,1), sendsize, MPI_RPREC,                            &
        rank_of_coord(coord), 2, interComm, ierr )
    call mpi_send( w_p(1,1,1), sendsize, MPI_RPREC,                            &
        rank_of_coord(coord), 3, interComm, ierr )
#ifdef PPSCALARS
    call mpi_send( theta_p(1,1,1), sendsize, MPI_RPREC,                        &
        rank_of_coord(coord), 4, interComm, ierr )
#endif
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_red_send)
if (coriolis_forcing>0) then
    if (cps_timing) call cps_timer_start(cps_t_sub)
    call mpi_send(G, 1, MPI_RPREC, rank_of_coord(coord), 5, interComm, ierr)
    call mpi_send(alpha, 1, MPI_RPREC, rank_of_coord(coord), 6, interComm, ierr)
    if (cps_timing) call cps_timer_accum(cps_t_sub, cps_time_sync_coriolis)
end if

else
   call error(sub_name, 'Erroneous color specification')
endif

nullify(u_p, v_p, w_p)
#ifdef PPSCALARS
nullify(theta_p)
#endif
if (cps_timing) call cps_timer_accum(cps_t0, cps_time_sync)

end subroutine synchronize_cps

!*******************************************************************************
subroutine inflow_cps ()
!*******************************************************************************
!
!  Enforces prescribed inflow condition from an inlet velocity field
!  generated from a precursor simulation. The inflow condition is
!  enforced by direct modulation on the velocity in the fringe region.
!
use types, only : rprec
use param, only : nx, ny, nz, coord, wbase
use sim_param, only : u, v, w
#ifdef PPSCALARS
use scalars, only : theta
#endif
use messages, only : error
implicit none

character (*), parameter :: sub_name = 'inflow_cond_cps'
    integer :: i, i_w, j, k
real(rprec), pointer, dimension(:,:,:) :: u_p, v_p, w_p
#ifdef PPSCALARS
real(rprec), pointer, dimension(:,:,:) :: theta_p
#endif
logical :: cps_timing
real(rprec) :: cps_t0

nullify( u_p, v_p, w_p )
#ifdef PPSCALARS
nullify( theta_p )
#endif

u_p => vel_sample_t%u
v_p => vel_sample_t%v
w_p => vel_sample_t%w
#ifdef PPSCALARS
theta_p => vel_sample_t%theta
#endif
cps_timing = cps_stage_timing_enabled()
if (cps_timing) call cps_timer_start(cps_t0)

#ifdef PPLES_GPU
! Active LES GPU path: velocity sample buffers are already device-resident from
! synchronize_cps().  Apply the fringe on device.  Scalar CPU fallback remains
! host-side unless PPSCALARS_GPU is also enabled.
#ifdef PPSCALARS_GPU
!$acc parallel loop collapse(3) present(u, v, w, theta, u_p, v_p, w_p, theta_p, cps_iwrap_acc, cps_alpha_acc, cps_beta_acc) async(1)
#else
!$acc parallel loop collapse(3) present(u, v, w, u_p, v_p, w_p, cps_iwrap_acc, cps_alpha_acc, cps_beta_acc) async(1)
#endif
do k = 1, nz
do j = 1, ny
do i = 1, cps_fringe%nx
    i_w = cps_iwrap_acc(i)
    u(i_w,j,k) = cps_alpha_acc(i) * u(i_w,j,k)                                &
        + cps_beta_acc(i) * u_p(i,j,k)
    v(i_w,j,k) = cps_alpha_acc(i) * v(i_w,j,k)                                &
        + cps_beta_acc(i) * v_p(i,j,k)
    w(i_w,j,k) = cps_alpha_acc(i) * w(i_w,j,k)                                &
        + cps_beta_acc(i) * w_p(i,j,k)
#ifdef PPSCALARS_GPU
    theta(i_w,j,k) = cps_alpha_acc(i) * theta(i_w,j,k)                        &
        + cps_beta_acc(i) * theta_p(i,j,k)
#endif
end do
end do
end do
!$acc wait(1)

#if defined(PPSCALARS) && !defined(PPSCALARS_GPU)
do i = 1, cps_fringe%nx
    i_w = cps_fringe%iwrap(i)
    theta(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * theta(i_w,1:ny,1:nz)         &
        + cps_fringe%beta(i) * theta_p(i,1:ny,1:nz)
end do
#endif

nullify(u_p, v_p, w_p)
#ifdef PPSCALARS
nullify(theta_p)
#endif
if (cps_timing) then
    call cps_timer_accum(cps_t0, cps_time_inflow)
    cps_stage_count = cps_stage_count + 1
    call cps_stage_report(coord, wbase)
end if
return
#endif

#ifdef ENABLE_CUDA
    if (cps_cuda_enabled()) then
        !$cuf kernel do(3) <<<*,*>>>
        do k = 1, nz
        do j = 1, ny
        do i = 1, cps_fringe%nx
            i_w = cps_fringe%iwrap(i)
            u(i_w,j,k) = cps_fringe%alpha(i) * u(i_w,j,k)                         &
                + cps_fringe%beta(i) * vel_sample_t%u(i,j,k)
            v(i_w,j,k) = cps_fringe%alpha(i) * v(i_w,j,k)                         &
                + cps_fringe%beta(i) * vel_sample_t%v(i,j,k)
            w(i_w,j,k) = cps_fringe%alpha(i) * w(i_w,j,k)                         &
                + cps_fringe%beta(i) * vel_sample_t%w(i,j,k)
#ifdef PPSCALARS
            theta(i_w,j,k) = cps_fringe%alpha(i) * theta(i_w,j,k)                 &
                + cps_fringe%beta(i) * vel_sample_t%theta(i,j,k)
#endif
        end do
        end do
        end do
    else
#endif
    do i = 1, cps_fringe%nx
        i_w = cps_fringe%iwrap(i)
        u(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * u(i_w,1:ny,1:nz)                  &
            + cps_fringe%beta(i) * u_p(i,1:ny,1:nz)
        v(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * v(i_w,1:ny,1:nz)                  &
            + cps_fringe%beta(i) * v_p(i,1:ny,1:nz)
        w(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * w(i_w,1:ny,1:nz)                  &
            + cps_fringe%beta(i) * w_p(i,1:ny,1:nz)
#ifdef PPSCALARS
        theta(i_w,1:ny,1:nz) = cps_fringe%alpha(i) * theta(i_w,1:ny,1:nz)          &
            + cps_fringe%beta(i) * theta_p(i,1:ny,1:nz)
#endif
    end do
#ifdef ENABLE_CUDA
    endif
#endif

nullify(u_p, v_p, w_p)
#ifdef PPSCALARS
nullify(theta_p)
#endif
if (cps_timing) then
    call cps_timer_accum(cps_t0, cps_time_inflow)
    cps_stage_count = cps_stage_count + 1
    call cps_stage_report(coord, wbase)
end if

end subroutine inflow_cps

end module concurrent_precursor
