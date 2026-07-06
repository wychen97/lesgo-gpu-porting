!!
!!  Copyright (C) 2020  Johns Hopkins University
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
module shifted_inflow
!*******************************************************************************
use types, only : rprec
use fringe

implicit none

private
public shifted_inflow_init, inflow_shifted

type(fringe_t) :: sample_fringe
type(fringe_t) :: apply_fringe

#if defined(PPLES_GPU)
integer, allocatable, dimension(:) :: sample_iwrap_cuda
integer, allocatable, dimension(:) :: apply_iwrap_cuda
real(rprec), allocatable, dimension(:) :: apply_alpha_cuda
real(rprec), allocatable, dimension(:) :: apply_beta_cuda
real(rprec), allocatable, dimension(:,:,:) :: u_s, v_s, w_s
!$acc declare create(sample_iwrap_cuda, apply_iwrap_cuda)
!$acc declare create(apply_alpha_cuda, apply_beta_cuda)
!$acc declare create(u_s, v_s, w_s)
#else
real(rprec), allocatable, dimension(:,:,:) :: u_s, v_s, w_s
#endif
#ifdef PPSCALARS
real(rprec), allocatable, dimension(:,:,:) :: theta_s
#ifdef PPLES_GPU
!$acc declare create(theta_s)
#endif
#endif

contains

#if defined(PPLES_GPU)
!*******************************************************************************
logical function shifted_inflow_cuda_enabled()
!*******************************************************************************
implicit none

shifted_inflow_cuda_enabled = .true.

end function shifted_inflow_cuda_enabled

#endif

!*******************************************************************************
subroutine shifted_inflow_init
!*******************************************************************************
use param, only : fringe_region_end, fringe_region_len, sampling_region_end
use param, only : shift_n, ny, nz

sample_fringe = fringe_t(sampling_region_end, fringe_region_len)
apply_fringe = fringe_t(fringe_region_end, fringe_region_len)

#if defined(PPLES_GPU)
if (allocated(sample_iwrap_cuda)) deallocate(sample_iwrap_cuda)
if (allocated(apply_iwrap_cuda)) deallocate(apply_iwrap_cuda)
if (allocated(apply_alpha_cuda)) deallocate(apply_alpha_cuda)
if (allocated(apply_beta_cuda)) deallocate(apply_beta_cuda)
allocate(sample_iwrap_cuda(sample_fringe%nx))
allocate(apply_iwrap_cuda(apply_fringe%nx))
allocate(apply_alpha_cuda(apply_fringe%nx))
allocate(apply_beta_cuda(apply_fringe%nx))
sample_iwrap_cuda = sample_fringe%iwrap
apply_iwrap_cuda = apply_fringe%iwrap
apply_alpha_cuda = apply_fringe%alpha
apply_beta_cuda = apply_fringe%beta
#ifdef PPLES_GPU
!$acc update device(sample_iwrap_cuda, apply_iwrap_cuda, apply_alpha_cuda,      &
!$acc& apply_beta_cuda)
#endif
#endif

! Allocate the sample block
allocate(u_s(sample_fringe%nx, ny, nz ))
allocate(v_s(sample_fringe%nx, ny, nz ))
allocate(w_s(sample_fringe%nx, ny, nz ))
#ifdef PPSCALARS
allocate(theta_s(sample_fringe%nx, ny, nz))
#endif

! Only allow positive shifts than are less than ny/2
shift_n = modulo(abs(shift_n), ny)
if (shift_n > ny/2) shift_n = ny-shift_n
if (shift_n == 0) shift_n = 1

end subroutine shifted_inflow_init

!*******************************************************************************
subroutine inflow_shifted
!*******************************************************************************
use param, only : shift_n, ny, nz
use sim_param, only : u, v, w
#ifdef PPSCALARS
use scalars, only : theta
#endif
integer :: i, i_w, sample_nx, apply_nx
#if defined(PPLES_GPU)
integer :: j, k
#endif

sample_nx = sample_fringe%nx
apply_nx = apply_fringe%nx

#ifdef PPLES_GPU
if (shifted_inflow_cuda_enabled()) then
    !$acc parallel loop collapse(3) default(present)
    do k = 1, nz
    do j = shift_n+1, ny
    do i = 1, sample_nx
        u_s(i,j,k) = u(sample_iwrap_cuda(i),j-shift_n,k)
        v_s(i,j,k) = v(sample_iwrap_cuda(i),j-shift_n,k)
        w_s(i,j,k) = w(sample_iwrap_cuda(i),j-shift_n,k)
#ifdef PPSCALARS
        theta_s(i,j,k) = theta(sample_iwrap_cuda(i),j-shift_n,k)
#endif
    end do
    end do
    end do

    !$acc parallel loop collapse(3) default(present)
    do k = 1, nz
    do j = 1, shift_n
    do i = 1, sample_nx
        u_s(i,j,k) = u(sample_iwrap_cuda(i),ny-shift_n+j,k)
        v_s(i,j,k) = v(sample_iwrap_cuda(i),ny-shift_n+j,k)
        w_s(i,j,k) = w(sample_iwrap_cuda(i),ny-shift_n+j,k)
#ifdef PPSCALARS
        theta_s(i,j,k) = theta(sample_iwrap_cuda(i),ny-shift_n+j,k)
#endif
    end do
    end do
    end do

    !$acc parallel loop collapse(3) default(present) private(i_w)
    do k = 1, nz
    do j = 1, ny
    do i = 1, apply_nx
        i_w = apply_iwrap_cuda(i)
        u(i_w,j,k) = apply_alpha_cuda(i) * u(i_w,j,k)                         &
            + apply_beta_cuda(i) * u_s(i,j,k)
        v(i_w,j,k) = apply_alpha_cuda(i) * v(i_w,j,k)                         &
            + apply_beta_cuda(i) * v_s(i,j,k)
        w(i_w,j,k) = apply_alpha_cuda(i) * w(i_w,j,k)                         &
            + apply_beta_cuda(i) * w_s(i,j,k)
#ifdef PPSCALARS
        theta(i_w,j,k) = apply_alpha_cuda(i) * theta(i_w,j,k)                 &
            + apply_beta_cuda(i) * theta_s(i,j,k)
#endif
    end do
    end do
    end do
    return
end if
#endif

! Sample and shift velocity
u_s(:,shift_n+1:ny,:) = u(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
u_s(:,1:shift_n,:) = u(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
v_s(:,shift_n+1:ny,:) = v(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
v_s(:,1:shift_n,:) = v(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
w_s(:,shift_n+1:ny,:) = w(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
w_s(:,1:shift_n,:) = w(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
#ifdef PPSCALARS
theta_s(:,shift_n+1:ny,:) = theta(sample_fringe%iwrap(:),1:ny-shift_n,1:nz)
theta_s(:,1:shift_n,:) = theta(sample_fringe%iwrap(:),ny-shift_n+1:ny,1:nz)
#endif

! Apply inflow conditions
do i = 1, apply_nx
    i_w = apply_fringe%iwrap(i)
    u(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * u(i_w,1:ny,1:nz)                  &
        + apply_fringe%beta(i) * u_s(i,1:ny,1:nz)
    v(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * v(i_w,1:ny,1:nz)                  &
        + apply_fringe%beta(i) * v_s(i,1:ny,1:nz)
    w(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * w(i_w,1:ny,1:nz)                  &
        + apply_fringe%beta(i) * w_s(i,1:ny,1:nz)
#ifdef PPSCALARS
    theta(i_w,1:ny,1:nz) = apply_fringe%alpha(i) * theta(i_w,1:ny,1:nz)          &
        + apply_fringe%beta(i) * theta_s(i,1:ny,1:nz)
#endif
end do

end subroutine inflow_shifted

end module shifted_inflow
