!!
!!  Copyright (C) 2016-2020  Johns Hopkins University
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
module inflow
!*******************************************************************************
use types, only : rprec
use param, only : inflow_type
use fringe
#ifdef PPCPS
use concurrent_precursor, only : synchronize_cps, inflow_cps
#endif
#ifdef PPHIT
use hit_inflow
#endif
use shifted_inflow
implicit none

private
public :: inflow_init, apply_inflow
#if defined(PPLES_GPU)
public :: inflow_cuda_enabled
#endif

type (fringe_t) :: uniform_fringe
#if defined(PPLES_GPU)
integer, allocatable, dimension(:) :: uniform_iwrap_cuda
real(rprec), allocatable, dimension(:) :: uniform_alpha_cuda
real(rprec), allocatable, dimension(:) :: uniform_beta_cuda
!$acc declare create(uniform_iwrap_cuda, uniform_alpha_cuda, uniform_beta_cuda)
#endif

contains

#if defined(PPLES_GPU)
!*******************************************************************************
logical function inflow_cuda_enabled()
!*******************************************************************************
implicit none

inflow_cuda_enabled = .true.

end function inflow_cuda_enabled

#endif

!*******************************************************************************
subroutine inflow_init
!*******************************************************************************
use param, only : fringe_region_end, fringe_region_len

select case (inflow_type)
    ! uniform
    case (1)
        uniform_fringe = fringe_t(fringe_region_end, fringe_region_len)
#if defined(PPLES_GPU)
        if (allocated(uniform_iwrap_cuda)) deallocate(uniform_iwrap_cuda)
        if (allocated(uniform_alpha_cuda)) deallocate(uniform_alpha_cuda)
        if (allocated(uniform_beta_cuda)) deallocate(uniform_beta_cuda)
        allocate(uniform_iwrap_cuda(uniform_fringe%nx))
        allocate(uniform_alpha_cuda(uniform_fringe%nx))
        allocate(uniform_beta_cuda(uniform_fringe%nx))
        uniform_iwrap_cuda = uniform_fringe%iwrap
        uniform_alpha_cuda = uniform_fringe%alpha
        uniform_beta_cuda = uniform_fringe%beta
#ifdef PPLES_GPU
        !$acc update device(uniform_iwrap_cuda, uniform_alpha_cuda,             &
        !$acc& uniform_beta_cuda)
#endif
#endif
#ifdef PPHIT
    ! HIT
    case (2)
        call inflow_HIT()
#endif
    ! shifted
    case (3)
        call shifted_inflow_init()
#ifdef PPCPS
    ! CPS
    case (4)
        ! CPS sample buffers are allocated later by initialize_cps(), after
        ! the velocity/scalar initial conditions exist.  The actual CPS fringe
        ! application happens in apply_inflow().
        continue
#endif
end select

end subroutine inflow_init

!*******************************************************************************
subroutine apply_inflow
!*******************************************************************************

select case (inflow_type)
    ! uniform
    case (1)
        call inflow_uniform()
#ifdef PPHIT
    ! HIT
    case (2)
        call inflow_HIT()
#endif
    ! shifted
    case (3)
        call inflow_shifted()
#ifdef PPCPS
    ! CPS
    case (4)
        call synchronize_cps()
        call inflow_cps()
#endif
end select

end subroutine apply_inflow

!*******************************************************************************
subroutine inflow_uniform ()
!*******************************************************************************
!  Enforces prescribed inflow condition based on an uniform inflow
!  velocity.
use param, only : nx, ny, nz, inflow_velocity
use sim_param, only : u, v, w
integer :: i, i_w, j, k, uniform_nx
real(rprec) :: alpha_i, beta_i

uniform_nx = uniform_fringe%nx

#ifdef PPLES_GPU
if (inflow_cuda_enabled()) then
    !$acc parallel loop collapse(3) default(present) private(i_w, alpha_i, beta_i)
    do k = 1, nz
    do j = 1, ny
    do i = 1, uniform_nx
        i_w = uniform_iwrap_cuda(i)
        alpha_i = uniform_alpha_cuda(i)
        beta_i = uniform_beta_cuda(i)
        u(i_w,j,k) = alpha_i*u(i_w,j,k) + beta_i*inflow_velocity
        v(i_w,j,k) = alpha_i*v(i_w,j,k)
        w(i_w,j,k) = alpha_i*w(i_w,j,k)
    end do
    end do
    end do
else
#endif
!--skip istart since we know vel at istart, iend already
do i = 1, uniform_fringe%nx
    i_w = uniform_fringe%iwrap(i)
    u(i_w,1:ny,1:nz) = uniform_fringe%alpha(i) * u(i_w,1:ny,1:nz)              &
        + uniform_fringe%beta(i) * inflow_velocity
    v(i_w,1:ny,1:nz) = uniform_fringe%alpha(i) * v(i_w,1:ny,1:nz)
    w(i_w,1:ny,1:nz) = uniform_fringe%alpha(i) * w(i_w,1:ny,1:nz)
end do
#if defined(PPLES_GPU)
end if
#endif

end subroutine inflow_uniform

end module inflow
