!!
!!  Copyright (C) 2019  Johns Hopkins University
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
module sponge
!*******************************************************************************
! This module contains all of the subroutines associated with scalar transport
use types, only : rprec
#ifdef ENABLE_CUDA
use cudafor
#endif
implicit none

private
public :: sponge_init, sponge_force

! Sponge layer properties
logical, public :: use_sponge = .false.
real(rprec), public :: sponge_frequency = 3.9_rprec
real(rprec), public :: sponge_height = 0.75_rprec
#ifdef ENABLE_CUDA
real (rprec), managed, dimension (:), allocatable :: sp
#else
real (rprec), dimension (:), allocatable :: sp
#endif
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
!$acc declare create(sp)
#endif

contains

#ifdef ENABLE_CUDA
!******************************************************************************
logical function sponge_cuda_enabled()
!*******************************************************************************
implicit none

sponge_cuda_enabled = .true.

end function sponge_cuda_enabled

!******************************************************************************
subroutine sponge_cuda_sync(where)
!******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'sponge CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'sponge CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine sponge_cuda_sync
#endif

!******************************************************************************
subroutine sponge_init()
!******************************************************************************
use param, only : nz, lbz, pi, L_z
use types, only : rprec
use grid_m, only : grid

integer :: k

allocate (sp(lbz:nz)); sp = 0._rprec

do k = lbz, nz
    if (grid%z(k) > sponge_height) then
        sp(k) = 0.5_rprec*sponge_frequency*(1._rprec                           &
            - cos(pi*(grid%z(k) -sponge_height)/(L_z - sponge_height)))
    end if
end do

#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
!$acc update device(sp)
#endif

end subroutine sponge_init

!*******************************************************************************
subroutine sponge_force
!*******************************************************************************
! This subroutine calculates the sponge force term
use param, only : nx, ny, nz
use sim_param, only :  RHSx, RHSy, RHSz, u, v, w

integer :: k
#if defined(ENABLE_CUDA) || (defined(PPLES_GPU) && !defined(ENABLE_CUDA))
integer :: jx, jy
real(rprec) :: usum, vsum, wsum, umean, vmean, wmean, scale_uv, scale_w
#endif

if (use_sponge) then
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
    ! Direct OpenACC sponge path: wait for queued RHS/velocity work, then
    ! compute plane means and apply sponge forcing on the device.
    !$acc wait(1)
    do k = 1, nz-1
        usum = 0._rprec
        vsum = 0._rprec
        wsum = 0._rprec

        !$acc parallel loop collapse(2) default(present)                       &
        !$acc& reduction(+:usum, vsum, wsum)
        do jy = 1, ny
        do jx = 1, nx
            usum = usum + u(jx,jy,k)
            vsum = vsum + v(jx,jy,k)
            wsum = wsum + w(jx,jy,k)
        end do
        end do

        umean = usum / real(nx*ny, rprec)
        vmean = vsum / real(nx*ny, rprec)
        wmean = wsum / real(nx*ny, rprec)
        scale_uv = 0.5_rprec*(sp(k) + sp(k+1))
        scale_w = sp(k)

        !$acc parallel loop collapse(2) default(present)
        do jy = 1, ny
        do jx = 1, nx
            RHSx(jx,jy,k) = RHSx(jx,jy,k) - scale_uv                          &
                * (u(jx,jy,k) - umean)
            RHSy(jx,jy,k) = RHSy(jx,jy,k) - scale_uv                          &
                * (v(jx,jy,k) - vmean)
            RHSz(jx,jy,k) = RHSz(jx,jy,k) - scale_w                           &
                * (w(jx,jy,k) - wmean)
        end do
        end do
    end do
    return
#endif
#ifdef ENABLE_CUDA
    if (sponge_cuda_enabled()) then
        do k = 1, nz-1
            usum = 0._rprec
            vsum = 0._rprec
            wsum = 0._rprec

            !$cuf kernel do(2) <<<*,*>>> reduction(+:usum, vsum, wsum)
            do jy = 1, ny
            do jx = 1, nx
                usum = usum + u(jx,jy,k)
                vsum = vsum + v(jx,jy,k)
                wsum = wsum + w(jx,jy,k)
            end do
            end do
            call sponge_cuda_sync('plane means')

            umean = usum / real(nx*ny, rprec)
            vmean = vsum / real(nx*ny, rprec)
            wmean = wsum / real(nx*ny, rprec)
            scale_uv = 0.5_rprec*(sp(k) + sp(k+1))
            scale_w = sp(k)

            !$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, nx
                RHSx(jx,jy,k) = RHSx(jx,jy,k) - scale_uv                      &
                    * (u(jx,jy,k) - umean)
                RHSy(jx,jy,k) = RHSy(jx,jy,k) - scale_uv                      &
                    * (v(jx,jy,k) - vmean)
                RHSz(jx,jy,k) = RHSz(jx,jy,k) - scale_w                       &
                    * (w(jx,jy,k) - wmean)
            end do
            end do
            call sponge_cuda_sync('forcing apply')
        end do
        return
    end if
#endif
    do k = 1, nz-1
        RHSx(1:nx,1:ny,k) = RHSx(1:nx,1:ny,k) - 0.5_rprec*(sp(k) + sp(k+1))    &
            * (u(1:nx,1:ny,k) - sum(u(1:nx,1:ny,k))/(nx*ny))
        RHSy(1:nx,1:ny,k) = RHSy(1:nx,1:ny,k) - 0.5_rprec*(sp(k) + sp(k+1))    &
            * (v(1:nx,1:ny,k) - sum(v(1:nx,1:ny,k))/(nx*ny))
        RHSz(1:nx,1:ny,k) = RHSz(1:nx,1:ny,k) - sp(k)                          &
            * (w(1:nx,1:ny,k) - sum(w(1:nx,1:ny,k))/(nx*ny))
    end do
end if

end subroutine sponge_force


end module sponge
