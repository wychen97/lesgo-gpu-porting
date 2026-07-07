!!
!!  Copyright (C) 2010-2017  Johns Hopkins University
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
module derivatives
!*******************************************************************************
!
! This module contains all of the major subroutines used for computing
! derivatives.
!
! Navigation map:
!   - CUDA FFT helper paths: filt_da_cuda, xy_derivs_cuda, stress_*_cuda
!   - CUFFT plan management: ensure_*_cuda_plan routines
!   - synchronization/error helpers: derivatives_cuda_sync and require_* checks
!   - public CPU wrappers: ddx, ddy, ddxy, filt_da, filt_da_vel
!   - vertical derivatives: ddz_vel, ddz_uv, ddz_w
!
! The separate `derivatives_gpu.f90` module owns the optimized OpenACC/CUDA
! production path used by `USE_LES_GPU`.  Keep wrapper semantics identical.
implicit none

save
private

public ddx, ddy, ddxy, filt_da, filt_da_vel, ddz_vel, ddz_uv, ddz_w
public stress_uv_xy_derivs, stress_w_xy_derivs


contains


!*******************************************************************************
subroutine stress_uv_xy_derivs(txx, tyy, txy, dtxdx, dtydy2, dtxdx2, dtydy,   &
    lbz)
!*******************************************************************************
!
! Compute the horizontal stress derivatives used by divstress_uv.
!
use types, only : rprec
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: txx, tyy, txy
real(rprec), dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy2
real(rprec), dimension(:,:,lbz:), intent(inout) :: dtxdx2, dtydy

call ddx(txx, dtxdx, lbz)
call ddy(tyy, dtydy2, lbz)
call ddxy(txy, dtxdx2, dtydy, lbz)

end subroutine stress_uv_xy_derivs

!*******************************************************************************
subroutine stress_w_xy_derivs(tx, ty, dtxdx, dtydy, lbz)
!*******************************************************************************
!
! Compute the horizontal stress derivatives used by divstress_w.
!
use types, only : rprec
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: tx, ty
real(rprec), dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy

call ddx(tx, dtxdx, lbz)
call ddy(ty, dtydy, lbz)

end subroutine stress_w_xy_derivs

!*******************************************************************************
subroutine ddx(f,dfdx,lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz
use fft, only : back, forw, kx
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx
real(rprec) :: const
integer :: jz


const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lbz, nz
    !  Use dfdx to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = const*f(:,:,jz)
    call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz),dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(ld-1:ld,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdx to perform complex multiplication
    ! Optimized version for real(eye*kx)=0
    ! only passing imaginary part of eye*kx
    dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
enddo

end subroutine ddx

!*******************************************************************************
subroutine ddy(f,dfdy, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! y using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz
use fft, only : back, forw, ky
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
real(rprec) :: const
integer :: jz


const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lbz, nz
    !  Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdy(:,:,jz) = const * f(:,:,jz)
    call dfftw_execute_dft_r2c(forw, dfdy(:,:,jz), dfdy(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdy(ld-1:ld,:,jz) = 0._rprec
    dfdy(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdy(:,:,jz) .MULI. ky

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine ddy

!*******************************************************************************
subroutine ddxy (f, dfdx, dfdy, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to
! x and y using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz
use fft, only : back, forw, kx, ky
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
real(rprec) :: const
integer :: jz


const = 1._rprec / ( nx * ny )

! Loop through horizontal slices
do jz = lbz, nz
    ! Use dfdy to hold f; since we are doing in place FFTs this is required
    dfdx(:,:,jz) = const*f(:,:,jz)
    call dfftw_execute_dft_r2c(forw, dfdx(:,:,jz), dfdx(:,:,jz))

    ! Zero padded region and Nyquist frequency
    dfdx(ld-1:ld,:,jz) = 0._rprec
    dfdx(:,ny/2+1,jz) = 0._rprec

    ! Derivatives: must to y's first here, because we're using dfdx as storage
    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdy(:,:,jz) = dfdx(:,:,jz) .MULI. ky
    dfdx(:,:,jz) = dfdx(:,:,jz) .MULI. kx

    ! Perform inverse transform to get pseudospectral derivative
    call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine ddxy

!*******************************************************************************
subroutine filt_da(f,dfdx,dfdy, lbz)
!*******************************************************************************
!
! This subroutine kills the oddball components in f and computes the partial
! derivative of f with respect to x and y using spectral decomposition.
!
use types, only : rprec
use param, only : ld, nx, ny, nz
use fft, only : back, forw, kx, ky
use emul_complex, only : OPERATOR(.MULI.)
implicit none


integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(inout) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
real(rprec) :: const
integer :: jz


const = 1._rprec/(nx*ny)

! loop through horizontal slices
do jz = lbz, nz
    ! Calculate FFT in place
    f(:,:,jz) = const*f(:,:,jz)
    call dfftw_execute_dft_r2c(forw, f(:,:,jz), f(:,:,jz))

    ! Kill oddballs in zero padded region and Nyquist frequency
    f(ld-1:ld,:,jz) = 0._rprec
    f(:,ny/2+1,jz) = 0._rprec

    ! Use complex emulation of dfdy to perform complex multiplication
    ! Optimized version for real(eye*ky)=0
    ! only passing imaginary part of eye*ky
    dfdx(:,:,jz) = f(:,:,jz) .MULI. kx
    dfdy(:,:,jz) = f(:,:,jz) .MULI. ky

    ! Perform inverse transform to get pseudospectral derivative
    ! The oddballs for derivatives should already be dead, since they are for f
    ! inverse transform
    call dfftw_execute_dft_c2r(back, f(:,:,jz), f(:,:,jz))
    call dfftw_execute_dft_c2r(back, dfdx(:,:,jz), dfdx(:,:,jz))
    call dfftw_execute_dft_c2r(back, dfdy(:,:,jz), dfdy(:,:,jz))
end do

end subroutine filt_da

!*******************************************************************************
subroutine filt_da_vel(u, v, w, dudx, dudy, dvdx, dvdy, dwdx, dwdy, lbz)
!*******************************************************************************
!
! Combined velocity derivative/filter wrapper.  CUDA builds use a single
! 3-component batched cuFFT path; CPU/fallback builds keep the original calls.
!
use types, only : rprec
use param, only : ld, nx, ny, nz
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(inout) :: u, v, w
real(rprec), dimension(:,:,lbz:), intent(inout) :: dudx, dudy
real(rprec), dimension(:,:,lbz:), intent(inout) :: dvdx, dvdy
real(rprec), dimension(:,:,lbz:), intent(inout) :: dwdx, dwdy


call filt_da(u, dudx, dudy, lbz)
call filt_da(v, dvdx, dvdy, lbz)
call filt_da(w, dwdx, dwdy, lbz)


end subroutine filt_da_vel

!*******************************************************************************
subroutine ddz_vel(u, v, w, dudz, dvdz, dwdz, lbz)
!*******************************************************************************
!
! Combined velocity z-derivative wrapper. CUDA builds use one set of kernels for
! all three velocity components to avoid repeated launch/sync overhead.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, dz, BOGUS
#ifdef PPSAFETYMODE
use param, only : nproc, coord
#endif
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: u, v, w
real(rprec), dimension(:,:,lbz:), intent(inout) :: dudz, dvdz, dwdz
integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec/dz


call ddz_uv(u, dudz, lbz)
call ddz_uv(v, dvdz, lbz)
call ddz_w(w, dwdz, lbz)

end subroutine ddz_vel

!*******************************************************************************
subroutine ddz_uv(f, dfdz, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing. f is on the uv grid and dfdz is on the w grid.
! The serial version provides dfdz(:,:,2:nz), and the value at jz=1 is not
! touched. The MPI version provides dfdz(:,:,1:nz), except at the bottom
! process it only supplies 2:nz
!
use types, only : rprec
use param, only : ld, nx, ny, nz, dz, BOGUS
#ifdef PPSAFETYMODE
use param, only : nproc, coord
#endif
implicit none

integer, intent(in) :: lbz
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec/dz


#if defined(PPMPI) && defined(PPSAFETYMODE)
dfdz(:,:,0) = BOGUS
#endif

! Calculate derivative.
! The ghost node information is available here
! if coord == 0, dudz(1) will be set in wallstress
do jz = lbz+1, nz
do jy = 1, ny
do jx = 1, nx
    dfdz(jx,jy,jz) = const*(f(jx,jy,jz)-f(jx,jy,jz-1))
end do
end do
end do

! Not necessarily accurate at top and bottom boundary
! Set to BOGUS just to be safe
#ifdef PPSAFETYMODE
if (coord == 0) then
    dfdz(:,:,1) = BOGUS
end if
if (coord == nproc-1) then
    dfdz(:,:,nz) = BOGUS
end if
#endif

end subroutine ddz_uv

!*******************************************************************************
subroutine ddz_w(f, dfdz, lbz)
!*******************************************************************************
!
! This subroutine computes the partial derivative of f with respect to z using
! 2nd order finite differencing. f is on the w grid and dfdz is on the uv grid.
! The serial version provides dfdz(:,:,1:nz-1), and the value at jz=1 is not
! touched. The MPI version provides dfdz(:,:,0:nz-1), except at the top and
! bottom processes, which each has has 0:nz, and 1:nz-1, respectively.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, dz, BOGUS
#ifdef PPSAFETYMODE
#ifdef PPMPI
use param, only : coord
#endif
#endif
implicit none

real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
integer, intent(in) :: lbz
real(rprec)::const
integer :: jx, jy, jz

const = 1._rprec/dz


do jz = lbz, nz-1
do jy = 1, ny
do jx = 1, nx
    dfdz(jx,jy,jz) = const*(f(jx,jy,jz+1)-f(jx,jy,jz))
end do
end do
end do

#ifdef PPSAFETYMODE
#ifdef PPMPI
! bottom process cannot calculate dfdz(jz=0)
if (coord == 0) then
    dfdz(:,:,lbz) = BOGUS
endif
#endif
! All processes cannot calculate dfdz(jz=nz)
dfdz(:,:,nz) = BOGUS
#endif

end subroutine ddz_w

end module derivatives
