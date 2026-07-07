#if !defined(PPSGS_GPU) || (defined(PPLVLSET) && defined(PPLES_GPU))
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
subroutine divstress_w(divt, tx, ty, tz)
!*******************************************************************************
!
! This subroutine provides divt for 1:nz. MPI provides 1:nz-1,
! except at top, where 1:nz is provided
!
use types, only : rprec
use param, only : ld, nx, ny, nz, coord, BOGUS, lbz, nproc
use derivatives, only : ddz_uv, stress_w_xy_derivs
implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: divt
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: tx, ty, tz
real(rprec), dimension(ld,ny,lbz:nz) :: dtxdx, dtydy, dtzdz
integer :: jx, jy, jz

! compute stress gradients
! tx, ty 1:nz => horizontal gradients 1:nz
call stress_w_xy_derivs(tx, ty, dtxdx, dtydy, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
dtxdx(:, :, 0) = BOGUS
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef PPMPI
dtydy(:, :, 0) = BOGUS
#endif
#endif

! tz 0:nz-1 (special case) => dtzdz 1:nz-1 (default), 2:nz-1 (bottom),
! and 1:nz (top)
call ddz_uv(tz, dtzdz, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
dtzdz(:, :, 0) = BOGUS
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef PPMPI
divt(:, :, 0) = BOGUS
#endif
#endif


if (coord == 0) then
    ! Physical bottom boundary: no resolved tzz gradient is available at the
    ! boundary node, so the wall-node divergence uses horizontal stress terms.
    do jy = 1, ny
    do jx = 1, nx
        divt(jx,jy,1) = dtxdx(jx,jy,1)+dtydy(jx,jy,1)
    end do
    end do
else
    do jy = 1, ny
    do jx = 1, nx
        divt(jx,jy,1) = dtxdx(jx,jy,1)+dtydy(jx,jy,1)+dtzdz(jx,jy,1)
    end do
    end do
end if

! Channel
if (coord == nproc-1) then
    ! Physical top boundary: mirror the bottom-wall closure and omit dtzdz at
    ! the boundary node.
    do jy = 1, ny
    do jx = 1, nx
        divt(jx,jy,nz) = dtxdx(jx,jy,nz)+dtydy(jx,jy,nz)
    end do
    end do
else
    do jy = 1, ny
    do jx = 1, nx
        divt(jx,jy,nz) = dtxdx(jx,jy,nz) + dtydy(jx,jy,nz) + dtzdz(jx,jy,nz)
    end do
    end do
end if

do jz = 2, nz-1
do jy = 1, ny
do jx = 1, nx
    divt(jx,jy,jz) = dtxdx(jx,jy,jz)+dtydy(jx,jy,jz)+dtzdz(jx,jy,jz)
end do
end do
end do

! Clear FFT-aligned padding columns outside the active x range.
divt(ld-1:ld, :, 1:nz-1) = 0._rprec


end subroutine divstress_w
#endif
