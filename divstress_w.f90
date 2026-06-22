#if !defined(PPSGS_GPU) || (defined(PPLVLSET) && defined(PPLES_GPU) && !defined(ENABLE_CUDA))
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
use param, only : ld, nx, ny, nz, coord, BOGUS, lbz, nproc, coord
use derivatives, only : ddz_uv, stress_w_xy_derivs
#ifdef ENABLE_CUDA
use cudafor
use derivatives, only : stress_w_div_cuda
#endif
implicit none

#ifdef ENABLE_CUDA
real(rprec), managed, dimension(ld,ny,lbz:nz), intent(out) :: divt
real(rprec), managed, dimension(ld,ny,lbz:nz), intent(in) :: tx, ty, tz
real(rprec), managed, allocatable, save, dimension(:,:,:) :: dtxdx, dtydy
real(rprec), managed, allocatable, save, dimension(:,:,:) :: dtzdz
#else
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: divt
real(rprec), dimension(ld,ny,lbz:nz), intent(in) :: tx, ty, tz
real(rprec), dimension(ld,ny,lbz:nz) :: dtxdx, dtydy, dtzdz
#endif
integer :: jx, jy, jz
#ifdef ENABLE_CUDA
character(len=16) :: cuda_setting
integer :: cuda_stat, istat
logical, save :: cuda_initialized = .false.
logical, save :: cuda_enabled = .true.
logical, save :: cuda_extra_sync = .false.
integer, save :: scratch_ld = -1, scratch_ny = -1, scratch_lbz = -999999
integer, save :: scratch_nz = -1

if (.not. cuda_initialized) then
    cuda_enabled = .true.
    cuda_extra_sync = .false.
    cuda_initialized = .true.
end if

if (cuda_enabled) then
    call stress_w_div_cuda(tx, ty, tz, divt, lbz)
    return
end if

if (.not. allocated(dtxdx) .or. scratch_ld /= ld .or. scratch_ny /= ny .or.  &
    scratch_lbz /= lbz .or. scratch_nz /= nz) then
    if (allocated(dtxdx)) then
        deallocate(dtxdx, dtydy, dtzdz)
    end if
    allocate(dtxdx(ld,ny,lbz:nz), dtydy(ld,ny,lbz:nz),                       &
        dtzdz(ld,ny,lbz:nz))
    scratch_ld = ld
    scratch_ny = ny
    scratch_lbz = lbz
    scratch_nz = nz
end if
#endif

! compute stress gradients
! tx, ty 1:nz => horizontal gradients 1:nz
call stress_w_xy_derivs(tx, ty, dtxdx, dtydy, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
#ifdef ENABLE_CUDA
if (cuda_enabled) then
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dtxdx(jx,jy,0) = BOGUS
    end do
    end do
else
#endif
dtxdx(:, :, 0) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef PPMPI
#ifdef ENABLE_CUDA
if (cuda_enabled) then
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dtydy(jx,jy,0) = BOGUS
    end do
    end do
else
#endif
dtydy(:, :, 0) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif
#endif

! tz 0:nz-1 (special case) => dtzdz 1:nz-1 (default), 2:nz-1 (bottom),
! and 1:nz (top)
call ddz_uv(tz, dtzdz, lbz)
#ifdef PPSAFETYMODE
#ifdef PPMPI
#ifdef ENABLE_CUDA
if (cuda_enabled) then
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dtzdz(jx,jy,0) = BOGUS
    end do
    end do
else
#endif
dtzdz(:, :, 0) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef PPMPI
#ifdef ENABLE_CUDA
if (cuda_enabled) then
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        divt(jx,jy,0) = BOGUS
    end do
    end do
else
#endif
divt(:, :, 0) = BOGUS
#ifdef ENABLE_CUDA
end if
#endif
#endif
#endif

#ifdef ENABLE_CUDA
if (cuda_enabled) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, nx
        if ((jz == 1 .and. coord == 0) .or.                                  &
            (jz == nz .and. coord == nproc-1)) then
            divt(jx,jy,jz) = dtxdx(jx,jy,jz) + dtydy(jx,jy,jz)
        else
            divt(jx,jy,jz) = dtxdx(jx,jy,jz) + dtydy(jx,jy,jz) +              &
                dtzdz(jx,jy,jz)
        end if
    end do
    end do
    end do

    !$cuf kernel do(2) <<<*,*>>>
    do jz = 1, nz - 1
    do jy = 1, ny
        divt(ld-1,jy,jz) = 0._rprec
        divt(ld,jy,jz) = 0._rprec
    end do
    end do

    if (cuda_extra_sync) then
        istat = cudaDeviceSynchronize()
        if (istat /= 0) then
            print *, 'divstress_w CUDA sync failure: ', istat
            stop
        end if
    end if
    istat = cudaGetLastError()
    if (istat /= 0) then
        print *, 'divstress_w CUDA kernel failure: ', istat
        stop
    end if
else
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

#ifdef ENABLE_CUDA
end if
#endif

end subroutine divstress_w
#endif
