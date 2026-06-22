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
subroutine divstress_uv (divtx, divty, txx, txy, txz, tyy, tyz)
!*******************************************************************************
!
! This subroutine provides divt for 1:nz. MPI provides 1:nz-1,
! except at top, where 1:nz is provided
!
use types, only : rprec
use param, only : ld, ny, nz, BOGUS, lbz
use derivatives, only : ddz_w, stress_uv_xy_derivs
#ifdef ENABLE_CUDA
use cudafor
use derivatives, only : stress_uv_div_cuda
#endif
implicit none

#ifdef ENABLE_CUDA
real(rprec), managed, dimension(ld,ny,lbz:nz), intent(out) :: divtx, divty
real(rprec), managed, dimension(ld, ny, lbz:nz), intent(in) :: txx, txy, txz
real(rprec), managed, dimension(ld, ny, lbz:nz), intent(in) :: tyy, tyz
real(rprec), managed, allocatable, save, dimension(:,:,:) :: dtxdx, dtydy
real(rprec), managed, allocatable, save, dimension(:,:,:) :: dtzdz
real(rprec), managed, allocatable, save, dimension(:,:,:) :: dtxdx2, dtydy2
real(rprec), managed, allocatable, save, dimension(:,:,:) :: dtzdz2
#else
real(rprec), dimension(ld,ny,lbz:nz), intent(out) :: divtx, divty
real(rprec), dimension(ld, ny, lbz:nz), intent(in) :: txx, txy, txz, tyy, tyz
real(rprec), dimension(ld,ny,lbz:nz) :: dtxdx, dtydy, dtzdz
real(rprec), dimension(ld,ny,lbz:nz) :: dtxdx2, dtydy2, dtzdz2
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
    call stress_uv_div_cuda(txx, tyy, txy, txz, tyz, divtx, divty, lbz)
    return
end if

if (.not. allocated(dtxdx) .or. scratch_ld /= ld .or. scratch_ny /= ny .or.  &
    scratch_lbz /= lbz .or. scratch_nz /= nz) then
    if (allocated(dtxdx)) then
        deallocate(dtxdx, dtydy, dtzdz, dtxdx2, dtydy2, dtzdz2)
    end if
    allocate(dtxdx(ld,ny,lbz:nz), dtydy(ld,ny,lbz:nz),                       &
        dtzdz(ld,ny,lbz:nz), dtxdx2(ld,ny,lbz:nz),                           &
        dtydy2(ld,ny,lbz:nz), dtzdz2(ld,ny,lbz:nz))
    scratch_ld = ld
    scratch_ny = ny
    scratch_lbz = lbz
    scratch_nz = nz
end if
#endif

! compute stress gradients
! MPI: txx, tyy, and txy 1:nz-1 => horizontal gradients 1:nz-1
call stress_uv_xy_derivs(txx, tyy, txy, dtxdx, dtydy2, dtxdx2, dtydy, lbz)

! MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(txz, dtzdz, lbz)

! MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(tyz, dtzdz2, lbz)

! Historical wall-node correction reference.  The active path uses ddz_w above;
! applying this extra first-level correction changed the wall stress balance.
!      dtzdz(:,:,1) = (tz(:,:,2)-tz(:,:,1))/(0.5*dz)

#ifdef ENABLE_CUDA
if (cuda_enabled) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, ld
        if (jx >= ld - 1) then
            divtx(jx,jy,jz) = 0._rprec
            divty(jx,jy,jz) = 0._rprec
        else
            divtx(jx,jy,jz) = dtxdx(jx,jy,jz) + dtydy(jx,jy,jz) +             &
                dtzdz(jx,jy,jz)
            divty(jx,jy,jz) = dtxdx2(jx,jy,jz) + dtydy2(jx,jy,jz) +           &
                dtzdz2(jx,jy,jz)
        end if
    end do
    end do
    end do

#ifdef PPSAFETYMODE
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
#ifdef PPMPI
        divtx(jx,jy,0) = BOGUS
        divty(jx,jy,0) = BOGUS
#endif
        divtx(jx,jy,nz) = BOGUS
        divty(jx,jy,nz) = BOGUS
    end do
    end do
#endif

    if (cuda_extra_sync) then
        istat = cudaDeviceSynchronize()
        if (istat /= 0) then
            print *, 'divstress_uv CUDA sync failure: ', istat
            stop
        end if
    end if
    istat = cudaGetLastError()
    if (istat /= 0) then
        print *, 'divstress_uv CUDA kernel failure: ', istat
        stop
    end if
else
#endif
! only 1:nz-1 are valid
divtx(:,:,1:nz-1) = dtxdx(:,:,1:nz-1) + dtydy(:,:,1:nz-1) + dtzdz(:,:,1:nz-1)

! Set ld-1, ld to 0 (or could do BOGUS)
divtx(ld-1:ld, :, 1:nz-1) = 0._rprec

#ifdef PPSAFETYMODE
#ifdef PPMPI
divtx(:,:,0) = BOGUS
#endif
divtx(:,:,nz) = BOGUS
#endif

! only 1:nz-1 are valid
divty(:,:,1:nz-1) = dtxdx2(:,:,1:nz-1) + dtydy2(:,:,1:nz-1) + dtzdz2(:,:,1:nz-1)

! Set ld-1, ld to 0 (or could do BOGUS)
divty(ld-1:ld,:,1:nz-1) = 0._rprec

#ifdef PPSAFETYMODE
#ifdef PPMPI
divty(:,:,0) = BOGUS
#endif
divty(:,:,nz) = BOGUS
#endif

#ifdef ENABLE_CUDA
end if
#endif

end subroutine divstress_uv
#endif
