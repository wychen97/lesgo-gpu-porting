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
#ifdef ENABLE_CUDA
use types, only : rprec
use cudafor
#endif
implicit none

save
private

public ddx, ddy, ddxy, filt_da, filt_da_vel, ddz_vel, ddz_uv, ddz_w
public stress_uv_xy_derivs, stress_w_xy_derivs
#ifdef ENABLE_CUDA
public stress_uv_div_cuda, stress_w_div_cuda
#endif

#ifdef ENABLE_CUDA
logical, save :: filt_da_cuda_initialized = .false.
integer, save :: filt_da_fw_plan = 0
integer, save :: filt_da_bk_plan = 0
integer, save :: filt_da_plan_lbz = -999999
integer, save :: filt_da_plan_batch = -1
logical, save :: filt_da_vel_cuda_initialized = .false.
integer, save :: filt_da_vel_fw_plan = 0
integer, save :: filt_da_vel_bk_plan = 0
integer, save :: filt_da_vel_plan_lbz = -999999
integer, save :: filt_da_vel_plan_batch = -1
logical, save :: stress_uv_xy_cuda_initialized = .false.
integer, save :: stress_uv_xy_fw_plan = 0
integer, save :: stress_uv_xy_bk_plan = 0
integer, save :: stress_uv_xy_plan_lbz = -999999
integer, save :: stress_uv_xy_plan_batch = -1
logical, save :: stress_w_xy_cuda_initialized = .false.
integer, save :: stress_w_xy_fw_plan = 0
integer, save :: stress_w_xy_bk_plan = 0
integer, save :: stress_w_xy_plan_lbz = -999999
integer, save :: stress_w_xy_plan_batch = -1
real(rprec), device, allocatable, save :: filt_da_f_d(:,:,:)
real(rprec), device, allocatable, save :: filt_da_dfdx_d(:,:,:)
real(rprec), device, allocatable, save :: filt_da_dfdy_d(:,:,:)
real(rprec), device, allocatable, save :: filt_da_vel_f_d(:,:,:,:)
real(rprec), device, allocatable, save :: filt_da_vel_dfdx_d(:,:,:,:)
real(rprec), device, allocatable, save :: filt_da_vel_dfdy_d(:,:,:,:)
real(rprec), device, allocatable, save :: stress_uv_xy_f_d(:,:,:,:)
real(rprec), device, allocatable, save :: stress_uv_xy_deriv_d(:,:,:,:)
real(rprec), device, allocatable, save :: stress_w_xy_f_d(:,:,:,:)
real(rprec), device, allocatable, save :: stress_w_xy_deriv_d(:,:,:,:)
real(rprec), device, allocatable, save :: filt_da_kx_d(:,:)
real(rprec), device, allocatable, save :: filt_da_ky_d(:,:)
complex(rprec), device, allocatable, save :: filt_da_fh_d(:,:,:)
complex(rprec), device, allocatable, save :: filt_da_dfdxh_d(:,:,:)
complex(rprec), device, allocatable, save :: filt_da_dfdyh_d(:,:,:)
complex(rprec), device, allocatable, save :: filt_da_vel_fh_d(:,:,:,:)
complex(rprec), device, allocatable, save :: filt_da_vel_dfdxh_d(:,:,:,:)
complex(rprec), device, allocatable, save :: filt_da_vel_dfdyh_d(:,:,:,:)
complex(rprec), device, allocatable, save :: stress_uv_xy_fh_d(:,:,:,:)
complex(rprec), device, allocatable, save :: stress_uv_xy_dh_d(:,:,:,:)
complex(rprec), device, allocatable, save :: stress_w_xy_fh_d(:,:,:,:)
complex(rprec), device, allocatable, save :: stress_w_xy_dh_d(:,:,:,:)
#endif

contains

#ifdef ENABLE_CUDA
!*******************************************************************************
subroutine filt_da_cuda(f, dfdx, dfdy, lbz)
!*******************************************************************************
!
! Batched cuFFT implementation of filt_da.  The PlanMany layout is reversed for
! Fortran column-major storage: real slices are f(ld,ny), so cuFFT sees ny by nx.
!
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use cufft
implicit none

integer, intent(in) :: lbz
    real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: f
    real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy

integer :: istat
integer :: i, j, k, kk
real(rprec) :: const, fr, fi

const = 1._rprec / (nx * ny)

call ensure_filt_da_cuda_plan(lbz)

    !$cuf kernel do(3) <<<*,*>>>
    do kk = 1, filt_da_plan_batch
    do j = 1, ny
    do i = 1, ld
        k = lbz + kk - 1
        filt_da_f_d(i,j,kk) = const * f(i,j,k)
    end do
    end do
    end do

istat = cufftExecD2Z(filt_da_fw_plan, filt_da_f_d, filt_da_fh_d)
call require_filt_da_cufft_success('forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, filt_da_plan_batch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        filt_da_fh_d(i,j,kk) = cmplx(0._rprec, 0._rprec, kind=rprec)
        filt_da_dfdxh_d(i,j,kk) = cmplx(0._rprec, 0._rprec, kind=rprec)
        filt_da_dfdyh_d(i,j,kk) = cmplx(0._rprec, 0._rprec, kind=rprec)
    else
        fr = real(filt_da_fh_d(i,j,kk), kind=rprec)
        fi = aimag(filt_da_fh_d(i,j,kk))
        filt_da_dfdxh_d(i,j,kk) = cmplx(-fi * filt_da_kx_d(i,j),             &
            fr * filt_da_kx_d(i,j), kind=rprec)
        filt_da_dfdyh_d(i,j,kk) = cmplx(-fi * filt_da_ky_d(i,j),             &
            fr * filt_da_ky_d(i,j), kind=rprec)
    end if
end do
end do
end do

istat = cudaGetLastError()
call require_filt_da_cuda_success('derivative kernel', istat)

istat = cufftExecZ2D(filt_da_bk_plan, filt_da_fh_d, filt_da_f_d)
call require_filt_da_cufft_success('inverse filtered exec', istat)
istat = cufftExecZ2D(filt_da_bk_plan, filt_da_dfdxh_d, filt_da_dfdx_d)
call require_filt_da_cufft_success('inverse dfdx exec', istat)
istat = cufftExecZ2D(filt_da_bk_plan, filt_da_dfdyh_d, filt_da_dfdy_d)
call require_filt_da_cufft_success('inverse dfdy exec', istat)
    !$cuf kernel do(3) <<<*,*>>>
    do kk = 1, filt_da_plan_batch
    do j = 1, ny
    do i = 1, ld
        k = lbz + kk - 1
        f(i,j,k) = filt_da_f_d(i,j,kk)
        dfdx(i,j,k) = filt_da_dfdx_d(i,j,kk)
        dfdy(i,j,k) = filt_da_dfdy_d(i,j,kk)
    end do
    end do
    end do

end subroutine filt_da_cuda

!*******************************************************************************
subroutine xy_derivs_cuda(f, lbz, dfdx, dfdy)
!*******************************************************************************
!
! Batched cuFFT implementation for ddx, ddy, and ddxy.  Inputs and outputs may
! be host or managed arrays; only the internal work arrays are used in kernels.
!
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use cufft
implicit none

integer, intent(in) :: lbz
    real(rprec), managed, dimension(:,:,lbz:), intent(in) :: f
    real(rprec), managed, dimension(:,:,lbz:), intent(inout), optional :: dfdx, dfdy

integer :: istat
integer :: i, j, k, kk
real(rprec) :: const, fr, fi

if (.not. present(dfdx) .and. .not. present(dfdy)) return

const = 1._rprec / (nx * ny)

call ensure_filt_da_cuda_plan(lbz)

    !$cuf kernel do(3) <<<*,*>>>
    do kk = 1, filt_da_plan_batch
    do j = 1, ny
    do i = 1, ld
        k = lbz + kk - 1
        filt_da_f_d(i,j,kk) = const * f(i,j,k)
    end do
    end do
    end do

istat = cufftExecD2Z(filt_da_fw_plan, filt_da_f_d, filt_da_fh_d)
call require_filt_da_cufft_success('xy forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, filt_da_plan_batch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        filt_da_dfdxh_d(i,j,kk) = cmplx(0._rprec, 0._rprec, kind=rprec)
        filt_da_dfdyh_d(i,j,kk) = cmplx(0._rprec, 0._rprec, kind=rprec)
    else
        fr = real(filt_da_fh_d(i,j,kk), kind=rprec)
        fi = aimag(filt_da_fh_d(i,j,kk))
        filt_da_dfdxh_d(i,j,kk) = cmplx(-fi * filt_da_kx_d(i,j),             &
            fr * filt_da_kx_d(i,j), kind=rprec)
        filt_da_dfdyh_d(i,j,kk) = cmplx(-fi * filt_da_ky_d(i,j),             &
            fr * filt_da_ky_d(i,j), kind=rprec)
    end if
end do
end do
end do

if (derivatives_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    call require_filt_da_cuda_success('xy derivative kernel sync', istat)
end if
istat = cudaGetLastError()
call require_filt_da_cuda_success('xy derivative kernel', istat)

if (present(dfdx)) then
    istat = cufftExecZ2D(filt_da_bk_plan, filt_da_dfdxh_d, filt_da_dfdx_d)
    call require_filt_da_cufft_success('xy inverse dfdx exec', istat)
end if
if (present(dfdy)) then
    istat = cufftExecZ2D(filt_da_bk_plan, filt_da_dfdyh_d, filt_da_dfdy_d)
    call require_filt_da_cufft_success('xy inverse dfdy exec', istat)
end if
if (derivatives_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    call require_filt_da_cuda_success('xy inverse sync', istat)
end if

    if (present(dfdx) .and. present(dfdy)) then
        !$cuf kernel do(3) <<<*,*>>>
        do kk = 1, filt_da_plan_batch
        do j = 1, ny
        do i = 1, ld
            k = lbz + kk - 1
            dfdx(i,j,k) = filt_da_dfdx_d(i,j,kk)
            dfdy(i,j,k) = filt_da_dfdy_d(i,j,kk)
        end do
        end do
        end do
    else if (present(dfdx)) then
        !$cuf kernel do(3) <<<*,*>>>
        do kk = 1, filt_da_plan_batch
        do j = 1, ny
        do i = 1, ld
            k = lbz + kk - 1
            dfdx(i,j,k) = filt_da_dfdx_d(i,j,kk)
        end do
        end do
        end do
    else if (present(dfdy)) then
        !$cuf kernel do(3) <<<*,*>>>
        do kk = 1, filt_da_plan_batch
        do j = 1, ny
        do i = 1, ld
            k = lbz + kk - 1
            dfdy(i,j,k) = filt_da_dfdy_d(i,j,kk)
        end do
        end do
        end do
    end if

call derivatives_cuda_sync('xy_derivs_cuda')

end subroutine xy_derivs_cuda

!*******************************************************************************
subroutine stress_uv_xy_derivs_cuda(txx, tyy, txy, dtxdx, dtydy2, dtxdx2,      &
    dtydy, lbz)
!*******************************************************************************
!
! Batched horizontal derivatives for divstress_uv.  This keeps the original
! spectral operations but groups txx, tyy, and txy into one forward cuFFT and
! the four required derivative fields into one inverse cuFFT.
!
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use cufft
implicit none

integer, intent(in) :: lbz
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: txx, tyy, txy
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy2
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dtxdx2, dtydy

integer :: istat
integer :: i, j, k, kk
real(rprec) :: const, fr1, fi1, fr2, fi2, fr3, fi3

const = 1._rprec / (nx * ny)

call ensure_stress_uv_xy_cuda_plan(lbz)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_uv_xy_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    stress_uv_xy_f_d(i,j,kk,1) = const * txx(i,j,k)
    stress_uv_xy_f_d(i,j,kk,2) = const * tyy(i,j,k)
    stress_uv_xy_f_d(i,j,kk,3) = const * txy(i,j,k)
end do
end do
end do

istat = cufftExecD2Z(stress_uv_xy_fw_plan, stress_uv_xy_f_d(1,1,1,1),         &
    stress_uv_xy_fh_d(1,1,1,1))
call require_filt_da_cufft_success('stress uv xy forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_uv_xy_plan_batch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        stress_uv_xy_dh_d(i,j,kk,1) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,2) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,3) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,4) = cmplx(0._rprec, 0._rprec, kind=rprec)
    else
        fr1 = real(stress_uv_xy_fh_d(i,j,kk,1), kind=rprec)
        fi1 = aimag(stress_uv_xy_fh_d(i,j,kk,1))
        fr2 = real(stress_uv_xy_fh_d(i,j,kk,2), kind=rprec)
        fi2 = aimag(stress_uv_xy_fh_d(i,j,kk,2))
        fr3 = real(stress_uv_xy_fh_d(i,j,kk,3), kind=rprec)
        fi3 = aimag(stress_uv_xy_fh_d(i,j,kk,3))
        stress_uv_xy_dh_d(i,j,kk,1) = cmplx(-fi1 * filt_da_kx_d(i,j),          &
            fr1 * filt_da_kx_d(i,j), kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,2) = cmplx(-fi2 * filt_da_ky_d(i,j),          &
            fr2 * filt_da_ky_d(i,j), kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,3) = cmplx(-fi3 * filt_da_kx_d(i,j),          &
            fr3 * filt_da_kx_d(i,j), kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,4) = cmplx(-fi3 * filt_da_ky_d(i,j),          &
            fr3 * filt_da_ky_d(i,j), kind=rprec)
    end if
end do
end do
end do

if (derivatives_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    call require_filt_da_cuda_success('stress uv xy derivative sync', istat)
end if
istat = cudaGetLastError()
call require_filt_da_cuda_success('stress uv xy derivative kernel', istat)

istat = cufftExecZ2D(stress_uv_xy_bk_plan, stress_uv_xy_dh_d(1,1,1,1),        &
    stress_uv_xy_deriv_d(1,1,1,1))
call require_filt_da_cufft_success('stress uv xy inverse exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_uv_xy_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    dtxdx(i,j,k) = stress_uv_xy_deriv_d(i,j,kk,1)
    dtydy2(i,j,k) = stress_uv_xy_deriv_d(i,j,kk,2)
    dtxdx2(i,j,k) = stress_uv_xy_deriv_d(i,j,kk,3)
    dtydy(i,j,k) = stress_uv_xy_deriv_d(i,j,kk,4)
end do
end do
end do

call derivatives_cuda_sync('stress_uv_xy_derivs_cuda')

end subroutine stress_uv_xy_derivs_cuda

!*******************************************************************************
subroutine stress_w_xy_derivs_cuda(tx, ty, dtxdx, dtydy, lbz)
!*******************************************************************************
!
! Batched horizontal derivatives for divstress_w.
!
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use cufft
implicit none

integer, intent(in) :: lbz
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: tx, ty
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy

integer :: istat
integer :: i, j, k, kk
real(rprec) :: const, fr1, fi1, fr2, fi2

const = 1._rprec / (nx * ny)

call ensure_stress_w_xy_cuda_plan(lbz)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_w_xy_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    stress_w_xy_f_d(i,j,kk,1) = const * tx(i,j,k)
    stress_w_xy_f_d(i,j,kk,2) = const * ty(i,j,k)
end do
end do
end do

istat = cufftExecD2Z(stress_w_xy_fw_plan, stress_w_xy_f_d(1,1,1,1),           &
    stress_w_xy_fh_d(1,1,1,1))
call require_filt_da_cufft_success('stress w xy forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_w_xy_plan_batch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        stress_w_xy_dh_d(i,j,kk,1) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_w_xy_dh_d(i,j,kk,2) = cmplx(0._rprec, 0._rprec, kind=rprec)
    else
        fr1 = real(stress_w_xy_fh_d(i,j,kk,1), kind=rprec)
        fi1 = aimag(stress_w_xy_fh_d(i,j,kk,1))
        fr2 = real(stress_w_xy_fh_d(i,j,kk,2), kind=rprec)
        fi2 = aimag(stress_w_xy_fh_d(i,j,kk,2))
        stress_w_xy_dh_d(i,j,kk,1) = cmplx(-fi1 * filt_da_kx_d(i,j),           &
            fr1 * filt_da_kx_d(i,j), kind=rprec)
        stress_w_xy_dh_d(i,j,kk,2) = cmplx(-fi2 * filt_da_ky_d(i,j),           &
            fr2 * filt_da_ky_d(i,j), kind=rprec)
    end if
end do
end do
end do

if (derivatives_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    call require_filt_da_cuda_success('stress w xy derivative sync', istat)
end if
istat = cudaGetLastError()
call require_filt_da_cuda_success('stress w xy derivative kernel', istat)

istat = cufftExecZ2D(stress_w_xy_bk_plan, stress_w_xy_dh_d(1,1,1,1),          &
    stress_w_xy_deriv_d(1,1,1,1))
call require_filt_da_cufft_success('stress w xy inverse exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_w_xy_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    dtxdx(i,j,k) = stress_w_xy_deriv_d(i,j,kk,1)
    dtydy(i,j,k) = stress_w_xy_deriv_d(i,j,kk,2)
end do
end do
end do

call derivatives_cuda_sync('stress_w_xy_derivs_cuda')

end subroutine stress_w_xy_derivs_cuda

!*******************************************************************************
subroutine stress_uv_div_cuda(txx, tyy, txy, txz, tyz, divtx, divty, lbz)
!*******************************************************************************
!
! Complete GPU divstress_uv path.  The horizontal derivatives remain in the
! module device work array from the batched stress derivative transform, and the
! vertical stress derivatives are formed directly in the final combine kernel.
!
use types, only : rprec
use param, only : ld, lh, ny, nz, nx, dz, BOGUS
use cufft
implicit none

integer, intent(in) :: lbz
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: txx, tyy, txy
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: txz, tyz
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: divtx, divty

integer :: istat
integer :: i, j, k, kk
real(rprec) :: const_xy, const_z, dz_txz, dz_tyz
real(rprec) :: fr1, fi1, fr2, fi2, fr3, fi3

const_xy = 1._rprec / (nx * ny)
const_z = 1._rprec / dz

call ensure_stress_uv_xy_cuda_plan(lbz)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_uv_xy_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    stress_uv_xy_f_d(i,j,kk,1) = const_xy * txx(i,j,k)
    stress_uv_xy_f_d(i,j,kk,2) = const_xy * tyy(i,j,k)
    stress_uv_xy_f_d(i,j,kk,3) = const_xy * txy(i,j,k)
end do
end do
end do

istat = cufftExecD2Z(stress_uv_xy_fw_plan, stress_uv_xy_f_d(1,1,1,1),         &
    stress_uv_xy_fh_d(1,1,1,1))
call require_filt_da_cufft_success('stress uv div forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_uv_xy_plan_batch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        stress_uv_xy_dh_d(i,j,kk,1) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,2) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,3) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,4) = cmplx(0._rprec, 0._rprec, kind=rprec)
    else
        fr1 = real(stress_uv_xy_fh_d(i,j,kk,1), kind=rprec)
        fi1 = aimag(stress_uv_xy_fh_d(i,j,kk,1))
        fr2 = real(stress_uv_xy_fh_d(i,j,kk,2), kind=rprec)
        fi2 = aimag(stress_uv_xy_fh_d(i,j,kk,2))
        fr3 = real(stress_uv_xy_fh_d(i,j,kk,3), kind=rprec)
        fi3 = aimag(stress_uv_xy_fh_d(i,j,kk,3))
        stress_uv_xy_dh_d(i,j,kk,1) = cmplx(-fi1 * filt_da_kx_d(i,j),          &
            fr1 * filt_da_kx_d(i,j), kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,2) = cmplx(-fi2 * filt_da_ky_d(i,j),          &
            fr2 * filt_da_ky_d(i,j), kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,3) = cmplx(-fi3 * filt_da_kx_d(i,j),          &
            fr3 * filt_da_kx_d(i,j), kind=rprec)
        stress_uv_xy_dh_d(i,j,kk,4) = cmplx(-fi3 * filt_da_ky_d(i,j),          &
            fr3 * filt_da_ky_d(i,j), kind=rprec)
    end if
end do
end do
end do

if (derivatives_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    call require_filt_da_cuda_success('stress uv div derivative sync', istat)
end if
istat = cudaGetLastError()
call require_filt_da_cuda_success('stress uv div derivative kernel', istat)

istat = cufftExecZ2D(stress_uv_xy_bk_plan, stress_uv_xy_dh_d(1,1,1,1),        &
    stress_uv_xy_deriv_d(1,1,1,1))
call require_filt_da_cufft_success('stress uv div inverse exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do k = 1, nz - 1
do j = 1, ny
do i = 1, ld
    kk = k - lbz + 1
    if (i >= ld - 1) then
        divtx(i,j,k) = 0._rprec
        divty(i,j,k) = 0._rprec
    else
        dz_txz = const_z * (txz(i,j,k+1) - txz(i,j,k))
        dz_tyz = const_z * (tyz(i,j,k+1) - tyz(i,j,k))
        divtx(i,j,k) = stress_uv_xy_deriv_d(i,j,kk,1) +                      &
            stress_uv_xy_deriv_d(i,j,kk,4) + dz_txz
        divty(i,j,k) = stress_uv_xy_deriv_d(i,j,kk,3) +                      &
            stress_uv_xy_deriv_d(i,j,kk,2) + dz_tyz
    end if
end do
end do
end do

#ifdef PPSAFETYMODE
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
#ifdef PPMPI
    divtx(i,j,0) = BOGUS
    divty(i,j,0) = BOGUS
#endif
    divtx(i,j,nz) = BOGUS
    divty(i,j,nz) = BOGUS
end do
end do
#endif

istat = cudaDeviceSynchronize()
call require_filt_da_cuda_success('stress uv div final sync', istat)
istat = cudaGetLastError()
call require_filt_da_cuda_success('stress uv div final kernel', istat)

end subroutine stress_uv_div_cuda

!*******************************************************************************
subroutine stress_w_div_cuda(tx, ty, tz, divt, lbz)
!*******************************************************************************
!
! Complete GPU divstress_w path with the same scratch-elimination strategy used
! for stress_uv_div_cuda.
!
use types, only : rprec
use param, only : ld, lh, nx, ny, nz, dz, coord, nproc, BOGUS
use cufft
implicit none

integer, intent(in) :: lbz
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: tx, ty, tz
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: divt

integer :: istat
integer :: i, j, k, kk
real(rprec) :: const_xy, const_z, dz_tz
real(rprec) :: fr1, fi1, fr2, fi2

const_xy = 1._rprec / (nx * ny)
const_z = 1._rprec / dz

call ensure_stress_w_xy_cuda_plan(lbz)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_w_xy_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    stress_w_xy_f_d(i,j,kk,1) = const_xy * tx(i,j,k)
    stress_w_xy_f_d(i,j,kk,2) = const_xy * ty(i,j,k)
end do
end do
end do

istat = cufftExecD2Z(stress_w_xy_fw_plan, stress_w_xy_f_d(1,1,1,1),           &
    stress_w_xy_fh_d(1,1,1,1))
call require_filt_da_cufft_success('stress w div forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, stress_w_xy_plan_batch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        stress_w_xy_dh_d(i,j,kk,1) = cmplx(0._rprec, 0._rprec, kind=rprec)
        stress_w_xy_dh_d(i,j,kk,2) = cmplx(0._rprec, 0._rprec, kind=rprec)
    else
        fr1 = real(stress_w_xy_fh_d(i,j,kk,1), kind=rprec)
        fi1 = aimag(stress_w_xy_fh_d(i,j,kk,1))
        fr2 = real(stress_w_xy_fh_d(i,j,kk,2), kind=rprec)
        fi2 = aimag(stress_w_xy_fh_d(i,j,kk,2))
        stress_w_xy_dh_d(i,j,kk,1) = cmplx(-fi1 * filt_da_kx_d(i,j),           &
            fr1 * filt_da_kx_d(i,j), kind=rprec)
        stress_w_xy_dh_d(i,j,kk,2) = cmplx(-fi2 * filt_da_ky_d(i,j),           &
            fr2 * filt_da_ky_d(i,j), kind=rprec)
    end if
end do
end do
end do

if (derivatives_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    call require_filt_da_cuda_success('stress w div derivative sync', istat)
end if
istat = cudaGetLastError()
call require_filt_da_cuda_success('stress w div derivative kernel', istat)

istat = cufftExecZ2D(stress_w_xy_bk_plan, stress_w_xy_dh_d(1,1,1,1),          &
    stress_w_xy_deriv_d(1,1,1,1))
call require_filt_da_cufft_success('stress w div inverse exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do k = 1, nz
do j = 1, ny
do i = 1, ld
    kk = k - lbz + 1
    if (i > nx) then
        divt(i,j,k) = 0._rprec
    else if ((k == 1 .and. coord == 0) .or.                                  &
        (k == nz .and. coord == nproc - 1)) then
        divt(i,j,k) = stress_w_xy_deriv_d(i,j,kk,1) +                         &
            stress_w_xy_deriv_d(i,j,kk,2)
    else
        dz_tz = const_z * (tz(i,j,k) - tz(i,j,k-1))
        divt(i,j,k) = stress_w_xy_deriv_d(i,j,kk,1) +                         &
            stress_w_xy_deriv_d(i,j,kk,2) + dz_tz
    end if
end do
end do
end do

#ifdef PPSAFETYMODE
#ifdef PPMPI
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    divt(i,j,0) = BOGUS
end do
end do
#endif
#endif

istat = cudaDeviceSynchronize()
call require_filt_da_cuda_success('stress w div final sync', istat)
istat = cudaGetLastError()
call require_filt_da_cuda_success('stress w div final kernel', istat)

end subroutine stress_w_div_cuda

!*******************************************************************************
subroutine filt_da_vel_cuda(u, v, w, dudx, dudy, dvdx, dvdy, dwdx, dwdy, lbz)
!*******************************************************************************
!
! Batched velocity version of filt_da.  It uses one 3-component cuFFT batch
! instead of calling filt_da separately for u, v, and w.
!
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use cufft
implicit none

integer, intent(in) :: lbz
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: u, v, w
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dudx, dudy
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dvdx, dvdy
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dwdx, dwdy

integer :: istat
integer :: i, j, k, kk, c
real(rprec) :: const, fr, fi

const = 1._rprec / (nx * ny)

call ensure_filt_da_vel_cuda_plan(lbz)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, filt_da_vel_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    filt_da_vel_f_d(i,j,kk,1) = const * u(i,j,k)
    filt_da_vel_f_d(i,j,kk,2) = const * v(i,j,k)
    filt_da_vel_f_d(i,j,kk,3) = const * w(i,j,k)
end do
end do
end do

istat = cufftExecD2Z(filt_da_vel_fw_plan, filt_da_vel_f_d(1,1,1,1),            &
    filt_da_vel_fh_d(1,1,1,1))
call require_filt_da_cufft_success('velocity forward exec', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, filt_da_vel_plan_batch
do j = 1, ny
do i = 1, lh
    do c = 1, 3
        if (i == lh .or. j == ny/2 + 1) then
            filt_da_vel_fh_d(i,j,kk,c) = cmplx(0._rprec, 0._rprec, kind=rprec)
            filt_da_vel_dfdxh_d(i,j,kk,c) = cmplx(0._rprec, 0._rprec,          &
                kind=rprec)
            filt_da_vel_dfdyh_d(i,j,kk,c) = cmplx(0._rprec, 0._rprec,          &
                kind=rprec)
        else
            fr = real(filt_da_vel_fh_d(i,j,kk,c), kind=rprec)
            fi = aimag(filt_da_vel_fh_d(i,j,kk,c))
            filt_da_vel_dfdxh_d(i,j,kk,c) = cmplx(-fi * filt_da_kx_d(i,j),     &
                fr * filt_da_kx_d(i,j), kind=rprec)
            filt_da_vel_dfdyh_d(i,j,kk,c) = cmplx(-fi * filt_da_ky_d(i,j),     &
                fr * filt_da_ky_d(i,j), kind=rprec)
        end if
    end do
end do
end do
end do

istat = cudaGetLastError()
call require_filt_da_cuda_success('velocity derivative kernel', istat)

istat = cufftExecZ2D(filt_da_vel_bk_plan, filt_da_vel_fh_d(1,1,1,1),           &
    filt_da_vel_f_d(1,1,1,1))
call require_filt_da_cufft_success('velocity inverse filtered exec', istat)
istat = cufftExecZ2D(filt_da_vel_bk_plan, filt_da_vel_dfdxh_d(1,1,1,1),        &
    filt_da_vel_dfdx_d(1,1,1,1))
call require_filt_da_cufft_success('velocity inverse dfdx exec', istat)
istat = cufftExecZ2D(filt_da_vel_bk_plan, filt_da_vel_dfdyh_d(1,1,1,1),        &
    filt_da_vel_dfdy_d(1,1,1,1))
call require_filt_da_cufft_success('velocity inverse dfdy exec', istat)
istat = cudaDeviceSynchronize()
call require_filt_da_cuda_success('velocity inverse sync', istat)

!$cuf kernel do(3) <<<*,*>>>
do kk = 1, filt_da_vel_plan_batch
do j = 1, ny
do i = 1, ld
    k = lbz + kk - 1
    u(i,j,k) = filt_da_vel_f_d(i,j,kk,1)
    dudx(i,j,k) = filt_da_vel_dfdx_d(i,j,kk,1)
    dudy(i,j,k) = filt_da_vel_dfdy_d(i,j,kk,1)
    v(i,j,k) = filt_da_vel_f_d(i,j,kk,2)
    dvdx(i,j,k) = filt_da_vel_dfdx_d(i,j,kk,2)
    dvdy(i,j,k) = filt_da_vel_dfdy_d(i,j,kk,2)
    w(i,j,k) = filt_da_vel_f_d(i,j,kk,3)
    dwdx(i,j,k) = filt_da_vel_dfdx_d(i,j,kk,3)
    dwdy(i,j,k) = filt_da_vel_dfdy_d(i,j,kk,3)
end do
end do
end do

end subroutine filt_da_vel_cuda

!*******************************************************************************
subroutine ensure_filt_da_cuda_plan(lbz)
!*******************************************************************************
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use fft, only : kx, ky
use cufft
implicit none

integer, intent(in) :: lbz
integer :: batch, istat
integer :: n_s(2), inem_s(2), onem_s(2)

batch = nz - lbz + 1

if (filt_da_cuda_initialized .and. filt_da_plan_lbz == lbz .and.              &
    filt_da_plan_batch == batch) then
    return
end if

if (filt_da_cuda_initialized) then
    istat = cufftDestroy(filt_da_fw_plan)
    call require_filt_da_cufft_success('destroy forward plan', istat)
    istat = cufftDestroy(filt_da_bk_plan)
    call require_filt_da_cufft_success('destroy inverse plan', istat)
    filt_da_cuda_initialized = .false.
end if

if (allocated(filt_da_f_d)) deallocate(filt_da_f_d)
if (allocated(filt_da_dfdx_d)) deallocate(filt_da_dfdx_d)
if (allocated(filt_da_dfdy_d)) deallocate(filt_da_dfdy_d)
if (allocated(filt_da_fh_d)) deallocate(filt_da_fh_d)
if (allocated(filt_da_dfdxh_d)) deallocate(filt_da_dfdxh_d)
if (allocated(filt_da_dfdyh_d)) deallocate(filt_da_dfdyh_d)

allocate(filt_da_f_d(ld, ny, batch))
allocate(filt_da_dfdx_d(ld, ny, batch))
allocate(filt_da_dfdy_d(ld, ny, batch))
allocate(filt_da_fh_d(lh, ny, batch))
allocate(filt_da_dfdxh_d(lh, ny, batch))
allocate(filt_da_dfdyh_d(lh, ny, batch))

if (.not. allocated(filt_da_kx_d)) allocate(filt_da_kx_d(lh, ny))
if (.not. allocated(filt_da_ky_d)) allocate(filt_da_ky_d(lh, ny))
filt_da_kx_d = kx
filt_da_ky_d = ky

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(filt_da_fw_plan, 2, n_s, inem_s, 1, ld*ny,              &
    onem_s, 1, lh*ny, CUFFT_D2Z, batch)
call require_filt_da_cufft_success('forward plan', istat)

istat = cufftPlanMany(filt_da_bk_plan, 2, n_s, onem_s, 1, lh*ny,              &
    inem_s, 1, ld*ny, CUFFT_Z2D, batch)
call require_filt_da_cufft_success('inverse plan', istat)

filt_da_plan_lbz = lbz
filt_da_plan_batch = batch
filt_da_cuda_initialized = .true.

end subroutine ensure_filt_da_cuda_plan

!*******************************************************************************
subroutine ensure_filt_da_vel_cuda_plan(lbz)
!*******************************************************************************
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use fft, only : kx, ky
use cufft
implicit none

integer, intent(in) :: lbz
integer :: batch, istat
integer :: n_s(2), inem_s(2), onem_s(2)

batch = nz - lbz + 1

if (filt_da_vel_cuda_initialized .and. filt_da_vel_plan_lbz == lbz .and.       &
    filt_da_vel_plan_batch == batch) then
    return
end if

if (filt_da_vel_cuda_initialized) then
    istat = cufftDestroy(filt_da_vel_fw_plan)
    call require_filt_da_cufft_success('destroy velocity forward plan', istat)
    istat = cufftDestroy(filt_da_vel_bk_plan)
    call require_filt_da_cufft_success('destroy velocity inverse plan', istat)
    filt_da_vel_cuda_initialized = .false.
end if

if (allocated(filt_da_vel_f_d)) deallocate(filt_da_vel_f_d)
if (allocated(filt_da_vel_dfdx_d)) deallocate(filt_da_vel_dfdx_d)
if (allocated(filt_da_vel_dfdy_d)) deallocate(filt_da_vel_dfdy_d)
if (allocated(filt_da_vel_fh_d)) deallocate(filt_da_vel_fh_d)
if (allocated(filt_da_vel_dfdxh_d)) deallocate(filt_da_vel_dfdxh_d)
if (allocated(filt_da_vel_dfdyh_d)) deallocate(filt_da_vel_dfdyh_d)

allocate(filt_da_vel_f_d(ld, ny, batch, 3))
allocate(filt_da_vel_dfdx_d(ld, ny, batch, 3))
allocate(filt_da_vel_dfdy_d(ld, ny, batch, 3))
allocate(filt_da_vel_fh_d(lh, ny, batch, 3))
allocate(filt_da_vel_dfdxh_d(lh, ny, batch, 3))
allocate(filt_da_vel_dfdyh_d(lh, ny, batch, 3))

if (.not. allocated(filt_da_kx_d)) allocate(filt_da_kx_d(lh, ny))
if (.not. allocated(filt_da_ky_d)) allocate(filt_da_ky_d(lh, ny))
filt_da_kx_d = kx
filt_da_ky_d = ky

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(filt_da_vel_fw_plan, 2, n_s, inem_s, 1, ld*ny,           &
    onem_s, 1, lh*ny, CUFFT_D2Z, 3*batch)
call require_filt_da_cufft_success('velocity forward plan', istat)

istat = cufftPlanMany(filt_da_vel_bk_plan, 2, n_s, onem_s, 1, lh*ny,           &
    inem_s, 1, ld*ny, CUFFT_Z2D, 3*batch)
call require_filt_da_cufft_success('velocity inverse plan', istat)

filt_da_vel_plan_lbz = lbz
filt_da_vel_plan_batch = batch
filt_da_vel_cuda_initialized = .true.

end subroutine ensure_filt_da_vel_cuda_plan

!*******************************************************************************
subroutine ensure_stress_uv_xy_cuda_plan(lbz)
!*******************************************************************************
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use fft, only : kx, ky
use cufft
implicit none

integer, intent(in) :: lbz
integer :: batch, istat
integer :: n_s(2), inem_s(2), onem_s(2)

batch = nz - lbz + 1

if (stress_uv_xy_cuda_initialized .and. stress_uv_xy_plan_lbz == lbz .and.    &
    stress_uv_xy_plan_batch == batch) then
    return
end if

if (stress_uv_xy_cuda_initialized) then
    istat = cufftDestroy(stress_uv_xy_fw_plan)
    call require_filt_da_cufft_success('destroy stress uv xy forward plan',    &
        istat)
    istat = cufftDestroy(stress_uv_xy_bk_plan)
    call require_filt_da_cufft_success('destroy stress uv xy inverse plan',    &
        istat)
    stress_uv_xy_cuda_initialized = .false.
end if

if (allocated(stress_uv_xy_f_d)) deallocate(stress_uv_xy_f_d)
if (allocated(stress_uv_xy_deriv_d)) deallocate(stress_uv_xy_deriv_d)
if (allocated(stress_uv_xy_fh_d)) deallocate(stress_uv_xy_fh_d)
if (allocated(stress_uv_xy_dh_d)) deallocate(stress_uv_xy_dh_d)

allocate(stress_uv_xy_f_d(ld, ny, batch, 3))
allocate(stress_uv_xy_deriv_d(ld, ny, batch, 4))
allocate(stress_uv_xy_fh_d(lh, ny, batch, 3))
allocate(stress_uv_xy_dh_d(lh, ny, batch, 4))

if (.not. allocated(filt_da_kx_d)) allocate(filt_da_kx_d(lh, ny))
if (.not. allocated(filt_da_ky_d)) allocate(filt_da_ky_d(lh, ny))
filt_da_kx_d = kx
filt_da_ky_d = ky

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(stress_uv_xy_fw_plan, 2, n_s, inem_s, 1, ld*ny,        &
    onem_s, 1, lh*ny, CUFFT_D2Z, 3*batch)
call require_filt_da_cufft_success('stress uv xy forward plan', istat)

istat = cufftPlanMany(stress_uv_xy_bk_plan, 2, n_s, onem_s, 1, lh*ny,        &
    inem_s, 1, ld*ny, CUFFT_Z2D, 4*batch)
call require_filt_da_cufft_success('stress uv xy inverse plan', istat)

stress_uv_xy_plan_lbz = lbz
stress_uv_xy_plan_batch = batch
stress_uv_xy_cuda_initialized = .true.

end subroutine ensure_stress_uv_xy_cuda_plan

!*******************************************************************************
subroutine ensure_stress_w_xy_cuda_plan(lbz)
!*******************************************************************************
use types, only : rprec
use param, only : ld, lh, nx, ny, nz
use fft, only : kx, ky
use cufft
implicit none

integer, intent(in) :: lbz
integer :: batch, istat
integer :: n_s(2), inem_s(2), onem_s(2)

batch = nz - lbz + 1

if (stress_w_xy_cuda_initialized .and. stress_w_xy_plan_lbz == lbz .and.      &
    stress_w_xy_plan_batch == batch) then
    return
end if

if (stress_w_xy_cuda_initialized) then
    istat = cufftDestroy(stress_w_xy_fw_plan)
    call require_filt_da_cufft_success('destroy stress w xy forward plan',     &
        istat)
    istat = cufftDestroy(stress_w_xy_bk_plan)
    call require_filt_da_cufft_success('destroy stress w xy inverse plan',     &
        istat)
    stress_w_xy_cuda_initialized = .false.
end if

if (allocated(stress_w_xy_f_d)) deallocate(stress_w_xy_f_d)
if (allocated(stress_w_xy_deriv_d)) deallocate(stress_w_xy_deriv_d)
if (allocated(stress_w_xy_fh_d)) deallocate(stress_w_xy_fh_d)
if (allocated(stress_w_xy_dh_d)) deallocate(stress_w_xy_dh_d)

allocate(stress_w_xy_f_d(ld, ny, batch, 2))
allocate(stress_w_xy_deriv_d(ld, ny, batch, 2))
allocate(stress_w_xy_fh_d(lh, ny, batch, 2))
allocate(stress_w_xy_dh_d(lh, ny, batch, 2))

if (.not. allocated(filt_da_kx_d)) allocate(filt_da_kx_d(lh, ny))
if (.not. allocated(filt_da_ky_d)) allocate(filt_da_ky_d(lh, ny))
filt_da_kx_d = kx
filt_da_ky_d = ky

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(stress_w_xy_fw_plan, 2, n_s, inem_s, 1, ld*ny,         &
    onem_s, 1, lh*ny, CUFFT_D2Z, 2*batch)
call require_filt_da_cufft_success('stress w xy forward plan', istat)

istat = cufftPlanMany(stress_w_xy_bk_plan, 2, n_s, onem_s, 1, lh*ny,         &
    inem_s, 1, ld*ny, CUFFT_Z2D, 2*batch)
call require_filt_da_cufft_success('stress w xy inverse plan', istat)

stress_w_xy_plan_lbz = lbz
stress_w_xy_plan_batch = batch
stress_w_xy_cuda_initialized = .true.

end subroutine ensure_stress_w_xy_cuda_plan

!*******************************************************************************
logical function filt_da_cuda_enabled()
!*******************************************************************************
implicit none

filt_da_cuda_enabled = .true.

end function filt_da_cuda_enabled

!*******************************************************************************
logical function filt_da_vel_cuda_enabled()
!*******************************************************************************
implicit none

filt_da_vel_cuda_enabled = .false.

end function filt_da_vel_cuda_enabled

!*******************************************************************************
logical function spectral_derivs_cuda_enabled()
!*******************************************************************************
implicit none

spectral_derivs_cuda_enabled = .true.

end function spectral_derivs_cuda_enabled

!*******************************************************************************
logical function stress_uv_xy_cuda_enabled()
!*******************************************************************************
implicit none

stress_uv_xy_cuda_enabled = .true.

end function stress_uv_xy_cuda_enabled

!*******************************************************************************
logical function stress_w_xy_cuda_enabled()
!*******************************************************************************
implicit none

stress_w_xy_cuda_enabled = .true.

end function stress_w_xy_cuda_enabled

!*******************************************************************************
logical function vertical_derivs_cuda_enabled()
!*******************************************************************************
implicit none

vertical_derivs_cuda_enabled = .true.

end function vertical_derivs_cuda_enabled

!*******************************************************************************
logical function derivatives_extra_sync_enabled()
!*******************************************************************************
implicit none

derivatives_extra_sync_enabled = .false.

end function derivatives_extra_sync_enabled

!*******************************************************************************
subroutine derivatives_cuda_sync(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'derivatives CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'derivatives CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine derivatives_cuda_sync

!*******************************************************************************
subroutine require_filt_da_cufft_success(where, istat)
!*******************************************************************************
use cufft
implicit none

character(len=*), intent(in) :: where
integer, intent(in) :: istat

if (istat /= CUFFT_SUCCESS) then
    print *, 'filt_da_cuda cuFFT failure at ', trim(where), ': ', istat
    stop
end if

end subroutine require_filt_da_cufft_success

!*******************************************************************************
subroutine require_filt_da_cuda_success(where, istat)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer, intent(in) :: istat

if (istat /= 0) then
    print *, 'filt_da_cuda CUDA failure at ', trim(where), ': ', istat
    stop
end if

end subroutine require_filt_da_cuda_success
#endif

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
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: txx, tyy, txy
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy2
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dtxdx2, dtydy

if (stress_uv_xy_cuda_enabled()) then
    call stress_uv_xy_derivs_cuda(txx, tyy, txy, dtxdx, dtydy2, dtxdx2,        &
        dtydy, lbz)
    return
end if
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: txx, tyy, txy
real(rprec), dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy2
real(rprec), dimension(:,:,lbz:), intent(inout) :: dtxdx2, dtydy
#endif

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
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: tx, ty
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy

if (stress_w_xy_cuda_enabled()) then
    call stress_w_xy_derivs_cuda(tx, ty, dtxdx, dtydy, lbz)
    return
end if
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: tx, ty
real(rprec), dimension(:,:,lbz:), intent(inout) :: dtxdx, dtydy
#endif

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
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: f
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdx
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx
#endif
real(rprec) :: const
integer :: jz

#ifdef ENABLE_CUDA
if (spectral_derivs_cuda_enabled()) then
    call xy_derivs_cuda(f, lbz, dfdx=dfdx)
    return
end if
#endif

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
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: f
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdy
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdy
#endif
real(rprec) :: const
integer :: jz

#ifdef ENABLE_CUDA
if (spectral_derivs_cuda_enabled()) then
    call xy_derivs_cuda(f, lbz, dfdy=dfdy)
    return
end if
#endif

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
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none

integer, intent(in) :: lbz
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: f
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
#endif
real(rprec) :: const
integer :: jz

#ifdef ENABLE_CUDA
if (spectral_derivs_cuda_enabled()) then
    call xy_derivs_cuda(f, lbz, dfdx=dfdx, dfdy=dfdy)
    return
end if
#endif

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
use fft
use emul_complex, only : OPERATOR(.MULI.)
implicit none


integer, intent(in) :: lbz
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: f
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
#else
real(rprec), dimension(:,:,lbz:), intent(inout) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdx, dfdy
#endif
real(rprec) :: const
integer :: jz

#ifdef ENABLE_CUDA
if (filt_da_cuda_enabled()) then
    call filt_da_cuda(f, dfdx, dfdy, lbz)
    return
end if
#endif

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
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: u, v, w
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dudx, dudy
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dvdx, dvdy
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dwdx, dwdy
#else
real(rprec), dimension(:,:,lbz:), intent(inout) :: u, v, w
real(rprec), dimension(:,:,lbz:), intent(inout) :: dudx, dudy
real(rprec), dimension(:,:,lbz:), intent(inout) :: dvdx, dvdy
real(rprec), dimension(:,:,lbz:), intent(inout) :: dwdx, dwdy
#endif

#ifdef ENABLE_CUDA
if (filt_da_vel_cuda_enabled()) then
    call filt_da_vel_cuda(u, v, w, dudx, dudy, dvdx, dvdy, dwdx, dwdy, lbz)
    call derivatives_cuda_sync('filt_da_vel')
    return
end if
#endif

call filt_da(u, dudx, dudy, lbz)
call filt_da(v, dvdx, dvdy, lbz)
call filt_da(w, dwdx, dwdy, lbz)

#ifdef ENABLE_CUDA
call derivatives_cuda_sync('filt_da_vel')
#endif

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
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: u, v, w
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dudz, dvdz, dwdz
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: u, v, w
real(rprec), dimension(:,:,lbz:), intent(inout) :: dudz, dvdz, dwdz
#endif
integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec/dz

#ifdef ENABLE_CUDA
if (vertical_derivs_cuda_enabled()) then
#if defined(PPMPI) && defined(PPSAFETYMODE)
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dudz(jx,jy,0) = BOGUS
        dvdz(jx,jy,0) = BOGUS
    end do
    end do
#endif

    !$cuf kernel do(3) <<<*,*>>>
    do jz = lbz + 1, nz
    do jy = 1, ny
    do jx = 1, nx
        dudz(jx,jy,jz) = const * (u(jx,jy,jz) - u(jx,jy,jz-1))
        dvdz(jx,jy,jz) = const * (v(jx,jy,jz) - v(jx,jy,jz-1))
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do jz = lbz, nz - 1
    do jy = 1, ny
    do jx = 1, nx
        dwdz(jx,jy,jz) = const * (w(jx,jy,jz+1) - w(jx,jy,jz))
    end do
    end do
    end do

#ifdef PPSAFETYMODE
    if (coord == 0) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            dudz(jx,jy,1) = BOGUS
            dvdz(jx,jy,1) = BOGUS
#ifdef PPMPI
            dwdz(jx,jy,lbz) = BOGUS
#endif
        end do
        end do
    end if
    if (coord == nproc-1) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            dudz(jx,jy,nz) = BOGUS
            dvdz(jx,jy,nz) = BOGUS
        end do
        end do
    end if

    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dwdz(jx,jy,nz) = BOGUS
    end do
    end do
#endif

    call derivatives_cuda_sync('ddz_vel')
    return
end if
#endif

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
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: f
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdz
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
#endif
integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec/dz

#ifdef ENABLE_CUDA
if (vertical_derivs_cuda_enabled()) then
#if defined(PPMPI) && defined(PPSAFETYMODE)
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dfdz(jx,jy,0) = BOGUS
    end do
    end do
#endif

    !$cuf kernel do(3) <<<*,*>>>
    do jz = lbz + 1, nz
    do jy = 1, ny
    do jx = 1, nx
        dfdz(jx,jy,jz) = const * (f(jx,jy,jz) - f(jx,jy,jz-1))
    end do
    end do
    end do

#ifdef PPSAFETYMODE
    if (coord == 0) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            dfdz(jx,jy,1) = BOGUS
        end do
        end do
    end if
    if (coord == nproc-1) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            dfdz(jx,jy,nz) = BOGUS
        end do
        end do
    end if
#endif

    call derivatives_cuda_sync('ddz_uv')
    return
end if
#endif

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

#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,lbz:), intent(in) :: f
real(rprec), managed, dimension(:,:,lbz:), intent(inout) :: dfdz
#else
real(rprec), dimension(:,:,lbz:), intent(in) :: f
real(rprec), dimension(:,:,lbz:), intent(inout) :: dfdz
#endif
integer, intent(in) :: lbz
real(rprec)::const
integer :: jx, jy, jz

const = 1._rprec/dz

#ifdef ENABLE_CUDA
if (vertical_derivs_cuda_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = lbz, nz - 1
    do jy = 1, ny
    do jx = 1, nx
        dfdz(jx,jy,jz) = const * (f(jx,jy,jz+1) - f(jx,jy,jz))
    end do
    end do
    end do

#ifdef PPSAFETYMODE
#ifdef PPMPI
    if (coord == 0) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            dfdz(jx,jy,lbz) = BOGUS
        end do
        end do
    endif
#endif
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, ld
        dfdz(jx,jy,nz) = BOGUS
    end do
    end do
#endif

    call derivatives_cuda_sync('ddz_w')
    return
end if
#endif

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
