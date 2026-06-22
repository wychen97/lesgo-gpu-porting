#ifndef PPCONVEC_GPU
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
#ifdef ENABLE_CUDA
module convec_cuda_m
!*******************************************************************************
!
! FP64 CUDA implementation of convec.  This keeps the original rotation-form
! algorithm, but batches the 2D FFTs and uses CUF kernels for pointwise work.
!
! Navigation map:
!   - CUDA feature gates: convec_cuda_enabled and related helper switches
!   - CUDA implementation: convec_cuda_impl
!   - CUFFT plan setup: ensure_convec_cuda
!   - padding helpers: padd_3d_dp and unpadd_3d_dp
!   - diagnostics: check_convec_cuda, convec_cuda_sync, require_convec_cufft
!   - public fallback entry point: convec near the file end
!
! The optimized production `USE_LES_GPU` path uses convec_gpu.f90.  Keep this
! legacy/CUDA Fortran path correct as a fallback and validation reference.
use types, only : rprec
use param
use sim_param, only : u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy
use sim_param, only : RHSx, RHSy, RHSz
use cudafor
use cufft
implicit none

private
public :: convec_cuda_enabled, convec_cuda_impl

real(rprec), device, save, allocatable, dimension(:,:,:) :: cc_big_d
real(rprec), device, save, allocatable, dimension(:,:,:) :: u_big_d
real(rprec), device, save, allocatable, dimension(:,:,:) :: v_big_d
real(rprec), device, save, allocatable, dimension(:,:,:) :: w_big_d
real(rprec), device, save, allocatable, dimension(:,:,:) :: vort1_big_d
real(rprec), device, save, allocatable, dimension(:,:,:) :: vort2_big_d
real(rprec), device, save, allocatable, dimension(:,:,:) :: vort3_big_d

integer, save :: plan_fw_s_lbz = 0
integer, save :: plan_bk_b_lbz = 0
integer, save :: plan_fw_s_nz = 0
integer, save :: plan_bk_b_nz = 0
integer, save :: plan_fw_b_nz = 0
integer, save :: plan_bk_s_nz = 0
integer, save :: plan_fw_b_nz1 = 0
integer, save :: plan_bk_s_nz1 = 0
logical, save :: convec_cuda_initialized = .false.

contains

!*******************************************************************************
logical function convec_cuda_enabled()
!*******************************************************************************
implicit none

convec_cuda_enabled = .true.

end function convec_cuda_enabled

!*******************************************************************************
logical function convec_extra_sync_enabled()
!*******************************************************************************
implicit none

convec_extra_sync_enabled = .false.

end function convec_extra_sync_enabled

!*******************************************************************************
logical function convec_fused_pad_enabled()
!*******************************************************************************
implicit none

convec_fused_pad_enabled = .true.

end function convec_fused_pad_enabled

!*******************************************************************************
subroutine convec_cuda_impl()
!*******************************************************************************
implicit none

integer :: jx, jy, jz
integer :: jz_min, jz_max, jzLo, jzHi
integer :: istat
real(rprec) :: const

call ensure_convec_cuda()

if (sgs) then
    jzLo = 2
    jzHi = nz - 1
else
    jzLo = 1
    jzHi = nz - 1
end if

! Velocity fields: small physical -> small spectral -> padded big physical.
const = 1._rprec / (nx * ny)
!$cuf kernel do(3) <<<*,*>>>
do jz = lbz, nz
do jy = 1, ny
do jx = 1, ld
    RHSx(jx,jy,jz) = const * u(jx,jy,jz)
    RHSy(jx,jy,jz) = const * v(jx,jy,jz)
    RHSz(jx,jy,jz) = const * w(jx,jy,jz)
end do
end do
end do
call check_convec_cuda('velocity scale')

istat = cufftExecD2Z(plan_fw_s_lbz, RHSx(:,:,lbz), RHSx(:,:,lbz))
call require_convec_cufft('velocity u forward', istat)
istat = cufftExecD2Z(plan_fw_s_lbz, RHSy(:,:,lbz), RHSy(:,:,lbz))
call require_convec_cufft('velocity v forward', istat)
istat = cufftExecD2Z(plan_fw_s_lbz, RHSz(:,:,lbz), RHSz(:,:,lbz))
call require_convec_cufft('velocity w forward', istat)

call padd_3d_dp(u_big_d, RHSx, lbz, lbz, lbz, nz)
call padd_3d_dp(v_big_d, RHSy, lbz, lbz, lbz, nz)
call padd_3d_dp(w_big_d, RHSz, lbz, lbz, lbz, nz)

istat = cufftExecZ2D(plan_bk_b_lbz, u_big_d(:,:,lbz), u_big_d(:,:,lbz))
call require_convec_cufft('velocity u inverse big', istat)
istat = cufftExecZ2D(plan_bk_b_lbz, v_big_d(:,:,lbz), v_big_d(:,:,lbz))
call require_convec_cufft('velocity v inverse big', istat)
istat = cufftExecZ2D(plan_bk_b_lbz, w_big_d(:,:,lbz), w_big_d(:,:,lbz))
call require_convec_cufft('velocity w inverse big', istat)
call check_convec_cuda('velocity inverse big')

! Vorticity fields: small physical -> small spectral -> padded big physical.
!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    if (coord == 0 .and. jz == 1) then
        if (lbc_mom == 0) then
            RHSx(jx,jy,1) = 0._rprec
            RHSy(jx,jy,1) = 0._rprec
        else
            RHSx(jx,jy,1) = const * (0.5_rprec * (dwdy(jx,jy,1) +            &
                dwdy(jx,jy,2)) - dvdz(jx,jy,1))
            RHSy(jx,jy,1) = const * (dudz(jx,jy,1) - 0.5_rprec *             &
                (dwdx(jx,jy,1) + dwdx(jx,jy,2)))
        end if
    else if (coord == nproc-1 .and. jz == nz) then
        if (ubc_mom == 0) then
            RHSx(jx,jy,nz) = 0._rprec
            RHSy(jx,jy,nz) = 0._rprec
        else
            RHSx(jx,jy,nz) = const * (0.5_rprec * (dwdy(jx,jy,nz-1) +        &
                dwdy(jx,jy,nz)) - dvdz(jx,jy,nz-1))
            RHSy(jx,jy,nz) = const * (dudz(jx,jy,nz-1) - 0.5_rprec *         &
                (dwdx(jx,jy,nz-1) + dwdx(jx,jy,nz)))
        end if
    else
        RHSx(jx,jy,jz) = const * (dwdy(jx,jy,jz) - dvdz(jx,jy,jz))
        RHSy(jx,jy,jz) = const * (dudz(jx,jy,jz) - dwdx(jx,jy,jz))
    end if
    RHSz(jx,jy,jz) = const * (dvdx(jx,jy,jz) - dudy(jx,jy,jz))
end do
end do
end do
call check_convec_cuda('vorticity scale')

istat = cufftExecD2Z(plan_fw_s_nz, RHSx(:,:,1), RHSx(:,:,1))
call require_convec_cufft('vorticity x forward', istat)
istat = cufftExecD2Z(plan_fw_s_nz, RHSy(:,:,1), RHSy(:,:,1))
call require_convec_cufft('vorticity y forward', istat)
istat = cufftExecD2Z(plan_fw_s_nz, RHSz(:,:,1), RHSz(:,:,1))
call require_convec_cufft('vorticity z forward', istat)

call padd_3d_dp(vort1_big_d, RHSx, 1, lbz, 1, nz)
call padd_3d_dp(vort2_big_d, RHSy, 1, lbz, 1, nz)
call padd_3d_dp(vort3_big_d, RHSz, 1, lbz, 1, nz)

istat = cufftExecZ2D(plan_bk_b_nz, vort1_big_d(:,:,1), vort1_big_d(:,:,1))
call require_convec_cufft('vorticity x inverse big', istat)
istat = cufftExecZ2D(plan_bk_b_nz, vort2_big_d(:,:,1), vort2_big_d(:,:,1))
call require_convec_cufft('vorticity y inverse big', istat)
istat = cufftExecZ2D(plan_bk_b_nz, vort3_big_d(:,:,1), vort3_big_d(:,:,1))
call require_convec_cufft('vorticity z inverse big', istat)
call check_convec_cuda('vorticity inverse big')

! RHSx = v*w3 - w*w2 on the padded grid, then truncate.
const = 1._rprec / (nx2 * ny2)
!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz-1
do jy = 1, ny2
do jx = 1, ld_big
    if (coord == 0 .and. jz == 1) then
        cc_big_d(jx,jy,1) = const * (v_big_d(jx,jy,1) *                      &
            (-vort3_big_d(jx,jy,1)) + 0.5_rprec * w_big_d(jx,jy,2) *          &
            vort2_big_d(jx,jy,jzLo))
    else if (coord == nproc-1 .and. jz == nz-1) then
        cc_big_d(jx,jy,nz-1) = const * (v_big_d(jx,jy,nz-1) *                 &
            (-vort3_big_d(jx,jy,nz-1)) + 0.5_rprec * w_big_d(jx,jy,nz-1) *    &
            vort2_big_d(jx,jy,jzHi))
    else
        cc_big_d(jx,jy,jz) = const * (v_big_d(jx,jy,jz) *                     &
            (-vort3_big_d(jx,jy,jz)) + 0.5_rprec *                            &
            (w_big_d(jx,jy,jz+1) * vort2_big_d(jx,jy,jz+1) +                 &
            w_big_d(jx,jy,jz) * vort2_big_d(jx,jy,jz)))
    end if
end do
end do
end do
call check_convec_cuda('rhsx product')

istat = cufftExecD2Z(plan_fw_b_nz1, cc_big_d(:,:,1), cc_big_d(:,:,1))
call require_convec_cufft('rhsx forward big', istat)
call unpadd_3d_dp(RHSx, cc_big_d, lbz, 1, 1, nz-1)
istat = cufftExecZ2D(plan_bk_s_nz1, RHSx(:,:,1), RHSx(:,:,1))
call require_convec_cufft('rhsx inverse small', istat)

! RHSy = u*w3 - w*w1 on the padded grid, then truncate.
!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz-1
do jy = 1, ny2
do jx = 1, ld_big
    if (coord == 0 .and. jz == 1) then
        cc_big_d(jx,jy,1) = const * (u_big_d(jx,jy,1) *                       &
            vort3_big_d(jx,jy,1) + 0.5_rprec * w_big_d(jx,jy,2) *             &
            (-vort1_big_d(jx,jy,jzLo)))
    else if (coord == nproc-1 .and. jz == nz-1) then
        cc_big_d(jx,jy,nz-1) = const * (u_big_d(jx,jy,nz-1) *                 &
            vort3_big_d(jx,jy,nz-1) + 0.5_rprec * w_big_d(jx,jy,nz-1) *       &
            (-vort1_big_d(jx,jy,jzHi)))
    else
        cc_big_d(jx,jy,jz) = const * (u_big_d(jx,jy,jz) *                     &
            vort3_big_d(jx,jy,jz) + 0.5_rprec *                               &
            (w_big_d(jx,jy,jz+1) * (-vort1_big_d(jx,jy,jz+1)) +              &
            w_big_d(jx,jy,jz) * (-vort1_big_d(jx,jy,jz))))
    end if
end do
end do
end do
call check_convec_cuda('rhsy product')

istat = cufftExecD2Z(plan_fw_b_nz1, cc_big_d(:,:,1), cc_big_d(:,:,1))
call require_convec_cufft('rhsy forward big', istat)
call unpadd_3d_dp(RHSy, cc_big_d, lbz, 1, 1, nz-1)
istat = cufftExecZ2D(plan_bk_s_nz1, RHSy(:,:,1), RHSy(:,:,1))
call require_convec_cufft('rhsy inverse small', istat)

! RHSz = u*w2 - v*w1 on the padded grid, then truncate.
if (coord == 0) then
    jz_min = 2
else
    jz_min = 1
end if
jz_max = nz - 1

!$cuf kernel do(3) <<<*,*>>>
do jz = 1, nz
do jy = 1, ny2
do jx = 1, ld_big
    if (coord == 0 .and. jz == 1) then
        cc_big_d(jx,jy,1) = 0._rprec
    else if (coord == nproc-1 .and. jz == nz) then
        cc_big_d(jx,jy,nz) = 0._rprec
    else if (jz >= jz_min .and. jz <= jz_max) then
        cc_big_d(jx,jy,jz) = const * 0.5_rprec *                              &
            ((u_big_d(jx,jy,jz) + u_big_d(jx,jy,jz-1)) *                     &
            (-vort2_big_d(jx,jy,jz)) +                                        &
            (v_big_d(jx,jy,jz) + v_big_d(jx,jy,jz-1)) *                      &
            vort1_big_d(jx,jy,jz))
    else
        cc_big_d(jx,jy,jz) = 0._rprec
    end if
end do
end do
end do
call check_convec_cuda('rhsz product')

istat = cufftExecD2Z(plan_fw_b_nz, cc_big_d(:,:,1), cc_big_d(:,:,1))
call require_convec_cufft('rhsz forward big', istat)
call unpadd_3d_dp(RHSz, cc_big_d, lbz, 1, 1, nz)
istat = cufftExecZ2D(plan_bk_s_nz, RHSz(:,:,1), RHSz(:,:,1))
call require_convec_cufft('rhsz inverse small', istat)
call check_convec_cuda('rhs inverse small')

#ifdef PPMPI
#ifdef PPSAFETYMODE
!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, ld
    RHSx(jx, jy, 0) = BOGUS
    RHSy(jx, jy, 0) = BOGUS
    RHSz(jx, jy, 0) = BOGUS
end do
end do
#endif
#endif

#ifdef PPSAFETYMODE
!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, ld
    RHSx(jx, jy, nz) = BOGUS
    RHSy(jx, jy, nz) = BOGUS
    if (coord < nproc-1) RHSz(jx, jy, nz) = BOGUS
end do
end do
#endif

call convec_cuda_sync('convec final')

end subroutine convec_cuda_impl

!*******************************************************************************
subroutine ensure_convec_cuda()
!*******************************************************************************
implicit none

integer :: istat
integer :: n_s(2), inem_s(2), onem_s(2)
integer :: n_b(2), inem_b(2), onem_b(2)

if (.not. allocated(cc_big_d)) then
    allocate(cc_big_d(ld_big, ny2, nz))
    allocate(u_big_d(ld_big, ny2, lbz:nz))
    allocate(v_big_d(ld_big, ny2, lbz:nz))
    allocate(w_big_d(ld_big, ny2, lbz:nz))
    allocate(vort1_big_d(ld_big, ny2, nz))
    allocate(vort2_big_d(ld_big, ny2, nz))
    allocate(vort3_big_d(ld_big, ny2, nz))
end if

if (convec_cuda_initialized) return

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = ld / 2

n_b(1) = ny2
n_b(2) = nx2
inem_b(1) = ny2
inem_b(2) = ld_big
onem_b(1) = ny2
onem_b(2) = ld_big / 2

istat = cufftPlanMany(plan_fw_s_lbz, 2, n_s, inem_s, 1, ld*ny,               &
    onem_s, 1, (ld/2)*ny, CUFFT_D2Z, nz-lbz+1)
call require_convec_cufft('small forward lbz plan', istat)

istat = cufftPlanMany(plan_bk_b_lbz, 2, n_b, onem_b, 1, (ld_big/2)*ny2,       &
    inem_b, 1, ld_big*ny2, CUFFT_Z2D, nz-lbz+1)
call require_convec_cufft('big inverse lbz plan', istat)

istat = cufftPlanMany(plan_fw_s_nz, 2, n_s, inem_s, 1, ld*ny,                &
    onem_s, 1, (ld/2)*ny, CUFFT_D2Z, nz)
call require_convec_cufft('small forward nz plan', istat)

istat = cufftPlanMany(plan_bk_b_nz, 2, n_b, onem_b, 1, (ld_big/2)*ny2,        &
    inem_b, 1, ld_big*ny2, CUFFT_Z2D, nz)
call require_convec_cufft('big inverse nz plan', istat)

istat = cufftPlanMany(plan_fw_b_nz, 2, n_b, inem_b, 1, ld_big*ny2,            &
    onem_b, 1, (ld_big/2)*ny2, CUFFT_D2Z, nz)
call require_convec_cufft('big forward nz plan', istat)

istat = cufftPlanMany(plan_bk_s_nz, 2, n_s, onem_s, 1, (ld/2)*ny,             &
    inem_s, 1, ld*ny, CUFFT_Z2D, nz)
call require_convec_cufft('small inverse nz plan', istat)

istat = cufftPlanMany(plan_fw_b_nz1, 2, n_b, inem_b, 1, ld_big*ny2,           &
    onem_b, 1, (ld_big/2)*ny2, CUFFT_D2Z, nz-1)
call require_convec_cufft('big forward nz-1 plan', istat)

istat = cufftPlanMany(plan_bk_s_nz1, 2, n_s, onem_s, 1, (ld/2)*ny,            &
    inem_s, 1, ld*ny, CUFFT_Z2D, nz-1)
call require_convec_cufft('small inverse nz-1 plan', istat)

convec_cuda_initialized = .true.

end subroutine ensure_convec_cuda

!*******************************************************************************
subroutine padd_3d_dp(u_big, u_small, kb_big, kb_small, k_start, k_end)
!*******************************************************************************
implicit none

integer, intent(in) :: kb_big, kb_small, k_start, k_end
real(rprec), intent(inout) :: u_big(ld_big, ny2, kb_big:*)
real(rprec), intent(in) :: u_small(ld, ny, kb_small:*)
attributes(device) :: u_big, u_small
integer :: i, j, k
integer :: ny_h, j_s, j_big_s

ny_h = ny / 2
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

if (convec_fused_pad_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 1, ny2
    do i = 1, ld_big
        if (i <= nx .and. j <= ny_h) then
            u_big(i,j,k) = u_small(i,j,k)
        else if (i <= nx .and. j >= j_big_s) then
            u_big(i,j,k) = u_small(i,j - j_big_s + j_s,k)
        else
            u_big(i,j,k) = 0._rprec
        end if
    end do
    end do
    end do
else
    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 1, ny2
    do i = 1, ld_big
        u_big(i,j,k) = 0._rprec
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 1, ny_h
    do i = 1, nx
        u_big(i,j,k) = u_small(i,j,k)
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 0, ny - j_s
    do i = 1, nx
        u_big(i,j_big_s+j,k) = u_small(i,j_s+j,k)
    end do
    end do
    end do
end if

call check_convec_cuda('padd')

end subroutine padd_3d_dp

!*******************************************************************************
subroutine unpadd_3d_dp(cc, cc_big, kb_cc, kb_big, k_start, k_end)
!*******************************************************************************
implicit none

integer, intent(in) :: kb_cc, kb_big, k_start, k_end
real(rprec), intent(inout) :: cc(ld, ny, kb_cc:*)
real(rprec), intent(in) :: cc_big(ld_big, ny2, kb_big:*)
attributes(device) :: cc, cc_big
integer :: i, j, k
integer :: ny_h, j_s, j_big_s

ny_h = ny / 2
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

if (convec_fused_pad_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 1, ny
    do i = 1, ld
        if (i >= ld - 1 .or. j == ny_h + 1) then
            cc(i,j,k) = 0._rprec
        else if (i <= nx .and. j <= ny_h) then
            cc(i,j,k) = cc_big(i,j,k)
        else if (i <= nx .and. j >= j_s) then
            cc(i,j,k) = cc_big(i,j_big_s + j - j_s,k)
        else
            cc(i,j,k) = 0._rprec
        end if
    end do
    end do
    end do
else
    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 1, ny_h
    do i = 1, nx
        cc(i,j,k) = cc_big(i,j,k)
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 0, ny - j_s
    do i = 1, nx
        cc(i,j_s+j,k) = cc_big(i,j_big_s+j,k)
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do k = k_start, k_end
    do j = 1, ny
    do i = ld-1, ld
        cc(i,j,k) = 0._rprec
    end do
    end do
    end do

    !$cuf kernel do(2) <<<*,*>>>
    do k = k_start, k_end
    do i = 1, ld
        cc(i,ny_h+1,k) = 0._rprec
    end do
    end do
end if

call check_convec_cuda('unpadd')

end subroutine unpadd_3d_dp

!*******************************************************************************
subroutine check_convec_cuda(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

if (convec_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    if (istat /= 0) then
        print *, 'convec_cuda CUDA failure at ', trim(where), ': ', istat
        stop
    end if
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'convec_cuda kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine check_convec_cuda

!*******************************************************************************
subroutine convec_cuda_sync(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'convec_cuda CUDA failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'convec_cuda kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine convec_cuda_sync

!*******************************************************************************
subroutine require_convec_cufft(where, istat)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer, intent(in) :: istat

if (istat /= CUFFT_SUCCESS) then
    print *, 'convec_cuda cuFFT failure at ', trim(where), ': ', istat
    stop
end if

end subroutine require_convec_cufft

end module convec_cuda_m
#endif

!*******************************************************************************
subroutine convec
!*******************************************************************************
!
! Computes the rotation convective term in physical space
!       c = - (u X vort)
! Uses 3/2-rule for dealiasing for more info see Canuto 1991 Spectral Methods,
! chapter 7
!
use types, only : rprec
use param
use sim_param, only : u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy
use sim_param, only : RHSx, RHSy, RHSz
use fft
#ifdef ENABLE_CUDA
use convec_cuda_m, only : convec_cuda_enabled, convec_cuda_impl
#endif

implicit none

integer :: jz
integer :: jz_min
integer :: jzLo, jzHi, jz_max  ! added for full channel capabilities

real(rprec), save, allocatable, dimension(:,:,:) :: cc_big,                    &
    u_big, v_big, w_big, vort1_big, vort2_big, vort3_big
logical, save :: arrays_allocated = .false.

real(rprec) :: const

#ifdef ENABLE_CUDA
if (convec_cuda_enabled()) then
    call convec_cuda_impl()
    return
end if
#endif

! Boundary vorticity index used in the wall cross-product terms below.
! LES/SGS uses the first interior uvp slab at the lower wall; DNS/no-SGS keeps
! the boundary slab.  The upper-wall term uses nz-1 in both cases.
if (sgs) then
    jzLo = 2
    jzHi = nz-1
else
    jzLo = 1
    jzHi = nz-1
endif

if( .not. arrays_allocated ) then
   allocate( cc_big( ld_big,ny2,nz ) )
   allocate( u_big(ld_big, ny2, lbz:nz) )
   allocate( v_big(ld_big, ny2, lbz:nz) )
   allocate( w_big(ld_big, ny2, lbz:nz) )
   allocate( vort1_big( ld_big,ny2,nz ) )
   allocate( vort2_big( ld_big,ny2,nz ) )
   allocate( vort3_big( ld_big,ny2,nz ) )
   arrays_allocated = .true.
endif

! Recall dudz, and dvdz are on UVP node for k=1 only
! So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! Loop through horizontal slices
! MPI: u_big, v_big needed at jz = 0, w_big not needed though
! MPI: could get u{1,2}_big
const = 1._rprec/(nx*ny)
do jz = lbz, nz
    ! use RHSx,RHSy,RHSz for temp storage
    RHSx(:,:,jz)=const*u(:,:,jz)
    RHSy(:,:,jz)=const*v(:,:,jz)
    RHSz(:,:,jz)=const*w(:,:,jz)

    ! do forward fft on normal-size arrays
    call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))

    ! zero pad: padd takes care of the oddballs
    call padd(u_big(:,:,jz), RHSx(:,:,jz))
    call padd(v_big(:,:,jz), RHSy(:,:,jz))
    call padd(w_big(:,:,jz), RHSz(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back_big, u_big(:,:,jz), u_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, v_big(:,:,jz), v_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, w_big(:,:,jz), w_big(:,:,jz))
end do


! Do the same thing with the vorticity
do jz = 1, nz
    ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
    ! definition of the vorticity in that case which also interpolates
    ! dwdx, dwdy to the u-node at jz=1
    if ( (coord == 0) .and. (jz == 1) ) then

        select case (lbc_mom)
        ! Stress free
        case (0)
            RHSx(:, :, 1) = 0._rprec
            RHSy(:, :, 1) = 0._rprec

        ! Wall (all cases >= 1)
        case (1:)
            ! Wall-node vorticity: interpolate dwdy to the first uvp node.
            RHSx(:, :, 1) = const * ( 0.5_rprec * (dwdy(:, :, 1) +             &
                dwdy(:, :, 2))  - dvdz(:, :, 1) )
            ! Wall-node vorticity: interpolate dwdx to the first uvp node.
            RHSy(:, :, 1) = const * ( dudz(:, :, 1) -                          &
                0.5_rprec * (dwdx(:, :, 1) + dwdx(:, :, 2)) )

        end select
  endif

  if ( (coord == nproc-1) .and. (jz == nz) ) then

     select case (ubc_mom)

     ! Stress free
     case (0)

         RHSx(:, :, nz) = 0._rprec
         RHSy(:, :, nz) = 0._rprec

      ! No-slip and wall model
      case (1:)

         ! RHSx = vort1 at uvp nz-1, stored in the w-node nz work slot.
         RHSx(:, :, nz) = const * ( 0.5_rprec * (dwdy(:, :, nz-1) +            &
            dwdy(:, :, nz)) - dvdz(:, :, nz-1) )
         ! RHSy = vort2 at uvp nz-1, stored in the w-node nz work slot.
         RHSy(:, :, nz) = const * ( dudz(:, :, nz-1) -                         &
            0.5_rprec * (dwdx(:, :, nz-1) + dwdx(:, :, nz)) )

      end select
   endif

    ! Boundary slabs with explicit wall vorticity formulas were handled above;
    ! do not overwrite them with the general interior expression.
    if (.not.(coord==0 .and. jz==1) .and. .not. (ubc_mom>0 .and.               &
        coord==nproc-1 .and. jz==nz)  ) then
        RHSx(:,:,jz)=const*(dwdy(:,:,jz)-dvdz(:,:,jz))
        RHSy(:,:,jz)=const*(dudz(:,:,jz)-dwdx(:,:,jz))
    end if

    RHSz(:,:,jz)=const*(dvdx(:,:,jz)-dudy(:,:,jz))

    ! do forward fft on normal-size arrays
    call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    call padd(vort1_big(:,:,jz), RHSx(:,:,jz))
    call padd(vort2_big(:,:,jz), RHSy(:,:,jz))
    call padd(vort3_big(:,:,jz), RHSz(:,:,jz))

    ! Back to physical space
    ! FFTW backward transforms are unnormalized; const is applied when the
    ! padded-grid products are assembled below.
    call dfftw_execute_dft_c2r(back_big, vort1_big(:,:,jz), vort1_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, vort2_big(:,:,jz), vort2_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, vort3_big(:,:,jz), vort3_big(:,:,jz))
end do

! RHSx
! redefinition of const
const=1._rprec/(nx2*ny2)

if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,1)=const*(v_big(:,:,1)*(-vort3_big(:,:,1))&
       +0.5_rprec*w_big(:,:,2)*(vort2_big(:,:,jzLo)))
    ! jzLo selects the LES/DNS lower-wall vorticity slab described above.
    ! The 0.5 factor interpolates w(:,:,2) to the first uvp node above the
    ! wall; changing it is a wall-discretization change that needs validation.
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1 ) then  ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,nz-1)=const*(v_big(:,:,nz-1)*(-vort3_big(:,:,nz-1))&
        +0.5_rprec*w_big(:,:,nz-1)*(vort2_big(:,:,jzHi)))
    ! jzHi selects the upper-wall vorticity slab described above.
    ! The 0.5 factor interpolates w(:,:,nz-1) to the uvp node below the wall;
    ! changing it is a wall-discretization change that needs validation.

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max    !nz-1   ! channel
    cc_big(:,:,jz)=const*(v_big(:,:,jz)*(-vort3_big(:,:,jz))&
        +0.5_rprec*(w_big(:,:,jz+1)*(vort2_big(:,:,jz+1))&
        +w_big(:,:,jz)*(vort2_big(:,:,jz))))
end do

! Loop through horizontal slices
do jz=1,nz-1
    call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz),cc_big(:,:,jz))
    ! un-zero pad
    ! note: cc_big is going into RHSx
    call unpadd(RHSx(:,:,jz),cc_big(:,:,jz))
    ! Back to physical space
    call dfftw_execute_dft_c2r(back, RHSx(:,:,jz), RHSx(:,:,jz))
end do

! RHSy
! const remains the padded-grid inverse transform normalization.
if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,1)=const*(u_big(:,:,1)*(vort3_big(:,:,1))&
        +0.5_rprec*w_big(:,:,2)*(-vort1_big(:,:,jzLo)))
    ! jzLo selects the LES/DNS lower-wall vorticity slab described above.
    ! The 0.5 factor interpolates w(:,:,2) to the first uvp node above the
    ! wall; changing it is a wall-discretization change that needs validation.
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then   ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,nz-1)=const*(u_big(:,:,nz-1)*(vort3_big(:,:,nz-1))&
        +0.5_rprec*w_big(:,:,nz-1)*(-vort1_big(:,:,jzHi)))
    ! jzHi selects the upper-wall vorticity slab described above.
    !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    !  below the wall

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max  !nz - 1   ! channel
   cc_big(:,:,jz)=const*(u_big(:,:,jz)*(vort3_big(:,:,jz))&
        +0.5_rprec*(w_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
        +w_big(:,:,jz)*(-vort1_big(:,:,jz))))
end do

do jz=1,nz-1
    call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz), cc_big(:,:,jz))
    ! un-zero pad
    ! note: cc_big is going into RHSy
    call unpadd(RHSy(:,:,jz), cc_big(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back, RHSy(:,:,jz), RHSy(:,:,jz))
end do

! RHSz

if (coord == 0) then
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
    !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
    ! slip BC
    cc_big(:,:,1)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then     ! channel
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--but now we do solve z-eqn (w-eqn) at top wall --pj
    !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
    ! slip BC
    cc_big(:,:,nz)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_max = nz-1
else
    jz_max = nz-1   !! or nz ?       ! channel
end if

!#ifdef PPMPI
!  if (coord == nproc-1) then
!    cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!    jz_max = nz - 1
!  else
!    jz_max = nz
!  endif
!#else
!  cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!  jz_max = nz - 1
!#endif

! channel
do jz = jz_min, jz_max    !nz - 1
    cc_big(:,:,jz) = const*0.5_rprec*(                                         &
        (u_big(:,:,jz)+u_big(:,:,jz-1))*(-vort2_big(:,:,jz))                   &
        +(v_big(:,:,jz)+v_big(:,:,jz-1))*(vort1_big(:,:,jz)))
end do

! Loop through horizontal slices
do jz=1,nz !nz - 1
    call dfftw_execute_dft_r2c(forw_big,cc_big(:,:,jz),cc_big(:,:,jz))

    ! un-zero pad
    ! note: cc_big is going into RHSz!!!!
    call unpadd(RHSz(:,:,jz),cc_big(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back,RHSz(:,:,jz),   RHSz(:,:,jz))
end do

#ifdef PPMPI
#ifdef PPSAFETYMODE
RHSx(:, :, 0) = BOGUS
RHSy(:, :, 0) = BOGUS
RHSz(: ,:, 0) = BOGUS
#endif
#endif

!--top level is not valid
#ifdef PPSAFETYMODE
RHSx(:, :, nz) = BOGUS
RHSy(:, :, nz) = BOGUS
if(coord<nproc-1) RHSz(:, :, nz) = BOGUS
#endif

end subroutine convec
#endif
