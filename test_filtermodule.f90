!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!7
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
module test_filtermodule
!*******************************************************************************
! Navigation map:
!   - filter state: G_test/G_test_test and optional CUDA/OpenACC work buffers
!   - CUDA Fortran helper path: ensure_*_cuda_plan and apply_*_cuda routines
!   - initialization: test_filter_init
!   - CPU/managed wrappers: test_filter*, test_test_filter*
!   - OpenACC production wrappers: test_filter_*_gpu and test_test_filter_*_gpu
!
! This module is shared by SGS, derivatives, and wall-model paths.  Keep the
! CPU, CUDA Fortran, and OpenACC wrappers semantically paired when changing
! filter coefficients or padding behavior.
use types, only : rprec
use param, only : lh, ny
#ifdef ENABLE_CUDA
use cudafor
#endif

private lh, ny

! the implicit filter (1=grid size)
integer, parameter :: filter_size=1
! alpha is ratio of test filter to grid filter widths
real(rprec) :: alpha_test = 2.0_rprec * filter_size
real(rprec) :: alpha_test_test = 4.0_rprec * filter_size
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:), allocatable :: G_test, G_test_test
real(rprec), device, dimension(:,:), allocatable :: G_test_d, G_test_test_d
real(rprec), device, dimension(:,:), allocatable :: test_filter_work
integer, save :: test_filter_fw_plan = 0
integer, save :: test_filter_bk_plan = 0
integer, save :: test_filter_fw_plan_many(6) = 0
integer, save :: test_filter_bk_plan_many(6) = 0
integer, save :: test_filter_fw_plan_12 = 0
integer, save :: test_filter_bk_plan_12 = 0
logical, save :: test_filter_cuda_initialized = .false.
logical, save :: test_filter_cuda_many_initialized(6) = .false.
logical, save :: test_filter_cuda_12_initialized = .false.
real(rprec), device, dimension(:,:,:), allocatable :: test_filter_work_many
real(rprec), device, dimension(:,:,:), allocatable :: test_filter_work_many12
#else
real(rprec), dimension(:,:), allocatable :: G_test, G_test_test
#endif

#ifdef PPSGS_GPU
! OpenACC explicit-residency filter kernels used by the optimized LES GPU
! route. The existing CUDA Fortran managed filter path remains available when
! PPSGS_GPU is not enabled.
!$acc declare create(G_test, G_test_test)
#endif

contains

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function test_filter_cuda_enabled()
!*******************************************************************************
implicit none

test_filter_cuda_enabled = .true.

end function test_filter_cuda_enabled

!*******************************************************************************
logical function test_filter_extra_sync_enabled()
!*******************************************************************************
implicit none

test_filter_extra_sync_enabled = .false.

end function test_filter_extra_sync_enabled

!*******************************************************************************
subroutine require_test_filter_cufft_success(where, istat)
!*******************************************************************************
use cufft
implicit none

character(len=*), intent(in) :: where
integer, intent(in) :: istat

if (istat /= CUFFT_SUCCESS) then
    print *, 'test_filter cuFFT failure at ', trim(where), ': ', istat
    stop
end if

end subroutine require_test_filter_cufft_success

!*******************************************************************************
subroutine test_filter_cuda_sync(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

if (test_filter_extra_sync_enabled()) then
    istat = cudaDeviceSynchronize()
    if (istat /= 0) then
        print *, 'test_filter CUDA sync failure at ', trim(where), ': ', istat
        stop
    end if
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'test_filter CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine test_filter_cuda_sync

!*******************************************************************************
subroutine test_filter_cuda_barrier(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'test_filter CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'test_filter CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine test_filter_cuda_barrier

!*******************************************************************************
subroutine ensure_test_filter_cuda_plan()
!*******************************************************************************
use param, only : ld, nx, ny
use cufft
implicit none

integer :: istat
integer :: n_s(2), inem_s(2), onem_s(2)

if (test_filter_cuda_initialized) return

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(test_filter_fw_plan, 2, n_s, inem_s, 1, ld*ny,          &
    onem_s, 1, lh*ny, CUFFT_D2Z, 1)
call require_test_filter_cufft_success('forward plan', istat)

istat = cufftPlanMany(test_filter_bk_plan, 2, n_s, onem_s, 1, lh*ny,          &
    inem_s, 1, ld*ny, CUFFT_Z2D, 1)
call require_test_filter_cufft_success('inverse plan', istat)

test_filter_cuda_initialized = .true.

end subroutine ensure_test_filter_cuda_plan

!*******************************************************************************
subroutine ensure_test_filter_cuda_many_plan(nbatch)
!*******************************************************************************
use param, only : ld, nx, ny
use cufft
implicit none

integer, intent(in) :: nbatch
integer :: istat
integer :: n_s(2), inem_s(2), onem_s(2)

if (test_filter_cuda_many_initialized(nbatch)) return

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(test_filter_fw_plan_many(nbatch), 2, n_s, inem_s, 1,    &
    ld*ny, onem_s, 1, lh*ny, CUFFT_D2Z, nbatch)
call require_test_filter_cufft_success('batched forward plan', istat)

istat = cufftPlanMany(test_filter_bk_plan_many(nbatch), 2, n_s, onem_s, 1,    &
    lh*ny, inem_s, 1, ld*ny, CUFFT_Z2D, nbatch)
call require_test_filter_cufft_success('batched inverse plan', istat)

test_filter_cuda_many_initialized(nbatch) = .true.

end subroutine ensure_test_filter_cuda_many_plan

!*******************************************************************************
subroutine ensure_test_filter_cuda_12_plan()
!*******************************************************************************
use param, only : ld, nx, ny
use cufft
implicit none

integer :: istat
integer :: n_s(2), inem_s(2), onem_s(2)

if (test_filter_cuda_12_initialized) return

n_s(1) = ny
n_s(2) = nx
inem_s(1) = ny
inem_s(2) = ld
onem_s(1) = ny
onem_s(2) = lh

istat = cufftPlanMany(test_filter_fw_plan_12, 2, n_s, inem_s, 1,             &
    ld*ny, onem_s, 1, lh*ny, CUFFT_D2Z, 12)
call require_test_filter_cufft_success('batched-12 forward plan', istat)

istat = cufftPlanMany(test_filter_bk_plan_12, 2, n_s, onem_s, 1,             &
    lh*ny, inem_s, 1, ld*ny, CUFFT_Z2D, 12)
call require_test_filter_cufft_success('batched-12 inverse plan', istat)

test_filter_cuda_12_initialized = .true.

end subroutine ensure_test_filter_cuda_12_plan

!*******************************************************************************
    subroutine apply_test_filter_cuda(f, G_d, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), dimension(:,:), intent(inout) :: f
    real(rprec), device, dimension(:,:), intent(in) :: G_d
character(len=*), intent(in) :: where
integer :: i, j, istat

call ensure_test_filter_cuda_plan()
if (.not. allocated(test_filter_work)) allocate(test_filter_work(ld, ny))
test_filter_work = f

istat = cufftExecD2Z(test_filter_fw_plan, test_filter_work, test_filter_work)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, lh
        test_filter_work(2*i-1,j) = test_filter_work(2*i-1,j) * G_d(i,j)
        test_filter_work(2*i,j) = test_filter_work(2*i,j) * G_d(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan, test_filter_work, test_filter_work)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_barrier(trim(where)//' inverse')
f = test_filter_work

end subroutine apply_test_filter_cuda

!*******************************************************************************
    subroutine apply_test_filter_cuda_3(f1, f2, f3, G_d, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), dimension(:,:), intent(inout) :: f1, f2, f3
    real(rprec), device, dimension(:,:), intent(in) :: G_d
character(len=*), intent(in) :: where
integer :: i, j, b, istat

call ensure_test_filter_cuda_many_plan(3)
if (.not. allocated(test_filter_work_many)) allocate(test_filter_work_many(ld, ny, 6))
test_filter_work_many(:,:,1) = f1
test_filter_work_many(:,:,2) = f2
test_filter_work_many(:,:,3) = f3
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan_many(3), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(3) <<<*,*>>>
do b = 1, 3
do j = 1, ny
do i = 1, lh
        test_filter_work_many(2*i-1,j,b) = test_filter_work_many(2*i-1,j,b) * G_d(i,j)
        test_filter_work_many(2*i,j,b) = test_filter_work_many(2*i,j,b) * G_d(i,j)
end do
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan_many(3), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_barrier(trim(where)//' inverse')

f1 = test_filter_work_many(:,:,1)
f2 = test_filter_work_many(:,:,2)
f3 = test_filter_work_many(:,:,3)
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_3

!*******************************************************************************
    subroutine apply_test_filter_cuda_6(f1, f2, f3, f4, f5, f6, G_d, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6
    real(rprec), device, dimension(:,:), intent(in) :: G_d
character(len=*), intent(in) :: where
integer :: i, j, b, istat

call ensure_test_filter_cuda_many_plan(6)
if (.not. allocated(test_filter_work_many)) allocate(test_filter_work_many(ld, ny, 6))
test_filter_work_many(:,:,1) = f1
test_filter_work_many(:,:,2) = f2
test_filter_work_many(:,:,3) = f3
test_filter_work_many(:,:,4) = f4
test_filter_work_many(:,:,5) = f5
test_filter_work_many(:,:,6) = f6
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan_many(6), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(3) <<<*,*>>>
do b = 1, 6
do j = 1, ny
do i = 1, lh
        test_filter_work_many(2*i-1,j,b) = test_filter_work_many(2*i-1,j,b) * G_d(i,j)
        test_filter_work_many(2*i,j,b) = test_filter_work_many(2*i,j,b) * G_d(i,j)
end do
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan_many(6), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_barrier(trim(where)//' inverse')

f1 = test_filter_work_many(:,:,1)
f2 = test_filter_work_many(:,:,2)
f3 = test_filter_work_many(:,:,3)
f4 = test_filter_work_many(:,:,4)
f5 = test_filter_work_many(:,:,5)
f6 = test_filter_work_many(:,:,6)
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_6

!*******************************************************************************
    subroutine apply_test_filter_cuda_managed(f, G_d, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f
    real(rprec), device, dimension(:,:), intent(in) :: G_d
character(len=*), intent(in) :: where
integer :: i, j, istat

call ensure_test_filter_cuda_plan()
if (.not. allocated(test_filter_work)) allocate(test_filter_work(ld, ny))
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    test_filter_work(i,j) = f(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan, test_filter_work, test_filter_work)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, lh
        test_filter_work(2*i-1,j) = test_filter_work(2*i-1,j) * G_d(i,j)
        test_filter_work(2*i,j) = test_filter_work(2*i,j) * G_d(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan, test_filter_work, test_filter_work)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_sync(trim(where)//' inverse')

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    f(i,j) = test_filter_work(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_managed

!*******************************************************************************
    subroutine apply_test_filter_cuda_3_managed(f1, f2, f3, G_d, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3
    real(rprec), device, dimension(:,:), intent(in) :: G_d
character(len=*), intent(in) :: where
integer :: i, j, b, istat

call ensure_test_filter_cuda_many_plan(3)
if (.not. allocated(test_filter_work_many)) allocate(test_filter_work_many(ld, ny, 6))
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    test_filter_work_many(i,j,1) = f1(i,j)
    test_filter_work_many(i,j,2) = f2(i,j)
    test_filter_work_many(i,j,3) = f3(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan_many(3), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(3) <<<*,*>>>
do b = 1, 3
do j = 1, ny
do i = 1, lh
        test_filter_work_many(2*i-1,j,b) = test_filter_work_many(2*i-1,j,b) * G_d(i,j)
        test_filter_work_many(2*i,j,b) = test_filter_work_many(2*i,j,b) * G_d(i,j)
end do
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan_many(3), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_sync(trim(where)//' inverse')

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    f1(i,j) = test_filter_work_many(i,j,1)
    f2(i,j) = test_filter_work_many(i,j,2)
    f3(i,j) = test_filter_work_many(i,j,3)
end do
end do
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_3_managed

!*******************************************************************************
    subroutine apply_test_filter_cuda_6_managed(f1, f2, f3, f4, f5, f6, G_d, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6
    real(rprec), device, dimension(:,:), intent(in) :: G_d
character(len=*), intent(in) :: where
integer :: i, j, b, istat

call ensure_test_filter_cuda_many_plan(6)
if (.not. allocated(test_filter_work_many)) allocate(test_filter_work_many(ld, ny, 6))
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    test_filter_work_many(i,j,1) = f1(i,j)
    test_filter_work_many(i,j,2) = f2(i,j)
    test_filter_work_many(i,j,3) = f3(i,j)
    test_filter_work_many(i,j,4) = f4(i,j)
    test_filter_work_many(i,j,5) = f5(i,j)
    test_filter_work_many(i,j,6) = f6(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan_many(6), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(3) <<<*,*>>>
do b = 1, 6
do j = 1, ny
do i = 1, lh
        test_filter_work_many(2*i-1,j,b) = test_filter_work_many(2*i-1,j,b) * G_d(i,j)
        test_filter_work_many(2*i,j,b) = test_filter_work_many(2*i,j,b) * G_d(i,j)
end do
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan_many(6), test_filter_work_many,       &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_sync(trim(where)//' inverse')

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    f1(i,j) = test_filter_work_many(i,j,1)
    f2(i,j) = test_filter_work_many(i,j,2)
    f3(i,j) = test_filter_work_many(i,j,3)
    f4(i,j) = test_filter_work_many(i,j,4)
    f5(i,j) = test_filter_work_many(i,j,5)
    f6(i,j) = test_filter_work_many(i,j,6)
end do
end do
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_6_managed

!*******************************************************************************
subroutine apply_test_filter_cuda_3_dual_managed(f1, f2, f3, g1, g2, g3, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3
real(rprec), managed, dimension(:,:), intent(inout) :: g1, g2, g3
character(len=*), intent(in) :: where
integer :: i, j, b, istat

call ensure_test_filter_cuda_many_plan(6)
if (.not. allocated(test_filter_work_many)) allocate(test_filter_work_many(ld, ny, 6))
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    test_filter_work_many(i,j,1) = f1(i,j)
    test_filter_work_many(i,j,2) = f2(i,j)
    test_filter_work_many(i,j,3) = f3(i,j)
    test_filter_work_many(i,j,4) = g1(i,j)
    test_filter_work_many(i,j,5) = g2(i,j)
    test_filter_work_many(i,j,6) = g3(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan_many(6), test_filter_work_many,      &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(3) <<<*,*>>>
do b = 1, 6
do j = 1, ny
do i = 1, lh
    if (b <= 3) then
        test_filter_work_many(2*i-1,j,b) = test_filter_work_many(2*i-1,j,b)   &
            * G_test_d(i,j)
        test_filter_work_many(2*i,j,b) = test_filter_work_many(2*i,j,b)       &
            * G_test_d(i,j)
    else
        test_filter_work_many(2*i-1,j,b) = test_filter_work_many(2*i-1,j,b)   &
            * G_test_test_d(i,j)
        test_filter_work_many(2*i,j,b) = test_filter_work_many(2*i,j,b)       &
            * G_test_test_d(i,j)
    end if
end do
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan_many(6), test_filter_work_many,      &
    test_filter_work_many)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_sync(trim(where)//' inverse')

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    f1(i,j) = test_filter_work_many(i,j,1)
    f2(i,j) = test_filter_work_many(i,j,2)
    f3(i,j) = test_filter_work_many(i,j,3)
    g1(i,j) = test_filter_work_many(i,j,4)
    g2(i,j) = test_filter_work_many(i,j,5)
    g3(i,j) = test_filter_work_many(i,j,6)
end do
end do
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_3_dual_managed

!*******************************************************************************
subroutine apply_test_filter_cuda_6_dual_managed(f1, f2, f3, f4, f5, f6,      &
    g1, g2, g3, g4, g5, g6, where)
!*******************************************************************************
use param, only : ld, ny
use cufft
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6
real(rprec), managed, dimension(:,:), intent(inout) :: g1, g2, g3, g4, g5, g6
character(len=*), intent(in) :: where
integer :: i, j, b, istat

call ensure_test_filter_cuda_12_plan()
if (.not. allocated(test_filter_work_many12))                            &
    allocate(test_filter_work_many12(ld, ny, 12))
!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    test_filter_work_many12(i,j,1) = f1(i,j)
    test_filter_work_many12(i,j,2) = f2(i,j)
    test_filter_work_many12(i,j,3) = f3(i,j)
    test_filter_work_many12(i,j,4) = f4(i,j)
    test_filter_work_many12(i,j,5) = f5(i,j)
    test_filter_work_many12(i,j,6) = f6(i,j)
    test_filter_work_many12(i,j,7) = g1(i,j)
    test_filter_work_many12(i,j,8) = g2(i,j)
    test_filter_work_many12(i,j,9) = g3(i,j)
    test_filter_work_many12(i,j,10) = g4(i,j)
    test_filter_work_many12(i,j,11) = g5(i,j)
    test_filter_work_many12(i,j,12) = g6(i,j)
end do
end do
call test_filter_cuda_sync(trim(where)//' pack')

istat = cufftExecD2Z(test_filter_fw_plan_12, test_filter_work_many12,         &
    test_filter_work_many12)
call require_test_filter_cufft_success(trim(where)//' forward', istat)

!$cuf kernel do(3) <<<*,*>>>
do b = 1, 12
do j = 1, ny
do i = 1, lh
    if (b <= 6) then
        test_filter_work_many12(2*i-1,j,b) =                              &
            test_filter_work_many12(2*i-1,j,b) * G_test_d(i,j)
        test_filter_work_many12(2*i,j,b) =                                  &
            test_filter_work_many12(2*i,j,b) * G_test_d(i,j)
    else
        test_filter_work_many12(2*i-1,j,b) =                              &
            test_filter_work_many12(2*i-1,j,b) * G_test_test_d(i,j)
        test_filter_work_many12(2*i,j,b) =                                  &
            test_filter_work_many12(2*i,j,b) * G_test_test_d(i,j)
    end if
end do
end do
end do
call test_filter_cuda_sync(trim(where)//' multiply')

istat = cufftExecZ2D(test_filter_bk_plan_12, test_filter_work_many12,         &
    test_filter_work_many12)
call require_test_filter_cufft_success(trim(where)//' inverse', istat)
call test_filter_cuda_sync(trim(where)//' inverse')

!$cuf kernel do(2) <<<*,*>>>
do j = 1, ny
do i = 1, ld
    f1(i,j) = test_filter_work_many12(i,j,1)
    f2(i,j) = test_filter_work_many12(i,j,2)
    f3(i,j) = test_filter_work_many12(i,j,3)
    f4(i,j) = test_filter_work_many12(i,j,4)
    f5(i,j) = test_filter_work_many12(i,j,5)
    f6(i,j) = test_filter_work_many12(i,j,6)
    g1(i,j) = test_filter_work_many12(i,j,7)
    g2(i,j) = test_filter_work_many12(i,j,8)
    g3(i,j) = test_filter_work_many12(i,j,9)
    g4(i,j) = test_filter_work_many12(i,j,10)
    g5(i,j) = test_filter_work_many12(i,j,11)
    g6(i,j) = test_filter_work_many12(i,j,12)
end do
end do
call test_filter_cuda_sync(trim(where)//' unpack')

end subroutine apply_test_filter_cuda_6_dual_managed
#endif

!*******************************************************************************
subroutine test_filter_init()
!*******************************************************************************
! Creates the kernels which will be used for filtering the field
use types, only : rprec
use param, only : lh, nx, ny, dx, dy, pi, ifilter, sgs_model
use fft
implicit none

real(rprec) :: delta_test, kc2_test, delta_test_test, kc2_test_test

! Allocate the arrays
allocate( G_test(lh,ny) )

! Include the normalization for the forward FFT
G_test = 1._rprec/(nx*ny)

! Filter characteristic width
! "2d-delta", not full 3d one
delta_test = alpha_test * sqrt(dx*dy)

! Calculate the kernel
! spectral cutoff filter
if(ifilter==1) then
    if (sgs_model==6.OR.sgs_model==7) then
        print *, 'Use Gaussian or Top-hat filter for mixed models'
        stop
    endif
    kc2_test = (pi/(delta_test))**2
    where (real(k2, rprec) >= kc2_test) G_test = 0._rprec

! Gaussian filter
else if(ifilter==2) then
    G_test=exp(-(delta_test)**2*k2/(4._rprec*6._rprec))*G_test

! Top-hat (Box) filter
else if(ifilter==3) then
    G_test = (sin(kx*delta_test/2._rprec)*sin(ky*delta_test/2._rprec)+1E-8)/   &
        (kx*delta_test/2._rprec*ky*delta_test/2._rprec+1E-8)*G_test
endif

! since our k2 has zero at Nyquist, we have to do this by hand
G_test(lh,:) = 0._rprec
G_test(:,ny/2+1) = 0._rprec
#ifdef ENABLE_CUDA
allocate(G_test_d(lh,ny))
G_test_d = G_test
#endif
#ifdef PPSGS_GPU
!$acc update device(G_test)
#endif

! Second test filter, if necessary (with scale dependent dynamic)
if ((sgs_model == 3) .or. (sgs_model == 5)) then
    ! Allocate the arrays
    allocate ( G_test_test(lh,ny) )

    ! Include the normalization
    G_test_test = 1._rprec/(nx*ny)

    ! Filter characteristic width
    delta_test_test = alpha_test_test * sqrt(dx*dy)

    ! Calculate the kernel
    ! spectral cutoff filter
    if (ifilter==1) then
        if (sgs_model==6.OR.sgs_model==7) then
            print *, 'Use Gaussian or Top-hat filter for mixed models'
            stop
        endif

        kc2_test_test = (pi/(delta_test_test))**2
        where (real(k2, rprec) >= kc2_test_test) G_test_test = 0._rprec

    ! Gaussian filter
    else if(ifilter==2) then
        G_test_test=exp(-(delta_test_test)**2*k2/(4._rprec*6._rprec))          &
            * G_test_test

    ! Top-hat (Box) filter
    else if(ifilter==3) then
        G_test_test= (sin(kx*delta_test_test/2._rprec)                         &
            * sin(ky*delta_test_test/2._rprec)+1E-8)                           &
            / (kx*delta_test_test/2._rprec*ky*delta_test_test/2._rprec+1E-8)   &
            * G_test_test
    endif

    ! since our k2 has zero at Nyquist, we have to do this by hand
        G_test_test(lh,:) = 0._rprec
        G_test_test(:,ny/2+1) = 0._rprec
#ifdef ENABLE_CUDA
    allocate(G_test_test_d(lh,ny))
    G_test_test_d = G_test_test
#endif
#ifdef PPSGS_GPU
    !$acc update device(G_test_test)
#endif

endif

end subroutine test_filter_init

!*******************************************************************************
subroutine test_filter(f)
!*******************************************************************************
! note: this filters in-place, so input is ruined
use types, only : rprec
use fft
use param, only : ny
use emul_complex, only : OPERATOR(.MULR.)
implicit none

real(rprec), dimension(:,:), intent(inout) :: f

#ifdef ENABLE_CUDA
if (test_filter_cuda_enabled()) then
        call apply_test_filter_cuda(f, G_test_d, 'test_filter')
    return
end if
#endif

!  Perform in-place FFT
call dfftw_execute_dft_r2c(forw, f, f)

! Perform f = G_test*f, emulating f as complex
! Nyquist frequency and normalization is taken care of with G_test
f = f .MULR. G_test

call dfftw_execute_dft_c2r(back, f, f)

end subroutine test_filter

!*******************************************************************************
subroutine test_filter_3(f1, f2, f3)
!*******************************************************************************
use types, only : rprec
implicit none

real(rprec), dimension(:,:), intent(inout) :: f1, f2, f3

#ifdef ENABLE_CUDA
if (test_filter_cuda_enabled()) then
        call apply_test_filter_cuda_3(f1, f2, f3, G_test_d, 'test_filter_3')
    return
end if
#endif

call test_filter(f1)
call test_filter(f2)
call test_filter(f3)

end subroutine test_filter_3

!*******************************************************************************
subroutine test_filter_6(f1, f2, f3, f4, f5, f6)
!*******************************************************************************
use types, only : rprec
implicit none

real(rprec), dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6

#ifdef ENABLE_CUDA
if (test_filter_cuda_enabled()) then
        call apply_test_filter_cuda_6(f1, f2, f3, f4, f5, f6, G_test_d,             &
        'test_filter_6')
    return
end if
#endif

call test_filter(f1)
call test_filter(f2)
call test_filter(f3)
call test_filter(f4)
call test_filter(f5)
call test_filter(f6)

end subroutine test_filter_6

!*******************************************************************************
subroutine test_test_filter(f)
!*******************************************************************************
! note: this filters in-place, so input is ruined
use types, only : rprec
use fft
use param, only : ny
use emul_complex, only : OPERATOR(.MULR.)
implicit none

real(rprec), dimension(:,:), intent(inout) :: f

#ifdef ENABLE_CUDA
if (test_filter_cuda_enabled()) then
        call apply_test_filter_cuda(f, G_test_test_d, 'test_test_filter')
    return
end if
#endif

!  Perform in-place FFT
call dfftw_execute_dft_r2c(forw, f, f)

! Perform f = G_test*f, emulating f as complex
! Nyquist frequency and normalization is taken care of with G_test_test
f = f .MULR. G_test_test

call dfftw_execute_dft_c2r(back, f, f)

end subroutine test_test_filter

!*******************************************************************************
subroutine test_test_filter_3(f1, f2, f3)
!*******************************************************************************
use types, only : rprec
implicit none

real(rprec), dimension(:,:), intent(inout) :: f1, f2, f3

#ifdef ENABLE_CUDA
if (test_filter_cuda_enabled()) then
        call apply_test_filter_cuda_3(f1, f2, f3, G_test_test_d, 'test_test_filter_3')
    return
end if
#endif

call test_test_filter(f1)
call test_test_filter(f2)
call test_test_filter(f3)

end subroutine test_test_filter_3

!*******************************************************************************
subroutine test_test_filter_6(f1, f2, f3, f4, f5, f6)
!*******************************************************************************
use types, only : rprec
implicit none

real(rprec), dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6

#ifdef ENABLE_CUDA
if (test_filter_cuda_enabled()) then
        call apply_test_filter_cuda_6(f1, f2, f3, f4, f5, f6, G_test_test_d,        &
        'test_test_filter_6')
    return
end if
#endif

call test_test_filter(f1)
call test_test_filter(f2)
call test_test_filter(f3)
call test_test_filter(f4)
call test_test_filter(f5)
call test_test_filter(f6)

end subroutine test_test_filter_6

#ifdef PPSGS_GPU
!*******************************************************************************
subroutine test_filter_b_gpu(f, n)
!*******************************************************************************
! Batched OpenACC test filter for explicit-residency Lagrangian SGS data.
! f has shape (ld,ny,n); the current route uses n == nz.
use param, only : ny, nz
use fft_gpu, only : plan_forw_small_nz, plan_back_small_nz,                  &
                    fft_gpu_exec_d2z, fft_gpu_exec_z2d
implicit none

integer, intent(in) :: n
real(rprec), dimension(:,:,:), intent(inout) :: f
integer :: i, j, k

if (n /= nz) then
    write(*,*) 'test_filter_b_gpu: only batch=nz is wired; got n=', n
    error stop
end if

call fft_gpu_exec_d2z(plan_forw_small_nz, f)

!$acc parallel loop collapse(3) present(f, G_test) default(present) async(1)
do k = 1, n
do j = 1, ny
do i = 1, lh
    f(2*i-1, j, k) = f(2*i-1, j, k) * G_test(i, j)
    f(2*i,   j, k) = f(2*i,   j, k) * G_test(i, j)
end do
end do
end do

call fft_gpu_exec_z2d(plan_back_small_nz, f)

end subroutine test_filter_b_gpu

!*******************************************************************************
subroutine test_test_filter_b_gpu(f, n)
!*******************************************************************************
! Batched second-level OpenACC test filter for explicit-residency Lagrangian SGS.
use param, only : ny, nz
use fft_gpu, only : plan_forw_small_nz, plan_back_small_nz,                  &
                    fft_gpu_exec_d2z, fft_gpu_exec_z2d
implicit none

integer, intent(in) :: n
real(rprec), dimension(:,:,:), intent(inout) :: f
integer :: i, j, k

if (n /= nz) then
    write(*,*) 'test_test_filter_b_gpu: only batch=nz is wired; got n=', n
    error stop
end if

call fft_gpu_exec_d2z(plan_forw_small_nz, f)

!$acc parallel loop collapse(3) present(f, G_test_test) default(present) async(1)
do k = 1, n
do j = 1, ny
do i = 1, lh
    f(2*i-1, j, k) = f(2*i-1, j, k) * G_test_test(i, j)
    f(2*i,   j, k) = f(2*i,   j, k) * G_test_test(i, j)
end do
end do
end do

call fft_gpu_exec_z2d(plan_back_small_nz, f)

end subroutine test_test_filter_b_gpu

!*******************************************************************************
subroutine test_filter_plane_gpu(f)
!*******************************************************************************
! Single-plane (ld,ny) OpenACC test filter - used by the device-resident
! wallstress wall plane. f must be device-present.
use param, only : ny
use fft_gpu, only : plan_forw_small_one, plan_back_small_one,                 &
                    fft_gpu_exec_d2z, fft_gpu_exec_z2d
implicit none

real(rprec), dimension(:,:), intent(inout) :: f
integer :: i, j

call fft_gpu_exec_d2z(plan_forw_small_one, f)

!$acc parallel loop collapse(2) present(f, G_test) default(present) async(1)
do j = 1, ny
do i = 1, lh
    f(2*i-1, j) = f(2*i-1, j) * G_test(i, j)
    f(2*i,   j) = f(2*i,   j) * G_test(i, j)
end do
end do

call fft_gpu_exec_z2d(plan_back_small_one, f)

end subroutine test_filter_plane_gpu

!*******************************************************************************
subroutine test_test_filter_plane_gpu(f)
!*******************************************************************************
! Single-plane second-level OpenACC test filter. f must be device-present.
use param, only : ny
use fft_gpu, only : plan_forw_small_one, plan_back_small_one,                 &
                    fft_gpu_exec_d2z, fft_gpu_exec_z2d
implicit none

real(rprec), dimension(:,:), intent(inout) :: f
integer :: i, j

call fft_gpu_exec_d2z(plan_forw_small_one, f)

!$acc parallel loop collapse(2) present(f, G_test_test) default(present) async(1)
do j = 1, ny
do i = 1, lh
    f(2*i-1, j) = f(2*i-1, j) * G_test_test(i, j)
    f(2*i,   j) = f(2*i,   j) * G_test_test(i, j)
end do
end do

call fft_gpu_exec_z2d(plan_back_small_one, f)

end subroutine test_test_filter_plane_gpu
#endif

#ifdef ENABLE_CUDA
!*******************************************************************************
subroutine test_filter_gpu(f)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f

call apply_test_filter_cuda_managed(f, G_test_d, 'test_filter_gpu')

end subroutine test_filter_gpu

!*******************************************************************************
subroutine test_filter_3_gpu(f1, f2, f3)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3

call apply_test_filter_cuda_3_managed(f1, f2, f3, G_test_d, 'test_filter_3_gpu')

end subroutine test_filter_3_gpu

!*******************************************************************************
subroutine test_filter_6_gpu(f1, f2, f3, f4, f5, f6)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6

call apply_test_filter_cuda_6_managed(f1, f2, f3, f4, f5, f6, G_test_d,        &
    'test_filter_6_gpu')

end subroutine test_filter_6_gpu

!*******************************************************************************
subroutine test_test_filter_gpu(f)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f

call apply_test_filter_cuda_managed(f, G_test_test_d, 'test_test_filter_gpu')

end subroutine test_test_filter_gpu

!*******************************************************************************
subroutine test_test_filter_3_gpu(f1, f2, f3)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3

call apply_test_filter_cuda_3_managed(f1, f2, f3, G_test_test_d,               &
    'test_test_filter_3_gpu')

end subroutine test_test_filter_3_gpu

!*******************************************************************************
subroutine test_test_filter_6_gpu(f1, f2, f3, f4, f5, f6)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6

call apply_test_filter_cuda_6_managed(f1, f2, f3, f4, f5, f6, G_test_test_d,   &
    'test_test_filter_6_gpu')

end subroutine test_test_filter_6_gpu

!*******************************************************************************
subroutine test_filter_3_dual_gpu(f1, f2, f3, g1, g2, g3)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3
real(rprec), managed, dimension(:,:), intent(inout) :: g1, g2, g3

call apply_test_filter_cuda_3_dual_managed(f1, f2, f3, g1, g2, g3,            &
    'test_filter_3_dual_gpu')

end subroutine test_filter_3_dual_gpu

!*******************************************************************************
subroutine test_filter_6_dual_gpu(f1, f2, f3, f4, f5, f6, g1, g2, g3, g4,     &
    g5, g6)
!*******************************************************************************
implicit none

real(rprec), managed, dimension(:,:), intent(inout) :: f1, f2, f3, f4, f5, f6
real(rprec), managed, dimension(:,:), intent(inout) :: g1, g2, g3, g4, g5, g6

call apply_test_filter_cuda_6_dual_managed(f1, f2, f3, f4, f5, f6,            &
    g1, g2, g3, g4, g5, g6, 'test_filter_6_dual_gpu')

end subroutine test_filter_6_dual_gpu
#endif

end module test_filtermodule
