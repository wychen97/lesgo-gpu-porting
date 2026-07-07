!*******************************************************************************
! GPU FFT layer (cuFFT batched 2D plans).
!
! Mirrors the FFTW plan handles in `fft.f90` but as cuFFT batched plans that
! transform all z-slabs of a (ld, ny, batch) array in a single launch.
!
! Only compiled and linked when USE_LES_GPU=ON (CMake flag) -> CPP macro
! PPLES_GPU is defined. On gfortran/intel/nvfortran-CPU builds this module is
! never seen.
!
! Plan inventory (created once in init_fft_gpu, destroyed in finalize):
!
!                       small (nx x ny)              big (nx2 x ny2)
!                       D2Z forward, Z2D backward    D2Z forward, Z2D backward
!   batch = nz_full     plan_*_small_full            plan_*_big_full
!     (= nz - lbz + 1)  used by Phase A
!   batch = nz          plan_*_small_nz              plan_*_big_nz
!                       used by Phase B and E
!   batch = nz - 1      plan_*_small_nzm1            plan_*_big_nzm1
!                       used by Phases C and D
!
! 4 sizes/dirs x 3 batch sizes = 12 plans. cuFFT plans are cheap (~microseconds to
! create, KB of state); enumerating them removes runtime branching from the
! exec path.
!
! Layout convention (matches FFTW in-place r2c on a (ld, ny) Fortran array
! with ld = 2*(nx/2+1)):
!   - real input  : (ld, ny, batch)  with jx fastest, batch slowest
!   - complex out : same memory, viewed as (lh, ny, batch) double-complex
!   - cuFFT seen as C-row-major: outer dim ny, inner dim ld (real) / lh (cmplx)
!
! Stream: cuFFT calls are bound to OpenACC async queue 1 via
! cufftSetStream(plan, acc_get_cuda_stream(1)) in init_fft_gpu. Every
! !$acc parallel loop in *_gpu.f90 uses `async(1)`, so cuFFT and OpenACC
! kernels serialize through CUDA stream-1 ordering without any !$acc wait
! between them. Explicit `!$acc wait(1)` is only needed before host-visible
! operations (acc update self/device, MPI sendrecv via host_data,
! end-of-routine). A numbered queue is used (not acc_async_noval) so that
! wait(1) maps to cudaStreamSynchronize, not the ~10x costlier cuCtxSynchronize.
!*******************************************************************************
module fft_gpu
#ifdef PPLES_GPU
use types, only : rprec
use param, only : ld, ny, ld_big, ny2, nx, nx2, lbz, nz
use cufft, only : CUFFT_D2Z, CUFFT_SUCCESS, CUFFT_Z2D, cufftDestroy,          &
                  cufftExecD2Z, cufftExecZ2D, cufftPlanMany, cufftSetStream
use openacc, only : acc_get_cuda_stream, acc_handle_kind
implicit none
save

private
public :: init_fft_gpu, finalize_fft_gpu
public :: fft_gpu_exec_d2z, fft_gpu_exec_z2d
public :: plan_forw_small_full, plan_back_small_full
public :: plan_forw_big_full,   plan_back_big_full
public :: plan_forw_small_nz,   plan_back_small_nz
public :: plan_forw_big_nz,     plan_back_big_nz
public :: plan_forw_small_nzm1, plan_back_small_nzm1
public :: plan_forw_big_nzm1,   plan_back_big_nzm1
public :: plan_forw_small_one,  plan_back_small_one
public :: nz_full

! Batch sizes used by convec / derivs / press (will grow as we port more)
integer :: nz_full = 0    ! = nz - lbz + 1  (covers slabs lbz..nz, Phase A)

! Plan handles
integer :: plan_forw_small_full = -1, plan_back_small_full = -1
integer :: plan_forw_big_full   = -1, plan_back_big_full   = -1
integer :: plan_forw_small_nz   = -1, plan_back_small_nz   = -1
integer :: plan_forw_big_nz     = -1, plan_back_big_nz     = -1
integer :: plan_forw_small_nzm1 = -1, plan_back_small_nzm1 = -1
integer :: plan_forw_big_nzm1   = -1, plan_back_big_nzm1   = -1
! Single-slab plans (batch=1) - used for the rbottomw/rtopw and the boundary
! slabs in press_gpu (Phase A's nz-1 and Phase E's slab 0/nz inverse FFTs).
integer :: plan_forw_small_one  = -1, plan_back_small_one  = -1

contains

!*******************************************************************************
subroutine init_fft_gpu()
!*******************************************************************************
! Build all batched cuFFT plans. Must be called AFTER param.f90's
! ld/lh/nx/ny/ld_big/ny2 have been initialized (i.e. after read_input_conf)
! and AFTER the CUDA context exists (i.e. after the first OpenACC region).
!*******************************************************************************
implicit none
integer :: lh, lh_big_local
integer(acc_handle_kind) :: stream

nz_full = nz - lbz + 1
lh           = nx /2 + 1
lh_big_local = nx2/2 + 1

! Create the 12 plans
call make_plan(plan_forw_small_full, plan_back_small_full,                     &
               nx, ny, ld, lh, nz_full)
call make_plan(plan_forw_big_full,   plan_back_big_full,                       &
               nx2, ny2, ld_big, lh_big_local, nz_full)

call make_plan(plan_forw_small_nz, plan_back_small_nz,                         &
               nx, ny, ld, lh, nz)
call make_plan(plan_forw_big_nz,   plan_back_big_nz,                           &
               nx2, ny2, ld_big, lh_big_local, nz)

call make_plan(plan_forw_small_nzm1, plan_back_small_nzm1,                     &
               nx, ny, ld, lh, nz - 1)
call make_plan(plan_forw_big_nzm1,   plan_back_big_nzm1,                       &
               nx2, ny2, ld_big, lh_big_local, nz - 1)

call make_plan(plan_forw_small_one,  plan_back_small_one,                      &
               nx, ny, ld, lh, 1)

! Bind cuFFT calls to OpenACC async queue 1. Every !$acc parallel loop in
! *_gpu.f90 uses `async(1)` to match this queue. We use a numbered queue
! (vs acc_async_noval) so that `!$acc wait(1)` maps to cudaStreamSynchronize
! on this specific stream rather than `!$acc wait` -> cuCtxSynchronize
! (which is roughly 10x more expensive in NVHPC 25.5).
stream = acc_get_cuda_stream(1)
call set_stream(plan_forw_small_full, stream)
call set_stream(plan_back_small_full, stream)
call set_stream(plan_forw_big_full,   stream)
call set_stream(plan_back_big_full,   stream)
call set_stream(plan_forw_small_nz,   stream)
call set_stream(plan_back_small_nz,   stream)
call set_stream(plan_forw_big_nz,     stream)
call set_stream(plan_back_big_nz,     stream)
call set_stream(plan_forw_small_nzm1, stream)
call set_stream(plan_back_small_nzm1, stream)
call set_stream(plan_forw_big_nzm1,   stream)
call set_stream(plan_back_big_nzm1,   stream)
call set_stream(plan_forw_small_one,  stream)
call set_stream(plan_back_small_one,  stream)

end subroutine init_fft_gpu

!*******************************************************************************
subroutine make_plan(plan_forw, plan_back, n1, n2, ld_real, ld_cmplx, batch)
!*******************************************************************************
! Build a forward (D2Z) and backward (Z2D) batched 2D plan. The real-space
! layout is (ld_real, n2) per slab with batch slabs; the complex-space layout
! is (ld_cmplx, n2) per slab with the same batch.
!*******************************************************************************
implicit none
integer, intent(out) :: plan_forw, plan_back
integer, intent(in)  :: n1, n2, ld_real, ld_cmplx, batch
integer :: istat
integer :: n(2), inembed(2), onembed(2)

! C-order: outer dim = n2, inner dim = n1 (the dim being halved for r2c)
n(1)       = n2
n(2)       = n1
inembed(1) = n2
inembed(2) = ld_real
onembed(1) = n2
onembed(2) = ld_cmplx

istat = cufftPlanMany(plan_forw, 2, n,                                         &
                      inembed, 1, ld_real*n2,                                  &
                      onembed, 1, ld_cmplx*n2,                                 &
                      CUFFT_D2Z, batch)
call check_cufft(istat, 'cufftPlanMany forw')

istat = cufftPlanMany(plan_back, 2, n,                                         &
                      onembed, 1, ld_cmplx*n2,                                 &
                      inembed, 1, ld_real*n2,                                  &
                      CUFFT_Z2D, batch)
call check_cufft(istat, 'cufftPlanMany back')

end subroutine make_plan

!*******************************************************************************
subroutine set_stream(plan, stream)
!*******************************************************************************
implicit none
integer, intent(in) :: plan
integer(acc_handle_kind), intent(in) :: stream
integer :: istat
istat = cufftSetStream(plan, stream)
call check_cufft(istat, 'cufftSetStream')
end subroutine set_stream

!*******************************************************************************
subroutine finalize_fft_gpu()
!*******************************************************************************
implicit none
call destroy_plan(plan_forw_small_full); call destroy_plan(plan_back_small_full)
call destroy_plan(plan_forw_big_full);   call destroy_plan(plan_back_big_full)
call destroy_plan(plan_forw_small_nz);   call destroy_plan(plan_back_small_nz)
call destroy_plan(plan_forw_big_nz);     call destroy_plan(plan_back_big_nz)
call destroy_plan(plan_forw_small_nzm1); call destroy_plan(plan_back_small_nzm1)
call destroy_plan(plan_forw_big_nzm1);   call destroy_plan(plan_back_big_nzm1)
call destroy_plan(plan_forw_small_one);  call destroy_plan(plan_back_small_one)
end subroutine finalize_fft_gpu

subroutine destroy_plan(plan)
integer, intent(inout) :: plan
integer :: istat
if (plan >= 0) then
    istat = cufftDestroy(plan)
    plan = -1
end if
end subroutine destroy_plan

!*******************************************************************************
! Thin exec wrappers. Caller passes the plan handle and an array pointer; cuFFT
! is invoked in-place. Caller is responsible for choosing the plan whose batch
! count matches the array slab count.
!
! The `data` argument is declared as a 1D array (no shape contract) so callers
! can pass any (ld, ny, batch) slice - e.g. `RHSx(:,:,lbz)` or `cc_big(:,:,1)`.
!*******************************************************************************
subroutine fft_gpu_exec_d2z(plan, data)
implicit none
integer, intent(in) :: plan
real(rprec), intent(inout) :: data(*)
integer :: istat
!$acc host_data use_device(data)
istat = cufftExecD2Z(plan, data, data)
!$acc end host_data
call check_cufft(istat, 'cufftExecD2Z')
end subroutine fft_gpu_exec_d2z

subroutine fft_gpu_exec_z2d(plan, data)
implicit none
integer, intent(in) :: plan
real(rprec), intent(inout) :: data(*)
integer :: istat
!$acc host_data use_device(data)
istat = cufftExecZ2D(plan, data, data)
!$acc end host_data
call check_cufft(istat, 'cufftExecZ2D')
end subroutine fft_gpu_exec_z2d

!*******************************************************************************
subroutine check_cufft(istat, where)
!*******************************************************************************
implicit none
integer, intent(in) :: istat
character(len=*), intent(in) :: where
if (istat /= CUFFT_SUCCESS) then
    write(*,'(3a,i0)') 'cuFFT error at ', trim(where), ': code=', istat
    error stop
end if
end subroutine check_cufft

#endif
end module fft_gpu
