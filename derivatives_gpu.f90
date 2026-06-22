!!
!!  Copyright (C) 2010-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!  GPU port: see docs/gpu_module_contracts.md for current ownership context.
!!

!*******************************************************************************
module derivatives_gpu_m
!*******************************************************************************
! GPU implementation of the three derivative routines called per time step from
! main.f90:
!
!   filt_da_gpu  - pseudospectral df/dx and df/dy of f, plus de-aliased f.
!                  CPU equivalent: filt_da() in derivatives.f90.
!   ddz_uv_gpu   - 2nd-order finite-difference d/dz, uv-grid -> w-grid.
!                  CPU equivalent: ddz_uv().
!   ddz_w_gpu    - 2nd-order finite-difference d/dz, w-grid -> uv-grid.
!                  CPU equivalent: ddz_w().
!
! NOT ported (still CPU):
!   ddx, ddy, ddxy - only used by divstress_uv/w which run during the SGS
!                    block. Will be done when sgs_gpu lands.
!
! Data movement contract (post sim_param refactor for u, v, w):
!   filt_da_gpu:   present(f, dfdx, dfdy)         no PCIe
!   ddz_uv_gpu:    present(f, dfdz)               no PCIe
!   ddz_w_gpu:     present(f, dfdz)               no PCIe
!
! All inputs and outputs are device-resident via sim_param's `!$acc declare
! create(u, v, w, dudx..dwdz)`. The caller (main.f90) issues a single
! `!$acc update self(u, v, w, dudx..dwdz)` after the derivatives + convec
! block before host code (wallstress, sgs_stag) needs them.
!
! filt_da_gpu modifies f in-place (zeros the spectral oddballs after the
! forward FFT, then inverse-transforms). Because f is device-resident the
! de-aliased version stays on the device for free; the corresponding host
! copy is refreshed via the next `update self`.
!
! kx, ky live permanently on the device (push happens once in
! fft.f90::init_wavenumber).
!
! cuFFT plans are batched over nz_full = nz - lbz + 1 slabs, identical to
! the plans used by Phase A of convec_gpu.
!*******************************************************************************
#ifdef PPDERIVS_GPU
use types, only : rprec
use param, only : ld, nx, ny, nz, lbz, dz, BOGUS
#ifdef PPSAFETYMODE
use param, only : nproc, coord
#endif
use fft, only : kx, ky
use fft_gpu, only : fft_gpu_exec_d2z, fft_gpu_exec_z2d,                        &
                    plan_forw_small_full, plan_back_small_full
implicit none
save
private
public :: filt_da_gpu, ddz_uv_gpu, ddz_w_gpu
public :: ddx_gpu, ddy_gpu, ddxy_gpu

contains

!*******************************************************************************
subroutine filt_da_gpu(f, dfdx, dfdy, lbz_arg)
!*******************************************************************************
! GPU equivalent of filt_da. Computes df/dx and df/dy in spectral space and
! also writes back the de-aliased f (oddballs zeroed). All three transforms
! are batched over the full lbz..nz slab range.
!*******************************************************************************
implicit none
integer, intent(in) :: lbz_arg
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: f
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: dfdx, dfdy

integer  :: jz, jy, i, ir, ii
integer  :: lh, ny_h
real(rprec) :: const, ar, ai, kxv, kyv

const = 1._rprec / real(nx*ny, rprec)
lh    = nx/2 + 1
ny_h  = ny/2

! All inputs/outputs are device-resident via sim_param's `!$acc declare
! create`:
!   f       - u, v, or w. Modified in-place (de-aliased) on device; the
!             host copy stays stale until the caller's next update_self.
!   dfdx, dfdy - derivative outputs.
!   kx, ky  - wavenumber arrays (init_wavenumber did the one-time push).
!$acc data present(f, dfdx, dfdy, kx, ky)

! ----- 1. const-multiply f *= 1/(nx*ny) -----
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg, nz
    do jy = 1, ny
        do i = 1, ld
            f(i, jy, jz) = const * f(i, jy, jz)
        end do
    end do
end do

! ----- 2. forward FFT in-place (batched, all nz_full slabs) -----
call fft_gpu_exec_d2z(plan_forw_small_full, f(1, 1, lbz_arg))

! ----- 3. zero oddballs in f (the (ld-1, ld) Nyquist row and the ny_h+1 col) -----
!$acc parallel loop collapse(2) default(present) async(1)
do jz = lbz_arg, nz
    do jy = 1, ny
        f(ld-1, jy, jz) = 0._rprec
        f(ld,   jy, jz) = 0._rprec
    end do
end do
!$acc parallel loop collapse(2) default(present) async(1)
do jz = lbz_arg, nz
    do i = 1, ld
        f(i, ny_h+1, jz) = 0._rprec
    end do
end do

! ----- 4. dfdx = f .MULI. kx ; dfdy = f .MULI. ky  (inlined) -----
! MULI rule (real array f viewed as interleaved complex):
!   for each complex index i (1..lh):
!     ir = 2i - 1       ! real part of f
!     ii = 2i           ! imag part of f
!     b(ir) = -f(ii) * k(i)
!     b(ii) =  f(ir) * k(i)
! Both dfdx and dfdy are computed in one kernel, sharing the loaded f values.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg, nz
    do jy = 1, ny
        do i = 1, lh
            ir  = 2*i - 1
            ii  = 2*i
            ar  = f(ir, jy, jz)
            ai  = f(ii, jy, jz)
            kxv = kx(i, jy)
            kyv = ky(i, jy)
            dfdx(ir, jy, jz) = -ai * kxv
            dfdx(ii, jy, jz) =  ar * kxv
            dfdy(ir, jy, jz) = -ai * kyv
            dfdy(ii, jy, jz) =  ar * kyv
        end do
    end do
end do

! ----- 5. Inverse FFT all three arrays back to physical -----
call fft_gpu_exec_z2d(plan_back_small_full, f(1,    1, lbz_arg))
call fft_gpu_exec_z2d(plan_back_small_full, dfdx(1, 1, lbz_arg))
call fft_gpu_exec_z2d(plan_back_small_full, dfdy(1, 1, lbz_arg))

! No terminal wait: all work is queued on stream 1 and all arrays remain
! device-resident. Callers already synchronize before host-visible operations.
!$acc end data    ! no PCIe - all arrays device-resident

end subroutine filt_da_gpu

!*******************************************************************************
subroutine ddz_uv_gpu(f, dfdz, lbz_arg)
!*******************************************************************************
! GPU equivalent of ddz_uv. f is on uv grid, dfdz is on w grid.
! Pure stencil: dfdz(jz) = (f(jz) - f(jz-1)) / dz
!*******************************************************************************
implicit none
integer, intent(in) :: lbz_arg
real(rprec), dimension(:,:,lbz_arg:), intent(in)    :: f
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: dfdz

integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec / dz

! Both f (u/v) and dfdz (dudz/dvdz) are device-resident via sim_param.
!$acc data present(f, dfdz)

#if defined(PPMPI) && defined(PPSAFETYMODE)
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, nx
        dfdz(jx, jy, 0) = BOGUS
    end do
end do
#endif

!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg+1, nz
    do jy = 1, ny
        do jx = 1, nx
            dfdz(jx, jy, jz) = const * (f(jx, jy, jz) - f(jx, jy, jz-1))
        end do
    end do
end do

#ifdef PPSAFETYMODE
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            dfdz(jx, jy, 1) = BOGUS
        end do
    end do
end if
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            dfdz(jx, jy, nz) = BOGUS
        end do
    end do
end if
#endif

! No terminal wait: queue-1 ordering preserves dependencies for downstream GPU
! kernels, and callers synchronize before host-visible operations.
!$acc end data

end subroutine ddz_uv_gpu

!*******************************************************************************
subroutine ddz_w_gpu(f, dfdz, lbz_arg)
!*******************************************************************************
! GPU equivalent of ddz_w. f is on w grid, dfdz is on uv grid.
! Pure stencil: dfdz(jz) = (f(jz+1) - f(jz)) / dz
!*******************************************************************************
implicit none
integer, intent(in) :: lbz_arg
real(rprec), dimension(:,:,lbz_arg:), intent(in)    :: f
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: dfdz

integer :: jx, jy, jz
real(rprec) :: const

const = 1._rprec / dz

! Both f (w) and dfdz (dwdz) are device-resident via sim_param.
!$acc data present(f, dfdz)

!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg, nz-1
    do jy = 1, ny
        do jx = 1, nx
            dfdz(jx, jy, jz) = const * (f(jx, jy, jz+1) - f(jx, jy, jz))
        end do
    end do
end do

#ifdef PPSAFETYMODE
#ifdef PPMPI
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            dfdz(jx, jy, lbz_arg) = BOGUS
        end do
    end do
end if
#endif
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, nx
        dfdz(jx, jy, nz) = BOGUS
    end do
end do
#endif

! No terminal wait: queue-1 ordering preserves dependencies for downstream GPU
! kernels, and callers synchronize before host-visible operations.
!$acc end data

end subroutine ddz_w_gpu

!*******************************************************************************
subroutine ddx_gpu(f, dfdx, lbz_arg)
!*******************************************************************************
! GPU equivalent of ddx (CPU derivatives.f90). Computes df/dx in spectral
! space across slabs lbz_arg..nz, batched in a single cuFFT call.
!
! Unlike filt_da_gpu, f is read-only - we copy const*f into dfdx, do the
! forward FFT in dfdx, zero the oddballs in dfdx, multiply by kx in place,
! then inverse FFT. f is never modified.
!
! Both f and dfdx must be device-resident in the caller's data scope.
!*******************************************************************************
implicit none
integer, intent(in) :: lbz_arg
real(rprec), dimension(:,:,lbz_arg:), intent(in)    :: f
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: dfdx

integer  :: jz, jy, i, ir, ii, lh, ny_h
real(rprec) :: const, ar, ai, kxv

const = 1._rprec / real(nx*ny, rprec)
lh    = nx/2 + 1
ny_h  = ny/2

!$acc data present(f, dfdx, kx)

! 1. dfdx = const * f
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg, nz
    do jy = 1, ny
        do i = 1, ld
            dfdx(i, jy, jz) = const * f(i, jy, jz)
        end do
    end do
end do

! 2. forward FFT (batched on slabs lbz_arg..nz)
call fft_gpu_exec_d2z(plan_forw_small_full, dfdx(1, 1, lbz_arg))

    ! 3. dfdx = dfdx .MULI. kx (in place using local registers).
    !    The old separate oddball-zero kernels are fused here: the x-Nyquist
    !    complex mode and y-Nyquist row become exactly zero derivative output.
    !$acc parallel loop collapse(3) default(present) private(ar, ai, kxv, ir, ii) async(1)
    do jz = lbz_arg, nz
        do jy = 1, ny
            do i = 1, lh
                ir  = 2*i - 1
                ii  = 2*i
                if (i == lh .or. jy == ny_h+1) then
                    dfdx(ir, jy, jz) = 0._rprec
                    dfdx(ii, jy, jz) = 0._rprec
                else
                    ar  = dfdx(ir, jy, jz)
                    ai  = dfdx(ii, jy, jz)
                    kxv = kx(i, jy)
                    dfdx(ir, jy, jz) = -ai * kxv
                    dfdx(ii, jy, jz) =  ar * kxv
                end if
            end do
        end do
    end do

    ! 4. inverse FFT (batched)
call fft_gpu_exec_z2d(plan_back_small_full, dfdx(1, 1, lbz_arg))

! No terminal wait: queue-1 ordering preserves dependencies for downstream GPU
! kernels, and callers synchronize before host-visible operations.
!$acc end data

end subroutine ddx_gpu

!*******************************************************************************
subroutine ddy_gpu(f, dfdy, lbz_arg)
!*******************************************************************************
! GPU equivalent of ddy. Same structure as ddx_gpu but multiplies by ky.
!*******************************************************************************
implicit none
integer, intent(in) :: lbz_arg
real(rprec), dimension(:,:,lbz_arg:), intent(in)    :: f
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: dfdy

integer  :: jz, jy, i, ir, ii, lh, ny_h
real(rprec) :: const, ar, ai, kyv

const = 1._rprec / real(nx*ny, rprec)
lh    = nx/2 + 1
ny_h  = ny/2

!$acc data present(f, dfdy, ky)

! 1. dfdy = const * f
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg, nz
    do jy = 1, ny
        do i = 1, ld
            dfdy(i, jy, jz) = const * f(i, jy, jz)
        end do
    end do
end do

! 2. forward FFT
call fft_gpu_exec_d2z(plan_forw_small_full, dfdy(1, 1, lbz_arg))

    ! 3. dfdy = dfdy .MULI. ky.  Fuse oddball zeroing into the multiply
    !    kernel to avoid two extra small launches per derivative call.
    !$acc parallel loop collapse(3) default(present) private(ar, ai, kyv, ir, ii) async(1)
    do jz = lbz_arg, nz
        do jy = 1, ny
            do i = 1, lh
                ir  = 2*i - 1
                ii  = 2*i
                if (i == lh .or. jy == ny_h+1) then
                    dfdy(ir, jy, jz) = 0._rprec
                    dfdy(ii, jy, jz) = 0._rprec
                else
                    ar  = dfdy(ir, jy, jz)
                    ai  = dfdy(ii, jy, jz)
                    kyv = ky(i, jy)
                    dfdy(ir, jy, jz) = -ai * kyv
                    dfdy(ii, jy, jz) =  ar * kyv
                end if
            end do
        end do
    end do

    ! 4. inverse FFT
call fft_gpu_exec_z2d(plan_back_small_full, dfdy(1, 1, lbz_arg))

! No terminal wait: queue-1 ordering preserves dependencies for downstream GPU
! kernels, and callers synchronize before host-visible operations.
!$acc end data

end subroutine ddy_gpu

!*******************************************************************************
subroutine ddxy_gpu(f, dfdx, dfdy, lbz_arg)
!*******************************************************************************
! GPU equivalent of ddxy: computes both df/dx and df/dy in one fused pass -
! single forward FFT of f, two .MULI. multiplies, two inverse FFTs.
!*******************************************************************************
implicit none
integer, intent(in) :: lbz_arg
real(rprec), dimension(:,:,lbz_arg:), intent(in)    :: f
real(rprec), dimension(:,:,lbz_arg:), intent(inout) :: dfdx, dfdy

integer  :: jz, jy, i, ir, ii, lh, ny_h
real(rprec) :: const, ar, ai, kxv, kyv

const = 1._rprec / real(nx*ny, rprec)
lh    = nx/2 + 1
ny_h  = ny/2

!$acc data present(f, dfdx, dfdy, kx, ky)

! 1. dfdx = const * f  (we use dfdx as the FFT scratch, then copy to dfdy
!    after multiply since each output diverges)
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz_arg, nz
    do jy = 1, ny
        do i = 1, ld
            dfdx(i, jy, jz) = const * f(i, jy, jz)
        end do
    end do
end do

! 2. forward FFT
call fft_gpu_exec_d2z(plan_forw_small_full, dfdx(1, 1, lbz_arg))

    ! 3. dfdy = dfdx .MULI. ky ; dfdx = dfdx .MULI. kx.
    !    Critical ordering: write dfdy first (it reads dfdx's pre-multiply state),
    !    then overwrite dfdx in-place.  Oddball zeroing is fused into this kernel.
    !$acc parallel loop collapse(3) default(present) private(ar, ai, kxv, kyv, ir, ii) async(1)
    do jz = lbz_arg, nz
        do jy = 1, ny
            do i = 1, lh
                ir  = 2*i - 1
                ii  = 2*i
                if (i == lh .or. jy == ny_h+1) then
                    dfdy(ir, jy, jz) = 0._rprec
                    dfdy(ii, jy, jz) = 0._rprec
                    dfdx(ir, jy, jz) = 0._rprec
                    dfdx(ii, jy, jz) = 0._rprec
                else
                    ar  = dfdx(ir, jy, jz)
                    ai  = dfdx(ii, jy, jz)
                    kxv = kx(i, jy)
                    kyv = ky(i, jy)
                    dfdy(ir, jy, jz) = -ai * kyv
                    dfdy(ii, jy, jz) =  ar * kyv
                    dfdx(ir, jy, jz) = -ai * kxv
                    dfdx(ii, jy, jz) =  ar * kxv
                end if
            end do
        end do
    end do

    ! 4. inverse FFTs
call fft_gpu_exec_z2d(plan_back_small_full, dfdx(1, 1, lbz_arg))
call fft_gpu_exec_z2d(plan_back_small_full, dfdy(1, 1, lbz_arg))

! No terminal wait: queue-1 ordering preserves dependencies for downstream GPU
! kernels, and callers synchronize before host-visible operations.
!$acc end data

end subroutine ddxy_gpu

#endif
end module derivatives_gpu_m
