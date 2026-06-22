!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!  GPU port: see docs/gpu_module_contracts.md for current ownership context.
!!

!*******************************************************************************
module convec_gpu_m
!*******************************************************************************
! GPU implementation of the convective term c = -(u cross omega) using
! batched cuFFT plans plus OpenACC kernels for const-multiply, padd,
! vorticity, and cross-product arithmetic.
!
! Drop-in replacement for the CPU `convec()` subroutine in convec.f90. Selected
! at compile time via PPCONVEC_GPU (set by CMake when USE_LES_GPU=ON). Dispatch
! happens in main.f90.
!
! Algorithm and array semantics are kept identical to the CPU version
! (5 phases A-E, same boundary handling, same 3/2 dealiasing). The only
! differences from the CPU path are:
!   1. FFTs are batched along z (one cuFFT launch per direction-per-array,
!      replacing one FFTW call per (jz, array)).
!   2. The 7 padded big arrays (u_big, v_big, w_big, vort1/2/3_big, cc_big)
!      are persistent device-resident allocatables, never round-tripped to
!      host between phases.
!   3. The arithmetic kernels are `!$acc parallel loop collapse(3)` over
!      (jz, jy, jx). Wall BCs are issued as small separate kernels guarded by
!      host-side `if(coord==...)` so we don't get branch divergence inside the
!      hot loops.
!
! Lifecycle:
!   convec_gpu_init      - called once on first use; allocates the 7 device
!                          arrays and registers them with OpenACC.
!   convec_gpu           - replaces convec(); called every time step.
!   convec_gpu_finalize  - called from finalize.f90 via main exit path.
!
! Data movement contract: u, v, w, dudy/dudz/dvdx/dvdz/dwdx/dwdy and RHSx/y/z
! are all device-resident (sim_param `!$acc declare create`). convec reads the
! velocity gradients and writes RHSx/y/z entirely on the device; no PCIe. The
! `!$acc data present(...)` region just asserts residency; RHSx/y/z are read on
! device by main.f90's RHS = -RHS - divt kernel.
!*******************************************************************************
#ifdef PPCONVEC_GPU
use types, only : rprec
use param, only : ld, ny, nz, lbz, ld_big, ny2, nx, nx2,                       &
                  nproc, coord, sgs, lbc_mom, ubc_mom, BOGUS
use sim_param, only : u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy,             &
                      RHSx, RHSy, RHSz
use fft_gpu, only : fft_gpu_exec_d2z, fft_gpu_exec_z2d,                        &
                    plan_forw_small_full, plan_back_big_full,                  &
                    plan_forw_small_nz,   plan_back_big_nz,                    &
                    plan_forw_big_nzm1,   plan_back_small_nzm1,                &
                    plan_forw_big_nz,     plan_back_small_nz
implicit none
save
private
public :: convec_gpu, convec_gpu_finalize, convec_gpu_big_available
public :: u_big_d, v_big_d, w_big_d

! Persistent device-resident padded big arrays (3/2 dealiasing grid).
! Allocated lazily on first call (see convec_gpu_init).
real(rprec), allocatable, dimension(:,:,:) :: cc_big_d
real(rprec), allocatable, dimension(:,:,:) :: u_big_d, v_big_d, w_big_d
real(rprec), allocatable, dimension(:,:,:) :: vort1_big_d, vort2_big_d, vort3_big_d
logical :: initialized = .false.

contains

!*******************************************************************************
subroutine convec_gpu_init()
!*******************************************************************************
! First-touch allocation + OpenACC enter-data for the persistent big arrays.
!*******************************************************************************
implicit none
if (initialized) return

! u_big, v_big, w_big carry a halo at jz=lbz (filled in Phase A loop).
allocate(u_big_d(ld_big, ny2, lbz:nz)); u_big_d = 0._rprec
allocate(v_big_d(ld_big, ny2, lbz:nz)); v_big_d = 0._rprec
allocate(w_big_d(ld_big, ny2, lbz:nz)); w_big_d = 0._rprec
! vort and cc don't carry a halo (used at jz = 1..nz only).
allocate(vort1_big_d(ld_big, ny2, nz)); vort1_big_d = 0._rprec
allocate(vort2_big_d(ld_big, ny2, nz)); vort2_big_d = 0._rprec
allocate(vort3_big_d(ld_big, ny2, nz)); vort3_big_d = 0._rprec
allocate(cc_big_d(ld_big, ny2, nz));     cc_big_d   = 0._rprec

!$acc enter data copyin(u_big_d, v_big_d, w_big_d,                             &
!$acc                   vort1_big_d, vort2_big_d, vort3_big_d, cc_big_d)

initialized = .true.
end subroutine convec_gpu_init

!*******************************************************************************
logical function convec_gpu_big_available()
!*******************************************************************************
! True after convec_gpu_init has allocated and entered the persistent padded
! velocity fields.  Scalar transport can then reuse u_big_d/v_big_d/w_big_d
! instead of redoing the same small-to-big transforms.
!*******************************************************************************
implicit none

convec_gpu_big_available = initialized

end function convec_gpu_big_available

!*******************************************************************************
subroutine convec_gpu_finalize()
!*******************************************************************************
implicit none
if (.not. initialized) return
!$acc exit data delete(u_big_d, v_big_d, w_big_d,                              &
!$acc                  vort1_big_d, vort2_big_d, vort3_big_d, cc_big_d)
deallocate(u_big_d, v_big_d, w_big_d,                                          &
           vort1_big_d, vort2_big_d, vort3_big_d, cc_big_d)
initialized = .false.
end subroutine convec_gpu_finalize

!*******************************************************************************
subroutine convec_gpu()
!*******************************************************************************
! GPU replacement for convec(). Computes RHSx, RHSy, RHSz = -(u cross omega) using
! 3/2-rule dealiasing.
!*******************************************************************************
implicit none

integer :: jz, jy, jx
integer :: ny_h, j_s, j_big_s
integer :: jzLo, jzHi, jz_min, jz_max
real(rprec) :: const, const_big
real(rprec) :: half

call convec_gpu_init()

const     = 1._rprec / real(nx*ny,   rprec)
const_big = 1._rprec / real(nx2*ny2, rprec)
half      = 0.5_rprec

ny_h    = ny/2
j_s     = ny_h + 2
j_big_s = ny2 - ny_h + 2

! sgs flag changes the BC vorticity index used at the boundary slabs.
if (sgs) then
    jzLo = 2
    jzHi = nz - 1
else
    jzLo = 1
    jzHi = nz - 1
end if

! All inputs/outputs are device-resident via sim_param's `!$acc declare
! create`: u, v, w + dudy/dudz/dvdx/dvdz/dwdx/dwdy are inputs; RHSx/y/z are
! outputs read on device by main.f90's RHS = -RHS - divt kernel. They are listed
! `present` (not copyout); the declare-create device copy is
! authoritative; the host copy is never read on the steady-state path. Every
! RHS element is overwritten by the const-multiply kernel before any read.
!$acc data present(RHSx, RHSy, RHSz,                                           &
!$acc              u, v, w,                                                    &
!$acc              dudy, dudz, dvdx, dvdz, dwdx, dwdy,                         &
!$acc              u_big_d, v_big_d, w_big_d,                                  &
!$acc              vort1_big_d, vort2_big_d, vort3_big_d, cc_big_d)

! ============================================================================
! PHASE A: u, v, w -> small-grid spectral -> pad to big -> big physical
! ============================================================================

! A.1: const-multiply  RHSx,y,z = const * u,v,w  (slabs lbz..nz)
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz, nz
    do jy = 1, ny
        do jx = 1, ld
            RHSx(jx, jy, jz) = const * u(jx, jy, jz)
            RHSy(jx, jy, jz) = const * v(jx, jy, jz)
            RHSz(jx, jy, jz) = const * w(jx, jy, jz)
        end do
    end do
end do

! A.2: forward FFT (small grid, batch=nz_full) on each of RHSx, RHSy, RHSz
call fft_gpu_exec_d2z(plan_forw_small_full, RHSx(1,1,lbz))
call fft_gpu_exec_d2z(plan_forw_small_full, RHSy(1,1,lbz))
call fft_gpu_exec_d2z(plan_forw_small_full, RHSz(1,1,lbz))

! A.3: pad small spectral -> big spectral. Three steps to mirror padd():
!   (1) zero the big array entirely
!   (2) copy lower-frequency block [1..nx, 1..ny_h]
!   (3) copy upper-frequency block [1..nx, j_s..ny] -> [1..nx, j_big_s..ny2]

! Step (1): zero u_big, v_big, w_big over lbz..nz, all (jy, jx)
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz, nz
    do jy = 1, ny2
        do jx = 1, ld_big
            u_big_d(jx, jy, jz) = 0._rprec
            v_big_d(jx, jy, jz) = 0._rprec
            w_big_d(jx, jy, jz) = 0._rprec
        end do
    end do
end do

! Step (2): lower freq block
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz, nz
    do jy = 1, ny_h
        do jx = 1, nx
            u_big_d(jx, jy, jz) = RHSx(jx, jy, jz)
            v_big_d(jx, jy, jz) = RHSy(jx, jy, jz)
            w_big_d(jx, jy, jz) = RHSz(jx, jy, jz)
        end do
    end do
end do

! Step (3): upper freq block (wrapped to high end of ny2)
!$acc parallel loop collapse(3) default(present) async(1)
do jz = lbz, nz
    do jy = j_s, ny    !  ny_h - 1 iterations  (jy = ny_h+2 .. ny)
        do jx = 1, nx
            u_big_d(jx, j_big_s + (jy - j_s), jz) = RHSx(jx, jy, jz)
            v_big_d(jx, j_big_s + (jy - j_s), jz) = RHSy(jx, jy, jz)
            w_big_d(jx, j_big_s + (jy - j_s), jz) = RHSz(jx, jy, jz)
        end do
    end do
end do

! A.4: backward FFT (big grid, batch=nz_full) -> physical big velocity
call fft_gpu_exec_z2d(plan_back_big_full, u_big_d(1,1,lbz))
call fft_gpu_exec_z2d(plan_back_big_full, v_big_d(1,1,lbz))
call fft_gpu_exec_z2d(plan_back_big_full, w_big_d(1,1,lbz))

! ============================================================================
! PHASE B: vorticity (omega_x, omega_y, omega_z) -> small spectral -> pad
! to big -> big physical
! ============================================================================

! B.1: General interior formula for ALL slabs jz=1..nz. Wall BCs may overwrite.
!      RHSx = const * (dwdy - dvdz)
!      RHSy = const * (dudz - dwdx)
!      RHSz = const * (dvdx - dudy)
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
    do jy = 1, ny
        do jx = 1, ld
            RHSx(jx, jy, jz) = const * (dwdy(jx, jy, jz) - dvdz(jx, jy, jz))
            RHSy(jx, jy, jz) = const * (dudz(jx, jy, jz) - dwdx(jx, jy, jz))
            RHSz(jx, jy, jz) = const * (dvdx(jx, jy, jz) - dudy(jx, jy, jz))
        end do
    end do
end do

! B.2: bottom-wall BC for vort1, vort2 (overwrites general formula at jz=1)
if (coord == 0) then
    select case (lbc_mom)
    case (0)    ! stress-free
        !$acc parallel loop collapse(2) default(present) async(1)
        do jy = 1, ny
            do jx = 1, ld
                RHSx(jx, jy, 1) = 0._rprec
                RHSy(jx, jy, 1) = 0._rprec
            end do
        end do
    case (1:)   ! wall (no-slip / wall-model)
        !$acc parallel loop collapse(2) default(present) async(1)
        do jy = 1, ny
            do jx = 1, ld
                RHSx(jx, jy, 1) = const * (                                    &
                    half * (dwdy(jx, jy, 1) + dwdy(jx, jy, 2))                 &
                    - dvdz(jx, jy, 1) )
                RHSy(jx, jy, 1) = const * (                                    &
                    dudz(jx, jy, 1)                                            &
                    - half * (dwdx(jx, jy, 1) + dwdx(jx, jy, 2)) )
            end do
        end do
    end select
end if

! B.3: top-wall BC for vort1, vort2 (overwrites general formula at jz=nz).
!      Note: when ubc_mom == 0 the CPU path falls through to general formula
!      (the case-0 zero is overwritten by step 3 of the original loop). We
!      mirror that behaviour by NOT issuing a top BC kernel for ubc_mom == 0.
if (coord == nproc-1 .and. ubc_mom > 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            RHSx(jx, jy, nz) = const * (                                       &
                half * (dwdy(jx, jy, nz-1) + dwdy(jx, jy, nz))                 &
                - dvdz(jx, jy, nz-1) )
            RHSy(jx, jy, nz) = const * (                                       &
                dudz(jx, jy, nz-1)                                             &
                - half * (dwdx(jx, jy, nz-1) + dwdx(jx, jy, nz)) )
        end do
    end do
end if

! B.4: forward FFT (small, batch=nz) on RHSx, RHSy, RHSz
call fft_gpu_exec_d2z(plan_forw_small_nz, RHSx(1,1,1))
call fft_gpu_exec_d2z(plan_forw_small_nz, RHSy(1,1,1))
call fft_gpu_exec_d2z(plan_forw_small_nz, RHSz(1,1,1))

! B.5: pad small -> big (vort1, vort2, vort3). Same 3-step pattern as A.3 but
!      with destination arrays vort1/2/3_big_d and source arrays RHSx/y/z.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
    do jy = 1, ny2
        do jx = 1, ld_big
            vort1_big_d(jx, jy, jz) = 0._rprec
            vort2_big_d(jx, jy, jz) = 0._rprec
            vort3_big_d(jx, jy, jz) = 0._rprec
        end do
    end do
end do

!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
    do jy = 1, ny_h
        do jx = 1, nx
            vort1_big_d(jx, jy, jz) = RHSx(jx, jy, jz)
            vort2_big_d(jx, jy, jz) = RHSy(jx, jy, jz)
            vort3_big_d(jx, jy, jz) = RHSz(jx, jy, jz)
        end do
    end do
end do

!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
    do jy = j_s, ny
        do jx = 1, nx
            vort1_big_d(jx, j_big_s + (jy - j_s), jz) = RHSx(jx, jy, jz)
            vort2_big_d(jx, j_big_s + (jy - j_s), jz) = RHSy(jx, jy, jz)
            vort3_big_d(jx, j_big_s + (jy - j_s), jz) = RHSz(jx, jy, jz)
        end do
    end do
end do

! B.6: backward FFT (big, batch=nz) -> physical big vorticity
call fft_gpu_exec_z2d(plan_back_big_nz, vort1_big_d(1,1,1))
call fft_gpu_exec_z2d(plan_back_big_nz, vort2_big_d(1,1,1))
call fft_gpu_exec_z2d(plan_back_big_nz, vort3_big_d(1,1,1))

! ============================================================================
! PHASE C: RHSx convective product:
!          cc_big = -v*omega_z + 0.5*(w*omega_y stencil)
! ============================================================================
! Per-coord interior limits (boundaries handled by separate kernels):
jz_min = 1; if (coord == 0)        jz_min = 2
jz_max = nz-1; if (coord == nproc-1) jz_max = nz-2

! C.1: interior cc_big for jz_min..jz_max
!$acc parallel loop collapse(3) default(present) async(1)
do jz = jz_min, jz_max
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, jz) = const_big * (                               &
                v_big_d(jx, jy, jz) * (-vort3_big_d(jx, jy, jz))               &
                + half * ( w_big_d(jx, jy, jz+1) * vort2_big_d(jx, jy, jz+1)   &
                         + w_big_d(jx, jy, jz)   * vort2_big_d(jx, jy, jz)) )
        end do
    end do
end do

! C.2: bottom wall slab (jz=1) for coord==0
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, 1) = const_big * (                                &
                v_big_d(jx, jy, 1) * (-vort3_big_d(jx, jy, 1))                 &
                + half * w_big_d(jx, jy, 2) * vort2_big_d(jx, jy, jzLo) )
        end do
    end do
end if

! C.3: top channel slab (jz=nz-1) for coord==nproc-1
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, nz-1) = const_big * (                             &
                v_big_d(jx, jy, nz-1) * (-vort3_big_d(jx, jy, nz-1))           &
                + half * w_big_d(jx, jy, nz-1) * vort2_big_d(jx, jy, jzHi) )
        end do
    end do
end if

! C.4: forward FFT (big, batch=nz-1) on cc_big slabs 1..nz-1
call fft_gpu_exec_d2z(plan_forw_big_nzm1, cc_big_d(1,1,1))

! C.5: unpadd big -> small -> RHSx[1..nz-1]
call unpadd_gpu(RHSx, cc_big_d, nz-1)

! C.6: backward FFT (small, batch=nz-1) on RHSx
call fft_gpu_exec_z2d(plan_back_small_nzm1, RHSx(1,1,1))

! ============================================================================
! PHASE D: RHSy convective product:
!          cc_big = u*omega_z - 0.5*(w*omega_x stencil)
! ============================================================================
jz_min = 1; if (coord == 0)        jz_min = 2
jz_max = nz-1; if (coord == nproc-1) jz_max = nz-2

! D.1: interior
!$acc parallel loop collapse(3) default(present) async(1)
do jz = jz_min, jz_max
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, jz) = const_big * (                               &
                u_big_d(jx, jy, jz) * vort3_big_d(jx, jy, jz)                  &
                + half * ( w_big_d(jx, jy, jz+1) * (-vort1_big_d(jx, jy, jz+1))&
                         + w_big_d(jx, jy, jz)   * (-vort1_big_d(jx, jy, jz))) )
        end do
    end do
end do

! D.2: bottom wall slab
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, 1) = const_big * (                                &
                u_big_d(jx, jy, 1) * vort3_big_d(jx, jy, 1)                    &
                + half * w_big_d(jx, jy, 2) * (-vort1_big_d(jx, jy, jzLo)) )
        end do
    end do
end if

! D.3: top channel slab
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, nz-1) = const_big * (                             &
                u_big_d(jx, jy, nz-1) * vort3_big_d(jx, jy, nz-1)              &
                + half * w_big_d(jx, jy, nz-1) * (-vort1_big_d(jx, jy, jzHi)) )
        end do
    end do
end if

! D.4: forward FFT, unpadd, backward FFT -> RHSy[1..nz-1]
call fft_gpu_exec_d2z(plan_forw_big_nzm1, cc_big_d(1,1,1))
call unpadd_gpu(RHSy, cc_big_d, nz-1)
call fft_gpu_exec_z2d(plan_back_small_nzm1, RHSy(1,1,1))

! ============================================================================
! PHASE E: RHSz convective product:
!          cc_big = 0.5*((u+u_below)*(-omega_y) + (v+v_below)*omega_x)
! ============================================================================
! Note: this stencil reads u_big(:,:,jz-1). At interior ranks, jz=1 reads
! u_big(:,:,0) which is the halo slab populated in Phase A. At coord==0 we
! force jz=1 to zero (no convective accel of w at the wall) so we never
! actually read u_big(:,:,0) on rank 0.
jz_min = 1; if (coord == 0)        jz_min = 2
jz_max = nz-1   ! same upper limit on all ranks for E (channel)

! E.1: interior
!$acc parallel loop collapse(3) default(present) async(1)
do jz = jz_min, jz_max
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, jz) = const_big * half * (                        &
                (u_big_d(jx, jy, jz) + u_big_d(jx, jy, jz-1))                  &
                    * (-vort2_big_d(jx, jy, jz))                               &
                + (v_big_d(jx, jy, jz) + v_big_d(jx, jy, jz-1))                &
                    * vort1_big_d(jx, jy, jz) )
        end do
    end do
end do

! E.2: bottom wall slab - no convective accel of w at wall
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, 1) = 0._rprec
        end do
    end do
end if

! E.3: top channel slab - no convective accel of w at top (cc_big(:,:,nz)=0)
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny2
        do jx = 1, ld_big
            cc_big_d(jx, jy, nz) = 0._rprec
        end do
    end do
end if

! E.4: forward FFT, unpadd, backward FFT -> RHSz[1..nz]
!      This phase covers slabs 1..nz (one extra over C/D).
call fft_gpu_exec_d2z(plan_forw_big_nz, cc_big_d(1,1,1))
call unpadd_gpu(RHSz, cc_big_d, nz)
call fft_gpu_exec_z2d(plan_back_small_nz, RHSz(1,1,1))

! ============================================================================
! Trailing safety-mode boundary scrubs (mirror CPU PPSAFETYMODE block)
! ============================================================================
#ifdef PPMPI
#ifdef PPSAFETYMODE
!$acc parallel loop collapse(2) default(present) async(1)
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
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        RHSx(jx, jy, nz) = BOGUS
        RHSy(jx, jy, nz) = BOGUS
    end do
end do
if (coord < nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            RHSz(jx, jy, nz) = BOGUS
        end do
    end do
end if
#endif

! NOTE (overlap): no terminal `wait(1)` - the data region above is
! present-only (no create/copy), so nothing is freed at `end data` and the
! async(1) kernels keep draining after we return. The convection backlog then
! overlaps the ATM phase-1 host work (blade-force model + MPI gather) that
! main.f90 runs right after this call; queue-1 FIFO order keeps every
! downstream consumer (RHS accumulation, velocity update, pressure) correct,
! and the pressure solver's pre-halo wait(1) is the absorber.
!$acc end data    ! RHSx/y/z stay device-resident (declare create); no copyout

end subroutine convec_gpu

!*******************************************************************************
subroutine unpadd_gpu(cc, cc_big, nz_use)
!*******************************************************************************
! Inverse of padd: copy big spectral cc_big[1..nz_use] back into the (ld, ny)
! small spectral cc array, zero the oddballs. Mirrors the CPU `unpadd()` in
! fft.f90 line by line, just collapsed over jz.
!*******************************************************************************
implicit none
real(rprec), intent(inout) :: cc(ld, ny, lbz:nz)
real(rprec), intent(in)    :: cc_big(ld_big, ny2, nz)
integer,     intent(in)    :: nz_use

integer :: jz, jy, jx
integer :: ny_h, j_s, j_big_s

ny_h    = ny/2
j_s     = ny_h + 2
j_big_s = ny2 - ny_h + 2

! Lower freq block: cc[1..nx, 1..ny_h] = cc_big[1..nx, 1..ny_h]
!$acc parallel loop collapse(3) present(cc, cc_big) async(1)
do jz = 1, nz_use
    do jy = 1, ny_h
        do jx = 1, nx
            cc(jx, jy, jz) = cc_big(jx, jy, jz)
        end do
    end do
end do

! Oddballs: cc[ld-1..ld, :] = 0  and  cc[:, ny_h+1] = 0
!$acc parallel loop collapse(2) present(cc) async(1)
do jz = 1, nz_use
    do jy = 1, ny
        cc(ld-1, jy, jz) = 0._rprec
        cc(ld,   jy, jz) = 0._rprec
    end do
end do
!$acc parallel loop collapse(2) present(cc) async(1)
do jz = 1, nz_use
    do jx = 1, ld
        cc(jx, ny_h+1, jz) = 0._rprec
    end do
end do

! Upper freq block: cc[1..nx, j_s..ny] = cc_big[1..nx, j_big_s..ny2]
!$acc parallel loop collapse(3) present(cc, cc_big) async(1)
do jz = 1, nz_use
    do jy = j_s, ny
        do jx = 1, nx
            cc(jx, jy, jz) = cc_big(jx, j_big_s + (jy - j_s), jz)
        end do
    end do
end do

end subroutine unpadd_gpu

#endif
end module convec_gpu_m
