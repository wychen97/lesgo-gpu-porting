!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!  GPU port: see docs/gpu_module_contracts.md for current ownership context.
!!

!*******************************************************************************
module press_gpu_m
!*******************************************************************************
! GPU implementation of press_stag_array.
!
! Drop-in replacement for the CPU `press_stag_array()` subroutine in
! press_stag_array.f90. Selected at compile time via PPPRESS_GPU (set by
! CMake when USE_LES_GPU=ON). Dispatch happens in main.f90.
!
! Algorithm and array semantics are identical to the CPU version:
!   Phase A - build spectral RHS (rH_x, rH_y, rH_z), boundary terms (rtopw,
!             rbottomw), zero oddballs.
!   Phase B - build the per-(jx,jy) tridiagonal system coefficients
!             (a, b, c, RHS_col) on the device.
!   Phase C - call `tridag_array_gpu` on device-resident coefficients and RHS.
!             tridag_gpu_m owns the MPI-chain policy, including the
!             GPU-aware path and host-staged fallback.
!   Phase D - DC-mode (jx=1, jy=1) integration. Tiny, stays on host.
!   Phase E - zero p oddballs; build dpdx, dpdy from spectral p; batched
!             inverse FFTs of (p, dpdx, dpdy).
!   Phase F - dpdz = (p(jz) - p(jz-1)) / dz finite difference.
!
! Ownership map:
!   - press_gpu_m owns pressure RHS/coefficient assembly and pressure-gradient
!     output assembly on persistent device scratch.
!   - tridag_gpu_m owns the distributed tridiagonal solve and decides whether
!     its rank-chain messages use GPU-aware MPI or the host-staged fallback.
!   - fft_gpu owns cuFFT plan lifecycle and spectral transforms.
!   - sim_param owns u/v/w/divtz inputs and p/dpdx/dpdy/dpdz outputs.
!
! Data movement contract:
!   copyin:  u, v, w, divtz
!   copyout: p, dpdx, dpdy, dpdz
!   Internal scratch arrays (rH_x_d, rH_y_d, rH_z_d, RHS_col_d, a_d, b_d,
!   c_d, rtopw_d, rbottomw_d) are persistent on device - allocated lazily on
!   first call and freed by press_gpu_finalize.
!
! Host-visible boundaries per call:
!   1. MPI halo exchange of rH_x/y/z slabs between Phase A and B.
!   2. The small DC-mode column loop and zero-mode MPI send/recv in Phase D.
!   3. Any host-staged fallback inside tridag_gpu_m when GPU-aware MPI is not
!      compiled or is disabled for A/B debugging.
!*******************************************************************************
#ifdef PPPRESS_GPU
use types, only : rprec
use param, only : ld, lh, ny, nz, lbz, nx, dz, dt, tadv1, BOGUS,               &
                  nproc, coord, jt_total,                                      &
                  domain_calc, domain_nstart, domain_nend, domain_nskip
#ifdef PPMPI
use param, only : MPI_RPREC, comm, ierr, status, up, down
use mpi
#endif
use sim_param, only : u, v, w, divtz, p, dpdx, dpdy, dpdz
use fft, only : kx, ky
use fft_gpu, only : fft_gpu_exec_d2z, fft_gpu_exec_z2d,                        &
                    plan_forw_small_nz,   plan_back_small_nz,                  &
                    plan_forw_small_nzm1, plan_back_small_nzm1,                &
                    plan_forw_small_one,  plan_back_small_one
use tridag_gpu_m, only : tridag_array_gpu
implicit none
save
private
public :: press_stag_array_gpu, press_gpu_finalize

! Persistent device-resident scratch.
real(rprec), allocatable, dimension(:,:,:) :: rH_x_d, rH_y_d, rH_z_d
real(rprec), allocatable, dimension(:,:,:) :: RHS_col_d
real(rprec), allocatable, dimension(:,:,:) :: a_d, b_d, c_d
real(rprec), allocatable, dimension(:,:)   :: rtopw_d, rbottomw_d
logical :: initialized = .false.

contains

!*******************************************************************************
subroutine press_gpu_init()
!*******************************************************************************
implicit none
if (initialized) return

allocate(rH_x_d(ld, ny, lbz:nz));      rH_x_d = 0._rprec
allocate(rH_y_d(ld, ny, lbz:nz));      rH_y_d = 0._rprec
allocate(rH_z_d(ld, ny, lbz:nz));      rH_z_d = 0._rprec
allocate(rtopw_d(ld, ny));             rtopw_d = 0._rprec
allocate(rbottomw_d(ld, ny));          rbottomw_d = 0._rprec
allocate(RHS_col_d(ld, ny, nz+1));     RHS_col_d = 0._rprec
allocate(a_d(lh, ny, nz+1));           a_d = 0._rprec
allocate(b_d(lh, ny, nz+1));           b_d = 0._rprec
allocate(c_d(lh, ny, nz+1));           c_d = 0._rprec

!$acc enter data copyin(rH_x_d, rH_y_d, rH_z_d,                                &
!$acc                   rtopw_d, rbottomw_d,                                   &
!$acc                   RHS_col_d, a_d, b_d, c_d)

initialized = .true.
end subroutine press_gpu_init

!*******************************************************************************
subroutine press_gpu_finalize()
!*******************************************************************************
implicit none
if (.not. initialized) return
!$acc exit data delete(rH_x_d, rH_y_d, rH_z_d,                                 &
!$acc                  rtopw_d, rbottomw_d,                                    &
!$acc                  RHS_col_d, a_d, b_d, c_d)
deallocate(rH_x_d, rH_y_d, rH_z_d, rtopw_d, rbottomw_d,                        &
           RHS_col_d, a_d, b_d, c_d)
initialized = .false.
end subroutine press_gpu_finalize

!*******************************************************************************
subroutine press_stag_array_gpu()
!*******************************************************************************
implicit none

real(rprec) :: const, const2, const3, const4
integer :: jx, jy, jz, ir, ii
integer :: jz_min

call press_gpu_init()

! Cache constants (mirror CPU)
const  = 1._rprec / real(nx*ny, rprec)
const2 = const / tadv1 / dt
const3 = 1._rprec / (dz**2)
const4 = 1._rprec / dz

! ---- Set up p(:,:,0) on host (CPU lines 66-72); seeded to device after the
!      data region opens (see slab-only update device below). ----
if (coord == 0) then
    p(:,:,0) = 0._rprec
#ifdef PPSAFETYMODE
else
    p(:,:,0) = BOGUS
#endif
end if

! u, v, w are device-resident via sim_param's `!$acc declare create`.
! divtz and dpdx/y/z are also device-resident (sim_param declare create).
! divstress_w_gpu writes divtz on device; the RHS update reads dpdx/y/z on
! device via the GPU AB kernel in main.f90 + project on device.
!
! p is a press-internal scratch. We `create` (not copyin) it on device: the
! only host-set value at entry is the BC ghost slab p(:,:,0) (above), which we
! seed with a slab-only H2D below. Slabs 1..nz are written by tridag (the
! Thomas solve fills u(:,:,1..) as p, with the index shift u(:,:,k)=p(:,:,k-1)
! / the DC-mode column / the 0<->nz-1 halo before any read in Phase E/F - so
! they never need a host seed. (Equivalent to the old copyin, which pushed a
! fresh p(:,:,0) plus STALE slabs 1..nz that were overwritten anyway; this just
! drops the full 34 MB H2D/iter to a 0.5 MB slab.)
!$acc data create(p)                                                           &
!$acc      present(u, v, w, divtz, dpdx, dpdy, dpdz,                           &
!$acc              rH_x_d, rH_y_d, rH_z_d,                                     &
!$acc              rtopw_d, rbottomw_d,                                        &
!$acc              RHS_col_d, a_d, b_d, c_d, kx, ky)

! NOTE (round 3): the old per-step `update device(p(:,:,0))` seed here was a
! dead store - rank 0's tridag boundary row (j=1, which is p slab 0 through
! the u(:,:,k)=p(:,:,k-1) index shift) writes it, and non-bottom ranks
! receive it in the tridag forward sweep and again via the p halo, all
! before the Phase E read. Dropped: saves a 9.4 MB H2D + a queue drain
! every step.

! ============================================================================
! PHASE A: build spectral RHS
! ============================================================================

! A.1: const-multiply rH_x/y/z = const2 * u/v/w  for jz = 1..nz-1
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        do jx = 1, ld
            rH_x_d(jx, jy, jz) = const2 * u(jx, jy, jz)
            rH_y_d(jx, jy, jz) = const2 * v(jx, jy, jz)
            rH_z_d(jx, jy, jz) = const2 * w(jx, jy, jz)
        end do
    end do
end do

! A.2: forward FFT (small, batch=nz-1) on each of rH_x, rH_y, rH_z
call fft_gpu_exec_d2z(plan_forw_small_nzm1, rH_x_d(1, 1, 1))
call fft_gpu_exec_d2z(plan_forw_small_nzm1, rH_y_d(1, 1, 1))
call fft_gpu_exec_d2z(plan_forw_small_nzm1, rH_z_d(1, 1, 1))

! A.3: top-rank rH_z(nz) - single-slab const-mult + FFT
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            rH_z_d(jx, jy, nz) = const2 * w(jx, jy, nz)
        end do
    end do
    call fft_gpu_exec_d2z(plan_forw_small_one, rH_z_d(1, 1, nz))
end if

! A.4: bottom-rank rbottomw = const * divtz(:,:,1), single-slab FFT
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            rbottomw_d(jx, jy) = const * divtz(jx, jy, 1)
        end do
    end do
    call fft_gpu_exec_d2z(plan_forw_small_one, rbottomw_d(1, 1))
end if

! A.5: top-rank rtopw = const * divtz(:,:,nz), single-slab FFT
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            rtopw_d(jx, jy) = const * divtz(jx, jy, nz)
        end do
    end do
    call fft_gpu_exec_d2z(plan_forw_small_one, rtopw_d(1, 1))
end if

! A.6: zero oddballs (Nyquist row + ny_h+1 column) for rH_x/y/z over jz=1..nz-1
!$acc parallel loop collapse(2) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        rH_x_d(ld-1, jy, jz) = 0._rprec; rH_x_d(ld, jy, jz) = 0._rprec
        rH_y_d(ld-1, jy, jz) = 0._rprec; rH_y_d(ld, jy, jz) = 0._rprec
        rH_z_d(ld-1, jy, jz) = 0._rprec; rH_z_d(ld, jy, jz) = 0._rprec
    end do
end do
!$acc parallel loop collapse(2) default(present) async(1)
do jz = 1, nz-1
    do jx = 1, ld
        rH_x_d(jx, ny/2+1, jz) = 0._rprec
        rH_y_d(jx, ny/2+1, jz) = 0._rprec
        rH_z_d(jx, ny/2+1, jz) = 0._rprec
    end do
end do

if (coord == nproc-1) then
    !$acc parallel loop default(present) async(1)
    do jy = 1, ny
        rH_z_d(ld-1, jy, nz) = 0._rprec
        rH_z_d(ld,   jy, nz) = 0._rprec
    end do
    !$acc parallel loop default(present) async(1)
    do jx = 1, ld
        rH_z_d(jx, ny/2+1, nz) = 0._rprec
    end do
end if

! Zero oddballs in rtopw/rbottomw (always do - value is 0 if not coord==0/nproc-1)
!$acc parallel loop default(present) async(1)
do jy = 1, ny
    rtopw_d(ld-1, jy)    = 0._rprec
    rtopw_d(ld,   jy)    = 0._rprec
    rbottomw_d(ld-1, jy) = 0._rprec
    rbottomw_d(ld,   jy) = 0._rprec
end do
!$acc parallel loop default(present) async(1)
do jx = 1, ld
    rtopw_d(jx,    ny/2+1) = 0._rprec
    rbottomw_d(jx, ny/2+1) = 0._rprec
end do

! ============================================================================
! MPI halo exchange of rH_x/y/z (need rH(:,:,0) and rH_z(:,:,nz)).
! GPU-aware builds pass device pointers; non-GPU-aware builds stage slabs
! through host memory in the fallback branch below.
! ============================================================================
#ifdef PPMPI
#ifdef PPGPU_AWARE_MPI
! rH_x_d/y_d/z_d are device-resident (declare-create in fft_gpu); send GPU
! pointers straight to MPICH.
!$acc wait(1)
!$acc host_data use_device(rH_x_d, rH_y_d, rH_z_d)
call mpi_sendrecv(rH_x_d(1, 1, nz-1), ld*ny, MPI_RPREC, up,   1,               &
                  rH_x_d(1, 1, 0),    ld*ny, MPI_RPREC, down, 1,               &
                  comm, status, ierr)
call mpi_sendrecv(rH_y_d(1, 1, nz-1), ld*ny, MPI_RPREC, up,   2,               &
                  rH_y_d(1, 1, 0),    ld*ny, MPI_RPREC, down, 2,               &
                  comm, status, ierr)
call mpi_sendrecv(rH_z_d(1, 1, nz-1), ld*ny, MPI_RPREC, up,   3,               &
                  rH_z_d(1, 1, 0),    ld*ny, MPI_RPREC, down, 3,               &
                  comm, status, ierr)
call mpi_sendrecv(rH_z_d(1, 1, 1),    ld*ny, MPI_RPREC, down, 6,               &
                  rH_z_d(1, 1, nz),   ld*ny, MPI_RPREC, up,   6,               &
                  comm, status, ierr)
!$acc end host_data
#else
! Pull the slabs we need to send back to host
!$acc wait(1)
!$acc update self(rH_x_d(:,:,nz-1), rH_y_d(:,:,nz-1), rH_z_d(:,:,nz-1),        &
!$acc             rH_z_d(:,:,1))

call mpi_sendrecv(rH_x_d(1, 1, nz-1), ld*ny, MPI_RPREC, up,   1,               &
                  rH_x_d(1, 1, 0),    ld*ny, MPI_RPREC, down, 1,               &
                  comm, status, ierr)
call mpi_sendrecv(rH_y_d(1, 1, nz-1), ld*ny, MPI_RPREC, up,   2,               &
                  rH_y_d(1, 1, 0),    ld*ny, MPI_RPREC, down, 2,               &
                  comm, status, ierr)
call mpi_sendrecv(rH_z_d(1, 1, nz-1), ld*ny, MPI_RPREC, up,   3,               &
                  rH_z_d(1, 1, 0),    ld*ny, MPI_RPREC, down, 3,               &
                  comm, status, ierr)
call mpi_sendrecv(rH_z_d(1, 1, 1),    ld*ny, MPI_RPREC, down, 6,               &
                  rH_z_d(1, 1, nz),   ld*ny, MPI_RPREC, up,   6,               &
                  comm, status, ierr)

! Push the received halos back to device
!$acc wait(1)
!$acc update device(rH_x_d(:,:,0), rH_y_d(:,:,0), rH_z_d(:,:,0),               &
!$acc               rH_z_d(:,:,nz))
#endif
#endif

! ============================================================================
! PHASE B: build a, b, c, RHS_col on device
! ============================================================================

! Boundary slab at jz=1 (coord==0). Mirrors CPU lines 149-158. Writes for ALL
! (jx, jy) including jx=jy=1 (DC mode is overwritten later).
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, lh
#ifdef PPSAFETYMODE
            a_d(jx, jy, 1) = BOGUS
#endif
            b_d(jx, jy, 1) = -1._rprec
            c_d(jx, jy, 1) =  1._rprec
        end do
    end do
    ! RHS_col(:,:,1) = -dz * rbottomw(:,:)   (interleaved-real layout)
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            RHS_col_d(jx, jy, 1) = -dz * rbottomw_d(jx, jy)
        end do
    end do
    jz_min = 2
else
    jz_min = 1
end if

! Boundary slab at jz=nz+1 (coord==nproc-1). Mirrors CPU lines 167-175.
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, lh
            a_d(jx, jy, nz+1) = -1._rprec
            b_d(jx, jy, nz+1) =  1._rprec
#ifdef PPSAFETYMODE
            c_d(jx, jy, nz+1) = BOGUS
#endif
        end do
    end do
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            RHS_col_d(jx, jy, nz+1) = -dz * rtopw_d(jx, jy)
        end do
    end do
end if

! Interior slabs jz_min..nz: build a, b, c (real-valued) and RHS_col.
! We compute for ALL (jy, jx) - the CPU code skips (jy=ny/2+1) and (jx*jy=1)
! but tridag_array also skips them, so the unused values are don't-care.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = jz_min, nz
    do jy = 1, ny
        do jx = 1, lh
            ii = 2*jx
            ir = ii - 1
            ! Coefficients (real-valued)
            a_d(jx, jy, jz) = const3
            b_d(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2 + 2._rprec*const3)
            c_d(jx, jy, jz) = const3
            ! RHS_col (interleaved-real); inlined MULI of rH_x and rH_y at jz-1.
            !   real:  -rH_x(ii, jz-1)*kx - rH_y(ii, jz-1)*ky + (rH_z(ir,jz)-rH_z(ir,jz-1))/dz
            !   imag:   rH_x(ir, jz-1)*kx + rH_y(ir, jz-1)*ky + (rH_z(ii,jz)-rH_z(ii,jz-1))/dz
            RHS_col_d(ir, jy, jz) =                                            &
                  - rH_x_d(ii, jy, jz-1) * kx(jx, jy)                          &
                  - rH_y_d(ii, jy, jz-1) * ky(jx, jy)                          &
                  + (rH_z_d(ir, jy, jz) - rH_z_d(ir, jy, jz-1)) * const4
            RHS_col_d(ii, jy, jz) =                                            &
                    rH_x_d(ir, jy, jz-1) * kx(jx, jy)                          &
                  + rH_y_d(ir, jy, jz-1) * ky(jx, jy)                          &
                  + (rH_z_d(ii, jy, jz) - rH_z_d(ii, jy, jz-1)) * const4
        end do
    end do
end do

! ============================================================================
! PHASE C: tridiagonal solve on GPU.
! tridag_array_gpu uses the device-resident a/b/c/RHS_col directly. Its MPI
! chain uses GPU-aware rank-boundary messages when PPGPU_AWARE_MPI is enabled,
! otherwise it falls back to the host-staged path in tridag_gpu_m. cuSPARSE is
! not used because each global tridiagonal system spans the z-decomposed ranks.
! ============================================================================
call tridag_array_gpu(a_d, b_d, c_d, RHS_col_d, p)

! ============================================================================
! PHASE D: zero-wavenumber (DC mode). The recurrence is sequential along jz and
! touches ONLY the DC column (jx=1:2, jy=1), so we keep it on the host but move
! just that column across PCIe instead of the whole (ld,ny,0:nz) p array. The
! oddball-zero and the 0<->nz-1 halo are done on the device. This eliminates
! three full-p transfers/iter (~3 x 34 MB at 256^3/4) - see OPTIMIZATION_PLAN.
!
! tridag's DC-column result is singular and discarded - the DC mode recomputes
! p(1:2,1,:) from rH_z / rbottomw, so we never need tridag's p on the host.
! ============================================================================
! Pull only the DC-column inputs (jx=1:2, jy=1) - ~1 KB instead of 34 MB.
!$acc wait(1)
!$acc update self(rH_z_d(1:2,1,:), rbottomw_d(1:2,1))

#ifdef PPMPI
! wait for p(1, 1, 1) from "down"
call mpi_recv(p(1:2, 1, 1), 2, MPI_RPREC, down, 8, comm, status, ierr)
#endif

if (coord == 0) then
    p(1:2, 1, 0) = 0._rprec
    p(1:2, 1, 1) = p(1:2, 1, 0) - dz * rbottomw_d(1:2, 1)
end if

do jz = 2, nz
    p(1:2, 1, jz) = p(1:2, 1, jz-1) + rH_z_d(1:2, 1, jz) * dz
end do

#ifdef PPMPI
! send p(1, 1, nz) to "up"
call mpi_send(p(1:2, 1, nz), 2, MPI_RPREC, up, 8, comm, ierr)
#endif

! Push the DC column back to device so the halo exchange + Phase E/F see the
! DC-mode-corrected p. The rest of p is already device-resident from tridag.
!$acc wait(1)
!$acc update device(p(1:2,1,:))

#ifdef PPMPI
! sync 0 <-> nz-1 (ghost slab read by Phase F dpdz at jz=1) - slab-only, on
! device. Mirrors the host sendrecv that used to ride on the full-p transfer.
#ifdef PPGPU_AWARE_MPI
!$acc wait(1)
!$acc host_data use_device(p)
call mpi_sendrecv(p(1, 1, nz-1), ld*ny, MPI_RPREC, up,   2,                    &
                  p(1, 1, 0),    ld*ny, MPI_RPREC, down, 2,                    &
                  comm, status, ierr)
!$acc end host_data
#else
!$acc wait(1)
!$acc update self(p(:,:,nz-1))
call mpi_sendrecv(p(1, 1, nz-1), ld*ny, MPI_RPREC, up,   2,                    &
                  p(1, 1, 0),    ld*ny, MPI_RPREC, down, 2,                    &
                  comm, status, ierr)
!$acc wait(1)
!$acc update device(p(:,:,0))
#endif
#endif

! Zero p oddballs on device over all z slabs (mirrors host p(ld-1:ld,:,:)=0
! and p(:,ny/2+1,:)=0 over the full 0:nz range, done AFTER the halo so the
! received p(:,:,0) oddballs get zeroed too - same order as the CPU path).
!$acc parallel loop collapse(2) default(present) async(1)
do jz = 0, nz
    do jy = 1, ny
        p(ld-1, jy, jz) = 0._rprec
        p(ld,   jy, jz) = 0._rprec
    end do
end do
!$acc parallel loop collapse(2) default(present) async(1)
do jz = 0, nz
    do jx = 1, ld
        p(jx, ny/2+1, jz) = 0._rprec
    end do
end do

! ============================================================================
! PHASE E: build dpdx, dpdy from spectral p; inverse FFT all three.
! ============================================================================
! Build dpdx, dpdy via inlined MULI for jz = 1..nz-1 (mirror CPU lines 256-265).
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        do jx = 1, lh
            ii = 2*jx
            ir = ii - 1
            dpdx(ir, jy, jz) = -p(ii, jy, jz) * kx(jx, jy)
            dpdx(ii, jy, jz) =  p(ir, jy, jz) * kx(jx, jy)
            dpdy(ir, jy, jz) = -p(ii, jy, jz) * ky(jx, jy)
            dpdy(ii, jy, jz) =  p(ir, jy, jz) * ky(jx, jy)
        end do
    end do
end do

! Inverse FFT p(:,:,0) - single slab.
call fft_gpu_exec_z2d(plan_back_small_one, p(1, 1, 0))

! Batched inverse FFT for jz=1..nz-1 of (p, dpdx, dpdy)
call fft_gpu_exec_z2d(plan_back_small_nzm1, p(1,    1, 1))
call fft_gpu_exec_z2d(plan_back_small_nzm1, dpdx(1, 1, 1))
call fft_gpu_exec_z2d(plan_back_small_nzm1, dpdy(1, 1, 1))

! Top rank: also inverse FFT p(:,:,nz)
if (coord == nproc-1) then
    call fft_gpu_exec_z2d(plan_back_small_one, p(1, 1, nz))
end if

#ifdef PPSAFETYMODE
! dpdx(:,:,nz), dpdy(:,:,nz) are not needed elsewhere - set to BOGUS for safety.
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        dpdx(jx, jy, nz) = BOGUS
        dpdy(jx, jy, nz) = BOGUS
    end do
end do
if (coord < nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            p(jx, jy, nz) = BOGUS
        end do
    end do
end if
#endif

! ============================================================================
! PHASE F: dpdz = (p(jz) - p(jz-1)) / dz
! ============================================================================
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        do jx = 1, nx
            dpdz(jx, jy, jz) = (p(jx, jy, jz) - p(jx, jy, jz-1)) / dz
        end do
    end do
end do

! The lower pressure boundary is homogeneous Neumann.  The host pressure path
! produces p(:,:,1) == p(:,:,0), so its bottom w-grid gradient is exactly zero.
! Enforce the same public dpdz contract after the device inverse transforms;
! RHSz(:,:,1) is a boundary value, but diagnostics and optional consumers still
! observe this array.
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            dpdz(jx, jy, 1) = 0._rprec
        end do
    end do
end if

#ifdef PPSAFETYMODE
if (coord < nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            dpdz(jx, jy, nz) = BOGUS
        end do
    end do
end if
#endif

if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            dpdz(jx, jy, nz) = (p(jx, jy, nz) - p(jx, jy, nz-1)) / dz
        end do
    end do
end if

! Refresh host p ONLY when domain output (io.f90 inst_write itype=2,
! pres_real = p - 0.5*(u^2+v^2+w^2)) will read it this step. That domain branch
! is the sole consumer of host-resident p: checkpoint does not write p, and
! point/x/y/z-plane output do not use it. Condition matches output_loop's
! domain gate exactly. Skipping it on non-output steps avoids a full
! (ld,ny,0:nz) D2H every iteration (in this config domain_nstart > nsteps, so
! it never fires). p is not device-resident outside this data region, so this
! is the only place the host copy can be refreshed.
if (domain_calc .and. jt_total >= domain_nstart .and.                          &
    jt_total <= domain_nend .and.                                              &
    mod(jt_total - domain_nstart, domain_nskip) == 0) then
    !$acc wait(1)
    !$acc update self(p)
end if

!$acc wait(1)
!$acc end data    ! create(p) - no copyout; dpdx/y/z are declare-create (present)

end subroutine press_stag_array_gpu

#endif
end module press_gpu_m
