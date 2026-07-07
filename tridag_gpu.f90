!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!  GPU port: see docs/gpu_module_contracts.md for current ownership context.
!!

!*******************************************************************************
module tridag_gpu_m
!*******************************************************************************
! GPU implementation of tridag_array.
!
! Algorithmically identical to tridag_array.f90's MPI version: distributed
! Thomas algorithm with MPI-chained forward sweep and backward substitution
! across the z-decomposed ranks. The local sweeps run on the GPU.
!
! Why no cuSPARSE: the (jx, jy) systems are NOT independent across MPI ranks
! along z (each global system spans all nproc slabs). cuSPARSE only solves
! complete independent systems, so it would require gather-solve-scatter to
! a single rank - which is *slower* than the chained Thomas due to data
! movement (~270 MB to gather a/b/c/RHS to rank 0 without GPU-aware MPI).
!
! Parallelism strategy: one CUDA thread per (jx, jy) system, each thread
! sequential along j.
!
! MPI-chain strategy (PPGPU_AWARE_MPI): the (jx, jy) systems are mutually
! independent, so the jy range is split into nchunk_pipeline blocks that move
! through the rank chain in a pipeline - rank k computes chunk m while rank
! k-1 computes chunk m+1. Boundary data (c, bet, u at slab edges) travels as
! GPU-aware MPI messages (device pointers straight to MPICH, no host bounce).
! This shrinks the serial pipeline bubble from nproc*(full sweep + 18.6 MB
! host-staged message) to roughly (nproc + nchunk)*(chunk sweep + chunk
! message): ~10-20 ms instead of ~150-250 ms at nproc=16 for 3072x384x400.
! Per-system arithmetic and its ordering are untouched - results are
! bit-identical to the monolithic sweep.
!
! Runtime controls (read once, must be identical on all ranks; srun
! --export=ALL takes care of that):
!   LESGO_TRIDAG_NCHUNK   pipeline chunk count (default 8)
!   LESGO_TRIDAG_GPU_MPI  set to 0 to fall back to the monolithic
!                         host-staged chain at runtime (A/B debugging)
!
! Without PPGPU_AWARE_MPI the monolithic host-staged chain is the only path.
!*******************************************************************************
#ifdef PPPRESS_GPU
use types, only : rprec
use param, only : ld, lh, ny, nz, nproc, coord
#ifdef PPMPI
use param, only : MPI_RPREC, comm, ierr, status, up, down
use mpi, only : MPI_STATUSES_IGNORE, mpi_irecv, mpi_isend, mpi_recv, mpi_send, &
    mpi_wait, mpi_waitall
#endif
implicit none
save
private
public :: tridag_array_gpu, tridag_gpu_finalize

! Persistent device-resident scratch:
!   bet_d holds the running diagonal in the modified-Thomas recursion.
!   gam_d holds the multiplier from forward sweep, needed in backward sub.
real(rprec), allocatable, dimension(:,:)   :: bet_d
real(rprec), allocatable, dimension(:,:,:) :: gam_d
logical :: initialized = .false.

! Pipeline chunking over jy (GPU-aware path only)
integer, parameter :: MAXCHUNK = 64
integer :: nchunk_pipeline = 8
logical :: pipeline_enabled = .true.

contains

!*******************************************************************************
integer function tridag_gpu_env_positive_int(name, default_value)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
integer, intent(in) :: default_value
character(64) :: envstr
integer :: envval, envios

tridag_gpu_env_positive_int = default_value
envstr = ''
call get_environment_variable(name, value=envstr)
if (len_trim(envstr) > 0) then
    read(envstr, *, iostat=envios) envval
    if (envios == 0 .and. envval >= 1)                                &
        tridag_gpu_env_positive_int = envval
end if

end function tridag_gpu_env_positive_int

!*******************************************************************************
logical function tridag_gpu_env_exact_zero(name)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
character(64) :: envstr

envstr = ''
call get_environment_variable(name, value=envstr)
tridag_gpu_env_exact_zero = trim(envstr) == '0'

end function tridag_gpu_env_exact_zero

!*******************************************************************************
subroutine tridag_gpu_init()
!*******************************************************************************
implicit none
if (initialized) return
allocate(bet_d(lh, ny));         bet_d = 0._rprec
allocate(gam_d(lh, ny, nz+1));   gam_d = 0._rprec
!$acc enter data copyin(bet_d, gam_d)

nchunk_pipeline = tridag_gpu_env_positive_int('LESGO_TRIDAG_NCHUNK',          &
    nchunk_pipeline)
nchunk_pipeline = max(1, min(nchunk_pipeline, ny, MAXCHUNK))

if (tridag_gpu_env_exact_zero('LESGO_TRIDAG_GPU_MPI'))                       &
    pipeline_enabled = .false.

initialized = .true.
end subroutine tridag_gpu_init

!*******************************************************************************
subroutine tridag_gpu_finalize()
!*******************************************************************************
implicit none
if (.not. initialized) return
!$acc exit data delete(bet_d, gam_d)
deallocate(bet_d, gam_d)
initialized = .false.
end subroutine tridag_gpu_finalize

!*******************************************************************************
subroutine tridag_array_gpu(a, b, c, r, u)
!*******************************************************************************
! Drop-in replacement for tridag_array. Inputs a, b, c, r and the output u
! must already be device-resident (caller's responsibility).
!
! Note: c is INPUT but we receive INTO c(:,:,1) on non-bottom ranks (mirrors
! CPU tridag_array which also overwrites the local c(:,:,1) with the value
! from below). This is part of the algorithm - c then re-reads the received
! value as c[j-1] in the first forward step.
!*******************************************************************************
implicit none
real(rprec), dimension(lh, ny, nz+1), intent(inout) :: a, b, c    ! intent in fact, but c gets recv-overwritten
real(rprec), dimension(ld, ny, nz+1), intent(in)    :: r
real(rprec), dimension(ld, ny, nz+1), intent(inout) :: u

integer :: jx, jy, j, ii, ir
integer :: j_max, j_min
real(rprec) :: bet_local
#if defined(PPMPI) && defined(PPGPU_AWARE_MPI)
integer :: kchunk, jy0, jy1, clen, nck, nsb
integer :: rreq_f(3, MAXCHUNK), rreq_b(2, MAXCHUNK)
integer :: sreq_f(3, MAXCHUNK), sreq_b(2*MAXCHUNK)
#endif

call tridag_gpu_init()

! Per-rank j range
if (coord == nproc-1) then
    j_max = nz + 1
else
    j_max = nz
end if
if (coord == 0) then
    j_min = 1
else
    j_min = 2
end if

#if defined(PPMPI) && defined(PPGPU_AWARE_MPI)
if (pipeline_enabled) then
! ============================================================================
! GPU-AWARE CHUNK-PIPELINED THOMAS
! ============================================================================
nck = nchunk_pipeline
nsb = 0

! All prior async(1) work (the a/b/c/RHS build kernels) must be complete
! before MPI is allowed to land data in the device buffers.
!$acc wait(1)

! Pre-post every per-chunk receive up front so the network can deliver chunk
! m+1 while chunk m's sweep kernel runs. Same-tag messages between a rank
! pair match in posting order, which equals chunk order on both sides.
!$acc host_data use_device(c, bet_d, u, gam_d)
if (coord /= 0) then
    do kchunk = 1, nck
        jy0 = 1 + ((kchunk-1)*ny)/nck
        jy1 = (kchunk*ny)/nck
        clen = jy1 - jy0 + 1
        call mpi_irecv(c(1, jy0, 1),  lh*clen, MPI_RPREC, down, 1, comm,       &
                       rreq_f(1, kchunk), ierr)
        call mpi_irecv(bet_d(1, jy0), lh*clen, MPI_RPREC, down, 2, comm,       &
                       rreq_f(2, kchunk), ierr)
        call mpi_irecv(u(1, jy0, 1),  ld*clen, MPI_RPREC, down, 3, comm,       &
                       rreq_f(3, kchunk), ierr)
    end do
end if
if (coord /= nproc-1) then
    do kchunk = 1, nck
        jy0 = 1 + ((kchunk-1)*ny)/nck
        jy1 = (kchunk*ny)/nck
        clen = jy1 - jy0 + 1
        call mpi_irecv(u(1, jy0, nz+1),     ld*clen, MPI_RPREC, up, 4, comm,   &
                       rreq_b(1, kchunk), ierr)
        call mpi_irecv(gam_d(1, jy0, nz+1), lh*clen, MPI_RPREC, up, 5, comm,   &
                       rreq_b(2, kchunk), ierr)
    end do
end if
!$acc end host_data

! ---- forward sweep, chunk by chunk ----
! The top rank needs no backward receive (it owns j=nz+1), so its backward
! sweep is fused into the same chunk loop: the backward wave starts down the
! chain as soon as chunk 1 is done, overlapping the remaining forward chunks.
do kchunk = 1, nck
    jy0 = 1 + ((kchunk-1)*ny)/nck
    jy1 = (kchunk*ny)/nck
    clen = jy1 - jy0 + 1

    if (coord == 0) then
        ! Initialize j=1: bet = b(:,:,1), u = r/b. Mirrors CPU lines 52-66.
        !$acc parallel loop collapse(2) default(present) async(1)
        do jy = jy0, jy1
            do jx = 1, lh
                ii = 2*jx
                ir = ii - 1
                bet_d(jx, jy)    = b(jx, jy, 1)
                u(ir, jy, 1)     = r(ir, jy, 1) / bet_d(jx, jy)
                u(ii, jy, 1)     = r(ii, jy, 1) / bet_d(jx, jy)
            end do
        end do
    else
        call mpi_wait(rreq_f(1, kchunk), status, ierr)
        call mpi_wait(rreq_f(2, kchunk), status, ierr)
        call mpi_wait(rreq_f(3, kchunk), status, ierr)
    end if

    ! Forward sweep j = 2 .. j_max for this chunk's (jx, jy) systems. Each
    ! thread holds bet in a private register through the inner j loop and
    ! stores back to bet_d at the end (so the send can pick it up).
    !$acc parallel loop collapse(2) private(bet_local, ii, ir, j) default(present) async(1)
    do jy = jy0, jy1
        do jx = 1, lh
            ii = 2*jx
            ir = ii - 1
            bet_local = bet_d(jx, jy)
            do j = 2, j_max
                gam_d(jx, jy, j) = c(jx, jy, j-1) / bet_local
                bet_local = b(jx, jy, j) - a(jx, jy, j) * gam_d(jx, jy, j)
                u(ir, jy, j) = (r(ir, jy, j) - a(jx, jy, j) * u(ir, jy, j-1)) / bet_local
                u(ii, jy, j) = (r(ii, jy, j) - a(jx, jy, j) * u(ii, jy, j-1)) / bet_local
            end do
            bet_d(jx, jy) = bet_local
        end do
    end do

    if (coord == nproc-1) then
        ! Fused backward sub for the top rank (no receive needed).
        !$acc parallel loop collapse(2) private(ii, ir, j) default(present) async(1)
        do jy = jy0, jy1
            do jx = 1, lh
                ii = 2*jx
                ir = ii - 1
                do j = nz, j_min, -1
                    u(ir, jy, j) = u(ir, jy, j) - gam_d(jx, jy, j+1) * u(ir, jy, j+1)
                    u(ii, jy, j) = u(ii, jy, j) - gam_d(jx, jy, j+1) * u(ii, jy, j+1)
                end do
            end do
        end do
        if (coord /= 0) then
            !$acc wait(1)
            !$acc host_data use_device(u, gam_d)
            call mpi_isend(u(1, jy0, 2),     ld*clen, MPI_RPREC, down, 4, comm, &
                           sreq_b(nsb+1), ierr)
            call mpi_isend(gam_d(1, jy0, 2), lh*clen, MPI_RPREC, down, 5, comm, &
                           sreq_b(nsb+2), ierr)
            !$acc end host_data
            nsb = nsb + 2
        end if
    else
        !$acc wait(1)
        !$acc host_data use_device(c, bet_d, u)
        call mpi_isend(c(1, jy0, nz),  lh*clen, MPI_RPREC, up, 1, comm,        &
                       sreq_f(1, kchunk), ierr)
        call mpi_isend(bet_d(1, jy0),  lh*clen, MPI_RPREC, up, 2, comm,        &
                       sreq_f(2, kchunk), ierr)
        call mpi_isend(u(1, jy0, nz),  ld*clen, MPI_RPREC, up, 3, comm,        &
                       sreq_f(3, kchunk), ierr)
        !$acc end host_data
    end if
end do

! ---- backward substitution for non-top ranks, chunk by chunk ----
if (coord /= nproc-1) then
    do kchunk = 1, nck
        jy0 = 1 + ((kchunk-1)*ny)/nck
        jy1 = (kchunk*ny)/nck
        clen = jy1 - jy0 + 1

        ! The backward kernel below overwrites u(:,jy chunk,nz), which the
        ! forward isend of this chunk read. By causality those sends have
        ! already been consumed upstream (the backward data we are about to
        ! receive depends on them), so these waits return immediately - they
        ! exist to satisfy the MPI rule that a send buffer may not be
        ! modified before its request completes.
        call mpi_wait(sreq_f(1, kchunk), status, ierr)
        call mpi_wait(sreq_f(2, kchunk), status, ierr)
        call mpi_wait(sreq_f(3, kchunk), status, ierr)

        call mpi_wait(rreq_b(1, kchunk), status, ierr)
        call mpi_wait(rreq_b(2, kchunk), status, ierr)

        !$acc parallel loop collapse(2) private(ii, ir, j) default(present) async(1)
        do jy = jy0, jy1
            do jx = 1, lh
                ii = 2*jx
                ir = ii - 1
                do j = nz, j_min, -1
                    u(ir, jy, j) = u(ir, jy, j) - gam_d(jx, jy, j+1) * u(ir, jy, j+1)
                    u(ii, jy, j) = u(ii, jy, j) - gam_d(jx, jy, j+1) * u(ii, jy, j+1)
                end do
            end do
        end do

        if (coord /= 0) then
            !$acc wait(1)
            !$acc host_data use_device(u, gam_d)
            call mpi_isend(u(1, jy0, 2),     ld*clen, MPI_RPREC, down, 4, comm, &
                           sreq_b(nsb+1), ierr)
            call mpi_isend(gam_d(1, jy0, 2), lh*clen, MPI_RPREC, down, 5, comm, &
                           sreq_b(nsb+2), ierr)
            !$acc end host_data
            nsb = nsb + 2
        end if
    end do
end if

if (nsb > 0) call mpi_waitall(nsb, sreq_b, MPI_STATUSES_IGNORE, ierr)

return
end if
#endif

! ============================================================================
! FALLBACK: monolithic host-staged chain (original implementation)
! ============================================================================
! FORWARD SWEEP
if (coord == 0) then
    ! Initialize j=1: bet = b(:,:,1), u = r/b. Mirrors CPU lines 52-66.
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, lh
            ii = 2*jx
            ir = ii - 1
            bet_d(jx, jy)    = b(jx, jy, 1)
            u(ir, jy, 1)     = r(ir, jy, 1) / bet_d(jx, jy)
            u(ii, jy, 1)     = r(ii, jy, 1) / bet_d(jx, jy)
        end do
    end do
#ifdef PPMPI
else
    ! Receive c(:,:,1), bet, u(:,:,1) from down rank (overwrites this rank's
    ! c(:,:,1) - that's the algorithm). Receive into HOST buffers, then push
    ! to device.
    call mpi_recv(c(1, 1, 1),  lh*ny, MPI_RPREC, down, 1, comm, status, ierr)
    call mpi_recv(bet_d(1, 1), lh*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_recv(u(1, 1, 1),  ld*ny, MPI_RPREC, down, 3, comm, status, ierr)
    !$acc wait(1)
    !$acc update device(c(:,:,1), bet_d, u(:,:,1))
#endif
end if

! Forward sweep j = 2 .. j_max for ALL (jx, jy). Each thread holds bet in a
! private register through the inner j loop and stores back to bet_d at the
! end (so the MPI send can pick it up).
!$acc parallel loop collapse(2) private(bet_local, ii, ir, j) default(present) async(1)
do jy = 1, ny
    do jx = 1, lh
        ii = 2*jx
        ir = ii - 1
        bet_local = bet_d(jx, jy)
        do j = 2, j_max
            gam_d(jx, jy, j) = c(jx, jy, j-1) / bet_local
            bet_local = b(jx, jy, j) - a(jx, jy, j) * gam_d(jx, jy, j)
            u(ir, jy, j) = (r(ir, jy, j) - a(jx, jy, j) * u(ir, jy, j-1)) / bet_local
            u(ii, jy, j) = (r(ii, jy, j) - a(jx, jy, j) * u(ii, jy, j-1)) / bet_local
        end do
        bet_d(jx, jy) = bet_local
    end do
end do

! Send up to next rank (no-op when coord==nproc-1 because `up` is MPI_PROC_NULL)
#ifdef PPMPI
if (coord /= nproc-1) then
    !$acc wait(1)
    !$acc update self(c(:,:,nz), bet_d, u(:,:,nz))
    call mpi_send(c(1, 1, nz),  lh*ny, MPI_RPREC, up, 1, comm, ierr)
    call mpi_send(bet_d(1, 1),  lh*ny, MPI_RPREC, up, 2, comm, ierr)
    call mpi_send(u(1, 1, nz),  ld*ny, MPI_RPREC, up, 3, comm, ierr)
end if
#endif

! BACKWARD SUBSTITUTION
#ifdef PPMPI
if (coord /= nproc-1) then
    ! Receive u(:,:,nz+1), gam(:,:,nz+1) from up rank
    call mpi_recv(u(1, 1, nz+1),     ld*ny, MPI_RPREC, up, 4, comm, status, ierr)
    call mpi_recv(gam_d(1, 1, nz+1), lh*ny, MPI_RPREC, up, 5, comm, status, ierr)
    !$acc wait(1)
    !$acc update device(u(:,:,nz+1), gam_d(:,:,nz+1))
end if
#endif

! Backward sub j = nz .. j_min for ALL (jx, jy)
!$acc parallel loop collapse(2) private(ii, ir, j) default(present) async(1)
do jy = 1, ny
    do jx = 1, lh
        ii = 2*jx
        ir = ii - 1
        do j = nz, j_min, -1
            u(ir, jy, j) = u(ir, jy, j) - gam_d(jx, jy, j+1) * u(ir, jy, j+1)
            u(ii, jy, j) = u(ii, jy, j) - gam_d(jx, jy, j+1) * u(ii, jy, j+1)
        end do
    end do
end do

! Send down to previous rank (no-op when coord==0 because `down` is MPI_PROC_NULL)
#ifdef PPMPI
if (coord /= 0) then
    !$acc wait(1)
    !$acc update self(u(:,:,2), gam_d(:,:,2))
    call mpi_send(u(1, 1, 2),     ld*ny, MPI_RPREC, down, 4, comm, ierr)
    call mpi_send(gam_d(1, 1, 2), lh*ny, MPI_RPREC, down, 5, comm, ierr)
end if
#endif

end subroutine tridag_array_gpu

#endif
end module tridag_gpu_m
