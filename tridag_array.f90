#if !defined(PPPRESS_GPU) || (defined(PPLVLSET) && defined(PPLES_GPU) && !defined(ENABLE_CUDA))
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

! Navigation map:
!   - tridag_transpose_gpu_mod: CUDA packing kernels for transpose-based solve
!   - tridag_array: public pressure tridiagonal entry point
!   - environment helpers: tridag_apply_env_* diagnostics/policy parsing
!   - GPU solver variants: transpose Thomas, spike2, and replicated paths
!   - CPU fallback: final tridag_array implementation in the non-GPU branch
!
! The pressure solver depends on MPI decomposition and GPU-aware MPI behavior.
! Keep host-staged fallbacks correct, but document any change that alters the
! default GPU solver selection or data movement pattern.
!
! Ownership map:
!   - press_stag_array.f90 owns pressure RHS construction and calls this solver.
!   - tridag_array owns CPU/fallback orchestration, transpose staging, and
!     policy selection among GPU-aware, host-staged, and replicated variants.
!   - tridag_gpu.f90 owns lower-level GPU tridiagonal helper kernels.
!   - MPI buffers here are pressure-solver-private; do not reuse them as
!     general halo storage outside the pressure path.

#ifdef ENABLE_CUDA
module tridag_transpose_gpu_mod
use types, only : rprec
use cudafor
implicit none
contains

attributes(global) subroutine tr_pack_in_kernel(r, sol, sbuf, ldv, nyv,       &
    nymv, send_jlo, own_jlo, rem_jlo, local_glo, total)
implicit none
integer, value :: ldv, nyv, nymv, send_jlo, own_jlo, rem_jlo, local_glo
integer, value :: total
real(rprec), device :: r(ldv,nyv,*)
real(rprec), device :: sol(ldv,nymv,*)
real(rprec), device :: sbuf(ldv,nymv,*)
integer :: tid, k, jyl, j

tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
if (tid < 0 .or. tid >= total) return

k = mod(tid, ldv) + 1
jyl = mod(tid / ldv, nymv) + 1
j = tid / (ldv * nymv) + 1
sol(k,jyl,local_glo + j - 1) = r(k,own_jlo + jyl - 1, send_jlo + j - 1)
sbuf(k,jyl,j) = r(k,rem_jlo + jyl - 1, send_jlo + j - 1)

end subroutine tr_pack_in_kernel

attributes(global) subroutine tr_unpack_in_kernel(rbuf, sol, ldv, nymv,       &
    rem_glo, total)
implicit none
integer, value :: ldv, nymv, rem_glo, total
real(rprec), device :: rbuf(ldv,nymv,*)
real(rprec), device :: sol(ldv,nymv,*)
integer :: tid, k, jyl, j

tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
if (tid < 0 .or. tid >= total) return

k = mod(tid, ldv) + 1
jyl = mod(tid / ldv, nymv) + 1
j = tid / (ldv * nymv) + 1
sol(k,jyl,rem_glo + j - 1) = rbuf(k,jyl,j)

end subroutine tr_unpack_in_kernel

attributes(global) subroutine tr_pack_out_kernel(sol, u, sbuf, ldv, nyv,      &
    nymv, nv, nzv, coordv, own_jlo, rem_coord, total)
implicit none
integer, value :: ldv, nyv, nymv, nv, nzv, coordv, own_jlo, rem_coord, total
real(rprec), device :: sol(ldv,nymv,*)
real(rprec), device :: u(ldv,nyv,*)
real(rprec), device :: sbuf(ldv,nymv,*)
integer :: tid, k, jyl, j, jy, g

tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
if (tid < 0 .or. tid >= total) return

k = mod(tid, ldv) + 1
jyl = mod(tid / ldv, nymv) + 1
j = tid / (ldv * nymv) + 1

g = coordv * (nzv - 1) + j
jy = own_jlo + jyl - 1
if ((jy == 1) .and. (k <= 2)) then
    ! Zero mode is handled by the dedicated pressure zero-mode pass.
else if ((k >= ldv - 1) .or. (jy == nyv/2 + 1)) then
    u(k,jy,j) = 0._rprec
else
    u(k,jy,j) = sol(k,jyl,g)
end if

g = rem_coord * (nzv - 1) + j
sbuf(k,jyl,j) = sol(k,jyl,g)

end subroutine tr_pack_out_kernel

attributes(global) subroutine tr_unpack_out_kernel(rbuf, u, ldv, nyv, nymv,   &
    rem_jlo, total)
implicit none
integer, value :: ldv, nyv, nymv, rem_jlo, total
real(rprec), device :: rbuf(ldv,nymv,*)
real(rprec), device :: u(ldv,nyv,*)
integer :: tid, k, jyl, j, jy

tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
if (tid < 0 .or. tid >= total) return

k = mod(tid, ldv) + 1
jyl = mod(tid / ldv, nymv) + 1
j = tid / (ldv * nymv) + 1
jy = rem_jlo + jyl - 1

if ((jy == 1) .and. (k <= 2)) then
    ! Zero mode is handled by the dedicated pressure zero-mode pass.
else if ((k >= ldv - 1) .or. (jy == nyv/2 + 1)) then
    u(k,jy,j) = 0._rprec
else
    u(k,jy,j) = rbuf(k,jyl,j)
end if

end subroutine tr_unpack_out_kernel

attributes(global) subroutine tr_self_copy_in_2_kernel(r, sol, ldv, nyv,     &
    nymv, send_jlo, own_jlo, local_glo, nplanes)
implicit none
integer, value :: ldv, nyv, nymv, send_jlo, own_jlo, local_glo, nplanes
real(rprec), device :: r(ldv,nyv,*)
real(rprec), device :: sol(ldv,nymv,*)
integer :: k, jyl, j

k = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jyl = (blockIdx%y - 1) * blockDim%y + threadIdx%y
j = (blockIdx%z - 1) * blockDim%z + threadIdx%z
if (k > ldv .or. jyl > nymv .or. j > nplanes) return

sol(k,jyl,local_glo + j - 1) = r(k,own_jlo + jyl - 1, send_jlo + j - 1)

end subroutine tr_self_copy_in_2_kernel

attributes(global) subroutine tr_pack_remote_in_2_kernel(r, sbuf, ldv, nyv,  &
    nymv, send_jlo, rem_jlo, nplanes)
implicit none
integer, value :: ldv, nyv, nymv, send_jlo, rem_jlo, nplanes
real(rprec), device :: r(ldv,nyv,*)
real(rprec), device :: sbuf(ldv,nymv,*)
integer :: k, jyl, j

k = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jyl = (blockIdx%y - 1) * blockDim%y + threadIdx%y
j = (blockIdx%z - 1) * blockDim%z + threadIdx%z
if (k > ldv .or. jyl > nymv .or. j > nplanes) return

sbuf(k,jyl,j) = r(k,rem_jlo + jyl - 1, send_jlo + j - 1)

end subroutine tr_pack_remote_in_2_kernel

attributes(global) subroutine tr_unpack_in_2_kernel(rbuf, sol, ldv, nymv,    &
    rem_glo, nplanes)
implicit none
integer, value :: ldv, nymv, rem_glo, nplanes
real(rprec), device :: rbuf(ldv,nymv,*)
real(rprec), device :: sol(ldv,nymv,*)
integer :: k, jyl, j

k = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jyl = (blockIdx%y - 1) * blockDim%y + threadIdx%y
j = (blockIdx%z - 1) * blockDim%z + threadIdx%z
if (k > ldv .or. jyl > nymv .or. j > nplanes) return

sol(k,jyl,rem_glo + j - 1) = rbuf(k,jyl,j)

end subroutine tr_unpack_in_2_kernel

attributes(global) subroutine tr_pack_out_2_kernel(sol, u, sbuf, ldv, nyv,   &
    nymv, nv, nzv, coordv, own_jlo, rem_coord)
implicit none
integer, value :: ldv, nyv, nymv, nv, nzv, coordv, own_jlo, rem_coord
real(rprec), device :: sol(ldv,nymv,*)
real(rprec), device :: u(ldv,nyv,*)
real(rprec), device :: sbuf(ldv,nymv,*)
integer :: k, jyl, j, jy, g

k = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jyl = (blockIdx%y - 1) * blockDim%y + threadIdx%y
j = (blockIdx%z - 1) * blockDim%z + threadIdx%z
if (k > ldv .or. jyl > nymv .or. j > nv) return

g = coordv * (nzv - 1) + j
jy = own_jlo + jyl - 1
if ((jy == 1) .and. (k <= 2)) then
    ! Zero mode is handled by the dedicated pressure zero-mode pass.
else if ((k >= ldv - 1) .or. (jy == nyv/2 + 1)) then
    u(k,jy,j) = 0._rprec
else
    u(k,jy,j) = sol(k,jyl,g)
end if

g = rem_coord * (nzv - 1) + j
sbuf(k,jyl,j) = sol(k,jyl,g)

end subroutine tr_pack_out_2_kernel

attributes(global) subroutine tr_unpack_out_2_kernel(rbuf, u, ldv, nyv,      &
    nymv, rem_jlo, nv)
implicit none
integer, value :: ldv, nyv, nymv, rem_jlo, nv
real(rprec), device :: rbuf(ldv,nymv,*)
real(rprec), device :: u(ldv,nyv,*)
integer :: k, jyl, j, jy

k = (blockIdx%x - 1) * blockDim%x + threadIdx%x
jyl = (blockIdx%y - 1) * blockDim%y + threadIdx%y
j = (blockIdx%z - 1) * blockDim%z + threadIdx%z
if (k > ldv .or. jyl > nymv .or. j > nv) return

jy = rem_jlo + jyl - 1
if ((jy == 1) .and. (k <= 2)) then
    ! Zero mode is handled by the dedicated pressure zero-mode pass.
else if ((k >= ldv - 1) .or. (jy == nyv/2 + 1)) then
    u(k,jy,j) = 0._rprec
else
    u(k,jy,j) = rbuf(k,jyl,j)
end if

end subroutine tr_unpack_out_2_kernel

end module tridag_transpose_gpu_mod
#endif

#ifdef PPMPI
!*******************************************************************************
subroutine tridag_array (a, b, c, r, u, p_halo_req, p_halo_nreq)
!*******************************************************************************
use types, only : rprec
use param
use cuda_mpi_debug, only : mpi_dbg_send_r, mpi_dbg_recv_r
use mpi
#ifdef ENABLE_CUDA
use cudafor
#endif
implicit none

#ifdef ENABLE_CUDA
real(rprec), managed, dimension(lh,ny,nz+1), intent(in) :: a, b, c

!  u and r are interleaved as complex arrays
real(rprec), managed, dimension(ld,ny,nz+1), intent(in) :: r
real(rprec), managed, dimension(ld,ny,nz+1), intent(out) :: u
#else
real(rprec), dimension(lh,ny,nz+1), intent(in) :: a, b, c

!  u and r are interleaved as complex arrays
real(rprec), dimension(ld,ny,nz+1), intent(in) :: r
real(rprec), dimension(ld,ny,nz+1), intent(out) :: u
#endif
integer, intent(out) :: p_halo_req(2)
integer, intent(out) :: p_halo_nreq

integer :: n, nchunks
integer :: chunksize
integer :: cstart, cend
integer :: jx, jy, j, j_min, j_max
integer :: tag0
integer :: q
integer :: ir, ii
#ifdef ENABLE_CUDA
integer :: cuda_istat
logical, save :: tridag_cuda_checked = .false.
logical, save :: tridag_cuda_enabled = .true.
logical, save :: tridag_fused_enabled = .true.
logical, save :: tridag_const_cache_enabled = .true.
logical, save :: tridag_min_sync_enabled = .true.
logical, save :: tridag_local_cache_enabled = .true.
logical, save :: tridag_p_halo_enabled = .true.
! Replicated global solve is retained as a validation/benchmark variant.  It is
! not the default production pressure path because it gathers each full global
! z-system on every rank before solving.
logical, save :: tridag_replicated_enabled = .false.
logical, save :: tridag_transpose_enabled = .true.
! Default-on for validated two-rank pressure solves; env=0 restores chasing.
logical, save :: tridag_spike2_enabled = .true.
logical, save :: tridag_const_cache_ready = .false.
logical, save :: tridag_local_cache_ready = .false.
logical, save :: tridag_cuda_allocated = .false.
integer :: tridag_use_const_cache
integer :: tridag_use_local_cache
integer :: tridag_back_start
real(rprec), managed, save, allocatable, dimension(:,:) :: bet_gpu
real(rprec), managed, save, allocatable, dimension(:,:) :: bet_in_gpu
real(rprec), managed, save, allocatable, dimension(:,:) :: c_in_gpu
real(rprec), managed, save, allocatable, dimension(:,:) :: gam_top_gpu
real(rprec), managed, save, allocatable, dimension(:,:,:) :: inv_bet_gpu
real(rprec), managed, save, allocatable, dimension(:,:,:) :: gam_gpu
#endif

real(rprec) :: bet(lh, ny)
real(rprec), dimension(lh,ny,nz+1) :: gam

! Initialize variables
n = nz + 1
nchunks = ny

! make sure nchunks divides ny evenly
chunksize = ny / nchunks

#ifdef ENABLE_CUDA
if (.not. tridag_cuda_checked) then
    tridag_cuda_enabled = .true.
    tridag_fused_enabled = .true.
    tridag_const_cache_enabled = .true.
    tridag_min_sync_enabled = .true.
    tridag_local_cache_enabled = .true.
    tridag_p_halo_enabled = .true.
    tridag_replicated_enabled = .false.
    tridag_transpose_enabled = .true.
    tridag_spike2_enabled = .false.
    tridag_cuda_checked = .true.
end if

if (tridag_cuda_enabled) then
    if (tridag_transpose_enabled .and. nproc == 2) then
        call tridag_array_transpose_thomas_cuda(a, b, c, r, u,                 &
            p_halo_req, p_halo_nreq)
        return
    end if

    if (tridag_spike2_enabled .and. nproc == 2) then
        call tridag_array_spike2_cuda(a, b, c, r, u, p_halo_req, p_halo_nreq)
        return
    end if

    if (tridag_replicated_enabled .and. nproc > 1) then
        call tridag_array_replicated_cuda(a, b, c, r, u, p_halo_req,           &
            p_halo_nreq)
        return
    end if

    p_halo_req = MPI_REQUEST_NULL
    p_halo_nreq = 0
    tridag_back_start = n - 1

    if (.not. tridag_cuda_allocated) then
        allocate(bet_gpu(lh, ny), bet_in_gpu(lh, ny), c_in_gpu(lh, ny),        &
            gam_top_gpu(lh, ny), inv_bet_gpu(lh, ny, nz+1),                   &
            gam_gpu(lh, ny, nz+1))
        tridag_cuda_allocated = .true.
    end if

    tridag_use_const_cache = 0
    if (tridag_const_cache_enabled .and. tridag_const_cache_ready)             &
        tridag_use_const_cache = 1
    tridag_use_local_cache = 0
    if (tridag_local_cache_enabled .and. tridag_local_cache_ready)             &
        tridag_use_local_cache = 1

    tag0 = 0
    if (coord == 0) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, lh-1
            ii = 2*jx
            ir = ii - 1
            u(ir,jy,1) = r(ir,jy,1) / b(jx,jy,1)
            u(ii,jy,1) = r(ii,jy,1) / b(jx,jy,1)
            bet_gpu(jx,jy) = b(jx,jy,1)
            if (tridag_use_local_cache == 0)                                  &
                inv_bet_gpu(jx,jy,1) = 1._rprec / b(jx,jy,1)
        end do
        end do
        if (.not. tridag_min_sync_enabled) then
            cuda_istat = cudaDeviceSynchronize()
            if (cuda_istat /= 0) stop 'tridag_array GPU init sync failed'
        end if
        cuda_istat = cudaGetLastError()
        if (cuda_istat /= 0) stop 'tridag_array GPU init kernel failed'
        j_min = 1
    else
        if (tridag_use_const_cache == 1) then
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                bet_gpu(jx,jy) = bet_in_gpu(jx,jy)
            end do
            end do
        else
            call mpi_dbg_recv_r(c(1, 1, 1), lh*ny, MPI_RPREC, down,           &
                           tag0+1, comm, status, ierr, 'tridag_fwd_c_recv')
            call mpi_dbg_recv_r(bet_gpu(1, 1), lh*ny, MPI_RPREC, down,        &
                           tag0+2, comm, status, ierr, 'tridag_fwd_bet_recv')
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                c_in_gpu(jx,jy) = c(jx,jy,1)
                bet_in_gpu(jx,jy) = bet_gpu(jx,jy)
            end do
            end do
        end if
        call mpi_dbg_recv_r(u(1,1,1), ld*ny, MPI_RPREC, down,                 &
                       tag0+3, comm, status, ierr, 'tridag_fwd_u_recv')
        j_min = 2
    end if

    if (coord == nproc-1) then
        j_max = n
    else
        j_max = n-1
    end if

    if (tridag_fused_enabled) then
        if (tridag_use_local_cache == 1) then
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
                    ii = 2*jx
                    ir = ii - 1
                    do j = 2, j_max
                        u(ir, jy, j) = (r(ir, jy, j) - a(jx, jy, j) *         &
                            u(ir, jy, j-1)) * inv_bet_gpu(jx, jy, j)
                        u(ii, jy, j) = (r(ii, jy, j) - a(jx, jy, j) *         &
                            u(ii, jy, j-1)) * inv_bet_gpu(jx, jy, j)
                    end do
                end if
            end do
            end do
        else
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
                    ii = 2*jx
                    ir = ii - 1
                    do j = 2, j_max
                        if (j == 2 .and. coord /= 0 .and.                     &
                            tridag_use_const_cache == 1) then
                            gam_gpu(jx, jy, j) = c_in_gpu(jx, jy) /           &
                                bet_gpu(jx,jy)
                        else
                            gam_gpu(jx, jy, j) = c(jx, jy, j-1) /             &
                                bet_gpu(jx, jy)
                        end if
                        bet_gpu(jx, jy) = b(jx, jy, j) - a(jx, jy, j) *       &
                            gam_gpu(jx, jy, j)
                        inv_bet_gpu(jx, jy, j) = 1._rprec / bet_gpu(jx, jy)
                        u(ir, jy, j) = (r(ir, jy, j) - a(jx, jy, j) *         &
                            u(ir, jy, j-1)) * inv_bet_gpu(jx, jy, j)
                        u(ii, jy, j) = (r(ii, jy, j) - a(jx, jy, j) *         &
                            u(ii, jy, j-1)) * inv_bet_gpu(jx, jy, j)
                    end do
                end if
            end do
            end do
        end if
    else
        do j = 2, j_max
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
                    ii = 2*jx
                    ir = ii - 1
                    if (j == 2 .and. coord /= 0 .and.                         &
                        tridag_use_const_cache == 1) then
                        gam_gpu(jx, jy, j) = c_in_gpu(jx, jy) / bet_gpu(jx,jy)
                    else
                        gam_gpu(jx, jy, j) = c(jx, jy, j-1) / bet_gpu(jx, jy)
                    end if
                    bet_gpu(jx, jy) = b(jx, jy, j) - a(jx, jy, j) *           &
                        gam_gpu(jx, jy, j)
                    inv_bet_gpu(jx, jy, j) = 1._rprec / bet_gpu(jx, jy)
                    u(ir, jy, j) = (r(ir, jy, j) - a(jx, jy, j) *             &
                        u(ir, jy, j-1)) * inv_bet_gpu(jx, jy, j)
                    u(ii, jy, j) = (r(ii, jy, j) - a(jx, jy, j) *             &
                        u(ii, jy, j-1)) * inv_bet_gpu(jx, jy, j)
                end if
            end do
            end do
        end do
    end if
    if ((.not. tridag_min_sync_enabled) .or. (coord /= nproc - 1)) then
        cuda_istat = cudaDeviceSynchronize()
        if (cuda_istat /= 0) stop 'tridag_array GPU forward sync failed'
    end if
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'tridag_array GPU forward kernel failed'

    if (coord /= nproc - 1) then
        if (tridag_use_const_cache == 0) then
            call mpi_dbg_send_r (c(1, 1, n-1), lh*ny, MPI_RPREC, up,          &
                           tag0+1, comm, ierr, 'tridag_fwd_c_send')
            call mpi_dbg_send_r (bet_gpu(1, 1), lh*ny, MPI_RPREC, up,         &
                           tag0+2, comm, ierr, 'tridag_fwd_bet_send')
        end if
        call mpi_dbg_send_r (u(1, 1, n-1), ld*ny, MPI_RPREC, up,              &
                       tag0+3, comm, ierr, 'tridag_fwd_u_send')
    end if

    if (tridag_p_halo_enabled .and. coord /= 0) then
        p_halo_nreq = p_halo_nreq + 1
        call mpi_irecv(u(1, 1, 1), ld*ny, MPI_RPREC, down, 2, comm,           &
            p_halo_req(p_halo_nreq), ierr)
    end if

    if (coord /= nproc - 1) then
        call mpi_dbg_recv_r (u(1, 1, n), ld*ny, MPI_RPREC, up,                &
                       tag0+4, comm, status, ierr, 'tridag_back_u_recv')
        if (tridag_use_const_cache == 1) then
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                gam_gpu(jx,jy,n) = gam_top_gpu(jx,jy)
            end do
            end do
        else
            call mpi_dbg_recv_r (gam_gpu(1, 1, n), lh*ny, MPI_RPREC, up,      &
                           tag0+5, comm, status, ierr, 'tridag_back_gam_recv')
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                gam_top_gpu(jx,jy) = gam_gpu(jx,jy,n)
            end do
            end do
        end if
    end if

    if (tridag_p_halo_enabled .and. coord /= nproc - 1) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, lh-1
            if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
                ii = 2*jx
                ir = ii - 1
                u(ir, jy, n-1) = u(ir, jy, n-1) - gam_gpu(jx, jy, n) *       &
                    u(ir, jy, n)
                u(ii, jy, n-1) = u(ii, jy, n-1) - gam_gpu(jx, jy, n) *       &
                    u(ii, jy, n)
            end if
        end do
        end do
        cuda_istat = cudaDeviceSynchronize()
        if (cuda_istat /= 0) stop 'tridag_array GPU p halo sync failed'
        cuda_istat = cudaGetLastError()
        if (cuda_istat /= 0) stop 'tridag_array GPU p halo kernel failed'
        p_halo_nreq = p_halo_nreq + 1
        call mpi_isend(u(1, 1, n-1), ld*ny, MPI_RPREC, up, 2, comm,          &
            p_halo_req(p_halo_nreq), ierr)
        tridag_back_start = n - 2
    end if

    if (tridag_fused_enabled) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, lh-1
            if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
                ii = 2*jx
                ir = ii - 1
                do j = tridag_back_start, j_min, -1
                    u(ir, jy, j) = u(ir, jy, j) - gam_gpu(jx, jy, j+1) *      &
                        u(ir, jy, j+1)
                    u(ii, jy, j) = u(ii, jy, j) - gam_gpu(jx, jy, j+1) *      &
                        u(ii, jy, j+1)
                end do
            end if
        end do
        end do
    else
        do j = tridag_back_start, j_min, -1
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, lh-1
                if ((jy /= ny/2+1) .and. (jx*jy /= 1)) then
                    ii = 2*jx
                    ir = ii - 1
                    u(ir, jy, j) = u(ir, jy, j) - gam_gpu(jx, jy, j+1) *      &
                        u(ir, jy, j+1)
                    u(ii, jy, j) = u(ii, jy, j) - gam_gpu(jx, jy, j+1) *      &
                        u(ii, jy, j+1)
                end if
            end do
            end do
        end do
    end if
    if ((.not. tridag_min_sync_enabled) .or. (coord /= 0)) then
        cuda_istat = cudaDeviceSynchronize()
        if (cuda_istat /= 0) stop 'tridag_array GPU backward sync failed'
    end if
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'tridag_array GPU backward kernel failed'

    if (coord /= 0) then
        call mpi_dbg_send_r (u(1, 1, 2), ld*ny, MPI_RPREC, down,              &
                       tag0+4, comm, ierr, 'tridag_back_u_send')
        if (tridag_use_const_cache == 0) then
            call mpi_dbg_send_r (gam_gpu(1, 1, 2), lh*ny, MPI_RPREC, down,    &
                       tag0+5, comm, ierr, 'tridag_back_gam_send')
        end if
    end if

    if (tridag_const_cache_enabled .and. .not. tridag_const_cache_ready)       &
        tridag_const_cache_ready = .true.
    if (tridag_local_cache_enabled .and. .not. tridag_local_cache_ready)       &
        tridag_local_cache_ready = .true.

    return
end if
#endif

p_halo_req = MPI_REQUEST_NULL
p_halo_nreq = 0

if (coord == 0) then
    do jy = 1, ny
        do jx = 1, lh-1
#ifdef PPSAFETYMODE
            if (b(jx, jy, 1) == 0._rprec) then
                write (*, *) 'tridag_array: rewrite eqs, jx, jy= ', jx, jy
                stop
            end if
#endif
            ii = 2*jx
            ir = ii - 1
            u(ir:ii,jy,1) = r(ir:ii,jy,1) / b(jx,jy,1)
        end do
    end do
    bet = b(:, :, 1)
    j_min = 1  ! this is only for backward pass
else
    j_min = 2  ! this is only for backward pass
end if

if (coord == nproc-1) then
    j_max = n
else
    j_max = n-1
end if

do q = 1, nchunks
    cstart = 1 + (q - 1) * chunksize
    cend = cstart + chunksize - 1
    tag0 = 0 + 10 * (q - 1)

    if (coord /= 0) then
        ! wait for c(:,:,1), bet(:,:), u(:,:,1) from "down"
        call mpi_dbg_recv_r(c(1, cstart, 1), lh*chunksize, MPI_RPREC, down,    &
                       tag0+1, comm, status, ierr, 'tridag_fwd_c_recv')
        call mpi_dbg_recv_r(bet(1, cstart), lh*chunksize, MPI_RPREC, down,     &
                       tag0+2, comm, status, ierr, 'tridag_fwd_bet_recv')
        call mpi_dbg_recv_r(u(1,cstart,1), ld*chunksize, MPI_RPREC, down,      &
                       tag0+3, comm, status, ierr, 'tridag_fwd_u_recv')
    end if

    do j = 2, j_max
        do jy = cstart, cend
            if (jy == ny/2+1) cycle
            do jx = 1, lh-1
                if (jx*jy == 1) cycle
                gam(jx, jy, j) = c(jx, jy, j-1) / bet(jx, jy)
                bet(jx, jy) = b(jx, jy, j) - a(jx, jy, j)*gam(jx, jy, j)

#ifdef PPSAFETYMODE
                if (bet(jx, jy) == 0._rprec) then
                    write (*, *) 'tridag_array failed at jx,jy,j=', jx, jy, j
                    write (*, *) 'a,b,c,gam,bet=', a(jx, jy, j), b(jx, jy, j), &
                        c(jx, jy, j), gam(jx, jy, j), bet(jx, jy)
                    stop
                end if
#endif
                ii = 2*jx
                ir = ii - 1
                u(ir:ii, jy, j) = (r(ir:ii, jy, j) - a(jx, jy, j) *            &
                u(ir:ii, jy, j-1)) / bet(jx, jy)
            end do
        end do
    end do

    if (coord /= nproc - 1) then
        ! send c(n-1), bet, u(n-1) to "up"
        call mpi_dbg_send_r (c(1, cstart, n-1), lh*chunksize, MPI_RPREC, up,   &
                       tag0+1, comm, ierr, 'tridag_fwd_c_send')
        call mpi_dbg_send_r (bet(1, cstart), lh*chunksize, MPI_RPREC, up,      &
                       tag0+2, comm, ierr, 'tridag_fwd_bet_send')
        call mpi_dbg_send_r (u(1, cstart, n-1), ld*chunksize, MPI_RPREC, up,   &
                       tag0+3, comm, ierr, 'tridag_fwd_u_send')
    end if
end do

do q = 1, nchunks
    cstart = 1 + (q - 1) * chunksize
    cend = cstart + chunksize - 1
    tag0 = 0 + 10 * (q - 1)

    if (coord /= nproc - 1) then
        ! wait for u(n), gam(n) from "up"
        call mpi_dbg_recv_r (u(1, cstart, n), ld*chunksize, MPI_RPREC, up,     &
                       tag0+4, comm, status, ierr, 'tridag_back_u_recv')
        call mpi_dbg_recv_r (gam(1, cstart, n), lh*chunksize, MPI_RPREC, up,   &
                       tag0+5, comm, status, ierr, 'tridag_back_gam_recv')
    end if

    do j = n-1, j_min, -1
        do jy = cstart, cend
            if (jy == ny/2+1) cycle
            do jx = 1, lh-1
                if (jx*jy == 1) cycle
                ii = 2*jx
                ir = ii - 1
                u(ir:ii, jy, j) = u(ir:ii, jy, j) - gam(jx, jy, j+1) *         &
                    u(ir:ii, jy, j+1)
            end do
        end do
    end do

    ! send u(2), gam(2) to "down"
    call mpi_dbg_send_r (u(1, cstart, 2), ld*chunksize, MPI_RPREC, down,       &
                   tag0+4, comm, ierr, 'tridag_back_u_send')
    call mpi_dbg_send_r (gam(1, cstart, 2), lh*chunksize, MPI_RPREC, down,     &
                   tag0+5, comm, ierr, 'tridag_back_gam_send')

end do

end subroutine tridag_array

!*******************************************************************************
subroutine tridag_apply_env_enabled_unless_false(name, enabled)
!*******************************************************************************
!
! Preserve pressure/tridiagonal switches where an unset variable keeps the
! current default, explicit false tokens disable the path, and any other set
! value enables it.
!
implicit none

character(len=*), intent(in) :: name
logical, intent(inout) :: enabled
character(len=32) :: setting
integer :: stat

call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('0', 'false', 'FALSE', 'False', 'off', 'OFF', 'Off', 'no', 'NO',    &
        'No')
        enabled = .false.
    case default
        enabled = .true.
    end select
end if

end subroutine tridag_apply_env_enabled_unless_false

!*******************************************************************************
subroutine tridag_apply_env_true_token(name, enabled)
!*******************************************************************************
!
! Preserve benchmark switches where the variable only takes effect when set;
! once set, only explicit true tokens enable the option.
!
implicit none

character(len=*), intent(in) :: name
logical, intent(inout) :: enabled
character(len=32) :: setting
integer :: stat

call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('1', 'true', 'TRUE', 'True', 'on', 'ON', 'On', 'yes', 'YES', 'Yes')
        enabled = .true.
    case default
        enabled = .false.
    end select
end if

end subroutine tridag_apply_env_true_token

#ifdef ENABLE_CUDA
!*******************************************************************************
subroutine tridag_array_transpose_thomas_cuda(a, b, c, r, u,                  &
    p_halo_req, p_halo_nreq)
!*******************************************************************************
use types, only : rprec
use param
use mpi
use cudafor
use tridag_transpose_gpu_mod
implicit none

real(rprec), managed, dimension(lh,ny,nz+1), intent(in) :: a, b, c
real(rprec), managed, dimension(ld,ny,nz+1), intent(in) :: r
real(rprec), managed, dimension(ld,ny,nz+1), intent(out) :: u
integer, intent(out) :: p_halo_req(2)
integer, intent(out) :: p_halo_nreq

integer :: n, ng, nym
integer :: own_jlo, rem_jlo, rem_coord, rem_rank
integer :: send_jlo, send_planes, local_glo
integer :: rem_send_jlo, rem_glo
integer :: jx, jyl, jy, j, g, k, ir, ii, idx
integer :: total_in, total_out
integer :: threads, blocks
integer :: block_x, block_y, block_z
integer :: cuda_istat, env_stat
integer, save :: trans_cached_n = 0
integer, save :: trans_cached_ng = 0
integer, save :: trans_cached_nym = 0
integer :: send_stat, recv_stat, last, dev, dev_stat
type(dim3) :: block3, grid3
type(cudaPointerAttributes) :: send_attr, recv_attr
type(cudaEvent), save :: trans_evt_start, trans_evt_mid, trans_evt_stop
real(rprec) :: bet_val, gam_val, a_val, b_val, c_prev
real(rprec) :: t0, t1, t_total0
real(rprec) :: t_self_copy_in, t_pack_remote_in
real(rprec) :: t_pack_in, t_exchange_in, t_unpack_in, t_thomas
real(rprec) :: t_pack_out, t_exchange_out, t_unpack_out, t_total
real(rprec) :: t_direct_out, t_unpack_in_gpu, t_thomas_gpu
real(rprec) :: t_direct_out_gpu, t_pack_out_gpu, t_unpack_out_gpu
real(rprec) :: t_tmp, t_self_copy_rank, t_pack_remote_rank
real(rprec) :: t_exchange_in_rank, t_unpack_in_rank, t_thomas_rank
real(rprec) :: t_pack_out_rank, t_exchange_out_rank, t_unpack_out_rank
real(rprec) :: t_total_rank
real(rprec) :: t_pack_in_debt, t_ready_in, t_ready_rel_in, t_ready_min_in
real(rprec) :: t_ready_max_in, t_arrival_skew_in, t_barrier_in
real(rprec) :: t_pack_out_debt, t_ready_out, t_ready_rel_out
real(rprec) :: t_ready_min_out, t_ready_max_out, t_arrival_skew_out
real(rprec) :: t_barrier_out
real(rprec) :: pnorm_local, pnorm_global
real(rprec) :: pack_in_bytes, exchange_in_bytes, pack_out_bytes
real(rprec) :: exchange_out_bytes, gbs_self, gbs_pack, gbs_exchange
real(rprec) :: gbs_unpack
real :: elapsed_ms
character(len=64) :: mpich_gpu, visible, send_type_name, recv_type_name
integer :: env_len
logical, save :: trans_allocated = .false.
logical, save :: trans_coeff_ready = .false.
logical, save :: trans_printed = .false.
logical, save :: trans_timing_checked = .false.
logical, save :: trans_timing_enabled = .false.
logical, save :: trans_timing_header_printed = .false.
logical, save :: trans_timing2_header_printed = .false.
logical, save :: trans_direct_header_printed = .false.
logical, save :: trans_checksum_header_printed = .false.
logical, save :: trans_queue_header_printed = .false.
logical, save :: trans_debug_sync_enabled = .false.
logical, save :: trans_explicit_enabled = .false.
logical, save :: trans_generic_enabled = .false.
logical, save :: trans_direct_out_enabled = .true.
logical, save :: trans_checksum_enabled = .false.
logical, save :: trans_queue_attrib_enabled = .false.
logical, save :: trans_barrier_in_enabled = .false.
logical, save :: trans_barrier_out_enabled = .false.
logical, save :: trans_pointer_audit_printed = .false.
logical, save :: trans_events_ready = .false.
real(rprec), device, save, allocatable, dimension(:,:,:) :: trans_sol
real(rprec), device, save, allocatable, dimension(:,:,:) :: trans_send
real(rprec), device, save, allocatable, dimension(:,:,:) :: trans_recv
real(rprec), device, save, allocatable, dimension(:,:,:) :: trans_inv_bet
real(rprec), device, save, allocatable, dimension(:,:,:) :: trans_gam

n = nz + 1
ng = nz_tot + 1
p_halo_req = MPI_REQUEST_NULL
p_halo_nreq = 0

if (nproc /= 2) stop 'transpose Thomas pressure path requires nproc=2'
if (mod(ny, 2) /= 0) stop 'transpose Thomas pressure path requires even ny'

nym = ny / 2
rem_coord = 1 - coord
rem_rank = rank_of_coord(rem_coord)
if (coord == 0) then
    own_jlo = 1
    rem_jlo = nym + 1
    send_jlo = 1
else
    own_jlo = nym + 1
    rem_jlo = 1
    send_jlo = 2
end if
send_planes = nz
total_in = ld * nym * send_planes
total_out = ld * nym * n
local_glo = coord * (nz - 1) + send_jlo
if (rem_coord == 0) then
    rem_send_jlo = 1
else
    rem_send_jlo = 2
end if
rem_glo = rem_coord * (nz - 1) + rem_send_jlo

if (.not. trans_timing_checked) then
    call tridag_apply_env_enabled_unless_false(                               &
        'LESGO_PRESS_TRANSPOSE_TIMING', trans_timing_enabled)
    trans_debug_sync_enabled = .false.
    trans_explicit_enabled = .false.
    call tridag_apply_env_true_token('LESGO_PRESS_TRANSPOSE_GENERIC',         &
        trans_generic_enabled)
    call tridag_apply_env_true_token('LESGO_PRESS_DIRECT_THOMAS_OUT',         &
        trans_direct_out_enabled)
    trans_checksum_enabled = .false.
    trans_queue_attrib_enabled = .false.
    trans_barrier_in_enabled = .false.
    trans_barrier_out_enabled = .false.
    trans_timing_checked = .true.
end if

if (coord == 0 .and. .not. trans_printed) then
    write(*,*) 'Pressure tridiagonal: transpose-owner full-z Thomas ON',      &
        ' explicit kernels=', trans_explicit_enabled,                         &
        ' generic=', trans_generic_enabled,                                    &
        ' direct_out=', trans_direct_out_enabled
    trans_printed = .true.
end if

if (trans_allocated) then
    if ((trans_cached_n /= n) .or. (trans_cached_ng /= ng) .or.               &
        (trans_cached_nym /= nym)) then
        deallocate(trans_sol, trans_send, trans_recv, trans_inv_bet,          &
            trans_gam)
        trans_allocated = .false.
        trans_coeff_ready = .false.
    end if
end if

if (.not. trans_allocated) then
    allocate(trans_sol(ld,nym,ng), trans_send(ld,nym,n),                      &
        trans_recv(ld,nym,n), trans_inv_bet(lh,nym,ng),                       &
        trans_gam(lh,nym,ng))
    trans_cached_n = n
    trans_cached_ng = ng
    trans_cached_nym = nym
    trans_allocated = .true.
    trans_coeff_ready = .false.
end if

t_pack_in = 0._rprec
t_self_copy_in = 0._rprec
t_pack_remote_in = 0._rprec
t_exchange_in = 0._rprec
t_unpack_in = 0._rprec
t_thomas = 0._rprec
t_direct_out = 0._rprec
t_pack_out = 0._rprec
t_exchange_out = 0._rprec
t_unpack_out = 0._rprec
t_unpack_in_gpu = 0._rprec
t_thomas_gpu = 0._rprec
t_direct_out_gpu = 0._rprec
t_pack_out_gpu = 0._rprec
t_unpack_out_gpu = 0._rprec
t_pack_in_debt = 0._rprec
t_ready_rel_in = 0._rprec
t_arrival_skew_in = 0._rprec
t_barrier_in = 0._rprec
t_pack_out_debt = 0._rprec
t_ready_rel_out = 0._rprec
t_arrival_skew_out = 0._rprec
t_barrier_out = 0._rprec
pack_in_bytes = real(2 * int(ld,8) * int(nym,8) * int(send_planes,8) * 8_8,  &
    rprec)
exchange_in_bytes = real(int(ld,8) * int(nym,8) * int(send_planes,8) * 8_8,  &
    rprec)
pack_out_bytes = real(2 * int(ld,8) * int(nym,8) * int(n,8) * 8_8, rprec)
exchange_out_bytes = real(int(ld,8) * int(nym,8) * int(n,8) * 8_8, rprec)
block_x = 32
block_y = 4
block_z = 2
block3 = dim3(block_x, block_y, block_z)
if ((trans_timing_enabled .or. trans_queue_attrib_enabled) .and.              &
    (.not. trans_events_ready)) then
    cuda_istat = cudaEventCreate(trans_evt_start)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event create failed'
    cuda_istat = cudaEventCreate(trans_evt_mid)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event create failed'
    cuda_istat = cudaEventCreate(trans_evt_stop)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event create failed'
    trans_events_ready = .true.
end if
t_total0 = mpi_wtime()

t0 = mpi_wtime()
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_start, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas unpack event start failed'
end if
if (trans_generic_enabled) then
    if (trans_explicit_enabled) then
        threads = 256
        blocks = (total_in + threads - 1) / threads
        call tr_pack_in_kernel<<<blocks,threads>>>(r, trans_sol, trans_send,  &
            ld, ny, nym, send_jlo, own_jlo, rem_jlo, local_glo, total_in)
    else
!$cuf kernel do(1) <<<*,*>>>
    do idx = 1, total_in
        k = mod(idx - 1, ld) + 1
        jyl = mod((idx - 1) / ld, nym) + 1
        j = (idx - 1) / (ld * nym) + 1
        trans_sol(k,jyl,local_glo + j - 1) = r(k,own_jlo + jyl - 1,          &
            send_jlo + j - 1)
        trans_send(k,jyl,j) = r(k,rem_jlo + jyl - 1, send_jlo + j - 1)
    end do
    end if
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas pack-in sync failed'
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'transpose Thomas pack-in kernel failed'
    t1 = mpi_wtime()
    t_pack_in = t1 - t0
else
    grid3 = dim3((ld + block_x - 1) / block_x,                               &
        (nym + block_y - 1) / block_y,                                       &
        (send_planes + block_z - 1) / block_z)
    if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
        cuda_istat = cudaEventRecord(trans_evt_start, 0)
        if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event start failed'
    end if
    call tr_self_copy_in_2_kernel<<<grid3,block3>>>(r, trans_sol, ld, ny,    &
        nym, send_jlo, own_jlo, local_glo, send_planes)
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'transpose Thomas self-copy-in kernel failed'
    if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
        cuda_istat = cudaEventRecord(trans_evt_mid, 0)
        if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event mid failed'
    end if
    call tr_pack_remote_in_2_kernel<<<grid3,block3>>>(r, trans_send, ld, ny, &
        nym, send_jlo, rem_jlo, send_planes)
    if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
        cuda_istat = cudaEventRecord(trans_evt_stop, 0)
        if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event stop failed'
    end if
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas pack-remote-in sync failed'
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'transpose Thomas pack-remote-in kernel failed'
    t1 = mpi_wtime()
    t_pack_in = t1 - t0
    if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
        cuda_istat = cudaEventElapsedTime(elapsed_ms, trans_evt_start,        &
            trans_evt_mid)
        if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event self failed'
        t_self_copy_in = real(elapsed_ms, rprec) * 1.0e-3_rprec
        cuda_istat = cudaEventElapsedTime(elapsed_ms, trans_evt_mid,          &
            trans_evt_stop)
        if (cuda_istat /= cudaSuccess) stop 'transpose Thomas event pack failed'
        t_pack_remote_in = real(elapsed_ms, rprec) * 1.0e-3_rprec
    else
        t_self_copy_in = 0._rprec
        t_pack_remote_in = t_pack_in
    end if
end if

t_pack_in_debt = max(0._rprec, t_pack_in - t_self_copy_in - t_pack_remote_in)
t_ready_in = mpi_wtime()
if (trans_queue_attrib_enabled .or. trans_barrier_in_enabled) then
    call mpi_allreduce(t_ready_in, t_ready_min_in, 1, MPI_RPREC, MPI_MIN,     &
        comm, ierr)
    call mpi_allreduce(t_ready_in, t_ready_max_in, 1, MPI_RPREC, MPI_MAX,     &
        comm, ierr)
    t_ready_rel_in = t_ready_in - t_ready_min_in
    t_arrival_skew_in = t_ready_max_in - t_ready_min_in
end if

if (.not. trans_pointer_audit_printed) then
    trans_pointer_audit_printed = .true.
    call get_environment_variable('MPICH_GPU_SUPPORT_ENABLED', mpich_gpu,     &
        length=env_len, status=env_stat)
    if (env_stat /= 0 .or. env_len <= 0) mpich_gpu = 'unset'
    call get_environment_variable('CUDA_VISIBLE_DEVICES', visible,            &
        length=env_len, status=env_stat)
    if (env_stat /= 0 .or. env_len <= 0) visible = 'unset'
    dev_stat = cudaGetDevice(dev)
    send_attr%type = -1
    send_attr%device = -1
    recv_attr%type = -1
    recv_attr%device = -1
    send_stat = cudaPointerGetAttributes(send_attr, trans_send)
    last = cudaGetLastError()
    recv_stat = cudaPointerGetAttributes(recv_attr, trans_recv)
    last = cudaGetLastError()
    select case (send_attr%type)
    case (1)
        send_type_name = 'host'
    case (2)
        send_type_name = 'device'
    case (3)
        send_type_name = 'managed'
    case default
        write(send_type_name,'(a,i0)') 'type', send_attr%type
    end select
    select case (recv_attr%type)
    case (1)
        recv_type_name = 'host'
    case (2)
        recv_type_name = 'device'
    case (3)
        recv_type_name = 'managed'
    case default
        write(recv_type_name,'(a,i0)') 'type', recv_attr%type
    end select
    write(*,*) 'PRESS_TRANSPOSE_POINTER_AUDIT rank=', coord,                 &
        ' path_generic=', trans_generic_enabled,                              &
        ' selected_gpu=', dev, ' cudaGetDevice_status=', dev_stat,           &
        ' MPICH_GPU_SUPPORT=', trim(mpich_gpu),                              &
        ' CUDA_VISIBLE_DEVICES=', trim(visible),                             &
        ' send_status=', send_stat, ' send_type=', send_attr%type,           &
        ' send_type_name=', trim(send_type_name),                            &
        ' send_device=', send_attr%device,                                   &
        ' recv_status=', recv_stat, ' recv_type=', recv_attr%type,           &
        ' recv_type_name=', trim(recv_type_name),                            &
        ' recv_device=', recv_attr%device
    flush(6)
end if

t0 = mpi_wtime()
if (trans_barrier_in_enabled) then
    call mpi_barrier(comm, ierr)
    t1 = mpi_wtime()
    t_barrier_in = t1 - t0
    t0 = t1
end if
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas pre exchange-in cuda error'
call mpi_sendrecv(trans_send(1,1,1), ld*nym*send_planes, MPI_RPREC,          &
    rem_rank, 231, trans_recv(1,1,1), ld*nym*send_planes, MPI_RPREC,         &
    rem_rank, 231, comm, status, ierr)
if (ierr /= 0) stop 'transpose Thomas exchange-in failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas post exchange-in cuda error'
t1 = mpi_wtime()
t_exchange_in = t1 - t0

t0 = mpi_wtime()
if (trans_generic_enabled) then
if (trans_explicit_enabled) then
    threads = 256
    blocks = (total_in + threads - 1) / threads
    call tr_unpack_in_kernel<<<blocks,threads>>>(trans_recv, trans_sol, ld,    &
        nym, rem_glo, total_in)
else
!$cuf kernel do(1) <<<*,*>>>
do idx = 1, total_in
    k = mod(idx - 1, ld) + 1
    jyl = mod((idx - 1) / ld, nym) + 1
    j = (idx - 1) / (ld * nym) + 1
    trans_sol(k,jyl,rem_glo + j - 1) = trans_recv(k,jyl,j)
end do
end if
else
    grid3 = dim3((ld + block_x - 1) / block_x,                               &
        (nym + block_y - 1) / block_y,                                       &
        (send_planes + block_z - 1) / block_z)
    call tr_unpack_in_2_kernel<<<grid3,block3>>>(trans_recv, trans_sol, ld,   &
        nym, rem_glo, send_planes)
end if
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_stop, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas unpack event stop failed'
end if
if (trans_timing_enabled .or. trans_debug_sync_enabled .or.                  &
    (.not. trans_generic_enabled)) then
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas unpack-in sync failed'
end if
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas unpack-in kernel failed'
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventElapsedTime(elapsed_ms, trans_evt_start,            &
        trans_evt_stop)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas unpack event elapsed failed'
    t_unpack_in_gpu = real(elapsed_ms, rprec) * 1.0e-3_rprec
end if
t1 = mpi_wtime()
t_unpack_in = t1 - t0

if (.not. trans_coeff_ready) then
!$cuf kernel do(2) <<<*,*>>>
    do jyl = 1, nym
    do jx = 1, lh - 1
        jy = own_jlo + jyl - 1
        if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
            bet_val = -1._rprec
            trans_gam(jx,jyl,1) = 0._rprec
            trans_inv_bet(jx,jyl,1) = 1._rprec / bet_val
            do g = 2, ng
                if (g == 2) then
                    c_prev = 1._rprec
                else
                    c_prev = c(jx,jy,2)
                end if
                if (g == ng) then
                    a_val = -1._rprec
                    b_val = 1._rprec
                else
                    a_val = a(jx,jy,2)
                    b_val = b(jx,jy,2)
                end if
                gam_val = c_prev / bet_val
                trans_gam(jx,jyl,g) = gam_val
                bet_val = b_val - a_val * gam_val
                trans_inv_bet(jx,jyl,g) = 1._rprec / bet_val
            end do
        end if
    end do
    end do
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas coeff sync failed'
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'transpose Thomas coeff kernel failed'
    trans_coeff_ready = .true.
end if

t0 = mpi_wtime()
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_start, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas solve event start failed'
end if
!$cuf kernel do(2) <<<*,*>>>
do jyl = 1, nym
do jx = 1, lh - 1
    jy = own_jlo + jyl - 1
    if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
        ii = 2*jx
        ir = ii - 1
        trans_sol(ir,jyl,1) = trans_sol(ir,jyl,1) *                           &
            trans_inv_bet(jx,jyl,1)
        trans_sol(ii,jyl,1) = trans_sol(ii,jyl,1) *                           &
            trans_inv_bet(jx,jyl,1)
        do g = 2, ng
            if (g == ng) then
                a_val = -1._rprec
            else
                a_val = a(jx,jy,2)
            end if
            trans_sol(ir,jyl,g) = (trans_sol(ir,jyl,g) - a_val *              &
                trans_sol(ir,jyl,g-1)) * trans_inv_bet(jx,jyl,g)
            trans_sol(ii,jyl,g) = (trans_sol(ii,jyl,g) - a_val *              &
                trans_sol(ii,jyl,g-1)) * trans_inv_bet(jx,jyl,g)
        end do
        if (trans_direct_out_enabled .and. (.not. trans_generic_enabled)) then
            j = ng - coord * (nz - 1)
            if ((j >= 1) .and. (j <= n)) then
                u(ir,jy,j) = trans_sol(ir,jyl,ng)
                u(ii,jy,j) = trans_sol(ii,jyl,ng)
            end if
            j = ng - rem_coord * (nz - 1)
            if ((j >= 1) .and. (j <= n)) then
                trans_send(ir,jyl,j) = trans_sol(ir,jyl,ng)
                trans_send(ii,jyl,j) = trans_sol(ii,jyl,ng)
            end if
        end if
        do g = ng - 1, 1, -1
            trans_sol(ir,jyl,g) = trans_sol(ir,jyl,g) -                       &
                trans_gam(jx,jyl,g+1) * trans_sol(ir,jyl,g+1)
            trans_sol(ii,jyl,g) = trans_sol(ii,jyl,g) -                       &
                trans_gam(jx,jyl,g+1) * trans_sol(ii,jyl,g+1)
            if (trans_direct_out_enabled .and.                                &
                (.not. trans_generic_enabled)) then
                j = g - coord * (nz - 1)
                if ((j >= 1) .and. (j <= n)) then
                    u(ir,jy,j) = trans_sol(ir,jyl,g)
                    u(ii,jy,j) = trans_sol(ii,jyl,g)
                end if
                j = g - rem_coord * (nz - 1)
                if ((j >= 1) .and. (j <= n)) then
                    trans_send(ir,jyl,j) = trans_sol(ir,jyl,g)
                    trans_send(ii,jyl,j) = trans_sol(ii,jyl,g)
                end if
            end if
        end do
    end if
end do
end do
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_stop, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas solve event stop failed'
end if
if (trans_timing_enabled .or. trans_debug_sync_enabled) then
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas solve sync failed'
end if
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas solve kernel failed'
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventElapsedTime(elapsed_ms, trans_evt_start,            &
        trans_evt_stop)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas solve event elapsed failed'
    t_thomas_gpu = real(elapsed_ms, rprec) * 1.0e-3_rprec
end if
t1 = mpi_wtime()
t_thomas = t1 - t0

t0 = mpi_wtime()
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_start, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas pack-out event start failed'
end if
if (trans_direct_out_enabled .and. (.not. trans_generic_enabled)) then
!$cuf kernel do(3) <<<*,*>>>
    do j = 1, n
    do jyl = 1, nym
    do k = ld - 1, ld
        jy = own_jlo + jyl - 1
        if (.not. ((jy == 1) .and. (k <= 2))) u(k,jy,j) = 0._rprec
    end do
    end do
    end do
    if ((own_jlo <= ny/2 + 1) .and.                                      &
        ((own_jlo + nym - 1) >= ny/2 + 1)) then
        jy = ny/2 + 1
!$cuf kernel do(2) <<<*,*>>>
        do j = 1, n
        do k = 1, ld
            if (.not. ((jy == 1) .and. (k <= 2))) u(k,jy,j) = 0._rprec
        end do
        end do
    end if
else
    if (trans_generic_enabled) then
    if (trans_explicit_enabled) then
        threads = 256
        blocks = (total_out + threads - 1) / threads
        call tr_pack_out_kernel<<<blocks,threads>>>(trans_sol, u, trans_send, &
            ld, ny, nym, n, nz, coord, own_jlo, rem_coord, total_out)
    else
!$cuf kernel do(1) <<<*,*>>>
    do idx = 1, total_out
        k = mod(idx - 1, ld) + 1
        jyl = mod((idx - 1) / ld, nym) + 1
        j = (idx - 1) / (ld * nym) + 1
        g = coord * (nz - 1) + j
        jy = own_jlo + jyl - 1
        if ((jy == 1) .and. (k <= 2)) then
            ! Zero mode is handled by the dedicated pressure zero-mode pass.
        else if ((k >= ld - 1) .or. (jy == ny/2 + 1)) then
            u(k,jy,j) = 0._rprec
        else
            u(k,jy,j) = trans_sol(k,jyl,g)
        end if
        g = rem_coord * (nz - 1) + j
        trans_send(k,jyl,j) = trans_sol(k,jyl,g)
    end do
    end if
    else
        grid3 = dim3((ld + block_x - 1) / block_x,                           &
            (nym + block_y - 1) / block_y, (n + block_z - 1) / block_z)
        call tr_pack_out_2_kernel<<<grid3,block3>>>(trans_sol, u, trans_send, &
            ld, ny, nym, n, nz, coord, own_jlo, rem_coord)
    end if
end if
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_stop, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas pack-out event stop failed'
end if
cuda_istat = cudaDeviceSynchronize()
if (cuda_istat /= 0) stop 'transpose Thomas pack-out sync failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas pack-out kernel failed'
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventElapsedTime(elapsed_ms, trans_evt_start,            &
        trans_evt_stop)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas pack-out event elapsed failed'
    t_pack_out_gpu = real(elapsed_ms, rprec) * 1.0e-3_rprec
end if
t1 = mpi_wtime()
if (trans_direct_out_enabled .and. (.not. trans_generic_enabled)) then
    t_direct_out = t1 - t0
    t_direct_out_gpu = t_pack_out_gpu
    t_pack_out = 0._rprec
    t_pack_out_gpu = 0._rprec
    t_pack_out_debt = t_direct_out
else
    t_pack_out = t1 - t0
    t_pack_out_debt = t_pack_out
end if
t_ready_out = mpi_wtime()
if (trans_queue_attrib_enabled .or. trans_barrier_out_enabled) then
    call mpi_allreduce(t_ready_out, t_ready_min_out, 1, MPI_RPREC, MPI_MIN,   &
        comm, ierr)
    call mpi_allreduce(t_ready_out, t_ready_max_out, 1, MPI_RPREC, MPI_MAX,   &
        comm, ierr)
    t_ready_rel_out = t_ready_out - t_ready_min_out
    t_arrival_skew_out = t_ready_max_out - t_ready_min_out
end if

t0 = mpi_wtime()
if (trans_barrier_out_enabled) then
    call mpi_barrier(comm, ierr)
    t1 = mpi_wtime()
    t_barrier_out = t1 - t0
    t0 = t1
end if
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas pre exchange-out cuda error'
call mpi_sendrecv(trans_send(1,1,1), ld*nym*n, MPI_RPREC, rem_rank, 232,      &
    trans_recv(1,1,1), ld*nym*n, MPI_RPREC, rem_rank, 232, comm, status,      &
    ierr)
if (ierr /= 0) stop 'transpose Thomas exchange-out failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas post exchange-out cuda error'
t1 = mpi_wtime()
t_exchange_out = t1 - t0

t0 = mpi_wtime()
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_start, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas unpack-out event start failed'
end if
if (trans_generic_enabled) then
if (trans_explicit_enabled) then
    threads = 256
    blocks = (total_out + threads - 1) / threads
    call tr_unpack_out_kernel<<<blocks,threads>>>(trans_recv, u, ld, ny,      &
        nym, rem_jlo, total_out)
else
!$cuf kernel do(1) <<<*,*>>>
do idx = 1, total_out
    k = mod(idx - 1, ld) + 1
    jyl = mod((idx - 1) / ld, nym) + 1
    j = (idx - 1) / (ld * nym) + 1
    jy = rem_jlo + jyl - 1
    if ((jy == 1) .and. (k <= 2)) then
        ! Zero mode is handled by the dedicated pressure zero-mode pass.
    else if ((k >= ld - 1) .or. (jy == ny/2 + 1)) then
        u(k,jy,j) = 0._rprec
    else
        u(k,jy,j) = trans_recv(k,jyl,j)
    end if
end do
end if
else
    grid3 = dim3((ld + block_x - 1) / block_x,                               &
        (nym + block_y - 1) / block_y, (n + block_z - 1) / block_z)
    call tr_unpack_out_2_kernel<<<grid3,block3>>>(trans_recv, u, ld, ny,      &
        nym, rem_jlo, n)
end if
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventRecord(trans_evt_stop, 0)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas unpack-out event stop failed'
end if
if (trans_timing_enabled .or. trans_debug_sync_enabled .or.                  &
    (.not. trans_generic_enabled)) then
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas unpack-out sync failed'
end if
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'transpose Thomas unpack-out kernel failed'
if (trans_timing_enabled .or. trans_queue_attrib_enabled) then
    cuda_istat = cudaEventElapsedTime(elapsed_ms, trans_evt_start,            &
        trans_evt_stop)
    if (cuda_istat /= cudaSuccess) stop 'transpose Thomas unpack-out event elapsed failed'
    t_unpack_out_gpu = real(elapsed_ms, rprec) * 1.0e-3_rprec
end if
t1 = mpi_wtime()
t_unpack_out = t1 - t0
t_total = t1 - t_total0

if (trans_checksum_enabled .and. (jt == nsteps)) then
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'transpose Thomas checksum sync failed'
    pnorm_local = 0._rprec
    do j = 1, n
    do jy = 1, ny
    do k = 1, ld
        pnorm_local = pnorm_local + u(k,jy,j) * u(k,jy,j)
    end do
    end do
    end do
    call mpi_allreduce(pnorm_local, pnorm_global, 1, MPI_RPREC, MPI_SUM,       &
        comm, ierr)
    if (coord == 0) then
        if (.not. trans_checksum_header_printed) then
            write(*,'(A)') 'PRESS_TRANSPOSE_P_NORM2 fields: jt direct_out norm2'
            trans_checksum_header_printed = .true.
        end if
        write(*,'(A,I8,1X,L1,1X,ES24.16)') 'PRESS_TRANSPOSE_P_NORM2 ', jt,    &
            (trans_direct_out_enabled .and. (.not. trans_generic_enabled)),    &
            pnorm_global
    end if
end if

if (trans_timing_enabled) then
    t_self_copy_rank = t_self_copy_in
    t_pack_remote_rank = t_pack_remote_in
    t_exchange_in_rank = t_exchange_in
    t_unpack_in_rank = t_unpack_in
    t_thomas_rank = t_thomas
    t_pack_out_rank = t_pack_out
    t_exchange_out_rank = t_exchange_out
    t_unpack_out_rank = t_unpack_out
    t_total_rank = t_total

    call mpi_allreduce(t_self_copy_in, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,    &
        ierr)
    t_self_copy_in = t_tmp
    call mpi_allreduce(t_pack_remote_in, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,  &
        ierr)
    t_pack_remote_in = t_tmp
    call mpi_allreduce(t_pack_in, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_pack_in = t_tmp
    call mpi_allreduce(t_exchange_in, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,     &
        ierr)
    t_exchange_in = t_tmp
    call mpi_allreduce(t_unpack_in, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_unpack_in = t_tmp
    call mpi_allreduce(t_thomas, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_thomas = t_tmp
    call mpi_allreduce(t_direct_out, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_direct_out = t_tmp
    call mpi_allreduce(t_pack_out, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_pack_out = t_tmp
    call mpi_allreduce(t_exchange_out, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,    &
        ierr)
    t_exchange_out = t_tmp
    call mpi_allreduce(t_unpack_out, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,      &
        ierr)
    t_unpack_out = t_tmp
    call mpi_allreduce(t_total, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_total = t_tmp
    call mpi_allreduce(t_unpack_in_gpu, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,   &
        ierr)
    t_unpack_in_gpu = t_tmp
    call mpi_allreduce(t_thomas_gpu, t_tmp, 1, MPI_RPREC, MPI_MAX, comm, ierr)
    t_thomas_gpu = t_tmp
    call mpi_allreduce(t_direct_out_gpu, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,   &
        ierr)
    t_direct_out_gpu = t_tmp
    call mpi_allreduce(t_pack_out_gpu, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,    &
        ierr)
    t_pack_out_gpu = t_tmp
    call mpi_allreduce(t_unpack_out_gpu, t_tmp, 1, MPI_RPREC, MPI_MAX, comm,  &
        ierr)
    t_unpack_out_gpu = t_tmp

    if (coord == 0) then
        if (.not. trans_timing_header_printed) then
            write(*,'(A)') 'PRESS_TRANSPOSE_THOMAS fields: pack_in ' //       &
                'exchange_in unpack_in Thomas pack_out exchange_out ' //      &
                'unpack_out total'
            trans_timing_header_printed = .true.
        end if
        write(*,'(A,8(1X,ES12.5))') 'PRESS_TRANSPOSE_THOMAS ',                &
            t_pack_in, t_exchange_in, t_unpack_in, t_thomas,                  &
            t_pack_out, t_exchange_out, t_unpack_out, t_total
        if (t_self_copy_in > 0._rprec)                                        &
            gbs_self = exchange_in_bytes / t_self_copy_in / 1.0e9_rprec
        if (t_self_copy_in <= 0._rprec) gbs_self = 0._rprec
        if (t_pack_remote_in > 0._rprec)                                      &
            gbs_pack = exchange_in_bytes / t_pack_remote_in / 1.0e9_rprec
        if (t_pack_remote_in <= 0._rprec) gbs_pack = 0._rprec
        if (t_exchange_in > 0._rprec)                                         &
            gbs_exchange = exchange_in_bytes / t_exchange_in / 1.0e9_rprec
        if (t_exchange_in <= 0._rprec) gbs_exchange = 0._rprec
        if (t_unpack_in > 0._rprec)                                           &
            gbs_unpack = exchange_in_bytes / t_unpack_in / 1.0e9_rprec
        if (t_unpack_in <= 0._rprec) gbs_unpack = 0._rprec
        if (.not. trans_timing2_header_printed) then
            write(*,'(A)') 'PRESS_TRANSPOSE_THOMAS2 fields: generic ' //      &
                'self_copy_in pack_remote_in exchange_in unpack_in Thomas '   &
                // 'pack_out exchange_out unpack_out pack_in_bytes '          &
                // 'exchange_in_bytes pack_out_bytes exchange_out_bytes '     &
                // 'self_gbs pack_gbs exchange_gbs unpack_gbs'
            write(*,'(A)') 'PRESS_TRANSPOSE_RANK fields: generic rank ' //    &
                'self_copy_in pack_remote_in exchange_in unpack_in Thomas '   &
                // 'pack_out exchange_out unpack_out pack_in_bytes '          &
                // 'exchange_in_bytes pack_out_bytes exchange_out_bytes '     &
                // 'self_gbs pack_gbs exchange_gbs unpack_gbs'
            trans_timing2_header_printed = .true.
        end if
        write(*,'(A,L1,8(1X,ES12.5),8(1X,ES12.5))')                           &
            'PRESS_TRANSPOSE_THOMAS2 generic=', trans_generic_enabled,        &
            t_self_copy_in, t_pack_remote_in, t_exchange_in, t_unpack_in,      &
            t_thomas, t_pack_out, t_exchange_out, t_unpack_out,               &
            pack_in_bytes, exchange_in_bytes, pack_out_bytes,                 &
            exchange_out_bytes, gbs_self, gbs_pack, gbs_exchange, gbs_unpack
        write(*,'(A)') 'PRESS_TRANSPOSE_HELPER_GPU fields: self_copy_in ' //  &
            'pack_remote_in unpack_in Thomas pack_out unpack_out'
        write(*,'(A,6(1X,ES12.5))') 'PRESS_TRANSPOSE_HELPER_GPU ',            &
            t_self_copy_in, t_pack_remote_in, t_unpack_in_gpu, t_thomas_gpu,  &
            t_pack_out_gpu, t_unpack_out_gpu
        if (.not. trans_direct_header_printed) then
            write(*,'(A)') 'PRESS_TRANSPOSE_DIRECT_OUT fields: enabled ' //    &
                'direct_out direct_out_gpu Thomas Thomas_gpu pack_out ' //     &
                'pack_out_gpu exchange_out unpack_out total'
            trans_direct_header_printed = .true.
        end if
        write(*,'(A,L1,9(1X,ES12.5))') 'PRESS_TRANSPOSE_DIRECT_OUT enabled=', &
            (trans_direct_out_enabled .and. (.not. trans_generic_enabled)),    &
            t_direct_out, t_direct_out_gpu, t_thomas, t_thomas_gpu,            &
            t_pack_out, t_pack_out_gpu, t_exchange_out, t_unpack_out, t_total
    end if
    if (t_self_copy_rank > 0._rprec)                                          &
        gbs_self = exchange_in_bytes / t_self_copy_rank / 1.0e9_rprec
    if (t_self_copy_rank <= 0._rprec) gbs_self = 0._rprec
    if (t_pack_remote_rank > 0._rprec)                                        &
        gbs_pack = exchange_in_bytes / t_pack_remote_rank / 1.0e9_rprec
    if (t_pack_remote_rank <= 0._rprec) gbs_pack = 0._rprec
    if (t_exchange_in_rank > 0._rprec)                                        &
        gbs_exchange = exchange_in_bytes / t_exchange_in_rank / 1.0e9_rprec
    if (t_exchange_in_rank <= 0._rprec) gbs_exchange = 0._rprec
    if (t_unpack_in_rank > 0._rprec)                                          &
        gbs_unpack = exchange_in_bytes / t_unpack_in_rank / 1.0e9_rprec
    if (t_unpack_in_rank <= 0._rprec) gbs_unpack = 0._rprec
    call mpi_barrier(comm, ierr)
    if (coord == 0) then
        write(*,'(A,L1,I4,8(1X,ES12.5),8(1X,ES12.5))')                        &
            'PRESS_TRANSPOSE_RANK generic=', trans_generic_enabled, coord,    &
            t_self_copy_rank, t_pack_remote_rank, t_exchange_in_rank,          &
            t_unpack_in_rank, t_thomas_rank, t_pack_out_rank,                  &
            t_exchange_out_rank, t_unpack_out_rank, pack_in_bytes,             &
            exchange_in_bytes, pack_out_bytes, exchange_out_bytes,             &
            gbs_self, gbs_pack, gbs_exchange, gbs_unpack
        flush(6)
    end if
    call mpi_barrier(comm, ierr)
    if (coord == 1) then
        write(*,'(A,L1,I4,8(1X,ES12.5),8(1X,ES12.5))')                        &
            'PRESS_TRANSPOSE_RANK generic=', trans_generic_enabled, coord,    &
            t_self_copy_rank, t_pack_remote_rank, t_exchange_in_rank,          &
            t_unpack_in_rank, t_thomas_rank, t_pack_out_rank,                  &
            t_exchange_out_rank, t_unpack_out_rank, pack_in_bytes,             &
            exchange_in_bytes, pack_out_bytes, exchange_out_bytes,             &
            gbs_self, gbs_pack, gbs_exchange, gbs_unpack
        flush(6)
    end if
    call mpi_barrier(comm, ierr)
end if

if (trans_queue_attrib_enabled) then
    if (coord == 0 .and. .not. trans_queue_header_printed) then
        write(*,'(A)') 'PRESS_EXCHANGE_ATTRIB fields: rank stage pre_sync_debt'&
            // ' ready_rel arrival_skew mpi_wait barrier bytes'
        trans_queue_header_printed = .true.
    end if
    call mpi_barrier(comm, ierr)
    if (coord == 0) then
        write(*,'(A,I4,1X,A,6(1X,ES12.5))') 'PRESS_EXCHANGE_ATTRIB ', coord, &
            'in ', t_pack_in_debt, t_ready_rel_in, t_arrival_skew_in,         &
            t_exchange_in, t_barrier_in, exchange_in_bytes
        write(*,'(A,I4,1X,A,6(1X,ES12.5))') 'PRESS_EXCHANGE_ATTRIB ', coord, &
            'out', t_pack_out_debt, t_ready_rel_out, t_arrival_skew_out,      &
            t_exchange_out, t_barrier_out, exchange_out_bytes
        flush(6)
    end if
    call mpi_barrier(comm, ierr)
    if (coord == 1) then
        write(*,'(A,I4,1X,A,6(1X,ES12.5))') 'PRESS_EXCHANGE_ATTRIB ', coord, &
            'in ', t_pack_in_debt, t_ready_rel_in, t_arrival_skew_in,         &
            t_exchange_in, t_barrier_in, exchange_in_bytes
        write(*,'(A,I4,1X,A,6(1X,ES12.5))') 'PRESS_EXCHANGE_ATTRIB ', coord, &
            'out', t_pack_out_debt, t_ready_rel_out, t_arrival_skew_out,      &
            t_exchange_out, t_barrier_out, exchange_out_bytes
        flush(6)
    end if
    call mpi_barrier(comm, ierr)
end if

end subroutine tridag_array_transpose_thomas_cuda

!*******************************************************************************
subroutine tridag_array_spike2_cuda(a, b, c, r, u, p_halo_req, p_halo_nreq)
!*******************************************************************************
use types, only : rprec
use param
use mpi
use cudafor
implicit none

real(rprec), managed, dimension(lh,ny,nz+1), intent(in) :: a, b, c
real(rprec), managed, dimension(ld,ny,nz+1), intent(in) :: r
real(rprec), managed, dimension(ld,ny,nz+1), intent(out) :: u
integer, intent(out) :: p_halo_req(2)
integer, intent(out) :: p_halo_nreq

integer :: n, jlo, jhi, jiface
integer :: jx, jy, j, ir, ii
integer :: cuda_istat
real(rprec) :: bet_val, z_rhs, denom
real(rprec) :: y0r, y0i, y1r, y1i, z0, z1, x0r, x0i, x1r, x1i
logical, save :: spike_allocated = .false.
logical, save :: spike_coeff_ready = .false.
logical, save :: spike_printed = .false.
integer, save :: spike_cached_n = 0
real(rprec), managed, save, allocatable, dimension(:,:,:) :: spike_inv_bet
real(rprec), managed, save, allocatable, dimension(:,:,:) :: spike_gam
real(rprec), managed, save, allocatable, dimension(:,:,:) :: spike_z
real(rprec), managed, save, allocatable, dimension(:,:) :: spike_y_recv
real(rprec), managed, save, allocatable, dimension(:,:) :: spike_z_recv

n = nz + 1
p_halo_req = MPI_REQUEST_NULL
p_halo_nreq = 0

if (coord == 0 .and. .not. spike_printed) then
    write(*,*) 'Pressure tridiagonal: 2-rank SPIKE solve ON'
    spike_printed = .true.
end if

if (coord == 0) then
    jlo = 1
else
    jlo = 2
end if
if (coord == nproc - 1) then
    jhi = n
else
    jhi = n - 1
end if
if (coord == 0) then
    jiface = jhi
else
    jiface = jlo
end if

if (spike_allocated) then
    if (spike_cached_n /= n) then
        deallocate(spike_inv_bet, spike_gam, spike_z, spike_y_recv,           &
            spike_z_recv)
        spike_allocated = .false.
        spike_coeff_ready = .false.
    end if
end if

if (.not. spike_allocated) then
    allocate(spike_inv_bet(lh,ny,n), spike_gam(lh,ny,n), spike_z(lh,ny,n),    &
        spike_y_recv(ld,ny), spike_z_recv(lh,ny))
    spike_cached_n = n
    spike_allocated = .true.
    spike_coeff_ready = .false.
end if

if (.not. spike_coeff_ready) then
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, lh - 1
        if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
            bet_val = b(jx,jy,jlo)
            spike_gam(jx,jy,jlo) = 0._rprec
            spike_inv_bet(jx,jy,jlo) = 1._rprec / bet_val
            do j = jlo + 1, jhi
                spike_gam(jx,jy,j) = c(jx,jy,j-1) / bet_val
                bet_val = b(jx,jy,j) - a(jx,jy,j) * spike_gam(jx,jy,j)
                spike_inv_bet(jx,jy,j) = 1._rprec / bet_val
            end do

            do j = jlo, jhi
                if ((coord == 0 .and. j == jhi) .or.                          &
                    (coord == nproc - 1 .and. j == jlo)) then
                    if (coord == 0) then
                        z_rhs = c(jx,jy,jhi)
                    else
                        z_rhs = a(jx,jy,jlo)
                    end if
                else
                    z_rhs = 0._rprec
                end if
                if (j == jlo) then
                    spike_z(jx,jy,j) = z_rhs * spike_inv_bet(jx,jy,j)
                else
                    spike_z(jx,jy,j) = (z_rhs - a(jx,jy,j) *                  &
                        spike_z(jx,jy,j-1)) * spike_inv_bet(jx,jy,j)
                end if
            end do
            do j = jhi - 1, jlo, -1
                spike_z(jx,jy,j) = spike_z(jx,jy,j) -                         &
                    spike_gam(jx,jy,j+1) * spike_z(jx,jy,j+1)
            end do
        end if
    end do
    end do
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'spike2 tridag coeff sync failed'
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'spike2 tridag coeff kernel failed'
    call mpi_sendrecv(spike_z(1,1,jiface), lh*ny, MPI_RPREC, 1 - coord, 22,   &
        spike_z_recv(1,1), lh*ny, MPI_RPREC, 1 - coord, 22, comm, status,     &
        ierr)
    if (ierr /= 0) stop 'spike2 tridag z interface exchange failed'
    spike_coeff_ready = .true.
end if

!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, lh - 1
    if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
        ii = 2*jx
        ir = ii - 1
        u(ir,jy,jlo) = r(ir,jy,jlo) * spike_inv_bet(jx,jy,jlo)
        u(ii,jy,jlo) = r(ii,jy,jlo) * spike_inv_bet(jx,jy,jlo)
        do j = jlo + 1, jhi
            u(ir,jy,j) = (r(ir,jy,j) - a(jx,jy,j) * u(ir,jy,j-1)) *           &
                spike_inv_bet(jx,jy,j)
            u(ii,jy,j) = (r(ii,jy,j) - a(jx,jy,j) * u(ii,jy,j-1)) *           &
                spike_inv_bet(jx,jy,j)
        end do
        do j = jhi - 1, jlo, -1
            u(ir,jy,j) = u(ir,jy,j) - spike_gam(jx,jy,j+1) * u(ir,jy,j+1)
            u(ii,jy,j) = u(ii,jy,j) - spike_gam(jx,jy,j+1) * u(ii,jy,j+1)
        end do
    end if
end do
end do

cuda_istat = cudaDeviceSynchronize()
if (cuda_istat /= 0) stop 'spike2 tridag local solve sync failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'spike2 tridag local solve kernel failed'

call mpi_sendrecv(u(1,1,jiface), ld*ny, MPI_RPREC, 1 - coord, 21,             &
    spike_y_recv(1,1), ld*ny, MPI_RPREC, 1 - coord, 21, comm, status, ierr)
if (ierr /= 0) stop 'spike2 tridag y interface exchange failed'

!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, lh - 1
    ii = 2*jx
    ir = ii - 1
    if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
        if (coord == 0) then
            y0r = u(ir,jy,jhi)
            y0i = u(ii,jy,jhi)
            z0 = spike_z(jx,jy,jhi)
            y1r = spike_y_recv(ir,jy)
            y1i = spike_y_recv(ii,jy)
            z1 = spike_z_recv(jx,jy)
        else
            y0r = spike_y_recv(ir,jy)
            y0i = spike_y_recv(ii,jy)
            z0 = spike_z_recv(jx,jy)
            y1r = u(ir,jy,jlo)
            y1i = u(ii,jy,jlo)
            z1 = spike_z(jx,jy,jlo)
        end if
        denom = 1._rprec - z0 * z1
        x0r = (y0r - z0 * y1r) / denom
        x0i = (y0i - z0 * y1i) / denom
        x1r = y1r - z1 * x0r
        x1i = y1i - z1 * x0i
        if (coord == 0) then
            do j = jlo, jhi
                u(ir,jy,j) = u(ir,jy,j) - spike_z(jx,jy,j) * x1r
                u(ii,jy,j) = u(ii,jy,j) - spike_z(jx,jy,j) * x1i
            end do
            u(ir,jy,n) = x1r
            u(ii,jy,n) = x1i
        else
            u(ir,jy,1) = x0r
            u(ii,jy,1) = x0i
            do j = jlo, jhi
                u(ir,jy,j) = u(ir,jy,j) - spike_z(jx,jy,j) * x0r
                u(ii,jy,j) = u(ii,jy,j) - spike_z(jx,jy,j) * x0i
            end do
        end if
    else
        do j = 1, n
            u(ir,jy,j) = 0._rprec
            u(ii,jy,j) = 0._rprec
        end do
    end if
end do
end do
cuda_istat = cudaDeviceSynchronize()
if (cuda_istat /= 0) stop 'spike2 tridag correction sync failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'spike2 tridag correction kernel failed'

end subroutine tridag_array_spike2_cuda

!*******************************************************************************
subroutine tridag_array_replicated_cuda(a, b, c, r, u, p_halo_req,             &
    p_halo_nreq)
!*******************************************************************************
use types, only : rprec
use param
use mpi
use cudafor
implicit none

real(rprec), managed, dimension(lh,ny,nz+1), intent(in) :: a, b, c
real(rprec), managed, dimension(ld,ny,nz+1), intent(in) :: r
real(rprec), managed, dimension(ld,ny,nz+1), intent(out) :: u
integer, intent(out) :: p_halo_req(2)
integer, intent(out) :: p_halo_nreq

integer :: n, ng
integer :: send_jlo, send_planes, send_count
integer :: ip, ip_rank, ip_jlo, ip_planes, ip_glo
integer :: jx, jy, j, g, ir, ii
integer :: cuda_istat
real(rprec) :: bet_val, gam_val, a_val, b_val, c_prev
logical, save :: rep_allocated = .false.
logical, save :: rep_coeff_ready = .false.
logical, save :: rep_printed = .false.
integer, save :: rep_cached_ng = 0
integer, save :: rep_cached_nproc = 0
integer, save :: rep_cached_nz = 0
real(rprec), managed, save, allocatable, dimension(:,:,:) :: rep_sol
real(rprec), managed, save, allocatable, dimension(:,:,:) :: rep_inv_bet
real(rprec), managed, save, allocatable, dimension(:,:,:) :: rep_gam
integer, save, allocatable, dimension(:) :: rep_counts, rep_displs

n = nz + 1
ng = nz_tot + 1
p_halo_req = MPI_REQUEST_NULL
p_halo_nreq = 0

if (coord == 0 .and. .not. rep_printed) then
    write(*,*) 'Pressure tridiagonal: validation replicated global solve ON'
    rep_printed = .true.
end if

if (rep_allocated) then
    if ((rep_cached_ng /= ng) .or. (rep_cached_nproc /= nproc) .or.            &
        (rep_cached_nz /= nz)) then
        deallocate(rep_sol, rep_inv_bet, rep_gam, rep_counts, rep_displs)
        rep_allocated = .false.
        rep_coeff_ready = .false.
    end if
end if

if (.not. rep_allocated) then
    allocate(rep_sol(ld, ny, ng), rep_inv_bet(lh, ny, ng),                     &
        rep_gam(lh, ny, ng), rep_counts(nproc), rep_displs(nproc))
    rep_cached_ng = ng
    rep_cached_nproc = nproc
    rep_cached_nz = nz
    rep_allocated = .true.
    rep_coeff_ready = .false.
end if

rep_counts = 0
rep_displs = 0
do ip = 0, nproc - 1
    if (ip == 0) then
        ip_jlo = 1
    else
        ip_jlo = 2
    end if
    if ((ip == 0) .or. (ip == nproc - 1)) then
        ip_planes = nz
    else
        ip_planes = nz - 1
    end if
    ip_glo = ip * (nz - 1) + ip_jlo
    ip_rank = rank_of_coord(ip)
    rep_counts(ip_rank + 1) = ld * ny * ip_planes
    rep_displs(ip_rank + 1) = ld * ny * (ip_glo - 1)
end do

if (coord == 0) then
    send_jlo = 1
else
    send_jlo = 2
end if
if ((coord == 0) .or. (coord == nproc - 1)) then
    send_planes = nz
else
    send_planes = nz - 1
end if
send_count = ld * ny * send_planes

cuda_istat = cudaDeviceSynchronize()
if (cuda_istat /= 0) stop 'replicated tridag pre-allgather sync failed'
call mpi_allgatherv(r(1,1,send_jlo), send_count, MPI_RPREC, rep_sol(1,1,1),   &
    rep_counts, rep_displs, MPI_RPREC, comm, ierr)
if (ierr /= 0) stop 'replicated tridag allgatherv failed'

if (.not. rep_coeff_ready) then
!$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, lh - 1
        if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
            bet_val = -1._rprec
            rep_gam(jx,jy,1) = 0._rprec
            rep_inv_bet(jx,jy,1) = 1._rprec / bet_val
            do g = 2, ng
                if (g == 2) then
                    c_prev = 1._rprec
                else
                    c_prev = c(jx,jy,2)
                end if
                if (g == ng) then
                    a_val = -1._rprec
                    b_val = 1._rprec
                else
                    a_val = a(jx,jy,2)
                    b_val = b(jx,jy,2)
                end if
                gam_val = c_prev / bet_val
                rep_gam(jx,jy,g) = gam_val
                bet_val = b_val - a_val * gam_val
                rep_inv_bet(jx,jy,g) = 1._rprec / bet_val
            end do
        end if
    end do
    end do
    cuda_istat = cudaDeviceSynchronize()
    if (cuda_istat /= 0) stop 'replicated tridag coeff sync failed'
    cuda_istat = cudaGetLastError()
    if (cuda_istat /= 0) stop 'replicated tridag coeff kernel failed'
    rep_coeff_ready = .true.
end if

!$cuf kernel do(2) <<<*,*>>>
do jy = 1, ny
do jx = 1, lh - 1
    if ((jy /= ny/2 + 1) .and. (jx*jy /= 1)) then
        ii = 2*jx
        ir = ii - 1
        rep_sol(ir,jy,1) = rep_sol(ir,jy,1) * rep_inv_bet(jx,jy,1)
        rep_sol(ii,jy,1) = rep_sol(ii,jy,1) * rep_inv_bet(jx,jy,1)
        do g = 2, ng
            if (g == ng) then
                a_val = -1._rprec
            else
                a_val = a(jx,jy,2)
            end if
            rep_sol(ir,jy,g) = (rep_sol(ir,jy,g) - a_val *                    &
                rep_sol(ir,jy,g-1)) * rep_inv_bet(jx,jy,g)
            rep_sol(ii,jy,g) = (rep_sol(ii,jy,g) - a_val *                    &
                rep_sol(ii,jy,g-1)) * rep_inv_bet(jx,jy,g)
        end do
        do g = ng - 1, 1, -1
            rep_sol(ir,jy,g) = rep_sol(ir,jy,g) - rep_gam(jx,jy,g+1) *         &
                rep_sol(ir,jy,g+1)
            rep_sol(ii,jy,g) = rep_sol(ii,jy,g) - rep_gam(jx,jy,g+1) *         &
                rep_sol(ii,jy,g+1)
        end do
    end if
end do
end do
cuda_istat = cudaDeviceSynchronize()
if (cuda_istat /= 0) stop 'replicated tridag solve sync failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'replicated tridag solve kernel failed'

!$cuf kernel do(3) <<<*,*>>>
do j = 1, n
do jy = 1, ny
do jx = 1, ld
    g = coord * (nz - 1) + j
    if ((jx >= ld - 1) .or. (jy == ny/2 + 1) .or.                             &
        ((jy == 1) .and. (jx <= 2))) then
        u(jx,jy,j) = 0._rprec
    else
        u(jx,jy,j) = rep_sol(jx,jy,g)
    end if
end do
end do
end do
cuda_istat = cudaDeviceSynchronize()
if (cuda_istat /= 0) stop 'replicated tridag copy sync failed'
cuda_istat = cudaGetLastError()
if (cuda_istat /= 0) stop 'replicated tridag copy kernel failed'

end subroutine tridag_array_replicated_cuda
#endif

#else
!*******************************************************************************
subroutine tridag_array(a, b, c, r, u)
!*******************************************************************************
use types, only : rprec
use param
implicit none

real(rprec),dimension(lh,ny,nz+1), intent(in) :: a, b, c

!  u and r are interleaved as complex arrays
real(rprec), dimension(ld,ny,nz+1), intent(in) :: r
real(rprec), dimension(ld,ny,nz+1), intent(out) :: u

integer :: n
integer :: jx, jy, j, j_min, j_max
real(rprec) :: bet(lh, ny)
real(rprec), dimension(lh,ny,nz+1) :: gam
integer :: ir, ii

n = nz+1

if (coord == 0) then
    do jy = 1, ny
        do jx = 1, lh-1

        if (b(jx, jy, 1) == 0._rprec) then
            write (*, *) 'tridag_array: rewrite eqs, jx, jy= ', jx, jy
            stop
        end if

        ii = 2*jx
        ir = ii - 1
        u(ir:ii,jy,1) = r(ir:ii,jy,1) / b(jx,jy,1)
        end do
    end do
    bet = b(:, :, 1)
    j_min = 1  ! this is only for backward pass
else
    j_min = 2  ! this is only for backward pass
end if

j_max = n

do j = 2, j_max
    do jy = 1, ny
        if (jy == ny/2+1) cycle
        do jx = 1, lh-1
            if (jx*jy == 1) cycle

            gam(jx, jy, j) = c(jx, jy, j-1) / bet(jx, jy)
            bet(jx, jy) = b(jx, jy, j) - a(jx, jy, j)*gam(jx, jy, j)

            if (bet(jx, jy) == 0._rprec) then
                write (*, *) 'tridag_array failed at jx,jy,j=', jx, jy, j
                write (*, *) 'a,b,c,gam,bet=', a(jx, jy, j), b(jx, jy, j),     &
                    c(jx, jy, j), gam(jx, jy, j), bet(jx, jy)
                stop
            end if

            !  u and r are interleaved
            ii = 2*jx
            ir = ii - 1
            u(ir:ii, jy, j) = (r(ir:ii, jy, j) - a(jx, jy, j) *                &
                u(ir:ii, jy, j-1)) /  bet(jx, jy)
        end do
    end do
end do

do j = n-1, j_min, -1
    do jy = 1, ny
        if (jy == ny/2+1) cycle
        do jx = 1, lh-1
            if (jx*jy == 1) cycle
            ii = 2*jx
            ir = ii - 1
            u(ir:ii, jy, j) = u(ir:ii, jy, j) - gam(jx, jy, j+1) *             &
                u(ir:ii, jy, j+1)
        end do
    end do
end do

end subroutine tridag_array

#endif
#endif
