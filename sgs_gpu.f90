!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!  GPU port: see docs/gpu_module_contracts.md for current ownership context.
!!

!*******************************************************************************
module sgs_gpu_m
!*******************************************************************************
! GPU implementation of the steady-state subgrid stress computation, i.e. the
! parts of sgs_stag()/calc_Sij() that run *every* time step.  Runtime dynamic
! SGS updates are also available on the GPU for the supported production SGS
! values: standard dynamic, scale-dependent dynamic, similarity Lagrangian,
! and scale-dependent Lagrangian.
!
! What sgs_stag_gpu does:
!   1. Update lagran_dt scalar (host bookkeeping).
!   2. Compute S11..S33 from device-resident dudx..dwdz (calc_Sij_gpu).
!   3. Initialize Cs_opt2 = 0.03 on device at jt=1 if inilag (matches CPU
!      sgs_stag's one-shot init for dynamic models).
!   4. Compute Nu_t = |S| * Cs_opt2 * delta^2 fused into a single kernel
!      (no intermediate S(:,:) scratch buffer needed).
!   5. Wall BCs for txx,txy,tyy,tzz at jz=1 (coord 0) and jz=nz-1 (top coord),
!      including the w-grid stress slots that the wall path writes.
!   6. Bulk tau loop for jz_min..jz_max.
!   7. Apply the immersed-surface stress treatment on the device when Level
!      Set GPU support is enabled.
!   8. MPI sync of txz, tyz, directly from device memory when GPU-aware MPI
!      is available or through boundary slabs on host otherwise.
!   9. BOGUS lines for safety mode.
!
! Navigation map:
!   - dynamic coefficient paths: std_dynamic_pples_gpu,
!     scaledep_dynamic_pples_gpu, lagrange_Ssim_gpu, lagrange_Sdep_gpu
!   - timestep SGS entry point: sgs_stag_gpu
!   - strain-rate assembly: calc_Sij_gpu
!   - divstress output: divstress_uv_gpu and divstress_w_gpu
!
! Device residency contract (declared in sgs_param):
!   S11..S33, Nu_t, Cs_opt2 are `!$acc declare create`. S11..S33 and Nu_t are
!   internal scratch - always overwritten before read - so they don't need an
!   initial host->device push. Cs_opt2 is pushed once in initialize.f90 right
!   after `call initial()` (matches the u/v/w pattern). If a Lagrangian step
!   runs on the CPU (host sgs_stag -> lagrange_Sdep), main.f90 issues
!   `!$acc update device(Cs_opt2)` afterward to refresh the device copy.
!
! Inputs from sim_param (device-resident via sim_param's `!$acc declare create`):
!   dudx..dwdz - only the dudz/dvdz at jz=1 (and nz) were touched on host by
!                wallstress; main.f90 issues update_device on those slabs.
!
! Outputs to sim_param (device-resident):
!   txx, txy, txz, tyy, tyz, tzz remain on the device for divstress and the
!   pressure pipeline. Only MPI boundary slabs are staged when required.
!*******************************************************************************
#ifdef PPSGS_GPU
use types, only : rprec
use param, only : ld, nx, ny, nz, lbz, dz, BOGUS,                              &
                  coord, nproc, lbc_mom, ubc_mom,                              &
                  jt, jt_total, dt, DYN_init, cs_count, sgs, sgs_model,        &
                  use_cfl_dt, inilag, Co, wall_damp_exp, vonk
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
use sim_param, only : u, v, w
use sim_param, only : txx, txy, txz, tyy, tyz, tzz
use sgs_param, only : S11, S12, S13, S22, S23, S33, Nu_t, Cs_opt2,             &
                      delta, nu, lagran_dt, S, ee_now,                         &
                      L11, L12, L13, L22, L23, L33,                             &
                      M11, M12, M13, M22, M23, M33,                             &
                      Q11, Q12, Q13, Q22, Q23, Q33,                             &
                      S11_bar, S12_bar, S13_bar, S22_bar, S23_bar, S33_bar,     &
                      S11_hat, S12_hat, S13_hat, S22_hat, S23_hat, S33_hat,     &
                      S_S11_bar, S_S12_bar, S_S13_bar,                         &
                      S_S22_bar, S_S23_bar, S_S33_bar,                          &
                      S_S11_hat, S_S12_hat, S_S13_hat,                         &
                      S_S22_hat, S_S23_hat, S_S33_hat,                          &
                      u_bar, v_bar, w_bar, u_hat, v_hat, w_hat, S_bar, S_hat,  &
                      SGS_MODEL_SMAGORINSKY, SGS_MODEL_STANDARD_DYNAMIC,        &
                      SGS_MODEL_SCALE_DEP_DYNAMIC,                              &
                      SGS_MODEL_LAGRANGE_SIMILARITY,                            &
                      SGS_MODEL_LAGRANGE_SCALE_DEP
use derivatives_gpu_m, only : ddx_gpu, ddy_gpu, ddxy_gpu, ddz_uv_gpu, ddz_w_gpu
use lagrange_Sdep_gpu_m, only : lagrange_Sdep_gpu, lagrange_Ssim_gpu
use test_filtermodule, only : test_filter_plane_gpu, test_test_filter_plane_gpu
use sgs_stag_util, only : rtnewt
#if defined(PPLVLSET) && defined(PPLVLSET_GPU)
use level_set, only : level_set_BC, level_set_Cs
#endif

#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
use param, only : MPI_RPREC, down, up, comm, status, ierr
#endif

implicit none
save
private
public :: sgs_stag_gpu, divstress_uv_gpu, divstress_w_gpu

! Persistent device-resident scratch for divstress_uv_gpu / divstress_w_gpu.
! These used to be per-call automatic arrays bracketed by `!$acc data create`,
! so every timestep allocated/freed ~4 GB of device scratch (+ first-touch
! page-in under managed memory). At 3072x384x400 that churn dominated the SGS
! cost. Hoisted to module scope (allocated once, reused every step), mirroring
! the press_gpu scratch pattern. divstress_w reuses dtxdx/dtydy/dtzdz.
real(rprec), allocatable, dimension(:,:,:) :: dtxdx, dtydy, dtzdz
real(rprec), allocatable, dimension(:,:,:) :: dtxdx2, dtydy2, dtzdz2
logical :: divstress_initialized = .false.

contains

!*******************************************************************************
subroutine divstress_gpu_init()
!*******************************************************************************
! Lazily allocate the persistent divstress scratch on host + device, once.
!*******************************************************************************
implicit none
if (divstress_initialized) return
allocate(dtxdx (ld, ny, lbz:nz));  dtxdx  = 0._rprec
allocate(dtydy (ld, ny, lbz:nz));  dtydy  = 0._rprec
allocate(dtzdz (ld, ny, lbz:nz));  dtzdz  = 0._rprec
allocate(dtxdx2(ld, ny, lbz:nz));  dtxdx2 = 0._rprec
allocate(dtydy2(ld, ny, lbz:nz));  dtydy2 = 0._rprec
allocate(dtzdz2(ld, ny, lbz:nz));  dtzdz2 = 0._rprec
!$acc enter data copyin(dtxdx, dtydy, dtzdz, dtxdx2, dtydy2, dtzdz2)
divstress_initialized = .true.
end subroutine divstress_gpu_init

!*******************************************************************************
subroutine std_dynamic_pples_gpu()
!*******************************************************************************
! PPLES/OpenACC standard dynamic SGS coefficient update (sgs_model=2).
! Mirrors std_dynamic.f90 layer-by-layer, but leaves all plane work arrays and
! Cs_opt2(:,:,jz) on the device.  The later Nu_t kernel consumes Cs_opt2.
!*******************************************************************************
implicit none

integer :: jx, jy, jz
real(rprec) :: const, lm_sum, mm_sum, cs_val

do jz = 1, nz
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            if ((coord == 0) .and. (jz == 1)) then
                u_bar(jx,jy) = u(jx,jy,1)
                v_bar(jx,jy) = v(jx,jy,1)
                w_bar(jx,jy) = 0.25_rprec*w(jx,jy,2)
            else
                u_bar(jx,jy) = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
                v_bar(jx,jy) = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
                w_bar(jx,jy) = w(jx,jy,jz)
            end if
            L11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
            L23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
            L33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
        end do
    end do

    call test_filter_plane_gpu(u_bar)
    call test_filter_plane_gpu(v_bar)
    call test_filter_plane_gpu(w_bar)

    call test_filter_plane_gpu(L11)
    call test_filter_plane_gpu(L12)
    call test_filter_plane_gpu(L13)
    call test_filter_plane_gpu(L22)
    call test_filter_plane_gpu(L23)
    call test_filter_plane_gpu(L33)

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            L11(jx,jy) = L11(jx,jy) - u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = L12(jx,jy) - u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = L13(jx,jy) - u_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = L22(jx,jy) - v_bar(jx,jy)*v_bar(jx,jy)
            L23(jx,jy) = L23(jx,jy) - v_bar(jx,jy)*w_bar(jx,jy)
            L33(jx,jy) = L33(jx,jy) - w_bar(jx,jy)*w_bar(jx,jy)

            S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +  &
                S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +               &
                S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
            S11_bar(jx,jy) = S11(jx,jy,jz)
            S12_bar(jx,jy) = S12(jx,jy,jz)
            S13_bar(jx,jy) = S13(jx,jy,jz)
            S22_bar(jx,jy) = S22(jx,jy,jz)
            S23_bar(jx,jy) = S23(jx,jy,jz)
            S33_bar(jx,jy) = S33(jx,jy,jz)
            S_S11_bar(jx,jy) = S(jx,jy)*S11(jx,jy,jz)
            S_S12_bar(jx,jy) = S(jx,jy)*S12(jx,jy,jz)
            S_S13_bar(jx,jy) = S(jx,jy)*S13(jx,jy,jz)
            S_S22_bar(jx,jy) = S(jx,jy)*S22(jx,jy,jz)
            S_S23_bar(jx,jy) = S(jx,jy)*S23(jx,jy,jz)
            S_S33_bar(jx,jy) = S(jx,jy)*S33(jx,jy,jz)
        end do
    end do

    call test_filter_plane_gpu(S11_bar)
    call test_filter_plane_gpu(S12_bar)
    call test_filter_plane_gpu(S13_bar)
    call test_filter_plane_gpu(S22_bar)
    call test_filter_plane_gpu(S23_bar)
    call test_filter_plane_gpu(S33_bar)

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            S_bar(jx,jy) = sqrt(2._rprec*(S11_bar(jx,jy)**2 +                 &
                S22_bar(jx,jy)**2 + S33_bar(jx,jy)**2 + 2._rprec*(            &
                S12_bar(jx,jy)**2 + S13_bar(jx,jy)**2 + S23_bar(jx,jy)**2)))
        end do
    end do

    call test_filter_plane_gpu(S_S11_bar)
    call test_filter_plane_gpu(S_S12_bar)
    call test_filter_plane_gpu(S_S13_bar)
    call test_filter_plane_gpu(S_S22_bar)
    call test_filter_plane_gpu(S_S23_bar)
    call test_filter_plane_gpu(S_S33_bar)

    const = 2._rprec*delta**2
    lm_sum = 0._rprec
    mm_sum = 0._rprec
    !$acc parallel loop collapse(2) default(present)                           &
    !$acc          reduction(+:lm_sum,mm_sum) async(1)
    do jy = 1, ny
        do jx = 1, ld
            M11(jx,jy) = const*(S_S11_bar(jx,jy) -                            &
                4._rprec*S_bar(jx,jy)*S11_bar(jx,jy))
            M12(jx,jy) = const*(S_S12_bar(jx,jy) -                            &
                4._rprec*S_bar(jx,jy)*S12_bar(jx,jy))
            M13(jx,jy) = const*(S_S13_bar(jx,jy) -                            &
                4._rprec*S_bar(jx,jy)*S13_bar(jx,jy))
            M22(jx,jy) = const*(S_S22_bar(jx,jy) -                            &
                4._rprec*S_bar(jx,jy)*S22_bar(jx,jy))
            M23(jx,jy) = const*(S_S23_bar(jx,jy) -                            &
                4._rprec*S_bar(jx,jy)*S23_bar(jx,jy))
            M33(jx,jy) = const*(S_S33_bar(jx,jy) -                            &
                4._rprec*S_bar(jx,jy)*S33_bar(jx,jy))
            lm_sum = lm_sum + L11(jx,jy)*M11(jx,jy) +                         &
                L22(jx,jy)*M22(jx,jy) + L33(jx,jy)*M33(jx,jy) +               &
                2._rprec*(L12(jx,jy)*M12(jx,jy) +                             &
                L13(jx,jy)*M13(jx,jy) + L23(jx,jy)*M23(jx,jy))
            mm_sum = mm_sum + M11(jx,jy)**2 + M22(jx,jy)**2 +                 &
                M33(jx,jy)**2 + 2._rprec*(M12(jx,jy)**2 +                     &
                M13(jx,jy)**2 + M23(jx,jy)**2)
        end do
    end do
    !$acc wait(1)

    if (mm_sum > 0._rprec) then
        cs_val = max(0._rprec, lm_sum/mm_sum)
    else
        cs_val = 0._rprec
    end if

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            Cs_opt2(jx,jy,jz) = cs_val
            ee_now(jx,jy,jz) = L11(jx,jy)**2 + L22(jx,jy)**2 +               &
                L33(jx,jy)**2 + 2._rprec*(L12(jx,jy)**2 +                    &
                L13(jx,jy)**2 + L23(jx,jy)**2) -                             &
                2._rprec*(L11(jx,jy)*M11(jx,jy) +                            &
                L22(jx,jy)*M22(jx,jy) + L33(jx,jy)*M33(jx,jy) +              &
                2._rprec*(L12(jx,jy)*M12(jx,jy) +                            &
                L13(jx,jy)*M13(jx,jy) + L23(jx,jy)*M23(jx,jy)))*cs_val +     &
                (M11(jx,jy)**2 + M22(jx,jy)**2 + M33(jx,jy)**2 +             &
                2._rprec*(M12(jx,jy)**2 + M13(jx,jy)**2 +                    &
                M23(jx,jy)**2))*cs_val**2
        end do
    end do
end do

end subroutine std_dynamic_pples_gpu

!*******************************************************************************
subroutine scaledep_dynamic_pples_gpu()
!*******************************************************************************
! PPLES/OpenACC scale-dependent dynamic SGS coefficient update (sgs_model=3).
! Mirrors scaledep_dynamic.f90, with per-plane reductions returning only the
! scalar polynomial coefficients to the host for rtnewt().
!*******************************************************************************
implicit none

integer :: jx, jy, jz
real(rprec) :: const, beta_val, cs_val, lm_sum, mm_sum
real(rprec) :: a1, b1, c1, d1, e1, a2, b2, c2, d2, e2
real(rprec), dimension(0:5) :: A

do jz = 1, nz
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            if ((coord == 0) .and. (jz == 1)) then
                u_bar(jx,jy) = u(jx,jy,1)
                v_bar(jx,jy) = v(jx,jy,1)
                w_bar(jx,jy) = 0.25_rprec*w(jx,jy,2)
            else
                u_bar(jx,jy) = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
                v_bar(jx,jy) = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
                w_bar(jx,jy) = w(jx,jy,jz)
            end if
            L11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
            L23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
            L33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
            u_hat(jx,jy) = u_bar(jx,jy)
            v_hat(jx,jy) = v_bar(jx,jy)
            w_hat(jx,jy) = w_bar(jx,jy)
        end do
    end do

    call test_filter_plane_gpu(u_bar)
    call test_filter_plane_gpu(v_bar)
    call test_filter_plane_gpu(w_bar)
    call test_filter_plane_gpu(L11)
    call test_filter_plane_gpu(L12)
    call test_filter_plane_gpu(L13)
    call test_filter_plane_gpu(L22)
    call test_filter_plane_gpu(L23)
    call test_filter_plane_gpu(L33)

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            L11(jx,jy) = L11(jx,jy) - u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = L12(jx,jy) - u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = L13(jx,jy) - u_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = L22(jx,jy) - v_bar(jx,jy)*v_bar(jx,jy)
            L23(jx,jy) = L23(jx,jy) - v_bar(jx,jy)*w_bar(jx,jy)
            L33(jx,jy) = L33(jx,jy) - w_bar(jx,jy)*w_bar(jx,jy)
            Q11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
            Q12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
            Q13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
            Q22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
            Q23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
            Q33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
        end do
    end do

    call test_test_filter_plane_gpu(u_hat)
    call test_test_filter_plane_gpu(v_hat)
    call test_test_filter_plane_gpu(w_hat)
    call test_test_filter_plane_gpu(Q11)
    call test_test_filter_plane_gpu(Q12)
    call test_test_filter_plane_gpu(Q13)
    call test_test_filter_plane_gpu(Q22)
    call test_test_filter_plane_gpu(Q23)
    call test_test_filter_plane_gpu(Q33)

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            Q11(jx,jy) = Q11(jx,jy) - u_hat(jx,jy)*u_hat(jx,jy)
            Q12(jx,jy) = Q12(jx,jy) - u_hat(jx,jy)*v_hat(jx,jy)
            Q13(jx,jy) = Q13(jx,jy) - u_hat(jx,jy)*w_hat(jx,jy)
            Q22(jx,jy) = Q22(jx,jy) - v_hat(jx,jy)*v_hat(jx,jy)
            Q23(jx,jy) = Q23(jx,jy) - v_hat(jx,jy)*w_hat(jx,jy)
            Q33(jx,jy) = Q33(jx,jy) - w_hat(jx,jy)*w_hat(jx,jy)

            S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +  &
                S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +               &
                S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
            S11_bar(jx,jy) = S11(jx,jy,jz)
            S12_bar(jx,jy) = S12(jx,jy,jz)
            S13_bar(jx,jy) = S13(jx,jy,jz)
            S22_bar(jx,jy) = S22(jx,jy,jz)
            S23_bar(jx,jy) = S23(jx,jy,jz)
            S33_bar(jx,jy) = S33(jx,jy,jz)
            S11_hat(jx,jy) = S11_bar(jx,jy)
            S12_hat(jx,jy) = S12_bar(jx,jy)
            S13_hat(jx,jy) = S13_bar(jx,jy)
            S22_hat(jx,jy) = S22_bar(jx,jy)
            S23_hat(jx,jy) = S23_bar(jx,jy)
            S33_hat(jx,jy) = S33_bar(jx,jy)
            S_S11_bar(jx,jy) = S(jx,jy)*S11(jx,jy,jz)
            S_S12_bar(jx,jy) = S(jx,jy)*S12(jx,jy,jz)
            S_S13_bar(jx,jy) = S(jx,jy)*S13(jx,jy,jz)
            S_S22_bar(jx,jy) = S(jx,jy)*S22(jx,jy,jz)
            S_S23_bar(jx,jy) = S(jx,jy)*S23(jx,jy,jz)
            S_S33_bar(jx,jy) = S(jx,jy)*S33(jx,jy,jz)
            S_S11_hat(jx,jy) = S_S11_bar(jx,jy)
            S_S12_hat(jx,jy) = S_S12_bar(jx,jy)
            S_S13_hat(jx,jy) = S_S13_bar(jx,jy)
            S_S22_hat(jx,jy) = S_S22_bar(jx,jy)
            S_S23_hat(jx,jy) = S_S23_bar(jx,jy)
            S_S33_hat(jx,jy) = S_S33_bar(jx,jy)
        end do
    end do

    call test_filter_plane_gpu(S11_bar)
    call test_filter_plane_gpu(S12_bar)
    call test_filter_plane_gpu(S13_bar)
    call test_filter_plane_gpu(S22_bar)
    call test_filter_plane_gpu(S23_bar)
    call test_filter_plane_gpu(S33_bar)
    call test_test_filter_plane_gpu(S11_hat)
    call test_test_filter_plane_gpu(S12_hat)
    call test_test_filter_plane_gpu(S13_hat)
    call test_test_filter_plane_gpu(S22_hat)
    call test_test_filter_plane_gpu(S23_hat)
    call test_test_filter_plane_gpu(S33_hat)

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            S_bar(jx,jy) = sqrt(2._rprec*(S11_bar(jx,jy)**2 +                 &
                S22_bar(jx,jy)**2 + S33_bar(jx,jy)**2 + 2._rprec*(            &
                S12_bar(jx,jy)**2 + S13_bar(jx,jy)**2 + S23_bar(jx,jy)**2)))
            S_hat(jx,jy) = sqrt(2._rprec*(S11_hat(jx,jy)**2 +                 &
                S22_hat(jx,jy)**2 + S33_hat(jx,jy)**2 + 2._rprec*(            &
                S12_hat(jx,jy)**2 + S13_hat(jx,jy)**2 + S23_hat(jx,jy)**2)))
        end do
    end do

    call test_filter_plane_gpu(S_S11_bar)
    call test_filter_plane_gpu(S_S12_bar)
    call test_filter_plane_gpu(S_S13_bar)
    call test_filter_plane_gpu(S_S22_bar)
    call test_filter_plane_gpu(S_S23_bar)
    call test_filter_plane_gpu(S_S33_bar)
    call test_test_filter_plane_gpu(S_S11_hat)
    call test_test_filter_plane_gpu(S_S12_hat)
    call test_test_filter_plane_gpu(S_S13_hat)
    call test_test_filter_plane_gpu(S_S22_hat)
    call test_test_filter_plane_gpu(S_S23_hat)
    call test_test_filter_plane_gpu(S_S33_hat)

    a1 = 0._rprec; b1 = 0._rprec; c1 = 0._rprec; d1 = 0._rprec; e1 = 0._rprec
    a2 = 0._rprec; b2 = 0._rprec; c2 = 0._rprec; d2 = 0._rprec; e2 = 0._rprec
    !$acc parallel loop collapse(2) default(present)                           &
    !$acc reduction(+:a1,b1,c1,d1,e1,a2,b2,c2,d2,e2) async(1)
    do jy = 1, ny
        do jx = 1, ld
            a1 = a1 + S_bar(jx,jy)*(S11_bar(jx,jy)*L11(jx,jy) +              &
                S22_bar(jx,jy)*L22(jx,jy) + S33_bar(jx,jy)*L33(jx,jy) +      &
                2._rprec*(S12_bar(jx,jy)*L12(jx,jy) +                       &
                S13_bar(jx,jy)*L13(jx,jy) + S23_bar(jx,jy)*L23(jx,jy)))
            b1 = b1 + S_S11_bar(jx,jy)*L11(jx,jy) +                         &
                S_S22_bar(jx,jy)*L22(jx,jy) + S_S33_bar(jx,jy)*L33(jx,jy) +  &
                2._rprec*(S_S12_bar(jx,jy)*L12(jx,jy) +                     &
                S_S13_bar(jx,jy)*L13(jx,jy) + S_S23_bar(jx,jy)*L23(jx,jy))
            c1 = c1 + S_S11_bar(jx,jy)**2 + S_S22_bar(jx,jy)**2 +            &
                S_S33_bar(jx,jy)**2 + 2._rprec*(S_S12_bar(jx,jy)**2 +        &
                S_S13_bar(jx,jy)**2 + S_S23_bar(jx,jy)**2)
            d1 = d1 + 0.5_rprec*S_bar(jx,jy)**4
            e1 = e1 + S_bar(jx,jy)*(S11_bar(jx,jy)*S_S11_bar(jx,jy) +        &
                S22_bar(jx,jy)*S_S22_bar(jx,jy) +                            &
                S33_bar(jx,jy)*S_S33_bar(jx,jy) + 2._rprec*(                &
                S12_bar(jx,jy)*S_S12_bar(jx,jy) +                            &
                S13_bar(jx,jy)*S_S13_bar(jx,jy) +                            &
                S23_bar(jx,jy)*S_S23_bar(jx,jy)))
            a2 = a2 + S_hat(jx,jy)*(S11_hat(jx,jy)*Q11(jx,jy) +              &
                S22_hat(jx,jy)*Q22(jx,jy) + S33_hat(jx,jy)*Q33(jx,jy) +      &
                2._rprec*(S12_hat(jx,jy)*Q12(jx,jy) +                       &
                S13_hat(jx,jy)*Q13(jx,jy) + S23_hat(jx,jy)*Q23(jx,jy)))
            b2 = b2 + S_S11_hat(jx,jy)*Q11(jx,jy) +                         &
                S_S22_hat(jx,jy)*Q22(jx,jy) + S_S33_hat(jx,jy)*Q33(jx,jy) +  &
                2._rprec*(S_S12_hat(jx,jy)*Q12(jx,jy) +                     &
                S_S13_hat(jx,jy)*Q13(jx,jy) + S_S23_hat(jx,jy)*Q23(jx,jy))
            c2 = c2 + S_S11_hat(jx,jy)**2 + S_S22_hat(jx,jy)**2 +            &
                S_S33_hat(jx,jy)**2 + 2._rprec*(S_S12_hat(jx,jy)**2 +        &
                S_S13_hat(jx,jy)**2 + S_S23_hat(jx,jy)**2)
            d2 = d2 + 0.5_rprec*S_hat(jx,jy)**4
            e2 = e2 + S_hat(jx,jy)*(S11_hat(jx,jy)*S_S11_hat(jx,jy) +        &
                S22_hat(jx,jy)*S_S22_hat(jx,jy) +                            &
                S33_hat(jx,jy)*S_S33_hat(jx,jy) + 2._rprec*(                &
                S12_hat(jx,jy)*S_S12_hat(jx,jy) +                            &
                S13_hat(jx,jy)*S_S13_hat(jx,jy) +                            &
                S23_hat(jx,jy)*S_S23_hat(jx,jy)))
        end do
    end do
    !$acc wait(1)

    a1 = -2._rprec*(delta**2)*4._rprec*a1/(nx*ny)
    b1 = -2._rprec*(delta**2)*b1/(nx*ny)
    c1 = (2._rprec*delta**2)**2*c1/(nx*ny)
    d1 = (2._rprec*delta**2)**2*16._rprec*d1/(nx*ny)
    e1 = 2._rprec*(2._rprec*delta**2)**2*4._rprec*e1/(nx*ny)
    a2 = -2._rprec*(delta**2)*16._rprec*a2/(nx*ny)
    b2 = -2._rprec*(delta**2)*b2/(nx*ny)
    c2 = (2._rprec*delta**2)**2*c2/(nx*ny)
    d2 = (2._rprec*delta**2)**2*256._rprec*d2/(nx*ny)
    e2 = 2._rprec*(2._rprec*delta**2)**2*16._rprec*e2/(nx*ny)

    A(0) = b2*c1 - b1*c2
    A(1) = a1*c2 - b2*e1
    A(2) = b2*d1 + b1*e2 - a2*c1
    A(3) = a2*e1 - a1*e2
    A(4) = -a2*d1 - b1*d2
    A(5) = a1*d2
    beta_val = rtnewt(A, jz)

    const = 2._rprec*delta**2
    lm_sum = 0._rprec
    mm_sum = 0._rprec
    !$acc parallel loop collapse(2) default(present)                           &
    !$acc reduction(+:lm_sum,mm_sum) async(1)
    do jy = 1, ny
        do jx = 1, ld
            M11(jx,jy) = const*(S_S11_bar(jx,jy) -                            &
                4._rprec*beta_val*S_bar(jx,jy)*S11_bar(jx,jy))
            M12(jx,jy) = const*(S_S12_bar(jx,jy) -                            &
                4._rprec*beta_val*S_bar(jx,jy)*S12_bar(jx,jy))
            M13(jx,jy) = const*(S_S13_bar(jx,jy) -                            &
                4._rprec*beta_val*S_bar(jx,jy)*S13_bar(jx,jy))
            M22(jx,jy) = const*(S_S22_bar(jx,jy) -                            &
                4._rprec*beta_val*S_bar(jx,jy)*S22_bar(jx,jy))
            M23(jx,jy) = const*(S_S23_bar(jx,jy) -                            &
                4._rprec*beta_val*S_bar(jx,jy)*S23_bar(jx,jy))
            M33(jx,jy) = const*(S_S33_bar(jx,jy) -                            &
                4._rprec*beta_val*S_bar(jx,jy)*S33_bar(jx,jy))
            lm_sum = lm_sum + L11(jx,jy)*M11(jx,jy) +                         &
                L22(jx,jy)*M22(jx,jy) + L33(jx,jy)*M33(jx,jy) +               &
                2._rprec*(L12(jx,jy)*M12(jx,jy) +                             &
                L13(jx,jy)*M13(jx,jy) + L23(jx,jy)*M23(jx,jy))
            mm_sum = mm_sum + M11(jx,jy)**2 + M22(jx,jy)**2 +                 &
                M33(jx,jy)**2 + 2._rprec*(M12(jx,jy)**2 +                     &
                M13(jx,jy)**2 + M23(jx,jy)**2)
        end do
    end do
    !$acc wait(1)

    if (mm_sum > 0._rprec) then
        cs_val = max(0._rprec, lm_sum/mm_sum)
    else
        cs_val = 0._rprec
    end if

    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, ld
            Cs_opt2(jx,jy,jz) = cs_val
            ee_now(jx,jy,jz) = L11(jx,jy)**2 + L22(jx,jy)**2 +               &
                L33(jx,jy)**2 + 2._rprec*(L12(jx,jy)**2 +                    &
                L13(jx,jy)**2 + L23(jx,jy)**2) -                             &
                2._rprec*(L11(jx,jy)*M11(jx,jy) +                            &
                L22(jx,jy)*M22(jx,jy) + L33(jx,jy)*M33(jx,jy) +              &
                2._rprec*(L12(jx,jy)*M12(jx,jy) +                            &
                L13(jx,jy)*M13(jx,jy) + L23(jx,jy)*M23(jx,jy)))*cs_val +     &
                (M11(jx,jy)**2 + M22(jx,jy)**2 + M33(jx,jy)**2 +             &
                2._rprec*(M12(jx,jy)**2 + M13(jx,jy)**2 +                    &
                M23(jx,jy)**2))*cs_val**2
        end do
    end do
end do

end subroutine scaledep_dynamic_pples_gpu

!*******************************************************************************
subroutine sgs_stag_gpu()
!*******************************************************************************
! Always-on SGS path (sgs_model 2-5, jt < DYN_init OR mod(jt_total,cs_count)/=0).
! Matches the CPU sgs_stag() flow, but with all stencil work on the device.
!*******************************************************************************
implicit none

integer  :: jx, jy, jz
integer  :: jz_min, jz_max
real(rprec) :: const, const2, const3, const4
real(rprec) :: dsq, wall_exp, zloc, lsgs, lsq, domain_height

#if defined(PPLVLSET) && defined(PPLVLSET_GPU)
! Level Set wall-stress and derivative inputs can be produced by a mixture of
! default-stream and queue-1 kernels. Close that ownership boundary before the
! device SGS construction consumes them.
!$acc wait
#endif

! ----- 1. Host scalar bookkeeping (Lagrangian dt accumulator). -----
if (use_cfl_dt) then
    if (sgs_model == SGS_MODEL_LAGRANGE_SIMILARITY .or.                       &
        sgs_model == SGS_MODEL_LAGRANGE_SCALE_DEP) then
        if (jt_total >= DYN_init - cs_count + 1) then
            lagran_dt = lagran_dt + dt
        end if
    end if
else
    lagran_dt = cs_count * dt
end if

! ----- 2. Strain-rate tensor on device. -----
call calc_Sij_gpu()

! ----- 3. One-shot Cs_opt2 init at jt=1 (mirrors CPU sgs_stag line 187-189). -----
if (sgs .and. sgs_model /= SGS_MODEL_SMAGORINSKY .and. jt == 1 .and.          &
    inilag) then
    if (coord == 0) write(*,*) 'CS_opt2 initialiazed (GPU)'
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz
        do jy = 1, ny
            do jx = 1, ld
                Cs_opt2(jx, jy, jz) = 0.03_rprec
            end do
        end do
    end do
end if

! ----- 3b. Dynamic-coefficient updates.
! Mirrors the CPU sgs_stag dispatch EXACTLY (sgs_stag_util.f90:192):
!   jt_total >= DYN_init .AND. mod(jt_total, cs_count) == 0
! The cumulative counter preserves the dynamic-model schedule across a restart;
! the local jt counter resets to one for every resumed segment.
if (sgs .and. jt_total >= DYN_init .and.                                       &
    mod(jt_total, cs_count) == 0) then
    select case (sgs_model)
    case (SGS_MODEL_STANDARD_DYNAMIC)
        if (jt_total == DYN_init .and. coord == 0) write(*,*)                   &
            'running dynamic sgs_model = ', sgs_model, ' (GPU)'
        call std_dynamic_pples_gpu()
    case (SGS_MODEL_SCALE_DEP_DYNAMIC)
        if (jt_total == DYN_init .and. coord == 0) write(*,*)                   &
            'running dynamic sgs_model = ', sgs_model, ' (GPU)'
        call scaledep_dynamic_pples_gpu()
    case (SGS_MODEL_LAGRANGE_SIMILARITY)
        call lagrange_Ssim_gpu()
    case (SGS_MODEL_LAGRANGE_SCALE_DEP)
        call lagrange_Sdep_gpu()
    end select
end if

! ----- 4. tau region: txx..tzz are device-resident (sim_param declare create).
!   wallstress' jz=1 / jz=nz writes to txz/tyz were synced device-side by
!   main.f90 (update_device after wallstress block).
!$acc data present(txx, txy, txz, tyy, tyz, tzz,                               &
!$acc              S11, S12, S13, S22, S23, S33, Nu_t, Cs_opt2)

! ----- 5. Nu_t = sqrt(2*Sij^2) * coefficient * filter_length^2.
!   For sgs_model=5, the filter size is constant delta and Cs_opt2 is updated
!   by lagrange_Sdep_gpu.  For sgs_model=1, mirror the CPU Smagorinsky branch:
!   Cs_opt2 = Co^2 and l(z) follows Mason wall damping.  When SGS is disabled,
!   skip Nu_t entirely; the later stress kernels use the molecular branches and
!   do not read Nu_t.
dsq = delta * delta
#if defined(PPLVLSET) && defined(PPLVLSET_GPU)
if (sgs .and. sgs_model == SGS_MODEL_SMAGORINSKY) call level_set_Cs(delta)
#endif
if (sgs) then
    wall_exp = real(wall_damp_exp, rprec)
    domain_height = real((nz - 1) * nproc, rprec) * dz
    !$acc parallel loop collapse(3) default(present)                            &
    !$acc          private(const, zloc, lsgs, lsq) async(1)
    do jz = 1, nz
        do jy = 1, ny
            do jx = 1, nx
                const = sqrt(2._rprec * (S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +&
                    S33(jx,jy,jz)**2 + 2._rprec * (S12(jx,jy,jz)**2 +          &
                    S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))

                if (sgs_model == SGS_MODEL_SMAGORINSKY) then
#if defined(PPLVLSET) && defined(PPLVLSET_GPU)
                    Nu_t(jx,jy,jz) = const * Cs_opt2(jx,jy,jz) * dsq
#else
                    if (lbc_mom == 0 .and. ubc_mom == 0) then
                        lsq = dsq
                    else if (lbc_mom > 0 .and. ubc_mom == 0) then
                        if (coord == 0 .and. jz == 1) then
                            zloc = 0.5_rprec * dz
                        else
                            zloc = real((jz - 1) + coord * (nz - 1), rprec) * dz
                        end if
                        lsgs = (Co**wall_exp * (vonk*zloc)**(-wall_exp) +      &
                            delta**(-wall_exp))**(-1._rprec/wall_exp)
                        lsq = lsgs * lsgs
                    else if (lbc_mom > 0 .and. ubc_mom > 0) then
                        if (coord == 0 .and. jz == 1) then
                            zloc = 0.5_rprec * dz
                        else if (coord == nproc - 1 .and. jz == nz) then
                            zloc = 0.5_rprec * dz
                        else
                            zloc = real((jz - 1) + coord * (nz - 1), rprec) * dz
                            zloc = min(zloc, domain_height - zloc)
                        end if
                        lsgs = (Co**wall_exp * (vonk*zloc)**(-wall_exp) +      &
                            delta**(-wall_exp))**(-1._rprec/wall_exp)
                        lsq = lsgs * lsgs
                    else
                        if (coord == nproc - 1 .and. jz == nz) then
                            zloc = 0.5_rprec * dz
                        else
                            zloc = real((nproc - coord) * (nz - 1) -           &
                                (jz - 1), rprec) * dz
                        end if
                        lsgs = (Co**wall_exp * (vonk*zloc)**(-wall_exp) +      &
                            delta**(-wall_exp))**(-1._rprec/wall_exp)
                        lsq = lsgs * lsgs
                    end if
                    Nu_t(jx,jy,jz) = const * Co * Co * lsq
#endif
                else
                    Nu_t(jx,jy,jz) = const * Cs_opt2(jx,jy,jz) * dsq
                end if
            end do
        end do
    end do
end if

! ----- 6. Wall BC for jz=1 (coord==0 only). -----
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free: txx, txy, tyy, tzz on w-nodes (averaged across jz=1,2)
        case (0)
            if (sgs) then
                !$acc parallel loop collapse(2) default(present) private(const) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        const = 0.5_rprec*(Nu_t(jx,jy,1) + Nu_t(jx,jy,2)) + nu
                        txx(jx,jy,1) = -const*(S11(jx,jy,1) + S11(jx,jy,2))
                        txy(jx,jy,1) = -const*(S12(jx,jy,1) + S12(jx,jy,2))
                        tyy(jx,jy,1) = -const*(S22(jx,jy,1) + S22(jx,jy,2))
                        tzz(jx,jy,1) = -const*(S33(jx,jy,1) + S33(jx,jy,2))
                    end do
                end do
            else
                !$acc parallel loop collapse(2) default(present) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,1) = -nu*(S11(jx,jy,1) + S11(jx,jy,2))
                        txy(jx,jy,1) = -nu*(S12(jx,jy,1) + S12(jx,jy,2))
                        tyy(jx,jy,1) = -nu*(S22(jx,jy,1) + S22(jx,jy,2))
                        tzz(jx,jy,1) = -nu*(S33(jx,jy,1) + S33(jx,jy,2))
                    end do
                end do
            end if

        ! Wall: Sij stored on uvp-nodes at jz=1.
        case (1:)
            if (sgs) then
                !$acc parallel loop collapse(2) default(present) private(const) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        const = -2._rprec*(Nu_t(jx,jy,1) + nu)
                        txx(jx,jy,1) = const*S11(jx,jy,1)
                        txy(jx,jy,1) = const*S12(jx,jy,1)
                        tyy(jx,jy,1) = const*S22(jx,jy,1)
                        tzz(jx,jy,1) = const*S33(jx,jy,1)
                    end do
                end do
            else
                !$acc parallel loop collapse(2) default(present) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,1) = -2._rprec*nu*S11(jx,jy,1)
                        txy(jx,jy,1) = -2._rprec*nu*S12(jx,jy,1)
                        tyy(jx,jy,1) = -2._rprec*nu*S22(jx,jy,1)
                        tzz(jx,jy,1) = -2._rprec*nu*S33(jx,jy,1)
                    end do
                end do
            end if
    end select
    jz_min = 2
else
    jz_min = 1
end if

! ----- 7. Wall BC for jz=nz (coord==nproc-1 only). -----
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        case (0)
            if (sgs) then
                !$acc parallel loop collapse(2) default(present)            &
                !$acc          private(const, const2) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        const  = 0.5_rprec*(Nu_t(jx,jy,nz-1) + Nu_t(jx,jy,nz)) + nu
                        const2 = 2._rprec*(Nu_t(jx,jy,nz-1) + nu)
                        txx(jx,jy,nz-1) = -const*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                        txy(jx,jy,nz-1) = -const*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                        tyy(jx,jy,nz-1) = -const*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                        tzz(jx,jy,nz-1) = -const*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                        txz(jx,jy,nz-1) = -const2*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = -const2*S23(jx,jy,nz-1)
                    end do
                end do
            else
                !$acc parallel loop collapse(2) default(present) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,nz-1) = -nu*(S11(jx,jy,nz-1) + S11(jx,jy,nz))
                        txy(jx,jy,nz-1) = -nu*(S12(jx,jy,nz-1) + S12(jx,jy,nz))
                        tyy(jx,jy,nz-1) = -nu*(S22(jx,jy,nz-1) + S22(jx,jy,nz))
                        tzz(jx,jy,nz-1) = -nu*(S33(jx,jy,nz-1) + S33(jx,jy,nz))
                        txz(jx,jy,nz-1) = -2._rprec*nu*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = -2._rprec*nu*S23(jx,jy,nz-1)
                    end do
                end do
            end if

        ! Wall
        case (1:)
            if (sgs) then
                !$acc parallel loop collapse(2) default(present)            &
                !$acc          private(const, const2) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        const  = -2._rprec*(Nu_t(jx,jy,nz) + nu)
                        const2 = -2._rprec*(Nu_t(jx,jy,nz-1) + nu)
                        txx(jx,jy,nz-1) = const*S11(jx,jy,nz)
                        txy(jx,jy,nz-1) = const*S12(jx,jy,nz)
                        tyy(jx,jy,nz-1) = const*S22(jx,jy,nz)
                        tzz(jx,jy,nz-1) = const*S33(jx,jy,nz)
                        txz(jx,jy,nz-1) = const2*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = const2*S23(jx,jy,nz-1)
                    end do
                end do
            else
                !$acc parallel loop collapse(2) default(present) async(1)
                do jy = 1, ny
                    do jx = 1, nx
                        txx(jx,jy,nz-1) = -2._rprec*nu*S11(jx,jy,nz-1)
                        txy(jx,jy,nz-1) = -2._rprec*nu*S12(jx,jy,nz-1)
                        tyy(jx,jy,nz-1) = -2._rprec*nu*S22(jx,jy,nz-1)
                        tzz(jx,jy,nz-1) = -2._rprec*nu*S33(jx,jy,nz-1)
                        txz(jx,jy,nz-1) = -2._rprec*nu*S13(jx,jy,nz-1)
                        tyz(jx,jy,nz-1) = -2._rprec*nu*S23(jx,jy,nz-1)
                    end do
                end do
            end if
    end select
    jz_max = nz - 2
else
    jz_max = nz - 1
end if

! ----- 8. Bulk tau (jz_min..jz_max). -----
if (sgs) then
    const3 = -2._rprec*nu*0.5_rprec
    const4 = -2._rprec*nu
    !$acc parallel loop collapse(3) default(present) private(const, const2) async(1)
    do jz = jz_min, jz_max
        do jy = 1, ny
            do jx = 1, nx
                const  = -0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))
                const2 = -2._rprec*Nu_t(jx,jy,jz)
                txx(jx,jy,jz) = (const + const3)*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
                txy(jx,jy,jz) = (const + const3)*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
                tyy(jx,jy,jz) = (const + const3)*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
                tzz(jx,jy,jz) = (const + const3)*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
                txz(jx,jy,jz) = (const2 + const4)*S13(jx,jy,jz)
                tyz(jx,jy,jz) = (const2 + const4)*S23(jx,jy,jz)
            end do
        end do
    end do
else
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = jz_min, jz_max
        do jy = 1, ny
            do jx = 1, nx
                txx(jx,jy,jz) = -nu*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
                txy(jx,jy,jz) = -nu*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
                tyy(jx,jy,jz) = -nu*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
                tzz(jx,jy,jz) = -nu*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
                txz(jx,jy,jz) = -2._rprec*nu*S13(jx,jy,jz)
                tyz(jx,jy,jz) = -2._rprec*nu*S23(jx,jy,jz)
            end do
        end do
    end do
end if

!$acc wait(1)
!$acc end data    ! txx..tzz remain device-resident

#if defined(PPLVLSET) && defined(PPLVLSET_GPU)
! Apply immersed-surface stress treatment after the dense stress construction.
call level_set_BC()
#endif

! ----- 9. MPI sync of txz, tyz across ranks (DOWN: send slab 1, recv slab nz).
!   txz, tyz are device-resident.
#ifdef PPMPI
#ifdef PPGPU_AWARE_MPI
! GPU-aware MPI path: send GPU pointers straight to MPICH. mpi_sync_real_array
! is host-only, so use explicit sendrecv calls with host_data use_device
! wrappers.
!$acc wait(1)
!$acc host_data use_device(txz, tyz)
call mpi_sendrecv (txz(:,:,1), ld*ny, MPI_RPREC, down, 21,                     &
                   txz(:,:,nz), ld*ny, MPI_RPREC, up, 21,                      &
                   comm, status, ierr)
call mpi_sendrecv (tyz(:,:,1), ld*ny, MPI_RPREC, down, 22,                     &
                   tyz(:,:,nz), ld*ny, MPI_RPREC, up, 22,                      &
                   comm, status, ierr)
!$acc end host_data
#else
! Sandwich fall-back: pull outbound slab, MPI on host, push inbound.
!$acc wait(1)
!$acc update self(txz(:,:,1), tyz(:,:,1))
call mpi_sync_real_array(txz, 0, MPI_SYNC_DOWN)
call mpi_sync_real_array(tyz, 0, MPI_SYNC_DOWN)
!$acc wait(1)
!$acc update device(txz(:,:,nz), tyz(:,:,nz))
#endif

#ifdef PPSAFETYMODE
! Write BOGUS into the lbz=0 ghost slabs on the device (matches CPU semantics).
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        txx(jx, jy, 0) = BOGUS
        txy(jx, jy, 0) = BOGUS
        txz(jx, jy, 0) = BOGUS
        tyy(jx, jy, 0) = BOGUS
        tyz(jx, jy, 0) = BOGUS
        tzz(jx, jy, 0) = BOGUS
    end do
end do
#endif
#endif

#ifdef PPSAFETYMODE
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        txx(jx, jy, nz) = BOGUS
        txy(jx, jy, nz) = BOGUS
        tyy(jx, jy, nz) = BOGUS
        tzz(jx, jy, nz) = BOGUS
    end do
end do
#endif

end subroutine sgs_stag_gpu


!*******************************************************************************
subroutine calc_Sij_gpu()
!*******************************************************************************
! GPU equivalent of calc_Sij. Builds the strain-rate tensor S11..S33 on w-nodes
! (jz=1..nz) from the device-resident velocity derivatives.
!
! Boundary slabs (jz=1 on coord 0, jz=nz on top coord) get special-case
! formulas depending on lbc_mom / ubc_mom. The bulk loop computes the rest by
! averaging slab jz with slab jz-1. A single MPI_SYNC_DOWN of dwdz fills the
! coord+1 ghost slab needed at jz=nz on non-top ranks.
!*******************************************************************************
implicit none

integer :: jx, jy, jz
integer :: jz_min, jz_max
real(rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz

! ----- 1. jz=1 BC (coord==0 only). -----
if (coord == 0) then
    select case (lbc_mom)

        ! Stress free: Sij on w-nodes.
        case (0)
            !$acc parallel loop collapse(2) default(present)                  &
            !$acc          private(ux, uy, uz, vx, vy, vz, wx, wy, wz) async(1)
            do jy = 1, ny
                do jx = 1, nx
                    ux = dudx(jx,jy,1)
                    uy = dudy(jx,jy,1)
                    uz = dudz(jx,jy,1)
                    vx = dvdx(jx,jy,1)
                    vy = dvdy(jx,jy,1)
                    vz = dvdz(jx,jy,1)
                    wx = dwdx(jx,jy,1)
                    wy = dwdy(jx,jy,1)
                    wz = 0.5_rprec*(dwdz(jx,jy,1) + 0._rprec)

                    S11(jx,jy,1) = ux
                    S12(jx,jy,1) = 0.5_rprec*(uy + vx)
                    S13(jx,jy,1) = 0.5_rprec*(uz + wx)
                    S22(jx,jy,1) = vy
                    S23(jx,jy,1) = 0.5_rprec*(vz + wy)
                    S33(jx,jy,1) = wz
                end do
            end do

        ! Wall: Sij on uvp-nodes at jz=1 (dudz, dvdz from wallstress).
        case (1:)
            !$acc parallel loop collapse(2) default(present) private(wx, wy) async(1)
            do jy = 1, ny
                do jx = 1, nx
                    S11(jx,jy,1) = dudx(jx,jy,1)
                    S12(jx,jy,1) = 0.5_rprec*(dudy(jx,jy,1) + dvdx(jx,jy,1))
                    wx = 0.5_rprec*(dwdx(jx,jy,1) + dwdx(jx,jy,2))
                    S13(jx,jy,1) = 0.5_rprec*(dudz(jx,jy,1) + wx)
                    S22(jx,jy,1) = dvdy(jx,jy,1)
                    wy = 0.5_rprec*(dwdy(jx,jy,1) + dwdy(jx,jy,2))
                    S23(jx,jy,1) = 0.5_rprec*(dvdz(jx,jy,1) + wy)
                    S33(jx,jy,1) = dwdz(jx,jy,1)
                end do
            end do
    end select
    jz_min = 2
else
    jz_min = 1
end if

! ----- 2. jz=nz BC (coord==nproc-1 only). -----
if (coord == nproc-1) then
    select case (ubc_mom)

        ! Stress free
        case (0)
            !$acc parallel loop collapse(2) default(present)                  &
            !$acc          private(ux, uy, uz, vx, vy, vz, wx, wy, wz) async(1)
            do jy = 1, ny
                do jx = 1, nx
                    ux = dudx(jx,jy,nz-1)
                    uy = dudy(jx,jy,nz-1)
                    uz = dudz(jx,jy,nz)   ! from wallstress (zero)
                    vx = dvdx(jx,jy,nz-1)
                    vy = dvdy(jx,jy,nz-1)
                    vz = dvdz(jx,jy,nz)   ! from wallstress (zero)
                    wx = dwdx(jx,jy,nz)
                    wy = dwdy(jx,jy,nz)
                    wz = 0.5_rprec*(dwdz(jx,jy,nz-1) + 0._rprec)

                    S11(jx,jy,nz) = ux
                    S12(jx,jy,nz) = 0.5_rprec*(uy + vx)
                    S13(jx,jy,nz) = 0.5_rprec*(uz + wx)
                    S22(jx,jy,nz) = vy
                    S23(jx,jy,nz) = 0.5_rprec*(vz + wy)
                    S33(jx,jy,nz) = wz
                end do
            end do

        ! Wall
        case (1:)
            !$acc parallel loop collapse(2) default(present) private(wx, wy) async(1)
            do jy = 1, ny
                do jx = 1, nx
                    S11(jx,jy,nz) = dudx(jx,jy,nz-1)
                    S12(jx,jy,nz) = 0.5_rprec*(dudy(jx,jy,nz-1) + dvdx(jx,jy,nz-1))
                    wx = 0.5_rprec*(dwdx(jx,jy,nz-1) + dwdx(jx,jy,nz))
                    S13(jx,jy,nz) = 0.5_rprec*(dudz(jx,jy,nz) + wx)
                    S22(jx,jy,nz) = dvdy(jx,jy,nz-1)
                    wy = 0.5_rprec*(dwdy(jx,jy,nz-1) + dwdy(jx,jy,nz))
                    S23(jx,jy,nz) = 0.5_rprec*(dvdz(jx,jy,nz) + wy)
                    S33(jx,jy,nz) = dwdz(jx,jy,nz-1)
                end do
            end do
    end select
    jz_max = nz - 1
else
    jz_max = nz
end if

! ----- 3. MPI sync of dwdz so the bulk loop can read dwdz(:,:,nz) on
!         non-top ranks. dwdz is device-resident.
#ifdef PPMPI
#ifdef PPGPU_AWARE_MPI
! GPU-aware MPI path: send the device-resident dwdz(:,:,1) slab directly.
! This is the production path because the SGS gradients remain device-resident
! after derivative assembly.
!$acc wait(1)
!$acc host_data use_device(dwdz)
call mpi_sendrecv (dwdz(:,:,1), ld*ny, MPI_RPREC, down, 23,                    &
                   dwdz(:,:,nz), ld*ny, MPI_RPREC, up, 23,                     &
                   comm, status, ierr)
!$acc end host_data
#else
! Host-staged fallback for non-GPU-aware MPI builds.  Refresh only the slab
! that mpi_sync_real_array sends, then push the received top ghost slab back to
! the device.  Keep this path correct for validation, but benchmark production
! runs with PPGPU_AWARE_MPI so the halo exchange stays device-resident.
!$acc wait(1)
!$acc update self(dwdz(:,:,1))
call mpi_sync_real_array(dwdz(:,:,1:), 1, MPI_SYNC_DOWN)
! On non-top ranks, dwdz(:,:,nz) was overwritten with coord+1's slab 1.
! Push it back to the device so the bulk Sij kernel sees it.
if (coord /= nproc - 1) then
    !$acc wait(1)
    !$acc update device(dwdz(:,:,nz))
end if
#endif
#endif

! ----- 4. Bulk Sij (jz_min..jz_max). -----
!$acc parallel loop collapse(3) default(present) private(uy, vx) async(1)
do jz = jz_min, jz_max
    do jy = 1, ny
        do jx = 1, nx
            S11(jx,jy,jz) = 0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))
            uy = (dudy(jx,jy,jz) + dudy(jx,jy,jz-1))
            vx = (dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))
            S12(jx,jy,jz) = 0.25_rprec*(uy + vx)
            S13(jx,jy,jz) = 0.5_rprec*(dudz(jx,jy,jz) + dwdx(jx,jy,jz))
            S22(jx,jy,jz) = 0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))
            S23(jx,jy,jz) = 0.5_rprec*(dvdz(jx,jy,jz) + dwdy(jx,jy,jz))
            S33(jx,jy,jz) = 0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))
        end do
    end do
end do

end subroutine calc_Sij_gpu


!*******************************************************************************
subroutine divstress_uv_gpu(divtx, divty)
!*******************************************************************************
! GPU equivalent of divstress_uv. Computes
!   divtx = dx(txx) + dy(txy) + dz(txz)
!   divty = dx(txy) + dy(tyy) + dz(tyz)
! using batched cuFFT for the spectral x/y derivatives and the existing
! ddz_w_gpu finite-difference kernel for z derivatives.
!
! Inputs (txx, txy, txz, tyy, tyz) are device-resident via sim_param's
! `!$acc declare create`. divtx/divty are also device-resident (sim_param
! declare create). The RHS update in main.f90 runs as a GPU kernel that reads
! divtx/divty directly on device - no copyout needed.
!*******************************************************************************
implicit none
real(rprec), dimension(:,:,lbz:), intent(out) :: divtx, divty

integer :: jx, jy, jz

! dtxdx..dtzdz2 are module-level persistent scratch (allocated once in
! divstress_gpu_init). The `!$acc data create` below is present-or-create, so
! after the enter-data it finds them already on the device - no per-step malloc.
call divstress_gpu_init()

!$acc data create(dtxdx, dtydy, dtzdz, dtxdx2, dtydy2, dtzdz2)                 &
!$acc      present(txx, txy, txz, tyy, tyz, divtx, divty)

! ----- Spectral derivatives (each call is fully on the device). -----
call ddx_gpu (txx, dtxdx, lbz)            ! dx(txx)
call ddy_gpu (tyy, dtydy2, lbz)           ! dy(tyy)
call ddxy_gpu(txy, dtxdx2, dtydy, lbz)    ! dx(txy) -> dtxdx2; dy(txy) -> dtydy

! ----- Finite-difference z derivatives (uv-grid output from w-grid input). -----
call ddz_w_gpu(txz, dtzdz,  lbz)          ! dz(txz)
call ddz_w_gpu(tyz, dtzdz2, lbz)          ! dz(tyz)

#ifdef PPSAFETYMODE
#ifdef PPMPI
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        dtxdx (jx, jy, 0) = BOGUS
        dtydy2(jx, jy, 0) = BOGUS
        dtxdx2(jx, jy, 0) = BOGUS
        dtydy (jx, jy, 0) = BOGUS
        dtzdz (jx, jy, 0) = BOGUS
        dtzdz2(jx, jy, 0) = BOGUS
    end do
end do
#endif
#endif

! ----- Combine: divtx, divty for jz = 1..nz-1.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        do jx = 1, nx
            divtx(jx, jy, jz) = dtxdx (jx, jy, jz) + dtydy (jx, jy, jz)        &
                                                    + dtzdz (jx, jy, jz)
            divty(jx, jy, jz) = dtxdx2(jx, jy, jz) + dtydy2(jx, jy, jz)        &
                                                    + dtzdz2(jx, jy, jz)
        end do
    end do
end do

! ----- Zero ld-1, ld (Nyquist padding). -----
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        do jx = ld-1, ld
            divtx(jx, jy, jz) = 0._rprec
            divty(jx, jy, jz) = 0._rprec
        end do
    end do
end do

! ----- BOGUS the unused ghost slabs on device, so the copyout to host
!       carries the safety-mode markers (matches CPU behavior).
#ifdef PPSAFETYMODE
#ifdef PPMPI
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        divtx(jx, jy, 0) = BOGUS
        divty(jx, jy, 0) = BOGUS
    end do
end do
#endif
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        divtx(jx, jy, nz) = BOGUS
        divty(jx, jy, nz) = BOGUS
    end do
end do
#endif

!$acc wait(1)
!$acc end data    ! divtx, divty stay device-resident (declare create)

end subroutine divstress_uv_gpu


!*******************************************************************************
subroutine divstress_w_gpu(divtz)
!*******************************************************************************
! GPU equivalent of divstress_w. Computes
!   divtz = dx(txz) + dy(tyz) + dz(tzz)   (with wall-BC adjustments at jz=1/nz)
!
! tx (=txz), ty (=tyz), tz (=tzz) are device-resident via sim_param.
! divtz is also device-resident (sim_param declare create); the RHS update in
! main.f90 reads divtz directly on device.
!*******************************************************************************
implicit none
real(rprec), dimension(:,:,lbz:), intent(out) :: divtz

integer :: jx, jy, jz

! dtxdx/dtydy/dtzdz are module-level persistent scratch shared with
! divstress_uv_gpu (allocated once in divstress_gpu_init). present-or-create.
call divstress_gpu_init()

!$acc data create(dtxdx, dtydy, dtzdz)                                         &
!$acc      present(txz, tyz, tzz, divtz)

! ----- dx(txz), dy(tyz), dz(tzz). -----
call ddx_gpu  (txz, dtxdx, lbz)
call ddy_gpu  (tyz, dtydy, lbz)
call ddz_uv_gpu(tzz, dtzdz, lbz)

#ifdef PPSAFETYMODE
#ifdef PPMPI
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        dtxdx(jx, jy, 0) = BOGUS
        dtydy(jx, jy, 0) = BOGUS
        dtzdz(jx, jy, 0) = BOGUS
    end do
end do
#endif
#endif

#ifdef PPSAFETYMODE
#ifdef PPMPI
!$acc parallel loop collapse(2) default(present) async(1)
do jy = 1, ny
    do jx = 1, ld
        divtz(jx, jy, 0) = BOGUS
    end do
end do
#endif
#endif

! ----- jz = 1: at bottom wall, drop the dz term (CPU does the same, since
!       d(tzz)/dz at the wall is assumed zero).
if (coord == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            divtz(jx, jy, 1) = dtxdx(jx, jy, 1) + dtydy(jx, jy, 1)
        end do
    end do
else
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            divtz(jx, jy, 1) = dtxdx(jx, jy, 1) + dtydy(jx, jy, 1)             &
                                                + dtzdz(jx, jy, 1)
        end do
    end do
end if

! ----- jz = nz: top-wall analogue; only top rank drops dz. -----
if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            divtz(jx, jy, nz) = dtxdx(jx, jy, nz) + dtydy(jx, jy, nz)
        end do
    end do
else
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
        do jx = 1, nx
            divtz(jx, jy, nz) = dtxdx(jx, jy, nz) + dtydy(jx, jy, nz)          &
                                                  + dtzdz(jx, jy, nz)
        end do
    end do
end if

! ----- Bulk: jz = 2..nz-1.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 2, nz-1
    do jy = 1, ny
        do jx = 1, nx
            divtz(jx, jy, jz) = dtxdx(jx, jy, jz) + dtydy(jx, jy, jz)          &
                                                  + dtzdz(jx, jy, jz)
        end do
    end do
end do

! ----- Zero ld-1, ld (Nyquist padding). -----
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz-1
    do jy = 1, ny
        do jx = ld-1, ld
            divtz(jx, jy, jz) = 0._rprec
        end do
    end do
end do

!$acc wait(1)
!$acc end data    ! divtz stays device-resident (declare create)

end subroutine divstress_w_gpu

#endif
end module sgs_gpu_m
