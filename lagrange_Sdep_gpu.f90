!*******************************************************************************
! GPU port of the Lagrangian Scale-Dependent (LASD, sgs_model=5) dynamic
! coefficient update. Mirrors the CPU `lagrange_Sdep` in lagrange_Sdep.f90,
! but uses a *batched-over-jz* layout: every per-jz scratch buffer (u_bar,
! L11, M11, ...) is allocated as a (ld, ny, nz) 3-D device array, and every
! call to test_filter_gpu is replaced by a single batched cuFFT over all jz
! slabs. The original (ld, ny) sgs_param scratch arrays remain in place for
! the CPU sgs_stag fall-back paths (sgs_model 1/4) and the level-set
! features that are not yet ported.
!
! Why batched: the naive per-jz translation (initial G2 implementation) ran
! ~10 000 GPU operations per Lagrangian step (42 filter calls x 65 jz x 3
! cuFFT/multiply ops + ~30 small kernels per jz). Each cuFFT invocation on a
! batch=1 plan carries ~5-10 ms of host-side launch + sync overhead, which
! dominated runtime (~105 s/step). Batching across jz collapses the per-step
! cost to a handful of large cuFFT launches (~150 GPU ops total) - same
! computational work, but cuFFT runs at full throughput.
!
! interpolag_Sdep is still on the host (G3 will port it). Velocity / F_*
! sync wrappers around it are the only PCIe round-trip in this routine.
!
! Compiled only when USE_LES_GPU=ON (CPP macro PPSGS_GPU defined).
!*******************************************************************************
module lagrange_Sdep_gpu_m
#ifdef PPSGS_GPU
use types, only : rprec
implicit none
save
private

public :: lagrange_Sdep_gpu, lagrange_Ssim_gpu, lagrange_Sdep_gpu_init,        &
          interpolag_Sdep_gpu, interpolag_Ssim_gpu

! Routine map:
!   - lagrange_Sdep_gpu_init: allocate and register batched device scratch
!   - lagrange_Ssim_gpu: similarity SGS coefficient path
!   - lagrange_Sdep_gpu: scale-dependent SGS coefficient path
!   - interpolag_*_gpu: GPU wrappers around interpolation/restart state
!   - sync_downup_F: vertical neighbor exchange for Lagrangian histories
!
! Ownership map:
!   - sgs_stag_util.f90 owns runtime SGS dispatch and calls this module.
!   - lagrange_Sdep_gpu_init owns allocation of 3-D GPU scratch and tiny grid
!     mirrors used by interpolation kernels.
!   - sgs_param owns the persistent F_* Lagrangian history arrays; this module
!     updates their GPU path and keeps backup tempF_* buffers for old values.
!   - CPU lagrange_Sdep.f90 remains the semantic reference for model algebra.

! ---------------------------------------------------------------------------
! 3-D scratch arrays, all (ld, ny, nz). Allocated once via
! lagrange_Sdep_gpu_init and made device-resident by `!$acc declare create`.
! Naming: same as CPU lagrange_Sdep variable but with `_3d` suffix to make
! the dimension change visible at the call site.
! ---------------------------------------------------------------------------
real(rprec), dimension(:,:,:), allocatable :: u_bar_3d, v_bar_3d, w_bar_3d
real(rprec), dimension(:,:,:), allocatable :: u_hat_3d, v_hat_3d, w_hat_3d
real(rprec), dimension(:,:,:), allocatable :: L11_3d, L12_3d, L13_3d,           &
                                               L22_3d, L23_3d, L33_3d
real(rprec), dimension(:,:,:), allocatable :: Q11_3d, Q12_3d, Q13_3d,           &
                                               Q22_3d, Q23_3d, Q33_3d
real(rprec), dimension(:,:,:), allocatable :: S11_bar_3d, S12_bar_3d,           &
                                               S13_bar_3d, S22_bar_3d,          &
                                               S23_bar_3d, S33_bar_3d
real(rprec), dimension(:,:,:), allocatable :: S11_hat_3d, S12_hat_3d,           &
                                               S13_hat_3d, S22_hat_3d,          &
                                               S23_hat_3d, S33_hat_3d
real(rprec), dimension(:,:,:), allocatable :: S_S11_bar_3d, S_S12_bar_3d,       &
                                               S_S13_bar_3d, S_S22_bar_3d,      &
                                               S_S23_bar_3d, S_S33_bar_3d
real(rprec), dimension(:,:,:), allocatable :: S_S11_hat_3d, S_S12_hat_3d,       &
                                               S_S13_hat_3d, S_S22_hat_3d,      &
                                               S_S23_hat_3d, S_S33_hat_3d
real(rprec), dimension(:,:,:), allocatable :: S_3d, S_bar_3d, S_hat_3d
! Per-jz contractions and Lagrangian timescale, persist across the F_LM/F_QN
! kernels so they can be reused by Cs_opt2 step.
real(rprec), dimension(:,:,:), allocatable :: LM_3d, MM_3d, QN_3d, NN_3d
real(rprec), dimension(:,:,:), allocatable :: Tn_3d
real(rprec), dimension(:,:,:), allocatable :: Cs_2d_3d, Cs_4d_3d

!$acc declare create(u_bar_3d, v_bar_3d, w_bar_3d,                              &
!$acc                u_hat_3d, v_hat_3d, w_hat_3d,                              &
!$acc                L11_3d, L12_3d, L13_3d, L22_3d, L23_3d, L33_3d,            &
!$acc                Q11_3d, Q12_3d, Q13_3d, Q22_3d, Q23_3d, Q33_3d,            &
!$acc                S11_bar_3d, S12_bar_3d, S13_bar_3d,                        &
!$acc                S22_bar_3d, S23_bar_3d, S33_bar_3d,                        &
!$acc                S11_hat_3d, S12_hat_3d, S13_hat_3d,                        &
!$acc                S22_hat_3d, S23_hat_3d, S33_hat_3d,                        &
!$acc                S_S11_bar_3d, S_S12_bar_3d, S_S13_bar_3d,                  &
!$acc                S_S22_bar_3d, S_S23_bar_3d, S_S33_bar_3d,                  &
!$acc                S_S11_hat_3d, S_S12_hat_3d, S_S13_hat_3d,                  &
!$acc                S_S22_hat_3d, S_S23_hat_3d, S_S33_hat_3d,                  &
!$acc                S_3d, S_bar_3d, S_hat_3d,                                  &
!$acc                LM_3d, MM_3d, QN_3d, NN_3d,                                &
!$acc                Tn_3d, Cs_2d_3d, Cs_4d_3d)

! Initialization-flag pair, mirroring the `save` logicals in the CPU
! lagrange_Sdep so the F_LM/MM/QN/NN one-shot init at jt==DYN_init happens
! exactly once per binary launch.
logical, save :: F_LM_MM_init = .false.
logical, save :: F_QN_NN_init = .false.

! ---------------------------------------------------------------------------
! Device-resident grid arrays used by interpolag_Sdep_gpu. Tiny: x_dev/y_dev
! are length nx+1/ny+1 (a few KB); z_dev/zw_dev are length nz (sub-KB);
! autowrap_i/j are length nx+2/ny+2 ints. Populated once by
! lagrange_Sdep_gpu_init from the grid_m %x/%y/%z/%zw/%autowrap_* pointers
! (which themselves are populated by grid%build()).
! ---------------------------------------------------------------------------
real(rprec), dimension(:), allocatable :: x_dev, y_dev, z_dev, zw_dev
integer,     dimension(:), allocatable :: autowrap_i_dev, autowrap_j_dev

!$acc declare create(x_dev, y_dev, z_dev, zw_dev,                               &
!$acc                autowrap_i_dev, autowrap_j_dev)

! ---------------------------------------------------------------------------
! Backup buffers for the F_* arrays read by trilinear interpolation. The
! CPU interpolag_Sdep uses `tempF_LM = F_LM` etc; we replicate that on the
! device so kernels write F_LM (new) by reading tempF_LM (old).
! ---------------------------------------------------------------------------
real(rprec), dimension(:,:,:), allocatable :: tempF_LM, tempF_MM,               &
                                               tempF_QN, tempF_NN
#ifdef PPDYN_TN
real(rprec), dimension(:,:,:), allocatable :: tempF_ee2, tempF_deedt2,          &
                                               tempee_past
#endif

!$acc declare create(tempF_LM, tempF_MM, tempF_QN, tempF_NN)
#ifdef PPDYN_TN
!$acc declare create(tempF_ee2, tempF_deedt2, tempee_past)
#endif

contains

!*******************************************************************************
subroutine lagrange_Sdep_gpu_init()
!*******************************************************************************
use param, only : ld, ny, nz, nx, lbz
use grid_m, only : grid
implicit none

if (allocated(u_bar_3d)) return

! Make sure the host-side grid arrays are populated before we copy them.
! cell_indx() and others call this lazily on first use; we trigger it here
! so the device copies below are valid.
if (.not. grid%built) call grid%build()

! Device-resident grid arrays.
allocate(x_dev(nx+1), y_dev(ny+1))
allocate(z_dev(lbz:nz), zw_dev(lbz:nz))
allocate(autowrap_i_dev(0:nx+1), autowrap_j_dev(0:ny+1))
x_dev          = grid%x
y_dev          = grid%y
z_dev          = grid%z
zw_dev         = grid%zw
autowrap_i_dev = grid%autowrap_i
autowrap_j_dev = grid%autowrap_j
!$acc wait(1)
!$acc update device(x_dev, y_dev, z_dev, zw_dev,                                &
!$acc               autowrap_i_dev, autowrap_j_dev)

! Backup F_* buffers for interpolag.
allocate(tempF_LM(ld, ny, lbz:nz), tempF_MM(ld, ny, lbz:nz),                    &
         tempF_QN(ld, ny, lbz:nz), tempF_NN(ld, ny, lbz:nz))
tempF_LM = 0._rprec; tempF_MM = 0._rprec
tempF_QN = 0._rprec; tempF_NN = 0._rprec
#ifdef PPDYN_TN
allocate(tempF_ee2(ld, ny, lbz:nz), tempF_deedt2(ld, ny, lbz:nz),               &
         tempee_past(ld, ny, lbz:nz))
tempF_ee2 = 0._rprec; tempF_deedt2 = 0._rprec; tempee_past = 0._rprec
#endif

allocate(u_bar_3d(ld,ny,nz), v_bar_3d(ld,ny,nz), w_bar_3d(ld,ny,nz))
allocate(u_hat_3d(ld,ny,nz), v_hat_3d(ld,ny,nz), w_hat_3d(ld,ny,nz))
allocate(L11_3d(ld,ny,nz), L12_3d(ld,ny,nz), L13_3d(ld,ny,nz),                  &
         L22_3d(ld,ny,nz), L23_3d(ld,ny,nz), L33_3d(ld,ny,nz))
allocate(Q11_3d(ld,ny,nz), Q12_3d(ld,ny,nz), Q13_3d(ld,ny,nz),                  &
         Q22_3d(ld,ny,nz), Q23_3d(ld,ny,nz), Q33_3d(ld,ny,nz))
allocate(S11_bar_3d(ld,ny,nz), S12_bar_3d(ld,ny,nz), S13_bar_3d(ld,ny,nz),      &
         S22_bar_3d(ld,ny,nz), S23_bar_3d(ld,ny,nz), S33_bar_3d(ld,ny,nz))
allocate(S11_hat_3d(ld,ny,nz), S12_hat_3d(ld,ny,nz), S13_hat_3d(ld,ny,nz),      &
         S22_hat_3d(ld,ny,nz), S23_hat_3d(ld,ny,nz), S33_hat_3d(ld,ny,nz))
allocate(S_S11_bar_3d(ld,ny,nz), S_S12_bar_3d(ld,ny,nz),                        &
         S_S13_bar_3d(ld,ny,nz), S_S22_bar_3d(ld,ny,nz),                        &
         S_S23_bar_3d(ld,ny,nz), S_S33_bar_3d(ld,ny,nz))
allocate(S_S11_hat_3d(ld,ny,nz), S_S12_hat_3d(ld,ny,nz),                        &
         S_S13_hat_3d(ld,ny,nz), S_S22_hat_3d(ld,ny,nz),                        &
         S_S23_hat_3d(ld,ny,nz), S_S33_hat_3d(ld,ny,nz))
allocate(S_3d(ld,ny,nz), S_bar_3d(ld,ny,nz), S_hat_3d(ld,ny,nz))
allocate(LM_3d(ld,ny,nz), MM_3d(ld,ny,nz), QN_3d(ld,ny,nz), NN_3d(ld,ny,nz))
allocate(Tn_3d(ld,ny,nz), Cs_2d_3d(ld,ny,nz), Cs_4d_3d(ld,ny,nz))

! Zero so device residency entries are well-defined before first use.
u_bar_3d = 0._rprec; v_bar_3d = 0._rprec; w_bar_3d = 0._rprec
u_hat_3d = 0._rprec; v_hat_3d = 0._rprec; w_hat_3d = 0._rprec
L11_3d = 0._rprec; L12_3d = 0._rprec; L13_3d = 0._rprec
L22_3d = 0._rprec; L23_3d = 0._rprec; L33_3d = 0._rprec
Q11_3d = 0._rprec; Q12_3d = 0._rprec; Q13_3d = 0._rprec
Q22_3d = 0._rprec; Q23_3d = 0._rprec; Q33_3d = 0._rprec
S11_bar_3d = 0._rprec; S12_bar_3d = 0._rprec; S13_bar_3d = 0._rprec
S22_bar_3d = 0._rprec; S23_bar_3d = 0._rprec; S33_bar_3d = 0._rprec
S11_hat_3d = 0._rprec; S12_hat_3d = 0._rprec; S13_hat_3d = 0._rprec
S22_hat_3d = 0._rprec; S23_hat_3d = 0._rprec; S33_hat_3d = 0._rprec
S_S11_bar_3d = 0._rprec; S_S12_bar_3d = 0._rprec; S_S13_bar_3d = 0._rprec
S_S22_bar_3d = 0._rprec; S_S23_bar_3d = 0._rprec; S_S33_bar_3d = 0._rprec
S_S11_hat_3d = 0._rprec; S_S12_hat_3d = 0._rprec; S_S13_hat_3d = 0._rprec
S_S22_hat_3d = 0._rprec; S_S23_hat_3d = 0._rprec; S_S33_hat_3d = 0._rprec
S_3d = 0._rprec; S_bar_3d = 0._rprec; S_hat_3d = 0._rprec
LM_3d = 0._rprec; MM_3d = 0._rprec; QN_3d = 0._rprec; NN_3d = 0._rprec
Tn_3d = 0._rprec; Cs_2d_3d = 0._rprec; Cs_4d_3d = 0._rprec

end subroutine lagrange_Sdep_gpu_init

!*******************************************************************************
subroutine lagrange_Ssim_gpu()
!*******************************************************************************
! Batched OpenACC Lagrangian scale-similarity SGS update (sgs_model=4).
! Mirrors lagrange_Ssim.f90, but filters all z planes as one batched operation.
!*******************************************************************************
use param
use sim_param, only : u, v, w
use sgs_param, only : F_LM, F_MM, Beta, Cs_opt2, opftime, lagran_dt,            &
                      S11, S12, S13, S22, S23, S33, delta, ee_now, Tn_all
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
use test_filtermodule, only : test_filter_b_gpu
implicit none

integer :: jx, jy, jz, istart, iend, iend_mod
real(rprec) :: const, opftdelta
real(rprec) :: ub, vb, wb
real(rprec) :: m11_loc, m12_loc, m13_loc, m22_loc, m23_loc, m33_loc
real(rprec) :: lm_loc, mm_loc, dumfac_loc, epsi_loc, tprod_loc
real(rprec), parameter :: eps = 1.e-32_rprec

opftdelta = opftime*delta
const = 2._rprec*delta*delta

! Model 4 uses beta=1 in the Mij term.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    Beta(jx,jy,jz) = 1._rprec
end do
end do
end do

call interpolag_Ssim_gpu()

! Interpolate u, v, w onto w nodes and build Lij first terms.
!$acc parallel loop collapse(3) default(present) private(ub, vb, wb) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    if ((coord == 0) .and. (jz == 1)) then
        ub = u(jx,jy,1)
        vb = v(jx,jy,1)
        wb = 0.25_rprec*w(jx,jy,2)
    else
        ub = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
        vb = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
        wb = w(jx,jy,jz)
    end if
    u_bar_3d(jx,jy,jz) = ub
    v_bar_3d(jx,jy,jz) = vb
    w_bar_3d(jx,jy,jz) = wb
    L11_3d(jx,jy,jz) = ub*ub
    L12_3d(jx,jy,jz) = ub*vb
    L13_3d(jx,jy,jz) = ub*wb
    L22_3d(jx,jy,jz) = vb*vb
    L23_3d(jx,jy,jz) = vb*wb
    L33_3d(jx,jy,jz) = wb*wb
end do
end do
end do

call test_filter_b_gpu(u_bar_3d, nz)
call test_filter_b_gpu(v_bar_3d, nz)
call test_filter_b_gpu(w_bar_3d, nz)
call test_filter_b_gpu(L11_3d,  nz)
call test_filter_b_gpu(L12_3d,  nz)
call test_filter_b_gpu(L13_3d,  nz)
call test_filter_b_gpu(L22_3d,  nz)
call test_filter_b_gpu(L23_3d,  nz)
call test_filter_b_gpu(L33_3d,  nz)

! Lij = filtered(ui uj) - filtered(ui) filtered(uj).
!$acc parallel loop collapse(3) default(present) private(ub, vb, wb) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    ub = u_bar_3d(jx,jy,jz)
    vb = v_bar_3d(jx,jy,jz)
    wb = w_bar_3d(jx,jy,jz)
    L11_3d(jx,jy,jz) = L11_3d(jx,jy,jz) - ub*ub
    L12_3d(jx,jy,jz) = L12_3d(jx,jy,jz) - ub*vb
    L13_3d(jx,jy,jz) = L13_3d(jx,jy,jz) - ub*wb
    L22_3d(jx,jy,jz) = L22_3d(jx,jy,jz) - vb*vb
    L23_3d(jx,jy,jz) = L23_3d(jx,jy,jz) - vb*wb
    L33_3d(jx,jy,jz) = L33_3d(jx,jy,jz) - wb*wb
end do
end do
end do

! Prepare Sij and |S|Sij terms.
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    S_3d(jx,jy,jz) = sqrt(2._rprec*(S11(jx,jy,jz)**2 +                       &
        S22(jx,jy,jz)**2 + S33(jx,jy,jz)**2 + 2._rprec*(                     &
        S12(jx,jy,jz)**2 + S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
    S11_bar_3d(jx,jy,jz) = S11(jx,jy,jz)
    S12_bar_3d(jx,jy,jz) = S12(jx,jy,jz)
    S13_bar_3d(jx,jy,jz) = S13(jx,jy,jz)
    S22_bar_3d(jx,jy,jz) = S22(jx,jy,jz)
    S23_bar_3d(jx,jy,jz) = S23(jx,jy,jz)
    S33_bar_3d(jx,jy,jz) = S33(jx,jy,jz)
    S_S11_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S11(jx,jy,jz)
    S_S12_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S12(jx,jy,jz)
    S_S13_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S13(jx,jy,jz)
    S_S22_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S22(jx,jy,jz)
    S_S23_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S23(jx,jy,jz)
    S_S33_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S33(jx,jy,jz)
end do
end do
end do

call test_filter_b_gpu(S11_bar_3d, nz)
call test_filter_b_gpu(S12_bar_3d, nz)
call test_filter_b_gpu(S13_bar_3d, nz)
call test_filter_b_gpu(S22_bar_3d, nz)
call test_filter_b_gpu(S23_bar_3d, nz)
call test_filter_b_gpu(S33_bar_3d, nz)

!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    S_bar_3d(jx,jy,jz) = sqrt(2._rprec*(S11_bar_3d(jx,jy,jz)**2 +             &
        S22_bar_3d(jx,jy,jz)**2 + S33_bar_3d(jx,jy,jz)**2 + 2._rprec*(       &
        S12_bar_3d(jx,jy,jz)**2 + S13_bar_3d(jx,jy,jz)**2 +                  &
        S23_bar_3d(jx,jy,jz)**2)))
end do
end do
end do

call test_filter_b_gpu(S_S11_bar_3d, nz)
call test_filter_b_gpu(S_S12_bar_3d, nz)
call test_filter_b_gpu(S_S13_bar_3d, nz)
call test_filter_b_gpu(S_S22_bar_3d, nz)
call test_filter_b_gpu(S_S23_bar_3d, nz)
call test_filter_b_gpu(S_S33_bar_3d, nz)

! Form Mij, contractions, and ee_now.
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(m11_loc, m12_loc, m13_loc, m22_loc, m23_loc, m33_loc,    &
!$acc                  lm_loc, mm_loc) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    m11_loc = const*(S_S11_bar_3d(jx,jy,jz) -                                 &
        4._rprec*Beta(jx,jy,jz)*S_bar_3d(jx,jy,jz)*S11_bar_3d(jx,jy,jz))
    m12_loc = const*(S_S12_bar_3d(jx,jy,jz) -                                 &
        4._rprec*Beta(jx,jy,jz)*S_bar_3d(jx,jy,jz)*S12_bar_3d(jx,jy,jz))
    m13_loc = const*(S_S13_bar_3d(jx,jy,jz) -                                 &
        4._rprec*Beta(jx,jy,jz)*S_bar_3d(jx,jy,jz)*S13_bar_3d(jx,jy,jz))
    m22_loc = const*(S_S22_bar_3d(jx,jy,jz) -                                 &
        4._rprec*Beta(jx,jy,jz)*S_bar_3d(jx,jy,jz)*S22_bar_3d(jx,jy,jz))
    m23_loc = const*(S_S23_bar_3d(jx,jy,jz) -                                 &
        4._rprec*Beta(jx,jy,jz)*S_bar_3d(jx,jy,jz)*S23_bar_3d(jx,jy,jz))
    m33_loc = const*(S_S33_bar_3d(jx,jy,jz) -                                 &
        4._rprec*Beta(jx,jy,jz)*S_bar_3d(jx,jy,jz)*S33_bar_3d(jx,jy,jz))
    lm_loc = L11_3d(jx,jy,jz)*m11_loc + L22_3d(jx,jy,jz)*m22_loc +            &
        L33_3d(jx,jy,jz)*m33_loc + 2._rprec*(L12_3d(jx,jy,jz)*m12_loc +       &
        L13_3d(jx,jy,jz)*m13_loc + L23_3d(jx,jy,jz)*m23_loc)
    mm_loc = m11_loc**2 + m22_loc**2 + m33_loc**2 +                           &
        2._rprec*(m12_loc**2 + m13_loc**2 + m23_loc**2)
    LM_3d(jx,jy,jz) = lm_loc
    MM_3d(jx,jy,jz) = mm_loc
    ee_now(jx,jy,jz) = L11_3d(jx,jy,jz)**2 + L22_3d(jx,jy,jz)**2 +            &
        L33_3d(jx,jy,jz)**2 + 2._rprec*(L12_3d(jx,jy,jz)**2 +                 &
        L13_3d(jx,jy,jz)**2 + L23_3d(jx,jy,jz)**2) -                          &
        2._rprec*lm_loc*Cs_opt2(jx,jy,jz) + mm_loc*Cs_opt2(jx,jy,jz)**2
end do
end do
end do

if (inilag .and. (.not. F_LM_MM_init) .and.                                   &
    (jt == cs_count .or. jt == DYN_init)) then
    if (coord == 0) print *, 'F_MM and F_LM initialized (GPU)'
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, ld
        F_MM(jx,jy,jz) = MM_3d(jx,jy,jz)
        F_LM(jx,jy,jz) = 0.025_rprec*MM_3d(jx,jy,jz)
    end do
    end do
    end do
    !$acc parallel loop collapse(2) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
        F_MM(ld-1,jy,jz) = 1._rprec
        F_MM(ld,  jy,jz) = 1._rprec
        F_LM(ld-1,jy,jz) = 1._rprec
        F_LM(ld,  jy,jz) = 1._rprec
    end do
    end do
    F_LM_MM_init = .true.
end if

if (inflow_type > 0) then
    iend = floor(fringe_region_end * nx + 1._rprec)
    istart = floor((fringe_region_end - fringe_region_len) * nx + 1._rprec)
    iend_mod = modulo(iend - 1, nx) + 1
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, ld
        MM_3d(jx,jy,jz) = max(MM_3d(jx,jy,jz),                                &
            0.1_rprec*const*S_3d(jx,jy,jz)**2)
        if ((iend <= nx .and. jx >= istart .and. jx <= iend) .or.             &
            (iend > nx .and. (jx >= istart .or. jx <= iend_mod))) then
            LM_3d(jx,jy,jz) = 0._rprec
            F_LM(jx,jy,jz) = 0._rprec
        end if
    end do
    end do
    end do
end if

! Running averages and Cs_opt2.
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(tprod_loc, dumfac_loc, epsi_loc) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
#ifdef PPDYN_TN
    Tn_3d(jx,jy,jz) = 4._rprec*pi*sqrt(F_ee2(jx,jy,jz)/F_deedt2(jx,jy,jz))
#else
    tprod_loc = max(F_LM(jx,jy,jz)*F_MM(jx,jy,jz), eps)
    Tn_3d(jx,jy,jz) = opftdelta/sqrt(sqrt(sqrt(tprod_loc)))
#endif
    dumfac_loc = lagran_dt/Tn_3d(jx,jy,jz)
    epsi_loc = dumfac_loc/(1._rprec + dumfac_loc)
    F_LM(jx,jy,jz) = epsi_loc*LM_3d(jx,jy,jz) +                               &
        (1._rprec - epsi_loc)*F_LM(jx,jy,jz)
    F_MM(jx,jy,jz) = epsi_loc*MM_3d(jx,jy,jz) +                               &
        (1._rprec - epsi_loc)*F_MM(jx,jy,jz)
    F_LM(jx,jy,jz) = max(eps, F_LM(jx,jy,jz))
#ifdef PPDYN_TN
    F_ee2(jx,jy,jz) = epsi_loc*ee_now(jx,jy,jz)**2 +                          &
        (1._rprec - epsi_loc)*F_ee2(jx,jy,jz)
    F_deedt2(jx,jy,jz) = epsi_loc*                                            &
        (((ee_now(jx,jy,jz) - ee_past(jx,jy,jz))/lagran_dt)**2) +             &
        (1._rprec - epsi_loc)*F_deedt2(jx,jy,jz)
    ee_past(jx,jy,jz) = ee_now(jx,jy,jz)
#endif
    if (jx >= ld-1) then
        Cs_opt2(jx,jy,jz) = eps
    else
        Cs_opt2(jx,jy,jz) = max(eps, F_LM(jx,jy,jz)/(F_MM(jx,jy,jz) + eps))
    end if
    Tn_all(jx,jy,jz) = Tn_3d(jx,jy,jz)
end do
end do
end do

#ifdef PPMPI
call sync_downup_F(F_LM)
call sync_downup_F(F_MM)
call sync_downup_F(Tn_all)
#ifdef PPDYN_TN
call sync_downup_F(F_ee2)
call sync_downup_F(F_deedt2)
call sync_downup_F(ee_past)
#endif
#endif

if (use_cfl_dt) lagran_dt = 0._rprec

end subroutine lagrange_Ssim_gpu

!*******************************************************************************
subroutine lagrange_Sdep_gpu()
!*******************************************************************************
use param
use sim_param, only : u, v, w
use sgs_param, only : F_LM, F_MM, F_QN, F_NN, beta, Cs_opt2, opftime,           &
                      lagran_dt, S11, S12, S13, S22, S23, S33, delta, ee_now,   &
                      Tn_all
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
use test_filtermodule, only : test_filter_b_gpu, test_test_filter_b_gpu
use nvtx, only : nvtxStartRange, nvtxEndRange
#ifdef PPMPI
use mpi
#endif
implicit none

integer :: jx, jy, jz
integer :: istart, iend
real(rprec) :: tf1, tf2, tf1_2, tf2_2
real(rprec) :: const, opftdelta, powcoeff
real(rprec) :: M11_loc, M12_loc, M13_loc, M22_loc, M23_loc, M33_loc
real(rprec) :: N11_loc, N12_loc, N13_loc, N22_loc, N23_loc, N33_loc
real(rprec) :: ub, vb, wb, uh, vh, wh
real(rprec) :: dumfac_loc, epsi_loc, betaclip
real(rprec), parameter :: zero = 1.e-24_rprec

opftdelta = opftime*delta
powcoeff  = -1._rprec/8._rprec
const     = 2._rprec*(delta**2)
tf1       = 2._rprec
tf2       = 4._rprec
tf1_2     = tf1**2
tf2_2     = tf2**2

call nvtxStartRange("lag_step")

! ---------------------------------------------------------------------------
! Lagrangian backwards interpolation on device. No host fall-back: F_LM,
! F_MM, F_QN, F_NN stay device-resident; u, v, w are already device-resident
! via sim_param declare-create. MPI ghost-slab sync is inside the routine.
! ---------------------------------------------------------------------------
call nvtxStartRange("interpolag")
call interpolag_Sdep_gpu()
call nvtxEndRange()

! ---------------------------------------------------------------------------
! Step 1: build u_bar/v_bar/w_bar for ALL jz in a single kernel.
! Save the same values into u_hat/v_hat/w_hat (input to test_test_filter).
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(ub, vb, wb) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    if (coord == 0 .and. jz == 1) then
        if (lbc_mom == 0) then
            ub = u(jx,jy,1)
            vb = v(jx,jy,1)
            wb = 0._rprec
        else
            ub = u(jx,jy,1)
            vb = v(jx,jy,1)
            wb = 0.25_rprec*w(jx,jy,2)
        end if
    else if (coord == nproc-1 .and. jz == nz) then
        if (ubc_mom == 0) then
            ub = u(jx,jy,nz-1)
            vb = v(jx,jy,nz-1)
            wb = 0._rprec
        else
            ub = u(jx,jy,nz-1)
            vb = v(jx,jy,nz-1)
            wb = 0.25_rprec*w(jx,jy,nz-1)
        end if
    else
        ub = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
        vb = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
        wb = w(jx,jy,jz)
    end if
    u_bar_3d(jx,jy,jz) = ub
    v_bar_3d(jx,jy,jz) = vb
    w_bar_3d(jx,jy,jz) = wb
    u_hat_3d(jx,jy,jz) = ub
    v_hat_3d(jx,jy,jz) = vb
    w_hat_3d(jx,jy,jz) = wb
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 2: Lij/Qij first term (unfiltered products). Lij = Bi*Bj.
! Qij = Lij at this point (will diverge after filtering with different
! kernel widths).
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present) private(ub, vb, wb) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    ub = u_bar_3d(jx,jy,jz)
    vb = v_bar_3d(jx,jy,jz)
    wb = w_bar_3d(jx,jy,jz)
    L11_3d(jx,jy,jz) = ub*ub
    L12_3d(jx,jy,jz) = ub*vb
    L13_3d(jx,jy,jz) = ub*wb
    L22_3d(jx,jy,jz) = vb*vb
    L23_3d(jx,jy,jz) = vb*wb
    L33_3d(jx,jy,jz) = wb*wb
    Q11_3d(jx,jy,jz) = L11_3d(jx,jy,jz)
    Q12_3d(jx,jy,jz) = L12_3d(jx,jy,jz)
    Q13_3d(jx,jy,jz) = L13_3d(jx,jy,jz)
    Q22_3d(jx,jy,jz) = L22_3d(jx,jy,jz)
    Q23_3d(jx,jy,jz) = L23_3d(jx,jy,jz)
    Q33_3d(jx,jy,jz) = L33_3d(jx,jy,jz)
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 3: test-filter (2-delta) the velocity stack and Lij stack.
! 9 batched cuFFT round-trips total.
! ---------------------------------------------------------------------------
call nvtxStartRange("filter_velL")
call test_filter_b_gpu(u_bar_3d, nz)
call test_filter_b_gpu(v_bar_3d, nz)
call test_filter_b_gpu(w_bar_3d, nz)
call test_filter_b_gpu(L11_3d,  nz)
call test_filter_b_gpu(L12_3d,  nz)
call test_filter_b_gpu(L13_3d,  nz)
call test_filter_b_gpu(L22_3d,  nz)
call test_filter_b_gpu(L23_3d,  nz)
call test_filter_b_gpu(L33_3d,  nz)
call nvtxEndRange()

! Lij = Lij_filtered - filtered_ui * filtered_uj
!$acc parallel loop collapse(3) default(present) private(ub, vb, wb) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    ub = u_bar_3d(jx,jy,jz)
    vb = v_bar_3d(jx,jy,jz)
    wb = w_bar_3d(jx,jy,jz)
    L11_3d(jx,jy,jz) = L11_3d(jx,jy,jz) - ub*ub
    L12_3d(jx,jy,jz) = L12_3d(jx,jy,jz) - ub*vb
    L13_3d(jx,jy,jz) = L13_3d(jx,jy,jz) - ub*wb
    L22_3d(jx,jy,jz) = L22_3d(jx,jy,jz) - vb*vb
    L23_3d(jx,jy,jz) = L23_3d(jx,jy,jz) - vb*wb
    L33_3d(jx,jy,jz) = L33_3d(jx,jy,jz) - wb*wb
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 4: test-test-filter (4-delta) the velocity-hat stack and Qij stack.
! ---------------------------------------------------------------------------
call nvtxStartRange("filter_velQ")
call test_test_filter_b_gpu(u_hat_3d, nz)
call test_test_filter_b_gpu(v_hat_3d, nz)
call test_test_filter_b_gpu(w_hat_3d, nz)
call test_test_filter_b_gpu(Q11_3d,  nz)
call test_test_filter_b_gpu(Q12_3d,  nz)
call test_test_filter_b_gpu(Q13_3d,  nz)
call test_test_filter_b_gpu(Q22_3d,  nz)
call test_test_filter_b_gpu(Q23_3d,  nz)
call test_test_filter_b_gpu(Q33_3d,  nz)
call nvtxEndRange()

!$acc parallel loop collapse(3) default(present) private(uh, vh, wh) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    uh = u_hat_3d(jx,jy,jz)
    vh = v_hat_3d(jx,jy,jz)
    wh = w_hat_3d(jx,jy,jz)
    Q11_3d(jx,jy,jz) = Q11_3d(jx,jy,jz) - uh*uh
    Q12_3d(jx,jy,jz) = Q12_3d(jx,jy,jz) - uh*vh
    Q13_3d(jx,jy,jz) = Q13_3d(jx,jy,jz) - uh*wh
    Q22_3d(jx,jy,jz) = Q22_3d(jx,jy,jz) - vh*vh
    Q23_3d(jx,jy,jz) = Q23_3d(jx,jy,jz) - vh*wh
    Q33_3d(jx,jy,jz) = Q33_3d(jx,jy,jz) - wh*wh
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 5: |S| and copy Sij(:,:,jz) into the per-jz _bar / _hat scratch
! ready for filtering.
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    S_3d(jx,jy,jz) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2         &
                       + S33(jx,jy,jz)**2                                       &
                       + 2._rprec*(S12(jx,jy,jz)**2 + S13(jx,jy,jz)**2          &
                                  + S23(jx,jy,jz)**2)))
    S11_bar_3d(jx,jy,jz) = S11(jx,jy,jz)
    S12_bar_3d(jx,jy,jz) = S12(jx,jy,jz)
    S13_bar_3d(jx,jy,jz) = S13(jx,jy,jz)
    S22_bar_3d(jx,jy,jz) = S22(jx,jy,jz)
    S23_bar_3d(jx,jy,jz) = S23(jx,jy,jz)
    S33_bar_3d(jx,jy,jz) = S33(jx,jy,jz)
    S11_hat_3d(jx,jy,jz) = S11(jx,jy,jz)
    S12_hat_3d(jx,jy,jz) = S12(jx,jy,jz)
    S13_hat_3d(jx,jy,jz) = S13(jx,jy,jz)
    S22_hat_3d(jx,jy,jz) = S22(jx,jy,jz)
    S23_hat_3d(jx,jy,jz) = S23(jx,jy,jz)
    S33_hat_3d(jx,jy,jz) = S33(jx,jy,jz)
end do
end do
end do

call nvtxStartRange("filter_Sbar")
call test_filter_b_gpu(S11_bar_3d, nz)
call test_filter_b_gpu(S12_bar_3d, nz)
call test_filter_b_gpu(S13_bar_3d, nz)
call test_filter_b_gpu(S22_bar_3d, nz)
call test_filter_b_gpu(S23_bar_3d, nz)
call test_filter_b_gpu(S33_bar_3d, nz)
call nvtxEndRange()

call nvtxStartRange("filter_Shat")
call test_test_filter_b_gpu(S11_hat_3d, nz)
call test_test_filter_b_gpu(S12_hat_3d, nz)
call test_test_filter_b_gpu(S13_hat_3d, nz)
call test_test_filter_b_gpu(S22_hat_3d, nz)
call test_test_filter_b_gpu(S23_hat_3d, nz)
call test_test_filter_b_gpu(S33_hat_3d, nz)
call nvtxEndRange()

! ---------------------------------------------------------------------------
! Step 6: |S_bar|, |S_hat|, and S_S{ij}_bar = |S| * S{ij},
!         S_S{ij}_hat = same (input to second filter)
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    S_bar_3d(jx,jy,jz) = sqrt(2._rprec*(S11_bar_3d(jx,jy,jz)**2                 &
                       + S22_bar_3d(jx,jy,jz)**2 + S33_bar_3d(jx,jy,jz)**2      &
                       + 2._rprec*(S12_bar_3d(jx,jy,jz)**2                      &
                                + S13_bar_3d(jx,jy,jz)**2                       &
                                + S23_bar_3d(jx,jy,jz)**2)))
    S_hat_3d(jx,jy,jz) = sqrt(2._rprec*(S11_hat_3d(jx,jy,jz)**2                 &
                       + S22_hat_3d(jx,jy,jz)**2 + S33_hat_3d(jx,jy,jz)**2      &
                       + 2._rprec*(S12_hat_3d(jx,jy,jz)**2                      &
                                + S13_hat_3d(jx,jy,jz)**2                       &
                                + S23_hat_3d(jx,jy,jz)**2)))

    S_S11_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S11(jx,jy,jz)
    S_S12_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S12(jx,jy,jz)
    S_S13_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S13(jx,jy,jz)
    S_S22_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S22(jx,jy,jz)
    S_S23_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S23(jx,jy,jz)
    S_S33_bar_3d(jx,jy,jz) = S_3d(jx,jy,jz)*S33(jx,jy,jz)

    S_S11_hat_3d(jx,jy,jz) = S_S11_bar_3d(jx,jy,jz)
    S_S12_hat_3d(jx,jy,jz) = S_S12_bar_3d(jx,jy,jz)
    S_S13_hat_3d(jx,jy,jz) = S_S13_bar_3d(jx,jy,jz)
    S_S22_hat_3d(jx,jy,jz) = S_S22_bar_3d(jx,jy,jz)
    S_S23_hat_3d(jx,jy,jz) = S_S23_bar_3d(jx,jy,jz)
    S_S33_hat_3d(jx,jy,jz) = S_S33_bar_3d(jx,jy,jz)
end do
end do
end do

call nvtxStartRange("filter_SSbar")
call test_filter_b_gpu(S_S11_bar_3d, nz)
call test_filter_b_gpu(S_S12_bar_3d, nz)
call test_filter_b_gpu(S_S13_bar_3d, nz)
call test_filter_b_gpu(S_S22_bar_3d, nz)
call test_filter_b_gpu(S_S23_bar_3d, nz)
call test_filter_b_gpu(S_S33_bar_3d, nz)
call nvtxEndRange()

call nvtxStartRange("filter_SShat")
call test_test_filter_b_gpu(S_S11_hat_3d, nz)
call test_test_filter_b_gpu(S_S12_hat_3d, nz)
call test_test_filter_b_gpu(S_S13_hat_3d, nz)
call test_test_filter_b_gpu(S_S22_hat_3d, nz)
call test_test_filter_b_gpu(S_S23_hat_3d, nz)
call test_test_filter_b_gpu(S_S33_hat_3d, nz)
call nvtxEndRange()

! ---------------------------------------------------------------------------
! Step 7: form Mij, Nij as per-(i,j,k) scalars; produce LM/MM/QN/NN/ee_now.
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(M11_loc, M12_loc, M13_loc, M22_loc, M23_loc, M33_loc,    &
!$acc                  N11_loc, N12_loc, N13_loc, N22_loc, N23_loc, N33_loc) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    M11_loc = const*(S_S11_bar_3d(jx,jy,jz) - tf1_2*S_bar_3d(jx,jy,jz)*S11_bar_3d(jx,jy,jz))
    M12_loc = const*(S_S12_bar_3d(jx,jy,jz) - tf1_2*S_bar_3d(jx,jy,jz)*S12_bar_3d(jx,jy,jz))
    M13_loc = const*(S_S13_bar_3d(jx,jy,jz) - tf1_2*S_bar_3d(jx,jy,jz)*S13_bar_3d(jx,jy,jz))
    M22_loc = const*(S_S22_bar_3d(jx,jy,jz) - tf1_2*S_bar_3d(jx,jy,jz)*S22_bar_3d(jx,jy,jz))
    M23_loc = const*(S_S23_bar_3d(jx,jy,jz) - tf1_2*S_bar_3d(jx,jy,jz)*S23_bar_3d(jx,jy,jz))
    M33_loc = const*(S_S33_bar_3d(jx,jy,jz) - tf1_2*S_bar_3d(jx,jy,jz)*S33_bar_3d(jx,jy,jz))

    N11_loc = const*(S_S11_hat_3d(jx,jy,jz) - tf2_2*S_hat_3d(jx,jy,jz)*S11_hat_3d(jx,jy,jz))
    N12_loc = const*(S_S12_hat_3d(jx,jy,jz) - tf2_2*S_hat_3d(jx,jy,jz)*S12_hat_3d(jx,jy,jz))
    N13_loc = const*(S_S13_hat_3d(jx,jy,jz) - tf2_2*S_hat_3d(jx,jy,jz)*S13_hat_3d(jx,jy,jz))
    N22_loc = const*(S_S22_hat_3d(jx,jy,jz) - tf2_2*S_hat_3d(jx,jy,jz)*S22_hat_3d(jx,jy,jz))
    N23_loc = const*(S_S23_hat_3d(jx,jy,jz) - tf2_2*S_hat_3d(jx,jy,jz)*S23_hat_3d(jx,jy,jz))
    N33_loc = const*(S_S33_hat_3d(jx,jy,jz) - tf2_2*S_hat_3d(jx,jy,jz)*S33_hat_3d(jx,jy,jz))

    LM_3d(jx,jy,jz) = L11_3d(jx,jy,jz)*M11_loc + L22_3d(jx,jy,jz)*M22_loc       &
                    + L33_3d(jx,jy,jz)*M33_loc                                  &
                    + 2._rprec*(L12_3d(jx,jy,jz)*M12_loc                        &
                              + L13_3d(jx,jy,jz)*M13_loc                        &
                              + L23_3d(jx,jy,jz)*M23_loc)
    MM_3d(jx,jy,jz) = M11_loc**2 + M22_loc**2 + M33_loc**2                      &
                    + 2._rprec*(M12_loc**2 + M13_loc**2 + M23_loc**2)
    QN_3d(jx,jy,jz) = Q11_3d(jx,jy,jz)*N11_loc + Q22_3d(jx,jy,jz)*N22_loc       &
                    + Q33_3d(jx,jy,jz)*N33_loc                                  &
                    + 2._rprec*(Q12_3d(jx,jy,jz)*N12_loc                        &
                              + Q13_3d(jx,jy,jz)*N13_loc                        &
                              + Q23_3d(jx,jy,jz)*N23_loc)
    NN_3d(jx,jy,jz) = N11_loc**2 + N22_loc**2 + N33_loc**2                      &
                    + 2._rprec*(N12_loc**2 + N13_loc**2 + N23_loc**2)

    ee_now(jx,jy,jz) = L11_3d(jx,jy,jz)**2 + L22_3d(jx,jy,jz)**2                &
                     + L33_3d(jx,jy,jz)**2                                      &
                     + 2._rprec*(L12_3d(jx,jy,jz)**2 + L13_3d(jx,jy,jz)**2      &
                               + L23_3d(jx,jy,jz)**2)                           &
                     - 2._rprec*LM_3d(jx,jy,jz)*Cs_opt2(jx,jy,jz)               &
                     + MM_3d(jx,jy,jz)*Cs_opt2(jx,jy,jz)**2
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 8: One-shot init of F_LM/F_MM/F_QN/F_NN at the first qualifying step
! (jt == cs_count or jt == DYN_init). Mirrors the CPU `inilag` branch but
! batched across all jz at once. F_LM_MM_init / F_QN_NN_init are
! module-private save logicals, set true after the first init.
! ---------------------------------------------------------------------------
if (inilag .and. (.not. F_LM_MM_init) .and.                                     &
    (jt == cs_count .or. jt == DYN_init)) then
    if (coord == 0) print *, 'F_MM and F_LM initialized (GPU)'
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, ld
        F_MM(jx,jy,jz) = MM_3d(jx,jy,jz)
        F_LM(jx,jy,jz) = 0.03_rprec*MM_3d(jx,jy,jz)
    end do
    end do
    end do
    !$acc parallel loop collapse(2) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
        F_MM(ld-1, jy, jz) = 1._rprec
        F_MM(ld,   jy, jz) = 1._rprec
        F_LM(ld-1, jy, jz) = 1._rprec
        F_LM(ld,   jy, jz) = 1._rprec
    end do
    end do
    F_LM_MM_init = .true.
end if

! ---------------------------------------------------------------------------
! Step 9: Inflow zero-out (rare). Done with two batched kernels.
! ---------------------------------------------------------------------------
if (inflow_type > 0) then
    iend = floor (fringe_region_end * nx + 1._rprec)
    iend = modulo (iend - 1, nx) + 1
    istart = floor ((fringe_region_end - fringe_region_len) * nx + 1._rprec)
    istart = modulo (istart - 1, nx) + 1

    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, ld
        if (MM_3d(jx,jy,jz) <= 0.1_rprec*const*S_3d(jx,jy,jz)**2) then
            MM_3d(jx,jy,jz) = 0.1_rprec*const*S_3d(jx,jy,jz)**2
        end if
        if (NN_3d(jx,jy,jz) <= 0.1_rprec*const*S_3d(jx,jy,jz)**2) then
            NN_3d(jx,jy,jz) = 0.1_rprec*const*S_3d(jx,jy,jz)**2
        end if
    end do
    end do
    end do

    ! Guard against empty range: when fringe_region_end=1.0 the modulo trick
    ! above wraps iend to 1 while istart stays at 169, giving jx=170..1. On
    ! CPU the equivalent array section LM(170:1,:) is a no-op zero-size slice;
    ! on GPU collapse(3) the negative trip count is reinterpreted as ~4 billion
    ! threads and the trivial kernel runs for ~1.6 s/Lag step (was 91% of total
    ! GPU time before this guard was added).
    if (iend >= istart + 1) then
        !$acc parallel loop collapse(3) default(present) async(1)
        do jz = 1, nz
        do jy = 1, ny
        do jx = istart+1, iend
            LM_3d(jx,jy,jz) = 0._rprec
            F_LM(jx,jy,jz)  = 0._rprec
            QN_3d(jx,jy,jz) = 0._rprec
            F_QN(jx,jy,jz)  = 0._rprec
        end do
        end do
        end do
    end if
end if

! ---------------------------------------------------------------------------
! Step 10: Lagrangian timescale Tn_3d (Meneveau, Lund, Cabot 1996),
! running-average update of F_LM/F_MM, and Cs_2d_3d.
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(dumfac_loc, epsi_loc) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
#ifdef PPDYN_TN
    Tn_3d(jx,jy,jz) = 4._rprec*pi*sqrt(F_ee2(jx,jy,jz)/F_deedt2(jx,jy,jz))
#else
    Tn_3d(jx,jy,jz) = max(F_LM(jx,jy,jz)*F_MM(jx,jy,jz), zero)
    Tn_3d(jx,jy,jz) = opftdelta*(Tn_3d(jx,jy,jz)**powcoeff)
    Tn_3d(jx,jy,jz) = max(zero, Tn_3d(jx,jy,jz))
#endif

    dumfac_loc = lagran_dt/Tn_3d(jx,jy,jz)
    epsi_loc   = dumfac_loc/(1._rprec + dumfac_loc)

    F_LM(jx,jy,jz) = epsi_loc*LM_3d(jx,jy,jz)                                   &
                   + (1._rprec - epsi_loc)*F_LM(jx,jy,jz)
    F_MM(jx,jy,jz) = epsi_loc*MM_3d(jx,jy,jz)                                   &
                   + (1._rprec - epsi_loc)*F_MM(jx,jy,jz)
    F_LM(jx,jy,jz) = max(zero, F_LM(jx,jy,jz))

    ! Fuse padding and clipping into the producer kernel to avoid two
    ! extra cleanup launches per dynamic SGS update.
    if (jx >= ld-1) then
        Cs_2d_3d(jx,jy,jz) = zero
    else
        Cs_2d_3d(jx,jy,jz) = max(zero, F_LM(jx,jy,jz)/(F_MM(jx,jy,jz) + zero))
    end if
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 11: One-shot init of F_QN/F_NN at the first qualifying step.
! ---------------------------------------------------------------------------
if (inilag .and. (.not. F_QN_NN_init) .and.                                     &
    (jt == cs_count .or. jt == DYN_init)) then
    if (coord == 0) print *, 'F_NN and F_QN initialized (GPU)'
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, ld
        F_NN(jx,jy,jz) = NN_3d(jx,jy,jz)
        F_QN(jx,jy,jz) = 0.03_rprec*NN_3d(jx,jy,jz)
    end do
    end do
    end do
    !$acc parallel loop collapse(2) default(present) async(1)
    do jz = 1, nz
    do jy = 1, ny
        F_NN(ld-1, jy, jz) = 1._rprec
        F_NN(ld,   jy, jz) = 1._rprec
        F_QN(ld-1, jy, jz) = 1._rprec
        F_QN(ld,   jy, jz) = 1._rprec
    end do
    end do
    F_QN_NN_init = .true.
end if

! ---------------------------------------------------------------------------
! Step 12: Tn for 4-delta, F_QN/F_NN running average, Cs_4d_3d
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(dumfac_loc, epsi_loc) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
#ifndef PPDYN_TN
    Tn_3d(jx,jy,jz) = max(F_QN(jx,jy,jz)*F_NN(jx,jy,jz), zero)
    Tn_3d(jx,jy,jz) = opftdelta*(Tn_3d(jx,jy,jz)**powcoeff)
    Tn_3d(jx,jy,jz) = max(zero, Tn_3d(jx,jy,jz))
#endif

    dumfac_loc = lagran_dt/Tn_3d(jx,jy,jz)
    epsi_loc   = dumfac_loc/(1._rprec + dumfac_loc)

    F_QN(jx,jy,jz) = epsi_loc*QN_3d(jx,jy,jz)                                   &
                   + (1._rprec - epsi_loc)*F_QN(jx,jy,jz)
    F_NN(jx,jy,jz) = epsi_loc*NN_3d(jx,jy,jz)                                   &
                   + (1._rprec - epsi_loc)*F_NN(jx,jy,jz)
    F_QN(jx,jy,jz) = max(zero, F_QN(jx,jy,jz))

#ifdef PPDYN_TN
    F_ee2(jx,jy,jz) = epsi_loc*ee_now(jx,jy,jz)**2                             &
                    + (1._rprec - epsi_loc)*F_ee2(jx,jy,jz)
    F_deedt2(jx,jy,jz) = epsi_loc*                                            &
        (((ee_now(jx,jy,jz) - ee_past(jx,jy,jz))/lagran_dt)**2)                &
                       + (1._rprec - epsi_loc)*F_deedt2(jx,jy,jz)
    ee_past(jx,jy,jz) = ee_now(jx,jy,jz)
#endif

    ! Same fused pad/clip as Cs_2d_3d.
    if (jx >= ld-1) then
        Cs_4d_3d(jx,jy,jz) = zero
    else
        Cs_4d_3d(jx,jy,jz) = max(zero, F_QN(jx,jy,jz)/(F_NN(jx,jy,jz) + zero))
    end if
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 13: Beta(:,:,jz) = (Cs_4d/Cs_2d)^(...) with BC overrides at top/bot
! walls (only for stress-free; lbc/ubc_mom == 0). Then Cs_opt2 from
! Cs_2d / max(Beta, 1/(tf1*tf2)).
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    Beta(jx,jy,jz) = (Cs_4d_3d(jx,jy,jz)/Cs_2d_3d(jx,jy,jz))                    &
                     **(log(tf1)/(log(tf2) - log(tf1)))
end do
end do
end do

#ifdef PPMPI
if (coord == nproc-1 .and. ubc_mom == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
    do jx = 1, ld
        Beta(jx,jy,nz) = 1._rprec
    end do
    end do
end if
#else
if (ubc_mom == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
    do jx = 1, ld
        Beta(jx,jy,nz) = 1._rprec
    end do
    end do
end if
#endif
if (coord == 0 .and. lbc_mom == 0) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
    do jx = 1, ld
        Beta(jx,jy,1) = 1._rprec
    end do
    end do
end if

!$acc parallel loop collapse(3) default(present) private(betaclip) async(1)
do jz = 1, nz
do jy = 1, ny
do jx = 1, ld
    if (jx >= ld-1) then
        Cs_opt2(jx,jy,jz) = zero
    else
        betaclip = max(Beta(jx,jy,jz), 1._rprec/(tf1*tf2))
        Cs_opt2(jx,jy,jz) = max(zero, Cs_2d_3d(jx,jy,jz)/betaclip)
    end if
    Tn_all(jx,jy,jz)  = Tn_3d(jx,jy,jz)
end do
end do
end do

! ---------------------------------------------------------------------------
! Cross-rank sync of F_LM/F_MM/F_QN/F_NN/Tn_all (ghost slabs k=0 and k=nz).
! ---------------------------------------------------------------------------
#ifdef PPMPI
call nvtxStartRange("sync_F")
call sync_downup_F(F_LM)
call sync_downup_F(F_MM)
call sync_downup_F(F_QN)
call sync_downup_F(F_NN)
call sync_downup_F(Tn_all)
#ifdef PPDYN_TN
call sync_downup_F(F_ee2)
call sync_downup_F(F_deedt2)
call sync_downup_F(ee_past)
#endif
call nvtxEndRange()
#endif

! Reset variable for use during next set of cs_count timesteps
if (use_cfl_dt) lagran_dt = 0._rprec

call nvtxEndRange()    ! pairs with nvtxStartRange("lag_step") at top

end subroutine lagrange_Sdep_gpu

!*******************************************************************************
subroutine interpolag_Ssim_gpu()
!*******************************************************************************
! GPU port of interpolag_Ssim. Updates only F_LM/F_MM and preserves the
! model-4 bottom/top boundary behavior from interpolag_Ssim.f90.
!*******************************************************************************
use param, only : ld, nx, ny, nz, lbz, dx, dy, dz, coord, nproc,               &
                  L_x, L_y, L_z, jt_total, lag_cfl_count, dt
use sgs_param, only : F_LM, F_MM, lagran_dt
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
use sim_param, only : u, v, w
use cfl_util, only : get_max_cfl
implicit none

integer :: i, j, k, i1, i2, j1, j2, k1, k2
real(rprec) :: px, py, pz, xdiff, ydiff, zdiff
real(rprec) :: wxgt, wygt, wzgt
real(rprec) :: c000, c100, c010, c110, c001, c101, c011, c111
real(rprec) :: zloc, ztop, lcfl
real(rprec), parameter :: thresh = 1.e-9_rprec

! Snapshot old averages before in-place interpolation.
!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, ld
    tempF_LM(i,j,k) = F_LM(i,j,k)
    tempF_MM(i,j,k) = F_MM(i,j,k)
#ifdef PPDYN_TN
    tempF_ee2(i,j,k) = F_ee2(i,j,k)
    tempF_deedt2(i,j,k) = F_deedt2(i,j,k)
    tempee_past(i,j,k) = ee_past(i,j,k)
#endif
end do
end do
end do

! Interpolate all model-4 levels that CPU interpolag_Ssim updates.
! coord==0,k==1 is intentionally left unchanged.
! Non-top ranks intentionally leave k==nz for MPI overlap.
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(px, py, pz, xdiff, ydiff, zdiff, i1, i2, j1, j2,        &
!$acc                  k1, k2, wxgt, wygt, wzgt, c000, c100, c010, c110,       &
!$acc                  c001, c101, c011, c111, zloc, ztop) async(1)
do k = 1, nz
do j = 1, ny
do i = 1, nx
    if (((coord /= 0) .or. (k >= 2)) .and.                                    &
        ((k <= nz-1) .or. ((coord == nproc-1) .and. (k == nz)))) then
        px = real(i - 1, rprec)*dx -                                          &
            0.5_rprec*(u(i,j,k-1) + u(i,j,k))*lagran_dt
        py = real(j - 1, rprec)*dy -                                          &
            0.5_rprec*(v(i,j,k-1) + v(i,j,k))*lagran_dt
        if ((coord == nproc-1) .and. (k == nz)) then
            pz = real(coord*(nz-1), rprec)*dz + real(k, rprec)*dz -           &
                0.5_rprec*dz - max(0._rprec, w(i,j,k))*lagran_dt
        else
            pz = real(coord*(nz-1), rprec)*dz + real(k, rprec)*dz -           &
                0.5_rprec*dz - w(i,j,k)*lagran_dt
        end if

        px = modulo(px, L_x)
        if (abs(px)/L_x < thresh) then
            i1 = 1
        else if (abs(px - L_x)/L_x < thresh) then
            i1 = nx
        else
            i1 = floor(px/dx) + 1
        end if
        i1 = min(max(i1, 1), nx)
        xdiff = px - real(i1 - 1, rprec)*dx

        py = modulo(py, L_y)
        if (abs(py)/L_y < thresh) then
            j1 = 1
        else if (abs(py - L_y)/L_y < thresh) then
            j1 = ny
        else
            j1 = floor(py/dy) + 1
        end if
        j1 = min(max(j1, 1), ny)
        ydiff = py - real(j1 - 1, rprec)*dy

        zloc = real(coord*(nz-1), rprec)*dz + 0.5_rprec*dz
        ztop = real(coord*(nz-1), rprec)*dz + real(nz, rprec)*dz -            &
            0.5_rprec*dz
        if (abs(pz - ztop)/L_z < thresh) then
            k1 = nz - 1
        else
            k1 = floor((pz - zloc)/dz) + 1
        end if
        k1 = min(max(k1, lbz), nz-1)
        zdiff = pz - (real(coord*(nz-1), rprec)*dz +                          &
            real(k1, rprec)*dz - 0.5_rprec*dz)

        i2 = i1 + 1
        if (i2 > nx) i2 = 1
        j2 = j1 + 1
        if (j2 > ny) j2 = 1
        k2 = k1 + 1

        wxgt = xdiff/dx
        wygt = ydiff/dy
        wzgt = zdiff/dz
        c000 = (1._rprec-wxgt)*(1._rprec-wygt)*(1._rprec-wzgt)
        c100 = wxgt*(1._rprec-wygt)*(1._rprec-wzgt)
        c010 = (1._rprec-wxgt)*wygt*(1._rprec-wzgt)
        c110 = wxgt*wygt*(1._rprec-wzgt)
        c001 = (1._rprec-wxgt)*(1._rprec-wygt)*wzgt
        c101 = wxgt*(1._rprec-wygt)*wzgt
        c011 = (1._rprec-wxgt)*wygt*wzgt
        c111 = wxgt*wygt*wzgt

        F_LM(i,j,k) = c000*tempF_LM(i1,j1,k1) + c100*tempF_LM(i2,j1,k1) +     &
            c010*tempF_LM(i1,j2,k1) + c110*tempF_LM(i2,j2,k1) +               &
            c001*tempF_LM(i1,j1,k2) + c101*tempF_LM(i2,j1,k2) +               &
            c011*tempF_LM(i1,j2,k2) + c111*tempF_LM(i2,j2,k2)
        F_MM(i,j,k) = c000*tempF_MM(i1,j1,k1) + c100*tempF_MM(i2,j1,k1) +     &
            c010*tempF_MM(i1,j2,k1) + c110*tempF_MM(i2,j2,k1) +               &
            c001*tempF_MM(i1,j1,k2) + c101*tempF_MM(i2,j1,k2) +               &
            c011*tempF_MM(i1,j2,k2) + c111*tempF_MM(i2,j2,k2)
#ifdef PPDYN_TN
        F_ee2(i,j,k) = c000*tempF_ee2(i1,j1,k1) +                             &
            c100*tempF_ee2(i2,j1,k1) + c010*tempF_ee2(i1,j2,k1) +             &
            c110*tempF_ee2(i2,j2,k1) + c001*tempF_ee2(i1,j1,k2) +             &
            c101*tempF_ee2(i2,j1,k2) + c011*tempF_ee2(i1,j2,k2) +             &
            c111*tempF_ee2(i2,j2,k2)
        F_deedt2(i,j,k) = c000*tempF_deedt2(i1,j1,k1) +                       &
            c100*tempF_deedt2(i2,j1,k1) + c010*tempF_deedt2(i1,j2,k1) +       &
            c110*tempF_deedt2(i2,j2,k1) + c001*tempF_deedt2(i1,j1,k2) +       &
            c101*tempF_deedt2(i2,j1,k2) + c011*tempF_deedt2(i1,j2,k2) +       &
            c111*tempF_deedt2(i2,j2,k2)
        ee_past(i,j,k) = c000*tempee_past(i1,j1,k1) +                         &
            c100*tempee_past(i2,j1,k1) + c010*tempee_past(i1,j2,k1) +         &
            c110*tempee_past(i2,j2,k1) + c001*tempee_past(i1,j1,k2) +         &
            c101*tempee_past(i2,j1,k2) + c011*tempee_past(i1,j2,k2) +         &
            c111*tempee_past(i2,j2,k2)
#endif
    end if
end do
end do
end do

#ifdef PPMPI
call sync_downup_F(F_LM)
call sync_downup_F(F_MM)
#ifdef PPDYN_TN
call sync_downup_F(F_ee2)
call sync_downup_F(F_deedt2)
call sync_downup_F(ee_past)
#endif
#endif

if (mod(jt_total, lag_cfl_count) == 0) then
    lcfl = get_max_cfl()
    lcfl = lcfl*lagran_dt/dt
#ifdef PPMPI
    if (coord == 0) print *, 'Lagrangian CFL condition= ', lcfl
#else
    print *, 'Lagrangian CFL condition= ', lcfl
#endif
end if

end subroutine interpolag_Ssim_gpu

!*******************************************************************************
subroutine interpolag_Sdep_gpu()
!*******************************************************************************
! GPU port of interpolag_Sdep. For every (i, j, k) in the F_* grid, compute
! the previous-timestep position xyz_past = (x - u*dt, y - v*dt, z - w*dt),
! locate it in the w-grid (cell_indx_w autowrap + diffs), and trilinearly
! interpolate the OLD F_* (snap-shotted into tempF_*) at xyz_past. The
! result overwrites F_* at (i, j, k). One kernel per region (bottom /
! intermediate / top) - the cell-index math is computed once per (i,j,k)
! and reused across all four F_* interpolations.
!
! Why inlined (no `!$acc routine seq`): NVHPC 25.5 won't generate a device
! variant of a function that references module-level scalar param vars
! (warning W-1054 about ld). Inlining avoids that and also fuses the four
! tempF_* reads into the same memory-locality region.
!
! Performance: replaces the host fall-back path in lagrange_Sdep_gpu which
! had been doing `!$acc update self(u, v, w, F_LM, F_MM, F_QN, F_NN)`
! (~250 MB pull) + CPU interpolag (~1-1.5 s) + push back. Now it's three
! collapse kernels on device data + a sync_downup_F per F_*.
!*******************************************************************************
use param, only : ld, nx, ny, nz, lbz, dx, dy, dz, coord, nproc, lbc_mom,       &
                  ubc_mom, L_x, L_y, L_z, jt_total, lag_cfl_count, dt
use sgs_param, only : F_LM, F_MM, F_QN, F_NN, lagran_dt
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
use sim_param, only : u, v, w
use cfl_util,  only : get_max_cfl
implicit none

integer :: i, j, k, kmin
integer :: istart, jstart, kstart, istart1, jstart1, kstart1
real(rprec) :: px, py, pz, xdiff, ydiff, zdiff
real(rprec) :: u1, u2, u3, u4, u5, u6
real(rprec) :: lcfl
real(rprec), parameter :: thresh = 1.e-9_rprec

! Step 1: snapshot F_* into tempF_* so trilinear-interp reads see the
! pre-update values while kernels write F_* in place.
!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, ld
    tempF_LM(i,j,k) = F_LM(i,j,k)
    tempF_MM(i,j,k) = F_MM(i,j,k)
    tempF_QN(i,j,k) = F_QN(i,j,k)
    tempF_NN(i,j,k) = F_NN(i,j,k)
#ifdef PPDYN_TN
    tempF_ee2(i,j,k)    = F_ee2(i,j,k)
    tempF_deedt2(i,j,k) = F_deedt2(i,j,k)
    tempee_past(i,j,k)  = ee_past(i,j,k)
#endif
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 2: bottom slab (coord==0, k=1).
! ---------------------------------------------------------------------------
if (coord == 0) then
    k = 1
    if (lbc_mom == 0) then
        ! Stress-free: u/v from neighbour, z stays on zw-node (no penetration)
        !$acc parallel loop collapse(2) default(present)                        &
        !$acc          private(px, py, pz, istart, jstart, kstart,              &
        !$acc                  istart1, jstart1, kstart1,                       &
        !$acc                  xdiff, ydiff, zdiff,                             &
        !$acc                  u1, u2, u3, u4, u5, u6) async(1)
        do j = 1, ny
        do i = 1, nx
            px = x_dev(i)  - u(i,j,k)*lagran_dt
            py = y_dev(j)  - v(i,j,k)*lagran_dt
            pz = zw_dev(k)
            ! ---- locate (istart, jstart) on uvp grid ----
            px = modulo(px, L_x)
            if (abs(px) / L_x < thresh) then
                istart = 1
            else if (abs(px - L_x) / L_x < thresh) then
                istart = nx
            else
                istart = floor(px / dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py) / L_y < thresh) then
                jstart = 1
            else if (abs(py - L_y) / L_y < thresh) then
                jstart = ny
            else
                jstart = floor(py / dy) + 1
            end if
            istart1 = autowrap_i_dev(istart + 1)
            jstart1 = autowrap_j_dev(jstart + 1)
            xdiff = px - x_dev(istart)
            ydiff = py - y_dev(jstart)
            ! ---- locate kstart on w-grid (special-case branches) ----
            ! lbc_mom == 0 here, so the "wall" branch (coord==0 + lbc_mom>0 +
            ! pz<zw(2)) cannot fire. Use the generic z-branch directly.
            if (abs(pz - zw_dev(nz)) / L_z < thresh) then
                kstart = nz-1
            else
                kstart = floor((pz - zw_dev(1)) / dz) + 1
            end if
            kstart1 = kstart + 1
            zdiff   = pz - zw_dev(kstart)
            ! ---- interp F_LM ----
            u1 = tempF_LM(istart,  jstart,  kstart)                             &
               + xdiff*(tempF_LM(istart1, jstart,  kstart)                      &
                      - tempF_LM(istart,  jstart,  kstart))/dx
            u2 = tempF_LM(istart,  jstart1, kstart)                             &
               + xdiff*(tempF_LM(istart1, jstart1, kstart)                      &
                      - tempF_LM(istart,  jstart1, kstart))/dx
            u3 = tempF_LM(istart,  jstart,  kstart1)                            &
               + xdiff*(tempF_LM(istart1, jstart,  kstart1)                     &
                      - tempF_LM(istart,  jstart,  kstart1))/dx
            u4 = tempF_LM(istart,  jstart1, kstart1)                            &
               + xdiff*(tempF_LM(istart1, jstart1, kstart1)                     &
                      - tempF_LM(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            F_LM(i,j,k) = u5 + zdiff*(u6 - u5)/dz
            ! ---- interp F_MM ----
            u1 = tempF_MM(istart,  jstart,  kstart)                             &
               + xdiff*(tempF_MM(istart1, jstart,  kstart)                      &
                      - tempF_MM(istart,  jstart,  kstart))/dx
            u2 = tempF_MM(istart,  jstart1, kstart)                             &
               + xdiff*(tempF_MM(istart1, jstart1, kstart)                      &
                      - tempF_MM(istart,  jstart1, kstart))/dx
            u3 = tempF_MM(istart,  jstart,  kstart1)                            &
               + xdiff*(tempF_MM(istart1, jstart,  kstart1)                     &
                      - tempF_MM(istart,  jstart,  kstart1))/dx
            u4 = tempF_MM(istart,  jstart1, kstart1)                            &
               + xdiff*(tempF_MM(istart1, jstart1, kstart1)                     &
                      - tempF_MM(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            F_MM(i,j,k) = u5 + zdiff*(u6 - u5)/dz
            ! ---- interp F_QN ----
            u1 = tempF_QN(istart,  jstart,  kstart)                             &
               + xdiff*(tempF_QN(istart1, jstart,  kstart)                      &
                      - tempF_QN(istart,  jstart,  kstart))/dx
            u2 = tempF_QN(istart,  jstart1, kstart)                             &
               + xdiff*(tempF_QN(istart1, jstart1, kstart)                      &
                      - tempF_QN(istart,  jstart1, kstart))/dx
            u3 = tempF_QN(istart,  jstart,  kstart1)                            &
               + xdiff*(tempF_QN(istart1, jstart,  kstart1)                     &
                      - tempF_QN(istart,  jstart,  kstart1))/dx
            u4 = tempF_QN(istart,  jstart1, kstart1)                            &
               + xdiff*(tempF_QN(istart1, jstart1, kstart1)                     &
                      - tempF_QN(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            F_QN(i,j,k) = u5 + zdiff*(u6 - u5)/dz
            ! ---- interp F_NN ----
            u1 = tempF_NN(istart,  jstart,  kstart)                             &
               + xdiff*(tempF_NN(istart1, jstart,  kstart)                      &
                      - tempF_NN(istart,  jstart,  kstart))/dx
            u2 = tempF_NN(istart,  jstart1, kstart)                             &
               + xdiff*(tempF_NN(istart1, jstart1, kstart)                      &
                      - tempF_NN(istart,  jstart1, kstart))/dx
            u3 = tempF_NN(istart,  jstart,  kstart1)                            &
               + xdiff*(tempF_NN(istart1, jstart,  kstart1)                     &
                      - tempF_NN(istart,  jstart,  kstart1))/dx
            u4 = tempF_NN(istart,  jstart1, kstart1)                            &
               + xdiff*(tempF_NN(istart1, jstart1, kstart1)                     &
                      - tempF_NN(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            F_NN(i,j,k) = u5 + zdiff*(u6 - u5)/dz
#ifdef PPDYN_TN
            ! ---- interp F_ee2 ----
            u1 = tempF_ee2(istart,  jstart,  kstart)                           &
               + xdiff*(tempF_ee2(istart1, jstart,  kstart)                    &
                      - tempF_ee2(istart,  jstart,  kstart))/dx
            u2 = tempF_ee2(istart,  jstart1, kstart)                           &
               + xdiff*(tempF_ee2(istart1, jstart1, kstart)                    &
                      - tempF_ee2(istart,  jstart1, kstart))/dx
            u3 = tempF_ee2(istart,  jstart,  kstart1)                          &
               + xdiff*(tempF_ee2(istart1, jstart,  kstart1)                   &
                      - tempF_ee2(istart,  jstart,  kstart1))/dx
            u4 = tempF_ee2(istart,  jstart1, kstart1)                          &
               + xdiff*(tempF_ee2(istart1, jstart1, kstart1)                   &
                      - tempF_ee2(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            F_ee2(i,j,k) = u5 + zdiff*(u6 - u5)/dz
            ! ---- interp F_deedt2 ----
            u1 = tempF_deedt2(istart,  jstart,  kstart)                        &
               + xdiff*(tempF_deedt2(istart1, jstart,  kstart)                 &
                      - tempF_deedt2(istart,  jstart,  kstart))/dx
            u2 = tempF_deedt2(istart,  jstart1, kstart)                        &
               + xdiff*(tempF_deedt2(istart1, jstart1, kstart)                 &
                      - tempF_deedt2(istart,  jstart1, kstart))/dx
            u3 = tempF_deedt2(istart,  jstart,  kstart1)                       &
               + xdiff*(tempF_deedt2(istart1, jstart,  kstart1)                &
                      - tempF_deedt2(istart,  jstart,  kstart1))/dx
            u4 = tempF_deedt2(istart,  jstart1, kstart1)                       &
               + xdiff*(tempF_deedt2(istart1, jstart1, kstart1)                &
                      - tempF_deedt2(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            F_deedt2(i,j,k) = u5 + zdiff*(u6 - u5)/dz
            ! ---- interp ee_past ----
            u1 = tempee_past(istart,  jstart,  kstart)                         &
               + xdiff*(tempee_past(istart1, jstart,  kstart)                  &
                      - tempee_past(istart,  jstart,  kstart))/dx
            u2 = tempee_past(istart,  jstart1, kstart)                         &
               + xdiff*(tempee_past(istart1, jstart1, kstart)                  &
                      - tempee_past(istart,  jstart1, kstart))/dx
            u3 = tempee_past(istart,  jstart,  kstart1)                        &
               + xdiff*(tempee_past(istart1, jstart,  kstart1)                 &
                      - tempee_past(istart,  jstart,  kstart1))/dx
            u4 = tempee_past(istart,  jstart1, kstart1)                        &
               + xdiff*(tempee_past(istart1, jstart1, kstart1)                 &
                      - tempee_past(istart,  jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2 - u1)/dy
            u6 = u3 + ydiff*(u4 - u3)/dy
            ee_past(i,j,k) = u5 + zdiff*(u6 - u5)/dz
#endif
        end do
        end do
    else
        ! Wall: u/v stay on uvp-node, z uses z_dev with 0.25*w(k+1) anomaly
        !$acc parallel loop collapse(2) default(present)                        &
        !$acc          private(px, py, pz, istart, jstart, kstart,              &
        !$acc                  istart1, jstart1, kstart1,                       &
        !$acc                  xdiff, ydiff, zdiff,                             &
        !$acc                  u1, u2, u3, u4, u5, u6) async(1)
        do j = 1, ny
        do i = 1, nx
            px = x_dev(i)  - u(i,j,k)*lagran_dt
            py = y_dev(j)  - v(i,j,k)*lagran_dt
            pz = z_dev(k)  - 0.25_rprec*w(i,j,k+1)*lagran_dt
            px = modulo(px, L_x)
            if (abs(px) / L_x < thresh) then
                istart = 1
            else if (abs(px - L_x) / L_x < thresh) then
                istart = nx
            else
                istart = floor(px / dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py) / L_y < thresh) then
                jstart = 1
            else if (abs(py - L_y) / L_y < thresh) then
                jstart = ny
            else
                jstart = floor(py / dy) + 1
            end if
            istart1 = autowrap_i_dev(istart + 1)
            jstart1 = autowrap_j_dev(jstart + 1)
            xdiff = px - x_dev(istart)
            ydiff = py - y_dev(jstart)
            ! lbc_mom > 0 + coord == 0: special wall branch
            if (pz < zw_dev(2)) then
                if (pz < z_dev(1)) then
                    kstart  = 1
                    kstart1 = 1
                    zdiff   = 0._rprec
                else
                    kstart  = 1
                    kstart1 = 2
                    zdiff   = 2._rprec*(pz - z_dev(kstart))
                end if
            else
                if (abs(pz - zw_dev(nz)) / L_z < thresh) then
                    kstart = nz-1
                else
                    kstart = floor((pz - zw_dev(1)) / dz) + 1
                end if
                kstart1 = kstart + 1
                zdiff   = pz - zw_dev(kstart)
            end if
            u1 = tempF_LM(istart, jstart, kstart) + xdiff*(tempF_LM(istart1, jstart, kstart) - tempF_LM(istart, jstart, kstart))/dx
            u2 = tempF_LM(istart, jstart1, kstart) + xdiff*(tempF_LM(istart1, jstart1, kstart) - tempF_LM(istart, jstart1, kstart))/dx
            u3 = tempF_LM(istart, jstart, kstart1) + xdiff*(tempF_LM(istart1, jstart, kstart1) - tempF_LM(istart, jstart, kstart1))/dx
            u4 = tempF_LM(istart, jstart1, kstart1) + xdiff*(tempF_LM(istart1, jstart1, kstart1) - tempF_LM(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_LM(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_MM(istart, jstart, kstart) + xdiff*(tempF_MM(istart1, jstart, kstart) - tempF_MM(istart, jstart, kstart))/dx
            u2 = tempF_MM(istart, jstart1, kstart) + xdiff*(tempF_MM(istart1, jstart1, kstart) - tempF_MM(istart, jstart1, kstart))/dx
            u3 = tempF_MM(istart, jstart, kstart1) + xdiff*(tempF_MM(istart1, jstart, kstart1) - tempF_MM(istart, jstart, kstart1))/dx
            u4 = tempF_MM(istart, jstart1, kstart1) + xdiff*(tempF_MM(istart1, jstart1, kstart1) - tempF_MM(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_MM(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_QN(istart, jstart, kstart) + xdiff*(tempF_QN(istart1, jstart, kstart) - tempF_QN(istart, jstart, kstart))/dx
            u2 = tempF_QN(istart, jstart1, kstart) + xdiff*(tempF_QN(istart1, jstart1, kstart) - tempF_QN(istart, jstart1, kstart))/dx
            u3 = tempF_QN(istart, jstart, kstart1) + xdiff*(tempF_QN(istart1, jstart, kstart1) - tempF_QN(istart, jstart, kstart1))/dx
            u4 = tempF_QN(istart, jstart1, kstart1) + xdiff*(tempF_QN(istart1, jstart1, kstart1) - tempF_QN(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_QN(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_NN(istart, jstart, kstart) + xdiff*(tempF_NN(istart1, jstart, kstart) - tempF_NN(istart, jstart, kstart))/dx
            u2 = tempF_NN(istart, jstart1, kstart) + xdiff*(tempF_NN(istart1, jstart1, kstart) - tempF_NN(istart, jstart1, kstart))/dx
            u3 = tempF_NN(istart, jstart, kstart1) + xdiff*(tempF_NN(istart1, jstart, kstart1) - tempF_NN(istart, jstart, kstart1))/dx
            u4 = tempF_NN(istart, jstart1, kstart1) + xdiff*(tempF_NN(istart1, jstart1, kstart1) - tempF_NN(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_NN(i,j,k) = u5 + zdiff*(u6-u5)/dz
#ifdef PPDYN_TN
            u1 = tempF_ee2(istart, jstart, kstart) + xdiff*(tempF_ee2(istart1, jstart, kstart) - tempF_ee2(istart, jstart, kstart))/dx
            u2 = tempF_ee2(istart, jstart1, kstart) + xdiff*(tempF_ee2(istart1, jstart1, kstart) - tempF_ee2(istart, jstart1, kstart))/dx
            u3 = tempF_ee2(istart, jstart, kstart1) + xdiff*(tempF_ee2(istart1, jstart, kstart1) - tempF_ee2(istart, jstart, kstart1))/dx
            u4 = tempF_ee2(istart, jstart1, kstart1) + xdiff*(tempF_ee2(istart1, jstart1, kstart1) - tempF_ee2(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_ee2(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_deedt2(istart, jstart, kstart) + xdiff*(tempF_deedt2(istart1, jstart, kstart) - tempF_deedt2(istart, jstart, kstart))/dx
            u2 = tempF_deedt2(istart, jstart1, kstart) + xdiff*(tempF_deedt2(istart1, jstart1, kstart) - tempF_deedt2(istart, jstart1, kstart))/dx
            u3 = tempF_deedt2(istart, jstart, kstart1) + xdiff*(tempF_deedt2(istart1, jstart, kstart1) - tempF_deedt2(istart, jstart, kstart1))/dx
            u4 = tempF_deedt2(istart, jstart1, kstart1) + xdiff*(tempF_deedt2(istart1, jstart1, kstart1) - tempF_deedt2(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_deedt2(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempee_past(istart, jstart, kstart) + xdiff*(tempee_past(istart1, jstart, kstart) - tempee_past(istart, jstart, kstart))/dx
            u2 = tempee_past(istart, jstart1, kstart) + xdiff*(tempee_past(istart1, jstart1, kstart) - tempee_past(istart, jstart1, kstart))/dx
            u3 = tempee_past(istart, jstart, kstart1) + xdiff*(tempee_past(istart1, jstart, kstart1) - tempee_past(istart, jstart, kstart1))/dx
            u4 = tempee_past(istart, jstart1, kstart1) + xdiff*(tempee_past(istart1, jstart1, kstart1) - tempee_past(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            ee_past(i,j,k) = u5 + zdiff*(u6-u5)/dz
#endif
        end do
        end do
    end if
    kmin = 2
else
    kmin = 1
end if

! ---------------------------------------------------------------------------
! Step 3: intermediate levels (kmin .. nz-1).
! ---------------------------------------------------------------------------
!$acc parallel loop collapse(3) default(present)                                &
!$acc          private(px, py, pz, istart, jstart, kstart,                      &
!$acc                  istart1, jstart1, kstart1,                               &
!$acc                  xdiff, ydiff, zdiff,                                     &
!$acc                  u1, u2, u3, u4, u5, u6) async(1)
do k = kmin, nz-1
do j = 1, ny
do i = 1, nx
    px = x_dev(i)  - 0.5_rprec*(u(i,j,k-1) + u(i,j,k))*lagran_dt
    py = y_dev(j)  - 0.5_rprec*(v(i,j,k-1) + v(i,j,k))*lagran_dt
    pz = zw_dev(k) - w(i,j,k)*lagran_dt
    px = modulo(px, L_x)
    if (abs(px) / L_x < thresh) then
        istart = 1
    else if (abs(px - L_x) / L_x < thresh) then
        istart = nx
    else
        istart = floor(px / dx) + 1
    end if
    py = modulo(py, L_y)
    if (abs(py) / L_y < thresh) then
        jstart = 1
    else if (abs(py - L_y) / L_y < thresh) then
        jstart = ny
    else
        jstart = floor(py / dy) + 1
    end if
    istart1 = autowrap_i_dev(istart + 1)
    jstart1 = autowrap_j_dev(jstart + 1)
    xdiff = px - x_dev(istart)
    ydiff = py - y_dev(jstart)
    ! generic z-branch with the two corner cases for stress-free walls
    if (coord == 0 .and. lbc_mom > 0 .and. pz < zw_dev(2)) then
        if (pz < z_dev(1)) then
            kstart = 1; kstart1 = 1; zdiff = 0._rprec
        else
            kstart = 1; kstart1 = 2; zdiff = 2._rprec*(pz - z_dev(kstart))
        end if
    else if (coord == nproc-1 .and. ubc_mom > 0 .and. pz > zw_dev(nz-1)) then
        if (pz > z_dev(nz-1)) then
            kstart = nz; kstart1 = nz; zdiff = 0._rprec
        else
            kstart = nz-1; kstart1 = nz; zdiff = 2._rprec*(pz - zw_dev(kstart))
        end if
    else
        if (abs(pz - zw_dev(nz)) / L_z < thresh) then
            kstart = nz-1
        else
            kstart = floor((pz - zw_dev(1)) / dz) + 1
        end if
        kstart1 = kstart + 1
        zdiff   = pz - zw_dev(kstart)
    end if
    u1 = tempF_LM(istart, jstart, kstart) + xdiff*(tempF_LM(istart1, jstart, kstart) - tempF_LM(istart, jstart, kstart))/dx
    u2 = tempF_LM(istart, jstart1, kstart) + xdiff*(tempF_LM(istart1, jstart1, kstart) - tempF_LM(istart, jstart1, kstart))/dx
    u3 = tempF_LM(istart, jstart, kstart1) + xdiff*(tempF_LM(istart1, jstart, kstart1) - tempF_LM(istart, jstart, kstart1))/dx
    u4 = tempF_LM(istart, jstart1, kstart1) + xdiff*(tempF_LM(istart1, jstart1, kstart1) - tempF_LM(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    F_LM(i,j,k) = u5 + zdiff*(u6-u5)/dz
    u1 = tempF_MM(istart, jstart, kstart) + xdiff*(tempF_MM(istart1, jstart, kstart) - tempF_MM(istart, jstart, kstart))/dx
    u2 = tempF_MM(istart, jstart1, kstart) + xdiff*(tempF_MM(istart1, jstart1, kstart) - tempF_MM(istart, jstart1, kstart))/dx
    u3 = tempF_MM(istart, jstart, kstart1) + xdiff*(tempF_MM(istart1, jstart, kstart1) - tempF_MM(istart, jstart, kstart1))/dx
    u4 = tempF_MM(istart, jstart1, kstart1) + xdiff*(tempF_MM(istart1, jstart1, kstart1) - tempF_MM(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    F_MM(i,j,k) = u5 + zdiff*(u6-u5)/dz
    u1 = tempF_QN(istart, jstart, kstart) + xdiff*(tempF_QN(istart1, jstart, kstart) - tempF_QN(istart, jstart, kstart))/dx
    u2 = tempF_QN(istart, jstart1, kstart) + xdiff*(tempF_QN(istart1, jstart1, kstart) - tempF_QN(istart, jstart1, kstart))/dx
    u3 = tempF_QN(istart, jstart, kstart1) + xdiff*(tempF_QN(istart1, jstart, kstart1) - tempF_QN(istart, jstart, kstart1))/dx
    u4 = tempF_QN(istart, jstart1, kstart1) + xdiff*(tempF_QN(istart1, jstart1, kstart1) - tempF_QN(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    F_QN(i,j,k) = u5 + zdiff*(u6-u5)/dz
    u1 = tempF_NN(istart, jstart, kstart) + xdiff*(tempF_NN(istart1, jstart, kstart) - tempF_NN(istart, jstart, kstart))/dx
    u2 = tempF_NN(istart, jstart1, kstart) + xdiff*(tempF_NN(istart1, jstart1, kstart) - tempF_NN(istart, jstart1, kstart))/dx
    u3 = tempF_NN(istart, jstart, kstart1) + xdiff*(tempF_NN(istart1, jstart, kstart1) - tempF_NN(istart, jstart, kstart1))/dx
    u4 = tempF_NN(istart, jstart1, kstart1) + xdiff*(tempF_NN(istart1, jstart1, kstart1) - tempF_NN(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    F_NN(i,j,k) = u5 + zdiff*(u6-u5)/dz
#ifdef PPDYN_TN
    u1 = tempF_ee2(istart, jstart, kstart) + xdiff*(tempF_ee2(istart1, jstart, kstart) - tempF_ee2(istart, jstart, kstart))/dx
    u2 = tempF_ee2(istart, jstart1, kstart) + xdiff*(tempF_ee2(istart1, jstart1, kstart) - tempF_ee2(istart, jstart1, kstart))/dx
    u3 = tempF_ee2(istart, jstart, kstart1) + xdiff*(tempF_ee2(istart1, jstart, kstart1) - tempF_ee2(istart, jstart, kstart1))/dx
    u4 = tempF_ee2(istart, jstart1, kstart1) + xdiff*(tempF_ee2(istart1, jstart1, kstart1) - tempF_ee2(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    F_ee2(i,j,k) = u5 + zdiff*(u6-u5)/dz
    u1 = tempF_deedt2(istart, jstart, kstart) + xdiff*(tempF_deedt2(istart1, jstart, kstart) - tempF_deedt2(istart, jstart, kstart))/dx
    u2 = tempF_deedt2(istart, jstart1, kstart) + xdiff*(tempF_deedt2(istart1, jstart1, kstart) - tempF_deedt2(istart, jstart1, kstart))/dx
    u3 = tempF_deedt2(istart, jstart, kstart1) + xdiff*(tempF_deedt2(istart1, jstart, kstart1) - tempF_deedt2(istart, jstart, kstart1))/dx
    u4 = tempF_deedt2(istart, jstart1, kstart1) + xdiff*(tempF_deedt2(istart1, jstart1, kstart1) - tempF_deedt2(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    F_deedt2(i,j,k) = u5 + zdiff*(u6-u5)/dz
    u1 = tempee_past(istart, jstart, kstart) + xdiff*(tempee_past(istart1, jstart, kstart) - tempee_past(istart, jstart, kstart))/dx
    u2 = tempee_past(istart, jstart1, kstart) + xdiff*(tempee_past(istart1, jstart1, kstart) - tempee_past(istart, jstart1, kstart))/dx
    u3 = tempee_past(istart, jstart, kstart1) + xdiff*(tempee_past(istart1, jstart, kstart1) - tempee_past(istart, jstart, kstart1))/dx
    u4 = tempee_past(istart, jstart1, kstart1) + xdiff*(tempee_past(istart1, jstart1, kstart1) - tempee_past(istart, jstart1, kstart1))/dx
    u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
    ee_past(i,j,k) = u5 + zdiff*(u6-u5)/dz
#endif
end do
end do
end do

! ---------------------------------------------------------------------------
! Step 4: top slab (coord==nproc-1, k=nz).
! ---------------------------------------------------------------------------
#ifdef PPMPI
if (coord == nproc-1) then
#endif
    k = nz
    if (ubc_mom == 0) then
        ! Stress-free at top: u/v use neighbor k-1, z is the wall zw(k)
        !$acc parallel loop collapse(2) default(present)                        &
        !$acc          private(px, py, pz, istart, jstart, kstart,              &
        !$acc                  istart1, jstart1, kstart1,                       &
        !$acc                  xdiff, ydiff, zdiff,                             &
        !$acc                  u1, u2, u3, u4, u5, u6) async(1)
        do j = 1, ny
        do i = 1, nx
            px = x_dev(i)  - u(i,j,k-1)*lagran_dt
            py = y_dev(j)  - v(i,j,k-1)*lagran_dt
            pz = zw_dev(k)
            px = modulo(px, L_x)
            if (abs(px) / L_x < thresh) then
                istart = 1
            else if (abs(px - L_x) / L_x < thresh) then
                istart = nx
            else
                istart = floor(px / dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py) / L_y < thresh) then
                jstart = 1
            else if (abs(py - L_y) / L_y < thresh) then
                jstart = ny
            else
                jstart = floor(py / dy) + 1
            end if
            istart1 = autowrap_i_dev(istart + 1)
            jstart1 = autowrap_j_dev(jstart + 1)
            xdiff = px - x_dev(istart)
            ydiff = py - y_dev(jstart)
            ! ubc_mom == 0 here, so the upper-wall branch can't fire.
            if (abs(pz - zw_dev(nz)) / L_z < thresh) then
                kstart = nz-1
            else
                kstart = floor((pz - zw_dev(1)) / dz) + 1
            end if
            kstart1 = kstart + 1
            zdiff   = pz - zw_dev(kstart)
            u1 = tempF_LM(istart, jstart, kstart) + xdiff*(tempF_LM(istart1, jstart, kstart) - tempF_LM(istart, jstart, kstart))/dx
            u2 = tempF_LM(istart, jstart1, kstart) + xdiff*(tempF_LM(istart1, jstart1, kstart) - tempF_LM(istart, jstart1, kstart))/dx
            u3 = tempF_LM(istart, jstart, kstart1) + xdiff*(tempF_LM(istart1, jstart, kstart1) - tempF_LM(istart, jstart, kstart1))/dx
            u4 = tempF_LM(istart, jstart1, kstart1) + xdiff*(tempF_LM(istart1, jstart1, kstart1) - tempF_LM(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_LM(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_MM(istart, jstart, kstart) + xdiff*(tempF_MM(istart1, jstart, kstart) - tempF_MM(istart, jstart, kstart))/dx
            u2 = tempF_MM(istart, jstart1, kstart) + xdiff*(tempF_MM(istart1, jstart1, kstart) - tempF_MM(istart, jstart1, kstart))/dx
            u3 = tempF_MM(istart, jstart, kstart1) + xdiff*(tempF_MM(istart1, jstart, kstart1) - tempF_MM(istart, jstart, kstart1))/dx
            u4 = tempF_MM(istart, jstart1, kstart1) + xdiff*(tempF_MM(istart1, jstart1, kstart1) - tempF_MM(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_MM(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_QN(istart, jstart, kstart) + xdiff*(tempF_QN(istart1, jstart, kstart) - tempF_QN(istart, jstart, kstart))/dx
            u2 = tempF_QN(istart, jstart1, kstart) + xdiff*(tempF_QN(istart1, jstart1, kstart) - tempF_QN(istart, jstart1, kstart))/dx
            u3 = tempF_QN(istart, jstart, kstart1) + xdiff*(tempF_QN(istart1, jstart, kstart1) - tempF_QN(istart, jstart, kstart1))/dx
            u4 = tempF_QN(istart, jstart1, kstart1) + xdiff*(tempF_QN(istart1, jstart1, kstart1) - tempF_QN(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_QN(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_NN(istart, jstart, kstart) + xdiff*(tempF_NN(istart1, jstart, kstart) - tempF_NN(istart, jstart, kstart))/dx
            u2 = tempF_NN(istart, jstart1, kstart) + xdiff*(tempF_NN(istart1, jstart1, kstart) - tempF_NN(istart, jstart1, kstart))/dx
            u3 = tempF_NN(istart, jstart, kstart1) + xdiff*(tempF_NN(istart1, jstart, kstart1) - tempF_NN(istart, jstart, kstart1))/dx
            u4 = tempF_NN(istart, jstart1, kstart1) + xdiff*(tempF_NN(istart1, jstart1, kstart1) - tempF_NN(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_NN(i,j,k) = u5 + zdiff*(u6-u5)/dz
#ifdef PPDYN_TN
            u1 = tempF_ee2(istart, jstart, kstart) + xdiff*(tempF_ee2(istart1, jstart, kstart) - tempF_ee2(istart, jstart, kstart))/dx
            u2 = tempF_ee2(istart, jstart1, kstart) + xdiff*(tempF_ee2(istart1, jstart1, kstart) - tempF_ee2(istart, jstart1, kstart))/dx
            u3 = tempF_ee2(istart, jstart, kstart1) + xdiff*(tempF_ee2(istart1, jstart, kstart1) - tempF_ee2(istart, jstart, kstart1))/dx
            u4 = tempF_ee2(istart, jstart1, kstart1) + xdiff*(tempF_ee2(istart1, jstart1, kstart1) - tempF_ee2(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_ee2(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_deedt2(istart, jstart, kstart) + xdiff*(tempF_deedt2(istart1, jstart, kstart) - tempF_deedt2(istart, jstart, kstart))/dx
            u2 = tempF_deedt2(istart, jstart1, kstart) + xdiff*(tempF_deedt2(istart1, jstart1, kstart) - tempF_deedt2(istart, jstart1, kstart))/dx
            u3 = tempF_deedt2(istart, jstart, kstart1) + xdiff*(tempF_deedt2(istart1, jstart, kstart1) - tempF_deedt2(istart, jstart, kstart1))/dx
            u4 = tempF_deedt2(istart, jstart1, kstart1) + xdiff*(tempF_deedt2(istart1, jstart1, kstart1) - tempF_deedt2(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_deedt2(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempee_past(istart, jstart, kstart) + xdiff*(tempee_past(istart1, jstart, kstart) - tempee_past(istart, jstart, kstart))/dx
            u2 = tempee_past(istart, jstart1, kstart) + xdiff*(tempee_past(istart1, jstart1, kstart) - tempee_past(istart, jstart1, kstart))/dx
            u3 = tempee_past(istart, jstart, kstart1) + xdiff*(tempee_past(istart1, jstart, kstart1) - tempee_past(istart, jstart, kstart1))/dx
            u4 = tempee_past(istart, jstart1, kstart1) + xdiff*(tempee_past(istart1, jstart1, kstart1) - tempee_past(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            ee_past(i,j,k) = u5 + zdiff*(u6-u5)/dz
#endif
        end do
        end do
    else
        ! Wall at top: similar to bottom wall (ubc_mom > 0 + special z-branch)
        !$acc parallel loop collapse(2) default(present)                        &
        !$acc          private(px, py, pz, istart, jstart, kstart,              &
        !$acc                  istart1, jstart1, kstart1,                       &
        !$acc                  xdiff, ydiff, zdiff,                             &
        !$acc                  u1, u2, u3, u4, u5, u6) async(1)
        do j = 1, ny
        do i = 1, nx
            px = x_dev(i)    - u(i,j,k-1)*lagran_dt
            py = y_dev(j)    - v(i,j,k-1)*lagran_dt
            pz = z_dev(k-1)  - 0.25_rprec*w(i,j,k-1)*lagran_dt
            px = modulo(px, L_x)
            if (abs(px) / L_x < thresh) then
                istart = 1
            else if (abs(px - L_x) / L_x < thresh) then
                istart = nx
            else
                istart = floor(px / dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py) / L_y < thresh) then
                jstart = 1
            else if (abs(py - L_y) / L_y < thresh) then
                jstart = ny
            else
                jstart = floor(py / dy) + 1
            end if
            istart1 = autowrap_i_dev(istart + 1)
            jstart1 = autowrap_j_dev(jstart + 1)
            xdiff = px - x_dev(istart)
            ydiff = py - y_dev(jstart)
            ! coord == nproc-1 + ubc_mom > 0: special top branch
            if (pz > zw_dev(nz-1)) then
                if (pz > z_dev(nz-1)) then
                    kstart  = nz; kstart1 = nz; zdiff = 0._rprec
                else
                    kstart  = nz-1; kstart1 = nz
                    zdiff   = 2._rprec*(pz - zw_dev(kstart))
                end if
            else
                if (abs(pz - zw_dev(nz)) / L_z < thresh) then
                    kstart = nz-1
                else
                    kstart = floor((pz - zw_dev(1)) / dz) + 1
                end if
                kstart1 = kstart + 1
                zdiff   = pz - zw_dev(kstart)
            end if
            u1 = tempF_LM(istart, jstart, kstart) + xdiff*(tempF_LM(istart1, jstart, kstart) - tempF_LM(istart, jstart, kstart))/dx
            u2 = tempF_LM(istart, jstart1, kstart) + xdiff*(tempF_LM(istart1, jstart1, kstart) - tempF_LM(istart, jstart1, kstart))/dx
            u3 = tempF_LM(istart, jstart, kstart1) + xdiff*(tempF_LM(istart1, jstart, kstart1) - tempF_LM(istart, jstart, kstart1))/dx
            u4 = tempF_LM(istart, jstart1, kstart1) + xdiff*(tempF_LM(istart1, jstart1, kstart1) - tempF_LM(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_LM(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_MM(istart, jstart, kstart) + xdiff*(tempF_MM(istart1, jstart, kstart) - tempF_MM(istart, jstart, kstart))/dx
            u2 = tempF_MM(istart, jstart1, kstart) + xdiff*(tempF_MM(istart1, jstart1, kstart) - tempF_MM(istart, jstart1, kstart))/dx
            u3 = tempF_MM(istart, jstart, kstart1) + xdiff*(tempF_MM(istart1, jstart, kstart1) - tempF_MM(istart, jstart, kstart1))/dx
            u4 = tempF_MM(istart, jstart1, kstart1) + xdiff*(tempF_MM(istart1, jstart1, kstart1) - tempF_MM(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_MM(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_QN(istart, jstart, kstart) + xdiff*(tempF_QN(istart1, jstart, kstart) - tempF_QN(istart, jstart, kstart))/dx
            u2 = tempF_QN(istart, jstart1, kstart) + xdiff*(tempF_QN(istart1, jstart1, kstart) - tempF_QN(istart, jstart1, kstart))/dx
            u3 = tempF_QN(istart, jstart, kstart1) + xdiff*(tempF_QN(istart1, jstart, kstart1) - tempF_QN(istart, jstart, kstart1))/dx
            u4 = tempF_QN(istart, jstart1, kstart1) + xdiff*(tempF_QN(istart1, jstart1, kstart1) - tempF_QN(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_QN(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_NN(istart, jstart, kstart) + xdiff*(tempF_NN(istart1, jstart, kstart) - tempF_NN(istart, jstart, kstart))/dx
            u2 = tempF_NN(istart, jstart1, kstart) + xdiff*(tempF_NN(istart1, jstart1, kstart) - tempF_NN(istart, jstart1, kstart))/dx
            u3 = tempF_NN(istart, jstart, kstart1) + xdiff*(tempF_NN(istart1, jstart, kstart1) - tempF_NN(istart, jstart, kstart1))/dx
            u4 = tempF_NN(istart, jstart1, kstart1) + xdiff*(tempF_NN(istart1, jstart1, kstart1) - tempF_NN(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_NN(i,j,k) = u5 + zdiff*(u6-u5)/dz
#ifdef PPDYN_TN
            u1 = tempF_ee2(istart, jstart, kstart) + xdiff*(tempF_ee2(istart1, jstart, kstart) - tempF_ee2(istart, jstart, kstart))/dx
            u2 = tempF_ee2(istart, jstart1, kstart) + xdiff*(tempF_ee2(istart1, jstart1, kstart) - tempF_ee2(istart, jstart1, kstart))/dx
            u3 = tempF_ee2(istart, jstart, kstart1) + xdiff*(tempF_ee2(istart1, jstart, kstart1) - tempF_ee2(istart, jstart, kstart1))/dx
            u4 = tempF_ee2(istart, jstart1, kstart1) + xdiff*(tempF_ee2(istart1, jstart1, kstart1) - tempF_ee2(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_ee2(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempF_deedt2(istart, jstart, kstart) + xdiff*(tempF_deedt2(istart1, jstart, kstart) - tempF_deedt2(istart, jstart, kstart))/dx
            u2 = tempF_deedt2(istart, jstart1, kstart) + xdiff*(tempF_deedt2(istart1, jstart1, kstart) - tempF_deedt2(istart, jstart1, kstart))/dx
            u3 = tempF_deedt2(istart, jstart, kstart1) + xdiff*(tempF_deedt2(istart1, jstart, kstart1) - tempF_deedt2(istart, jstart, kstart1))/dx
            u4 = tempF_deedt2(istart, jstart1, kstart1) + xdiff*(tempF_deedt2(istart1, jstart1, kstart1) - tempF_deedt2(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            F_deedt2(i,j,k) = u5 + zdiff*(u6-u5)/dz
            u1 = tempee_past(istart, jstart, kstart) + xdiff*(tempee_past(istart1, jstart, kstart) - tempee_past(istart, jstart, kstart))/dx
            u2 = tempee_past(istart, jstart1, kstart) + xdiff*(tempee_past(istart1, jstart1, kstart) - tempee_past(istart, jstart1, kstart))/dx
            u3 = tempee_past(istart, jstart, kstart1) + xdiff*(tempee_past(istart1, jstart, kstart1) - tempee_past(istart, jstart, kstart1))/dx
            u4 = tempee_past(istart, jstart1, kstart1) + xdiff*(tempee_past(istart1, jstart1, kstart1) - tempee_past(istart, jstart1, kstart1))/dx
            u5 = u1 + ydiff*(u2-u1)/dy; u6 = u3 + ydiff*(u4-u3)/dy
            ee_past(i,j,k) = u5 + zdiff*(u6-u5)/dz
#endif
        end do
        end do
    end if
#ifdef PPMPI
end if
#endif

! ---------------------------------------------------------------------------
! Step 5: MPI ghost-slab sync of new F_*.
! ---------------------------------------------------------------------------
#ifdef PPMPI
call sync_downup_F(F_LM)
call sync_downup_F(F_MM)
call sync_downup_F(F_QN)
call sync_downup_F(F_NN)
#ifdef PPDYN_TN
call sync_downup_F(F_ee2)
call sync_downup_F(F_deedt2)
call sync_downup_F(ee_past)
#endif
#endif

! Diagnostic CFL print (host-side, infrequent, doesn't need GPU data).
if (mod(jt_total, lag_cfl_count) == 0) then
    lcfl = get_max_cfl()
    lcfl = lcfl*lagran_dt/dt
    if (coord == 0) print*, 'Lagrangian CFL condition= ', lcfl
end if

end subroutine interpolag_Sdep_gpu

#ifdef PPMPI
!*******************************************************************************
subroutine sync_downup_F(F)
!*******************************************************************************
! GPU equivalent of mpi_sync_real_array(F, 0, MPI_SYNC_DOWNUP). Pull only
! the two slabs the sendrecv reads (k=1 for sync_down, k=nz-1 for sync_up),
! do the host-side MPI exchange, push only the two slabs the sendrecv wrote
! (k=nz for sync_down, k=0 for sync_up).
!*******************************************************************************
use mpi
use param, only : ld, ny, nz, MPI_RPREC, down, up, comm, status, ierr
implicit none
real(rprec), dimension(ld,ny,0:nz), intent(inout) :: F

#ifdef PPGPU_AWARE_MPI
!$acc wait(1)
!$acc host_data use_device(F)
call mpi_sendrecv(F(:,:,1),  ld*ny, MPI_RPREC, down, 1,                         &
                  F(:,:,nz), ld*ny, MPI_RPREC, up,   1, comm, status, ierr)
call mpi_sendrecv(F(:,:,nz-1), ld*ny, MPI_RPREC, up,   2,                       &
                  F(:,:,0),    ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
!$acc end host_data
#else
!$acc wait(1)
!$acc update self(F(:,:,1))
call mpi_sendrecv(F(:,:,1),  ld*ny, MPI_RPREC, down, 1,                         &
                  F(:,:,nz), ld*ny, MPI_RPREC, up,   1, comm, status, ierr)
!$acc wait(1)
!$acc update device(F(:,:,nz))

!$acc wait(1)
!$acc update self(F(:,:,nz-1))
call mpi_sendrecv(F(:,:,nz-1), ld*ny, MPI_RPREC, up,   2,                       &
                  F(:,:,0),    ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
!$acc wait(1)
!$acc update device(F(:,:,0))
#endif

end subroutine sync_downup_F
#endif

#endif /* PPSGS_GPU */
end module lagrange_Sdep_gpu_m
