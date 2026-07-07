!!
!!  Copyright (C) 2019  Johns Hopkins University
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
module scalars
!*******************************************************************************
! This module contains all of the subroutines associated with scalar transport
!
! Navigation map:
!   - module state: scalar fields, SGS scalar diffusivity, Monin-Obukhov arrays
!   - GPU/ACC helpers: scalars_*_acc and scalars_*_gpu routines near the top
!   - initialization: scalars_init and ic_scal* routines
!   - derivatives and stability: scalars_deriv and obukhov
!   - timestep update: scalars_transport
!   - output/restart: scalars_checkpoint
!
! `USE_SCALARS_GPU` is valid only with `USE_SCALARS` and `USE_LES_GPU`.
! Keep scalar CPU and GPU paths numerically paired when adding transport terms.
!
! Ownership map:
!   - this module owns theta/RHS_T, scalar gradients, scalar SGS diffusivity,
!     and Monin-Obukhov scalar arrays.
!   - velocity fields are borrowed from the LES core during scalar transport;
!     scalar code should not make them host-authoritative.
!   - fft_gpu owns GPU FFT plans used by the scalar GPU path.
!   - scalar checkpoint output is the explicit host-visible boundary for
!     theta/RHS_T/psi_m.
use types, only : rprec
use param, only : path
#ifdef PPSCALARS_GPU
use openacc
use fft_gpu, only : fft_gpu_exec_d2z, fft_gpu_exec_z2d,                       &
    plan_forw_small_full, plan_back_small_full, plan_back_big_full,            &
    plan_forw_big_nzm1, plan_back_small_nzm1
use fft, only : kx, ky
use derivatives_gpu_m, only : filt_da_gpu, ddz_uv_gpu, ddx_gpu, ddy_gpu,      &
    ddz_w_gpu
#endif
implicit none

save
private

public :: scalars_init, ic_scal, buoyancy_force, scalars_transport,            &
    scalars_checkpoint, obukhov, scalars_deriv

real(rprec), public, dimension(:,:,:), allocatable :: theta, dTdx, dTdy, dTdz, &
    RHS_T, RHS_Tf, u_big, v_big, w_big, dTdx_big, dTdy_big, dTdz_big, RHS_big, &
    pi_x, pi_y, pi_z, div_pi, temp_var

! SGS values
real(rprec), public :: Pr_sgs = 0.5
real(rprec), public, dimension(:,:,:), allocatable :: Kappa_t

! Monin-Obukhov BC
real(rprec), public, dimension(:,:), allocatable :: psi_m, phi_m, psi_h, phi_h,&
    L, tstar_lbc

! Gravitational acceleration
real(rprec), public :: g = 9.81_rprec
! Roughness length for scalars. typically zo/10
real(rprec), public :: zo_s = 0.00001_rprec
! Treat theta as passive_scalar (no buoyancy)
logical, public :: passive_scalar = .false.
! Whether to initialize theta field
logical, public :: inits = .true.
! Name of file for restarting
character(64) :: fname
! Reference temperature scale
real(rprec), public :: T_scale = 300._rprec

! Boundary conditions
! lbc: lower boundary condition
! ubc: upper boundary condition
!       0 - prescribed temperature, 1 - prescribed flux
integer, public :: lbc_scal = 0
real(rprec), public :: scal_bot = 300._rprec
real(rprec), public :: flux_bot = 0._rprec
real(rprec), dimension(:), allocatable, public :: ic_z, ic_theta
real(rprec), public :: ic_no_vel_noise_z
integer, public :: ic_nloc
real(rprec), public :: lapse_rate = 0._rprec
logical, public :: read_lbc_scal = .false.

! Interpolation of bottom boundary condition
real(rprec), dimension(:), allocatable :: t_interp, lbc_interp

#if defined(PPSCALARS_GPU)
integer, save :: scalars_stage_count = 0
real(rprec), save :: scalars_time_copy_rhs = 0._rprec
real(rprec), save :: scalars_time_to_big = 0._rprec
real(rprec), save :: scalars_time_advective = 0._rprec
real(rprec), save :: scalars_time_return = 0._rprec
real(rprec), save :: scalars_time_flux = 0._rprec
real(rprec), save :: scalars_time_div_update = 0._rprec
real(rprec), save :: scalars_time_divergence = 0._rprec
real(rprec), save :: scalars_time_rhs_update = 0._rprec
real(rprec), save :: scalars_time_halo = 0._rprec
real(rprec), save :: scalars_time_total = 0._rprec
logical, save :: scalar_deriv_big_ready = .false.
#endif

#ifdef PPSCALARS_GPU
real(rprec), dimension(:), allocatable :: theta_bar_acc
! Scalar GPU path for the active LES GPU build.  This follows the explicit
! OpenACC residency model used by sim_param/sgs_param instead of the older
! Legacy managed-memory branch removed.
!$acc declare create(theta, dTdx, dTdy, dTdz, RHS_T, RHS_Tf)
!$acc declare create(u_big, v_big, w_big, dTdx_big, dTdy_big, dTdz_big, RHS_big)
!$acc declare create(pi_x, pi_y, pi_z, div_pi, temp_var, Kappa_t)
!$acc declare create(psi_m, phi_m, psi_h, phi_h, L, tstar_lbc)
#endif

contains

!*******************************************************************************
logical function scalars_env_true_token_enabled(name)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
character(len=32) :: setting
integer :: stat

scalars_env_true_token_enabled = .false.
call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('1', 'true', 'TRUE', 'True', 'on', 'ON', 'On', 'yes', 'YES', 'Yes')
        scalars_env_true_token_enabled = .true.
    case default
        scalars_env_true_token_enabled = .false.
    end select
end if

end function scalars_env_true_token_enabled

#ifdef PPSCALARS_GPU
!*******************************************************************************
logical function scalars_acc_enabled()
!*******************************************************************************
implicit none

scalars_acc_enabled = .true.

end function scalars_acc_enabled

!*******************************************************************************
subroutine scalars_acc_sync()
!*******************************************************************************
implicit none

!$acc wait

end subroutine scalars_acc_sync

!*******************************************************************************
logical function scalars_stage_timing_enabled()
!*******************************************************************************
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = scalars_env_true_token_enabled('LESGO_SCALAR_STAGE_TIMING')
    initialized = .true.
end if

scalars_stage_timing_enabled = enabled

end function scalars_stage_timing_enabled

!*******************************************************************************
subroutine scalars_timer_start(t0)
!*******************************************************************************
implicit none

real(rprec), intent(out) :: t0

call scalars_acc_sync()
call cpu_time(t0)

end subroutine scalars_timer_start

!*******************************************************************************
subroutine scalars_timer_accum(t0, accum, where)
!*******************************************************************************
implicit none

real(rprec), intent(inout) :: t0
real(rprec), intent(inout) :: accum
character(len=*), intent(in) :: where
real(rprec) :: t1

call scalars_acc_sync()
call cpu_time(t1)
accum = accum + max(t1 - t0, 0._rprec)
t0 = t1

end subroutine scalars_timer_accum

!*******************************************************************************
subroutine scalars_stage_report(wbase_in)
!*******************************************************************************
use param, only : coord
implicit none

integer, intent(in) :: wbase_in

if (coord /= 0) return
if (wbase_in > 0) then
    if (mod(scalars_stage_count, wbase_in) /= 0) return
end if

write(*,'(a,i8)') 'Scalar stage timing (rank 0 cumulative), call ',           &
    scalars_stage_count
write(*,'(1a,E15.7)') '  copy RHS: ', scalars_time_copy_rhs
write(*,'(1a,E15.7)') '  to_big transforms: ', scalars_time_to_big
write(*,'(1a,E15.7)') '  advective product: ', scalars_time_advective
write(*,'(1a,E15.7)') '  return RHS transform: ', scalars_time_return
write(*,'(1a,E15.7)') '  flux build: ', scalars_time_flux
write(*,'(1a,E15.7)') '  div/update/halo: ', scalars_time_div_update
write(*,'(1a,E15.7)') '    divergence: ', scalars_time_divergence
write(*,'(1a,E15.7)') '    rhs/theta update: ', scalars_time_rhs_update
write(*,'(1a,E15.7)') '    halo/top bc: ', scalars_time_halo
write(*,'(1a,E15.7)') '  scalar transport total: ', scalars_time_total

end subroutine scalars_stage_report

!*******************************************************************************
subroutine scalars_deriv_xy_big_acc()
!*******************************************************************************
use param, only : lbz, ld, ld_big, nx, ny, ny2, nz, dz, coord, nproc
implicit none

integer :: i, j, k
integer :: ir, ii, lh, ny_h, j_s, j_big_s
real(rprec) :: const, const_z, ar, ai, kxv, kyv

const = 1._rprec / real(nx*ny, rprec)
const_z = 1._rprec / dz
lh = nx/2 + 1
ny_h = ny/2
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

!$acc data present(theta, dTdx, dTdy, dTdx_big, dTdy_big, dTdz_big, RHS_big, kx, ky)

!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, ld
    theta(i,j,k) = const*theta(i,j,k)
end do
end do
end do

call fft_gpu_exec_d2z(plan_forw_small_full, theta(1,1,lbz))

!$acc parallel loop collapse(2) default(present) async(1)
do k = lbz, nz
do j = 1, ny
    theta(ld-1,j,k) = 0._rprec
    theta(ld,j,k) = 0._rprec
end do
end do

!$acc parallel loop collapse(2) default(present) async(1)
do k = lbz, nz
do i = 1, ld
    theta(i,ny_h+1,k) = 0._rprec
end do
end do

!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, lh
    ir = 2*i - 1
    ii = 2*i
    ar = theta(ir,j,k)
    ai = theta(ii,j,k)
    kxv = kx(i,j)
    kyv = ky(i,j)
    dTdx(ir,j,k) = -ai*kxv
    dTdx(ii,j,k) = ar*kxv
    dTdy(ir,j,k) = -ai*kyv
    dTdy(ii,j,k) = ar*kyv
end do
end do
end do

! Pad spectral derivatives and theta directly to the 3/2 grid before the big
! inverse transforms.  RHS_big carries the padded theta field here; after the
! inverse, a vertical stencil gives dTdz_big directly and removes the later
! small-grid dTdz forward FFT.
!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny2
do i = 1, ld_big
    if (i <= nx .and. j <= ny_h) then
        dTdx_big(i,j,k) = dTdx(i,j,k)
        dTdy_big(i,j,k) = dTdy(i,j,k)
        RHS_big(i,j,k) = theta(i,j,k)
    else if (i <= nx .and. j >= j_big_s) then
        dTdx_big(i,j,k) = dTdx(i,j - j_big_s + j_s,k)
        dTdy_big(i,j,k) = dTdy(i,j - j_big_s + j_s,k)
        RHS_big(i,j,k) = theta(i,j - j_big_s + j_s,k)
    else
        dTdx_big(i,j,k) = 0._rprec
        dTdy_big(i,j,k) = 0._rprec
        RHS_big(i,j,k) = 0._rprec
    end if
end do
end do
end do

call fft_gpu_exec_z2d(plan_back_big_full, dTdx_big(1,1,lbz))
call fft_gpu_exec_z2d(plan_back_big_full, dTdy_big(1,1,lbz))
call fft_gpu_exec_z2d(plan_back_big_full, RHS_big(1,1,lbz))

! Horizontal interpolation/dealiasing and this vertical stencil commute.  Using
! theta_big in RHS_big avoids a separate to_big(dTdz) small-grid forward FFT.
!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz+1, nz
do j = 1, ny2
do i = 1, ld_big
    dTdz_big(i,j,k) = const_z*(RHS_big(i,j,k) - RHS_big(i,j,k-1))
end do
end do
end do

!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny2
do i = 1, ld_big
    dTdz_big(i,j,lbz) = 0._rprec
end do
end do

if (coord == nproc-1) then
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 1, ny2
    do i = 1, ld_big
        dTdz_big(i,j,nz) = lapse_rate
    end do
    end do
end if

call fft_gpu_exec_z2d(plan_back_small_full, theta(1,1,lbz))
call fft_gpu_exec_z2d(plan_back_small_full, dTdx(1,1,lbz))
call fft_gpu_exec_z2d(plan_back_small_full, dTdy(1,1,lbz))

! The caller waits before host-visible scalar halo/checkpoint work.  Keeping
! this routine queued on async(1) avoids an extra per-step stream sync.
!$acc end data

scalar_deriv_big_ready = .true.

end subroutine scalars_deriv_xy_big_acc

!*******************************************************************************
subroutine scalars_to_big_acc(a, a_big)
!*******************************************************************************
use param, only : lbz, ld, ld_big, nx, ny, ny2, nz
implicit none

real(rprec), intent(in) :: a(ld, ny, lbz:nz)
real(rprec), intent(inout) :: a_big(ld_big, ny2, lbz:nz)
integer :: i, j, k
integer :: ny_h, j_s, j_big_s
real(rprec) :: const

const = 1._rprec/real(nx*ny, rprec)
ny_h = ny / 2
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, ld
    temp_var(i,j,k) = const*a(i,j,k)
end do
end do
end do

call fft_gpu_exec_d2z(plan_forw_small_full, temp_var(1,1,lbz))

    ! Fuse zero-fill and spectral padding into one launch.  This preserves the
    ! CPU padd() layout while avoiding three kernels per scalar field.
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = lbz, nz
    do j = 1, ny2
    do i = 1, ld_big
        if (i <= nx .and. j <= ny_h) then
            a_big(i,j,k) = temp_var(i,j,k)
        else if (i <= nx .and. j >= j_big_s) then
            a_big(i,j,k) = temp_var(i,j - j_big_s + j_s,k)
        else
            a_big(i,j,k) = 0._rprec
        end if
    end do
    end do
    end do

call fft_gpu_exec_z2d(plan_back_big_full, a_big(1,1,lbz))

end subroutine scalars_to_big_acc

!*******************************************************************************
subroutine scalars_return_rhs_acc()
!*******************************************************************************
use param, only : ld, ld_big, nx, ny, ny2, nz
implicit none

integer :: i, j, k
integer :: ny_h, j_s, j_big_s
integer :: jb

ny_h = ny / 2
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

call fft_gpu_exec_d2z(plan_forw_big_nzm1, RHS_big(1,1,1))

! Fuse the old lower-half copy, oddball zero-fill, Nyquist-row zero-fill, and
! upper-half copy into one launch.  This preserves the CPU unpadd() layout:
! j=1:ny/2 copied directly, j=ny/2+1 zeroed, and the upper half copied from
! the high-y end of the 3/2 grid.
!$acc parallel loop collapse(3) default(present) async(1) private(jb)
do k = 1, nz-1
do j = 1, ny
do i = 1, ld
    if (i > nx .or. j == ny_h + 1) then
        RHS_T(i,j,k) = 0._rprec
    else if (j <= ny_h) then
        RHS_T(i,j,k) = RHS_big(i,j,k)
    else
        jb = j_big_s + (j - j_s)
        RHS_T(i,j,k) = RHS_big(i,jb,k)
    end if
end do
end do
end do

call fft_gpu_exec_z2d(plan_back_small_nzm1, RHS_T(1,1,1))

end subroutine scalars_return_rhs_acc

!*******************************************************************************
subroutine scalars_divergence_acc()
!*******************************************************************************
use param, only : lbz, ld, nx, ny, nz
implicit none

integer :: i, j, k, ir, ii, lh, ny_h
real(rprec) :: const, arx, aix, ary, aiy, kxv, kyv

const = 1._rprec / real(nx*ny, rprec)
lh = nx/2 + 1
ny_h = ny/2

! Compute d(pi_x)/dx + d(pi_y)/dy directly in spectral space.  This keeps the
! two required forward FFTs but replaces two inverse FFTs plus a physical-space
! xy sum with one fused spectral sum and one inverse FFT.
!$acc data present(pi_x, pi_y, div_pi, temp_var, kx, ky)

!$acc parallel loop collapse(3) default(present) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, ld
    div_pi(i,j,k) = const*pi_x(i,j,k)
    temp_var(i,j,k) = const*pi_y(i,j,k)
end do
end do
end do

call fft_gpu_exec_d2z(plan_forw_small_full, div_pi(1,1,lbz))
call fft_gpu_exec_d2z(plan_forw_small_full, temp_var(1,1,lbz))

!$acc parallel loop collapse(3) default(present) private(arx, aix, ary, aiy, kxv, kyv, ir, ii) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, lh
    ir = 2*i - 1
    ii = 2*i
    if (i == lh .or. j == ny_h+1) then
        div_pi(ir,j,k) = 0._rprec
        div_pi(ii,j,k) = 0._rprec
    else
        arx = div_pi(ir,j,k)
        aix = div_pi(ii,j,k)
        ary = temp_var(ir,j,k)
        aiy = temp_var(ii,j,k)
        kxv = kx(i,j)
        kyv = ky(i,j)
        div_pi(ir,j,k) = -aix*kxv - aiy*kyv
        div_pi(ii,j,k) =  arx*kxv + ary*kyv
    end if
end do
end do
end do

call fft_gpu_exec_z2d(plan_back_small_full, div_pi(1,1,lbz))

!$acc end data

end subroutine scalars_divergence_acc

!*******************************************************************************
subroutine scalars_update_device_state()
!*******************************************************************************
implicit none

!$acc update device(theta, dTdx, dTdy, dTdz, RHS_T, RHS_Tf)
!$acc update device(u_big, v_big, w_big, dTdx_big, dTdy_big, dTdz_big, RHS_big)
!$acc update device(pi_x, pi_y, pi_z, div_pi, temp_var, Kappa_t)
!$acc update device(psi_m, phi_m, psi_h, phi_h, L, tstar_lbc)

end subroutine scalars_update_device_state

!*******************************************************************************
subroutine scalars_copy_rhs_acc()
!*******************************************************************************
use param, only : lbz, ld, ny, nz
implicit none

integer :: i, j, k

!$acc parallel loop collapse(3) present(RHS_Tf, RHS_T) async(1)
do k = lbz, nz
do j = 1, ny
do i = 1, ld
    RHS_Tf(i,j,k) = RHS_T(i,j,k)
end do
end do
end do

end subroutine scalars_copy_rhs_acc

!*******************************************************************************
subroutine scalars_advective_acc(jz_min, jz_max, const,                       &
    u_adv_big, v_adv_big, w_adv_big)
!*******************************************************************************
use param, only : lbz, ld_big, ny2, nz, nproc, coord
implicit none

integer, intent(in) :: jz_min, jz_max
real(rprec), intent(in) :: const
real(rprec), intent(in) :: u_adv_big(ld_big, ny2, lbz:nz)
real(rprec), intent(in) :: v_adv_big(ld_big, ny2, lbz:nz)
real(rprec), intent(in) :: w_adv_big(ld_big, ny2, lbz:nz)
integer :: i, j, k

!$acc parallel loop collapse(3) present(RHS_big, u_adv_big, v_adv_big, w_adv_big, dTdx_big, dTdy_big, dTdz_big) async(1)
do k = jz_min, jz_max
do j = 1, ny2
do i = 1, ld_big
    RHS_big(i,j,k) = const*(u_adv_big(i,j,k)*dTdx_big(i,j,k)                   &
        + v_adv_big(i,j,k)*dTdy_big(i,j,k)                                     &
        + 0.5_rprec*w_adv_big(i,j,k+1)*dTdz_big(i,j,k+1)                       &
        + 0.5_rprec*w_adv_big(i,j,k)*dTdz_big(i,j,k))
end do
end do
end do

if (coord == 0) then
    !$acc parallel loop collapse(2) present(RHS_big, u_adv_big, v_adv_big, w_adv_big, dTdx_big, dTdy_big, dTdz_big) async(1)
    do j = 1, ny2
    do i = 1, ld_big
        RHS_big(i,j,1) = const*(u_adv_big(i,j,1)*dTdx_big(i,j,1)               &
            + v_adv_big(i,j,1)*dTdy_big(i,j,1)                                 &
            + 0.5_rprec*w_adv_big(i,j,2)*dTdz_big(i,j,2))
    end do
    end do
end if

if (coord == nproc-1) then
    !$acc parallel loop collapse(2) present(RHS_big, u_adv_big, v_adv_big, w_adv_big, dTdx_big, dTdy_big, dTdz_big) async(1)
    do j = 1, ny2
    do i = 1, ld_big
        RHS_big(i,j,nz-1) = const*(u_adv_big(i,j,nz-1)*dTdx_big(i,j,nz-1)      &
            + v_adv_big(i,j,nz-1)*dTdy_big(i,j,nz-1)                           &
            + 0.5_rprec*w_adv_big(i,j,nz-1)*dTdz_big(i,j,nz-1))
    end do
    end do
end if

end subroutine scalars_advective_acc

!*******************************************************************************
subroutine scalars_flux_acc(jz_min, jz_max)
!*******************************************************************************
use param, only : lbz, ld, ny, nz, nproc, coord, lbc_mom, ubc_mom
use sgs_param, only : Nu_t
implicit none

    integer, intent(in) :: jz_min, jz_max
    integer :: i, j, k
    real(rprec) :: inv_Pr_sgs

    inv_Pr_sgs = 1._rprec / Pr_sgs

    if (coord == 0) then
        select case (lbc_mom)
        case (0)
            !$acc parallel loop collapse(2) present(pi_x, pi_y, Nu_t, dTdx, dTdy) async(1)
            do j = 1, ny
            do i = 1, ld
                pi_x(i,j,1) = -0.5_rprec*inv_Pr_sgs                               &
                    *(Nu_t(i,j,1) + Nu_t(i,j,2))*dTdx(i,j,1)
                pi_y(i,j,1) = -0.5_rprec*inv_Pr_sgs                               &
                    *(Nu_t(i,j,1) + Nu_t(i,j,2))*dTdy(i,j,1)
            end do
            end do
        case (1:)
            !$acc parallel loop collapse(2) present(pi_x, pi_y, Nu_t, dTdx, dTdy) async(1)
            do j = 1, ny
            do i = 1, ld
                pi_x(i,j,1) = -inv_Pr_sgs*Nu_t(i,j,1)*dTdx(i,j,1)
                pi_y(i,j,1) = -inv_Pr_sgs*Nu_t(i,j,1)*dTdy(i,j,1)
            end do
            end do
        end select
    end if

    if (coord == nproc-1) then
        select case (ubc_mom)
        case (0)
            !$acc parallel loop collapse(2) present(pi_x, pi_y, pi_z, Nu_t, dTdx, dTdy, dTdz) async(1)
            do j = 1, ny
            do i = 1, ld
                pi_x(i,j,nz-1) = -0.5_rprec*inv_Pr_sgs                           &
                    *(Nu_t(i,j,nz-1) + Nu_t(i,j,nz))*dTdx(i,j,nz-1)
                pi_y(i,j,nz-1) = -0.5_rprec*inv_Pr_sgs                           &
                    *(Nu_t(i,j,nz-1) + Nu_t(i,j,nz))*dTdy(i,j,nz-1)
                pi_z(i,j,nz) = -inv_Pr_sgs*Nu_t(i,j,nz-1)*dTdz(i,j,nz)
            end do
            end do
        case (1:)
            !$acc parallel loop collapse(2) present(pi_x, pi_y, pi_z, Nu_t, dTdx, dTdy, dTdz) async(1)
            do j = 1, ny
            do i = 1, ld
                pi_x(i,j,nz-1) = -inv_Pr_sgs*Nu_t(i,j,nz-1)*dTdx(i,j,nz-1)
                pi_y(i,j,nz-1) = -inv_Pr_sgs*Nu_t(i,j,nz-1)*dTdy(i,j,nz-1)
                pi_z(i,j,nz) = -inv_Pr_sgs*Nu_t(i,j,nz-1)*dTdz(i,j,nz)
            end do
            end do
        end select
    end if

    !$acc parallel loop collapse(3) present(pi_x, pi_y, pi_z, Nu_t, dTdx, dTdy, dTdz) async(1)
    do k = jz_min, jz_max
    do j = 1, ny
    do i = 1, ld
        pi_x(i,j,k) = -0.5_rprec*inv_Pr_sgs                                     &
            *(Nu_t(i,j,k) + Nu_t(i,j,k+1))*dTdx(i,j,k)
        pi_y(i,j,k) = -0.5_rprec*inv_Pr_sgs                                     &
            *(Nu_t(i,j,k) + Nu_t(i,j,k+1))*dTdy(i,j,k)
        pi_z(i,j,k) = -inv_Pr_sgs*Nu_t(i,j,k)*dTdz(i,j,k)
    end do
    end do
    end do

    !$acc parallel loop collapse(2) present(pi_z, Nu_t, dTdz) async(1)
    do j = 1, ny
    do i = 1, ld
        pi_z(i,j,jz_max+1) = -inv_Pr_sgs*Nu_t(i,j,jz_max+1)                    &
            *dTdz(i,j,jz_max+1)
    end do
    end do

end subroutine scalars_flux_acc

!*******************************************************************************
subroutine scalars_rhs_theta_acc()
!*******************************************************************************
use param, only : lbz, ld, nx, ny, nz, coord, nproc, jt_total, dt, tadv1, tadv2, dz
implicit none

integer :: i, j, k
logical :: first_rhs_copy
real(rprec) :: rhs_new, rhs_old, inv_dz, divz, div_xy

! Fuse the final z-divergence contribution into the RHS update.  The scalar-GPU
! path stores d(pi_x)/dx + d(pi_y)/dy directly in div_pi; the pi_z stencil is
! evaluated here to avoid a separate ddz_w_gpu launch.
! Also fuse the theta update to avoid a second full-domain launch.  On the
! first initialized step the CPU path sets RHS_Tf = RHS_T after RHS_T is
! formed, so use the just-computed RHS for both AB terms in that case.
first_rhs_copy = (jt_total == 1) .and. (inits)
inv_dz = 1._rprec / dz
!$acc parallel loop collapse(3) present(theta, RHS_T, RHS_Tf, div_pi, pi_z) &
!$acc     async(1)                                                             &
!$acc     private(rhs_new, rhs_old, divz, div_xy)
do k = 1, nz-1
do j = 1, ny
do i = 1, nx
    div_xy = div_pi(i,j,k)
    divz = inv_dz*(pi_z(i,j,k+1) - pi_z(i,j,k))
    rhs_new = -RHS_T(i,j,k) - div_xy - divz
    RHS_T(i,j,k) = rhs_new
    if (first_rhs_copy) then
        RHS_Tf(i,j,k) = rhs_new
        rhs_old = rhs_new
    else
        rhs_old = RHS_Tf(i,j,k)
    end if
    theta(i,j,k) = theta(i,j,k) + dt*(tadv1*rhs_new + tadv2*rhs_old)
end do
end do
end do

if (coord == nproc-1) then
    !$acc parallel loop collapse(2) present(theta) async(1)
    do j = 1, ny
    do i = 1, ld
        theta(i,j,nz) = theta(i,j,nz-1) + lapse_rate*dz
    end do
    end do
end if

end subroutine scalars_rhs_theta_acc
#endif


    !*******************************************************************************
subroutine scalars_init
!*******************************************************************************
! This subroutine initializes the variables for the scalars module
use param, only : lbz, ld, ld_big, nx, ny, nz, ny2, u_star, z_i
use functions, only : count_lines
integer :: i, num_t, fid

! Allocate simulation variables
allocate ( theta(ld, ny, lbz:nz) ); theta = 0._rprec
allocate ( dTdx(ld, ny, lbz:nz) ); dTdx = 0._rprec
allocate ( dTdy(ld, ny, lbz:nz) ); dTdy = 0._rprec
allocate ( dTdz(ld, ny, lbz:nz) ); dTdz = 0._rprec
allocate ( RHS_T(ld, ny, lbz:nz) ); RHS_T = 0._rprec
allocate ( RHS_Tf(ld, ny, lbz:nz) ); RHS_Tf = 0._rprec
allocate ( u_big(ld_big, ny2, lbz:nz)); u_big = 0._rprec
allocate ( v_big(ld_big, ny2, lbz:nz)); v_big = 0._rprec
allocate ( w_big(ld_big, ny2, lbz:nz)); w_big = 0._rprec
allocate ( dTdx_big(ld_big, ny2, lbz:nz)); dTdx_big = 0._rprec
allocate ( dTdy_big(ld_big, ny2, lbz:nz)); dTdy_big = 0._rprec
allocate ( dTdz_big(ld_big, ny2, lbz:nz)); dTdz_big = 0._rprec
allocate ( RHS_big(ld_big, ny2, lbz:nz)); RHS_big = 0._rprec
allocate ( pi_x(ld, ny, lbz:nz) ); pi_x = 0._rprec
allocate ( pi_y(ld, ny, lbz:nz) ); pi_y = 0._rprec
allocate ( pi_z(ld, ny, lbz:nz) ); pi_z = 0._rprec
allocate ( div_pi(ld, ny, lbz:nz) ); div_pi = 0._rprec
allocate ( temp_var(ld, ny, lbz:nz) ); temp_var = 0._rprec
allocate ( Kappa_t(ld, ny, lbz:nz) ); Kappa_t = 0._rprec
#ifdef PPSCALARS_GPU
allocate ( theta_bar_acc(lbz:nz) ); theta_bar_acc = 0._rprec
!$acc enter data copyin(theta_bar_acc)
#endif

! Obukhov values (defaults for passive scalars)
allocate ( psi_m(nx, ny) ); psi_m = 0._rprec
allocate ( phi_m(nx, ny) ); phi_m = 1._rprec
allocate ( psi_h(nx, ny) ); psi_h = 0._rprec
allocate ( phi_h(nx, ny) ); phi_h = 1._rprec
allocate ( L(nx, ny) ); L = 0._rprec
    allocate ( tstar_lbc(nx, ny) ); tstar_lbc = 0._rprec

    ! Some GPU all-module smoke cases compile PPSCALARS without a SCALARS input
    ! block.  Keep those cases well-defined by using a neutral one-point profile;
    ! configured scalar cases still use the parsed IC_Z/IC_THETA vectors.
    if (.not. allocated(ic_z) .or. .not. allocated(ic_theta)) then
        if (allocated(ic_z)) deallocate(ic_z)
        if (allocated(ic_theta)) deallocate(ic_theta)
        ic_nloc = 1
        allocate(ic_z(ic_nloc), ic_theta(ic_nloc))
        ic_z = 0._rprec
        ic_theta = scal_bot
        ic_no_vel_noise_z = 0._rprec
    end if

    ! Nondimensionalize variables
    g = g*(z_i/(u_star**2))
flux_bot = flux_bot/u_star/T_scale
scal_bot = scal_bot/T_scale
lapse_rate = lapse_rate/T_scale*z_i
ic_theta = ic_theta/T_scale
ic_z = ic_z/z_i
ic_no_vel_noise_z = ic_no_vel_noise_z/z_i

! Read values from file
if (read_lbc_scal) then
    ! Count number of entries and allocate
    num_t = count_lines('lbc_scal.dat')
    allocate( t_interp(num_t) )
    allocate( lbc_interp(num_t) )

    ! Read entries
    open(newunit=fid, file='lbc_scal.dat', status='unknown', form='formatted', &
        position='rewind')
    do i = 1, num_t
        read(fid,*) t_interp(i), lbc_interp(i)
    end do
end if

#ifdef PPSCALARS_GPU
call scalars_update_device_state()
#endif

end subroutine scalars_init

!*******************************************************************************
subroutine ic_scal(interp_flag)
!*******************************************************************************
! Set initial profile for scalar
use param, only : coord
use string_util, only : string_concat
logical :: interp_flag

fname = path // 'scal.out'
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif
inquire (file=fname, exist=inits)

if (inits .and. .not.interp_flag) then
    inits = .true.
else
    inits = .false.
end if


if (inits) then
    write(*,*)  "--> Reading initial scalar field from file"
    call ic_scal_file
elseif (interp_flag) then
    write(*,*)  "--> Interpolating initial scalar field from file"
    call ic_scal_interp
else
    write(*,*)  "--> Creating initial boundary layer scalar field with LES BCs"
    call ic_scal_les
end if

#ifdef PPSCALARS_GPU
call scalars_update_device_state()
#endif

end subroutine ic_scal

!*******************************************************************************
subroutine ic_scal_file
!*******************************************************************************
! Read initial profile for scalar from file
use param, only : nx, nz, read_endian
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP

open(12, file=fname, form='unformatted', convert=read_endian)
read(12) theta(:, :, 1:nz), RHS_T(:, :, 1:nz), psi_m(1:nx, :)
close(12)

#ifdef PPMPI
call mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP)
call mpi_sync_real_array(RHS_T, 0, MPI_SYNC_DOWNUP)
#endif

end subroutine ic_scal_file

!*******************************************************************************
subroutine ic_scal_les
!*******************************************************************************
use param, only : lbz, nz
use grid_m, only : grid
use functions, only : linear_interp

integer :: i

do i = lbz, nz
    if (grid%z(i) < ic_z(ic_nloc)) then
        theta(:,:,i) = linear_interp(ic_z, ic_theta, grid%z(i))
    else
        theta(:,:,i) = ic_theta(ic_nloc) + lapse_rate*(grid%z(i) - ic_z(ic_nloc))
    end if
end do

end subroutine ic_scal_les

!*******************************************************************************
subroutine ic_scal_interp()
!*******************************************************************************
! This subroutine reads the initial conditions from a checkpoint file and
! interpolates onto the current grid
!
use param, only : nx, ny, nz, lbz, read_endian, path
use grid_m, only : grid
use functions, only : binary_search
integer :: nproc_f, Nx_f, Ny_f, Nz_f
real(rprec) :: Lx_f, Ly_f, Lz_f
integer :: i, j, k, z1, z2, ld_f, lh_f, Nz_tot_f
real(rprec) :: dx_f, dy_f, dz_f
integer :: i1, i2, j1, j2, k1, k2
real(rprec) :: ax, ay, az, bx, by, bz, xx, yy
real(rprec), allocatable, dimension(:) :: x_f, y_f!, z_f, zw_f
real(rprec), allocatable, dimension(:,:,:) :: theta_f
character(64) :: ff
integer :: npr1, npr2, nproc_r, Nz_tot_r
real(rprec), allocatable, dimension(:) :: z_r, zw_r
! real(rprec), allocatable, dimension(:,:,:) :: u_r, v_r, w_r

! Read grid information from file
open(12, file=path // 'grid.out', form='unformatted', convert=read_endian)
read(12) nproc_f, Nx_f, Ny_f, Nz_f, Lx_f, Ly_f, Lz_f
close(12)

! Compute intermediate values
lh_f = nx_f/2+1
ld_f = 2*lh_f
Nz_tot_f = ( nz_f - 1 ) * nproc_f + 1
dx_f = Lx_f / nx_f
dy_f = Ly_f / ny_f
dz_f = Lz_f / (nz_tot_f - 1)

! Figure out which processors actually need to be read
npr1 = max(floor(grid%z(lbz)/dz_f/Nz_f), 0)
npr2 = min(ceiling(grid%z(nz)/dz_f/Nz_f), nproc_f-1)
nproc_r = npr2-npr1+1
Nz_tot_r = ( nz_f - 1 ) * nproc_r + 1
! write(*,*) coord, npr1, npr2, nproc_r, Nz_tot_r

! Create file grid
allocate( z_r(nz_tot_r), zw_r(nz_tot_r))
do i = 1, nz_tot_r
    zw_r(i) = (i - 1 + npr1*(nz_f-1)) * dz_f
    z_r(i) = zw_r(i) + 0.5*dz_f
end do

! Create file grid
allocate( x_f(nx_f), y_f(ny_f))
do i = 1, nx_f
    x_f(i) = (i-1) * Lx_f/(nx_f)
end do
do i = 1, ny_f
    y_f(i) = (i-1) * Ly_f/ny_f
end do

! Read velocities from file
allocate( theta_f(ld_f, ny_f, nz_tot_r) )

! Loop through all levels
do i = 1, nproc_r
    ! Set level bounds
    z1 = nz_tot_f / nproc_f * (i-1) + 1
    z2 = nz_tot_f / nproc_f * i  + 1

    ! Read from file
    write(ff,*) i+npr1-1
    ff = path // "scal.out.c"//trim(adjustl(ff))
    open(12, file=ff,  action='read', form='unformatted')
    read(12) theta_f(:, :, z1:z2)
    close(12)
end do

! Calculate velocities
do i = 1, nx
    xx = grid%x(i) - floor(grid%x(i)/Lx_f) * Lx_f
    i1 = binary_search(x_f, xx)
    i2 = mod(i1, nx_f) + 1
    bx = (xx - x_f(i1)) / dx_f
    ax = 1._rprec - bx
    do j = 1, ny
        yy = grid%y(j) - floor(grid%y(j)/Ly_f) * Ly_f
        j1 = binary_search(y_f, yy)
        j2 = mod(j1, ny_f) + 1
        by = (yy - y_f(j1)) / dy_f
        ay = 1._rprec - by
        do k = 1, nz
            ! for u and v
            k1 = binary_search(z_r, grid%z(k))
            k2 = k1 + 1
            if (k1 == nz_tot_r) then
                theta(i,j,k) = ax*ay*theta_f(i1,j1,k1)                         &
                             + bx*ay*theta_f(i2,j1,k1)                         &
                             + ax*by*theta_f(i1,j2,k1)                         &
                             + bx*by*theta_f(i2,j2,k1)
            else if (k1 == 0) then
                theta(i,j,k) = ax*ay*theta_f(i1,j1,k2)                         &
                             + bx*ay*theta_f(i2,j1,k2)                         &
                             + ax*by*theta_f(i1,j2,k2)                         &
                             + bx*by*theta_f(i2,j2,k2)
            else
                bz = (grid%z(k) - z_r(k1)) / dz_f
                az = 1._rprec - bz
                theta(i,j,k) = ax*ay*az*theta_f(i1,j1,k1)                      &
                             + bx*ay*az*theta_f(i2,j1,k1)                      &
                             + ax*by*az*theta_f(i1,j2,k1)                      &
                             + bx*by*az*theta_f(i2,j2,k1)                      &
                             + ax*ay*bz*theta_f(i1,j1,k2)                      &
                             + bx*ay*bz*theta_f(i2,j1,k2)                      &
                             + ax*by*bz*theta_f(i1,j2,k2)                      &
                             + bx*by*bz*theta_f(i2,j2,k2)
            end if

        end do
    end do
end do

end subroutine ic_scal_interp

!*******************************************************************************
subroutine scalars_checkpoint
!*******************************************************************************
use param, only : nx, nz, write_endian

#ifdef PPSCALARS_GPU
! The scalar GPU path keeps theta/RHS_T resident during time stepping. Refresh
! host copies before writing Fortran checkpoint records.
!$acc wait(1)
!$acc update self(theta, RHS_T, psi_m)
#endif

!  Open scal.out (lun_default in io) for final output
open(11, file=fname, form='unformatted', convert=write_endian,                 &
    status='unknown', position='rewind')
write(11) theta(:, :, 1:nz), RHS_T(:, :, 1:nz), psi_m(1:nx, :)
close(11)

end subroutine scalars_checkpoint

!*******************************************************************************
subroutine scalars_deriv
!*******************************************************************************
use param, only : lbz, ld, ny, nz, coord, nproc
#if defined(PPSCALARS_GPU) && defined(PPGPU_AWARE_MPI)
use param, only : MPI_RPREC, down, up, comm, status, ierr
use mpi, only : mpi_sendrecv
#endif
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP
use derivatives, only : filt_da, ddz_uv

#ifdef PPSCALARS_GPU
integer :: i, j
logical :: scalar_acc

scalar_deriv_big_ready = .false.
scalar_acc = scalars_acc_enabled()

if (scalar_acc) then
    ! Reuse the accepted big-grid derivative path for decomposed runs too.
    ! It fills dTdx_big/dTdy_big/dTdz_big on device and avoids three later
    ! scalar-gradient to_big transforms in scalars_transport().  The small-grid
    ! dTdz halo below is still needed by the flux-divergence path.
    call scalars_deriv_xy_big_acc()
    call ddz_uv_gpu(theta, dTdz, lbz)

#ifdef PPMPI
    if (nproc > 1) then
        !$acc wait(1)
#ifdef PPGPU_AWARE_MPI
        !$acc host_data use_device(dTdz)
        call mpi_sendrecv(dTdz(1,1,1),    ld*ny, MPI_RPREC, down, 1,          &
                          dTdz(1,1,nz),   ld*ny, MPI_RPREC, up,   1, comm,   &
                          status, ierr)
        call mpi_sendrecv(dTdz(1,1,nz-1), ld*ny, MPI_RPREC, up,   2,          &
                          dTdz(1,1,0),    ld*ny, MPI_RPREC, down, 2, comm,   &
                          status, ierr)
        !$acc end host_data
#else
        !$acc update self(dTdz(:,:,1), dTdz(:,:,nz-1))
        call mpi_sync_real_array(dTdz, 0, MPI_SYNC_DOWNUP)
        !$acc update device(dTdz(:,:,0), dTdz(:,:,nz))
#endif
    end if
#endif

    if (coord == nproc-1) then
        !$acc parallel loop collapse(2) present(dTdz) async(1)
        do j = 1, ny
        do i = 1, ld
            dTdz(i,j,nz) = lapse_rate
        end do
        end do
    end if

    !$acc wait(1)
    return
end if
#endif

! Calculate derivatives of theta
call filt_da(theta, dTdx, dTdy, lbz)
call ddz_uv(theta, dTdz, lbz)

#ifdef PPMPI
call mpi_sync_real_array(dTdz, 0, MPI_SYNC_DOWNUP)
#endif

! Top boundary condition
if (coord == nproc-1) dTdz(:,:,nz) = lapse_rate

#ifdef PPSCALARS_GPU
if (scalar_acc) then
    !$acc update device(theta, dTdx, dTdy, dTdz)
end if
#endif

end subroutine scalars_deriv

!*******************************************************************************
subroutine obukhov(u_avg)
!*******************************************************************************
use param, only : vonk, dz, zo, nx, ny, ld, u_star, lbz, total_time_dim, pi
use sim_param, only : ustar_lbc
use coriolis, only : repeat_interval
use functions, only : linear_interp
use test_filtermodule, only : test_filter

real(rprec), dimension(nx, ny), intent(in) :: u_avg
real(rprec), dimension(ld, ny) :: theta1

integer :: i, j
#if defined(PPSCALARS_GPU)
real(rprec) :: zeta_acc, zeta_zo_acc, y_acc, yy_acc, xx_acc, psi_zero_acc
real(rprec) :: psi_m_half_acc, psi_m_zo_acc, psi_h_half_acc, psi_h_zo_acc
real(rprec), parameter :: am_acc=6.1_rprec, bm_acc=2.5_rprec
real(rprec), parameter :: ah_acc=5.3_rprec, bh_acc=1.1_rprec
real(rprec), parameter :: a_acc=0.33_rprec, b_acc=0.41_rprec
real(rprec), parameter :: c_acc=0.33_rprec, d_acc=0.057_rprec
real(rprec), parameter :: n_acc=0.78_rprec
#endif
#if defined(PPSCALARS_GPU)
logical, save :: obukhov_acc_allocated = .false.
real(rprec), save, allocatable, dimension(:,:) :: theta1_acc
#endif

#if defined(PPSCALARS_GPU)
if (scalars_acc_enabled()) then
    if (.not. obukhov_acc_allocated) then
        allocate(theta1_acc(ld, ny))
        !$acc enter data create(theta1_acc)
        obukhov_acc_allocated = .true.
    end if

    if (passive_scalar) then
        !$acc parallel loop collapse(2) default(present) async(1)
        do j = 1, ny
        do i = 1, nx
            ustar_lbc(i,j) = u_avg(i,j)*vonk/log(0.5_rprec*dz/zo)
        end do
        end do
        return
    end if

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 1, ny
    do i = 1, ld
        theta1_acc(i,j) = theta(i,j,1)
    end do
    end do

    call test_filter_plane_gpu(theta1_acc)

    if (read_lbc_scal) then
        if (lbc_scal == 0) then
            scal_bot = linear_interp(t_interp, lbc_interp,                    &
                mod(total_time_dim, repeat_interval))/T_scale
        else
            flux_bot = linear_interp(t_interp, lbc_interp,                    &
                mod(total_time_dim, repeat_interval))/T_scale/u_star
        end if
    end if

    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 1, ny
    do i = 1, nx
        ustar_lbc(i,j) = u_avg(i,j)*vonk/(log(0.5_rprec*dz/zo)                &
            + psi_m(i,j))
        tstar_lbc(i,j) = (theta1_acc(i,j) - scal_bot)*vonk                   &
            / (log(0.5_rprec*dz/zo_s) + psi_h(i,j))
        L(i,j) = ustar_lbc(i,j)*ustar_lbc(i,j)*theta1_acc(i,j)               &
            / (vonk*g*tstar_lbc(i,j))
        if (abs(L(i,j)) > 1._rprec/epsilon(0._rprec)) then
            phi_m(i,j) = 1._rprec
            phi_h(i,j) = 1._rprec
            psi_m(i,j) = 0._rprec
            psi_h(i,j) = 0._rprec
        else
            zeta_acc = 0.5_rprec*dz/L(i,j)
            if (zeta_acc < 0._rprec) then
                y_acc = -zeta_acc
                if (y_acc > b_acc**(-3._rprec)) then
                    yy_acc = b_acc**(-3._rprec)
                    xx_acc = (yy_acc/a_acc)**(1._rprec/3._rprec)
                    psi_zero_acc = -log(a_acc) + sqrt(3._rprec)*b_acc         &
                        *(a_acc**(1._rprec/3._rprec))*(pi/6._rprec)
                    psi_m_half_acc = log(a_acc+yy_acc) - 3._rprec*b_acc       &
                        *(yy_acc**(1._rprec/3._rprec)) + (b_acc              &
                        *(a_acc**(1._rprec/3._rprec))/2._rprec)              &
                        *log(((1._rprec+xx_acc)**2)/(1._rprec - xx_acc       &
                        + xx_acc*xx_acc)) + sqrt(3._rprec)*b_acc             &
                        *(a_acc**(1._rprec/3._rprec))*atan(((2._rprec        &
                        *xx_acc)-1._rprec)/sqrt(3._rprec)) + psi_zero_acc
                else
                    xx_acc = (y_acc/a_acc)**(1._rprec/3._rprec)
                    psi_zero_acc = -log(a_acc) + sqrt(3._rprec)*b_acc         &
                        *(a_acc**(1._rprec/3._rprec))*(pi/6._rprec)
                    psi_m_half_acc = log(a_acc+y_acc) - 3._rprec*b_acc        &
                        *(y_acc**(1._rprec/3._rprec)) + (b_acc               &
                        *(a_acc**(1._rprec/3._rprec))/2._rprec)              &
                        *log(((1._rprec+xx_acc)**2)/(1._rprec - xx_acc       &
                        + xx_acc*xx_acc)) + sqrt(3._rprec)*b_acc             &
                        *(a_acc**(1._rprec/3._rprec))*atan(((2._rprec        &
                        *xx_acc)-1._rprec)/sqrt(3._rprec)) + psi_zero_acc
                end if
                psi_h_half_acc = ((1._rprec-d_acc)/n_acc)                    &
                    * log((c_acc+y_acc**n_acc)/c_acc)
            else
                psi_m_half_acc = -am_acc*log(zeta_acc                        &
                    + (1._rprec + zeta_acc**2.5_rprec)**0.4_rprec)
                psi_h_half_acc = -ah_acc*log(zeta_acc                        &
                    + (1._rprec + zeta_acc**bh_acc)**(1._rprec/bh_acc))
            end if

            zeta_zo_acc = zo/L(i,j)
            if (zeta_zo_acc < 0._rprec) then
                y_acc = -zeta_zo_acc
                if (y_acc > b_acc**(-3._rprec)) then
                    yy_acc = b_acc**(-3._rprec)
                    xx_acc = (yy_acc/a_acc)**(1._rprec/3._rprec)
                    psi_zero_acc = -log(a_acc) + sqrt(3._rprec)*b_acc         &
                        *(a_acc**(1._rprec/3._rprec))*(pi/6._rprec)
                    psi_m_zo_acc = log(a_acc+yy_acc) - 3._rprec*b_acc         &
                        *(yy_acc**(1._rprec/3._rprec)) + (b_acc              &
                        *(a_acc**(1._rprec/3._rprec))/2._rprec)              &
                        *log(((1._rprec+xx_acc)**2)/(1._rprec - xx_acc       &
                        + xx_acc*xx_acc)) + sqrt(3._rprec)*b_acc             &
                        *(a_acc**(1._rprec/3._rprec))*atan(((2._rprec        &
                        *xx_acc)-1._rprec)/sqrt(3._rprec)) + psi_zero_acc
                else
                    xx_acc = (y_acc/a_acc)**(1._rprec/3._rprec)
                    psi_zero_acc = -log(a_acc) + sqrt(3._rprec)*b_acc         &
                        *(a_acc**(1._rprec/3._rprec))*(pi/6._rprec)
                    psi_m_zo_acc = log(a_acc+y_acc) - 3._rprec*b_acc          &
                        *(y_acc**(1._rprec/3._rprec)) + (b_acc               &
                        *(a_acc**(1._rprec/3._rprec))/2._rprec)              &
                        *log(((1._rprec+xx_acc)**2)/(1._rprec - xx_acc       &
                        + xx_acc*xx_acc)) + sqrt(3._rprec)*b_acc             &
                        *(a_acc**(1._rprec/3._rprec))*atan(((2._rprec        &
                        *xx_acc)-1._rprec)/sqrt(3._rprec)) + psi_zero_acc
                end if
                psi_h_zo_acc = ((1._rprec-d_acc)/n_acc)                      &
                    * log((c_acc+y_acc**n_acc)/c_acc)
            else
                psi_m_zo_acc = -am_acc*log(zeta_zo_acc                       &
                    + (1._rprec + zeta_zo_acc**2.5_rprec)**0.4_rprec)
                psi_h_zo_acc = -ah_acc*log(zeta_zo_acc                       &
                    + (1._rprec + zeta_zo_acc**bh_acc)**(1._rprec/bh_acc))
            end if

            psi_m(i,j) = -psi_m_half_acc + psi_m_zo_acc
            psi_h(i,j) = -psi_h_half_acc + psi_h_zo_acc
            if (zeta_acc < 0._rprec) then
                y_acc = -zeta_acc
                if (y_acc > b_acc**(-3._rprec)) then
                    phi_m(i,j) = 1._rprec
                else
                    phi_m(i,j) = (a_acc + b_acc                              &
                        *(y_acc**(4._rprec/3._rprec)))/(a_acc+y_acc)
                end if
                phi_h(i,j) = (c_acc + d_acc*(y_acc**n_acc))                  &
                    / (c_acc + y_acc**n_acc)
            else
                phi_m(i,j) = 1._rprec + am_acc*(zeta_acc + zeta_acc**bm_acc  &
                    *((1._rprec + zeta_acc**bm_acc)**(-1._rprec              &
                    + 1._rprec/bm_acc))) / (zeta_acc                         &
                    + ((1._rprec + zeta_acc**bm_acc)**(1._rprec/bm_acc)))
                phi_h(i,j) = 1._rprec + ah_acc*(zeta_acc + zeta_acc**bh_acc  &
                    *((1._rprec + zeta_acc**bh_acc)**(-1._rprec              &
                    + 1._rprec/bh_acc))) / (zeta_acc                         &
                    + ((1._rprec + zeta_acc**bh_acc)**(1._rprec/bh_acc)))
            end if
        end if
        ustar_lbc(i,j) = u_avg(i,j)*vonk/(log(0.5_rprec*dz/zo)                &
            + psi_m(i,j))
        if (lbc_scal == 0) then
            tstar_lbc(i,j) = (theta1_acc(i,j) - scal_bot)*vonk               &
                / (log(0.5_rprec*dz/zo_s) + psi_h(i,j))
        else
            tstar_lbc(i,j) = -flux_bot/ustar_lbc(i,j)
        end if
        dTdz(i,j,1) = tstar_lbc(i,j)/(vonk*dz*0.5_rprec)*phi_h(i,j)
        pi_z(i,j,1) = -tstar_lbc(i,j)*ustar_lbc(i,j)
    end do
    end do
    return
end if
#endif


! Use previous ustar_lbc to compute stability correction
if (passive_scalar) then
    ustar_lbc = u_avg*vonk/log(0.5_rprec*dz/zo)
    return
end if

theta1 = theta(:,:,1)
call test_filter(theta1)

! Using previous time step's psi_m and psi_h to calculate obukhov length and
! stability functions
ustar_lbc = u_avg*vonk/(log(0.5_rprec*dz/zo) + psi_m)
tstar_lbc = (theta1(1:nx,:) - scal_bot)*vonk                                   &
    / (log(0.5_rprec*dz/zo_s) + psi_h)

L = ustar_lbc**2*theta1(1:nx,:)/(vonk*g*tstar_lbc)
do i = 1, nx
    do j = 1, ny
        call stability(L(i,j), zo_s, phi_m(i,j), phi_h(i,j), psi_m(i,j),       &
            psi_h(i,j))
    end do
end do

! Recompute ustar_lbc using new values
ustar_lbc = u_avg*vonk/(log(0.5_rprec*dz/zo) + psi_m)

! Get boundary condition if reading from file
if (read_lbc_scal) then
    if (lbc_scal == 0) then
        scal_bot = linear_interp(t_interp, lbc_interp,                         &
            mod(total_time_dim, repeat_interval))/T_scale
    else
        flux_bot = linear_interp(t_interp, lbc_interp,                         &
            mod(total_time_dim, repeat_interval))/T_scale/u_star
    end if
end if

! Calculate tstar_lbc based on boundary condition
if (lbc_scal == 0) then
    tstar_lbc = (theta1(1:nx,:) - scal_bot)*vonk                               &
        / (log(0.5_rprec*dz/zo_s) + psi_h)
else
    tstar_lbc = -flux_bot/ustar_lbc
end if

! Calculate temperature gradient and flux
dTdz(1:nx,:,1) = tstar_lbc/(vonk*dz*0.5_rprec)*phi_h
pi_z(1:nx,:,1) = -tstar_lbc*ustar_lbc

end subroutine obukhov

!*******************************************************************************
subroutine scalars_transport()
!*******************************************************************************
use param, only : lbz, nx, nz, nx2, ny2, nproc, coord, dt, tadv1, tadv2,       &
    jt_total, wbase
use param, only : lbc_mom, ubc_mom, dz
#if defined(PPSCALARS_GPU) && defined(PPGPU_AWARE_MPI)
use param, only : ld, ny, MPI_RPREC, down, up, comm, status, ierr
use mpi, only : mpi_sendrecv
#endif
use sim_param, only : u, v, w
use sgs_param, only : Nu_t
use derivatives, only : filt_da, ddx, ddy, ddz_uv, ddz_w
use mpi_defs, only :  mpi_sync_real_array, MPI_SYNC_DOWNUP
use fft
use messages, only : error
#if defined(PPSCALARS_GPU) && defined(PPCONVEC_GPU)
use convec_gpu_m, only : convec_gpu_big_available,                            &
    convec_u_big => u_big_d, convec_v_big => v_big_d,                         &
    convec_w_big => w_big_d
#endif

    integer :: k, jz_min, jz_max
    real(rprec) :: const
#if defined(PPSCALARS_GPU)
    real(rprec) :: scalar_t0, scalar_t_stage, scalar_t_div0
    logical :: scalar_timing
#endif
#ifdef PPSCALARS_GPU
    logical :: scalar_acc
    logical :: reuse_convec_big

    scalar_acc = scalars_acc_enabled()
#ifdef PPCONVEC_GPU
    reuse_convec_big = scalar_acc .and. convec_gpu_big_available()
#else
    reuse_convec_big = .false.
#endif
    scalar_timing = scalar_acc .and. scalars_stage_timing_enabled()
    if (scalar_timing) then
        scalars_stage_count = scalars_stage_count + 1
        call scalars_timer_start(scalar_t0)
        scalar_t_stage = scalar_t0
    end if
#endif

! We do not advance the ground nodes, so start at k=2.
! For the MPI case, the means that we start from jz=2
! for coord=0 and jz=1 otherwise.
#ifdef PPMPI
    if (coord == 0) then
        jz_min = 2
    else
        jz_min = 1
    end if
    if (coord == nproc-1) then
        jz_max = nz-2
    else
        jz_max = nz-1
    end if
#else
    jz_max = nz-2
    jz_min = 2
#endif

    ! Save previous timestep's RHS
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
        call scalars_copy_rhs_acc()
    else
#endif
    RHS_Tf = RHS_T
#ifdef PPSCALARS_GPU
    endif
#endif
#if defined(PPSCALARS_GPU)
    if (scalar_timing) call scalars_timer_accum(scalar_t_stage,                &
        scalars_time_copy_rhs, 'scalar copy RHS')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advective term u_i d_i \theta computed using dealiasing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We could save memory and computional time by saving u_big, v_big, w_big from
! the convec subroutine.
! (put in sim_param, not as a saved variable in the subroutine)

! Set variables onto big domain for multiplication in physical space and
! dealiasing.  The CUDA path uses the same batched cuFFT/pad/unpad pattern as
! convec.f90; the original FFTW loop remains the CPU fallback.
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
        if (.not. reuse_convec_big) then
            call scalars_to_big_acc(u, u_big)
            call scalars_to_big_acc(v, v_big)
            call scalars_to_big_acc(w, w_big)
        end if
        if (.not. scalar_deriv_big_ready) then
            call scalars_to_big_acc(dTdx, dTdx_big)
            call scalars_to_big_acc(dTdy, dTdy_big)
            call scalars_to_big_acc(dTdz, dTdz_big)
        end if
    else
#endif
    call to_big(u, u_big)
    call to_big(v, v_big)
    call to_big(w, w_big)
    call to_big(dTdx, dTdx_big)
    call to_big(dTdy, dTdy_big)
    call to_big(dTdz, dTdz_big)
#ifdef PPSCALARS_GPU
    endif
#endif
#if defined(PPSCALARS_GPU)
    if (scalar_timing) call scalars_timer_accum(scalar_t_stage,                &
        scalars_time_to_big, 'scalar to_big transforms')
#endif

! Normalization for FFTs
const=1._rprec/(nx2*ny2)

    ! Interior and boundary planes of domain
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
#ifdef PPCONVEC_GPU
        if (reuse_convec_big) then
            call scalars_advective_acc(jz_min, jz_max, const,                  &
                convec_u_big, convec_v_big, convec_w_big)
        else
#endif
            call scalars_advective_acc(jz_min, jz_max, const, u_big, v_big,    &
                w_big)
#ifdef PPCONVEC_GPU
        end if
#endif
    else
#endif
    do k = jz_min, jz_max
        RHS_big(:,:,k) = const*(u_big(:,:,k)*dTdx_big(:,:,k)                       &
            + v_big(:,:,k)*dTdy_big(:,:,k)                                         &
            + 0.5_rprec*w_big(:,:,k+1)*dTdz_big(:,:,k+1)                           &
            + 0.5_rprec*w_big(:,:,k)*dTdz_big(:,:,k))
    end do

    if (coord == 0) then
        RHS_big(:,:,1) = const*(u_big(:,:,1)*dTdx_big(:,:,1)                       &
            + v_big(:,:,1)*dTdy_big(:,:,1)                                         &
            + 0.5_rprec*w_big(:,:,2)*dTdz_big(:,:,2))
    end if

    if (coord == nproc-1) then
        RHS_big(:,:,nz-1) = const*(u_big(:,:,nz-1)*dTdx_big(:,:,nz-1)              &
            + v_big(:,:,nz-1)*dTdy_big(:,:,nz-1)                                   &
            + 0.5_rprec*w_big(:,:,nz-1)*dTdz_big(:,:,nz-1))
    end if
#ifdef PPSCALARS_GPU
    endif
#endif
#if defined(PPSCALARS_GPU)
    if (scalar_timing) call scalars_timer_accum(scalar_t_stage,                &
        scalars_time_advective, 'scalar advective product')
#endif

! Put back on the smaller grid
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
        call scalars_return_rhs_acc()
    else
#endif
    do k = 1, nz-1
        call dfftw_execute_dft_r2c(forw_big, RHS_big(:,:,k), RHS_big(:,:,k))
        call unpadd(RHS_T(:,:,k), RHS_big(:,:,k))
        call dfftw_execute_dft_c2r(back, RHS_T(:,:,k), RHS_T(:,:,k))
    end do
#ifdef PPSCALARS_GPU
    endif
#endif
#if defined(PPSCALARS_GPU)
    if (scalar_timing) call scalars_timer_accum(scalar_t_stage,                &
        scalars_time_return, 'scalar return transform')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subgrid stress
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Calculate eddy diffusivity and heat fluxes
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
        call scalars_flux_acc(jz_min, jz_max)
    else
#endif
    Kappa_t = Nu_t/Pr_sgs

    if (coord == 0) then
        select case (lbc_mom)
        case (0)
            pi_x(:,:,1) = -0.5_rprec*(Kappa_t(:,:,1) + Kappa_t(:,:,2))*dTdx(:,:,1)
            pi_y(:,:,1) = -0.5_rprec*(Kappa_t(:,:,1) + Kappa_t(:,:,2))*dTdy(:,:,1)
        case (1:)
            pi_x(:,:,1) = -Kappa_t(:,:,1)*dTdx(:,:,1)
            pi_y(:,:,1) = -Kappa_t(:,:,1)*dTdy(:,:,1)
        end select
    end if

    if (coord == nproc-1) then
        select case (ubc_mom)
        case (0)
            pi_x(:,:,nz-1) = -0.5_rprec*(Kappa_t(:,:,nz-1) + Kappa_t(:,:,nz))*dTdx(:,:,nz-1)
            pi_y(:,:,nz-1) = -0.5_rprec*(Kappa_t(:,:,nz-1) + Kappa_t(:,:,nz))*dTdy(:,:,nz-1)
            pi_z(:,:,nz) = -Kappa_t(:,:,nz-1)*dTdz(:,:,nz)
        case (1:)
            pi_x(:,:,nz-1) = -Kappa_t(:,:,nz-1)*dTdx(:,:,nz-1)
            pi_y(:,:,nz-1) = -Kappa_t(:,:,nz-1)*dTdy(:,:,nz-1)
            pi_z(:,:,nz) = -Kappa_t(:,:,nz-1)*dTdz(:,:,nz)
        end select
    end if

    do k= jz_min, jz_max
        pi_x(:,:,k) = -0.5_rprec*(Kappa_t(:,:,k) + Kappa_t(:,:,k+1))*dTdx(:,:,k)
        pi_y(:,:,k) = -0.5_rprec*(Kappa_t(:,:,k) + Kappa_t(:,:,k+1))*dTdy(:,:,k)
        pi_z(:,:,k) = -Kappa_t(:,:,k)*dTdz(:,:,k)
    end do
    pi_z(:,:,jz_max+1) = -Kappa_t(:,:,jz_max+1)*dTdz(:,:,jz_max+1)
#ifdef PPSCALARS_GPU
    endif
#endif
#if defined(PPSCALARS_GPU)
    if (scalar_timing) then
        call scalars_timer_accum(scalar_t_stage, scalars_time_flux,             &
            'scalar flux build')
        scalar_t_div0 = scalar_t_stage
    end if
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Divergence of heat flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Store the derivatives in the stress values...can change if we need to output
! stuff
#ifdef PPSCALARS_GPU
if (scalar_acc) then
    call scalars_divergence_acc()
else
#endif
call ddx(pi_x, div_pi, lbz)
call ddy(pi_y, temp_var, lbz)
div_pi = div_pi + temp_var
call ddz_w(pi_z, temp_var, lbz)
#ifdef PPSCALARS_GPU
end if
#endif
#if defined(PPSCALARS_GPU)
    if (scalar_timing) call scalars_timer_accum(scalar_t_stage,                &
        scalars_time_divergence, 'scalar divergence')
#endif
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
        call scalars_rhs_theta_acc()
    else
#endif
    div_pi = div_pi + temp_var

    do k = 1, nz-1
        RHS_T(1:nx,:,k) = -RHS_T(1:nx,:,k) - div_pi(1:nx,:,k)
    end do

    if ((jt_total == 1) .and. (inits)) then
        RHS_Tf = RHS_T
    end if

    theta(1:nx,:,1:nz-1) = theta(1:nx,:,1:nz-1)                                    &
        + dt*(tadv1*RHS_T(1:nx,:,1:nz-1) + tadv2*RHS_Tf(1:nx,:,1:nz-1))
#ifdef PPSCALARS_GPU
    endif
#endif

#if defined(PPSCALARS_GPU)
    if (scalar_timing) call scalars_timer_accum(scalar_t_stage,                &
        scalars_time_rhs_update, 'scalar rhs/theta update')
#endif

#ifdef PPMPI
#ifdef PPSCALARS_GPU
if (scalar_acc .and. nproc == 1) then
    ! Single-rank scalar-GPU cases do not need a halo exchange.  Avoid the
    ! full theta device->host->device round trip that dominated 128^3 timing.
else
#ifdef PPGPU_AWARE_MPI
    if (scalar_acc) then
        ! GPU-aware scalar halo.  This matches
        ! mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP):
        !   tag 1 sends k=1 down and receives k=nz from up;
        !   tag 2 sends k=nz-1 up and receives k=0 from down.
        !$acc wait(1)
        !$acc host_data use_device(theta)
        call mpi_sendrecv(theta(1,1,1),    ld*ny, MPI_RPREC, down, 1,          &
                          theta(1,1,nz),   ld*ny, MPI_RPREC, up,   1, comm,   &
                          status, ierr)
        call mpi_sendrecv(theta(1,1,nz-1), ld*ny, MPI_RPREC, up,   2,          &
                          theta(1,1,0),    ld*ny, MPI_RPREC, down, 2, comm,   &
                          status, ierr)
        !$acc end host_data
    else
#endif
    if (scalar_acc) then
        !$acc wait(1)
        !$acc update self(theta)
    end if
#endif
call mpi_sync_real_array(theta, 0, MPI_SYNC_DOWNUP)
#ifdef PPSCALARS_GPU
    if (scalar_acc) then
        !$acc update device(theta)
    end if
#ifdef PPGPU_AWARE_MPI
    end if
#endif
end if
#endif
#endif

! Use gradient at top to project temperature above domain
#ifdef PPSCALARS_GPU
    if (.not. scalar_acc) then
#endif
    if (coord == nproc-1) then
            theta(:,:,nz) = theta(:,:,nz-1) + lapse_rate*dz
    end if
#ifdef PPSCALARS_GPU
    end if
#endif

#if defined(PPSCALARS_GPU)
    if (scalar_timing) then
        if (scalar_acc) then
            !$acc wait(1)
        end if
        call scalars_timer_accum(scalar_t_stage, scalars_time_halo,            &
            'scalar halo/top bc')
        scalars_time_div_update = scalars_time_div_update +                    &
            max(scalar_t_stage - scalar_t_div0, 0._rprec)
        scalars_time_total = scalars_time_total +                              &
            max(scalar_t_stage - scalar_t0, 0._rprec)
        call scalars_stage_report(wbase)
    end if
#endif

end subroutine scalars_transport

!*******************************************************************************
subroutine to_big(a, a_big)
!*******************************************************************************
use fft
use param, only : lbz, nx, ny, nz

real(rprec), dimension(ld, ny, lbz:nz), intent(inout) ::  a
real(rprec), dimension(ld_big, ny2, lbz:nz), intent(inout) :: a_big

integer :: jz
real(rprec) :: const

! Set variables onto big domain for multiplication in physical space and
! dealiasing
const = 1._rprec/(nx*ny)
do jz = lbz, nz
    temp_var(:,:,jz) = const*a(:,:,jz)
    call dfftw_execute_dft_r2c(forw, temp_var(:,:,jz), temp_var(:,:,jz))
    call padd(a_big(:,:,jz), temp_var(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, a_big(:,:,jz), a_big(:,:,jz))
end do

end subroutine to_big

!*******************************************************************************
subroutine buoyancy_force
!*******************************************************************************
! This subroutine calculates the buoyancy term due to temperature to be added to
! the RHS of the vertical momentum equation.
use param, only : coord, nx, ny, nz
use sim_param, only :  RHSz

integer :: i, j, k, jz_min
real(rprec) :: theta_bar
#if defined(PPSCALARS_GPU)
real(rprec) :: theta_sum
#endif

! We do not advance the ground nodes, so start at k=2.
! For the MPI case, the means that we start from jz=2
! for coord=0 and jz=1 otherwise.
#ifdef PPMPI
   if (coord == 0) then
      jz_min = 2
   else
      jz_min = 1
   end if
#else
   jz_min = 2
#endif

    ! Add to RHSz
    if ( .not.passive_scalar ) then
#if defined(PPSCALARS_GPU)
        if (scalars_acc_enabled()) then
            !$acc parallel loop default(present) async(1)
            do k = jz_min, nz-1
                theta_bar_acc(k) = 0._rprec
            end do
            !$acc parallel loop gang collapse(2) default(present)              &
            !$acc     private(theta_sum) async(1)
            do k = jz_min, nz-1
            do j = 1, ny
                theta_sum = 0._rprec
                !$acc loop vector reduction(+:theta_sum)
                do i = 1, nx
                    theta_sum = theta_sum +                                  &
                        0.5_rprec*(theta(i,j,k) + theta(i,j,k-1))
                end do
                !$acc atomic update
                theta_bar_acc(k) = theta_bar_acc(k) + theta_sum
            end do
            end do
            !$acc parallel loop default(present) async(1)
            do k = jz_min, nz-1
                theta_bar_acc(k) = theta_bar_acc(k)/nx/ny
            end do

            !$acc parallel loop collapse(3) default(present) async(1)
            do k = jz_min, nz-1
            do j = 1, ny
            do i = 1, nx
                RHSz(i,j,k) = RHSz(i,j,k) + g*(0.5_rprec*(theta(i,j,k)        &
                    + theta(i,j,k-1)) - theta_bar_acc(k))
            end do
            end do
            end do
        else
#endif
        do k = jz_min, nz-1
            theta_bar = sum(0.5_rprec*(theta(1:nx,:,k)+theta(1:nx,:,k-1)))/nx/ny
            RHSz(1:nx,:,k) = RHSz(1:nx,:,k) + g*(0.5_rprec*(theta(1:nx,:,k)+theta(1:nx,:,k-1)) - theta_bar)
        end do
#if defined(PPSCALARS_GPU)
        end if
#endif
    end if

end subroutine buoyancy_force

end module scalars
