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

!*******************************************************************************
subroutine wallstress
!*******************************************************************************
!
! This subroutine calculates the wall stress txz, tyz (w-nodes) and dudz,
! dvdz (w-nodes) at the first z-location k = 1. The wall stress is calculated
! depending on lower boundary condition lbc_mom. This subroutine should only
! be called after ensuring coord==0
!
! Options for lbc_mom:
!   0 - stress free
!       txz, tyz, dudz, and dvdz are all 0
!
!   1 - DNS wall boundary conditions
!       calculates wall stress values from the first grid point
!
!   2 - Equilibirum wall model
!       See John D. Albertson's dissertation, eqns (2.46)-(2.52)
!       Also see E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent
!           Lagrangian dynamic model for large eddy simulation of complex
!           turbulent flows" (2005) -- Appendix
!
!   3 - Integral wall model
!       See X.I.A. Yang, J. Sadique, R. Mittal & C. Meneveau, "Integral wall
!           model for large eddy simulations of wall-bounded turbulent flows." (2015)
!
use types, only : rprec
use param, only : lbc_mom
use param, only : ubc_mom, coord, nproc, nz ! these necessary only for upper bc
use messages, only : error
use iwmles, only : iwm_wallstress
use sim_param, only : txz, tyz, dudz, dvdz
implicit none
character(*), parameter :: sub_name = 'wallstress'

! Lower boundary condition
if (coord == 0) then
    select case (lbc_mom)
        ! Stress free
        case (0)
            call ws_free_lbc

        ! DNS wall
        case (1)
            call ws_dns_lbc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_lbc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call iwm_wallstress()

        ! Otherwise, invalid
        case default
            call error (sub_name, 'invalid lbc_mom')
    end select
end if

if (coord == nproc-1) then
    select case (ubc_mom)
        ! Stress free
        case (0)
            call ws_free_ubc

        ! DNS wall
        case (1)
            call ws_dns_ubc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_ubc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call error(sub_name, 'invalid ubc_mom')

        ! Otherwise, invalid
        case default
            call error(sub_name, 'invalid ubc_mom')
    end select
end if

contains


!*******************************************************************************
subroutine ws_free_lbc
!*******************************************************************************
use param, only : nx, ny
implicit none
integer :: i, j

#if defined(PPSGS_GPU)
! Device-resident wall plane (explicit-residency build): same arithmetic as
! the host loop below, on the async(1) queue so it orders with the SGS kernels.
!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny
do i = 1, nx
    txz(i,j,1) = 0._rprec
    tyz(i,j,1) = 0._rprec
    dudz(i,j,1) = 0._rprec
    dvdz(i,j,1) = 0._rprec
end do
end do
return
#endif

txz(:, :, 1) = 0._rprec
tyz(:, :, 1) = 0._rprec
dudz(:, :, 1) = 0._rprec
dvdz(:, :, 1) = 0._rprec

end subroutine ws_free_lbc

!*******************************************************************************
subroutine ws_free_ubc
!*******************************************************************************
use param, only : nx, ny
implicit none
integer :: i, j

#if defined(PPSGS_GPU)
!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny
do i = 1, nx
    txz(i,j,nz) = 0._rprec
    tyz(i,j,nz) = 0._rprec
    dudz(i,j,nz) = 0._rprec
    dvdz(i,j,nz) = 0._rprec
end do
end do
return
#endif

txz(:, :,nz) = 0._rprec
tyz(:, :,nz) = 0._rprec
dudz(:,:,nz) = 0._rprec
dvdz(:,:,nz) = 0._rprec

end subroutine ws_free_ubc

!*******************************************************************************
subroutine ws_dns_lbc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, dz
use param, only : ubot
use sim_param , only : u, v
implicit none
integer :: i, j

#if defined(PPSGS_GPU)
!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny
do i = 1, nx
    dudz(i,j,1) = ( u(i,j,1) - ubot ) / (0.5_rprec*dz)
    dvdz(i,j,1) = v(i,j,1) / (0.5_rprec*dz)
    txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
    tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)
end do
end do
return
#endif

do j = 1, ny
    do i = 1, nx
        dudz(i,j,1) = ( u(i,j,1) - ubot ) / (0.5_rprec*dz)
        dvdz(i,j,1) = v(i,j,1) / (0.5_rprec*dz)
        txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
        tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)
    end do
end do

end subroutine ws_dns_lbc

!*******************************************************************************
subroutine ws_dns_ubc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, dz
use param, only : utop
use sim_param , only : u, v
implicit none
integer :: i, j

#if defined(PPSGS_GPU)
!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny
do i = 1, nx
    dudz(i,j,nz) = ( utop - u(i,j,nz-1) ) / (0.5_rprec*dz)
    dvdz(i,j,nz) = -v(i,j,nz-1) / (0.5_rprec*dz)
    txz(i,j,nz) = -nu_molec/(z_i*u_star)*dudz(i,j,nz)
    tyz(i,j,nz) = -nu_molec/(z_i*u_star)*dvdz(i,j,nz)
end do
end do
return
#endif

do j = 1, ny
    do i = 1, nx
        dudz(i,j,nz) = ( utop - u(i,j,nz-1) ) / (0.5_rprec*dz)
        dvdz(i,j,nz) = -v(i,j,nz-1) / (0.5_rprec*dz)
        txz(i,j,nz) = -nu_molec/(z_i*u_star)*dudz(i,j,nz)
        tyz(i,j,nz) = -nu_molec/(z_i*u_star)*dvdz(i,j,nz)
    end do
end do

end subroutine ws_dns_ubc

!*******************************************************************************
subroutine ws_equilibrium_lbc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v, ustar_lbc
use test_filtermodule, only : test_filter
#ifdef PPSCALARS
use scalars, only : obukhov, phi_m, passive_scalar
#endif

implicit none

integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const
#if defined(PPSGS_GPU)
logical, save :: eq_lbc_acc_alloc = .false.
real(rprec), save, allocatable, dimension(:,:) :: u1_acc, v1_acc, u_avg_acc
real(rprec) :: u_avg_pt
#endif

#if defined(PPSGS_GPU)
! Device-resident equilibrium wall model for explicit-residency builds.
! In scalar builds, obukhov() handles both passive and active stability cases.
if (.not. eq_lbc_acc_alloc) then
    allocate(u1_acc(ld,ny), v1_acc(ld,ny), u_avg_acc(nx,ny))
    !$acc enter data create(u1_acc, v1_acc, u_avg_acc)
    eq_lbc_acc_alloc = .true.
end if

!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny
do i = 1, ld
    u1_acc(i,j) = u(i,j,1)
    v1_acc(i,j) = v(i,j,1)
end do
end do

call test_filter_plane_gpu(u1_acc)
call test_filter_plane_gpu(v1_acc)

!$acc parallel loop collapse(2) default(present) private(u_avg_pt) async(1)
do j = 1, ny
do i = 1, nx
    u_avg_pt = sqrt(u1_acc(i,j)**2 + v1_acc(i,j)**2)
    u_avg_acc(i,j) = u_avg_pt
#ifndef PPSCALARS
    ustar_lbc(i,j) = u_avg_pt*vonk/log(0.5_rprec*dz/zo)
#endif
end do
end do

#ifdef PPSCALARS
call obukhov(u_avg_acc)
#endif

!$acc parallel loop collapse(2) default(present) private(u_avg_pt, const) async(1)
do j = 1, ny
do i = 1, nx
    u_avg_pt = u_avg_acc(i,j)
    const = -(ustar_lbc(i,j)**2)/u_avg_pt
    txz(i,j,1) = const*u1_acc(i,j)
    tyz(i,j,1) = const*v1_acc(i,j)
#ifdef PPSCALARS
    dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg_pt        &
        * phi_m(i,j)
    dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg_pt        &
        * phi_m(i,j)
#else
    dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg_pt
    dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg_pt
#endif
    dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
    dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
end do
end do
return
#endif


u1 = u(:,:,1)
v1 = v(:,:,1)
call test_filter(u1)
call test_filter(v1)
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
#ifdef PPSCALARS
call obukhov(u_avg)
#else
ustar_lbc = u_avg*vonk/denom
#endif

do j = 1, ny
    do i = 1, nx
        const = -(ustar_lbc(i,j)**2)/u_avg(i,j)
        txz(i,j,1) = const*u1(i,j)
        tyz(i,j,1) = const*v1(i,j)
        !this is as in Moeng 84
#ifdef PPSCALARS
        dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)   &
            * phi_m(i,j)
        dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)   &
            * phi_m(i,j)
#else
        dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)
        dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)
#endif
        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_lbc

!*******************************************************************************
subroutine ws_equilibrium_ubc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v
use test_filtermodule, only : test_filter
implicit none
integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg, ustar
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const
#if defined(PPSGS_GPU)
logical, save :: eq_ubc_acc_alloc = .false.
real(rprec), save, allocatable, dimension(:,:) :: u1t_acc, v1t_acc
real(rprec) :: u_avg_pt, ustar_pt
#endif

#if defined(PPSGS_GPU)
! Device-resident equilibrium wall model for the top boundary.
if (.not. eq_ubc_acc_alloc) then
    allocate(u1t_acc(ld,ny), v1t_acc(ld,ny))
    !$acc enter data create(u1t_acc, v1t_acc)
    eq_ubc_acc_alloc = .true.
end if

!$acc parallel loop collapse(2) default(present) async(1)
do j = 1, ny
do i = 1, ld
    u1t_acc(i,j) = u(i,j,nz-1)
    v1t_acc(i,j) = v(i,j,nz-1)
end do
end do

call test_filter_plane_gpu(u1t_acc)
call test_filter_plane_gpu(v1t_acc)

!$acc parallel loop collapse(2) default(present) private(u_avg_pt, ustar_pt, const) async(1)
do j = 1, ny
do i = 1, nx
    u_avg_pt = sqrt(u1t_acc(i,j)**2 + v1t_acc(i,j)**2)
    ustar_pt = u_avg_pt*vonk/log(0.5_rprec*dz/zo)
    const = (ustar_pt**2)/u_avg_pt ! diff sign for upper b.c.
    txz(i,j,nz) = const*u1t_acc(i,j)
    tyz(i,j,nz) = const*v1t_acc(i,j)
    dudz(i,j,nz) = -ustar_pt/(0.5_rprec*dz*vonK)*u(i,j,nz-1)/u_avg_pt
    dvdz(i,j,nz) = -ustar_pt/(0.5_rprec*dz*vonK)*v(i,j,nz-1)/u_avg_pt
    dudz(i,j,nz) = merge(0._rprec,dudz(i,j,nz),u(i,j,nz-1).eq.0._rprec)
    dvdz(i,j,nz) = merge(0._rprec,dvdz(i,j,nz),v(i,j,nz-1).eq.0._rprec)
end do
end do
return
#endif


u1 = u(:,:,nz-1)
v1 = v(:,:,nz-1)
call test_filter(u1)
call test_filter(v1)
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar = u_avg*vonk/denom

do j = 1, ny
    do i = 1, nx
        const = (ustar(i,j)**2)/u_avg(i,j) ! diff sign for upper b.c.
        txz(i,j,nz) = const*u1(i,j)
        tyz(i,j,nz) = const*v1(i,j)
        !this is as in Moeng 84
        dudz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*u(i,j,nz-1)/u_avg(i,j)
        dvdz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*v(i,j,nz-1)/u_avg(i,j)
        dudz(i,j,nz) = merge(0._rprec,dudz(i,j,nz),u(i,j,nz-1).eq.0._rprec)
        dvdz(i,j,nz) = merge(0._rprec,dvdz(i,j,nz),v(i,j,nz-1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_ubc

end subroutine wallstress
