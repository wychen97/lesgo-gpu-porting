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

!*******************************************************************************
module sim_param
!*******************************************************************************
use types, only : rprec
use param, only : ld, nx, ny, nz, lbz, u_star
implicit none

save
public

logical :: sim_param_initialized = .false.

#if defined(ENABLE_CUDA) && !defined(PPLES_GPU)
  ! CUDA managed-memory path: arrays are declared managed for NVHPC/CUDA builds.
  real(rprec), managed, allocatable, dimension(:,:,:) :: u, v, w
  real(rprec), managed, allocatable, dimension(:,:,:) :: dudx, dudy, dudz
  real(rprec), managed, allocatable, dimension(:,:,:) :: dvdx, dvdy, dvdz
  real(rprec), managed, allocatable, dimension(:,:,:) :: dwdx, dwdy, dwdz
  real(rprec), managed, allocatable, dimension(:,:,:) :: RHSx, RHSy, RHSz
  real(rprec), managed, allocatable, dimension(:,:,:) :: RHSx_f, RHSy_f, RHSz_f
  real(rprec), managed, allocatable, dimension(:,:,:) :: dpdx, dpdy, dpdz
  real(rprec), managed, allocatable, dimension(:,:,:) :: txx, txy, tyy, txz, tyz, tzz
  real(rprec), managed, allocatable, dimension(:,:,:) :: divtx, divty, divtz
  real(rprec), managed, allocatable, dimension(:,:,:) :: fx, fy, fz, fxa, fya, fza

  real(rprec), managed, allocatable, target, dimension(:,:,:) :: p
  real(rprec), managed, allocatable, dimension(:,:) :: ustar_lbc
#else
  ! CPU and explicit-residency GPU builds use ordinary host allocatables.
  real(rprec), allocatable, dimension(:,:,:) :: u, v, w
  real(rprec), allocatable, dimension(:,:,:) :: dudx, dudy, dudz
  real(rprec), allocatable, dimension(:,:,:) :: dvdx, dvdy, dvdz
  real(rprec), allocatable, dimension(:,:,:) :: dwdx, dwdy, dwdz
  real(rprec), allocatable, dimension(:,:,:) :: RHSx, RHSy, RHSz
  real(rprec), allocatable, dimension(:,:,:) :: RHSx_f, RHSy_f, RHSz_f
  real(rprec), allocatable, dimension(:,:,:) :: dpdx, dpdy, dpdz
  real(rprec), allocatable, dimension(:,:,:) :: txx, txy, tyy, txz, tyz, tzz
  real(rprec), allocatable, dimension(:,:,:) :: divtx, divty, divtz
  real(rprec), allocatable, dimension(:,:,:) :: fx, fy, fz, fxa, fya, fza

  real(rprec), allocatable, target, dimension(:,:,:) :: p
  real(rprec), allocatable, dimension(:,:) :: ustar_lbc
#endif

#ifdef PPLES_GPU
! Experimental explicit-residency LES route. These arrays are ordinary host
! allocatables with persistent OpenACC device mirrors. Keep host/device movement
! explicit at call-site boundaries instead of relying on CUDA managed migration.
!$acc declare create(u, v, w, dudx, dudy, dudz, dvdx, dvdy, dvdz)
!$acc declare create(dwdx, dwdy, dwdz)
!$acc declare create(RHSx, RHSy, RHSz, RHSx_f, RHSy_f, RHSz_f)
!$acc declare create(dpdx, dpdy, dpdz)
!$acc declare create(txx, txy, tyy, txz, tyz, tzz)
!$acc declare create(divtx, divty, divtz)
!$acc declare create(p)
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
!$acc declare create(fxa, fya, fza)
#endif
#if defined(PPLVLSET) || defined(PPATM)
!$acc declare create(fx, fy, fz)
#endif
!$acc declare create(ustar_lbc)
#endif

contains

!*******************************************************************************
subroutine sim_param_init ()
!*******************************************************************************
!
! This subroutine initilizes all global arrays defined in the sim_param
! module. Here they are allocated and initialized to zero.
!
implicit none

allocate ( u(ld, ny, lbz:nz) ); u = 0.0_rprec
allocate ( v(ld, ny, lbz:nz) ); v = 0.0_rprec
allocate ( w(ld, ny, lbz:nz) ); w = 0.0_rprec
allocate( dudx(ld, ny, lbz:nz) ); dudx = 0.0_rprec
allocate( dudy(ld, ny, lbz:nz) ); dudy = 0.0_rprec
allocate( dudz(ld, ny, lbz:nz) ); dudz = 0.0_rprec
allocate( dvdx(ld, ny, lbz:nz) ); dvdx = 0.0_rprec
allocate( dvdy(ld, ny, lbz:nz) ); dvdy = 0.0_rprec
allocate( dvdz(ld, ny, lbz:nz) ); dvdz = 0.0_rprec
allocate( dwdx(ld, ny, lbz:nz) ); dwdx = 0.0_rprec
allocate( dwdy(ld, ny, lbz:nz) ); dwdy = 0.0_rprec
allocate( dwdz(ld, ny, lbz:nz) ); dwdz = 0.0_rprec
allocate( RHSx(ld, ny, lbz:nz) ); RHSx = 0.0_rprec
allocate( RHSy(ld, ny, lbz:nz) ); RHSy = 0.0_rprec
allocate( RHSz(ld, ny, lbz:nz) ); RHSz = 0.0_rprec
allocate( RHSx_f(ld, ny, lbz:nz) ); RHSx_f = 0.0_rprec
allocate( RHSy_f(ld, ny, lbz:nz) ); RHSy_f = 0.0_rprec
allocate( RHSz_f(ld, ny, lbz:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
allocate ( txx(ld, ny, lbz:nz) ); txx = 0.0_rprec
allocate ( txy(ld, ny, lbz:nz) ); txy = 0.0_rprec
allocate ( tyy(ld, ny, lbz:nz) ); tyy = 0.0_rprec
allocate ( txz(ld, ny, lbz:nz) ); txz = 0.0_rprec
allocate ( tyz(ld, ny, lbz:nz) ); tyz = 0.0_rprec
allocate ( tzz(ld, ny, lbz:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
allocate ( divtx(ld, ny, lbz:nz) ); divtx = 0.0_rprec
allocate ( divty(ld, ny, lbz:nz) ); divty = 0.0_rprec
allocate ( divtz(ld, ny, lbz:nz) ); divtz = 0.0_rprec

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
allocate ( fxa(ld, ny, lbz:nz) ); fxa = 0.0_rprec
allocate ( fya(ld, ny, lbz:nz) ); fya = 0.0_rprec
allocate ( fza(ld, ny, lbz:nz) ); fza = 0.0_rprec
#endif

#if defined(PPLVLSET) || defined(PPATM)
allocate ( fx(ld, ny, nz) ); fx = 0.0_rprec
allocate ( fy(ld, ny, nz) ); fy = 0.0_rprec
allocate ( fz(ld, ny, nz) ); fz = 0.0_rprec
#endif

allocate( ustar_lbc(nx, ny) ); ustar_lbc = u_star

#ifdef PPLES_GPU
! Initialize device mirrors for arrays whose host values are known here. The
! velocity field itself is pushed after initial-condition setup in initialize.f90.
!$acc update device(dudx, dudy, dudz, dvdx, dvdy, dvdz)
!$acc update device(dwdx, dwdy, dwdz)
!$acc update device(RHSx, RHSy, RHSz, RHSx_f, RHSy_f, RHSz_f)
!$acc update device(dpdx, dpdy, dpdz)
!$acc update device(txx, txy, tyy, txz, tyz, tzz)
!$acc update device(divtx, divty, divtz)
!$acc update device(p)
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
!$acc update device(fxa, fya, fza)
#endif
#if defined(PPLVLSET) || defined(PPATM)
!$acc update device(fx, fy, fz)
#endif
!$acc update device(ustar_lbc)
#endif

sim_param_initialized = .true.

end subroutine sim_param_init

end module sim_param
