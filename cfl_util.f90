!!
!!  Copyright (C) 2010-2017  Johns Hopkins University
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
module cfl_util
!*******************************************************************************
!
! This module provides the subroutines/functions for getting CFL related
! quantities
!

implicit none
save
private

public get_max_cfl, get_cfl_dt

contains


!*******************************************************************************
function get_max_cfl() result(cfl)
!*******************************************************************************
!
! This function provides the value of the maximum CFL in the entire
! domain
!
use types, only : rprec
use param, only : dt, dx, dy, dz, nx, ny, nz
use sim_param, only : u,v,w

#ifdef PPMPI
use mpi, only : MPI_COMM_WORLD, MPI_MAX, mpi_allreduce
use param, only : ierr, MPI_RPREC
#endif

implicit none
real(rprec) :: cfl

real(rprec) :: cfl_u, cfl_v, cfl_w

#ifdef PPMPI
real(rprec) :: cfl_buf
#endif
integer :: jx, jy, jz

#if defined(PPLES_GPU)
! Explicit-residency: reduce over the DEVICE velocity (host copy is stale here).
cfl_u = 0._rprec
cfl_v = 0._rprec
cfl_w = 0._rprec
!$acc wait(1)
!$acc parallel loop collapse(3) default(present) reduction(max:cfl_u,cfl_v,cfl_w)
do jz = 1, nz - 1
do jy = 1, ny
do jx = 1, nx
    cfl_u = max(cfl_u, abs(u(jx,jy,jz)))
    cfl_v = max(cfl_v, abs(v(jx,jy,jz)))
    cfl_w = max(cfl_w, abs(w(jx,jy,jz)))
end do
end do
end do
cfl_u = cfl_u / dx
cfl_v = cfl_v / dy
cfl_w = cfl_w / dz
#else
cfl_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
cfl_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
cfl_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz
#endif

cfl = dt * maxval( (/ cfl_u, cfl_v, cfl_w /) )

#ifdef PPMPI
call mpi_allreduce(cfl, cfl_buf, 1, MPI_RPREC, MPI_MAX, MPI_COMM_WORLD, ierr)
cfl = cfl_buf
#endif

end function get_max_cfl

!*******************************************************************************
function get_cfl_dt() result(dt)
!*******************************************************************************
!
! This functions determines the maximum allowable time step based on the CFL
! value specified in the param module
!
use types, only : rprec
use param, only : cfl, dx, dy, dz, nx, ny, nz
use sim_param, only : u,v,w

#ifdef PPMPI
use mpi, only : MPI_COMM_WORLD, MPI_MIN, mpi_allreduce
use param, only : ierr, MPI_RPREC
#endif

implicit none

real(rprec) :: dt

! dt inverse
real(rprec) :: dt_inv_u, dt_inv_v, dt_inv_w

#ifdef PPMPI
real(rprec) :: dt_buf
#endif
integer :: jx, jy, jz

! Avoid division by computing max dt^-1
#if defined(PPLES_GPU)
! Explicit-residency (mem:separate): reduce over the DEVICE velocity. The host
! copy of u,v,w is stale here (only synced every nenergy / at the ATM forcing),
! so the host maxval below would compute dt from a lagged field -> CFL blow-up.
! Mirrors the reference get_cfl_dt PPLES_GPU path.
dt_inv_u = 0._rprec
dt_inv_v = 0._rprec
dt_inv_w = 0._rprec
!$acc wait(1)
!$acc parallel loop collapse(3) default(present) reduction(max:dt_inv_u,dt_inv_v,dt_inv_w)
do jz = 1, nz - 1
do jy = 1, ny
do jx = 1, nx
    dt_inv_u = max(dt_inv_u, abs(u(jx,jy,jz)))
    dt_inv_v = max(dt_inv_v, abs(v(jx,jy,jz)))
    dt_inv_w = max(dt_inv_w, abs(w(jx,jy,jz)))
end do
end do
end do
dt_inv_u = dt_inv_u / dx
dt_inv_v = dt_inv_v / dy
dt_inv_w = dt_inv_w / dz
#else
dt_inv_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
dt_inv_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
dt_inv_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz
#endif

dt = cfl / maxval( (/ dt_inv_u, dt_inv_v, dt_inv_w /) )

#ifdef PPMPI
call mpi_allreduce(dt, dt_buf, 1, MPI_RPREC, MPI_MIN, MPI_COMM_WORLD, ierr)
dt = dt_buf
#endif

end function get_cfl_dt

end module cfl_util
