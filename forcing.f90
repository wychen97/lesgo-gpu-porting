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
module forcing
!*******************************************************************************
!
! Provides subroutines and functions for computing forcing terms on the
! velocity field. Provides driver routine for IBM forces
! (forcing_induced), driver routine for RNS, turbine, etc. forcing
! (forcing_applied), and for the projection step. Also included are
! routines for enforcing a uniform inflow and the fringe region
! treatment.
!
! Navigation map:
!   - environment and diagnostics: forcing_env_* helpers and project reports
!   - applied forcing: forcing_random and forcing_applied
!   - induced forcing and projection: forcing_induced and project
!   - velocity halo paths: project_sync_velocity_* and pack/unpack helpers
!   - optional bridges: HIT inflow, turbines, ATM, and LVLSET guarded paths
!
! Keep force application and projection ordering aligned with main.f90.  Moving
! work across that boundary changes timestep dependencies and validation scope.

#ifdef PPHIT
use hit_inflow, only : inflow_HIT
#endif

#ifdef ENABLE_CUDA
use cudafor
#endif
use types, only : rprec

implicit none

save

private

public :: forcing_random, forcing_applied, forcing_induced, project

#ifdef ENABLE_CUDA
real(rprec), device, allocatable, save :: project_halo_send_down(:)
real(rprec), device, allocatable, save :: project_halo_send_up(:)
real(rprec), device, allocatable, save :: project_halo_recv_down(:)
real(rprec), device, allocatable, save :: project_halo_recv_up(:)
integer, save :: project_stage_count = 0
#endif

contains

!*******************************************************************************
logical function forcing_env_true_token_enabled(name)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
character(len=32) :: setting
integer :: stat

forcing_env_true_token_enabled = .false.
call get_environment_variable(name, setting, status=stat)
if (stat == 0) then
    select case (trim(adjustl(setting)))
    case ('1', 'true', 'TRUE', 'True', 'on', 'ON', 'On', 'yes', 'YES', 'Yes')
        forcing_env_true_token_enabled = .true.
    case default
        forcing_env_true_token_enabled = .false.
    end select
end if

end function forcing_env_true_token_enabled

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function project_cuda_enabled()
!*******************************************************************************
implicit none

project_cuda_enabled = .true.

end function project_cuda_enabled

!*******************************************************************************
logical function project_packed_halo_enabled()
!*******************************************************************************
implicit none

project_packed_halo_enabled = .true.

end function project_packed_halo_enabled

!*******************************************************************************
logical function project_direct_halo_enabled()
!*******************************************************************************
implicit none

project_direct_halo_enabled = .false.

end function project_direct_halo_enabled

!*******************************************************************************
logical function project_overlap_boundary_enabled()
!*******************************************************************************
implicit none

project_overlap_boundary_enabled = .true.

end function project_overlap_boundary_enabled

!*******************************************************************************
logical function project_stage_timing_enabled()
!*******************************************************************************
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = forcing_env_true_token_enabled('LESGO_PROJECT_STAGE_TIMING')
    initialized = .true.
end if

project_stage_timing_enabled = enabled

end function project_stage_timing_enabled


!*******************************************************************************
logical function applied_force_reset_cuda_enabled()
!*******************************************************************************
implicit none

applied_force_reset_cuda_enabled = .true.

end function applied_force_reset_cuda_enabled

!*******************************************************************************
logical function random_force_cuda_enabled()
!*******************************************************************************
implicit none

random_force_cuda_enabled = .true.

end function random_force_cuda_enabled

!*******************************************************************************
subroutine forcing_cuda_sync(where)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'forcing CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'forcing CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine forcing_cuda_sync
#endif

#if defined(PPLVLSET) && defined(PPLES_GPU) && !defined(ENABLE_CUDA)
!*******************************************************************************
logical function lvlset_bridge_timing_enabled()
!*******************************************************************************
implicit none

logical, save :: initialized = .false.
logical, save :: enabled = .false.

if (.not. initialized) then
    enabled = forcing_env_true_token_enabled('LESGO_LVLSET_BRIDGE_TIMING')
    initialized = .true.
end if

lvlset_bridge_timing_enabled = enabled

end function lvlset_bridge_timing_enabled

!*******************************************************************************
subroutine lvlset_bridge_time(tnow)
!*******************************************************************************
#ifdef PPMPI
use mpi
#endif
implicit none

real(rprec), intent(out) :: tnow

#ifdef PPMPI
tnow = mpi_wtime()
#else
call cpu_time(tnow)
#endif

end subroutine lvlset_bridge_time

!*******************************************************************************
subroutine lvlset_bridge_report(call_count, wait_sync, update_self,            &
    host_force, update_device)
!*******************************************************************************
use param, only : coord
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC
use mpi
#endif
implicit none

integer, intent(in) :: call_count
real(rprec), intent(in) :: wait_sync, update_self, host_force, update_device
real(rprec) :: wait_max, self_max, force_max, device_max, total

#ifdef PPMPI
call mpi_allreduce(wait_sync, wait_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(update_self, self_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(host_force, force_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(update_device, device_max, 1, MPI_RPREC, MPI_MAX, comm,     &
    ierr)
#else
wait_max = wait_sync
self_max = update_self
force_max = host_force
device_max = update_device
#endif

total = wait_max + self_max + force_max + device_max

if (coord == 0) then
    write(*,'(a,i8)') 'LVLSET bridge timing (max rank), call ', call_count
    write(*,'(1a,E15.7)') '  Pre-bridge acc wait: ', wait_max
    write(*,'(1a,E15.7)') '  Device-to-host refresh: ', self_max
    write(*,'(1a,E15.7)') '  Host level_set_forcing: ', force_max
    write(*,'(1a,E15.7)') '  Host-to-device restore: ', device_max
    write(*,'(1a,E15.7)') '  Bridge sum: ', total
end if

end subroutine lvlset_bridge_report
#endif

!*******************************************************************************
subroutine forcing_random()
!*******************************************************************************
!
! This subroutine generates a random body force that is helpful to
! trigger transition at low Re DNS. The forces are applied to RHS in
! evaluation of u* (not at wall) so that mass conservation is preserved.
!
use types, only : rprec
use param, only : nx,ny,nz,rms_random_force
use sim_param, only : RHSy, RHSz

real(rprec) :: dummy_rand
integer :: jx,jy,jz
#if defined(ENABLE_CUDA) || (defined(PPLES_GPU) && !defined(ENABLE_CUDA))
    logical, save :: random_force_allocated = .false.
#ifdef ENABLE_CUDA
    real(rprec), managed, save, allocatable, dimension(:,:,:) :: rand_y, rand_z
#endif
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
    real(rprec), save, allocatable, dimension(:,:,:) :: rand_y, rand_z
#endif
#endif

    ! Note: the "default" rms of a unif variable is 0.289
call init_random_seed
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
if (.not. random_force_allocated) then
    allocate(rand_y(nx,ny,2:nz-1), rand_z(nx,ny,2:nz-1))
    !$acc enter data create(rand_y, rand_z)
    random_force_allocated = .true.
end if

do jz = 2, nz-1
do jy = 1, ny
do jx = 1, nx
    call random_number(rand_y(jx,jy,jz))
    call random_number(rand_z(jx,jy,jz))
end do
end do
end do

! Keep the intrinsic random sequence on the host, then apply the resulting
! forcing on the resident RHS fields.
!$acc update device(rand_y, rand_z) async(1)
!$acc parallel loop collapse(3) default(present) async(1)
do jz = 2, nz-1
do jy = 1, ny
do jx = 1, nx
    RHSy(jx,jy,jz) = RHSy(jx,jy,jz)                                              &
        + (rms_random_force/.289_rprec)*(rand_y(jx,jy,jz)-.5_rprec)
    RHSz(jx,jy,jz) = RHSz(jx,jy,jz)                                              &
        + (rms_random_force/.289_rprec)*(rand_z(jx,jy,jz)-.5_rprec)
end do
end do
end do
return
#endif
#ifdef ENABLE_CUDA
if (random_force_cuda_enabled()) then
    if (.not. random_force_allocated) then
        allocate(rand_y(nx,ny,2:nz-1), rand_z(nx,ny,2:nz-1))
        random_force_allocated = .true.
    end if

    do jz = 2, nz-1
    do jy = 1, ny
    do jx = 1, nx
        call random_number(rand_y(jx,jy,jz))
        call random_number(rand_z(jx,jy,jz))
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do jz = 2, nz-1
    do jy = 1, ny
    do jx = 1, nx
        RHSy(jx,jy,jz) = RHSy(jx,jy,jz)                                      &
            + (rms_random_force/.289_rprec)*(rand_y(jx,jy,jz)-.5_rprec)
        RHSz(jx,jy,jz) = RHSz(jx,jy,jz)                                      &
            + (rms_random_force/.289_rprec)*(rand_z(jx,jy,jz)-.5_rprec)
    end do
    end do
    end do
    call forcing_cuda_sync('forcing_random')
    return
end if
#endif
do jz = 2, nz-1 ! don't force too close to the wall
do jy = 1, ny
do jx = 1, nx
    call random_number(dummy_rand)
    RHSy(jx,jy,jz) = RHSy(jx,jy,jz) +                                          &
        (rms_random_force/.289_rprec)*(dummy_rand-.5_rprec)
    call random_number(dummy_rand)
    RHSz(jx,jy,jz) = RHSz(jx,jy,jz) +                                          &
        (rms_random_force/.289_rprec)*(dummy_rand-.5_rprec)
end do
end do
end do

end subroutine forcing_random

!*******************************************************************************
subroutine forcing_applied()
!*******************************************************************************
!
!  This subroutine acts as a driver for applying pointwise body forces
!  into the domain. Subroutines contained here should modify f{x,y,z}a
!  which are explicitly applied forces. These forces are applied to RHS
!  in the evaluation of u* so that mass conservation is preserved.
!
use types, only : rprec
use param, only : nx, ny, nz, lbz

#if defined(PPTURBINES) && !defined(PPATM)
#ifdef PPLES_GPU
use sim_param, only : u, v, w, fxa, fya, fza
#else
use sim_param, only : fxa, fya, fza
#endif
use turbines, only:turbines_forcing, turbines_acc_available, turbines_forcing_acc
#endif

#ifdef PPATM
use sim_param, only : fxa, fya, fza ! The body force components
use atm_lesgo_interface, only : atm_lesgo_forcing
#endif

implicit none
integer :: jx, jy, jz

#if defined(PPTURBINES) && !defined(PPATM)
! Reset applied force arrays for the legacy drag-disk turbine model.
! If PPATM is also compiled, ATM is the active turbine model for this
! executable and this legacy path is intentionally skipped.
#ifdef PPLES_GPU
if (turbines_acc_available()) then
    call turbines_forcing_acc()
else
    ! Fallback for dynamic or multi-rank legacy turbines. Preserve correctness
    ! under the explicit-residency LES path by bridging velocity and force
    ! arrays at the module boundary.
    !$acc wait(1)
    !$acc update self(u, v, w)
    fxa = 0._rprec
    fya = 0._rprec
    fza = 0._rprec
    call turbines_forcing ()
    !$acc update device(fxa, fya, fza)
end if
#elif defined(ENABLE_CUDA)
if (applied_force_reset_cuda_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = lbz, nz
    do jy = 1, ny
    do jx = 1, nx
        fxa(jx,jy,jz) = 0._rprec
        fya(jx,jy,jz) = 0._rprec
        fza(jx,jy,jz) = 0._rprec
    end do
    end do
    end do
    call forcing_cuda_sync('turbine applied-force reset')
else
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
end if
call turbines_forcing ()
#else
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
call turbines_forcing ()
#endif
#endif


#ifdef PPATM
#ifdef PPLES_GPU
! Explicit-residency, fully device-resident ATM: velocity sampling reads device
! u,v,w_uv (atm_sample_velocity_atpoint_gpu) and the convolution scatters forces
! straight into the device fxa/fya/fza, so NO host<->device velocity/force
! transfer is needed here -- just zero the force arrays on the device.
!$acc parallel loop collapse(3) default(present)
do jz = lbz, nz
do jy = 1, ny
do jx = 1, nx
    fxa(jx,jy,jz) = 0._rprec
    fya(jx,jy,jz) = 0._rprec
    fza(jx,jy,jz) = 0._rprec
end do
end do
end do
#elif defined(ENABLE_CUDA)
if (applied_force_reset_cuda_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = lbz, nz
    do jy = 1, ny
    do jx = 1, nx
        fxa(jx,jy,jz) = 0._rprec
        fya(jx,jy,jz) = 0._rprec
        fza(jx,jy,jz) = 0._rprec
    end do
    end do
    end do
    call forcing_cuda_sync('ATM applied-force reset')
else
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
end if
#else
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
#endif
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
! Phase 1 (sampling + blade force model + gather) already ran earlier in the
! step from main.f90, overlapped with the SGS/convection GPU kernels. Here we
! only convolute/apply the gathered forces.
call atm_lesgo_forcing (phase=2)
#else
call atm_lesgo_forcing ()
#endif
! Under PPLES_GPU, fxa/fya/fza were scattered into on the device by
! atm_convolute_atpoint_gpu, so no host->device force transfer is needed here.
#endif

end subroutine forcing_applied

!*******************************************************************************
subroutine forcing_induced()
!*******************************************************************************
!
!  These forces are designated as induced forces such that they are
!  chosen to obtain a desired velocity at time
!  step m+1. If this is not the case, care should be taken so that the forces
!  here are divergence free in order to preserve mass conservation. For
!  non-induced forces such as explicitly applied forces they should be
!  placed in forcing_applied.
!
use types, only : rprec
#ifdef PPLVLSET
use level_set, only : level_set_forcing
use sim_param, only : u, v, w, dpdx, dpdy, dpdz, fx, fy, fz
#endif
implicit none

#if defined(PPLVLSET) && defined(PPLES_GPU) && !defined(ENABLE_CUDA)
logical :: lvlset_timing
real(rprec) :: lvlset_t0, lvlset_t1, lvlset_t2, lvlset_t3, lvlset_t4
integer, save :: lvlset_bridge_count = 0
#endif

#ifdef PPLVLSET
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
lvlset_timing = lvlset_bridge_timing_enabled()
if (lvlset_timing) then
    lvlset_bridge_count = lvlset_bridge_count + 1
    call lvlset_bridge_time(lvlset_t0)
end if
! LVLSET forcing is still host code in the OpenACC LES build. Refresh the
! host inputs, then push the induced forces back for the GPU projection step.
!$acc wait(1)
if (lvlset_timing) call lvlset_bridge_time(lvlset_t1)
!$acc update self(u, v, w, dpdx, dpdy, dpdz)
if (lvlset_timing) call lvlset_bridge_time(lvlset_t2)
#endif
! Initialize
fx = 0._rprec
fy = 0._rprec
fz = 0._rprec
!  Compute the level set IBM forces
call level_set_forcing ()
#if defined(PPLES_GPU) && !defined(ENABLE_CUDA)
if (lvlset_timing) call lvlset_bridge_time(lvlset_t3)
!$acc update device(u, v, w, fx, fy, fz)
if (lvlset_timing) then
    call lvlset_bridge_time(lvlset_t4)
    call lvlset_bridge_report(lvlset_bridge_count, lvlset_t1 - lvlset_t0,      &
        lvlset_t2 - lvlset_t1, lvlset_t3 - lvlset_t2, lvlset_t4 - lvlset_t3)
end if
#endif
#endif

end subroutine forcing_induced

!*******************************************************************************
subroutine project ()
!*******************************************************************************
!
! provides u, v, w at 1:nz
!
use param
use sim_param
use messages
#ifdef ENABLE_CUDA
use inflow, only : apply_inflow, inflow_cuda_enabled
#else
use inflow, only : apply_inflow
#endif
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
use mpi

#endif
implicit none

integer :: jx, jy, jz
integer :: jz_min
real(rprec) :: RHS, tconst
#ifdef ENABLE_CUDA
logical :: project_boundary_done
logical :: project_stage_enabled
logical :: project_halo_already_synced
real(rprec) :: project_t0, project_t1
real(rprec) :: project_stage_update, project_stage_halo, project_stage_boundary
#endif

! Caching
tconst = tadv1 * dt
#ifdef PPLES_GPU
! Experimental explicit-residency projection path.  Keep this separate from the
! legacy CUDA Fortran path so PPLES_GPU arrays stay OpenACC-present instead of
! relying on CUDA managed dummy arguments.
if (coord == 0) then
    jz_min = 2
else
    jz_min = 1
end if

!$acc parallel loop collapse(3) default(present)
do jz = 1, nz - 1
do jy = 1, ny
do jx = 1, nx
#ifdef PPLVLSET
    u(jx, jy, jz) = u(jx, jy, jz) - tconst * dpdx(jx, jy, jz)                 &
        + dt * fx(jx, jy, jz)
    v(jx, jy, jz) = v(jx, jy, jz) - tconst * dpdy(jx, jy, jz)                 &
        + dt * fy(jx, jy, jz)
    if (jz >= jz_min) then
    w(jx, jy, jz) = w(jx, jy, jz) - tconst * dpdz(jx, jy, jz)                 &
        + dt * fz(jx, jy, jz)
    end if
#else
    u(jx, jy, jz) = u(jx, jy, jz) - tconst * dpdx(jx, jy, jz)
    v(jx, jy, jz) = v(jx, jy, jz) - tconst * dpdy(jx, jy, jz)
    if (jz >= jz_min) then
    w(jx, jy, jz) = w(jx, jy, jz) - tconst * dpdz(jx, jy, jz)
    end if
#endif
end do
end do
end do

call apply_inflow()

#ifdef PPMPI
if (nproc > 1) then
#ifdef PPGPU_AWARE_MPI
    ! GPU-aware halo exchange: device pointers straight to MPICH. Replicates
    ! mpi_sync_real_array(var, 0, MPI_SYNC_DOWNUP) exactly:
    !   sync_down (tag 1): send k=1 down, recv k=nz from up
    !   sync_up   (tag 2): send k=nz-1 up, recv k=0 from down
    !$acc wait(1)
    !$acc host_data use_device(u, v, w)
    call mpi_sendrecv (u(1,1,1),    ld*ny, MPI_RPREC, down, 1,                 &
                       u(1,1,nz),   ld*ny, MPI_RPREC, up,   1, comm, status, ierr)
    call mpi_sendrecv (v(1,1,1),    ld*ny, MPI_RPREC, down, 1,                 &
                       v(1,1,nz),   ld*ny, MPI_RPREC, up,   1, comm, status, ierr)
    call mpi_sendrecv (w(1,1,1),    ld*ny, MPI_RPREC, down, 1,                 &
                       w(1,1,nz),   ld*ny, MPI_RPREC, up,   1, comm, status, ierr)
    call mpi_sendrecv (u(1,1,nz-1), ld*ny, MPI_RPREC, up,   2,                 &
                       u(1,1,0),    ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv (v(1,1,nz-1), ld*ny, MPI_RPREC, up,   2,                 &
                       v(1,1,0),    ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
    call mpi_sendrecv (w(1,1,nz-1), ld*ny, MPI_RPREC, up,   2,                 &
                       w(1,1,0),    ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
    !$acc end host_data
#else
    ! Correctness-first halo staging.  Replace with host_data/use_device after
    ! the explicit-residency route passes multi-rank validation.
    !$acc update self(u(:,:,1), v(:,:,1), w(:,:,1))
    !$acc update self(u(:,:,nz-1), v(:,:,nz-1), w(:,:,nz-1))
    call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP )
    !$acc update device(u(:,:,0), v(:,:,0), w(:,:,0))
    !$acc update device(u(:,:,nz), v(:,:,nz), w(:,:,nz))
#endif
end if
#endif

#ifdef PPMPI
if (coord == nproc-1) then
#endif
    if (ubc_mom == 0) then
        !$acc parallel loop collapse(2) default(present)
        do jy = 1, ny
        do jx = 1, nx
            u(jx,jy,nz) = u(jx,jy,nz-1)
            v(jx,jy,nz) = v(jx,jy,nz-1)
            w(jx,jy,nz) = 0._rprec
        end do
        end do
    else
        !$acc parallel loop collapse(2) default(present)
        do jy = 1, ny
        do jx = 1, nx
            w(jx,jy,nz) = 0._rprec
        end do
        end do
    end if
#ifdef PPMPI
endif
#endif

if (coord == 0) then
    !$acc parallel loop collapse(2) default(present)
    do jy = 1, ny
    do jx = 1, nx
        w(jx,jy,1) = 0._rprec
    end do
    end do
end if

return
#endif
#ifdef ENABLE_CUDA
project_boundary_done = .false.
project_stage_update = 0._rprec
project_stage_halo = 0._rprec
project_stage_boundary = 0._rprec
project_stage_enabled = project_stage_timing_enabled()
project_halo_already_synced = .false.
if (project_stage_enabled) then
    project_stage_count = project_stage_count + 1
#ifdef PPMPI
    project_t0 = mpi_wtime()
#else
    call cpu_time(project_t0)
#endif
end if
#endif

#ifdef ENABLE_CUDA
if (project_cuda_enabled()) then
    !$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, nx
#ifdef PPLVLSET
        u(jx, jy, jz) = u(jx, jy, jz) - tconst * dpdx(jx, jy, jz)             &
            + dt * fx(jx, jy, jz)
        v(jx, jy, jz) = v(jx, jy, jz) - tconst * dpdy(jx, jy, jz)             &
            + dt * fy(jx, jy, jz)
#else
        u(jx, jy, jz) = u(jx, jy, jz) - tconst * dpdx(jx, jy, jz)
        v(jx, jy, jz) = v(jx, jy, jz) - tconst * dpdy(jx, jy, jz)
#endif
    end do
    end do
    end do

    if (coord == 0) then
        jz_min = 2
    else
        jz_min = 1
    end if

    !$cuf kernel do(3) <<<*,*>>>
    do jz = jz_min, nz - 1
    do jy = 1, ny
    do jx = 1, nx
#ifdef PPLVLSET
        w(jx, jy, jz) = w(jx, jy, jz) - tconst * dpdz(jx, jy, jz)             &
            + dt * fz(jx, jy, jz)
#else
        w(jx, jy, jz) = w(jx, jy, jz) - tconst * dpdz(jx, jy, jz)
#endif
    end do
    end do
    end do
    if (.not. (inflow_type == 1 .and. inflow_cuda_enabled())) then
        call forcing_cuda_sync('pressure-gradient update')
    end if
else
#endif
do jz = 1, nz - 1
do jy = 1, ny
do jx = 1, nx
#ifdef PPLVLSET
    RHS = -tadv1 * dpdx(jx, jy, jz)
    u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS + fx(jx, jy, jz)))
    RHS = -tadv1 * dpdy(jx, jy, jz)
    v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS + fy(jx, jy, jz)))
#else
    RHS = -tadv1 * dpdx(jx, jy, jz)
    u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS                 ))
    RHS = -tadv1 * dpdy(jx, jy, jz)
    v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS                 ))
#endif
end do
end do
end do

if (coord == 0) then
    jz_min = 2
else
    jz_min = 1
end if

do jz = jz_min, nz - 1
do jy = 1, ny
do jx = 1, nx
#ifdef PPLVLSET
    RHS = -tadv1 * dpdz(jx, jy, jz)
    w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS + fz(jx, jy, jz)))
#else
    RHS = -tadv1 * dpdz(jx, jy, jz)
    w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS                 ))
#endif
end do
end do
end do
#ifdef ENABLE_CUDA
end if
#endif

call apply_inflow()

#ifdef ENABLE_CUDA
if (project_stage_enabled) then
#ifdef PPMPI
    project_t1 = mpi_wtime()
#else
    call cpu_time(project_t1)
#endif
    project_stage_update = project_t1 - project_t0
    project_t0 = project_t1
end if
#endif

!--left this stuff last, so BCs are still enforced, no matter what
!  inflow_cond does

#ifdef PPMPI
    ! Exchange ghost node information (since coords overlap)
    if (nproc > 1) then
#ifdef ENABLE_CUDA
        if (project_cuda_enabled() .and. project_direct_halo_enabled() .and.   &
            project_overlap_boundary_enabled()) then
            call project_sync_velocity_direct_halos_overlap_cuda(u, v, w,      &
                project_halo_already_synced)
            project_boundary_done = .true.
        else if (project_cuda_enabled() .and. project_direct_halo_enabled()) then
            call project_sync_velocity_direct_halos_cuda(u, v, w)
        else if (project_cuda_enabled() .and. project_packed_halo_enabled()) then
            call project_sync_velocity_halos_cuda(u, v, w)
        else
#endif
    call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP )
#ifdef ENABLE_CUDA
        end if
#endif
    end if
#endif

#ifdef ENABLE_CUDA
if (project_stage_enabled) then
#ifdef PPMPI
    project_t1 = mpi_wtime()
#else
    call cpu_time(project_t1)
#endif
    project_stage_halo = project_t1 - project_t0
    project_t0 = project_t1
end if
#endif

!--enfore bc at top
#ifdef ENABLE_CUDA
if (.not. project_boundary_done) then
#endif
#ifdef PPMPI
if (coord == nproc-1) then
#endif
#ifdef ENABLE_CUDA
    if (project_cuda_enabled()) then
        if (ubc_mom == 0) then
            !$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, nx
                u(jx,jy,nz) = u(jx,jy,nz-1)
                v(jx,jy,nz) = v(jx,jy,nz-1)
            end do
            end do
        end if

        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, nx
            w(jx,jy,nz) = 0._rprec
        end do
        end do
    else
#endif
    ! Note: for ubc_mom > 0, u and v and nz will be written to output as BOGUS
    if (ubc_mom == 0) then    ! no-stress top
        u(:,:,nz) = u(:,:,nz-1)
        v(:,:,nz) = v(:,:,nz-1)
    endif
    ! no permeability
    w(:, :, nz)=0._rprec
#ifdef ENABLE_CUDA
    end if
#endif
#ifdef PPMPI
endif
#endif

if (coord == 0) then
  ! No modulation of u and v since if a stress free condition (lbc_mom=0) is
  ! applied, it is applied through the momentum equation.

  ! no permeability
#ifdef ENABLE_CUDA
  if (project_cuda_enabled()) then
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, nx
        w(jx,jy,1) = 0._rprec
    end do
    end do
  else
#endif
  w(:, :, 1)=0._rprec
#ifdef ENABLE_CUDA
  end if
#endif
end if

#ifdef ENABLE_CUDA
if (project_cuda_enabled()) call forcing_cuda_sync('boundary planes')
#endif

#ifdef ENABLE_CUDA
end if

if (project_stage_enabled) then
#ifdef PPMPI
    project_t1 = mpi_wtime()
#else
    call cpu_time(project_t1)
#endif
    project_stage_boundary = project_t1 - project_t0
#ifdef PPMPI
    call project_stage_report(project_stage_count, project_stage_update,       &
        project_stage_halo, project_stage_boundary)
#else
    write(*,'(a,i8)') 'Projection stage timing, call ', project_stage_count
    write(*,'(1a,E15.7)') '  Update + inflow: ', project_stage_update
    write(*,'(1a,E15.7)') '  Velocity halo: ', project_stage_halo
    write(*,'(1a,E15.7)') '  Boundary/final sync: ', project_stage_boundary
    write(*,'(1a,E15.7)') '  Stage sum: ', project_stage_update +             &
        project_stage_halo + project_stage_boundary
#endif
end if
#endif

end subroutine project

#ifdef ENABLE_CUDA
#ifdef PPMPI
!*******************************************************************************
subroutine project_sync_velocity_halos_cuda(u, v, w)
!*******************************************************************************
!
! Exchange projection halos for u/v/w using two packed GPU buffers.  This keeps
! the same slab-decomposition semantics as mpi_sync_real_array(...DOWNUP) while
! reducing MPI latency on multi-rank GPU runs.
!
use param, only : ld, nx, ny, nz, nproc, coord, comm, up, down, ierr,          &
    MPI_RPREC, MPI_STATUS_SIZE, MPI_REQUEST_NULL
implicit none

real(rprec), managed, intent(inout) :: u(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: v(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: w(ld,ny,0:nz)
integer :: plane_size
integer :: req(4)
integer :: statuses(MPI_STATUS_SIZE, 4)

plane_size = nx * ny
call project_ensure_halo_buffers(3 * plane_size)
call project_pack_velocity_halos_cuda(u, v, w, project_halo_send_down,         &
    project_halo_send_up)
call forcing_cuda_sync('project packed halo pack')

req = MPI_REQUEST_NULL
call mpi_irecv(project_halo_recv_up(1), 3*plane_size, MPI_RPREC, up, 171,      &
    comm, req(1), ierr)
call mpi_irecv(project_halo_recv_down(1), 3*plane_size, MPI_RPREC, down, 172,  &
    comm, req(2), ierr)
call mpi_isend(project_halo_send_down(1), 3*plane_size, MPI_RPREC, down, 171,  &
    comm, req(3), ierr)
call mpi_isend(project_halo_send_up(1), 3*plane_size, MPI_RPREC, up, 172,      &
    comm, req(4), ierr)
call mpi_waitall(4, req, statuses, ierr)

call project_unpack_velocity_halos_cuda(project_halo_recv_down,                &
    project_halo_recv_up, u, v, w)

end subroutine project_sync_velocity_halos_cuda

!*******************************************************************************
subroutine project_sync_velocity_direct_halos_cuda(u, v, w)
!*******************************************************************************
!
! Exchange projection halos for u/v/w by posting all contiguous halo-plane
! transfers together.  This avoids the pack/unpack kernels used by
! project_sync_velocity_halos_cuda while reducing the sequential latency of
! three separate mpi_sync_real_array calls.
!
use param, only : ld, ny, nz, nproc, coord, comm, up, down, ierr, MPI_RPREC,  &
    MPI_STATUS_SIZE, MPI_REQUEST_NULL
implicit none

real(rprec), managed, intent(inout) :: u(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: v(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: w(ld,ny,0:nz)
integer :: plane_size
integer :: req(12)
integer :: statuses(MPI_STATUS_SIZE, 12)
integer :: nreq

plane_size = ld * ny
call forcing_cuda_sync('project direct halo before MPI')

req = MPI_REQUEST_NULL
nreq = 0

if (coord < nproc - 1) then
    nreq = nreq + 1
    call mpi_irecv(u(1,1,nz), plane_size, MPI_RPREC, up, 181, comm, req(nreq),&
        ierr)
    nreq = nreq + 1
    call mpi_irecv(v(1,1,nz), plane_size, MPI_RPREC, up, 183, comm, req(nreq),&
        ierr)
    nreq = nreq + 1
    call mpi_irecv(w(1,1,nz), plane_size, MPI_RPREC, up, 185, comm, req(nreq),&
        ierr)
end if

if (coord > 0) then
    nreq = nreq + 1
    call mpi_irecv(u(1,1,0), plane_size, MPI_RPREC, down, 182, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_irecv(v(1,1,0), plane_size, MPI_RPREC, down, 184, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_irecv(w(1,1,0), plane_size, MPI_RPREC, down, 186, comm,           &
        req(nreq), ierr)
end if

if (coord > 0) then
    nreq = nreq + 1
    call mpi_isend(u(1,1,1), plane_size, MPI_RPREC, down, 181, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(v(1,1,1), plane_size, MPI_RPREC, down, 183, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(w(1,1,1), plane_size, MPI_RPREC, down, 185, comm,           &
        req(nreq), ierr)
end if

if (coord < nproc - 1) then
    nreq = nreq + 1
    call mpi_isend(u(1,1,nz-1), plane_size, MPI_RPREC, up, 182, comm,          &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(v(1,1,nz-1), plane_size, MPI_RPREC, up, 184, comm,          &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(w(1,1,nz-1), plane_size, MPI_RPREC, up, 186, comm,          &
        req(nreq), ierr)
end if

if (nreq > 0) then
    call mpi_waitall(nreq, req, statuses, ierr)
end if

end subroutine project_sync_velocity_direct_halos_cuda

!*******************************************************************************
subroutine project_sync_velocity_direct_halos_overlap_cuda(u, v, w,            &
    already_synced)
!*******************************************************************************
!
! Direct halo exchange with physical top/bottom boundary kernels launched while
! MPI is progressing.  The boundary kernels touch only physical domain planes,
! not the exchanged ghost planes or send planes.
!
use param, only : ld, nx, ny, nz, nproc, coord, comm, up, down, ierr,          &
    MPI_RPREC, MPI_STATUS_SIZE, MPI_REQUEST_NULL, ubc_mom
implicit none

real(rprec), managed, intent(inout) :: u(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: v(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: w(ld,ny,0:nz)
logical, intent(in) :: already_synced
integer :: plane_size
integer :: req(12)
integer :: statuses(MPI_STATUS_SIZE, 12)
integer :: nreq
integer :: jx, jy

plane_size = ld * ny
if (.not. already_synced) call forcing_cuda_sync('project direct overlap before MPI')

req = MPI_REQUEST_NULL
nreq = 0

if (coord < nproc - 1) then
    nreq = nreq + 1
    call mpi_irecv(u(1,1,nz), plane_size, MPI_RPREC, up, 181, comm, req(nreq),&
        ierr)
    nreq = nreq + 1
    call mpi_irecv(v(1,1,nz), plane_size, MPI_RPREC, up, 183, comm, req(nreq),&
        ierr)
    nreq = nreq + 1
    call mpi_irecv(w(1,1,nz), plane_size, MPI_RPREC, up, 185, comm, req(nreq),&
        ierr)
end if

if (coord > 0) then
    nreq = nreq + 1
    call mpi_irecv(u(1,1,0), plane_size, MPI_RPREC, down, 182, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_irecv(v(1,1,0), plane_size, MPI_RPREC, down, 184, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_irecv(w(1,1,0), plane_size, MPI_RPREC, down, 186, comm,           &
        req(nreq), ierr)
end if

if (coord > 0) then
    nreq = nreq + 1
    call mpi_isend(u(1,1,1), plane_size, MPI_RPREC, down, 181, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(v(1,1,1), plane_size, MPI_RPREC, down, 183, comm,           &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(w(1,1,1), plane_size, MPI_RPREC, down, 185, comm,           &
        req(nreq), ierr)
end if

if (coord < nproc - 1) then
    nreq = nreq + 1
    call mpi_isend(u(1,1,nz-1), plane_size, MPI_RPREC, up, 182, comm,          &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(v(1,1,nz-1), plane_size, MPI_RPREC, up, 184, comm,          &
        req(nreq), ierr)
    nreq = nreq + 1
    call mpi_isend(w(1,1,nz-1), plane_size, MPI_RPREC, up, 186, comm,          &
        req(nreq), ierr)
end if

if (coord == nproc-1) then
    if (ubc_mom == 0) then
        !$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, nx
            u(jx,jy,nz) = u(jx,jy,nz-1)
            v(jx,jy,nz) = v(jx,jy,nz-1)
        end do
        end do
    end if

    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, nx
        w(jx,jy,nz) = 0._rprec
    end do
    end do
end if

if (coord == 0) then
    !$cuf kernel do(2) <<<*,*>>>
    do jy = 1, ny
    do jx = 1, nx
        w(jx,jy,1) = 0._rprec
    end do
    end do
end if

if (nreq > 0) then
    call mpi_waitall(nreq, req, statuses, ierr)
end if

call forcing_cuda_sync('project overlap boundary')

end subroutine project_sync_velocity_direct_halos_overlap_cuda

!*******************************************************************************
subroutine project_stage_report(stage_count, update_inflow, halo, boundary)
!*******************************************************************************
use param, only : coord
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC
use mpi
#endif
implicit none

integer, intent(in) :: stage_count
real(rprec), intent(in) :: update_inflow, halo, boundary
real(rprec) :: update_max, halo_max, boundary_max, total

#ifdef PPMPI
call mpi_allreduce(update_inflow, update_max, 1, MPI_RPREC, MPI_MAX, comm,     &
    ierr)
call mpi_allreduce(halo, halo_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
call mpi_allreduce(boundary, boundary_max, 1, MPI_RPREC, MPI_MAX, comm, ierr)
#else
update_max = update_inflow
halo_max = halo
boundary_max = boundary
#endif

total = update_max + halo_max + boundary_max

if (coord == 0) then
    write(*,'(a,i8)') 'Projection stage timing (max rank), call ', stage_count
    write(*,'(1a,E15.7)') '  Update + inflow: ', update_max
    write(*,'(1a,E15.7)') '  Velocity halo: ', halo_max
    write(*,'(1a,E15.7)') '  Boundary/final sync: ', boundary_max
    write(*,'(1a,E15.7)') '  Stage sum: ', total
end if

end subroutine project_stage_report

!*******************************************************************************
subroutine project_ensure_halo_buffers(nitems)
!*******************************************************************************
implicit none

integer, intent(in) :: nitems

if (allocated(project_halo_send_down)) then
    if (size(project_halo_send_down) /= nitems) then
        deallocate(project_halo_send_down, project_halo_send_up,                &
            project_halo_recv_down, project_halo_recv_up)
    end if
end if

if (.not. allocated(project_halo_send_down)) then
    allocate(project_halo_send_down(nitems), project_halo_send_up(nitems),      &
        project_halo_recv_down(nitems), project_halo_recv_up(nitems))
end if

end subroutine project_ensure_halo_buffers

!*******************************************************************************
subroutine project_pack_velocity_halos_cuda(u, v, w, send_down, send_up)
!*******************************************************************************
use param, only : ld, nx, ny, nz
implicit none

real(rprec), managed, intent(in) :: u(ld,ny,0:nz)
real(rprec), managed, intent(in) :: v(ld,ny,0:nz)
real(rprec), managed, intent(in) :: w(ld,ny,0:nz)
real(rprec), device, intent(inout) :: send_down(3*nx*ny)
real(rprec), device, intent(inout) :: send_up(3*nx*ny)
integer :: idx, jx, jy, plane_size

plane_size = nx * ny
!$cuf kernel do(1) <<<*,*>>>
do idx = 1, plane_size
    jx = mod(idx - 1, nx) + 1
    jy = (idx - 1) / nx + 1

    send_down(idx) = u(jx,jy,1)
    send_down(plane_size + idx) = v(jx,jy,1)
    send_down(2*plane_size + idx) = w(jx,jy,1)

    send_up(idx) = u(jx,jy,nz-1)
    send_up(plane_size + idx) = v(jx,jy,nz-1)
    send_up(2*plane_size + idx) = w(jx,jy,nz-1)
end do

end subroutine project_pack_velocity_halos_cuda

!*******************************************************************************
subroutine project_unpack_velocity_halos_cuda(recv_down, recv_up, u, v, w)
!*******************************************************************************
use param, only : ld, nx, ny, nz, nproc, coord
implicit none

real(rprec), device, intent(in) :: recv_down(3*nx*ny)
real(rprec), device, intent(in) :: recv_up(3*nx*ny)
real(rprec), managed, intent(inout) :: u(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: v(ld,ny,0:nz)
real(rprec), managed, intent(inout) :: w(ld,ny,0:nz)
integer :: idx, jx, jy, plane_size

plane_size = nx * ny
!$cuf kernel do(1) <<<*,*>>>
do idx = 1, plane_size
    jx = mod(idx - 1, nx) + 1
    jy = (idx - 1) / nx + 1

    if (coord < nproc - 1) then
        u(jx,jy,nz) = recv_up(idx)
        v(jx,jy,nz) = recv_up(plane_size + idx)
        w(jx,jy,nz) = recv_up(2*plane_size + idx)
    end if

    if (coord > 0) then
        u(jx,jy,0) = recv_down(idx)
        v(jx,jy,0) = recv_down(plane_size + idx)
        w(jx,jy,0) = recv_down(2*plane_size + idx)
    end if
end do

end subroutine project_unpack_velocity_halos_cuda
#endif
#endif

end module forcing
