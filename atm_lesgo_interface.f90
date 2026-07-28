!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Written by:
!!
!!   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
!!
!!   Copyright (C) 2012-2013, Johns Hopkins University
!!
!!   This file is part of The Actuator Turbine Model Library.
!!
!!   The Actuator Turbine Model is free software: you can redistribute it
!!   and/or modify it under the terms of the GNU General Public License as
!!   published by the Free Software Foundation, either version 3 of the
!!   License, or (at your option) any later version.
!!
!!   The Actuator Turbine Model is distributed in the hope that it will be
!!   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details.
!!
!!   You should have received a copy of the GNU General Public License
!!   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************************
module atm_lesgo_interface
!*******************************************************************************
! Navigation map for this large interface module:
!   - declarations and policy flags: module header before `contains`
!   - diagnostics/timing helpers: atm_diag_*, atm_lesgo_report_timing
!   - GPU utility kernels: atm_prepare_direct_w and the batched atPoint routines
!   - lifecycle and geometry: initialize/finalize, diagnostics, and findCells
!   - timestep entry point: atm_lesgo_forcing(phase)
!   - legacy/CPU gather and force paths: atm_lesgo_mpi_gather*,
!     atm_lesgo_force, atm_lesgo_convolute_force, atm_lesgo_apply_force
!
! Keep CPU and GPU versions semantically paired.  When changing a GPU routine,
! update the matching CPU/legacy path or document why the path intentionally
! differs in docs/gpu_module_contracts.md.
! This module interfaces actuator turbine module with lesgo
! It is a lesgo specific module, unlike the atm module
! The MPI management is done only in this section of the code
! This is very code dependent and will have to be modified according to
! the code being used. In this case LESGO has its own MPI details
! Look into mpi_defs.f90 for the details

! Remember to always dimensionalize the variables from LESGO
! Length is non-dimensionalized by z_i

! Lesgo data used regarding the grid (LESGO)
use param, only : dt ,nx,ny,nz,nz_tot,dx,dy,dz,coord,nproc, z_i, u_star, lbz,  &
                  total_time, jt_total, L_x, L_y, L_z
! nx, ny, nz - nodes in every direction
! z_i - non-dimensionalizing length
! dt - time-step

! These are the forces, and velocities on x,y, and z respectively
use sim_param, only : fxa, fya, fza, u, v, w

! Grid definition (LESGO)
use grid_m, only : grid

! MPI implementation from LESGO
#ifdef PPMPI
  use mpi
  use param, only : ierr, mpi_rprec, comm, up, down
#endif

! Interpolating function for interpolating the velocity field to each
! actuator point
use functions, only : trilinear_interp, interp_to_uv_grid

use clock_m, only : clock_t

! Actuator Turbine Model module
use atm_base, only : distance
use actuator_turbine_model, only : atm_computeBladeForce,                     &
    atm_computeNacelleForce, atm_compute_cl_correction, atm_initialize,        &
    atm_initialize_output, atm_output, atm_structure_enabled,                  &
    atm_structure_timing_report, atm_update, atm_write_restart
use atm_input_util, only : numberOfTurbines, outputInterval, rprec,            &
                           turbineArray, turbineModel, updateInterval

! Used for testing time
! use clock_m

implicit none

! Variable for interpolating the velocity in w onto the uv grid
#ifdef PPLES_GPU
! PPLES GPU ownership map:
!   - LES fields u/v/w and fxa/fya/fza are owned by sim_param and are expected
!     to be device-resident during timestep forcing.
!   - w_uv is host scratch used only by CPU and legacy ATM sampling.
!   - batched atPoint tables below own the production sampling/convolution
!     copies of turbine geometry and forces.
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), allocatable, save, dimension(:,:) :: atm_wuv_send_down
real(rprec), allocatable, save, dimension(:,:) :: atm_wuv_recv_up
!$acc declare create(atm_wuv_send_down, atm_wuv_recv_up)
#else
real(rprec), allocatable, dimension(:,:,:) :: w_uv
#endif

private
public atm_lesgo_initialize, atm_lesgo_forcing, atm_lesgo_finalize,           &
    atm_lesgo_checkpoint

! This is a list that stores all the points in the domain with a body
! force due to the turbines.
type bodyForce_t
    integer :: c ! Number of cells
    ! i,j,k stores the index for the point in the domain
    integer, allocatable :: ijk(:,:)
    real(rprec), allocatable :: force(:,:) ! Force vector on uv grid
    real(rprec), allocatable :: location(:,:) ! Position vector on uv grid
end type bodyForce_t

#ifdef PPLES_GPU
! ---- Batched atPoint GPU state (round 3) ----
! Static concatenated tables for ALL turbines (built once): force-field cell
! locations/ijk/turbine-id, grid axes + autowrap, prefix offsets. Per force
! step only the flattened blade points/forces (~260 KB) and the per-turbine
! constants are re-uploaded; sampling and convolution then run as ONE kernel
! over all turbines instead of 60 data-region/kernel/sync rounds.
integer, parameter :: ATM_NTC = 14   ! per-turbine constant slots
integer, parameter :: ATM_SAMPLING_ATPOINT = 1
integer, parameter :: ATM_SAMPLING_SPALART = 2
logical :: atm_batch_ready   = .false.
logical :: atm_batch_sampled = .false.
logical :: atm_atpoint_present = .false.
logical :: atm_spalart_present = .false.
logical :: atm_host_velocity_bridge_required = .false.
integer :: atm_nbp_tot = 0, atm_cUV_tot = 0, atm_cW_tot = 0
integer,     allocatable :: atm_bp_off(:)                  ! (nTurb+1) blade-point prefix
real(rprec), allocatable :: atm_bp_all(:,:), atm_bf_all(:,:)   ! (3, nbp_tot)
real(rprec), allocatable :: atm_velbp_all(:,:)             ! (3, nbp_tot) sampled velocity
integer,     allocatable :: atm_inr_all(:)                 ! (nbp_tot) in-domain flag
real(rprec), allocatable :: atm_tconst(:,:)                ! (ATM_NTC, nTurb)
integer,     allocatable :: atm_sampling_mode(:)           ! static per-turbine mode
real(rprec), allocatable :: atm_locUV_all(:,:), atm_locW_all(:,:)  ! (3, c_tot) static
integer,     allocatable :: atm_ijkUV_all(:,:), atm_ijkW_all(:,:)  ! (3, c_tot) static
integer,     allocatable :: atm_tidUV(:), atm_tidW(:)      ! per-cell turbine id, static
real(rprec), allocatable :: atm_gx(:), atm_gy(:), atm_gz(:)
integer,     allocatable :: atm_awi(:), atm_awj(:)
! Optional compatibility buffers. They are allocated only when a Spalart
! turbine exists, and scatter its host-computed force into resident LES arrays
! without transferring or overwriting the complete fxa/fya/fza fields.
real(rprec), allocatable :: atm_spalart_forceUV(:,:), atm_spalart_forceW(:)

! ---- Batched Cl/tip correction (round 4) ----
! GPU port of atm_compute_cl_correction: the O(N^2) blade-to-blade induced
! velocity loop (~68 ms/step on the host) batched over all turbines. The
! model arithmetic is replicated verbatim (same mixed-precision literals,
! same accumulation order); all turbineArray outputs are copied back to the
! host each step so host state stays exactly as the host routine leaves it.
logical :: atm_clc_ready = .false.
integer,     allocatable :: atm_pt_turb(:), atm_pt_q(:), atm_pt_base(:), atm_pt_qq(:)
real(rprec), allocatable :: atm_chord_all(:), atm_brad_all(:)      ! static
real(rprec), allocatable :: atm_db_all(:)                      ! static blade-section width
real(rprec), allocatable :: atm_clc_tc(:,:)                       ! (3,nTurb): eps_s, optEpsChord, active
real(rprec), allocatable :: atm_wv_all(:,:)                       ! (3,nbp) windVectors in
real(rprec), allocatable :: atm_cl_all(:), atm_cd_all(:), atm_vmag_all(:)
real(rprec), allocatable :: atm_du_all(:,:)                       ! (3,nbp) state in/out
real(rprec), allocatable :: atm_uyopt_vec_all(:,:)                ! (3,nbp) state in/out
real(rprec), allocatable :: atm_uinf_all(:,:), atm_uxles_all(:,:) ! (3,nbp) out
real(rprec), allocatable :: atm_g_all(:), atm_dg_all(:), atm_epsopt_all(:)
real(rprec), allocatable :: atm_uyles_vec_all(:,:)                ! (3,nbp) out
real(rprec), allocatable :: atm_uyles_all(:), atm_uyopt_all(:)    ! magnitudes out
#endif

! Body force field
type(bodyForce_t), allocatable, target, dimension(:) :: forceFieldUV, forceFieldW

! The very crucial parameter pi
real(rprec), parameter :: pi=acos(-1._rprec)

type(clock_t), save :: atm_clock_interp_w, atm_clock_update, atm_clock_reset
type(clock_t), save :: atm_clock_sample, atm_clock_force, atm_clock_gather, atm_clock_convolve
type(clock_t), save :: atm_clock_clcorr, atm_clock_apply, atm_clock_output
real(rprec), save :: atm_time_interp_w = 0._rprec
real(rprec), save :: atm_time_update = 0._rprec
real(rprec), save :: atm_time_reset = 0._rprec
real(rprec), save :: atm_time_sample = 0._rprec
real(rprec), save :: atm_time_force = 0._rprec
real(rprec), save :: atm_time_gather = 0._rprec
real(rprec), save :: atm_time_convolve = 0._rprec
real(rprec), save :: atm_time_clcorr = 0._rprec
real(rprec), save :: atm_time_apply = 0._rprec
real(rprec), save :: atm_time_output = 0._rprec
integer, save :: atm_forcing_calls = 0
integer, save :: atm_last_checkpoint_step = -1
logical, save :: atm_diag_load_printed = .false.
logical, save :: atm_force_state_valid = .false.

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_interp_w_to_uv()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

end subroutine atm_interp_w_to_uv

#if defined(PPLES_GPU)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_prepare_direct_w()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer :: jx, jy

#ifdef PPMPI
integer :: status(MPI_STATUS_SIZE)
#endif

if (nproc <= 1) then
    ! No halo is needed on a single rank. The deferred LES queue only reads
    ! u/v/w at this point, so synchronizing it would unnecessarily serialize
    ! ATM phase 1 against convection. Keep one-element device-present buffers
    ! because the unified sampling kernel contains an unreachable MPI-halo
    ! branch and OpenACC still resolves its allocatable descriptors.
    if (.not. allocated(atm_wuv_send_down)) then
        allocate(atm_wuv_send_down(1,1), atm_wuv_recv_up(1,1))
    endif
    return
endif

if (.not. allocated(atm_wuv_send_down)) then
    allocate(atm_wuv_send_down(nx,ny), atm_wuv_recv_up(nx,ny))
    !$acc update device(atm_wuv_send_down, atm_wuv_recv_up)
endif

#ifdef PPMPI
! The direct sampler only needs the old w_uv(nz) halo from the rank above.
!$acc parallel loop collapse(2) default(present)
do jy = 1, ny
do jx = 1, nx
    atm_wuv_send_down(jx,jy) = 0.5_rprec * (w(jx,jy,1) + w(jx,jy,2))
end do
end do

! The boundary-pack loop has no async clause, so it is complete before the
! following GPU-aware MPI call.  Avoid a global OpenACC wait here because it can
! charge unrelated deferred LES work to the ATM direct-w timer.
!$acc host_data use_device(atm_wuv_send_down, atm_wuv_recv_up)
call mpi_sendrecv(atm_wuv_send_down(1,1), nx*ny, mpi_rprec, down, 991,        &
                  atm_wuv_recv_up(1,1), nx*ny, mpi_rprec, up, 991,            &
                  comm, status, ierr)
!$acc end host_data
if (ierr /= 0) stop 'ATM direct w device halo exchange failed'

! Blocking GPU-aware MPI sendrecv completes the device receive buffer before the
! following ATM sampling kernels are launched. Avoid a global OpenACC wait here:
! it can drain unrelated deferred LES work into the direct-w timer.
#endif

end subroutine atm_prepare_direct_w
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_initialize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

! Counter to establish number of points which are influenced by body forces
integer ::  m

! CPU sampling always uses host w_uv. The GPU build allocates it later only
! when Spalart or nacelle compatibility requires host interpolation.
#ifndef PPLES_GPU
allocate(w_uv(nx,ny,lbz:nz))
w_uv = 0._rprec
#endif

atm_last_checkpoint_step = -1
atm_force_state_valid = .false.
call atm_initialize(jt_total, total_time) ! Initialize the ATM state

#ifdef PPLES_GPU
! The optimized atPoint path remains fully device-resident. Spalart sampling
! and nacelle interpolation are established host algorithms, so GPU builds
! activate an explicit velocity bridge only when either feature is requested.
atm_atpoint_present = .false.
atm_spalart_present = .false.
atm_host_velocity_bridge_required = .false.
do m = 1, numberOfTurbines
    select case (trim(turbineArray(m) % sampling))
    case ('atPoint')
        atm_atpoint_present = .true.
    case ('Spalart')
        atm_spalart_present = .true.
        atm_host_velocity_bridge_required = .true.
    case default
        if (coord == 0) then
            write(*,'(1a,i0,2a)') 'Unsupported ATM sampling mode for turbine ', &
                m, ': ', trim(turbineArray(m) % sampling)
        endif
        error stop 'ATM sampling must be atPoint or Spalart'
    end select
    if (turbineArray(m) % nacelle) atm_host_velocity_bridge_required = .true.
enddo
if (coord == 0 .and. atm_host_velocity_bridge_required) then
    write(*,'(1a)') 'ATM GPU: enabling legacy host velocity compatibility bridge'
endif
if (atm_host_velocity_bridge_required) then
    allocate(w_uv(nx,ny,lbz:nz))
    w_uv = 0._rprec
endif
#endif

! Allocate the body force variables. It is an array with one per turbine.
allocate(forceFieldUV(numberOfTurbines))
allocate(forceFieldW(numberOfTurbines))

    do m=1, numberOfTurbines
        call atm_lesgo_findCells(m)
    enddo

    #ifdef PPMPI
        call mpi_barrier( comm, ierr )
    #endif


#ifdef PPMPI
    ! This will create the output files and write initialization to the screen
    if (coord == 0) then
        call atm_initialize_output()
    endif
#else
    call atm_initialize_output()
#endif


end subroutine atm_lesgo_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_checkpoint(checkpoint_step, checkpoint_time)
! Write ATM state at the same logical timestep as the LES field checkpoint.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: checkpoint_step
real(rprec), intent(in) :: checkpoint_time
integer :: i

if (checkpoint_step == atm_last_checkpoint_step) return

do i = 1, numberOfTurbines
    if (coord == turbineArray(i)%master) then
        call atm_write_restart(i, checkpoint_step, checkpoint_time)
    endif
enddo

#ifdef PPMPI
call mpi_barrier(comm, ierr)
#endif
atm_last_checkpoint_step = checkpoint_step

end subroutine atm_lesgo_checkpoint

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_finalize ()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialize the actuator turbine model
implicit none

call atm_lesgo_report_timing()
call atm_structure_timing_report()

! output_final normally wrote this timestep already. Keep finalization safe for
! callers that bypass the normal output driver; duplicate writes are skipped.
call atm_lesgo_checkpoint(jt_total, total_time)

! Write if on main node
if (coord == 0) then
    write(*,*) 'Finalizing ATM...'
endif

if (coord == 0) then
    write(*,*) 'Done finalizing ATM'
endif

end subroutine atm_lesgo_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_report_timing()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

real(rprec) :: vals(10), maxvals(10), measured_total

vals = (/ atm_time_interp_w, atm_time_update, atm_time_reset, atm_time_sample, &
          atm_time_force, atm_time_gather, atm_time_convolve, atm_time_clcorr, &
          atm_time_apply, atm_time_output /)

#ifdef PPMPI
call mpi_allreduce(vals, maxvals, size(vals), mpi_rprec, mpi_max, comm, ierr)
#else
maxvals = vals
#endif

measured_total = sum(maxvals)

if (coord == 0) then
    write(*,*) '==================================================='
    write(*,*) 'ATM Cumulative Times (s, max over ranks):'
    write(*,'(1a,I8)')    '  ATM forcing calls: ', atm_forcing_calls
    write(*,'(1a,E15.7)') '  w -> uv interpolation: ', maxvals(1)
    write(*,'(1a,E15.7)') '  turbine update/yaw/rotation: ', maxvals(2)
    write(*,'(1a,E15.7)') '  turbine/force reset: ', maxvals(3)
    write(*,'(1a,E15.7)') '  batched velocity sampling: ', maxvals(4)
    write(*,'(1a,E15.7)') '  blade/nacelle force: ', maxvals(5)
    write(*,'(1a,E15.7)') '  MPI gather: ', maxvals(6)
    write(*,'(1a,E15.7)') '  force convolution: ', maxvals(7)
    write(*,'(1a,E15.7)') '  tip correction: ', maxvals(8)
    write(*,'(1a,E15.7)') '  apply force to grid: ', maxvals(9)
    write(*,'(1a,E15.7)') '  ATM output: ', maxvals(10)
    write(*,'(1a,E15.7)') '  ATM measured subtotal: ', measured_total
    write(*,*) '==================================================='
end if


end subroutine atm_lesgo_report_timing


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_findCells (m)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine finds all the cells that surround the turbines

! The awkward if statements are to only consider points in front and behind
! the turbine without having to

implicit none

! The turbine number
integer, intent(in) :: m

! Counter to establish number of points which are influenced by body forces
integer :: cUV, cW  ! Counters for number of points affected on UV and W grids
integer :: i, j, k

! Vector used to store x, y, z locations
real(rprec), dimension(3) :: vector_point
! These are the pointers to the grid arrays
real(rprec), pointer, dimension(:) :: x,y,z,zw

! Variables for MPI implementation
#ifdef PPMPI
integer :: base_group ! The base group from comm --> MPI_COMM_WORLD (all processors)
integer :: local_group  ! The local group of processors
integer :: member !  (1 or 0) yes or no
integer :: num_of_members  ! total number of members
#endif

! List of all the cores that belong to this turbine
! This variable gets allocated for each turbine
integer, allocatable, dimension(:) :: ls_of_cores

nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

! Initialize internal counter to zero
forceFieldUV(m) % c = 0

! This will find all the locations that are influenced by each turbine
! It depends on a sphere centered on the rotor that extends beyond the blades
cUV=0  ! Initialize conuter
cW=0  ! Initialize conuter
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i

            ! Take into account the UV grid
            vector_point(3)=z(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cUV=cUV+1
!~ endif

                end if
                ! Take into account the W grid
                vector_point(3)=zw(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cW=cW+1
!~ endif
                end if
        enddo
    enddo
enddo

! Allocate space for the force fields in UV and W grids
forceFieldUV(m) % c = cUV  ! Counter
allocate(forceFieldUV(m) % force(3,cUV))
allocate(forceFieldUV(m) % location(3,cUV))
allocate(forceFieldUV(m) % ijk(3,cUV))

forceFieldW(m) % c = cW  ! Counter
allocate(forceFieldW(m) % force(3,cW))
allocate(forceFieldW(m) % location(3,cW))
allocate(forceFieldW(m) % ijk(3,cW))

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif
write(*,*) 'Number of cells being affected by ATM in turbine', m,              &
           ' cUV, cW = ', cUV, cW
#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

cUV=0
cW=0
! Run the same loop and save all variables
! The forceField arrays include all the forces which affect the domain
do i=1,nx ! Loop through grid points in x
    do j=1,ny ! Loop through grid points in y
        do k=1,nz ! Loop through grid points in z
            vector_point(1)=x(i)*z_i ! z_i used to dimensionalize LESGO
            vector_point(2)=y(j)*z_i
            vector_point(3)=z(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cUV=cUV+1
                    forceFieldUV(m) % ijk(1,cUV) = i
                    forceFieldUV(m) % ijk(2,cUV) = j
                    forceFieldUV(m) % ijk(3,cUV) = k
                    forceFieldUV(m) % location(1:3,cUV) = vector_point(1:3)
                    forceFieldUV(m) % force(1:3,cUV) = 0_rprec
                endif
!~ endif
            vector_point(3)=zw(k)*z_i
                if (distance(vector_point,turbineArray(m) %                    &
                    towerShaftIntersect)                                       &
                    .le. turbineArray(m) % sphereRadius ) then
!~ if ( ( (vector_point(1) - turbineArray(m) % towerShaftIntersect(1) )**2 ) <= ( turbineArray(m) % projectionRadius**2 )) then
                    cW=cW+1
                    forceFieldW(m) % ijk(1,cW) = i
                    forceFieldW(m) % ijk(2,cW) = j
                    forceFieldW(m) % ijk(3,cW) = k
                    forceFieldW(m) % location(1:3,cW) = vector_point(1:3)
                    forceFieldW(m) % force(:,cW) = 0_rprec
                endif
!~ endif
        enddo
    enddo
enddo


! MPI distribution
! This will create new communicator for each turbine
#ifdef PPMPI

! Store the base group from the global communicator mpi_comm_world
call MPI_COMM_GROUP(comm, base_group, ierr)

! Assign member
member = 0
! Flag to know if this turbine is operating or not
turbineArray(m) % operate = .FALSE.

! Assign proper values if turbine affects processors in this region
if (cUV > 0 .or. cW >0) then
member = 1
turbineArray(m) % operate = .TRUE.
endif

! Find the total number of processors for each turbine
call mpi_allreduce(member, num_of_members, 1, MPI_INTEGER , MPI_SUM, comm, ierr)

if (turbineArray(m) % operate) then
! Find the master processor for each turbine
    call mpi_allreduce(coord, turbineArray(m) % master, 1, MPI_INTEGER ,       &
                       MPI_MIN, comm, ierr)
else
    ! This is bogus since nz will always be less than number of processors
    ! This is done to ensure that the master is part of the processors
    ! that hold the turbine model
    call mpi_allreduce(nz_tot, turbineArray(m) % master, 1, MPI_INTEGER ,       &
                       MPI_MIN, comm, ierr)
endif

allocate(ls_of_cores(num_of_members))
ls_of_cores(1) = turbineArray(m) % master

! Notice this list is valid only for decomposition in 1 direction
do i = 2, num_of_members
    ls_of_cores(i) = ls_of_cores(i-1) + 1
enddo

! Write if this processor is the master
if (coord == turbineArray(m) % master) then
    write(*,*) 'Master for turbine',m, 'is processor', turbineArray(m) % master
endif

! Create the new communicator and group for this turbine
call MPI_GROUP_INCL(base_group, num_of_members, ls_of_cores, local_group, ierr)
call MPI_COMM_CREATE(comm, local_group, turbineArray(m) % TURBINE_COMM_WORLD, ierr)

if (turbineArray(m) % operate) then
    write(*,*) 'Processor', coord, 'has elements in turbine', m
else
    write(*,*) 'Processor', coord, 'does NOT have elements in turbine', m
endif

    call mpi_barrier( comm, ierr )

#endif

end subroutine atm_lesgo_findCells

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_forcing (phase)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutines calls the update function from the ATM Library
! and calculates the body forces needed in the domain
!
! Optional phase split (explicit-residency overlap): phase=1 prepares required
! velocity data, updates blades, samples velocity, evaluates the host airfoil
! model, and gathers; phase=2 deposits force, applies correction, and writes
! output. Called without phase, both portions run in order. The split lets
! phase 1 host work overlap queued SGS/convection kernels.
implicit none

integer, intent(in), optional :: phase
integer :: ph
integer :: i
logical :: atm_force_update_step, atm_output_step

!~ real(rprec) :: integrateNacelleForce, totForce
!~ integer :: c

!~ type(clock_t) :: myClock

ph = 0
if (present(phase)) ph = phase
! A fresh process has no cached grid force, including after a restart. Force
! one complete update before honoring updateInterval so the first continued
! timestep cannot apply an uninitialized/zero turbine field.
atm_force_update_step = mod(jt_total-1, updateInterval) == 0 .or.              &
                        .not. atm_force_state_valid

if (ph /= 2) then

atm_forcing_calls = atm_forcing_calls + 1

if (atm_force_update_step) then
    ! Prepare velocities only when the turbine force model will consume them.
    call atm_clock_interp_w%start()
#ifdef PPLES_GPU
    ! The GPU sampler reads w directly and only prepares a top-slab halo when
    ! MPI decomposition requires one. Spalart/nacelle compatibility is the
    ! only GPU configuration that refreshes host velocities and constructs
    ! host w_uv.
    if (atm_atpoint_present) call atm_prepare_direct_w()
    if (atm_host_velocity_bridge_required) then
        !$acc wait(1)
        !$acc update self(u, v, w)
        call atm_interp_w_to_uv()
    endif
#else
    call atm_interp_w_to_uv()
#endif
    call atm_clock_interp_w%stop()
    atm_time_interp_w = atm_time_interp_w + atm_clock_interp_w%time
endif


! Update the blade positions based on the time-step
! Time needs to be dimensionalized
! All processors carry the blade points
!~ call myCock%start_time();
!~ call atm_update(dt*z_i/u_star)

! Loop through all turbines and rotate the blades
call atm_clock_update%start()
do i = 1, numberOfTurbines
    ! If statement is for running code only with the processors on that turbine
        if (turbineArray(i) % operate) then
            ! Time is dimensionalize using velocity and length scale
            call atm_update(i, dt*z_i/u_star, total_time*z_i/u_star)
        endif
    enddo
call atm_clock_update%stop()
atm_time_update = atm_time_update + atm_clock_update%time

!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Update ', myClock % time

end if   ! ph /= 2  (head of phase 1)

! Only calculate new forces if interval is correct.
if (atm_force_update_step) then

    if (ph /= 2) then

#ifdef PPLES_GPU
    ! Batched device sampling for all atPoint turbines: one kernel + one D2H
    ! per force step instead of a data region + kernel + sync per turbine.
    ! Reads only bladePoints (already rotated by atm_update above) and the
    ! device u/v/w, so it is safe to run before the reset loop.
    if (atm_atpoint_present) then
        call atm_clock_sample%start()
        call atm_batch_sample_velocity_gpu()
        call atm_clock_sample%stop()
        atm_time_sample = atm_time_sample + atm_clock_sample%time
    endif
#endif

    ! Establish all turbine properties as zero
    ! This is essential for paralelization
    do i=1,numberOfTurbines
        call atm_clock_reset%start()
        turbineArray(i) % torqueRotor = 0._rprec
        turbineArray(i) % thrust = 0._rprec
        turbineArray(i) % nacelleForce = 0._rprec
        turbineArray(i) % VelNacelle_sampled = 0._rprec
        turbineArray(i) % VelNacelle_corrected = 0._rprec
        turbineArray(i) % bladeForces = 0._rprec
        turbineArray(i) % integratedBladeForces = 0._rprec
        turbineArray(i) % alpha = 0._rprec
        turbineArray(i) % Cd = 0._rprec
        turbineArray(i) % Cm = 0._rprec
        turbineArray(i) % Cl = 0._rprec
        turbineArray(i) % Cl_b = 0._rprec
        turbineArray(i) % G = 0._rprec
        turbineArray(i) % lift = 0._rprec
        turbineArray(i) % drag = 0._rprec
        turbineArray(i) % Vmag = 0._rprec
        turbineArray(i) % windVectors = 0._rprec
        turbineArray(i) % induction_a = 0._rprec
        turbineArray(i) % u_infinity = 0._rprec
        turbineArray(i) % bladeAlignedVectors = 0._rprec
        turbineArray(i) % axialForce = 0._rprec
        turbineArray(i) % tangentialForce = 0._rprec
        turbineArray(i) % pitchingMoment = 0._rprec

        ! Applied grid forces are overwritten by convolution below; clearing
        ! the large per-cell host arrays here would add avoidable memory traffic.
        call atm_clock_reset%stop()
        atm_time_reset = atm_time_reset + atm_clock_reset%time

        if (turbineArray(i) % operate) then
            ! Calculate forces for all turbines
            call atm_clock_force%start()
            call atm_lesgo_force(i)
            call atm_clock_force%stop()
            atm_time_force = atm_time_force + atm_clock_force%time

        endif

    enddo
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Forces ', myClock % time


!~  call myClock % start()
! This will gather all the blade forces from all processors
#ifdef PPMPI
    ! This will gather all values used in MPI
!~     call mpi_barrier( MPI_COMM_WORLD, ierr )
    if (nproc > 1) then
        call atm_clock_gather%start()
        call atm_lesgo_mpi_gather()
        call atm_clock_gather%stop()
        atm_time_gather = atm_time_gather + atm_clock_gather%time
    endif
!~     call mpi_barrier( MPI_COMM_WORLD, ierr )

#endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  MPI Gather ', myClock % time

    ! The force/basis state produced above is the load consumed by the next
    ! structural update. Fresh and legacy restarts skip structural advancement
    ! until this first complete load has been assembled.
    if (atm_structure_enabled()) then
        do i = 1, numberOfTurbines
            if (turbineArray(i)%operate) then
                turbineArray(i)%structure_load_history_valid = .true.
            endif
        enddo
    endif

    end if   ! ph /= 2  (sampling + blade force + gather)

    if (ph == 1) return

#ifdef PPLES_GPU
    ! Batched OpenACC Cl/tip correction updates the same induced-velocity
    ! state used by both rigid and structural turbine consumers.
    call atm_clock_clcorr%start()
    call atm_batch_cl_correction_gpu()
    call atm_clock_clcorr%stop()
    atm_time_clcorr = atm_time_clcorr + atm_clock_clcorr%time
#endif

    do i=1,numberOfTurbines
!~         if ( forceFieldUV(i) % c .gt. 0 .or. forceFieldW(i) % c .gt. 0) then

        ! Only perform is turbine is active in this processor
        if (turbineArray(i) % operate) then
#ifdef PPLES_GPU
            if (turbineArray(i) % sampling /= 'atPoint') then
#endif
            ! Convolute force onto the domain
            call atm_clock_convolve%start()
            call atm_lesgo_convolute_force(i)
            call atm_clock_convolve%stop()
            atm_time_convolve = atm_time_convolve + atm_clock_convolve%time
#ifdef PPLES_GPU
            endif
#endif

            ! Only do this if the correction is active
            if (turbineArray(i) % tipALMCorrection .eqv. .true.)  then
#ifdef PPLES_GPU
                if (turbineArray(i) % sampling /= 'atPoint') then
#endif
                ! Compute the correction for the Cl coefficient.  Structure-off
                ! and structure-on atPoint turbines are handled above by
                ! atm_batch_cl_correction_gpu().
                call atm_clock_clcorr%start()
                call atm_compute_cl_correction(i)
                call atm_clock_clcorr%stop()
                atm_time_clcorr = atm_time_clcorr + atm_clock_clcorr%time
#ifdef PPLES_GPU
                endif
#endif
            endif

        endif

!~         ! Sync the nacelle force
!~         integrateNacelleForce=0.
!~
!~         do c=1,forceFieldUV(i) % c
!~             if (turbineArray(i) % nacelle) then
!~                 integrateNacelleForce = integrateNacelleForce +  &
!~                     forceFieldUV(i) % force(1,c) * dx *dy * dz * z_i**2*u_star**2
!~             endif
!~         enddo


    enddo


!~         totForce=0.
!~         call mpi_allreduce( integrateNacelleForce,  totForce, 1,   &
!~                              mpi_rprec, mpi_sum, comm, ierr)

       !write(*,*) 'Integrated Nacelle Force is: ', integrateNacelleForce
!~         if (coord == 0) then
!~             write(*,*) 'Integrated Total Force is: ', totForce
!~         endif
endif
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Convolute force ', myClock % time

! Phase-1 calls stop here on non-update steps too (apply/output belong to
! phase 2 so they run exactly once per step).
if (ph == 1) return

#ifdef PPLES_GPU
! fxa/fya/fza are reset every timestep. Re-deposit the latest atPoint force
! field on non-update timesteps as well, matching the CPU path's held-force
! behavior when updateInterval > 1.
if (atm_atpoint_present) then
    call atm_clock_convolve%start()
    call atm_batch_convolute_force_gpu(atm_force_update_step)
    call atm_clock_convolve%stop()
    atm_time_convolve = atm_time_convolve + atm_clock_convolve%time
endif
#endif

if (atm_force_update_step) atm_force_state_valid = .true.

    ! This will apply body forces onto the flow field if there are forces within
    ! this domain
!~  call myClock % start()
    call atm_clock_apply%start()
    call atm_lesgo_apply_force()
    call atm_clock_apply%stop()
    atm_time_apply = atm_time_apply + atm_clock_apply%time
!~     call myClock % stop()
!~     write(*,*) 'coord ', coord, '  Apply force ', myClock % time

!!! Sync the integrated forces (used for debugging)
!do i=1,numberOfTurbines
!    j=turbineArray(i) % turbineTypeID ! The turbine type ID
!    ! Sync all the integrated blade forces
!    turbineArray(i) % bladeVectorDummy=turbineArray(i) % integratedBladeForces
!    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                   &
!                       turbineArray(i) % integratedBladeForces,              &
!                       size(turbineArray(i) % bladeVectorDummy),             &
!                       mpi_rprec, mpi_sum, comm, ierr)


!    if (coord==0) then

!    do q=1, turbineArray(i) % numBladePoints
!        do n=1, turbineArray(i) % numAnnulusSections
!            do m=1, turbineModel(j) % numBl
!                write(*,*) 'blade ',m,'section ',q, 'force ratio', &
!                turbineArray(i) % integratedBladeForces(m,n,q,1) /  &
!                turbineArray(i) % bladeForces(m,n,q,1) , &
!                turbineArray(i) % integratedBladeForces(m,n,q,2) /  &
!                turbineArray(i) % bladeForces(m,n,q,2) , &
!                turbineArray(i) % integratedBladeForces(m,n,q,3) /  &
!                turbineArray(i) % bladeForces(m,n,q,3)
!            enddo
!        enddo
!    enddo
!    endif

!enddo

atm_output_step = .false.
if (outputInterval > 0) then
    atm_output_step = mod(jt_total-1, outputInterval) == 0
endif
if (atm_output_step) then
    do i=1, numberOfTurbines
        if (coord == turbineArray(i) % master) then
        !~  call myClock % start()

            call atm_clock_output%start()
            call atm_output(i, jt_total, total_time*z_i/u_star)
            call atm_clock_output%stop()
            atm_time_output = atm_time_output + atm_clock_output%time
        !~     call myClock % stop()
        !~     write(*,*) 'coord ', coord, '  Output ', myClock % time
        endif
    enddo
endif

end subroutine atm_lesgo_forcing


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Complie this subroutines only if MPI will be used
#ifdef PPMPI

subroutine atm_lesgo_mpi_gather()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Gather all host ATM state with one packed reduction per operating turbine.
implicit none

if (nproc <= 1) return
! Structure-on runs add Cm/pitchingMoment to the same packed payload.
call atm_lesgo_mpi_gather_packed()

end subroutine atm_lesgo_mpi_gather

subroutine atm_lesgo_mpi_gather_packed()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Consolidate the ATM per-turbine reductions into one allreduce to reduce
! latency on multi-GPU runs.
implicit none

integer :: i, nitem, npack, pos
real(rprec), allocatable, save :: packed_send(:), packed_recv(:)
integer, pointer :: TURBINE_COMMUNICATOR
logical :: struct_active


struct_active = atm_structure_enabled()

do i=1,numberOfTurbines

    if (turbineArray(i) % operate) then

        TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

        npack = size(turbineArray(i) % bladeForces) +                         &
                size(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:)) +      &
                size(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:)) +      &
                size(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:)) +      &
                size(turbineArray(i) % alpha) + size(turbineArray(i) % lift) +&
                size(turbineArray(i) % drag) + size(turbineArray(i) % Cl) +   &
                size(turbineArray(i) % Cd) + size(turbineArray(i) % Vmag) +   &
                size(turbineArray(i) % axialForce) +                         &
                size(turbineArray(i) % tangentialForce) +                    &
                size(turbineArray(i) % windVectors(:,:,:,1:3)) +             &
                size(turbineArray(i) % induction_a(:,:,:)) +                 &
                size(turbineArray(i) % u_infinity(:,:,:)) + 7
        if (struct_active) then
            npack = npack + size(turbineArray(i) % Cm) +                      &
                size(turbineArray(i) % pitchingMoment)
        endif

        if (allocated(packed_send)) then
            if (size(packed_send) /= npack) then
                deallocate(packed_send, packed_recv)
            endif
        endif
        if (.not. allocated(packed_send)) then
            allocate(packed_send(npack), packed_recv(npack))
        endif

        pos = 1

        nitem = size(turbineArray(i) % bladeForces)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % bladeForces,  &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % bladeAlignedVectors(:,:,:,1,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % bladeAlignedVectors(:,:,:,2,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % bladeAlignedVectors(:,:,:,3,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % alpha)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % alpha,        &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % lift)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % lift,         &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % drag)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % drag,         &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cl)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cl,           &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cd)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cd,           &
            (/ nitem /))
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Cm,       &
                (/ nitem /))
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            packed_send(pos:pos+nitem-1) = reshape(                           &
                turbineArray(i) % pitchingMoment, (/ nitem /))
            pos = pos + nitem
        endif

        nitem = size(turbineArray(i) % Vmag)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % Vmag,         &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % axialForce)
        packed_send(pos:pos+nitem-1) = reshape(turbineArray(i) % axialForce,   &
            (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % tangentialForce)
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % tangentialForce, (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % windVectors(:,:,:,1:3))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % windVectors(:,:,:,1:3), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % induction_a(:,:,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % induction_a(:,:,:), (/ nitem /))
        pos = pos + nitem

        nitem = size(turbineArray(i) % u_infinity(:,:,:))
        packed_send(pos:pos+nitem-1) = reshape(                               &
            turbineArray(i) % u_infinity(:,:,:), (/ nitem /))
        pos = pos + nitem

        packed_send(pos) = turbineArray(i) % torqueRotor
        packed_send(pos+1) = turbineArray(i) % thrust
        packed_send(pos+2:pos+4) = turbineArray(i) % nacelleForce
        packed_send(pos+5) = turbineArray(i) % VelNacelle_sampled
        packed_send(pos+6) = turbineArray(i) % VelNacelle_corrected

        call mpi_allreduce(packed_send, packed_recv, npack, mpi_rprec,         &
                           mpi_sum, TURBINE_COMMUNICATOR, ierr)

        pos = 1

        nitem = size(turbineArray(i) % bladeForces)
        turbineArray(i) % bladeForces = reshape(packed_recv(pos:pos+nitem-1),  &
            shape(turbineArray(i) % bladeForces))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:))
        turbineArray(i) % bladeAlignedVectors(:,:,:,1,:) = reshape(           &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % bladeAlignedVectors(:,:,:,1,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:))
        turbineArray(i) % bladeAlignedVectors(:,:,:,2,:) = reshape(           &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % bladeAlignedVectors(:,:,:,2,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:))
        turbineArray(i) % bladeAlignedVectors(:,:,:,3,:) = reshape(           &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % bladeAlignedVectors(:,:,:,3,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % alpha)
        turbineArray(i) % alpha = reshape(packed_recv(pos:pos+nitem-1),        &
            shape(turbineArray(i) % alpha))
        pos = pos + nitem

        nitem = size(turbineArray(i) % lift)
        turbineArray(i) % lift = reshape(packed_recv(pos:pos+nitem-1),         &
            shape(turbineArray(i) % lift))
        pos = pos + nitem

        nitem = size(turbineArray(i) % drag)
        turbineArray(i) % drag = reshape(packed_recv(pos:pos+nitem-1),         &
            shape(turbineArray(i) % drag))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cl)
        turbineArray(i) % Cl = reshape(packed_recv(pos:pos+nitem-1),           &
            shape(turbineArray(i) % Cl))
        pos = pos + nitem

        nitem = size(turbineArray(i) % Cd)
        turbineArray(i) % Cd = reshape(packed_recv(pos:pos+nitem-1),           &
            shape(turbineArray(i) % Cd))
        pos = pos + nitem

        if (struct_active) then
            nitem = size(turbineArray(i) % Cm)
            turbineArray(i) % Cm = reshape(packed_recv(pos:pos+nitem-1),       &
                shape(turbineArray(i) % Cm))
            pos = pos + nitem

            nitem = size(turbineArray(i) % pitchingMoment)
            turbineArray(i) % pitchingMoment = reshape(                       &
                packed_recv(pos:pos+nitem-1),                                 &
                shape(turbineArray(i) % pitchingMoment))
            pos = pos + nitem
        endif

        nitem = size(turbineArray(i) % Vmag)
        turbineArray(i) % Vmag = reshape(packed_recv(pos:pos+nitem-1),         &
            shape(turbineArray(i) % Vmag))
        pos = pos + nitem

        nitem = size(turbineArray(i) % axialForce)
        turbineArray(i) % axialForce = reshape(                               &
            packed_recv(pos:pos+nitem-1), shape(turbineArray(i) % axialForce))
        pos = pos + nitem

        nitem = size(turbineArray(i) % tangentialForce)
        turbineArray(i) % tangentialForce = reshape(                          &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % tangentialForce))
        pos = pos + nitem

        nitem = size(turbineArray(i) % windVectors(:,:,:,1:3))
        turbineArray(i) % windVectors(:,:,:,1:3) = reshape(                   &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % windVectors(:,:,:,1:3)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % induction_a(:,:,:))
        turbineArray(i) % induction_a(:,:,:) = reshape(                       &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % induction_a(:,:,:)))
        pos = pos + nitem

        nitem = size(turbineArray(i) % u_infinity(:,:,:))
        turbineArray(i) % u_infinity(:,:,:) = reshape(                        &
            packed_recv(pos:pos+nitem-1),                                     &
            shape(turbineArray(i) % u_infinity(:,:,:)))
        pos = pos + nitem

        turbineArray(i) % torqueRotor = packed_recv(pos)
        turbineArray(i) % thrust = packed_recv(pos+1)
        turbineArray(i) % nacelleForce = packed_recv(pos+2:pos+4)
        turbineArray(i) % VelNacelle_sampled = packed_recv(pos+5)
        turbineArray(i) % VelNacelle_corrected = packed_recv(pos+6)

    endif

enddo

end subroutine atm_lesgo_mpi_gather_packed

#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will feed the velocity at all the actuator points into the atm
! This is done by using trilinear interpolation from lesgo
! Force will be calculated based on the velocities and stored on forceField
implicit none

integer, intent(in) :: i ! The turbine number
integer :: m,n,q,j
! mpi_velocity only used for Spalart method
real(rprec), dimension(3) :: velocity, mpi_velocity
real(rprec), dimension(3) :: xyz    ! Point onto which to interpolate velocity
real(rprec), pointer, dimension(:) :: x,y,z,zw
#ifdef PPLES_GPU
integer :: mm, nn, qq, p0
#endif

! The MPI turbine communcator
integer, pointer :: TURBINE_COMM

TURBINE_COMM => turbineArray(i) % TURBINE_COMM_WORLD

j=turbineArray(i) % turbineTypeID ! The turbine type ID


! Declare x, y, and z as pointers to the grid variables x, y, and z (LESGO)
nullify(x,y,z,zw)
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

if (turbineArray(i) % sampling == 'Spalart') then
    ! This loop goes through all the blade points and calculates the respective
    ! body forces then imposes it onto the force field
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl

                ! Actuator point onto which to interpolate the velocity
                xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

                velocity = 0._rprec
                mpi_velocity = 0._rprec

                call atm_lesgo_compute_spalart_u(i, xyz, velocity)

                mpi_velocity = velocity

                ! Complie this subroutines only if MPI will be used
#ifdef PPMPI
!~                     call mpi_barrier( TURBINE_COMM, ierr )
                    ! Sync all the blade forces
                    call mpi_allreduce(mpi_velocity, velocity, size(velocity), &
                           mpi_rprec, mpi_sum, TURBINE_COMM , ierr)
#endif

                ! This will compute the blade force for the specific point
                if (  z(1) <= xyz(3)/z_i .and. xyz(3)/z_i < z(nz) ) then
                    call atm_computeBladeForce(i,m,n,q,velocity)
                else
                    velocity = 0._rprec
                endif

            enddo
        enddo
    enddo


else if (turbineArray(i) % sampling == 'atPoint') then
    ! This loop goes through all the blade points and calculates the respective
    ! body forces then imposes it onto the force field
#ifdef PPLES_GPU
    ! Consume velocity sampled from resident u/v/w by the batched GPU kernel.
    ! The airfoil force model
    ! (atm_computeBladeForce) still runs on the host using the sampled velocity.
    mm = turbineModel(j) % numBl
    nn = turbineArray(i) % numAnnulusSections
    qq = turbineArray(i) % numBladePoints
    if (.not. atm_batch_sampled) stop 'ATM batched velocity sample is missing'
    do q = 1, qq
    do n = 1, nn
    do m = 1, mm
        p0 = atm_bp_off(i) + ((m-1)*nn + (n-1))*qq + q
        if (atm_inr_all(p0) == 1)                                              &
            call atm_computeBladeForce(i, m, n, q, atm_velbp_all(1:3,p0))
    end do
    end do
    end do
#else
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl

                ! Actuator point onto which to interpolate the velocity
                xyz=turbineArray(i) % bladePoints(m,n,q,1:3)

                ! Non-dimensionalizes the point location
                xyz=xyz/z_i

                ! Interpolate velocities if inside the domain
                if (  z(1) <= xyz(3) .and. xyz(3) < z(nz) ) then
                    velocity(1)=                                               &
                    trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                    velocity(2)=                                               &
                    trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
                    velocity(3)=                                               &
                    trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

                    ! This will compute the blade force for the specific point
                    call atm_computeBladeForce(i,m,n,q,velocity)

                endif

            enddo
        enddo
    enddo
#endif
endif

    ! Calculate Nacelle force
    if (turbineArray(i) % nacelle) then
        xyz=turbineArray(i) % nacelleLocation
        xyz=xyz/z_i
        if (  z(1) <= xyz(3) .and. xyz(3) < z(nz) ) then

            velocity(1)=                                                   &
            trilinear_interp(u(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
            velocity(2)=                                                   &
            trilinear_interp(v(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star
            velocity(3)=                                                   &
            trilinear_interp(w_uv(1:nx,1:ny,lbz:nz),lbz,xyz)*u_star

            call atm_computeNacelleForce(i,velocity)

        endif
    endif

end subroutine atm_lesgo_force

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_compute_Spalart_u(i, xyz, velocity)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will calculate the sampling velocity using the proposed method
! from Spalart
! n turbine number
! xyz actuator point position vector
! velocity reference velocity for computing lift and drag

implicit none

integer, intent(in) :: i
real(rprec), intent(in) :: xyz(3)
real(rprec), intent(inout) :: velocity(3)

integer :: c, m, n, q

! Pointers for mesh
real(rprec), pointer, dimension(:) :: z,zw

! Test for time optimization
real(rprec) :: dist, a(3), projectradius, epsilon
real(rprec) :: epsilon_sq, kernel_norm

nullify(z,zw)
z => grid % z
zw => grid % zw

! Value of epsilon
epsilon=turbineArray(i) % epsilon
epsilon_sq = epsilon * epsilon
kernel_norm = 1._rprec / ((epsilon * epsilon_sq) * (pi * sqrt(pi)))

! Projection radius
projectradius = turbineArray(i) % projectionRadius

! Set the velocity to zero
velocity = 0._rprec


do c=1,forceFieldUV(i) % c

    a = forceFieldUV(i) %  location(1:3, c)
    m = forceFieldUV(i) %  ijk(1, c)
    n = forceFieldUV(i) %  ijk(2, c)
    q = forceFieldUV(i) %  ijk(3, c)

    dist = sqrt((a(1)-xyz(1))**2 + (a(2)-xyz(2))**2 + (a(3)-xyz(3))**2)
    if (dist .le. projectradius * z_i) then
        if ( z(1) <= a(3)/z_i .and. a(3)/z_i < z(nz)) then

        ! The value of the kernel. This is the actual smoothing function
        velocity(1) = velocity(1) + u(m,n,q) * exp(-dist*dist/epsilon_sq)  &
                                 * kernel_norm
        velocity(2) = velocity(2) + v(m,n,q) * exp(-dist*dist/epsilon_sq)  &
                                 * kernel_norm
        endif
    endif
enddo

do c=1,forceFieldW(i) % c
    a = forceFieldW(i) %  location(1:3, c)
    m = forceFieldW(i) %  ijk(1, c)
    n = forceFieldW(i) %  ijk(2, c)
    q = forceFieldW(i) %  ijk(3, c)

    dist = sqrt((a(1)-xyz(1))**2 + (a(2)-xyz(2))**2 + (a(3)-xyz(3))**2)

    if (dist .le. projectradius) then
        if ( z(1) <= a(3)/z_i .and. a(3)/z_i < z(nz)) then

        ! The value of the kernel. This is the actual smoothing function
        velocity(3) = velocity(3) + w(m,n,q) * exp(-dist*dist/epsilon_sq)  &
                                 * kernel_norm
        endif
    endif
enddo

velocity = velocity * u_star * z_i * dx * z_i * dy *z_i * dz

end subroutine atm_lesgo_compute_Spalart_u

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_convolute_force(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will convolute the forces for each turbine

implicit none

!~ type(clock_t) :: myClock

integer, intent(in) :: i
integer :: j, m, n, q, c,mmend,nnend,qqend

integer :: ii, jj, kk  ! Indices for lesgo fields

! Test for time optimization
real(rprec) :: dist,a(3),b(3),projectradius,epsilon,const1,const2,const3
real(rprec) :: epsilon_sq, nacelle_epsilon_sq, nacelle_kernel_norm
real(rprec) :: nacelleEpsilon

! Variables for convolution force
real(rprec) :: kernel, force(3)

! Pointers for the turbineArray quantities
real(rprec), pointer, dimension(:,:,:,:) :: bladeForces, bladePoints

real(rprec), pointer, dimension(:,:) :: bodyForceUV, bodyForceW

nullify(bladeForces)
nullify(bladePoints)
nullify(bodyForceUV)
nullify(bodyForceW)

bladeForces => turbineArray(i) % bladeForces
bladePoints => turbineArray(i) % bladePoints

bodyForceUV => forceFieldUV(i) % force
bodyForceW =>  forceFieldW(i) % force

!real(rprec) :: dummyForce(3)  ! Debugging

j=turbineArray(i) % turbineTypeID ! The turbine type ID

! This will convolute the blade force onto the grid points
! affected by the turbines on both grids
! Only if the distance is less than specified value
mmend=turbineModel(j) % numBl
nnend=turbineArray(i) % numAnnulusSections
qqend=turbineArray(i) % numBladePoints
projectradius=turbineArray(i) % projectionRadius
epsilon=turbineArray(i) % epsilon
epsilon_sq = epsilon * epsilon
const1 = 1._rprec / ((epsilon * epsilon_sq) * (pi * sqrt(pi)))
const2 = z_i / (u_star*u_star)
const3=const1*const2
nacelle_epsilon_sq = 1._rprec
nacelle_kernel_norm = 0._rprec
if (turbineArray(i) % nacelle) then
    nacelleEpsilon = turbineArray(i) % nacelleEpsilon
    nacelle_epsilon_sq = nacelleEpsilon * nacelleEpsilon
    nacelle_kernel_norm = 1._rprec /                                        &
        ((nacelleEpsilon * nacelle_epsilon_sq) * (pi * sqrt(pi)))
endif

! Body Force implementation using velocity sampling at the actuator point
if (turbineArray(i) % sampling == 'atPoint') then

    !~  call myClock % start()
    do c=1,forceFieldUV(i) % c
        a= forceFieldUV(i) %  location(1:3,c)
        force=0._rprec

        ! Blade forces
        do m=1, mmend
            do n=1, nnend
               do q=1, qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2            &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                    ! The value of the kernel. This is the actual smoothing function
                     force(1:2) = force(1:2) + bladeForces(m,n,q,1:2)     &
                                  * exp(-dist*dist/epsilon_sq)
                    endif

                enddo
            enddo
        enddo
        force(1:2)=force(1:2)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
    !~         if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                &
                         * nacelle_kernel_norm
                !write(*,*) 'kernel Value= ', kernel
                force(1:2) = force(1:2)+turbineArray(i) % nacelleForce(1:2) *  &
                             kernel *const2
    !~          integrateNacelleForce=integrateNacelleForce+force(1) * dx *dy * dz * z_i**3

    !~         endif
        endif


        bodyForceUV(1:2,c) = force(1:2)
    !~     if (abs(bodyForceUV(1,c)) .gt. 0) then
    !~                 write(*,*) 'bodyForceUV is: ', bodyForceUV(1,c)
    !~     endif
    enddo


    do c=1,forceFieldW(i) % c
        a= forceFieldW(i) %  location(1:3,c)
        force=0._rprec

        ! Blade forces
        do m=1,mmend
            do n=1,nnend
               do q=1,qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2            &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                    ! The value of the kernel. This is the actual smoothing function
                    force(3) = force(3) +  bladeForces(m,n,q,3) &
                               * exp(-dist*dist/epsilon_sq)
                    endif

                enddo
            enddo
        enddo
        force(3)=force(3)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
            if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                &
                         * nacelle_kernel_norm
                force(3) = force(3)+turbineArray(i) % nacelleForce(3) *           &
                           kernel *const2
            endif
        endif

        bodyForceW(3,c) = force(3)
    enddo

! The Spalart method uses the local velocity field.
! For this reason it needs to be done explicitly in this module
! and cannot be generally coded from the actuator_turbine_model module
elseif (turbineArray(i) % sampling == 'Spalart') then

    !~  call myClock % start()
    do c=1,forceFieldUV(i) % c
        a= forceFieldUV(i) %  location(1:3,c)
        force=0._rprec
        ! Indices for velocity field
        ii = forceFieldUV(i) % ijk(1,c)
        jj = forceFieldUV(i) % ijk(2,c)
        kk = forceFieldUV(i) % ijk(3,c)

        ! Blade forces
        do m=1, mmend
            do n=1, nnend
               do q=1, qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2            &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                        ! The value of the kernel.
                        ! This is the actual smoothing function
                        ! Divide by velocity magnitude
                         force(1) = force(1) +  bladeForces(m,n,q,1) *         &
                                      exp(-dist*dist/epsilon_sq)              &
                         / (turbineArray(i) % Vmag(m,n,q)) *                     &
                         ( u(ii,jj,kk)  * u_star +                             &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))

                         force(2) = force(2) +  bladeForces(m,n,q,2) *         &
                                      exp(-dist*dist/epsilon_sq)              &
                         / turbineArray(i) % Vmag(m,n,q) *                     &
                         ( v(ii,jj,kk) * u_star +                              &
                         turbineArray(i) % bladeAlignedVectors(m,n,q,2,2) *    &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))

                    endif

                enddo
            enddo
        enddo
        force(1:2)=force(1:2)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
    !~         if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                   &
                         * nacelle_kernel_norm
                !write(*,*) 'kernel Value= ', kernel
                force(1:2) = force(1:2)+turbineArray(i) % nacelleForce(1:2) *  &
                             kernel *const2
    !~          integrateNacelleForce=integrateNacelleForce+force(1) * dx *dy * dz * z_i**3

    !~         endif
        endif


        bodyForceUV(1:2,c) = force(1:2)
    enddo


    do c=1,forceFieldW(i) % c
        a= forceFieldW(i) %  location(1:3,c)
        force=0._rprec
        ! Indices for velocity field
        ii = forceFieldW(i) % ijk(1,c)
        jj = forceFieldW(i) % ijk(2,c)
        kk = forceFieldW(i) % ijk(3,c)

        ! Blade forces
        do m=1,mmend
            do n=1,nnend
               do q=1,qqend

                    b= bladePoints(m,n,q,:)
                    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2               &
                               + (a(3)-b(3))**2)

                    if (dist .le. projectradius) then
                        ! The value of the kernel.
                        ! This is the actual smoothing function
                        force(3) = force(3) +  bladeForces(m,n,q,3) *          &
                                   exp(-dist*dist/epsilon_sq)                 &
                         / turbineArray(i) % Vmag(m,n,q) *                     &
                         ( w(ii,jj,kk) * u_star +                              &
                         turbineArray(i) % bladeAlignedVectors(m,n,q,2,3) *    &
                         turbineArray(i) % rotSpeed *                          &
                         turbineArray(i) % bladeRadius(m,n,q) *                &
                         cos(turbineModel(j) % PreCone))
                    endif

                enddo
            enddo
        enddo
        force(3)=force(3)* const3

        ! Nacelle force
        if (turbineArray(i) % nacelle) then
            b=turbineArray(i) % nacelleLocation
            dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
            if (dist .le. projectradius) then
                ! The value of the kernel. This is the actual smoothing function
                kernel = exp(-dist*dist/nacelle_epsilon_sq)                  &
                         * nacelle_kernel_norm
                force(3) = force(3)+turbineArray(i) % nacelleForce(3) *           &
                           kernel *const2
            endif
        endif

        bodyForceW(3,c) = force(3)
    enddo

endif

end subroutine atm_lesgo_convolute_force

#ifdef PPLES_GPU
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_atpoint_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Build the static device tables for the batched atPoint path: prefix offsets,
! concatenated force-field cell lists (location/ijk/turbine-id - static for
! fixed turbines), grid axes + autowrap, and the per-step staging buffers.
! Called once, lazily, on the first force step.
use grid_m, only : grid
use param,  only : nx, ny, nz, lbz
implicit none
integer :: i, j, c, p, nbp

if (atm_batch_ready) return

allocate(atm_bp_off(numberOfTurbines+1))
atm_bp_off(1) = 0
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    nbp = turbineModel(j) % numBl * turbineArray(i) % numAnnulusSections       &
                                  * turbineArray(i) % numBladePoints
    atm_bp_off(i+1) = atm_bp_off(i) + nbp
end do
atm_nbp_tot = atm_bp_off(numberOfTurbines+1)

atm_cUV_tot = 0
atm_cW_tot  = 0
do i = 1, numberOfTurbines
    atm_cUV_tot = atm_cUV_tot + forceFieldUV(i) % c
    atm_cW_tot  = atm_cW_tot  + forceFieldW(i)  % c
end do

allocate(atm_bp_all(3, max(atm_nbp_tot,1)), atm_bf_all(3, max(atm_nbp_tot,1)))
allocate(atm_velbp_all(3, max(atm_nbp_tot,1)), atm_inr_all(max(atm_nbp_tot,1)))
allocate(atm_tconst(ATM_NTC, numberOfTurbines))
allocate(atm_sampling_mode(numberOfTurbines))
atm_bp_all = 0._rprec; atm_bf_all = 0._rprec
atm_velbp_all = 0._rprec; atm_inr_all = 0; atm_tconst = 0._rprec
do i = 1, numberOfTurbines
    if (turbineArray(i) % sampling == 'atPoint') then
        atm_sampling_mode(i) = ATM_SAMPLING_ATPOINT
    else
        atm_sampling_mode(i) = ATM_SAMPLING_SPALART
    endif
enddo

allocate(atm_locUV_all(3, max(atm_cUV_tot,1)),                                 &
         atm_ijkUV_all(3, max(atm_cUV_tot,1)), atm_tidUV(max(atm_cUV_tot,1)))
allocate(atm_locW_all(3, max(atm_cW_tot,1)),                                   &
         atm_ijkW_all(3, max(atm_cW_tot,1)), atm_tidW(max(atm_cW_tot,1)))
atm_tidUV = 0; atm_tidW = 0
p = 0
do i = 1, numberOfTurbines
    do c = 1, forceFieldUV(i) % c
        p = p + 1
        atm_locUV_all(1:3,p) = forceFieldUV(i) % location(1:3,c)
        atm_ijkUV_all(1:3,p) = forceFieldUV(i) % ijk(1:3,c)
        atm_tidUV(p) = i
    end do
end do
p = 0
do i = 1, numberOfTurbines
    do c = 1, forceFieldW(i) % c
        p = p + 1
        atm_locW_all(1:3,p) = forceFieldW(i) % location(1:3,c)
        atm_ijkW_all(1:3,p) = forceFieldW(i) % ijk(1:3,c)
        atm_tidW(p) = i
    end do
end do

allocate(atm_gx(nx), atm_gy(ny), atm_gz(lbz:nz))
allocate(atm_awi(0:nx+1), atm_awj(0:ny+1))
atm_gx(1:nx)    = grid % x(1:nx)
atm_gy(1:ny)    = grid % y(1:ny)
atm_gz(lbz:nz)  = grid % z(lbz:nz)
atm_awi(0:nx+1) = grid % autowrap_i(0:nx+1)
atm_awj(0:ny+1) = grid % autowrap_j(0:ny+1)

!$acc enter data copyin(atm_bp_off, atm_bp_all, atm_bf_all,                    &
!$acc                   atm_velbp_all, atm_inr_all, atm_tconst,                &
!$acc                   atm_sampling_mode,                                     &
!$acc                   atm_locUV_all, atm_ijkUV_all, atm_tidUV,               &
!$acc                   atm_locW_all, atm_ijkW_all, atm_tidW,                  &
!$acc                   atm_gx, atm_gy, atm_gz, atm_awi, atm_awj)

if (atm_spalart_present) then
    allocate(atm_spalart_forceUV(2, max(atm_cUV_tot,1)))
    allocate(atm_spalart_forceW(max(atm_cW_tot,1)))
    atm_spalart_forceUV = 0._rprec
    atm_spalart_forceW = 0._rprec
    !$acc enter data create(atm_spalart_forceUV, atm_spalart_forceW)
endif

atm_batch_ready = .true.
end subroutine atm_batch_atpoint_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_sample_velocity_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched device velocity sampling at the actuator points of ALL turbines in
! one kernel. Identical per-point arithmetic to atm_sample_velocity_atpoint_gpu
! (cell_indx thresholds, autowrap, 8-corner trilinear of u, v, w_uv x u_star);
! only the per-turbine data regions/kernels/syncs are collapsed into one
! update device + one kernel + one update self per force step.
use sim_param, only : u, v, w
use param,     only : nx, ny, nz, lbz, dx, dy, dz, L_x, L_y, L_z, u_star, z_i
implicit none
integer :: i, j, m, n, q, p
integer :: is, js, ks, is1, js1, ks1
real(rprec) :: px, py, pz, xd, yd, zd
real(rprec) :: a1, a2, a3, a4, a5, a6, a7, a8
real(rprec), parameter :: thr = 1.e-9_rprec

call atm_batch_atpoint_init()
atm_batch_sampled = .false.
if (atm_nbp_tot == 0) return

! Flatten the (rotated) blade points: m-outer / q-fastest order -- the same
! order the convolution kernel sums over.
do i = 1, numberOfTurbines
    j = turbineArray(i) % turbineTypeID
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, turbineArray(i) % numAnnulusSections
    do q = 1, turbineArray(i) % numBladePoints
        p = p + 1
        atm_bp_all(1,p) = turbineArray(i) % bladePoints(m,n,q,1)
        atm_bp_all(2,p) = turbineArray(i) % bladePoints(m,n,q,2)
        atm_bp_all(3,p) = turbineArray(i) % bladePoints(m,n,q,3)
    end do
    end do
    end do
end do
!$acc update device(atm_bp_all)

!$acc parallel loop gang vector default(present)                              &
!$acc     private(px,py,pz,is,js,ks,is1,js1,ks1,xd,yd,zd,                    &
!$acc             a1,a2,a3,a4,a5,a6,a7,a8)
do p = 1, atm_nbp_tot
        px = atm_bp_all(1,p) / z_i
        py = atm_bp_all(2,p) / z_i
        pz = atm_bp_all(3,p) / z_i
        if (atm_gz(1) <= pz .and. pz < atm_gz(nz)) then
            atm_inr_all(p) = 1
            px = modulo(px, L_x)
            if (abs(px)/L_x < thr) then
                is = 1
            else if (abs(px-L_x)/L_x < thr) then
                is = nx
            else
                is = floor(px/dx) + 1
            end if
            py = modulo(py, L_y)
            if (abs(py)/L_y < thr) then
                js = 1
            else if (abs(py-L_y)/L_y < thr) then
                js = ny
            else
                js = floor(py/dy) + 1
            end if
            if (abs(pz - atm_gz(nz))/L_z < thr) then
                ks = nz - 1
            else
                ks = floor((pz - atm_gz(1))/dz) + 1
            end if
            is1 = atm_awi(is+1); js1 = atm_awj(js+1); ks1 = ks + 1
            xd  = px - atm_gx(is); yd = py - atm_gy(js); zd = pz - atm_gz(ks)
            a1 = u(is,js,ks)   + xd*(u(is1,js,ks)   - u(is,js,ks))  /dx
            a2 = u(is,js1,ks)  + xd*(u(is1,js1,ks)  - u(is,js1,ks)) /dx
            a3 = u(is,js,ks1)  + xd*(u(is1,js,ks1)  - u(is,js,ks1)) /dx
            a4 = u(is,js1,ks1) + xd*(u(is1,js1,ks1) - u(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(1,p) = (a5 + zd*(a6-a5)/dz) * u_star
            a1 = v(is,js,ks)   + xd*(v(is1,js,ks)   - v(is,js,ks))  /dx
            a2 = v(is,js1,ks)  + xd*(v(is1,js1,ks)  - v(is,js1,ks)) /dx
            a3 = v(is,js,ks1)  + xd*(v(is1,js,ks1)  - v(is,js,ks1)) /dx
            a4 = v(is,js1,ks1) + xd*(v(is1,js1,ks1) - v(is,js1,ks1))/dx
            a5 = a1 + yd*(a2-a1)/dy
            a6 = a3 + yd*(a4-a3)/dy
            atm_velbp_all(2,p) = (a5 + zd*(a6-a5)/dz) * u_star
            a1 = 0.5_rprec * (w(is,js,ks) + w(is,js,ks+1))
            a2 = 0.5_rprec * (w(is1,js,ks) + w(is1,js,ks+1))
            a3 = 0.5_rprec * (w(is,js1,ks) + w(is,js1,ks+1))
            a4 = 0.5_rprec * (w(is1,js1,ks) + w(is1,js1,ks+1))
            a5 = a1 + xd*(a2-a1)/dx
            a6 = a3 + xd*(a4-a3)/dx
            a7 = a5 + yd*(a6-a5)/dy
            if (ks1 == nz .and. coord < nproc - 1) then
                a1 = atm_wuv_recv_up(is,js)
                a2 = atm_wuv_recv_up(is1,js)
                a3 = atm_wuv_recv_up(is,js1)
                a4 = atm_wuv_recv_up(is1,js1)
            else if (ks1 == nz) then
                a1 = 0.5_rprec * (w(is,js,nz-1) + w(is,js,nz))
                a2 = 0.5_rprec * (w(is1,js,nz-1) + w(is1,js,nz))
                a3 = 0.5_rprec * (w(is,js1,nz-1) + w(is,js1,nz))
                a4 = 0.5_rprec * (w(is1,js1,nz-1) + w(is1,js1,nz))
            else
                a1 = 0.5_rprec * (w(is,js,ks1) + w(is,js,ks1+1))
                a2 = 0.5_rprec * (w(is1,js,ks1) + w(is1,js,ks1+1))
                a3 = 0.5_rprec * (w(is,js1,ks1) + w(is,js1,ks1+1))
                a4 = 0.5_rprec * (w(is1,js1,ks1) + w(is1,js1,ks1+1))
            end if
            a5 = a1 + xd*(a2-a1)/dx
            a6 = a3 + xd*(a4-a3)/dx
            a8 = a5 + yd*(a6-a5)/dy
            atm_velbp_all(3,p) = (a7 + zd*(a8-a7)/dz) * u_star
        else
            atm_inr_all(p)     = 0
            atm_velbp_all(1,p) = 0._rprec
            atm_velbp_all(2,p) = 0._rprec
            atm_velbp_all(3,p) = 0._rprec
        end if
end do

!$acc update self(atm_velbp_all, atm_inr_all)
atm_batch_sampled = .true.

end subroutine atm_batch_sample_velocity_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_convolute_force_gpu(refresh_state)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched atPoint Gaussian force convolution + apply for ALL turbines: one
! kernel per grid (UV, W) over the concatenated force-field cell lists,
! scattering straight into the device fxa/fya/fza. The blade-point summation
! order is unchanged, but squared-distance arithmetic avoids a redundant sqrt
! before the Gaussian radius test and exponent. The scatter uses atomics
! because cells of DIFFERENT turbines may coincide where
! projection regions overlap (within a turbine each cell is unique).
! The forceField%force host write-back of the per-turbine path is unnecessary:
! production convolution scatters directly into the resident LES force arrays.
use sim_param, only : fxa, fya, fza
use param,     only : z_i, u_star
implicit none
logical, intent(in) :: refresh_state
integer :: i, j, c, p, m, n, q
integer :: ii, jj, kk
real(rprec) :: eps, dist_sq, kw, fx1, fx2, fx3
real(rprec) :: pr, prsq, epsq, cc2, cc3

call atm_batch_atpoint_init()

if (refresh_state) then
    ! Cache per-turbine constants + flattened blade forces after the gather.
    ! On intervening updateInterval timesteps, reuse this state exactly as the
    ! CPU path reuses its previously convolved forceField.
    do i = 1, numberOfTurbines
        if (turbineArray(i) % operate .and.                                    &
            turbineArray(i) % sampling == 'atPoint') then
            atm_tconst(14,i) = 1._rprec
        else
            atm_tconst(14,i) = 0._rprec
        end if
        ! Keep the same const1*const2 operation order as the per-turbine path.
        ! The convolution kernels below use squared distance directly, so final
        ! fields are roundoff-equivalent rather than byte-identical.
        eps = turbineArray(i) % epsilon
        atm_tconst(1,i) = turbineArray(i) % projectionRadius
        atm_tconst(2,i) = eps * eps
        atm_tconst(4,i) = z_i / (u_star * u_star)
        atm_tconst(3,i) = (1._rprec /                                          &
            ((eps * atm_tconst(2,i)) * (pi * sqrt(pi)))) * atm_tconst(4,i)
        if (turbineArray(i) % nacelle) then
            atm_tconst(13,i) = 1._rprec
            eps = turbineArray(i) % nacelleEpsilon
            atm_tconst(5,i) = eps * eps
            atm_tconst(6,i) = 1._rprec /                                      &
                ((eps * eps * eps) * (pi * sqrt(pi)))
            atm_tconst(7:9,i) = turbineArray(i) % nacelleLocation(1:3)
            atm_tconst(10:12,i) = turbineArray(i) % nacelleForce(1:3)
        else
            atm_tconst(13,i) = 0._rprec
            atm_tconst(5,i) = 1._rprec
            atm_tconst(6,i) = 0._rprec
            atm_tconst(7:12,i) = 0._rprec
        end if

        j = turbineArray(i) % turbineTypeID
        p = atm_bp_off(i)
        do m = 1, turbineModel(j) % numBl
        do n = 1, turbineArray(i) % numAnnulusSections
        do q = 1, turbineArray(i) % numBladePoints
            p = p + 1
            atm_bf_all(1,p) = turbineArray(i) % bladeForces(m,n,q,1)
            atm_bf_all(2,p) = turbineArray(i) % bladeForces(m,n,q,2)
            atm_bf_all(3,p) = turbineArray(i) % bladeForces(m,n,q,3)
        end do
        end do
        end do
    end do
    !$acc update device(atm_bf_all, atm_tconst)
endif

! ---- UV grid (force components 1,2) ----
if (atm_cUV_tot > 0) then
    !$acc parallel loop gang vector default(present)                           &
    !$acc     private(i, pr, prsq, epsq, cc2, cc3, fx1, fx2, dist_sq, kw, ii, jj, kk, p)
    do c = 1, atm_cUV_tot
        i = atm_tidUV(c)
        if (atm_tconst(14,i) > 0.5_rprec) then
            pr   = atm_tconst(1,i)
            epsq = atm_tconst(2,i)
            prsq = pr * pr
            cc3  = atm_tconst(3,i)
            cc2  = atm_tconst(4,i)
            fx1 = 0._rprec
            fx2 = 0._rprec
            !$acc loop seq
            do p = atm_bp_off(i)+1, atm_bp_off(i+1)
                dist_sq = (atm_locUV_all(1,c)-atm_bp_all(1,p))**2             &
                        + (atm_locUV_all(2,c)-atm_bp_all(2,p))**2             &
                        + (atm_locUV_all(3,c)-atm_bp_all(3,p))**2
                if (dist_sq <= prsq) then
                    kw  = exp(-dist_sq/epsq)
                    fx1 = fx1 + atm_bf_all(1,p) * kw
                    fx2 = fx2 + atm_bf_all(2,p) * kw
                end if
            end do
            fx1 = fx1 * cc3
            fx2 = fx2 * cc3
            if (atm_tconst(13,i) > 0.5_rprec) then
                dist_sq = (atm_locUV_all(1,c)-atm_tconst(7,i))**2             &
                        + (atm_locUV_all(2,c)-atm_tconst(8,i))**2             &
                        + (atm_locUV_all(3,c)-atm_tconst(9,i))**2
                kw  = exp(-dist_sq/atm_tconst(5,i)) * atm_tconst(6,i)
                fx1 = fx1 + atm_tconst(10,i) * kw * cc2
                fx2 = fx2 + atm_tconst(11,i) * kw * cc2
            end if
            ii = atm_ijkUV_all(1,c); jj = atm_ijkUV_all(2,c); kk = atm_ijkUV_all(3,c)
            !$acc atomic update
            fxa(ii,jj,kk) = fxa(ii,jj,kk) + fx1
            !$acc atomic update
            fya(ii,jj,kk) = fya(ii,jj,kk) + fx2
        end if
    end do
end if

! ---- W grid (force component 3) ----
if (atm_cW_tot > 0) then
    !$acc parallel loop gang vector default(present)                           &
    !$acc     private(i, pr, prsq, epsq, cc2, cc3, fx3, dist_sq, kw, ii, jj, kk, p)
    do c = 1, atm_cW_tot
        i = atm_tidW(c)
        if (atm_tconst(14,i) > 0.5_rprec) then
            pr   = atm_tconst(1,i)
            epsq = atm_tconst(2,i)
            prsq = pr * pr
            cc3  = atm_tconst(3,i)
            cc2  = atm_tconst(4,i)
            fx3 = 0._rprec
            !$acc loop seq
            do p = atm_bp_off(i)+1, atm_bp_off(i+1)
                dist_sq = (atm_locW_all(1,c)-atm_bp_all(1,p))**2              &
                        + (atm_locW_all(2,c)-atm_bp_all(2,p))**2              &
                        + (atm_locW_all(3,c)-atm_bp_all(3,p))**2
                if (dist_sq <= prsq) then
                    fx3 = fx3 + atm_bf_all(3,p) * exp(-dist_sq/epsq)
                end if
            end do
            fx3 = fx3 * cc3
            if (atm_tconst(13,i) > 0.5_rprec) then
                dist_sq = (atm_locW_all(1,c)-atm_tconst(7,i))**2              &
                        + (atm_locW_all(2,c)-atm_tconst(8,i))**2              &
                        + (atm_locW_all(3,c)-atm_tconst(9,i))**2
                if (dist_sq <= prsq) then
                    kw  = exp(-dist_sq/atm_tconst(5,i)) * atm_tconst(6,i)
                    fx3 = fx3 + atm_tconst(12,i) * kw * cc2
                end if
            end if
            ii = atm_ijkW_all(1,c); jj = atm_ijkW_all(2,c); kk = atm_ijkW_all(3,c)
            !$acc atomic update
            fza(ii,jj,kk) = fza(ii,jj,kk) + fx3
        end if
    end do
end if

end subroutine atm_batch_convolute_force_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_apply_spalart_force_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Spalart sampling and convolution remain host algorithms. Pack only their
! completed per-cell forces and add them to the resident LES force arrays.
! This preserves mixed atPoint/Spalart farms without a full-domain force-field
! transfer that would be both expensive and capable of overwriting atPoint
! contributions already deposited on the device.
use sim_param, only : fxa, fya, fza
implicit none

integer :: i, c, p, ii, jj, kk

if (.not. atm_spalart_present) return
call atm_batch_atpoint_init()

p = 0
do i = 1, numberOfTurbines
    do c = 1, forceFieldUV(i) % c
        p = p + 1
        if (turbineArray(i) % operate .and.                                    &
            turbineArray(i) % sampling == 'Spalart') then
            atm_spalart_forceUV(1:2,p) = forceFieldUV(i) % force(1:2,c)
        else
            atm_spalart_forceUV(1:2,p) = 0._rprec
        endif
    enddo
enddo

p = 0
do i = 1, numberOfTurbines
    do c = 1, forceFieldW(i) % c
        p = p + 1
        if (turbineArray(i) % operate .and.                                    &
            turbineArray(i) % sampling == 'Spalart') then
            atm_spalart_forceW(p) = forceFieldW(i) % force(3,c)
        else
            atm_spalart_forceW(p) = 0._rprec
        endif
    enddo
enddo

!$acc update device(atm_spalart_forceUV, atm_spalart_forceW)

if (atm_cUV_tot > 0) then
    !$acc parallel loop gang vector default(present) private(i,ii,jj,kk)
    do c = 1, atm_cUV_tot
        i = atm_tidUV(c)
        if (atm_sampling_mode(i) == ATM_SAMPLING_SPALART) then
            ii = atm_ijkUV_all(1,c)
            jj = atm_ijkUV_all(2,c)
            kk = atm_ijkUV_all(3,c)
            !$acc atomic update
            fxa(ii,jj,kk) = fxa(ii,jj,kk) + atm_spalart_forceUV(1,c)
            !$acc atomic update
            fya(ii,jj,kk) = fya(ii,jj,kk) + atm_spalart_forceUV(2,c)
        endif
    enddo
endif

if (atm_cW_tot > 0) then
    !$acc parallel loop gang vector default(present) private(i,ii,jj,kk)
    do c = 1, atm_cW_tot
        i = atm_tidW(c)
        if (atm_sampling_mode(i) == ATM_SAMPLING_SPALART) then
            ii = atm_ijkW_all(1,c)
            jj = atm_ijkW_all(2,c)
            kk = atm_ijkW_all(3,c)
            !$acc atomic update
            fza(ii,jj,kk) = fza(ii,jj,kk) + atm_spalart_forceW(c)
        endif
    enddo
endif

end subroutine atm_apply_spalart_force_gpu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_clc_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Static tables for the batched Cl/tip correction: per-point decode
! (turbine id, q index, (m,n)-block base, block length) and the static
! chord/bladeRadius flattened in the same m-outer/q-fastest order as
! atm_bp_all.
implicit none
integer :: i, j, m, n, q, p, nn2, qq2

if (atm_clc_ready) return
call atm_batch_atpoint_init()

allocate(atm_pt_turb(max(atm_nbp_tot,1)), atm_pt_q(max(atm_nbp_tot,1)))
allocate(atm_pt_base(max(atm_nbp_tot,1)), atm_pt_qq(max(atm_nbp_tot,1)))
allocate(atm_chord_all(max(atm_nbp_tot,1)), atm_brad_all(max(atm_nbp_tot,1)))
allocate(atm_db_all(max(atm_nbp_tot,1)))
allocate(atm_clc_tc(3, numberOfTurbines))
allocate(atm_wv_all(3, max(atm_nbp_tot,1)))
allocate(atm_cl_all(max(atm_nbp_tot,1)), atm_cd_all(max(atm_nbp_tot,1)))
allocate(atm_vmag_all(max(atm_nbp_tot,1)))
allocate(atm_du_all(3, max(atm_nbp_tot,1)), atm_uyopt_vec_all(3, max(atm_nbp_tot,1)))
allocate(atm_uinf_all(3, max(atm_nbp_tot,1)), atm_uxles_all(3, max(atm_nbp_tot,1)))
allocate(atm_g_all(max(atm_nbp_tot,1)), atm_dg_all(max(atm_nbp_tot,1)))
allocate(atm_epsopt_all(max(atm_nbp_tot,1)))
allocate(atm_uyles_vec_all(3, max(atm_nbp_tot,1)))
allocate(atm_uyles_all(max(atm_nbp_tot,1)), atm_uyopt_all(max(atm_nbp_tot,1)))
atm_pt_turb = 0; atm_pt_q = 0; atm_pt_base = 0; atm_pt_qq = 0
atm_chord_all = 0._rprec; atm_brad_all = 0._rprec; atm_db_all = 0._rprec
atm_clc_tc = 0._rprec
atm_wv_all = 0._rprec; atm_cl_all = 0._rprec; atm_cd_all = 0._rprec
atm_vmag_all = 1._rprec; atm_du_all = 0._rprec; atm_uyopt_vec_all = 0._rprec
atm_uinf_all = 0._rprec; atm_uxles_all = 0._rprec
atm_g_all = 0._rprec; atm_dg_all = 0._rprec; atm_epsopt_all = 0._rprec
atm_uyles_vec_all = 0._rprec; atm_uyles_all = 0._rprec; atm_uyopt_all = 0._rprec

do i = 1, numberOfTurbines
    j   = turbineArray(i) % turbineTypeID
    nn2 = turbineArray(i) % numAnnulusSections
    qq2 = turbineArray(i) % numBladePoints
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, nn2
    do q = 1, qq2
        p = p + 1
        atm_pt_turb(p)  = i
        atm_pt_q(p)     = q
        atm_pt_base(p)  = atm_bp_off(i) + ((m-1)*nn2 + (n-1))*qq2
        atm_pt_qq(p)    = qq2
        atm_chord_all(p) = turbineArray(i) % chord(m,n,q)
        atm_brad_all(p)  = turbineArray(i) % bladeRadius(m,n,q)
        atm_db_all(p)    = turbineArray(i) % db(q)
    end do
    end do
    end do
end do

!$acc enter data copyin(atm_pt_turb, atm_pt_q, atm_pt_base, atm_pt_qq,         &
!$acc                   atm_chord_all, atm_brad_all, atm_db_all, atm_clc_tc,   &
!$acc                   atm_wv_all, atm_cl_all, atm_cd_all, atm_vmag_all,      &
!$acc                   atm_du_all, atm_uyopt_vec_all,                         &
!$acc                   atm_uinf_all, atm_uxles_all,                           &
!$acc                   atm_g_all, atm_dg_all, atm_epsopt_all,                 &
!$acc                   atm_uyles_vec_all, atm_uyles_all, atm_uyopt_all)

atm_clc_ready = .true.
end subroutine atm_batch_clc_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_batch_cl_correction_gpu()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Batched OpenACC port of the current exact-panel induced-velocity Cl/tip correction for
! all operating atPoint turbines with tipALMCorrection.  The LES panel face
! construction, optimal-epsilon panel integral, ux_LES contribution to du, and
! host-visible turbineArray outputs mirror atm_compute_cl_correction().
implicit none
integer :: i, j, m, n, q, p, k, pk, nn2, qq2
logical :: any_active
real(rprec) :: f
real(rprec) :: eps_s, opt_eps_chord, inv_eps_s, inv_2pi_eps_s
real(rprec) :: inv_4pi_sqrt, mag_uinf, inv_mag_tan
real(rprec) :: z_q, z_low, z_high, z_k, half_dzp, eps_opt_k
real(rprec) :: u1, u2, u3, uval, u2val, f_low, f_high, f_up, f_um
real(rprec) :: up, um, g_over_v_k, uy_les_k, uy_opt_k
real(rprec) :: acc_les1, acc_les2, acc_opt1, acc_opt2
real(rprec) :: tan1, tan2
real(rprec), parameter :: u_small = 1.0e-2_rprec

call atm_batch_clc_init()
if (atm_nbp_tot == 0) return

f = 0.1_rprec
inv_4pi_sqrt = 1._rprec / (4._rprec * sqrt(pi))

any_active = .false.
do i = 1, numberOfTurbines
    if (turbineArray(i) % operate .and.                                        &
        turbineArray(i) % sampling == 'atPoint' .and.                          &
        (turbineArray(i) % tipALMCorrection .eqv. .true.)) then
        atm_clc_tc(3,i) = 1._rprec
        any_active = .true.
    else
        atm_clc_tc(3,i) = 0._rprec
        cycle
    end if
    atm_clc_tc(1,i) = turbineArray(i) % epsilon
    atm_clc_tc(2,i) = turbineArray(i) % optimalEpsilonChord

    j   = turbineArray(i) % turbineTypeID
    nn2 = turbineArray(i) % numAnnulusSections
    qq2 = turbineArray(i) % numBladePoints
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, nn2
    do q = 1, qq2
        p = p + 1
        atm_wv_all(1,p) = turbineArray(i) % windVectors(m,n,q,1)
        atm_wv_all(2,p) = turbineArray(i) % windVectors(m,n,q,2)
        atm_wv_all(3,p) = turbineArray(i) % windVectors(m,n,q,3)
        atm_cl_all(p)   = turbineArray(i) % Cl(m,n,q)
        atm_cd_all(p)   = turbineArray(i) % Cd(m,n,q)
        atm_vmag_all(p) = turbineArray(i) % Vmag(m,n,q)
        atm_du_all(1,p) = turbineArray(i) % du(m,n,q,1)
        atm_du_all(2,p) = turbineArray(i) % du(m,n,q,2)
        atm_du_all(3,p) = turbineArray(i) % du(m,n,q,3)
        atm_uyopt_vec_all(1,p) = turbineArray(i) % uy_opt_vec(m,n,q,1)
        atm_uyopt_vec_all(2,p) = turbineArray(i) % uy_opt_vec(m,n,q,2)
        atm_uyopt_vec_all(3,p) = turbineArray(i) % uy_opt_vec(m,n,q,3)
    end do
    end do
    end do
end do
if (.not. any_active) return

!$acc update device(atm_wv_all, atm_cl_all, atm_cd_all, atm_vmag_all,          &
!$acc               atm_du_all, atm_uyopt_vec_all, atm_clc_tc)

! K1: point-local Uinf, ux_LES, G, epsilon_opt.
!$acc parallel loop gang vector default(present)                               &
!$acc     private(i, eps_s, opt_eps_chord, u1, u2, u3, mag_uinf)
do p = 1, atm_nbp_tot
    i = atm_pt_turb(p)
    if (atm_clc_tc(3,i) > 0.5_rprec) then
        eps_s = atm_clc_tc(1,i)
        opt_eps_chord = atm_clc_tc(2,i)
        u1 = atm_wv_all(1,p) - atm_uyopt_vec_all(1,p)
        u2 = atm_wv_all(2,p) - atm_uyopt_vec_all(2,p)
        u3 = atm_wv_all(3,p) - atm_uyopt_vec_all(3,p)
        atm_uinf_all(1,p) = u1
        atm_uinf_all(2,p) = u2
        atm_uinf_all(3,p) = u3
        mag_uinf = sqrt(u1*u1 + u2*u2 + u3*u3)
        atm_uxles_all(1,p) = atm_cd_all(p) * atm_chord_all(p) / eps_s          &
            * inv_4pi_sqrt * atm_vmag_all(p) * u1 / mag_uinf
        atm_uxles_all(2,p) = atm_cd_all(p) * atm_chord_all(p) / eps_s          &
            * inv_4pi_sqrt * atm_vmag_all(p) * u2 / mag_uinf
        atm_uxles_all(3,p) = atm_cd_all(p) * atm_chord_all(p) / eps_s          &
            * inv_4pi_sqrt * atm_vmag_all(p) * u3 / mag_uinf
        atm_g_all(p) = 0.5_rprec * atm_cl_all(p) * atm_chord_all(p)            &
            * atm_vmag_all(p) * atm_vmag_all(p)
        atm_epsopt_all(p) = atm_chord_all(p) * opt_eps_chord
    end if
end do

! K2: dG is retained for backward-compatible diagnostics/output.
!$acc parallel loop gang vector default(present) private(i, q)
do p = 1, atm_nbp_tot
    i = atm_pt_turb(p)
    if (atm_clc_tc(3,i) > 0.5_rprec) then
        q = atm_pt_q(p)
        if (q == 1) then
            atm_dg_all(p) = atm_g_all(p)
        else if (q == atm_pt_qq(p)) then
            atm_dg_all(p) = -atm_g_all(p)
        else
            atm_dg_all(p) = (atm_g_all(p+1) - atm_g_all(p-1)) * 0.5_rprec
        end if
    end if
end do

! K3: exact panel induced velocity, one thread per blade point, k loop in the
! same ascending order as the host routine.
!$acc parallel loop gang vector default(present)                               &
!$acc     private(i, q, k, pk, eps_s, inv_eps_s, inv_2pi_eps_s, z_q, z_low,    &
!$acc             z_high, z_k, half_dzp, eps_opt_k, uval, u2val, f_low,        &
!$acc             f_high, f_up, f_um, up, um, g_over_v_k, uy_les_k,            &
!$acc             uy_opt_k, inv_mag_tan, tan1, tan2, acc_les1, acc_les2,       &
!$acc             acc_opt1, acc_opt2)
do p = 1, atm_nbp_tot
    i = atm_pt_turb(p)
    if (atm_clc_tc(3,i) > 0.5_rprec) then
        q = atm_pt_q(p)
        eps_s = atm_clc_tc(1,i)
        inv_eps_s = 1._rprec / eps_s
        inv_2pi_eps_s = inv_eps_s / (2._rprec * pi)
        z_q = atm_brad_all(p)
        acc_les1 = 0._rprec
        acc_les2 = 0._rprec
        acc_opt1 = 0._rprec
        acc_opt2 = 0._rprec
        !$acc loop seq
        do k = 1, atm_pt_qq(p)
            pk = atm_pt_base(p) + k
            g_over_v_k = atm_g_all(pk) / atm_vmag_all(pk)

            if (k == 1) then
                z_low = atm_brad_all(pk) - 0.5_rprec * atm_db_all(pk)
            else
                z_low = atm_brad_all(pk-1) + 0.5_rprec * atm_db_all(pk-1)
            end if
            z_high = atm_brad_all(pk) + 0.5_rprec * atm_db_all(pk)

            uval = (z_q - z_low) * inv_eps_s
            u2val = uval * uval
            if (abs(uval) < u_small) then
                f_low = 0.5_rprec * uval *                                    &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_low = (1._rprec - exp(-u2val)) / (2._rprec * uval)
            endif
            uval = (z_q - z_high) * inv_eps_s
            u2val = uval * uval
            if (abs(uval) < u_small) then
                f_high = 0.5_rprec * uval *                                   &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_high = (1._rprec - exp(-u2val)) / (2._rprec * uval)
            endif
            uy_les_k = -g_over_v_k * inv_2pi_eps_s * (f_low - f_high)

            z_k = atm_brad_all(pk)
            half_dzp = 0.5_rprec * atm_db_all(pk)
            eps_opt_k = atm_epsopt_all(pk)
            up = (z_q - z_k + half_dzp) / eps_opt_k
            um = (z_q - z_k - half_dzp) / eps_opt_k
            u2val = up * up
            if (abs(up) < u_small) then
                f_up = 0.5_rprec * up *                                      &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_up = (1._rprec - exp(-u2val)) / (2._rprec * up)
            endif
            u2val = um * um
            if (abs(um) < u_small) then
                f_um = 0.5_rprec * um *                                      &
                    (1._rprec - 0.5_rprec * u2val *                           &
                    (1._rprec - (u2val / 3._rprec) *                          &
                    (1._rprec - 0.25_rprec * u2val)))
            else
                f_um = (1._rprec - exp(-u2val)) / (2._rprec * um)
            endif
            uy_opt_k = -g_over_v_k / (2._rprec * pi * eps_opt_k) *            &
                (f_up - f_um)

            inv_mag_tan = 1._rprec / sqrt(atm_uinf_all(1,pk)**2 +             &
                                          atm_uinf_all(2,pk)**2)
            tan1 =  atm_uinf_all(2,pk) * inv_mag_tan
            tan2 = -atm_uinf_all(1,pk) * inv_mag_tan
            acc_les1 = acc_les1 + uy_les_k * tan1
            acc_les2 = acc_les2 + uy_les_k * tan2
            acc_opt1 = acc_opt1 + uy_opt_k * tan1
            acc_opt2 = acc_opt2 + uy_opt_k * tan2
        end do
        atm_uyles_vec_all(1,p) = acc_les1
        atm_uyles_vec_all(2,p) = acc_les2
        atm_uyles_vec_all(3,p) = 0._rprec
        atm_uyopt_vec_all(1,p) = acc_opt1
        atm_uyopt_vec_all(2,p) = acc_opt2
        atm_uyopt_vec_all(3,p) = 0._rprec
        atm_uyles_all(p) = sqrt(acc_les1*acc_les1 + acc_les2*acc_les2)
        atm_uyopt_all(p) = sqrt(acc_opt1*acc_opt1 + acc_opt2*acc_opt2)
        atm_du_all(1,p) = atm_du_all(1,p) * (1._rprec - f) + f *              &
            (atm_uyopt_vec_all(1,p) - atm_uyles_vec_all(1,p) +                &
             atm_uxles_all(1,p))
        atm_du_all(2,p) = atm_du_all(2,p) * (1._rprec - f) + f *              &
            (atm_uyopt_vec_all(2,p) - atm_uyles_vec_all(2,p) +                &
             atm_uxles_all(2,p))
        atm_du_all(3,p) = atm_du_all(3,p) * (1._rprec - f) + f *              &
            (atm_uyopt_vec_all(3,p) - atm_uyles_vec_all(3,p) +                &
             atm_uxles_all(3,p))
    end if
end do

!$acc update self(atm_uinf_all, atm_uxles_all, atm_g_all, atm_dg_all,          &
!$acc             atm_epsopt_all, atm_uyles_vec_all, atm_uyopt_vec_all,        &
!$acc             atm_uyles_all, atm_uyopt_all, atm_du_all)

! Unpack everything the host routine writes, preserving host-visible state for
! downstream force, structural, restart, and output consumers.
do i = 1, numberOfTurbines
    if (atm_clc_tc(3,i) < 0.5_rprec) cycle
    j   = turbineArray(i) % turbineTypeID
    nn2 = turbineArray(i) % numAnnulusSections
    qq2 = turbineArray(i) % numBladePoints
    p = atm_bp_off(i)
    do m = 1, turbineModel(j) % numBl
    do n = 1, nn2
    do q = 1, qq2
        p = p + 1
        turbineArray(i) % Uinf_vec(m,n,q,1:3)   = atm_uinf_all(1:3,p)
        turbineArray(i) % ux_LES_vec(m,n,q,1:3) = atm_uxles_all(1:3,p)
        turbineArray(i) % G(m,n,q)              = atm_g_all(p)
        turbineArray(i) % dG(m,n,q)             = atm_dg_all(p)
        turbineArray(i) % epsilon_opt(m,n,q)    = atm_epsopt_all(p)
        turbineArray(i) % uy_LES_vec(m,n,q,1:3) = atm_uyles_vec_all(1:3,p)
        turbineArray(i) % uy_opt_vec(m,n,q,1:3) = atm_uyopt_vec_all(1:3,p)
        turbineArray(i) % uy_LES(m,n,q)         = atm_uyles_all(p)
        turbineArray(i) % uy_opt(m,n,q)         = atm_uyopt_all(p)
        turbineArray(i) % du(m,n,q,1:3)         = atm_du_all(1:3,p)
    end do
    end do
    end do
end do

end subroutine atm_batch_cl_correction_gpu
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_lesgo_apply_force()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This will apply the blade force onto the CFD grid by using the convolution
! function in the ATM library
implicit none

integer :: c,m
integer :: i,j,k

#ifdef PPLES_GPU
! atPoint forces were deposited by atm_batch_convolute_force_gpu. Optional
! Spalart forces are selectively packed and deposited without copying the
! complete LES force fields through the host.
call atm_apply_spalart_force_gpu()
return
#endif

do m=1, numberOfTurbines

    if (turbineArray(m) % operate) then
        ! Impose force field onto the flow field variables
        ! The forces are non-dimensionalized here as well
        do c=1,forceFieldUV(m) % c
            i=forceFieldUV(m) % ijk(1,c)
            j=forceFieldUV(m) % ijk(2,c)
            k=forceFieldUV(m) % ijk(3,c)

            fxa(i,j,k) = fxa(i,j,k) + forceFieldUV(m) % force(1,c)
            fya(i,j,k) = fya(i,j,k) + forceFieldUV(m) % force(2,c)

        enddo

        do c=1,forceFieldW(m) % c
            i=forceFieldW(m) % ijk(1,c)
            j=forceFieldW(m) % ijk(2,c)
            k=forceFieldW(m) % ijk(3,c)

            fza(i,j,k) = fza(i,j,k) + forceFieldW(m) % force(3,c)

        enddo
    endif
enddo

end subroutine atm_lesgo_apply_force


end module atm_lesgo_interface
