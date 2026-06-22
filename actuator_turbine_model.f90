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
module actuator_turbine_model
!*******************************************************************************
! This module has the subroutines to provide all calculations for use in the
! actuator turbine model (ATM)
!
! Navigation map:
!   - environment and diagnostics: atm_model_env_token, structure timing
!   - lifecycle and I/O: atm_initialize, restart, output initialization
!   - geometry and kinematics: create_points, update, yaw, rotor speed, rotate
!   - aerodynamic corrections: atm_compute_cl_correction and GPU companion
!   - blade/nacelle physics: calculate_variables, computeBladeForce,
!     computeNacelleForce, integrate_u, yawNacelle
!   - output and power: atm_output, atm_compute_power, process_output
!   - structural solver: atm_solve_structure and linear-system helpers
!
! This file owns turbine-side physics and structural state.  LESGO grid
! sampling and force deposition belong in atm_lesgo_interface.f90.
!
! Ownership map:
!   - turbineArray and turbineModel state are authoritative on the host here.
!   - structural-solver state and feedback gates are owned by this module.
!   - GPU companion routines mirror turbine-side computations but copy required
!     outputs back so host turbine state remains authoritative.
!   - turbine power, thrust, yaw, nacelle, and structural output files are
!     written here; LESGO field/checkpoint output remains in io.f90.

! Imported modules
use atm_base ! Include basic types and precision of real numbers

use atm_input_util ! Utilities to read input files
use param, only : coord

#ifdef ENABLE_CUDA
use cudafor
use cublas
use cusolverDn
#endif

implicit none

! Declare everything private except for subroutine which will be used
private
public :: atm_initialize, numberOfTurbines,                                    &
          atm_computeBladeForce, atm_update,                                   &
          vector_add, vector_divide, vector_mag, distance,                     &
          atm_output, atm_process_output,                                      &
          atm_initialize_output, atm_computeNacelleForce, atm_write_restart,   &
          atm_compute_cl_correction, atm_structure_enabled,                    &
          atm_structure_timing_report

! The very crucial parameter pi
real(rprec), parameter :: pi=acos(-1._rprec)

! These are used to do unit conversions
real(rprec) :: degRad = pi/180._rprec ! Degrees to radians conversion
real(rprec) :: rpmRadSec =  pi/30._rprec ! Set the revolutions/min to radians/s

logical :: pastFirstTimeStep ! Establishes if we are at the first time step

integer, save :: atm_structure_timing_calls = 0
real(rprec), save :: atm_structure_time_total = 0._rprec
real(rprec), save :: atm_structure_time_assembly = 0._rprec
real(rprec), save :: atm_structure_time_solve = 0._rprec
real(rprec), save :: atm_structure_time_update = 0._rprec
real(rprec), save :: atm_structure_time_other = 0._rprec

! Subroutines for the actuator turbine model
! All suboroutines names start with (atm_)
contains

! Environment-controlled ATM policy gates.  Keep this block synchronized with
! docs/environment_switches.md when adding, removing, or changing a switch.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_model_env_token(name, token, has_value)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(len=*), intent(in) :: name
character(len=*), intent(out) :: token
logical, intent(out) :: has_value
integer :: env_len, env_stat

token = ''
has_value = .false.
call get_environment_variable(name, token, env_len, env_stat)
if (env_stat == 0 .and. env_len > 0) then
    token = trim(adjustl(token(1:env_len)))
    has_value = .true.
endif

end subroutine atm_model_env_token

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_model_env_enabled(name)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(len=*), intent(in) :: name
character(len=32) :: env_value
logical :: has_value

atm_model_env_enabled = .false.
call atm_model_env_token(name, env_value, has_value)
if (has_value) then
    select case (env_value)
    case ('1','T','t','TRUE','true','True','Y','y','YES','yes','Yes')
        atm_model_env_enabled = .true.
    end select
endif

end function atm_model_env_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_STRUCTURE')
    initialized = .true.
endif

atm_structure_enabled = active

end function atm_structure_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_diag_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_STRUCTURE_DIAG')
    initialized = .true.
endif

atm_structure_diag_enabled = active

end function atm_structure_diag_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_timing_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_STRUCTURE_TIMING')
    initialized = .true.
endif

atm_structure_timing_enabled = active

end function atm_structure_timing_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_power_stdout_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Turbine power is written to turbineOutput/*/power.  Repeated stdout power
! lines are useful only for debugging and become serial log I/O in large farms.
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_POWER_STDOUT')
    initialized = .true.
endif

atm_power_stdout_enabled = active

end function atm_power_stdout_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_structure_timing_report()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

if (.not. atm_structure_timing_enabled()) return
if (atm_structure_timing_calls <= 0) return

write(*,'(a,i0,a,i0)') 'ATM structure timing rank=', coord, ' calls=',        &
    atm_structure_timing_calls
write(*,'(1a,E15.7)') '  structure total: ', atm_structure_time_total
write(*,'(1a,E15.7)') '  matrix/force assembly: ',                           &
    atm_structure_time_assembly
write(*,'(1a,E15.7)') '  linear solves: ', atm_structure_time_solve
write(*,'(1a,E15.7)') '  state update: ', atm_structure_time_update
write(*,'(1a,E15.7)') '  allocation/other: ', atm_structure_time_other

end subroutine atm_structure_timing_report

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_structure_diag_snapshot(i, stage, dt, time)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compact rank-local diagnostic for structure-on GPU-port work.  Disabled by
! default; when LESGO_ATM_STRUCTURE_DIAG=1 it records norms of the structural
! inputs and state immediately before/after atm_solve_structure().
implicit none

integer, intent(in) :: i
character(len=*), intent(in) :: stage
real(rprec), intent(in) :: dt, time

integer :: unit, ios
logical :: exists
character(len=512) :: diag_path
character(len=32) :: rank_tag
real(rprec) :: blade_force_l1, blade_force_linf
real(rprec) :: pitch_moment_l1, pitch_moment_linf
real(rprec) :: blade_aligned_l1, blade_aligned_linf
real(rprec) :: blade_points_l1, blade_points_linf
real(rprec) :: flap_disp_l1, flap_disp_linf
real(rprec) :: edge_disp_l1, edge_disp_linf
real(rprec) :: elastic_twist_l1, elastic_twist_linf
real(rprec) :: alpha_l1, alpha_linf, cl_l1, cl_linf, cd_l1, cd_linf, cm_l1, cm_linf

if (.not. atm_structure_diag_enabled()) return

write(rank_tag, '(i0)') coord
diag_path = './turbineOutput/'//trim(turbineArray(i) % turbineName)//          &
            '/structure_diag_rank'//trim(rank_tag)
inquire(file=trim(diag_path), exist=exists)

blade_force_l1 = sum(abs(turbineArray(i) % bladeForces))
blade_force_linf = maxval(abs(turbineArray(i) % bladeForces))
pitch_moment_l1 = sum(abs(turbineArray(i) % pitchingMoment))
pitch_moment_linf = maxval(abs(turbineArray(i) % pitchingMoment))
blade_aligned_l1 = sum(abs(turbineArray(i) % bladeAlignedVectors))
blade_aligned_linf = maxval(abs(turbineArray(i) % bladeAlignedVectors))
blade_points_l1 = sum(abs(turbineArray(i) % bladePoints))
blade_points_linf = maxval(abs(turbineArray(i) % bladePoints))
flap_disp_l1 = sum(abs(turbineArray(i) % flap_disp))
flap_disp_linf = maxval(abs(turbineArray(i) % flap_disp))
edge_disp_l1 = sum(abs(turbineArray(i) % edge_disp))
edge_disp_linf = maxval(abs(turbineArray(i) % edge_disp))
elastic_twist_l1 = sum(abs(turbineArray(i) % elastic_twist))
elastic_twist_linf = maxval(abs(turbineArray(i) % elastic_twist))
alpha_l1 = sum(abs(turbineArray(i) % alpha))
alpha_linf = maxval(abs(turbineArray(i) % alpha))
cl_l1 = sum(abs(turbineArray(i) % Cl))
cl_linf = maxval(abs(turbineArray(i) % Cl))
cd_l1 = sum(abs(turbineArray(i) % Cd))
cd_linf = maxval(abs(turbineArray(i) % Cd))
cm_l1 = sum(abs(turbineArray(i) % Cm))
cm_linf = maxval(abs(turbineArray(i) % Cm))

open(newunit=unit, file=trim(diag_path), status='unknown', position='append',  &
     action='write', iostat=ios)
if (ios /= 0) return
if (.not. exists) then
    write(unit,'(a)') '# stage time dt bf_l1 bf_linf pm_l1 pm_linf '//        &
        'bav_l1 bav_linf bp_l1 bp_linf flap_l1 flap_linf edge_l1 edge_linf '//&
        'twist_l1 twist_linf alpha_l1 alpha_linf cl_l1 cl_linf cd_l1 '//     &
        'cd_linf cm_l1 cm_linf'
endif
write(unit,'(a,1x,24(es25.16e3,1x))') trim(stage), time, dt,                &
    blade_force_l1, blade_force_linf, pitch_moment_l1, pitch_moment_linf,     &
    blade_aligned_l1, blade_aligned_linf, blade_points_l1, blade_points_linf, &
    flap_disp_l1, flap_disp_linf, edge_disp_l1, edge_disp_linf,               &
    elastic_twist_l1, elastic_twist_linf, alpha_l1, alpha_linf, cl_l1,        &
    cl_linf, cd_l1, cd_linf, cm_l1, cm_linf
close(unit)

end subroutine atm_structure_diag_snapshot

#ifdef ENABLE_CUDA
! GPU structural-solver validation and path-selection gates.  These switches are
! diagnostic/path controls; they should not be used to hide production defaults.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_gpu_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_STRUCTURE_GPU')
    initialized = .true.
endif

atm_structure_gpu_enabled = active

end function atm_structure_gpu_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_gpu_validate_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_STRUCTURE_GPU_VALIDATE')
    initialized = .true.
endif

atm_structure_gpu_validate_enabled = active

end function atm_structure_gpu_validate_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_gpu_direct_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

logical, save :: initialized = .false.
logical, save :: active = .false.

if (.not. initialized) then
    active = atm_model_env_enabled('LESGO_ATM_STRUCTURE_GPU_DIRECT')
    initialized = .true.
endif

atm_structure_gpu_direct_enabled = active

end function atm_structure_gpu_direct_enabled
#endif

#ifdef ENABLE_CUDA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_model_cuda_enabled()
!*******************************************************************************
implicit none

atm_model_cuda_enabled = .true.

end function atm_model_cuda_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_model_extra_sync_enabled()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

atm_model_extra_sync_enabled = .false.

end function atm_model_extra_sync_enabled

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_model_cuda_check(where)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: where
integer :: istat

if (atm_model_extra_sync_enabled()) then
    call atm_model_cuda_sync(where)
    return
end if

istat = cudaGetLastError()
if (istat /= cudaSuccess) then
    print *, 'ATM model CUDA kernel failure at ', trim(where), ': ', istat
    stop 1
end if

end subroutine atm_model_cuda_check

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_model_cuda_sync(where)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

character(*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= cudaSuccess) then
    print *, 'ATM model CUDA sync failure at ', trim(where), ': ', istat
    stop 1
end if
istat = cudaGetLastError()
if (istat /= cudaSuccess) then
    print *, 'ATM model CUDA kernel failure at ', trim(where), ': ', istat
    stop 1
end if

end subroutine atm_model_cuda_sync
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize()
! This subroutine initializes the ATM. It calls the subroutines in
! atm_input_util to read the input data and creates the initial geometry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
integer :: i
logical :: file_exists

pastFirstTimeStep=.false. ! The first time step not reached yet

write(*,*) 'Reading Actuator Turbine Model Input...'
call read_input_conf()  ! Read input data
write(*,*) 'Done Reading Actuator Turbine Model Input'
do i = 1,numberOfTurbines
    inquire(file = "./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/actuatorPoints", exist=file_exists)

    ! Creates the ATM points defining the geometry
    call atm_create_points(i)
    ! This will create the first yaw alignment
    turbineArray(i) % deltaNacYaw = turbineArray(i) % nacYaw
    call atm_yawNacelle(i)

    if (file_exists .eqv. .true.) then
        write(*,*) 'Reading bladePoints from Previous Simulation'
        call atm_read_actuator_points(i)
    endif

    call atm_calculate_variables(i) ! Calculates variables depending on input

    inquire(file = "./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/restart", exist=file_exists)

    if (file_exists .eqv. .true.) then
        write(*,*) 'Reading Turbine Properties from Previous Simulation'
        call atm_read_restart(i)
    endif

end do

pastFirstTimeStep=.true. ! Past the first time step

end subroutine atm_initialize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_read_actuator_points(i)
! This subroutine reads the location of the actuator points
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original location
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i ! Indicates the turbine number

integer :: j, m, n, q
integer :: rigidPointsFile=788
logical :: rigid_exists
character(len=512) :: rigid_points_path

j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=1, file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                  "/actuatorPoints", action='read')

do m=1, turbineModel(j) % numBl
    do n=1, turbineArray(i) %  numAnnulusSections
        do q=1, turbineArray(i) % numBladePoints
            read(1,*) turbineArray(i) % bladePoints(m,n,q,:)
        enddo
    enddo
enddo

close(1)

if (atm_structure_enabled()) then
    rigid_points_path = "./turbineOutput/"//                                  &
        trim(turbineArray(i) % turbineName)//"/actuatorPoints_rigid"
    inquire(file=trim(rigid_points_path), exist=rigid_exists)
    if (rigid_exists) then
        open(unit=rigidPointsFile, file=trim(rigid_points_path), action='read')
        do m=1, turbineModel(j) % numBl
            do n=1, turbineArray(i) %  numAnnulusSections
                do q=1, turbineArray(i) % numBladePoints
                    read(rigidPointsFile,*)                                   &
                        turbineArray(i) % bladePoints_rigid(m,n,q,:)
                enddo
            enddo
        enddo
        close(rigidPointsFile)
    else
        turbineArray(i) % bladePoints_rigid = turbineArray(i) % bladePoints
        write(*,*) 'Flexible restart has no actuatorPoints_rigid; ',          &
            'using actuatorPoints as the rigid reference.'
    endif
endif

end subroutine atm_read_actuator_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_read_restart(i)
! This subroutine reads the rotor speed
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original omega
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i  ! Indicates the turbine number
integer :: structureFile=789
logical :: structure_exists
character(len=512) :: structure_path

! Open the file at the last line (append)
open( unit=1, file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/restart", action='read') !, position='append')

! Bring the pointer to the last line
!~ backspace 1

! Read past the first line
read(1,*)

! Read the restart variables
read(1,*) turbineArray(i) % rotSpeed
read(1,*) turbineArray(i) % torqueGen
read(1,*) turbineArray(i) % torqueRotor
read(1,*) turbineArray(i) % u_infinity
read(1,*) turbineArray(i) % induction_a
read(1,*) turbineArray(i) % PitchControlAngle
read(1,*) turbineArray(i) % IntSpeedError
read(1,*) turbineArray(i) % nacYaw
read(1,*) turbineArray(i) % rotorApex
read(1,*) turbineArray(i) % uvShaft
close(1)

write(*,*) ' RotSpeed Value from previous simulation is ',                     &
                turbineArray(i) % rotSpeed
write(*,*) ' torqueGen Value from previous simulation is ',                    &
                turbineArray(i) % torqueGen
write(*,*) ' torqueRotor Value from previous simulation is ',                  &
                turbineArray(i) % torqueRotor
write(*,*) ' PitchControlAngle Value from previous simulation is ',            &
                turbineArray(i) % PitchControlAngle
write(*,*) ' IntSpeedError Value from previous simulation is ',                &
                turbineArray(i) % IntSpeedError
write(*,*) ' Yaw Value from previous simulation is ',                          &
                turbineArray(i) % nacYaw
write(*,*) ' Rotor Apex Value from previous simulation is ',                   &
                turbineArray(i) % rotorApex
write(*,*) ' uvShaft Value from previous simulation is ',                      &
                turbineArray(i) % uvShaft

if (atm_structure_enabled()) then
    structure_path = "./turbineOutput/"//trim(turbineArray(i) % turbineName)// &
        "/structure_restart"
    inquire(file=trim(structure_path), exist=structure_exists)
    if (structure_exists) then
        open(unit=structureFile, file=trim(structure_path), action='read')
        read(structureFile,*)
        read(structureFile,*) turbineArray(i) % flap_disp
        read(structureFile,*) turbineArray(i) % theta_disp
        read(structureFile,*) turbineArray(i) % flap_vel
        read(structureFile,*) turbineArray(i) % theta_vel
        read(structureFile,*) turbineArray(i) % flap_acc
        read(structureFile,*) turbineArray(i) % theta_acc
        read(structureFile,*) turbineArray(i) % edge_disp
        read(structureFile,*) turbineArray(i) % edge_theta_disp
        read(structureFile,*) turbineArray(i) % edge_vel
        read(structureFile,*) turbineArray(i) % edge_theta_vel
        read(structureFile,*) turbineArray(i) % edge_acc
        read(structureFile,*) turbineArray(i) % edge_theta_acc
        read(structureFile,*) turbineArray(i) % elastic_twist
        close(structureFile)
        write(*,*) ' Structural state restored from previous simulation.'
    else
        write(*,*) 'Flexible restart has no structure_restart; ',              &
            'starting structural DOFs from their initialized values.'
    endif
endif

end subroutine atm_read_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_write_restart(i)
! This subroutine reads the rotor speed
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original omega
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: pointsFile=787 ! File to write the actuator points
integer :: restartFile=21 ! File to write restart data
integer :: rigidPointsFile=788 ! File to write undeformed blade backbone points
integer :: structureFile=789 ! File to write structural restart data
integer j, m,n,q ! counters

! Open the file
open( unit=restartFile, file="./turbineOutput/"//                              &
            trim(turbineArray(i) % turbineName)//"/restart", status="replace")

write(restartFile,*) 'RotSpeed ', 'torqueGen ', 'torqueRotor ', 'u_infinity ', &
                     'induction_a ', 'PitchControlAngle ', 'IntSpeedError ',   &
                     'nacYaw ', 'rotorApex ', 'uvShaft'
! Store the rotSpeed value
write(restartFile,*) turbineArray(i) % rotSpeed
write(restartFile,*) turbineArray(i) % torqueGen
write(restartFile,*) turbineArray(i) % torqueRotor
write(restartFile,*) turbineArray(i) % u_infinity
write(restartFile,*) turbineArray(i) % induction_a
write(restartFile,*) turbineArray(i) % PitchControlAngle
write(restartFile,*) turbineArray(i) % IntSpeedError
write(restartFile,*) turbineArray(i) % nacYaw
write(restartFile,*) turbineArray(i) % rotorApex
write(restartFile,*) turbineArray(i) % uvShaft
close(restartFile)

! Write the actuator points at every time-step regardless
j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=pointsFile, status="replace", file="./turbineOutput/"//              &
                      trim(turbineArray(i) % turbineName)//"/actuatorPoints")

do m=1, turbineModel(j) % numBl
    do n=1, turbineArray(i) %  numAnnulusSections
        do q=1, turbineArray(i) % numBladePoints
            ! A new file will be created each time-step with the proper
            ! location of the blades
            write(pointsFile,*) turbineArray(i) % bladePoints(m,n,q,:)
        enddo
    enddo
enddo

close(pointsFile)

if (atm_structure_enabled()) then
    open(unit=rigidPointsFile, status="replace", file="./turbineOutput/"//      &
                          trim(turbineArray(i) % turbineName)//                &
                          "/actuatorPoints_rigid")

    do m=1, turbineModel(j) % numBl
        do n=1, turbineArray(i) %  numAnnulusSections
            do q=1, turbineArray(i) % numBladePoints
                write(rigidPointsFile,*)                                       &
                    turbineArray(i) % bladePoints_rigid(m,n,q,:)
            enddo
        enddo
    enddo

    close(rigidPointsFile)

    open(unit=structureFile, status="replace", file="./turbineOutput/"//       &
                          trim(turbineArray(i) % turbineName)//                &
                          "/structure_restart")
    write(structureFile,*) 'flap_disp ', 'theta_disp ', 'flap_vel ',           &
        'theta_vel ', 'flap_acc ', 'theta_acc ', 'edge_disp ',                 &
        'edge_theta_disp ', 'edge_vel ', 'edge_theta_vel ', 'edge_acc ',       &
        'edge_theta_acc ', 'elastic_twist'
    write(structureFile,*) turbineArray(i) % flap_disp
    write(structureFile,*) turbineArray(i) % theta_disp
    write(structureFile,*) turbineArray(i) % flap_vel
    write(structureFile,*) turbineArray(i) % theta_vel
    write(structureFile,*) turbineArray(i) % flap_acc
    write(structureFile,*) turbineArray(i) % theta_acc
    write(structureFile,*) turbineArray(i) % edge_disp
    write(structureFile,*) turbineArray(i) % edge_theta_disp
    write(structureFile,*) turbineArray(i) % edge_vel
    write(structureFile,*) turbineArray(i) % edge_theta_vel
    write(structureFile,*) turbineArray(i) % edge_acc
    write(structureFile,*) turbineArray(i) % edge_theta_acc
    write(structureFile,*) turbineArray(i) % elastic_twist
    close(structureFile)
endif

end subroutine atm_write_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_initialize_output()
! This subroutine initializes the output files for the ATM
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none
logical :: file_exists
integer :: i

! Write to the screen output start
call atm_print_initialize()

do i = 1,numberOfTurbines

    inquire(file="./turbineOutput/"//                                   &
                     trim(turbineArray(i) % turbineName),EXIST=file_exists)

    if (file_exists .eqv. .false.) then

        ! Create turbineOutput directory
        call system("mkdir -vp turbineOutput/"//                               &
                     trim(turbineArray(i) % turbineName))

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/power")
        write(1,*) 'time PowerRotor powerGen '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/thrust")
        write(1,*) 'time thrust '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/RotSpeed")
        write(1,*) 'time RotSpeed'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Yaw")
        write(1,*) 'time deltaNacYaw NacYaw'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/lift")
        write(1,*) 'turbineNumber bladeNumber '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/drag")
        write(1,*) 'turbineNumber bladeNumber '
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Cl")
        write(1,*) 'turbineNumber bladeNumber Cl'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Cd")
        write(1,*) 'turbineNumber bladeNumber Cd'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/alpha")
        write(1,*) 'turbineNumber bladeNumber alpha'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Vrel")
        write(1,*) 'turbineNumber bladeNumber Vrel'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Vaxial")
        write(1,*) 'turbineNumber bladeNumber Vaxial'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/Vtangential")
        write(1,*) 'turbineNumber bladeNumber Vtangential'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/tangentialForce")
        write(1,*) 'turbineNumber bladeNumber tangentialForce'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/axialForce")
        write(1,*) 'turbineNumber bladeNumber axialForce'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/nacelle")
        write(1,*) 'time Velocity-no-correction Velocity-w-correction'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/uy_LES")
        write(1,*) 'turbineNumber bladeNumber uy_LES'
        close(1)

        open(unit=1, file="./turbineOutput/"//                                 &
                     trim(turbineArray(i) % turbineName)//"/uy_opt")
        write(1,*) 'turbineNumber bladeNumber uy_opt'
        close(1)

        if (atm_structure_enabled()) then
            open(unit=1, file="./turbineOutput/"//                             &
                         trim(turbineArray(i) % turbineName)//"/elastic_twist")
            write(1,*) 'turbineNumber bladeNumber elastic_twist'
            close(1)

            open(unit=1, file="./turbineOutput/"//                             &
                         trim(turbineArray(i) % turbineName)//"/flap_disp")
            write(1,*) 'turbineNumber bladeNumber flap_disp'
            close(1)

            open(unit=1, file="./turbineOutput/"//                             &
                         trim(turbineArray(i) % turbineName)//"/edge_disp")
            write(1,*) 'turbineNumber bladeNumber edge_disp'
            close(1)

            open(unit=1, file="./turbineOutput/"//                             &
                         trim(turbineArray(i) % turbineName)//"/Cm")
            write(1,*) 'turbineNumber bladeNumber Cm'
            close(1)
        endif

    endif
enddo

end subroutine atm_initialize_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_create_points(i)
! This subroutine generate the set of blade points for each turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
integer :: m ! Indicates the blade point number
integer :: n ! Indicates number of actuator section
integer :: k
real(rprec), dimension (3) :: root ! Location of rotor apex
real(rprec) :: beta ! Difference between coning angle and shaft tilt
real(rprec) :: dist ! Distance from each actuator point
integer,     pointer :: numBladePoints
integer,     pointer :: numBl
integer,     pointer :: numAnnulusSections
real(rprec), pointer :: NacYaw
real(rprec),  pointer :: db(:)
real(rprec),  pointer :: bladePoints(:,:,:,:)
real(rprec),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: azimuth
real(rprec), pointer :: rotSpeed
real(rprec), pointer :: ShftTilt
real(rprec), pointer :: PreCone
real(rprec), pointer :: towerShaftIntersect(:)
real(rprec), pointer :: baseLocation(:)
real(rprec), pointer :: TowerHt
real(rprec), pointer :: Twr2Shft
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: uvShaftDir
real(rprec), pointer :: uvShaft(:)
real(rprec), pointer :: uvTower(:)
real(rprec), pointer :: TipRad
real(rprec), pointer :: HubRad
real(rprec), pointer :: annulusSectionAngle
real(rprec), pointer :: solidity(:,:,:)

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Variables to be used locally. They are stored in local variables within the
! subroutine for easier code following. The values are then passed to the
! proper type
numBladePoints => turbineArray(i) % numBladePoints
numBl=>turbineModel(j) % numBl
numAnnulusSections=>turbineArray(i) % numAnnulusSections

! Allocate variables depending on specific turbine properties and general
! turbine model properties
allocate(turbineArray(i) % db(numBladePoints))

allocate(turbineArray(i) % bladePoints(numBl, numAnnulusSections, &
         numBladePoints,3))

allocate(turbineArray(i) % bladeRadius(numBl,numAnnulusSections,numBladePoints))

allocate(turbineArray(i) % solidity(numBl,numAnnulusSections,numBladePoints))

! Assign Pointers turbineArray denpendent (i)
db=>turbineArray(i) % db
bladePoints=>turbineArray(i) % bladePoints
bladeRadius=>turbineArray(i) % bladeRadius
solidity=>turbineArray(i) % solidity
azimuth=>turbineArray(i) % azimuth
rotSpeed=>turbineArray(i) % rotSpeed
towerShaftIntersect=>turbineArray(i) % towerShaftIntersect
baseLocation=>turbineArray(i) % baseLocation
uvShaft=>turbineArray(i) % uvShaft
uvTower=>turbineArray(i) % uvTower
rotorApex=>turbineArray(i) % rotorApex
uvShaftDir=>turbineArray(i) % uvShaftDir
nacYaw=>turbineArray(i) % nacYaw
numAnnulusSections=>turbineArray(i) % numAnnulusSections
annulusSectionAngle=>turbineArray(i) % annulusSectionAngle


! Assign Pointers turbineModel (j)
ShftTilt=>turbineModel(j) % ShftTilt
preCone=>turbineModel(j) % preCone
TowerHt=>turbineModel(j) % TowerHt
Twr2Shft=> turbineModel(j) % Twr2Shft
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
TipRad=>turbineModel(j) % TipRad
HubRad=>turbineModel(j) % HubRad
PreCone=>turbineModel(j) %PreCone

!!-- Do all proper conversions for the required variables
! Convert nacelle yaw from compass directions to the standard convention
!~ call atm_compassToStandard(nacYaw)

! The nacelle Yaw is set to 0 deg in the streamwise direction
write(*,*) 'NacYaw is ', nacYaw
! Turbine specific
azimuth = degRad * azimuth
rotSpeed = rpmRadSec * rotSpeed
nacYaw = degRad * nacYaw

! Turbine model specific
shftTilt = degRad * shftTilt
preCone = degRad * preCone

! Calculate tower shaft intersection and rotor apex locations. (The i-index is
! at the turbine array level for each turbine and the j-index is for each type
! of turbine--if all turbines are the same, j- is always 0.)  The rotor apex is
! not yet rotated for initial yaw that is done below.
towerShaftIntersect = turbineArray(i) % baseLocation
towerShaftIntersect(3) = towerShaftIntersect(3) + TowerHt + Twr2Shft
rotorApex = towerShaftIntersect
rotorApex(1) = rotorApex(1) +  (OverHang + UndSling) * cos(ShftTilt)
rotorApex(3) = rotorApex(3) +  (OverHang + UndSling) * sin(ShftTilt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Create Nacelle Point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
turbineArray(i) % nacelleLocation = rotorApex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Create the first set of actuator points                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the shaft direction from tilt, not from the overhang offset.  The
! previous offset-based construction is singular for OverHang + UndSling = 0.
uvShaft(1) = cos(ShftTilt)
uvShaft(2) = 0.0_rprec
uvShaft(3) = sin(ShftTilt)

! Define vector aligned with the tower pointing from the ground to the nacelle
uvTower = vector_add(towerShaftIntersect, - baseLocation)
uvTower = vector_divide( uvTower, vector_mag(uvTower))

! Define thickness of each blade section
do k=1, numBladePoints
    db(k) = (TipRad - HubRad)/(numBladePoints)
enddo

! This creates the first set of points
do k=1, numBl
    root = rotorApex
    beta = PreCone - ShftTilt
    root(1)= root(1) + HubRad*sin(beta)
    root(3)= root(3) + HubRad*cos(beta)
!~     dist = HubRad
    dist = 0.

    ! Number of blade points for the first annular section
    do m=1, numBladePoints
        dist = dist + 0.5*(db(m))
        bladePoints(k,1,m,1) = root(1) + dist*sin(beta)
        bladePoints(k,1,m,2) = root(2)
        bladePoints(k,1,m,3) = root(3) + dist*cos(beta)
        do n=1,numAnnulusSections
!~             bladeRadius(k,n,m) = dist
            bladeRadius(k,n,m) = dist + HubRad
            solidity(k,n,m)=1./numAnnulusSections
        enddo
        dist = dist + 0.5*db(m)
    enddo

    ! If there are more than one blade create the points of other blades by
    ! rotating the points of the first blade
    if (k > 1) then
        do m=1, numBladePoints
            bladePoints(k,1,m,:)=rotatePoint(bladePoints(k,1,m,:), rotorApex,  &
            uvShaft,(360.0/NumBl)*(k-1)*degRad)
        enddo
    endif

    ! Rotate points for all the annular sections
    if (numAnnulusSections .lt. 2) cycle ! Cycle if only one section (ALM)
    do n=2, numAnnulusSections
        do m=1, numBladePoints
            bladePoints(k,n,m,:) =                                             &
            rotatePoint(bladePoints(k,1,m,:), rotorApex,                       &
            uvShaft,(annulusSectionAngle/(numAnnulusSections))*(n-1.)*degRad)
        enddo
    enddo
enddo

! Apply the first rotation
turbineArray(i) % bladePoints_rigid = turbineArray(i) % bladePoints
turbineArray(i) % deltaAzimuth = azimuth
call atm_rotateBlades(i)

end subroutine atm_create_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_update(i, dt, time)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
real(rprec), intent(in) :: dt                            ! Time step
real(rprec), intent(in) :: time                          ! Simulation time

! Rotate the blades
call atm_computeRotorSpeed(i,dt)
call atm_rotateBlades(i)

call atm_control_yaw(i, time)

if (atm_structure_enabled()) then
    call atm_structure_diag_snapshot(i, 'pre_solve', dt, time)
    call atm_solve_structure(i, dt)
    call atm_structure_diag_snapshot(i, 'post_solve', dt, time)
endif

!~ if(pastFirstTimeStep) then
    ! Compute the lift correction for this case
!~     call atm_compute_cl_correction(i)
!~ endif

end subroutine atm_update

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_control_yaw(i, time)
! This subroutine updates the model each time-step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
real(rprec), intent(in) :: time                          ! Simulation time

integer :: j                                             ! Turbine Type ID

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Will calculate the yaw angle  and yaw the nacelle (from degrees to radians)
if ( turbineModel(j) % YawControllerType == "timeYawTable" ) then
    turbineArray(i) % deltaNacYaw = interpolate(time,                          &
    turbineModel(j) % yaw_time(:), turbineModel(j) % yaw_angle(:)) * degRad -  &
    turbineArray(i) % NacYaw

    ! Yaw only if angle is greater than given tolerance
    if (abs(turbineArray(i) % deltaNacYaw) > 0.00000001) then
        call atm_yawNacelle(i)
    endif

!~     write(*,*) 'Delta Yaw is', turbineArray(i) % deltaNacYaw/degRad
!~     write(*,*) 'Nacelle Yaw is', turbineArray(i) % NacYaw/degRad
endif

end subroutine atm_control_yaw

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeRotorSpeed(i,dt)
! This subroutine rotates the turbine blades
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                     ! Turbine number
real(rprec), intent(in) :: dt                ! time step
integer :: j                                 ! Turbine type
!real(rprec) :: deltaAzimuth                  ! Angle of rotation

! Pointers to turbineArray(i)
real(rprec), pointer :: rotSpeed, torqueGen, torqueRotor, fluidDensity

! Pointers to turbineModel(j)
real(rprec), pointer :: GBRatio, CutInGenSpeed, RatedGenSpeed
real(rprec), pointer :: Region2StartGenSpeed, Region2EndGenSpeed
real(rprec), pointer :: CutInGenTorque,RateLimitGenTorque,RatedGenTorque
real(rprec), pointer :: KGen,TorqueControllerRelax, DriveTrainIner

! Other variables to be used
real(rprec) :: torqueGenOld, genSpeed, dGenSpeed, Region2StartGenTorque
real(rprec) :: torqueSlope, Region2EndGenTorque

! Pitch Controller values
real(rprec) :: GK, KP, KI, SpeedError
real(rprec), pointer :: IntSpeedError, PitchControlAngle

j=turbineArray(i) % turbineTypeID

rotSpeed=>turbineArray(i) % rotSpeed
torqueGen=>turbineArray(i) % torqueGen
torqueRotor => turbineArray(i) % torqueRotor
fluidDensity => turbineArray(i) % fluidDensity

GBRatio => turbineModel(j) % GBRatio
CutInGenSpeed => turbineModel(j) % CutInGenSpeed
CutInGenTorque => turbineModel(j) % CutInGenTorque
Region2StartGenSpeed => turbineModel(j) % Region2StartGenSpeed
KGen => turbineModel(j) % KGen
RatedGenSpeed => turbineModel(j) % RatedGenSpeed
Region2EndGenSpeed => turbineModel(j) % Region2EndGenSpeed
RatedGenTorque => turbineModel(j) % RatedGenTorque
RateLimitGenTorque => turbineModel(j) % RateLimitGenTorque
TorqueControllerRelax => turbineModel(j) % TorqueControllerRelax
DriveTrainIner => turbineModel(j) % DriveTrainIner

IntSpeedError => turbineArray(i) % IntSpeedError
PitchControlAngle => turbineArray(i) % PitchControlAngle

    ! No torque controller option
    if (turbineModel(j) % TorqueControllerType == "none") then

    elseif (turbineModel(j) % TorqueControllerType == "fiveRegion") then

        ! Get the generator speed.
        genSpeed = (rotSpeed/rpmRadSec)*GBRatio

        ! Save the generator torque from the last time step.
        torqueGenOld = torqueGen

        ! Region 1.
        if (genSpeed < CutInGenSpeed) then

            torqueGen = CutInGenTorque

        ! Region 1-1/2.
        elseif ((genSpeed >= CutInGenSpeed) .and.                              &
               (genSpeed < Region2StartGenSpeed)) then

        dGenSpeed = genSpeed - CutInGenSpeed
        Region2StartGenTorque = KGen * Region2StartGenSpeed *                  &
                                       Region2StartGenSpeed
        torqueSlope = (Region2StartGenTorque - CutInGenTorque) /               &
                      ( Region2StartGenSpeed - CutInGenSpeed )
        torqueGen = CutInGenTorque + torqueSlope*dGenSpeed

        ! Region 2.
        elseif ((genSpeed >= Region2StartGenSpeed) .and.                       &
                 (genSpeed < Region2EndGenSpeed)) then

                torqueGen = KGen * genSpeed * genSpeed

        ! Region 2-1/2.
        elseif ((genSpeed >= Region2EndGenSpeed) .and.                         &
                 (genSpeed < RatedGenSpeed)) then

                dGenSpeed = genSpeed - Region2EndGenSpeed
                Region2EndGenTorque = KGen * Region2EndGenSpeed *              &
                                             Region2EndGenSpeed
                torqueSlope = (RatedGenTorque - Region2EndGenTorque) /         &
                              ( RatedGenSpeed - Region2EndGenSpeed )
                torqueGen = Region2EndGenTorque + torqueSlope*dGenSpeed

        ! Region 3.
        elseif (genSpeed >= RatedGenSpeed) then

                torqueGen = RatedGenTorque
        endif

        ! Limit the change in generator torque if after first time step
        ! (otherwise it slowly ramps up from its zero initialized value--we
        ! want it to instantly be at its desired value on the first time
        ! step, but smoothly vary from there).
        if ((abs((torqueGen - torqueGenOld)/dt) > RateLimitGenTorque) &
              .and. (pastFirstTimeStep)) then

            if (torqueGen > torqueGenOld) then

                torqueGen = torqueGenOld + (RateLimitGenTorque * dt);

            elseif (torqueGen <= torqueGenOld) then

                torqueGen = torqueGenOld - (RateLimitGenTorque * dt);
            endif
        endif

        ! Update the rotor speed.
        rotSpeed = rotSpeed + TorqueControllerRelax * (dt/DriveTrainIner) *       &
                              (torqueRotor*fluidDensity - GBRatio*torqueGen)

        if (turbineModel(j) % PitchControllerType == "none") then
            ! Limit the rotor speed to be positive and such that the generator
            !does not turn faster than rated.
            rotSpeed = max(0.0_rprec,rotSpeed)
            rotSpeed = min(rotSpeed,(RatedGenSpeed*rpmRadSec)/GBRatio)
        endif

    ! Torque control for fixed tip speed ratio
    ! Note that this current method does NOT support Coning in the rotor
    elseif (turbineModel(j) % TorqueControllerType == "fixedTSR") then

        if (pastFirstTimeStep) then
            ! Integrate the velocity along all actuator points
            call atm_integrate_u(i)

            ! Match the rotor speed to a given TSR
            rotSpeed = turbineArray(i) % u_infinity_mean *     &
                       turbineArray(i) % TSR / turbineModel(j) % tipRad

            ! Important to get rid of negative values
            rotSpeed = max(0.0_rprec,rotSpeed)
        endif
    endif

    ! Pitch controllers (If there's no pitch controller, then don't do anything)
    if (turbineModel(j) % PitchControllerType == "gainScheduledPI") then

       ! Get the generator speed.
       genSpeed = (rotSpeed/rpmRadSec)*GBRatio

       ! Calculate the gain
       GK =  1.0/(1.0 + PitchControlAngle/turbineModel(j) % PitchControlAngleK)

       ! Calculate the Proportional and Integral terms
       KP = GK*turbineModel(j) % PitchControlKP0
       KI = GK*turbineModel(j) % PitchControlKI0

       ! Get speed error (generator in rpm) and update integral
       ! Integral is saturated to not push the angle beyond its limits
       SpeedError = genSpeed - RatedGenSpeed
       !write(*,*) 'Speed Error is: ', speedError

       IntSpeedError = IntSpeedError + SpeedError*dt
       IntSpeedError = min( max(IntSpeedError,                                 &
                       turbineModel(j) % PitchControlAngleMin/KI),             &
                       turbineModel(j) % PitchControlAngleMax/KI)

       ! Apply PI controller and saturate
       PitchControlAngle = KP*SpeedError + KI*IntSpeedError
       PitchControlAngle = min( max( PitchControlAngle,                        &
                           turbinemodel(j) % PitchControlAngleMin),            &
                           turbineModel(j) % PitchControlAngleMax)

    endif

    ! Compute the change in blade position at new rotor speed.
    turbineArray(i) % deltaAzimuth = rotSpeed * dt

end subroutine atm_computeRotorSpeed

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_rotateBlades(i)
! This subroutine rotates the turbine blades
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i                                 ! Turbine number
!real(rprec), intent(in) :: dt                            ! time step
integer :: j                                 ! Turbine type
integer :: m, n, q                           ! Counters tu be used in do loops
real(rprec) :: deltaAzimuth, deltaAzimuthI   ! Angle of rotation
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: rotSpeed
real(rprec), pointer :: uvShaft(:)
real(rprec), pointer :: azimuth
real(rprec) :: disp_ax, disp_tg, theta_tot
real(rprec) :: vec1_rot(3), vec2_rot(3), origin_zero(3)
#ifdef ENABLE_CUDA
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladePoints_d
real(rprec) :: rotorApex1, rotorApex2, rotorApex3
real(rprec) :: ax1, ax2, ax3, cang, sang
real(rprec) :: rm11, rm12, rm13, rm21, rm22, rm23, rm31, rm32, rm33
real(rprec) :: p1, p2, p3, r1, r2, r3
#endif

j=turbineArray(i) % turbineTypeID

! Variables which are used by pointers
rotorApex=> turbineArray(i) % rotorApex
rotSpeed=>turbineArray(i) % rotSpeed
uvShaft=>turbineArray(i) % uvShaft
azimuth=>turbineArray(i) % azimuth

! Angle of rotation
deltaAzimuth = turbineArray(i) % deltaAzimuth

! Check the rotation direction first and set the local delta azimuth
! variable accordingly.
if (turbineArray(i) % rotationDir == "cw") then
    deltaAzimuthI = deltaAzimuth
else if (turbineArray(i) % rotationDir == "ccw") then
    deltaAzimuthI = -deltaAzimuth
endif

#ifdef ENABLE_CUDA
if (atm_model_cuda_enabled() .and. .not. atm_structure_enabled()) then
    bladePoints_d => turbineArray(i) % bladePoints
    rotorApex1 = rotorApex(1)
    rotorApex2 = rotorApex(2)
    rotorApex3 = rotorApex(3)
    ax1 = uvShaft(1)
    ax2 = uvShaft(2)
    ax3 = uvShaft(3)
    cang = cos(deltaAzimuthI)
    sang = sin(deltaAzimuthI)

    rm11 = ax1*ax1 + (1._rprec - ax1*ax1) * cang
    rm12 = ax1*ax2 * (1._rprec - cang) - ax3 * sang
    rm13 = ax1*ax3 * (1._rprec - cang) + ax2 * sang
    rm21 = ax1*ax2 * (1._rprec - cang) + ax3 * sang
    rm22 = ax2*ax2 + (1._rprec - ax2*ax2) * cang
    rm23 = ax2*ax3 * (1._rprec - cang) - ax1 * sang
    rm31 = ax1*ax3 * (1._rprec - cang) - ax2 * sang
    rm32 = ax2*ax3 * (1._rprec - cang) + ax1 * sang
    rm33 = ax3*ax3 + (1._rprec - ax3*ax3) * cang

    !$cuf kernel do(3) <<<*,*>>>
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) % numAnnulusSections
            do m=1, turbineModel(j) % numBl
                p1 = bladePoints_d(m,n,q,1) - rotorApex1
                p2 = bladePoints_d(m,n,q,2) - rotorApex2
                p3 = bladePoints_d(m,n,q,3) - rotorApex3

                r1 = rm11*p1 + rm12*p2 + rm13*p3
                r2 = rm21*p1 + rm22*p2 + rm23*p3
                r3 = rm31*p1 + rm32*p2 + rm33*p3

                bladePoints_d(m,n,q,1) = r1 + rotorApex1
                bladePoints_d(m,n,q,2) = r2 + rotorApex2
                bladePoints_d(m,n,q,3) = r3 + rotorApex3
            enddo
        enddo
    enddo

    call atm_model_cuda_check('atm_rotateBlades')
else
#endif
origin_zero = (/ 0._rprec, 0._rprec, 0._rprec /)

! Loop through all the points and rotate them accordingly
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) % numAnnulusSections
        do m=1, turbineModel(j) % numBl
            if (atm_structure_enabled()) then
                turbineArray(i) % bladePoints_rigid(m,n,q,:) = rotatePoint(    &
                turbineArray(i) % bladePoints_rigid(m,n,q,:), rotorApex,       &
                uvShaft, deltaAzimuthI)

                vec1_rot = rotatePoint(                                        &
                    turbineArray(i) % bladeAlignedVectors(m,n,q,1,:),          &
                    origin_zero, uvShaft, deltaAzimuthI)
                vec2_rot = rotatePoint(                                        &
                    turbineArray(i) % bladeAlignedVectors(m,n,q,2,:),          &
                    origin_zero, uvShaft, deltaAzimuthI)

                theta_tot = (turbineArray(i)%twistAng(m,n,q) +                 &
                    turbineArray(i)%Pitch + turbineArray(i)%PitchControlAngle)  &
                    * degRad + turbineArray(i)%elastic_twist(m,1,q)
                disp_ax = turbineArray(i)%flap_disp(m,1,q) * cos(theta_tot) -  &
                          turbineArray(i)%edge_disp(m,1,q) * sin(theta_tot)
                disp_tg = turbineArray(i)%flap_disp(m,1,q) * sin(theta_tot) +  &
                          turbineArray(i)%edge_disp(m,1,q) * cos(theta_tot)

                turbineArray(i) % bladePoints(m,n,q,:) =                       &
                    turbineArray(i) % bladePoints_rigid(m,n,q,:) +             &
                    disp_ax * vec1_rot + disp_tg * vec2_rot
            else
                turbineArray(i) % bladePoints(m,n,q,:) = rotatePoint(          &
                turbineArray(i) % bladePoints(m,n,q,:), rotorApex, uvShaft,    &
                deltaAzimuthI)
            endif
        enddo
    enddo
enddo
#ifdef ENABLE_CUDA
endif
#endif

if (pastFirstTimeStep) then
    azimuth = azimuth + deltaAzimuth;
        if (azimuth .ge. 2.0 * pi) then
            azimuth =azimuth - 2.0 *pi;
        endif
endif

end subroutine atm_rotateBlades

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) function atm_panel_antiderivative(u) result(F)
! Smooth antiderivative used by the exact panel induced-velocity integral.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

real(rprec), intent(in) :: u
real(rprec) :: u2
real(rprec), parameter :: u_small = 1.0e-2_rprec

u2 = u * u
if (abs(u) < u_small) then
    F = 0.5_rprec * u *                                                     &
        (1._rprec - 0.5_rprec * u2 *                                        &
        (1._rprec - (u2 / 3._rprec) * (1._rprec - 0.25_rprec * u2)))
else
    F = (1._rprec - exp(-u2)) / (2._rprec * u)
endif

end function atm_panel_antiderivative

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compute_cl_correction(i)
! This subroutine computes the induced-velocity correction for the ALM kernel.
! The CPU path uses the exact panel integral used by the GPU 2024_panel path.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, intent(in) :: i             ! Turbine number
integer :: j                         ! Turbine type
integer :: m, n, q, k                ! Counters used in do loops
integer :: numBP, numBl_j, numAS     ! Cached loop bounds
real(rprec) :: eps_s                 ! The epsilon value in the simulation
real(rprec) :: f                     ! The relaxation factor
real(rprec) :: inv_eps_s, inv_2pi_eps_s, inv_4pi_sqrt
real(rprec) :: mag_uinf, inv_mag_tan
real(rprec) :: half_dzp, up, um, uy_les_k, uy_opt_k
real(rprec) :: acc_les1, acc_les2, acc_opt1, acc_opt2
real(rprec) :: z_q, g_over_v_k
real(rprec), allocatable :: z_loc(:), g_over_v_loc(:), eps_opt_loc(:)
real(rprec), allocatable :: db_loc(:), tan1_loc(:), tan2_loc(:)
real(rprec), allocatable :: inv_2pi_eps_opt_loc(:), z_face(:), f_face_les(:)

real(rprec), pointer :: uy_opt_vec(:,:,:,:)
real(rprec), pointer :: uy_LES_vec(:,:,:,:)
real(rprec), pointer :: ux_LES_vec(:,:,:,:)
real(rprec), pointer :: Uinf_vec(:,:,:,:)
real(rprec), pointer :: windVectors(:,:,:,:)
real(rprec), pointer :: cl(:,:,:)
real(rprec), pointer :: cd(:,:,:)
real(rprec), pointer :: Vmag(:,:,:)
real(rprec), pointer :: chord(:,:,:)

#ifdef ENABLE_CUDA
if (atm_model_cuda_enabled() .and. .not. atm_structure_enabled()) then
    call atm_compute_cl_correction_gpu(i)
    return
endif
#endif

! The wind vector
uy_opt_vec => turbineArray(i) % uy_opt_vec
uy_LES_vec => turbineArray(i) % uy_LES_vec
Uinf_vec => turbineArray(i) % Uinf_vec
windVectors => turbineArray(i) % windVectors
cl => turbineArray(i) % cl
cd => turbineArray(i) % cd
Vmag => turbineArray(i) % Vmag
chord => turbineArray(i) % chord
ux_LES_vec => turbineArray(i) % ux_LES_vec

j = turbineArray(i) % turbineTypeID
numBP = turbineArray(i) % numBladePoints
numAS = turbineArray(i) % numAnnulusSections
numBl_j = turbineModel(j) % numBl
eps_s = turbineArray(i) % epsilon
f = 0.1_rprec
inv_eps_s = 1._rprec / eps_s
inv_2pi_eps_s = inv_eps_s / (2._rprec * pi)
inv_4pi_sqrt = 1._rprec / (4._rprec * sqrt(pi))

allocate(z_loc(numBP), g_over_v_loc(numBP), eps_opt_loc(numBP))
allocate(db_loc(numBP), tan1_loc(numBP), tan2_loc(numBP))
allocate(inv_2pi_eps_opt_loc(numBP), z_face(numBP+1), f_face_les(numBP+1))

do q = 1, numBP
    do n = 1, numAS
        do m = 1, numBl_j
            Uinf_vec(m,n,q,:) = windVectors(m,n,q,:) - uy_opt_vec(m,n,q,:)
            mag_uinf = vector_mag(Uinf_vec(m,n,q,:))

            ux_LES_vec(m,n,q,:) = cd(m,n,q) * chord(m,n,q) / eps_s *           &
                inv_4pi_sqrt * Vmag(m,n,q) * Uinf_vec(m,n,q,:) / mag_uinf

            turbineArray(i) % G(m,n,q) = 0.5_rprec * cl(m,n,q) *              &
                chord(m,n,q) * Vmag(m,n,q)**2
            turbineArray(i) % epsilon_opt(m,n,q) = chord(m,n,q) *             &
                turbineArray(i) % optimalEpsilonChord
        enddo
    enddo
enddo

! Keep dG populated for diagnostics and backward-compatible output files.
do m = 1, numBl_j
    do n = 1, numAS
        turbineArray(i) % dG(m,n,1) = turbineArray(i) % G(m,n,1)
        do q = 2, numBP - 1
            turbineArray(i) % dG(m,n,q) =                                      &
                (turbineArray(i) % G(m,n,q+1) -                                &
                 turbineArray(i) % G(m,n,q-1)) * 0.5_rprec
        enddo
        turbineArray(i) % dG(m,n,numBP) = -turbineArray(i) % G(m,n,numBP)
    enddo
enddo

do m = 1, numBl_j
    do n = 1, numAS
        do k = 1, numBP
            z_loc(k) = turbineArray(i) % bladeRadius(m,n,k)
            g_over_v_loc(k) = turbineArray(i) % G(m,n,k) / Vmag(m,n,k)
            eps_opt_loc(k) = turbineArray(i) % epsilon_opt(m,n,k)
            db_loc(k) = turbineArray(i) % db(k)
            inv_2pi_eps_opt_loc(k) = 1._rprec / (2._rprec * pi * eps_opt_loc(k))

            inv_mag_tan = 1._rprec / sqrt(Uinf_vec(m,n,k,1)**2 +               &
                                           Uinf_vec(m,n,k,2)**2)
            tan1_loc(k) =  Uinf_vec(m,n,k,2) * inv_mag_tan
            tan2_loc(k) = -Uinf_vec(m,n,k,1) * inv_mag_tan
        enddo

        z_face(1) = z_loc(1) - 0.5_rprec * db_loc(1)
        do k = 1, numBP
            z_face(k+1) = z_loc(k) + 0.5_rprec * db_loc(k)
        enddo

        do q = 1, numBP
            acc_les1 = 0._rprec
            acc_les2 = 0._rprec
            acc_opt1 = 0._rprec
            acc_opt2 = 0._rprec
            z_q = z_loc(q)

            do k = 1, numBP + 1
                f_face_les(k) = atm_panel_antiderivative(                       &
                    (z_q - z_face(k)) * inv_eps_s)
            enddo

            do k = 1, numBP
                g_over_v_k = g_over_v_loc(k)
                uy_les_k = -g_over_v_k * inv_2pi_eps_s *                      &
                    (f_face_les(k) - f_face_les(k+1))

                half_dzp = 0.5_rprec * db_loc(k)
                up = (z_q - z_loc(k) + half_dzp) / eps_opt_loc(k)
                um = (z_q - z_loc(k) - half_dzp) / eps_opt_loc(k)
                uy_opt_k = -g_over_v_k * inv_2pi_eps_opt_loc(k) *             &
                    (atm_panel_antiderivative(up) - atm_panel_antiderivative(um))

                acc_les1 = acc_les1 + uy_les_k * tan1_loc(k)
                acc_les2 = acc_les2 + uy_les_k * tan2_loc(k)
                acc_opt1 = acc_opt1 + uy_opt_k * tan1_loc(k)
                acc_opt2 = acc_opt2 + uy_opt_k * tan2_loc(k)
            enddo

            uy_LES_vec(m,n,q,1) = acc_les1
            uy_LES_vec(m,n,q,2) = acc_les2
            uy_LES_vec(m,n,q,3) = 0._rprec
            uy_opt_vec(m,n,q,1) = acc_opt1
            uy_opt_vec(m,n,q,2) = acc_opt2
            uy_opt_vec(m,n,q,3) = 0._rprec

            turbineArray(i) % uy_LES(m,n,q) = sqrt(acc_les1*acc_les1 +         &
                                                    acc_les2*acc_les2)
            turbineArray(i) % uy_opt(m,n,q) = sqrt(acc_opt1*acc_opt1 +         &
                                                    acc_opt2*acc_opt2)

            turbineArray(i) % du(m,n,q,:) = turbineArray(i) % du(m,n,q,:) *    &
                (1._rprec - f) + f * (uy_opt_vec(m,n,q,:) -                   &
                uy_LES_vec(m,n,q,:) + ux_LES_vec(m,n,q,:))
        enddo
    enddo
enddo

deallocate(z_loc, g_over_v_loc, eps_opt_loc, db_loc, tan1_loc, tan2_loc)
deallocate(inv_2pi_eps_opt_loc, z_face, f_face_les)

end subroutine atm_compute_cl_correction

#ifdef ENABLE_CUDA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compute_cl_correction_gpu(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! GPU induced-velocity correction.
!
! Runtime selector:
!   LESGO_ATM_INDUCED_METHOD=legacy        old dG/dz filtered lifting-line path
!   LESGO_ATM_INDUCED_METHOD=2024_midpoint generalized 2024 midpoint quadrature
!   LESGO_ATM_INDUCED_METHOD=2024_panel    generalized 2024 exact panel integral
!
! The selector changes only the induced-velocity convolution.  Structural solver,
! force application, turbine controls, and LESGO flow numerics are unchanged.
! Update docs/environment_switches.md if these model switches change.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer, parameter :: induced_legacy = 0
integer, parameter :: induced_midpoint = 1
integer, parameter :: induced_panel = 2
integer, parameter :: option_disabled = 0
integer, parameter :: option_enabled = 1
integer :: j, m, n, q, k, mmend, nnend, qqend
integer, save :: induced_method = -1
integer, save :: du_include_ux = -1
integer, save :: ux_use_cd_eff = -1
logical, save :: induced_method_printed = .false.
character(len=32) :: env_value
logical :: env_has_value
real(rprec) :: eps_s, opt_eps_chord, relax, sqrt_pi
real(rprec) :: u1, u2, u3, o1, o2, o3, mag_uinf, mag_uyopt
real(rprec) :: dir1, dir2, dir3, cd_eff, fac, vmag_i
real(rprec) :: z_q, z_k, dzeta, dG_k, uy_les_s, uy_opt_s
real(rprec) :: tan1, tan2, tan3, mag_tan
real(rprec) :: les1, les2, les3, opt1, opt2, opt3, inv_vmag
real(rprec) :: eps_k, dzp_k, g_over_v, aexp, kernel_s
real(rprec) :: half_dzp, up, um, f_up, f_um, uval, u2val
real(rprec), parameter :: u_small = 1.0e-2_rprec
real(rprec), managed, pointer, dimension(:,:,:,:) :: uy_opt_vec
real(rprec), managed, pointer, dimension(:,:,:,:) :: uy_LES_vec
real(rprec), managed, pointer, dimension(:,:,:,:) :: ux_LES_vec
real(rprec), managed, pointer, dimension(:,:,:,:) :: Uinf_vec
real(rprec), managed, pointer, dimension(:,:,:,:) :: windVectors
real(rprec), managed, pointer, dimension(:,:,:) :: cl, cd, Vmag, chord
real(rprec), managed, pointer, dimension(:,:,:) :: G, dG
real(rprec), managed, pointer, dimension(:,:,:) :: epsilon_opt
real(rprec), managed, pointer, dimension(:,:,:) :: bladeRadius
real(rprec), managed, pointer, dimension(:,:,:) :: uy_LES, uy_opt
real(rprec), managed, pointer, dimension(:,:,:,:) :: du
real(rprec), managed, pointer, dimension(:) :: db

if (induced_method < 0) then
    induced_method = induced_panel
    call atm_model_env_token('LESGO_ATM_INDUCED_METHOD', env_value,            &
        env_has_value)
    if (env_has_value) then
        select case (env_value)
        case ('legacy','LEGACY','0','old','OLD')
            induced_method = induced_legacy
        case ('2024_uniform','2024-UNIFORM','uniform','UNIFORM',               &
              '2024_midpoint','2024-MIDPOINT','midpoint','MIDPOINT','1')
            induced_method = induced_midpoint
        ! Keep the personal-name token as a backward-compatible input alias;
        ! logs and documentation use the technical `2024_panel` method name.
        case ('2024_panel','2024-PANEL','panel','PANEL','atharva','ATHARVA',   &
              'exact','EXACT','2')
            induced_method = induced_panel
        case default
            induced_method = induced_panel
        end select
    endif
endif

if (du_include_ux < 0) then
    du_include_ux = option_enabled
    call atm_model_env_token('LESGO_ATM_DU_INCLUDE_UX', env_value,             &
        env_has_value)
    if (env_has_value) then
        if (env_value == '0') du_include_ux = option_disabled
        if (env_value == '1') du_include_ux = option_enabled
    endif
endif

if (ux_use_cd_eff < 0) then
    ux_use_cd_eff = option_disabled
    call atm_model_env_token('LESGO_ATM_UX_USE_CD_EFF', env_value,             &
        env_has_value)
    if (env_has_value) then
        if (env_value == '0') ux_use_cd_eff = option_disabled
        if (env_value == '1') ux_use_cd_eff = option_enabled
    endif
endif

if (.not. induced_method_printed) then
    if (induced_method == induced_legacy) then
        write(*,*) 'ATM induced velocity method: legacy dG/dz'
    elseif (induced_method == induced_midpoint) then
        write(*,*) 'ATM induced velocity method: 2024 uniform reference'
    else
        write(*,*) 'ATM induced velocity method: 2024 generalized exact panel'
    endif
    write(*,*) 'ATM induced options: du_include_ux=', du_include_ux,           &
        ' ux_use_cd_eff=', ux_use_cd_eff
    induced_method_printed = .true.
endif

j = turbineArray(i) % turbineTypeID
mmend = turbineModel(j) % numBl
nnend = turbineArray(i) % numAnnulusSections
qqend = turbineArray(i) % numBladePoints
eps_s = turbineArray(i) % epsilon
opt_eps_chord = turbineArray(i) % optimalEpsilonChord
relax = 0.1_rprec
sqrt_pi = sqrt(pi)

uy_opt_vec => turbineArray(i) % uy_opt_vec
uy_LES_vec => turbineArray(i) % uy_LES_vec
ux_LES_vec => turbineArray(i) % ux_LES_vec
Uinf_vec => turbineArray(i) % Uinf_vec
windVectors => turbineArray(i) % windVectors
cl => turbineArray(i) % cl
cd => turbineArray(i) % cd
Vmag => turbineArray(i) % Vmag
chord => turbineArray(i) % chord
G => turbineArray(i) % G
dG => turbineArray(i) % dG
epsilon_opt => turbineArray(i) % epsilon_opt
bladeRadius => turbineArray(i) % bladeRadius
uy_LES => turbineArray(i) % uy_LES
uy_opt => turbineArray(i) % uy_opt
du => turbineArray(i) % du
db => turbineArray(i) % db

!$cuf kernel do(3) <<<*,*>>>
do q=1, qqend
    do n=1, nnend
        do m=1, mmend
            Uinf_vec(m,n,q,1) = windVectors(m,n,q,1) - uy_opt_vec(m,n,q,1)
            Uinf_vec(m,n,q,2) = windVectors(m,n,q,2) - uy_opt_vec(m,n,q,2)
            Uinf_vec(m,n,q,3) = windVectors(m,n,q,3) - uy_opt_vec(m,n,q,3)

            u1 = Uinf_vec(m,n,q,1)
            u2 = Uinf_vec(m,n,q,2)
            u3 = Uinf_vec(m,n,q,3)
            mag_uinf = sqrt(u1*u1 + u2*u2 + u3*u3)
            dir1 = u1 / mag_uinf
            dir2 = u2 / mag_uinf
            dir3 = u3 / mag_uinf

            o1 = uy_opt_vec(m,n,q,1)
            o2 = uy_opt_vec(m,n,q,2)
            o3 = uy_opt_vec(m,n,q,3)
            mag_uyopt = sqrt(o1*o1 + o2*o2 + o3*o3)
            vmag_i = Vmag(m,n,q)
            cd_eff = cd(m,n,q) + cl(m,n,q) * mag_uyopt / vmag_i
            if (ux_use_cd_eff == option_enabled) then
                fac = cd_eff * chord(m,n,q) / eps_s /                         &
                    (4._rprec * sqrt_pi) * vmag_i
            else
                fac = cd(m,n,q) * chord(m,n,q) / eps_s /                      &
                    (4._rprec * sqrt_pi) * vmag_i
            endif
            ux_LES_vec(m,n,q,1) = fac * dir1
            ux_LES_vec(m,n,q,2) = fac * dir2
            ux_LES_vec(m,n,q,3) = fac * dir3

            G(m,n,q) = 0.5_rprec * cl(m,n,q) * chord(m,n,q) *                 &
                vmag_i * vmag_i
            epsilon_opt(m,n,q) = chord(m,n,q) * opt_eps_chord
        enddo
    enddo
enddo

! Legacy method still needs dG/dz.  It is harmless to compute it for all modes
! and keeps the kernel sequence stable across A/B/C tests.
!$cuf kernel do(3) <<<*,*>>>
do q=1, qqend
    do n=1, nnend
        do m=1, mmend
            if (q .eq. 1) then
                dG(m,n,q) = G(m,n,1)
            elseif (q .eq. qqend) then
                dG(m,n,q) = -G(m,n,qqend)
            else
                dG(m,n,q) = (G(m,n,q+1) - G(m,n,q-1)) / 2._rprec
            endif
        enddo
    enddo
enddo

!$cuf kernel do(3) <<<*,*>>>
do q=1, qqend
    do n=1, nnend
        do m=1, mmend
            uy_LES_vec(m,n,q,1) = 0._rprec
            uy_LES_vec(m,n,q,2) = 0._rprec
            uy_LES_vec(m,n,q,3) = 0._rprec
            uy_opt_vec(m,n,q,1) = 0._rprec
            uy_opt_vec(m,n,q,2) = 0._rprec
            uy_opt_vec(m,n,q,3) = 0._rprec

            z_q = bladeRadius(m,n,q)
            les1 = 0._rprec
            les2 = 0._rprec
            les3 = 0._rprec
            opt1 = 0._rprec
            opt2 = 0._rprec
            opt3 = 0._rprec

            do k=1, qqend
                z_k = bladeRadius(m,n,k)
                dzeta = z_q - z_k
                eps_k = epsilon_opt(m,n,k)
                dzp_k = db(k)
                if (induced_method == induced_midpoint) then
                    g_over_v = G(m,n,k)
                else
                    g_over_v = G(m,n,k) / Vmag(m,n,k)
                endif

                tan1 = Uinf_vec(m,n,k,2)
                tan2 = -Uinf_vec(m,n,k,1)
                tan3 = 0._rprec
                mag_tan = sqrt(tan1*tan1 + tan2*tan2 + tan3*tan3)
                tan1 = tan1 / mag_tan
                tan2 = tan2 / mag_tan
                tan3 = tan3 / mag_tan

                if (induced_method == induced_legacy) then
                    if (k /= q) then
                        dG_k = dG(m,n,k)
                        uy_les_s = -dG_k / (4._rprec * pi * dzeta) *          &
                            (1._rprec - exp(-((dzeta/eps_s) *                 &
                            (dzeta/eps_s))))
                        uy_opt_s = -dG_k / (4._rprec * pi * dzeta) *          &
                            (1._rprec - exp(-((dzeta/eps_k) *                 &
                            (dzeta/eps_k))))
                    else
                        uy_les_s = 0._rprec
                        uy_opt_s = 0._rprec
                    endif
                elseif (induced_method == induced_midpoint) then
                    if (abs(dzeta) > 1.0e-14_rprec) then
                        aexp = exp(-((dzeta/eps_s) * (dzeta/eps_s)))
                        kernel_s = 1._rprec / (2._rprec * pi * eps_s*eps_s) * &
                            (aexp - (1._rprec - aexp) * eps_s*eps_s /         &
                            (2._rprec * dzeta*dzeta))
                        uy_les_s = -g_over_v * kernel_s * dzp_k

                        aexp = exp(-((dzeta/eps_k) * (dzeta/eps_k)))
                        kernel_s = 1._rprec / (2._rprec * pi * eps_k*eps_k) * &
                            (aexp - (1._rprec - aexp) * eps_k*eps_k /         &
                            (2._rprec * dzeta*dzeta))
                        uy_opt_s = -g_over_v * kernel_s * dzp_k
                    else
                        uy_les_s = -g_over_v / (4._rprec * pi * eps_s*eps_s) *&
                            dzp_k
                        uy_opt_s = -g_over_v / (4._rprec * pi * eps_k*eps_k) *&
                            dzp_k
                    endif
                else
                    half_dzp = 0.5_rprec * dzp_k

                    uval = (dzeta + half_dzp) / eps_s
                    u2val = uval * uval
                    if (abs(uval) < u_small) then
                        f_up = 0.5_rprec * uval *                             &
                            (1._rprec - 0.5_rprec * u2val *                   &
                            (1._rprec - (u2val / 3._rprec) *                  &
                            (1._rprec - 0.25_rprec * u2val)))
                    else
                        f_up = (1._rprec - exp(-u2val)) / (2._rprec * uval)
                    endif
                    uval = (dzeta - half_dzp) / eps_s
                    u2val = uval * uval
                    if (abs(uval) < u_small) then
                        f_um = 0.5_rprec * uval *                             &
                            (1._rprec - 0.5_rprec * u2val *                   &
                            (1._rprec - (u2val / 3._rprec) *                  &
                            (1._rprec - 0.25_rprec * u2val)))
                    else
                        f_um = (1._rprec - exp(-u2val)) / (2._rprec * uval)
                    endif
                    uy_les_s = -g_over_v * (f_up - f_um) /                    &
                        (2._rprec * pi * eps_s)

                    up = (dzeta + half_dzp) / eps_k
                    um = (dzeta - half_dzp) / eps_k
                    u2val = up * up
                    if (abs(up) < u_small) then
                        f_up = 0.5_rprec * up *                               &
                            (1._rprec - 0.5_rprec * u2val *                   &
                            (1._rprec - (u2val / 3._rprec) *                  &
                            (1._rprec - 0.25_rprec * u2val)))
                    else
                        f_up = (1._rprec - exp(-u2val)) / (2._rprec * up)
                    endif
                    u2val = um * um
                    if (abs(um) < u_small) then
                        f_um = 0.5_rprec * um *                               &
                            (1._rprec - 0.5_rprec * u2val *                   &
                            (1._rprec - (u2val / 3._rprec) *                  &
                            (1._rprec - 0.25_rprec * u2val)))
                    else
                        f_um = (1._rprec - exp(-u2val)) / (2._rprec * um)
                    endif
                    uy_opt_s = -g_over_v * (f_up - f_um) /                    &
                        (2._rprec * pi * eps_k)
                endif

                les1 = les1 + uy_les_s * tan1
                les2 = les2 + uy_les_s * tan2
                les3 = les3 + uy_les_s * tan3
                opt1 = opt1 + uy_opt_s * tan1
                opt2 = opt2 + uy_opt_s * tan2
                opt3 = opt3 + uy_opt_s * tan3
            enddo

            if (induced_method /= induced_panel) then
                inv_vmag = 1._rprec / Vmag(m,n,q)
                uy_LES_vec(m,n,q,1) = les1 * inv_vmag
                uy_LES_vec(m,n,q,2) = les2 * inv_vmag
                uy_LES_vec(m,n,q,3) = les3 * inv_vmag
                uy_opt_vec(m,n,q,1) = opt1 * inv_vmag
                uy_opt_vec(m,n,q,2) = opt2 * inv_vmag
                uy_opt_vec(m,n,q,3) = opt3 * inv_vmag
            else
                uy_LES_vec(m,n,q,1) = les1
                uy_LES_vec(m,n,q,2) = les2
                uy_LES_vec(m,n,q,3) = les3
                uy_opt_vec(m,n,q,1) = opt1
                uy_opt_vec(m,n,q,2) = opt2
                uy_opt_vec(m,n,q,3) = opt3
            endif

            uy_LES(m,n,q) = sqrt(uy_LES_vec(m,n,q,1)*uy_LES_vec(m,n,q,1) +    &
                uy_LES_vec(m,n,q,2)*uy_LES_vec(m,n,q,2) +                    &
                uy_LES_vec(m,n,q,3)*uy_LES_vec(m,n,q,3))
            uy_opt(m,n,q) = sqrt(uy_opt_vec(m,n,q,1)*uy_opt_vec(m,n,q,1) +    &
                uy_opt_vec(m,n,q,2)*uy_opt_vec(m,n,q,2) +                    &
                uy_opt_vec(m,n,q,3)*uy_opt_vec(m,n,q,3))

            if (du_include_ux == option_enabled) then
                du(m,n,q,1) = du(m,n,q,1) * (1._rprec - relax) + relax *      &
                    (uy_opt_vec(m,n,q,1) - uy_LES_vec(m,n,q,1) +             &
                     ux_LES_vec(m,n,q,1))
                du(m,n,q,2) = du(m,n,q,2) * (1._rprec - relax) + relax *      &
                    (uy_opt_vec(m,n,q,2) - uy_LES_vec(m,n,q,2) +             &
                     ux_LES_vec(m,n,q,2))
                du(m,n,q,3) = du(m,n,q,3) * (1._rprec - relax) + relax *      &
                    (uy_opt_vec(m,n,q,3) - uy_LES_vec(m,n,q,3) +             &
                     ux_LES_vec(m,n,q,3))
            else
                du(m,n,q,1) = du(m,n,q,1) * (1._rprec - relax) + relax *      &
                    (uy_opt_vec(m,n,q,1) - uy_LES_vec(m,n,q,1))
                du(m,n,q,2) = du(m,n,q,2) * (1._rprec - relax) + relax *      &
                    (uy_opt_vec(m,n,q,2) - uy_LES_vec(m,n,q,2))
                du(m,n,q,3) = du(m,n,q,3) * (1._rprec - relax) + relax *      &
                    (uy_opt_vec(m,n,q,3) - uy_LES_vec(m,n,q,3))
            endif
        enddo
    enddo
enddo

call atm_model_cuda_sync('atm_compute_cl_correction')

end subroutine atm_compute_cl_correction_gpu
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function s_fit(psi, psipp, epsilon)
! This is the fit for the filtered lifting line theory
! Martinez-Tossas and Meneveau 2018
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

real(rprec), intent(in) :: psi, psipp, epsilon  ! Inputs for the fit
real(rprec) :: a,b,c                 ! Correction coefficients
real(rprec) :: one  ! Value used for 1 in the sign function
real(rprec) :: diff  ! The difference between psi and psipp
real(rprec) :: s_fit ! The function to be computed and returned
real(rprec) :: f ! Intermediate function

! One stored as a variable to be used in the sign function
one=1.

! Constants for the Fredholm integral equation solution fit
! These work well for epsilon >= 0.25
a=0.029
b=-2./3.
c=0.357

! The difference between psi and psipp (used many times in the formulas)
diff = psi-psipp + 0.000000000001

! The function with the first fit
f = 1./(4.*pi*abs(diff)) * (1. - exp(-diff**2)) -                      &
            a * epsilon**b / (diff**2) * (1. - exp(-c*abs(diff)**3))

! The final fit for the solution
s_fit = - dsign(one, diff) * (1. - 0.25 * exp(-epsilon) *               &
        (1. - exp(-.2 * abs(psipp)))) * f

return

end function s_fit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_calculate_variables(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Calculates the variables of the model that need information from the input
! files. It runs after reading input information.
implicit none

integer, intent(in) :: i ! Indicates the turbine number
integer :: j ! Indicates the turbine type
integer :: m, n, q ! Looping indices
integer,     pointer :: NumSec  ! Number of sections in lookup table
real(rprec),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: projectionRadius, projectionRadiusNacelle
real(rprec), pointer :: sphereRadius
real(rprec), pointer :: OverHang
real(rprec), pointer :: UndSling
real(rprec), pointer :: TipRad
real(rprec), pointer :: PreCone
real(rprec) :: cl1, cl2   ! Variables used to compute the slope of lift curve

! Identifies the turbineModel being used
j=turbineArray(i) % turbineTypeID ! The type of turbine (eg. NREL5MW)

! Number of sections in lookup table
NumSec => turbineModel(j) % NumSec

! Pointers dependent on turbineArray (i)
bladeRadius=>turbineArray(i) % bladeRadius
projectionRadius=>turbineArray(i) % projectionRadius
projectionRadiusNacelle=>turbineArray(i) % projectionRadiusNacelle
sphereRadius=>turbineArray(i) % sphereRadius

! Pointers dependent on turbineType (j)
OverHang=>turbineModel(j) % OverHang
UndSling=>turbineModel(j) % UndSling
TipRad=>turbineModel(j) % TipRad
PreCone=>turbineModel(j) %PreCone

! First compute the radius of the force projection (to the radius where the
! projection is only 0.001 its maximum value - this seems to recover 99.9% of
! the total forces when integrated
projectionRadius= turbineArray(i) % epsilon * sqrt(log(1.0/0.001))
projectionRadiusNacelle= turbineArray(i) % nacelleEpsilon*sqrt(log(1.0/0.001))

sphereRadius=sqrt(((OverHang + UndSling) + TipRad*sin(PreCone))**2 &
+ (TipRad*cos(PreCone))**2) + projectionRadius


! Compute the optimum value of epsilon for each blade section
! And compute the lift coefficient slope
do m=1, turbineModel(j) % numBl
    do n=1, turbineArray(i) %  numAnnulusSections
        do q=1, turbineArray(i) % numBladePoints

            ! Interpolate quantities through section
            turbineArray(i) % twistAng(m,n,q) =                                &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % twist(1:NumSec) )

            turbineArray(i) % chord(m,n,q) =                                   &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % chord(1:NumSec) )

            turbineArray(i) % sectionType(m,n,q) =                             &
                                   interpolate_i(bladeRadius(m,n,q),           &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % sectionType(1:NumSec))
            call atm_airfoil_blend_info(j, bladeRadius(m,n,q),                 &
                 turbineArray(i) % sectionTypeBlendLo(m,n,q),                  &
                 turbineArray(i) % sectionTypeBlendHi(m,n,q),                  &
                 turbineArray(i) % sectionTypeBlendW(m,n,q))

            turbineArray(i) % EI_blade(m,n,q) =                                &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % EI(1:NumSec) )
            turbineArray(i) % EI_edge_blade(m,n,q) =                           &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % EdgStff(1:NumSec) )
            turbineArray(i) % GJ_blade(m,n,q) =                                &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % GJStff(1:NumSec) )
            turbineArray(i) % rho_blade(m,n,q) =                               &
                                   interpolate(bladeRadius(m,n,q),             &
                                   turbineModel(j) % radius(1:NumSec),         &
                                   turbineModel(j) % rho(1:NumSec) )

            ! Compute lift coefficient slope
            ! 2 and 6 degrees
            cl1 = atm_blended_airfoil_coeff(j, bladeRadius(m,n,q),             &
                                            2._rprec, 1)
            cl2 = atm_blended_airfoil_coeff(j, bladeRadius(m,n,q),             &
                                            6._rprec, 1)

            ! Slope of the lift curve d Cl / d alpha
            ! Computed between 2 and 6 degrees
            turbineArray(i) % dCldalpha(m,n,q) = (cl2-cl1)/(4._rprec * pi/180.)

!~             write(*,*) 'Dcldalpha = ', turbineArray(i) % dCldalpha(m,n,q)

        enddo
    enddo
enddo


!~ ! Compute the lift correction for this case
!~ call atm_compute_cl_correction(i)

end subroutine atm_calculate_variables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) function atm_airfoil_component(foil, idx, coeff_kind) result(value)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(airfoilType_t), intent(in) :: foil
integer, intent(in) :: idx, coeff_kind

select case (coeff_kind)
case (1)
    value = foil % cl(idx)
case (2)
    value = foil % cd(idx)
case (3)
    value = foil % cm(idx)
case default
    value = 0._rprec
end select

end function atm_airfoil_component

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) function atm_airfoil_polar_value(foil, alpha_deg, coeff_kind)       &
    result(value)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

type(airfoilType_t), intent(in) :: foil
real(rprec), intent(in) :: alpha_deg
integer, intent(in) :: coeff_kind
integer :: idx, npts
real(rprec) :: xa, xb, ya, yb

npts = foil % n
if (npts <= 0) then
    value = 0._rprec
elseif (alpha_deg <= foil % AOA(1)) then
    value = atm_airfoil_component(foil, 1, coeff_kind)
elseif (alpha_deg >= foil % AOA(npts)) then
    value = atm_airfoil_component(foil, npts, coeff_kind)
else
    value = atm_airfoil_component(foil, npts, coeff_kind)
    do idx = 2, npts
        if (alpha_deg <= foil % AOA(idx)) then
            xa = foil % AOA(idx-1)
            xb = foil % AOA(idx)
            ya = atm_airfoil_component(foil, idx-1, coeff_kind)
            yb = atm_airfoil_component(foil, idx, coeff_kind)
            value = ya + (yb - ya) * (alpha_deg - xa) / (xb - xa)
            exit
        endif
    enddo
endif

end function atm_airfoil_polar_value

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_airfoil_blend_info(j, radius, foil_lo, foil_hi, weight)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: j
real(rprec), intent(in) :: radius
integer, intent(out) :: foil_lo, foil_hi
real(rprec), intent(out) :: weight
integer :: sec, nsec, nfoil
real(rprec) :: r_lo, r_hi

nsec = turbineModel(j) % NumSec
nfoil = size(turbineModel(j) % airfoilType)
if (nsec <= 0 .or. nfoil <= 0) then
    foil_lo = 1
    foil_hi = 1
    weight = 0._rprec
    return
endif

if (radius <= turbineModel(j) % radius(1)) then
    foil_lo = max(1, min(nfoil, turbineModel(j) % sectionType(1)))
    foil_hi = foil_lo
    weight = 0._rprec
    return
elseif (radius >= turbineModel(j) % radius(nsec)) then
    foil_lo = max(1, min(nfoil, turbineModel(j) % sectionType(nsec)))
    foil_hi = foil_lo
    weight = 0._rprec
    return
endif

do sec = 2, nsec
    if (radius <= turbineModel(j) % radius(sec)) then
        r_lo = turbineModel(j) % radius(sec-1)
        r_hi = turbineModel(j) % radius(sec)
        if (r_hi > r_lo) then
            weight = (radius - r_lo) / (r_hi - r_lo)
        else
            weight = 0._rprec
        endif
        foil_lo = max(1, min(nfoil, turbineModel(j) % sectionType(sec-1)))
        foil_hi = max(1, min(nfoil, turbineModel(j) % sectionType(sec)))
        return
    endif
enddo

foil_lo = max(1, min(nfoil, turbineModel(j) % sectionType(nsec)))
foil_hi = foil_lo
weight = 0._rprec

end subroutine atm_airfoil_blend_info

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real(rprec) function atm_blended_airfoil_coeff(j, radius, alpha_deg,            &
                                               coeff_kind) result(value)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: j, coeff_kind
real(rprec), intent(in) :: radius, alpha_deg
integer :: foil_lo, foil_hi
real(rprec) :: weight, coeff_lo, coeff_hi

call atm_airfoil_blend_info(j, radius, foil_lo, foil_hi, weight)
coeff_lo = atm_airfoil_polar_value(turbineModel(j) % airfoilType(foil_lo),     &
                                   alpha_deg, coeff_kind)
coeff_hi = atm_airfoil_polar_value(turbineModel(j) % airfoilType(foil_hi),     &
                                   alpha_deg, coeff_kind)
value = coeff_lo + weight * (coeff_hi - coeff_lo)

end function atm_blended_airfoil_coeff

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeBladeForce(i,m,n,q,U_local)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the wind vectors by projecting the velocity
! onto the transformed coordinates system
implicit none

integer, intent(in) :: i,m,n,q
! i - turbineTypeArray
! n - numAnnulusSections
! q - numBladePoints
! m - numBl
real(rprec), intent(in) :: U_local(3)    ! The local velocity at this point

! Local variables
integer :: j ! Use to identify turbine type
integer, save :: structure_feedback_init = 0
integer, save :: structure_vel_feedback = 1
integer, save :: structure_alpha_feedback = 1
character(len=16) :: env_value
logical :: env_has_value
real(rprec) :: twistAng_i, chord_i, windAng_i, db_i, sigma!, base_alpha
!real(rprec) :: solidity_i
real(rprec), dimension(3) :: dragVector, liftVector

! Pointers to be used
real(rprec), pointer :: rotorApex(:)
real(rprec), pointer :: bladeAlignedVectors(:,:,:,:,:)
real(rprec), pointer :: windVectors(:,:,:,:)
real(rprec),  pointer :: bladePoints(:,:,:,:)
real(rprec), pointer :: rotSpeed
real(rprec),  pointer :: bladeRadius(:,:,:)
real(rprec), pointer :: PreCone
real(rprec), pointer :: solidity(:,:,:),cl(:,:,:),cd(:,:,:),cm(:,:,:),alpha(:,:,:)
real(rprec), pointer :: Vmag(:,:,:)
real(rprec), pointer :: du(:,:,:,:)
real(rprec) :: theta_tot_vel, vel_ax, vel_tg

! Identifier for the turbine type
j= turbineArray(i) % turbineTypeID

! Pointers to trubineArray (i)
rotorApex => turbineArray(i) % rotorApex
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
windVectors => turbineArray(i) % windVectors
bladePoints => turbineArray(i) % bladePoints
rotSpeed => turbineArray(i) % rotSpeed
solidity=> turbineArray(i) % solidity
bladeRadius => turbineArray(i) % bladeRadius
cd => turbineArray(i) % cd       ! Drag coefficient
cl => turbineArray(i) % cl       ! Lift coefficient
cm => turbineArray(i) % cm       ! Pitching moment coefficient
alpha => turbineArray(i) % alpha ! Angle of attack
Vmag => turbineArray(i) % Vmag ! Velocity magnitude
du => turbineArray(i) % du ! Change in velocity

if (structure_feedback_init == 0) then
    structure_feedback_init = 1
    call atm_model_env_token('LESGO_ATM_STRUCTURE_VEL_FEEDBACK', env_value,    &
        env_has_value)
    if (env_has_value) then
        if (env_value == '0') structure_vel_feedback = 0
    endif
    call atm_model_env_token('LESGO_ATM_STRUCTURE_ALPHA_FEEDBACK', env_value,  &
        env_has_value)
    if (env_has_value) then
        if (env_value == '0') structure_alpha_feedback = 0
    endif
    write(*,*) 'ATM structure feedback options: velocity=',                    &
        structure_vel_feedback, ' alpha=', structure_alpha_feedback
endif

! Pointers for turbineModel (j)
PreCone => turbineModel(j) % PreCone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This will compute the vectors defining the local coordinate
! system of the actuator point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define vector in z'
! If clockwise rotating, this vector points along the blade toward the tip.
! If counter-clockwise rotating, this vector points along the blade towards
! the root.
if (turbineArray(i) % rotationDir == "cw")  then
    bladeAlignedVectors(m,n,q,3,:) =      &
                                     vector_add(bladePoints(m,n,q,:),-rotorApex)
elseif (turbineArray(i) % rotationDir == "ccw") then
    bladeAlignedVectors(m,n,q,3,:) =      &
                                     vector_add(-bladePoints(m,n,q,:),rotorApex)
endif

bladeAlignedVectors(m,n,q,3,:) =  &
                        vector_divide(bladeAlignedVectors(m,n,q,3,:),   &
                        vector_mag(bladeAlignedVectors(m,n,q,3,:)) )

! Define vector in y'
bladeAlignedVectors(m,n,q,2,:) = cross_product(bladeAlignedVectors(m,n,q,3,:), &
                                 turbineArray(i) % uvShaft)

bladeAlignedVectors(m,n,q,2,:) = vector_divide(bladeAlignedVectors(m,n,q,2,:), &
                                 vector_mag(bladeAlignedVectors(m,n,q,2,:)))

! Define vector in x'
bladeAlignedVectors(m,n,q,1,:) = cross_product(bladeAlignedVectors(m,n,q,2,:), &
                                 bladeAlignedVectors(m,n,q,3,:))

bladeAlignedVectors(m,n,q,1,:) = vector_divide(bladeAlignedVectors(m,n,q,1,:), &
                                 vector_mag(bladeAlignedVectors(m,n,q,1,:)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This concludes the definition of the local coordinate system


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now put the velocity in that cell into blade-oriented coordinates and add on
! the velocity due to blade rotation.
windVectors(m,n,q,1) = dot_product(bladeAlignedVectors(m,n,q,1,:) , U_local)
windVectors(m,n,q,2) = dot_product(bladeAlignedVectors(m,n,q,2,:), U_local) + &
                      (rotSpeed * bladeRadius(m,n,q) * cos(PreCone))
windVectors(m,n,q,3) = dot_product(bladeAlignedVectors(m,n,q,3,:), U_local)

! Interpolated quantities through section
twistAng_i = turbineArray(i) % twistAng(m,n,q)
chord_i = turbineArray(i) % chord(m,n,q)

if (atm_structure_enabled() .and. structure_vel_feedback == 1) then
    theta_tot_vel = (turbineArray(i)%twistAng(m,n,q) +                         &
        turbineArray(i)%Pitch + turbineArray(i)%PitchControlAngle) * degRad +  &
        turbineArray(i)%elastic_twist(m,1,q)
    vel_ax = turbineArray(i)%flap_vel(m,1,q) * cos(theta_tot_vel) -            &
             turbineArray(i)%edge_vel(m,1,q) * sin(theta_tot_vel)
    vel_tg = turbineArray(i)%flap_vel(m,1,q) * sin(theta_tot_vel) +            &
             turbineArray(i)%edge_vel(m,1,q) * cos(theta_tot_vel)
    windVectors(m,n,q,1) = windVectors(m,n,q,1) - vel_ax
    windVectors(m,n,q,2) = windVectors(m,n,q,2) - vel_tg
endif

! Correct the velocity
if (turbineArray(i) % tipALMCorrection .eqv. .true.)  then

    windVectors(m,n,q,1) = windVectors(m,n,q,1) + du(m,n,q,1)
    windVectors(m,n,q,2) = windVectors(m,n,q,2) + du(m,n,q,2)

endif

! Velocity magnitude must reflect the same velocity used for AoA and forces.
Vmag(m,n,q)=sqrt( windVectors(m,n,q,1)**2+windVectors(m,n,q,2)**2 )

! Angle between wind vector components
windAng_i = atan2( windVectors(m,n,q,1), windVectors(m,n,q,2) ) /degRad

! Local angle of attack
alpha(m,n,q) = windAng_i - twistAng_i - turbineArray(i) % Pitch - &
               turbineArray(i) % PitchControlAngle
if (atm_structure_enabled() .and. structure_alpha_feedback == 1) then
    alpha(m,n,q) = alpha(m,n,q) - turbineArray(i) % elastic_twist(m,1,q) / degRad
endif

! Lift coefficient
cl(m,n,q) = atm_blended_airfoil_coeff(j, bladeRadius(m,n,q),                  &
                                      alpha(m,n,q), 1)

! Drag coefficient
cd(m,n,q) = atm_blended_airfoil_coeff(j, bladeRadius(m,n,q),                  &
                                      alpha(m,n,q), 2)
if (atm_structure_enabled()) then
    cm(m,n,q) = atm_blended_airfoil_coeff(j, bladeRadius(m,n,q),              &
                                          alpha(m,n,q), 3)
else
    cm(m,n,q) = 0._rprec
endif

! The blade section width
db_i = turbineArray(i) % db(q)

! Lift force
turbineArray(i) % lift(m,n,q) = 0.5_rprec * cl(m,n,q) * (Vmag(m,n,q)**2) *     &
                                chord_i * db_i * solidity(m,n,q)

! Drag force
turbineArray(i) % drag(m,n,q) = 0.5_rprec * cd(m,n,q) * (Vmag(m,n,q)**2) *     &
                                chord_i * db_i * solidity(m,n,q)
if (atm_structure_enabled()) then
    turbineArray(i) % pitchingMoment(m,n,q) = 0.5_rprec * cm(m,n,q) *          &
        (Vmag(m,n,q)**2) * chord_i * chord_i * db_i * solidity(m,n,q)
else
    turbineArray(i) % pitchingMoment(m,n,q) = 0._rprec
endif

! This vector projects the drag onto the local coordinate system
dragVector = bladeAlignedVectors(m,n,q,1,:)*windVectors(m,n,q,1) +             &
             bladeAlignedVectors(m,n,q,2,:)*windVectors(m,n,q,2)

dragVector = vector_divide(dragVector,vector_mag(dragVector) )

! Lift vector
liftVector = cross_product(dragVector,bladeAlignedVectors(m,n,q,3,:) )
liftVector = liftVector/vector_mag(liftVector)

! Apply the lift and drag as vectors
liftVector = -turbineArray(i) % lift(m,n,q) * liftVector;
dragVector = -turbineArray(i) % drag(m,n,q) * dragVector;

! The blade force is the total lift and drag vectors
turbineArray(i) % bladeForces(m,n,q,:) = vector_add(liftVector, dragVector)

! Find the component of the blade element force/density in the axial
! (along the shaft) direction.
turbineArray(i) % axialForce(m,n,q) = dot_product(                           &
        -turbineArray(i) % bladeForces(m,n,q,:), turbineArray(i) % uvShaft)

! Find the component of the blade element force/density in the tangential
! (torque-creating) direction.
turbineArray(i) % tangentialForce(m,n,q) = dot_product(                      &
       turbineArray(i) % bladeForces(m,n,q,:), bladeAlignedVectors(m,n,q,2,:))

! Change this back to radians
windAng_i = windAng_i * degRad
! The solidity
sigma = chord_i * turbineModel(j) % NumBl/ (2.*pi * bladeRadius(m,n,q) )

! Calculate the induction factor
turbineArray(i) % induction_a(m,n,q) = 1. / ( 4. * sin(windAng_i)**2 /  &
                (sigma * ( Cl(m,n,q) * cos(windAng_i) +   &
                Cd(m,n,q) * sin(windAng_i))) + 1.)

!~ write(*,*) 'Induction ', turbineArray(i) % induction_a(m,n,q)
! Calculate u infinity
turbineArray(i) % u_infinity(m,n,q) = windVectors(m,n,q,1) !/    &
!~                              (1. - turbineArray(i) % induction_a(m,n,q))
!~ turbineArray(i) % u_infinity(m,n,q) = Vmag(m,n,q) * sin(windAng_i) / &
!~ (1-turbineArray(i) % induction_a(m,n,q))


!~             turbineArray(i) % u_infinity = turbineArray(i) % u_infinity  +     &
!~                              windVectors(m,n,q,1) /    &
!~                              (1. - turbineArray(i) % induction_a(m,n,q))

! Calculate output quantities based on each point
call atm_process_output(i,m,n,q)

end subroutine atm_computeBladeForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_computeNacelleForce(i,U_local)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the force from the nacelle
implicit none

integer, intent(in) :: i ! i - turbineTypeArray
real(rprec), intent(in), dimension(3) :: U_local ! Velocity input

integer :: j ! j - turbineModel
real(rprec) :: V ! Velocity projected
real(rprec), dimension(3) :: nacelleAlignedVector ! Nacelle vector
real(rprec) :: area, drag

! Identifier for the turbine type
j= turbineArray(i) % turbineTypeID

area = pi * turbineModel(j) % hubRad **2

nacelleAlignedVector = turbineArray(i) % uvShaft

! Velocity projected in the direction of the nacelle
V = dot_product( nacelleAlignedVector , U_local)
!~     write(*,*) 'Nacelle Velocity before correction is: ', V
!~     write(*,*) 'nacelleAlignedVector is: ', nacelleAlignedVector

! The sampled velocity (uncorrected)
turbineArray(i) % VelNacelle_sampled = V

! Apply the velocity correction
V = V / (1. - .25/ pi * turbineArray(i) % nacelleCd * area /                   &
                turbineArray(i) % nacelleEpsilon**2 )
!~ write(*,*) 'Nacelle Velocity after correction is: ', V

! The velocity (corrected)
turbineArray(i) % VelNacelle_corrected = V

if (V .ge. 0.) then
    ! Drag force
    drag = 0.5_rprec * turbineArray(i) % nacelleCd * (V*V) * area

    ! Drag Vector
    turbineArray(i) % nacelleForce = - drag * nacelleAlignedVector
!~     write(*,*) 'Nacelle Cd= ', turbineArray(i) % nacelleCd
!~     write(*,*) 'Nacelle Force is: ', turbineArray(i) % nacelleForce
endif

end subroutine atm_computeNacelleForce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_integrate_u(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will compute the induction factor a
! for each actuator point
implicit none

integer, intent(in) :: i ! i - turbineTypeArray
#ifdef ENABLE_CUDA
real(rprec), managed, pointer, dimension(:,:,:) :: u_infinity, induction_a
real(rprec) :: u_sum, induction_sum, npoints
integer :: m, n, q, mmend, nnend, qqend
#else
real(rprec), pointer, dimension(:,:,:) :: u_infinity, induction_a
#endif
real(rprec), pointer :: u_infinity_mean

induction_a => turbineArray(i) % induction_a
u_infinity_mean => turbineArray(i) % u_infinity_mean
u_infinity => turbineArray(i) % u_infinity

#ifdef ENABLE_CUDA
if (atm_model_cuda_enabled()) then
    mmend = size(u_infinity, 1)
    nnend = size(u_infinity, 2)
    qqend = size(u_infinity, 3)
    npoints = real(mmend*nnend*qqend, rprec)
    u_sum = 0._rprec
    induction_sum = 0._rprec

    !$cuf kernel do(3) <<<*,*>>> reduction(+:u_sum, induction_sum)
    do q = 1, qqend
    do n = 1, nnend
    do m = 1, mmend
        u_sum = u_sum + u_infinity(m,n,q)
        induction_sum = induction_sum + induction_a(m,n,q)
    end do
    end do
    end do
    call atm_model_cuda_sync('atm_integrate_u')

    u_infinity_mean = (u_sum / npoints) /                                     &
        (1._rprec - induction_sum / npoints)
else
#endif
u_infinity_mean = sum(u_infinity) / size(u_infinity) / &
(1. - sum(induction_a) / size(induction_a))
#ifdef ENABLE_CUDA
endif
#endif

!~ write(*,*) "U infinity is", u_infinity_mean, size(u_infinity)

end subroutine atm_integrate_u


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_yawNacelle(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine yaws the nacelle according to the yaw angle
implicit none

integer, intent(in) :: i
integer :: j
integer :: m,n,q
real(rprec) :: origin_zero(3)
#ifdef ENABLE_CUDA
real(rprec), managed, pointer, dimension(:,:,:,:) :: bladePoints_d
real(rprec) :: origin1, origin2, origin3
real(rprec) :: ax1, ax2, ax3, cang, sang
real(rprec) :: rm11, rm12, rm13, rm21, rm22, rm23, rm31, rm32, rm33
real(rprec) :: p1, p2, p3, r1, r2, r3
#endif
! Perform rotation for the turbine.
j=turbineArray(i) % turbineTypeID

! Rotate the rotor apex first.
turbineArray(i) % rotorApex = rotatePoint(turbineArray(i) % rotorApex,         &
                              turbineArray(i) % towerShaftIntersect,           &
                              turbineArray(i) % uvTower,                       &
                              turbineArray(i) % deltaNacYaw)

! Rotate the shaft direction itself. Reconstructing it from
! rotorApex - towerShaftIntersect is singular when the overhang is zero.
origin_zero = (/ 0._rprec, 0._rprec, 0._rprec /)
turbineArray(i) % uvShaft = rotatePoint(turbineArray(i) % uvShaft,             &
                             origin_zero, turbineArray(i) % uvTower,           &
                             turbineArray(i) % deltaNacYaw)
turbineArray(i) % uvShaft = vector_divide(turbineArray(i) % uvShaft,           &
                            vector_mag(turbineArray(i) % uvShaft))

! Rotate turbine blades, blade by blade, point by point.
#ifdef ENABLE_CUDA
if (atm_model_cuda_enabled() .and. .not. atm_structure_enabled()) then
    bladePoints_d => turbineArray(i) % bladePoints
    origin1 = turbineArray(i) % towerShaftIntersect(1)
    origin2 = turbineArray(i) % towerShaftIntersect(2)
    origin3 = turbineArray(i) % towerShaftIntersect(3)
    ax1 = turbineArray(i) % uvTower(1)
    ax2 = turbineArray(i) % uvTower(2)
    ax3 = turbineArray(i) % uvTower(3)
    cang = cos(turbineArray(i) % deltaNacYaw)
    sang = sin(turbineArray(i) % deltaNacYaw)

    rm11 = ax1*ax1 + (1._rprec - ax1*ax1) * cang
    rm12 = ax1*ax2 * (1._rprec - cang) - ax3 * sang
    rm13 = ax1*ax3 * (1._rprec - cang) + ax2 * sang
    rm21 = ax1*ax2 * (1._rprec - cang) + ax3 * sang
    rm22 = ax2*ax2 + (1._rprec - ax2*ax2) * cang
    rm23 = ax2*ax3 * (1._rprec - cang) - ax1 * sang
    rm31 = ax1*ax3 * (1._rprec - cang) - ax2 * sang
    rm32 = ax2*ax3 * (1._rprec - cang) + ax1 * sang
    rm33 = ax3*ax3 + (1._rprec - ax3*ax3) * cang

    !$cuf kernel do(3) <<<*,*>>>
    do q=1, turbineArray(i) % numBladePoints
        do n=1, turbineArray(i) %  numAnnulusSections
            do m=1, turbineModel(j) % numBl
                p1 = bladePoints_d(m,n,q,1) - origin1
                p2 = bladePoints_d(m,n,q,2) - origin2
                p3 = bladePoints_d(m,n,q,3) - origin3

                r1 = rm11*p1 + rm12*p2 + rm13*p3
                r2 = rm21*p1 + rm22*p2 + rm23*p3
                r3 = rm31*p1 + rm32*p2 + rm33*p3

                bladePoints_d(m,n,q,1) = r1 + origin1
                bladePoints_d(m,n,q,2) = r2 + origin2
                bladePoints_d(m,n,q,3) = r3 + origin3
            enddo
        enddo
    enddo

    call atm_model_cuda_sync('atm_yawNacelle')
else
#endif
do q=1, turbineArray(i) % numBladePoints
    do n=1, turbineArray(i) %  numAnnulusSections
        do m=1, turbineModel(j) % numBl
            if (atm_structure_enabled()) then
                turbineArray(i) % bladePoints_rigid(m,n,q,:) =                &
                rotatePoint(turbineArray(i) % bladePoints_rigid(m,n,q,:),     &
                turbineArray(i) % towerShaftIntersect,                        &
                turbineArray(i) % uvTower,                                    &
                turbineArray(i) % deltaNacYaw )
            endif
            turbineArray(i) % bladePoints(m,n,q,:) =                           &
            rotatePoint( turbineArray(i) % bladePoints(m,n,q,:),               &
            turbineArray(i) % towerShaftIntersect,                             &
            turbineArray(i) % uvTower,                                         &
            turbineArray(i) % deltaNacYaw )
        enddo
    enddo
enddo
#ifdef ENABLE_CUDA
endif
#endif

! Compute the new yaw angle and make sure it isn't bigger than 2*pi.
if (pastFirstTimeStep) then
    turbineArray(i) % nacYaw = turbineArray(i) % nacYaw +                      &
                               turbineArray(i) % deltaNacYaw
    if (turbineArray(i) % nacYaw .ge. 2.0 * pi) then
        turbineArray(i) % nacYaw = turbineArray(i) % nacYaw - 2.0 * pi
    endif
endif

end subroutine atm_yawNacelle

!~ !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!~ subroutine atm_compassToStandard(dir)
!~ ! This function converts nacelle yaw from compass directions to the standard
!~ ! convention of 0 degrees on the + x axis with positive degrees
!~ ! in the counter-clockwise direction.
!~ !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!~ implicit none
!~ real(rprec), intent(inout) :: dir
!~ dir = dir + 180.0
!~ if (dir .ge. 360.0) then
!~     dir = dir - 360.0
!~ endif
!~ dir = 90.0 - dir
!~ if (dir < 0.0) then
!~     dir = dir + 360.0
!~ endif
!~
!~ end subroutine atm_compassToStandard

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_output(i, jt_total, time)
! This subroutine will calculate and write the output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: jt_total ! Number of iteration fed in from solver
real(rprec), intent(in) :: time  ! time from simulation
integer, intent(in) :: i  ! The turbine number
integer :: j, m
integer :: powerFile=11, rotSpeedFile=12, bladeFile=13, liftFile=14, dragFile=15
integer :: ClFile=16, CdFile=17, alphaFile=18, VrelFile=19
integer :: VaxialFile=20, VtangentialFile=21, pitchFile=22, thrustFile=23
integer :: tangentialForceFile=24, axialForceFile=25, yawfile=26, nacelleFile=27
integer :: ClbFile=28
integer :: uy_LESFile=29
integer :: uy_optFile=30
integer :: elasticTwistFile=31, flapDispFile=32, edgeDispFile=33, CmFile=34
logical :: writeDiagnostics

! Output only if the number of intervals is right
if ( mod(jt_total-1, outputInterval) == 0) then

    writeDiagnostics = .false.
    if (diagnosticOutputInterval > 0) then
        writeDiagnostics = mod(jt_total-1, diagnosticOutputInterval) == 0
    endif

    ! Turbine output files carry the data.  Avoid per-output stdout banners;
    ! they become a measurable serial I/O cost in many-turbine output-heavy runs.

    j=turbineArray(i) % turbineTypeID ! The turbine type ID

    ! File for power output
    open(unit=powerFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/power")

    ! File for thrust
    open(unit=thrustFile,position="append",                                    &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/thrust")

    ! File for rotor speed
    open(unit=RotSpeedFile,position="append",                                  &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/RotSpeed")

    ! File for yaw
    open(unit=YawFile,position="append",                                  &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Yaw")

    open(unit=pitchFile,position="append",                                     &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/pitch")

    open(unit=nacelleFile,position="append",                                  &
    file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/nacelle")

    if (writeDiagnostics) then
        ! Heavy blade/airfoil diagnostics are useful for debugging, but are
        ! not required for scalar power/thrust convergence histories.
        open(unit=bladeFile,position="append",                                 &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/blade")

        open(unit=liftFile,position="append",                                  &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/lift")

        open(unit=dragFile,position="append",                                  &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/drag")

        open(unit=ClFile,position="append",                                    &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cl")

        open(unit=ClbFile,position="append",                                   &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Clb")

        open(unit=CdFile,position="append",                                    &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cd")

        open(unit=alphaFile,position="append",                                 &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/alpha")

        open(unit=VrelFile,position="append",                                  &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vrel")

        open(unit=VaxialFile,position="append",                                &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vaxial")

        open(unit=VtangentialFile,position="append",                           &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//         &
             "/Vtangential")

        open(unit=tangentialForceFile,position="append",                       &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//         &
             "/tangentialForce")

        open(unit=axialForceFile,position="append",                            &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//         &
             "/axialForce")

        open(unit=uy_LESFile,position="append",                                &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/uy_LES")

        open(unit=uy_optFile,position="append",                                &
        file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/uy_opt")

        if (atm_structure_enabled()) then
            open(unit=elasticTwistFile,position="append",                     &
            file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//     &
                 "/elastic_twist")
            open(unit=flapDispFile,position="append",                         &
            file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//     &
                 "/flap_disp")
            open(unit=edgeDispFile,position="append",                         &
            file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//     &
                 "/edge_disp")
            open(unit=CmFile,position="append",                               &
            file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cm")
        endif
    endif

    call atm_compute_power(i)
    write(powerFile,*) time, turbineArray(i) % powerRotor,                     &
                       turbineArray(i) % powerGen
    write(thrustFile,*) time, turbineArray(i) % thrust *                      &
                         turbineArray(i) % fluidDensity
    write(RotSpeedFile,*) time, turbineArray(i) % RotSpeed
    write(pitchFile,*) time, turbineArray(i) % PitchControlAngle,              &
                       turbineArray(i) % IntSpeedError
    write(YawFile,*) time, turbineArray(i) % deltaNacYaw,                      &
                           turbineArray(i) % NacYaw
    write(nacelleFile,*) time, turbineArray(i) % VelNacelle_sampled,           &
                           turbineArray(i) % VelNacelle_corrected

    if (writeDiagnostics) then
    ! Will write only the first actuator section of the blade
    do m=1, turbineModel(j) % numBl
        write(bladeFile,*) i, m, turbineArray(i) % bladeRadius(m,1,:)
        ! Internal blade loads are stored as section-integrated force per
        ! density.  Diagnostics are written as physical sectional loads.
        write(liftFile,*) i, m, turbineArray(i) % lift(m,1,:) *                &
                                turbineArray(i) % fluidDensity /               &
                                turbineArray(i) % db(:)
        write(dragFile,*) i, m, turbineArray(i) % drag(m,1,:) *                &
                                turbineArray(i) % fluidDensity /               &
                                turbineArray(i) % db(:)
        write(ClFile,*) i, m, turbineArray(i) % cl(m,1,:)
        write(ClbFile,*) i, m, turbineArray(i) % cl_b(m,1,:)
        write(CdFile,*) i, m, turbineArray(i) % cd(m,1,:)
        write(alphaFile,*) i, m, turbineArray(i) % alpha(m,1,:)
        write(VrelFile,*) i, m, turbineArray(i) % Vmag(m,1,:)
        write(VaxialFile,*) i, m, turbineArray(i) % windVectors(m,1,:,1)
        write(VtangentialFile,*) i, m, turbineArray(i) %                       &
                                       windVectors(m,1,:,2)
        write(tangentialForceFile,*) i, m, turbineArray(i) %                   &
                                            tangentialForce(m,1,:) *            &
                                            turbineArray(i) % fluidDensity /   &
                                            turbineArray(i) % db(:)
        write(axialForceFile,*) i, m, turbineArray(i) % axialForce(m,1,:) *    &
                                      turbineArray(i) % fluidDensity /         &
                                      turbineArray(i) % db(:)
        write(uy_LESFile,*) i, m, turbineArray(i) % uy_LES(m,1,:)
        write(uy_optFile,*) i, m, turbineArray(i) % uy_opt(m,1,:)
        if (atm_structure_enabled()) then
            write(elasticTwistFile,*) i, m,                                   &
                turbineArray(i) % elastic_twist(m,1,:)
            write(flapDispFile,*) i, m, turbineArray(i) % flap_disp(m,1,:)
            write(edgeDispFile,*) i, m, turbineArray(i) % edge_disp(m,1,:)
            write(CmFile,*) i, m, turbineArray(i) % cm(m,1,:)
        endif

    enddo
    endif

        ! Write blade points
!~         call atm_write_blade_points(i,jt_total)

    ! Close scalar files every output step.
    close(powerFile)
    close(thrustFile)
    close(rotSpeedFile)
    close(pitchFile)
    close(yawFile)
    close(nacelleFile)

    if (writeDiagnostics) then
        close(bladeFile)
        close(liftFile)
        close(dragFile)
        close(ClFile)
        close(ClbFile)
        close(CdFile)
        close(alphaFile)
        close(VrelFile)
        close(VaxialFile)
        close(VtangentialFile)
        close(tangentialForceFile)
        close(axialForceFile)
        close(uy_LESFile)
        close(uy_optFile)
        if (atm_structure_enabled()) then
            close(elasticTwistFile)
            close(flapDispFile)
            close(edgeDispFile)
            close(CmFile)
        endif
    endif

endif

end subroutine atm_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_compute_power(i)
! This subroutine will calculate the total power of the turbine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i
integer :: j

j = turbineArray(i) % turbineTypeID

turbineArray(i) % powerRotor = turbineArray(i) % torqueRotor *                 &
    turbineArray(i) % rotSpeed * turbineArray(i) % fluidDensity
if (turbineModel(j) % TorqueControllerType == "fiveRegion") then
    turbineArray(i) % powerGen = turbineArray(i) % torqueGen *                 &
        turbineArray(i) % rotSpeed * turbineModel(j) % GBRatio
else
    turbineArray(i) % powerGen = turbineArray(i) % powerRotor
endif

if (atm_power_stdout_enabled()) then
    write(*,*) 'Turbine ',i,' (Aerodynamic, Generator) Power is: ',            &
        turbineArray(i) % powerRotor, turbineArray(i) % powerGen
endif

end subroutine atm_compute_power

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_write_blade_points(i,time_counter)
! This subroutine writes the position of all the blades at each time step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i, time_counter
integer :: m, n, q, j

j=turbineArray(i) % turbineTypeID ! The turbine type ID

open(unit=231, file="./turbineOutput/"//trim(turbineArray(i) % turbineName)//  &
               '/blades'//trim(int2str(time_counter))//".vtk")

! Write the points to the blade file
do m=1, turbineModel(j) % numBl

    do n=1, turbineArray(i) %  numAnnulusSections

        do q=1, turbineArray(i) % numBladePoints

            write(231,*) turbineArray(i) % bladePoints(m,n,q,:)

        enddo

    enddo

enddo

close(231)

end subroutine atm_write_blade_points

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_process_output(i,m,n,q)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This subroutine will process the output for an individual actuator point
! Important quantities such as power, thrust, ect are calculated here
implicit none

integer, intent(in) :: i,m,n,q
integer :: j
! Identify the turbine type
j=turbineArray(i) % turbineTypeID

turbineArray(i) % axialForce(m,n,q) = dot_product(                             &
-turbineArray(i) % bladeForces(m,n,q,:), turbineArray(i) % uvShaft)

! Find the component of the blade element force/density in the
! tangential (torque-creating) direction.
turbineArray(i) % tangentialForce(m,n,q) = dot_product(                        &
                  turbineArray(i) % bladeForces(m,n,q,:) ,                     &
                  turbineArray(i) % bladeAlignedVectors(m,n,q,2,:))

! Add this blade element's contribution to thrust to the total
! turbine thrust.
turbineArray(i) % thrust = turbineArray(i) % thrust +                          &
                           turbineArray(i) % axialForce(m,n,q)

! Add this blade element's contribution to aerodynamic torque to
! the total turbine aerodynamic torque.
turbineArray(i) % torqueRotor = turbineArray(i) % torqueRotor +                &
                                turbineArray(i) % tangentialForce(m,n,q) *     &
                                turbineArray(i) % bladeRadius(m,n,q) *         &
                                cos(turbineModel(j) % PreCone)


end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical function atm_structure_has_data(i)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: i

atm_structure_has_data = allocated(turbineArray(i) % EI_blade) .and.           &
    allocated(turbineArray(i) % rho_blade) .and.                               &
    allocated(turbineArray(i) % EI_edge_blade) .and.                           &
    allocated(turbineArray(i) % GJ_blade)

if (atm_structure_has_data) then
    atm_structure_has_data = maxval(abs(turbineArray(i) % EI_blade)) >          &
        tiny(1._rprec) .and. maxval(abs(turbineArray(i) % rho_blade)) >         &
        tiny(1._rprec) .and. maxval(abs(turbineArray(i) % EI_edge_blade)) >     &
        tiny(1._rprec) .and. maxval(abs(turbineArray(i) % GJ_blade)) >          &
        tiny(1._rprec)
endif

end function atm_structure_has_data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine atm_solve_structure(i, dt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Two-way blade line-model update. This first integration keeps the solver on
! the host for traceability; the flag is off by default and GPU batching follows
! after this path is validated.
implicit none

integer, intent(in) :: i
real(rprec), intent(in) :: dt

integer :: j, m, q, s, N, DOFs
real(rprec), parameter :: beta = 0.25_rprec
real(rprec), parameter :: gamma = 0.5_rprec
real(rprec), parameter :: rayleigh_m = 0.03_rprec
real(rprec), parameter :: rayleigh_k = 0.0005_rprec
real(rprec) :: a0, a1, a2, a3, a4, a5
real(rprec) :: L, EI_f, EI_e, GJ_val, rho_val, c_kg, node_mass
real(rprec) :: theta_tot, F_aero_ax, F_aero_tg, F_grav_ax, F_grav_tg
real(rprec) :: old_acc, new_acc
logical :: structure_timing
real(rprec) :: structure_t0, structure_t1, structure_stage_t0
real(rprec) :: structure_assembly_local, structure_solve_local
real(rprec) :: structure_update_local, structure_total_local
real(rprec), allocatable :: Kf(:,:), Ke(:,:), Mmat(:,:), Cf(:,:), Ce(:,:)
real(rprec), allocatable :: Af(:,:), Ae(:,:), Kt(:,:)
real(rprec), allocatable :: Ff(:), Fe(:), Ft(:), Uf(:), Ue(:), Ut(:)
real(rprec), allocatable :: df(:), vf(:), accf(:), de(:), ve(:), acce(:)
real(rprec), allocatable :: Ttension(:)
logical, save :: missing_data_warned = .false.
#ifdef ENABLE_CUDA
logical :: use_gpu_solve
logical, save :: structure_gpu_path_printed = .false.
#endif

if (.not. atm_structure_enabled()) return
if (.not. atm_structure_has_data(i)) then
    if (.not. missing_data_warned) then
        write(*,*) 'ATM structure enabled, but turbine structural properties ', &
            'were not found in BladeData; skipping structural update.'
        missing_data_warned = .true.
    endif
    return
endif

j = turbineArray(i) % turbineTypeID
N = turbineArray(i) % numBladePoints
DOFs = 2 * N

structure_timing = atm_structure_timing_enabled()
structure_assembly_local = 0._rprec
structure_solve_local = 0._rprec
structure_update_local = 0._rprec
structure_total_local = 0._rprec
if (structure_timing) call cpu_time(structure_t0)

#ifdef ENABLE_CUDA
use_gpu_solve = atm_structure_gpu_enabled()
if (.not. structure_gpu_path_printed) then
    write(*,*) 'ATM structural FP64 GPU solve enabled=', use_gpu_solve
    structure_gpu_path_printed = .true.
endif
#endif

a0 = 1._rprec / (beta * dt * dt)
a1 = gamma / (beta * dt)
a2 = 1._rprec / (beta * dt)
a3 = 1._rprec / (2._rprec * beta) - 1._rprec
a4 = gamma / beta - 1._rprec
a5 = dt * (gamma / (2._rprec * beta) - 1._rprec)

allocate(Kf(DOFs,DOFs), Ke(DOFs,DOFs), Mmat(DOFs,DOFs))
allocate(Cf(DOFs,DOFs), Ce(DOFs,DOFs), Af(DOFs,DOFs), Ae(DOFs,DOFs))
allocate(Kt(N,N), Ff(DOFs), Fe(DOFs), Ft(N), Uf(DOFs), Ue(DOFs), Ut(N))
allocate(df(DOFs), vf(DOFs), accf(DOFs), de(DOFs), ve(DOFs), acce(DOFs))
allocate(Ttension(N))

do m = 1, turbineModel(j) % NumBl
    Kf = 0._rprec
    Ke = 0._rprec
    Mmat = 0._rprec
    Kt = 0._rprec
    Ff = 0._rprec
    Fe = 0._rprec
    Ft = 0._rprec
    Ttension = 0._rprec
    if (structure_timing) call cpu_time(structure_stage_t0)

    do q = N, 1, -1
        node_mass = turbineArray(i)%rho_blade(m,1,q) * turbineArray(i)%db(q)
        if (q == N) then
            Ttension(q) = node_mass * turbineArray(i)%rotSpeed**2 *             &
                turbineArray(i)%bladeRadius(m,1,q)
        else
            Ttension(q) = Ttension(q+1) + node_mass *                           &
                turbineArray(i)%rotSpeed**2 * turbineArray(i)%bladeRadius(m,1,q)
        endif
    enddo

    do q = 1, N - 1
        L = turbineArray(i)%bladeRadius(m,1,q+1) -                              &
            turbineArray(i)%bladeRadius(m,1,q)
        EI_f = 0.5_rprec * (turbineArray(i)%EI_blade(m,1,q) +                   &
                            turbineArray(i)%EI_blade(m,1,q+1))
        EI_e = 0.5_rprec * (turbineArray(i)%EI_edge_blade(m,1,q) +              &
                            turbineArray(i)%EI_edge_blade(m,1,q+1))
        GJ_val = 0.5_rprec * (turbineArray(i)%GJ_blade(m,1,q) +                 &
                              turbineArray(i)%GJ_blade(m,1,q+1))
        rho_val = 0.5_rprec * (turbineArray(i)%rho_blade(m,1,q) +               &
                               turbineArray(i)%rho_blade(m,1,q+1))
        c_kg = 0.5_rprec * (Ttension(q) + Ttension(q+1)) / (30._rprec * L)

        Kf(2*q-1,2*q-1) = Kf(2*q-1,2*q-1) + 12._rprec*EI_f/L**3 + c_kg*36._rprec
        Kf(2*q-1,2*q  ) = Kf(2*q-1,2*q  ) +  6._rprec*EI_f/L**2 + c_kg*3._rprec*L
        Kf(2*q-1,2*q+1) = Kf(2*q-1,2*q+1) - 12._rprec*EI_f/L**3 - c_kg*36._rprec
        Kf(2*q-1,2*q+2) = Kf(2*q-1,2*q+2) +  6._rprec*EI_f/L**2 + c_kg*3._rprec*L
        Kf(2*q,  2*q-1) = Kf(2*q,  2*q-1) +  6._rprec*EI_f/L**2 + c_kg*3._rprec*L
        Kf(2*q,  2*q  ) = Kf(2*q,  2*q  ) +  4._rprec*EI_f/L    + c_kg*4._rprec*L**2
        Kf(2*q,  2*q+1) = Kf(2*q,  2*q+1) -  6._rprec*EI_f/L**2 - c_kg*3._rprec*L
        Kf(2*q,  2*q+2) = Kf(2*q,  2*q+2) +  2._rprec*EI_f/L    - c_kg*L**2
        Kf(2*q+1,2*q-1) = Kf(2*q+1,2*q-1) - 12._rprec*EI_f/L**3 - c_kg*36._rprec
        Kf(2*q+1,2*q  ) = Kf(2*q+1,2*q  ) -  6._rprec*EI_f/L**2 - c_kg*3._rprec*L
        Kf(2*q+1,2*q+1) = Kf(2*q+1,2*q+1) + 12._rprec*EI_f/L**3 + c_kg*36._rprec
        Kf(2*q+1,2*q+2) = Kf(2*q+1,2*q+2) -  6._rprec*EI_f/L**2 - c_kg*3._rprec*L
        Kf(2*q+2,2*q-1) = Kf(2*q+2,2*q-1) +  6._rprec*EI_f/L**2 + c_kg*3._rprec*L
        Kf(2*q+2,2*q  ) = Kf(2*q+2,2*q  ) +  2._rprec*EI_f/L    - c_kg*L**2
        Kf(2*q+2,2*q+1) = Kf(2*q+2,2*q+1) -  6._rprec*EI_f/L**2 - c_kg*3._rprec*L
        Kf(2*q+2,2*q+2) = Kf(2*q+2,2*q+2) +  4._rprec*EI_f/L    + c_kg*4._rprec*L**2

        Ke(2*q-1,2*q-1) = Ke(2*q-1,2*q-1) + 12._rprec*EI_e/L**3 + c_kg*36._rprec
        Ke(2*q-1,2*q  ) = Ke(2*q-1,2*q  ) +  6._rprec*EI_e/L**2 + c_kg*3._rprec*L
        Ke(2*q-1,2*q+1) = Ke(2*q-1,2*q+1) - 12._rprec*EI_e/L**3 - c_kg*36._rprec
        Ke(2*q-1,2*q+2) = Ke(2*q-1,2*q+2) +  6._rprec*EI_e/L**2 + c_kg*3._rprec*L
        Ke(2*q,  2*q-1) = Ke(2*q,  2*q-1) +  6._rprec*EI_e/L**2 + c_kg*3._rprec*L
        Ke(2*q,  2*q  ) = Ke(2*q,  2*q  ) +  4._rprec*EI_e/L    + c_kg*4._rprec*L**2
        Ke(2*q,  2*q+1) = Ke(2*q,  2*q+1) -  6._rprec*EI_e/L**2 - c_kg*3._rprec*L
        Ke(2*q,  2*q+2) = Ke(2*q,  2*q+2) +  2._rprec*EI_e/L    - c_kg*L**2
        Ke(2*q+1,2*q-1) = Ke(2*q+1,2*q-1) - 12._rprec*EI_e/L**3 - c_kg*36._rprec
        Ke(2*q+1,2*q  ) = Ke(2*q+1,2*q  ) -  6._rprec*EI_e/L**2 - c_kg*3._rprec*L
        Ke(2*q+1,2*q+1) = Ke(2*q+1,2*q+1) + 12._rprec*EI_e/L**3 + c_kg*36._rprec
        Ke(2*q+1,2*q+2) = Ke(2*q+1,2*q+2) -  6._rprec*EI_e/L**2 - c_kg*3._rprec*L
        Ke(2*q+2,2*q-1) = Ke(2*q+2,2*q-1) +  6._rprec*EI_e/L**2 + c_kg*3._rprec*L
        Ke(2*q+2,2*q  ) = Ke(2*q+2,2*q  ) +  2._rprec*EI_e/L    - c_kg*L**2
        Ke(2*q+2,2*q+1) = Ke(2*q+2,2*q+1) -  6._rprec*EI_e/L**2 - c_kg*3._rprec*L
        Ke(2*q+2,2*q+2) = Ke(2*q+2,2*q+2) +  4._rprec*EI_e/L    + c_kg*4._rprec*L**2

        Mmat(2*q-1,2*q-1) = Mmat(2*q-1,2*q-1) + 156._rprec*rho_val*L/420._rprec
        Mmat(2*q-1,2*q  ) = Mmat(2*q-1,2*q  ) +  22._rprec*rho_val*L**2/420._rprec
        Mmat(2*q-1,2*q+1) = Mmat(2*q-1,2*q+1) +  54._rprec*rho_val*L/420._rprec
        Mmat(2*q-1,2*q+2) = Mmat(2*q-1,2*q+2) -  13._rprec*rho_val*L**2/420._rprec
        Mmat(2*q,  2*q-1) = Mmat(2*q,  2*q-1) +  22._rprec*rho_val*L**2/420._rprec
        Mmat(2*q,  2*q  ) = Mmat(2*q,  2*q  ) +   4._rprec*rho_val*L**3/420._rprec
        Mmat(2*q,  2*q+1) = Mmat(2*q,  2*q+1) +  13._rprec*rho_val*L**2/420._rprec
        Mmat(2*q,  2*q+2) = Mmat(2*q,  2*q+2) -   3._rprec*rho_val*L**3/420._rprec
        Mmat(2*q+1,2*q-1) = Mmat(2*q+1,2*q-1) +  54._rprec*rho_val*L/420._rprec
        Mmat(2*q+1,2*q  ) = Mmat(2*q+1,2*q  ) +  13._rprec*rho_val*L**2/420._rprec
        Mmat(2*q+1,2*q+1) = Mmat(2*q+1,2*q+1) + 156._rprec*rho_val*L/420._rprec
        Mmat(2*q+1,2*q+2) = Mmat(2*q+1,2*q+2) -  22._rprec*rho_val*L**2/420._rprec
        Mmat(2*q+2,2*q-1) = Mmat(2*q+2,2*q-1) -  13._rprec*rho_val*L**2/420._rprec
        Mmat(2*q+2,2*q  ) = Mmat(2*q+2,2*q  ) -   3._rprec*rho_val*L**3/420._rprec
        Mmat(2*q+2,2*q+1) = Mmat(2*q+2,2*q+1) -  22._rprec*rho_val*L**2/420._rprec
        Mmat(2*q+2,2*q+2) = Mmat(2*q+2,2*q+2) +   4._rprec*rho_val*L**3/420._rprec

        Kt(q,q) = Kt(q,q) + GJ_val / L
        Kt(q,q+1) = Kt(q,q+1) - GJ_val / L
        Kt(q+1,q) = Kt(q+1,q) - GJ_val / L
        Kt(q+1,q+1) = Kt(q+1,q+1) + GJ_val / L
    enddo

    Cf = rayleigh_m * Mmat + rayleigh_k * Kf
    Ce = rayleigh_m * Mmat + rayleigh_k * Ke

    do q = 1, N
        node_mass = turbineArray(i)%rho_blade(m,1,q) * turbineArray(i)%db(q)
        theta_tot = (turbineArray(i)%twistAng(m,1,q) + turbineArray(i)%Pitch + &
            turbineArray(i)%PitchControlAngle) * degRad +                      &
            turbineArray(i)%elastic_twist(m,1,q)
        F_aero_ax = -dot_product(turbineArray(i)%bladeForces(m,1,q,:),         &
            turbineArray(i)%bladeAlignedVectors(m,1,q,1,:)) *                  &
            turbineArray(i)%fluidDensity
        F_aero_tg = -dot_product(turbineArray(i)%bladeForces(m,1,q,:),         &
            turbineArray(i)%bladeAlignedVectors(m,1,q,2,:)) *                  &
            turbineArray(i)%fluidDensity
        F_grav_ax = node_mass * (-9.81_rprec *                                 &
            turbineArray(i)%bladeAlignedVectors(m,1,q,1,3))
        F_grav_tg = node_mass * (-9.81_rprec *                                 &
            turbineArray(i)%bladeAlignedVectors(m,1,q,2,3))
        Ff(2*q-1) = (F_aero_ax + F_grav_ax) * cos(theta_tot) +                 &
                    (F_aero_tg + F_grav_tg) * sin(theta_tot)
        Fe(2*q-1) = -(F_aero_ax + F_grav_ax) * sin(theta_tot) +                &
                    (F_aero_tg + F_grav_tg) * cos(theta_tot)
        ! Airfoil Cm uses the aerodynamic nose-down convention, opposite to the
        ! structural elastic-twist DOF used as positive pitch/twist in AoA.
        Ft(q) = -turbineArray(i)%pitchingMoment(m,1,q) *                       &
            turbineArray(i)%fluidDensity

        df(2*q-1) = turbineArray(i)%flap_disp(m,1,q)
        df(2*q) = turbineArray(i)%theta_disp(m,1,q)
        vf(2*q-1) = turbineArray(i)%flap_vel(m,1,q)
        vf(2*q) = turbineArray(i)%theta_vel(m,1,q)
        accf(2*q-1) = turbineArray(i)%flap_acc(m,1,q)
        accf(2*q) = turbineArray(i)%theta_acc(m,1,q)
        de(2*q-1) = turbineArray(i)%edge_disp(m,1,q)
        de(2*q) = turbineArray(i)%edge_theta_disp(m,1,q)
        ve(2*q-1) = turbineArray(i)%edge_vel(m,1,q)
        ve(2*q) = turbineArray(i)%edge_theta_vel(m,1,q)
        acce(2*q-1) = turbineArray(i)%edge_acc(m,1,q)
        acce(2*q) = turbineArray(i)%edge_theta_acc(m,1,q)
    enddo

    Af = Kf + a0 * Mmat + a1 * Cf
    Ae = Ke + a0 * Mmat + a1 * Ce
    Ff = Ff + matmul(Mmat, a0*df + a2*vf + a3*accf) +                          &
         matmul(Cf, a1*df + a4*vf + a5*accf)
    Fe = Fe + matmul(Mmat, a0*de + a2*ve + a3*acce) +                          &
         matmul(Ce, a1*de + a4*ve + a5*acce)

    Af(1,:) = 0._rprec; Af(1,1) = 1._rprec; Ff(1) = 0._rprec
    Af(2,:) = 0._rprec; Af(2,2) = 1._rprec; Ff(2) = 0._rprec
    Ae(1,:) = 0._rprec; Ae(1,1) = 1._rprec; Fe(1) = 0._rprec
    Ae(2,:) = 0._rprec; Ae(2,2) = 1._rprec; Fe(2) = 0._rprec
    Kt(1,:) = 0._rprec; Kt(1,1) = 1._rprec; Ft(1) = 0._rprec

    if (structure_timing) then
        call cpu_time(structure_t1)
        structure_assembly_local = structure_assembly_local +                 &
            max(structure_t1 - structure_stage_t0, 0._rprec)
        structure_stage_t0 = structure_t1
    endif

#ifdef ENABLE_CUDA
    if (use_gpu_solve) then
        call solve_linear_system_gpu_dp(DOFs, Af, Ff, Uf)
        call solve_linear_system_gpu_dp(DOFs, Ae, Fe, Ue)
        call solve_linear_system_gpu_dp(N, Kt, Ft, Ut)
    else
        call solve_linear_system_banded_dp(DOFs, 3, Af, Ff, Uf)
        call solve_linear_system_banded_dp(DOFs, 3, Ae, Fe, Ue)
        call solve_linear_system_banded_dp(N, 1, Kt, Ft, Ut)
    endif
#else
    call solve_linear_system_banded_dp(DOFs, 3, Af, Ff, Uf)
    call solve_linear_system_banded_dp(DOFs, 3, Ae, Fe, Ue)
    call solve_linear_system_banded_dp(N, 1, Kt, Ft, Ut)
#endif

    if (structure_timing) then
        call cpu_time(structure_t1)
        structure_solve_local = structure_solve_local +                       &
            max(structure_t1 - structure_stage_t0, 0._rprec)
        structure_stage_t0 = structure_t1
    endif

    do q = 1, N
        old_acc = turbineArray(i)%flap_acc(m,1,q)
        new_acc = a0 * (Uf(2*q-1) - turbineArray(i)%flap_disp(m,1,q)) -        &
                  a2 * turbineArray(i)%flap_vel(m,1,q) - a3 * old_acc
        turbineArray(i)%flap_vel(m,1,q) = turbineArray(i)%flap_vel(m,1,q) +    &
            dt * ((1._rprec - gamma) * old_acc + gamma * new_acc)
        turbineArray(i)%flap_acc(m,1,q) = new_acc
        turbineArray(i)%flap_disp(m,1,q) = Uf(2*q-1)

        old_acc = turbineArray(i)%theta_acc(m,1,q)
        new_acc = a0 * (Uf(2*q) - turbineArray(i)%theta_disp(m,1,q)) -         &
                  a2 * turbineArray(i)%theta_vel(m,1,q) - a3 * old_acc
        turbineArray(i)%theta_vel(m,1,q) = turbineArray(i)%theta_vel(m,1,q) +  &
            dt * ((1._rprec - gamma) * old_acc + gamma * new_acc)
        turbineArray(i)%theta_acc(m,1,q) = new_acc
        turbineArray(i)%theta_disp(m,1,q) = Uf(2*q)

        old_acc = turbineArray(i)%edge_acc(m,1,q)
        new_acc = a0 * (Ue(2*q-1) - turbineArray(i)%edge_disp(m,1,q)) -        &
                  a2 * turbineArray(i)%edge_vel(m,1,q) - a3 * old_acc
        turbineArray(i)%edge_vel(m,1,q) = turbineArray(i)%edge_vel(m,1,q) +    &
            dt * ((1._rprec - gamma) * old_acc + gamma * new_acc)
        turbineArray(i)%edge_acc(m,1,q) = new_acc
        turbineArray(i)%edge_disp(m,1,q) = Ue(2*q-1)

        old_acc = turbineArray(i)%edge_theta_acc(m,1,q)
        new_acc = a0 * (Ue(2*q) - turbineArray(i)%edge_theta_disp(m,1,q)) -    &
                  a2 * turbineArray(i)%edge_theta_vel(m,1,q) - a3 * old_acc
        turbineArray(i)%edge_theta_vel(m,1,q) =                                &
            turbineArray(i)%edge_theta_vel(m,1,q) +                            &
            dt * ((1._rprec - gamma) * old_acc + gamma * new_acc)
        turbineArray(i)%edge_theta_acc(m,1,q) = new_acc
        turbineArray(i)%edge_theta_disp(m,1,q) = Ue(2*q)

        turbineArray(i)%elastic_twist(m,1,q) = Ut(q)
    enddo
    if (structure_timing) then
        call cpu_time(structure_t1)
        structure_update_local = structure_update_local +                     &
            max(structure_t1 - structure_stage_t0, 0._rprec)
    endif
enddo

deallocate(Kf, Ke, Mmat, Cf, Ce, Af, Ae, Kt, Ff, Fe, Ft, Uf, Ue, Ut)
deallocate(df, vf, accf, de, ve, acce, Ttension)

if (structure_timing) then
    call cpu_time(structure_t1)
    structure_total_local = max(structure_t1 - structure_t0, 0._rprec)
    atm_structure_timing_calls = atm_structure_timing_calls + 1
    atm_structure_time_total = atm_structure_time_total + structure_total_local
    atm_structure_time_assembly = atm_structure_time_assembly +               &
        structure_assembly_local
    atm_structure_time_solve = atm_structure_time_solve + structure_solve_local
    atm_structure_time_update = atm_structure_time_update + structure_update_local
    atm_structure_time_other = atm_structure_time_other +                     &
        max(structure_total_local - structure_assembly_local -                &
        structure_solve_local - structure_update_local, 0._rprec)
endif

end subroutine atm_solve_structure

#ifdef ENABLE_CUDA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
attributes(global) subroutine solve_linear_system_gpu_kernel(N, lda, A, B, X,  &
    info)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, value :: N, lda
real(rprec), device :: A(lda,*), B(*), X(*)
integer, device :: info(*)

integer :: i, j, k, max_idx
real(rprec) :: factor, pivot, temp

if (blockIdx%x /= 1 .or. threadIdx%x /= 1) return

info(1) = 0
do k = 1, N - 1
    max_idx = k
    pivot = abs(A(k,k))
    do i = k + 1, N
        if (abs(A(i,k)) > pivot) then
            max_idx = i
            pivot = abs(A(i,k))
        endif
    enddo
    if (pivot <= tiny(1._rprec)) then
        info(1) = k
        return
    endif
    if (max_idx /= k) then
        do j = k, N
            temp = A(k,j)
            A(k,j) = A(max_idx,j)
            A(max_idx,j) = temp
        enddo
        temp = B(k)
        B(k) = B(max_idx)
        B(max_idx) = temp
    endif
    do i = k + 1, N
        factor = A(i,k) / A(k,k)
        do j = k, N
            A(i,j) = A(i,j) - factor * A(k,j)
        enddo
        B(i) = B(i) - factor * B(k)
    enddo
enddo

if (abs(A(N,N)) <= tiny(1._rprec)) then
    info(1) = N
    return
endif
X(N) = B(N) / A(N,N)
do i = N - 1, 1, -1
    temp = B(i)
    do j = i + 1, N
        temp = temp - A(i,j) * X(j)
    enddo
    X(i) = temp / A(i,i)
enddo

end subroutine solve_linear_system_gpu_kernel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine solve_linear_system_gpu_dp(N, A, B, X)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! FP64 structural linear solve on the GPU. This mirrors the validated host
! Gaussian-elimination interface exactly, including row pivoting.
implicit none

integer, intent(in) :: N
real(rprec), intent(inout) :: A(N,N)
real(rprec), intent(inout) :: B(N)
real(rprec), intent(out) :: X(N)

real(rprec), device, save, allocatable :: A_d(:,:), B_d(:), X_d(:), work_d(:)
integer, device, save, allocatable :: info_d(:), ipiv_d(:)
type(cusolverDnHandle), save :: solver_handle
integer, save :: capacity = 0
integer, save :: work_capacity = 0
integer, save :: validation_count = 0
logical, save :: solver_initialized = .false.
integer :: istat, lwork, info_h
real(rprec), allocatable :: A_ref(:,:), B_ref(:), X_ref(:)
real(rprec) :: max_diff, ref_norm

if (.not. solver_initialized) then
    istat = cusolverDnCreate(solver_handle)
    if (istat /= 0) then
        print *, 'ATM structural cuSOLVER create failed: ', istat
        stop 1
    endif
    solver_initialized = .true.
endif

if (capacity < N) then
    if (allocated(A_d)) deallocate(A_d, B_d, X_d, info_d, ipiv_d)
    allocate(A_d(N,N), B_d(N), X_d(N), info_d(1), ipiv_d(N))
    capacity = N

    istat = cusolverDnDgetrf_bufferSize(solver_handle, capacity, capacity,     &
                                        A_d, capacity, lwork)
    if (istat /= 0) then
        print *, 'ATM structural cuSOLVER buffer query failed: ', istat
        stop 1
    endif
    if (work_capacity < lwork) then
        if (allocated(work_d)) deallocate(work_d)
        allocate(work_d(lwork))
        work_capacity = lwork
    endif
endif

A_d(1:N,1:N) = A(1:N,1:N)
B_d(1:N) = B(1:N)
info_d(1) = 0

if (atm_structure_gpu_direct_enabled()) then
    call solve_linear_system_gpu_kernel<<<1,1>>>(N, capacity, A_d, B_d, X_d,   &
                                                info_d)
    istat = cudaGetLastError()
    if (istat /= 0) then
        print *, 'ATM structural GPU linear solve launch failed: ', istat
        stop 1
    endif
    call atm_model_cuda_sync('ATM structural GPU direct linear solve')
    info_h = info_d(1)
    if (info_h /= 0) then
        print *, 'ATM structural GPU direct singular matrix, info=', info_h
        stop 1
    endif
    X(1:N) = X_d(1:N)
else
    istat = cusolverDnDgetrf(solver_handle, N, N, A_d, capacity, work_d,       &
                             ipiv_d, info_d(1))
    if (istat /= 0) then
        print *, 'ATM structural cuSOLVER Dgetrf failed: ', istat
        stop 1
    endif

    istat = cusolverDnDgetrs(solver_handle, CUBLAS_OP_N, N, 1, A_d, capacity,  &
                             ipiv_d, B_d, N, info_d(1))
    if (istat /= 0) then
        print *, 'ATM structural cuSOLVER Dgetrs failed: ', istat
        stop 1
    endif
    call atm_model_cuda_sync('ATM structural cuSOLVER solve')
    info_h = info_d(1)
    if (info_h /= 0) then
        print *, 'ATM structural cuSOLVER solve failed, info=', info_h
        stop 1
    endif
    X(1:N) = B_d(1:N)
endif

if (atm_structure_gpu_validate_enabled()) then
    allocate(A_ref(N,N), B_ref(N), X_ref(N))
    A_ref = A
    B_ref = B
    call solve_linear_system_dp(N, A_ref, B_ref, X_ref)
    max_diff = maxval(abs(X(1:N) - X_ref(1:N)))
    ref_norm = max(maxval(abs(X_ref(1:N))), tiny(1._rprec))
    validation_count = validation_count + 1
    if (validation_count <= 12 .or. max_diff / ref_norm > 1.0e-10_rprec) then
        write(*,*) 'ATM structural GPU solve validation N=', N,                &
            ' maxabs=', max_diff, ' rel=', max_diff / ref_norm
    endif
    deallocate(A_ref, B_ref, X_ref)
endif

end subroutine solve_linear_system_gpu_dp
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine solve_linear_system_banded_dp(N, half_band, A, B, X)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Gaussian elimination for the structural banded matrices. The flap/edge beam
! matrices have half-bandwidth 3, and the torsion matrix is tridiagonal. The row
! work is limited to the physical lower band, but each active row update keeps
! the dense solver's full upper-row arithmetic for stricter reproducibility.
implicit none

integer, intent(in) :: N, half_band
real(rprec), intent(inout) :: A(N,N)
real(rprec), intent(inout) :: B(N)
real(rprec), intent(out) :: X(N)

integer :: i, j, k, max_idx, row_max
real(rprec) :: factor, pivot, temp

do k = 1, N - 1
    row_max = min(N, k + half_band)
    max_idx = k
    pivot = abs(A(k,k))
    do i = k + 1, row_max
        if (abs(A(i,k)) > pivot) then
            max_idx = i
            pivot = abs(A(i,k))
        endif
    enddo
    if (pivot <= tiny(1._rprec)) call error('singular structural matrix')
    if (max_idx /= k) then
        do j = k, N
            temp = A(k,j)
            A(k,j) = A(max_idx,j)
            A(max_idx,j) = temp
        enddo
        temp = B(k)
        B(k) = B(max_idx)
        B(max_idx) = temp
    endif

    do i = k + 1, row_max
        if (A(i,k) /= 0._rprec) then
            factor = A(i,k) / A(k,k)
            A(i,k) = 0._rprec
            do j = k + 1, N
                A(i,j) = A(i,j) - factor * A(k,j)
            enddo
            B(i) = B(i) - factor * B(k)
        endif
    enddo
enddo

if (abs(A(N,N)) <= tiny(1._rprec)) call error('singular structural matrix')
X(N) = B(N) / A(N,N)
do i = N - 1, 1, -1
    temp = B(i)
    do j = i + 1, N
        temp = temp - A(i,j) * X(j)
    enddo
    X(i) = temp / A(i,i)
enddo

end subroutine solve_linear_system_banded_dp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine solve_linear_system_dp(N, A, B, X)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent(in) :: N
real(rprec), intent(inout) :: A(N,N)
real(rprec), intent(inout) :: B(N)
real(rprec), intent(out) :: X(N)

integer :: i, j, k, max_idx
real(rprec) :: factor, pivot, temp

do k = 1, N - 1
    max_idx = k
    pivot = abs(A(k,k))
    do i = k + 1, N
        if (abs(A(i,k)) > pivot) then
            max_idx = i
            pivot = abs(A(i,k))
        endif
    enddo
    if (pivot <= tiny(1._rprec)) call error('singular structural matrix')
    if (max_idx /= k) then
        do j = k, N
            temp = A(k,j)
            A(k,j) = A(max_idx,j)
            A(max_idx,j) = temp
        enddo
        temp = B(k)
        B(k) = B(max_idx)
        B(max_idx) = temp
    endif
    do i = k + 1, N
        factor = A(i,k) / A(k,k)
        do j = k, N
            A(i,j) = A(i,j) - factor * A(k,j)
        enddo
        B(i) = B(i) - factor * B(k)
    enddo
enddo

if (abs(A(N,N)) <= tiny(1._rprec)) call error('singular structural matrix')
X(N) = B(N) / A(N,N)
do i = N - 1, 1, -1
    temp = B(i)
    do j = i + 1, N
        temp = temp - A(i,j) * X(j)
    enddo
    X(i) = temp / A(i,i)
enddo

end subroutine solve_linear_system_dp

!~ !-------------------------------------------------------------------------------
!~ function atm_convoluteForce(i,m,n,q,xyz)
!~ !-------------------------------------------------------------------------------
!~ ! This subroutine will convolute the body forces onto a point xyz
!~ integer, intent(in) :: i,m,n,q
!~ ! i - turbineTypeArray
!~ ! n - numAnnulusSections
!~ ! q - numBladePoints
!~ ! m - numBl
!~ real(rprec), intent(in) :: xyz(3)    ! Point onto which to convloute the force
!~ real(rprec) :: Force(3)   ! The blade force to be convoluted
!~ real(rprec) :: dis                ! Distance onto which convolute the force
!~ real(rprec) :: atm_convoluteForce(3)    ! The local velocity at this point
!~ real(rprec) :: kernel                ! Gaussian dsitribution value
!~
!~ ! Distance from the point of the force to the point where it is being convoluted
!~ dis=distance(xyz,turbineArray(i) % bladepoints(m,n,q,:))
!~
!~ ! The force which is being convoluted
!~ Force=turbineArray(i) % bladeForces(m,n,q,:)
!~
!~ ! The value of the kernel. This is the actual smoothing function
!~ kernel=exp(-(dis/turbineArray(i) % epsilon)**2._rprec) /                       &
!~ ((turbineArray(i) % epsilon**3._rprec)*(pi**1.5_rprec))
!~
!~ ! The force times the kernel will give the force/unitVolume
!~ atm_convoluteForce = Force * kernel
!~
!~ return
!~ end function atm_convoluteForce



end module actuator_turbine_model
