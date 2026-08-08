!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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

module level_set_base
use types, only : rp => rprec
use types, only : rprec
use param, only : ld, ny, nz, dx, lbz, nproc
use sgs_param, only : SGS_MODEL_LAGRANGE_SIMILARITY,                       &
                      SGS_MODEL_LAGRANGE_SCALE_DEP
use messages, only : error
implicit none

save

public
private :: rp, ld, ny, nz, dx, lbz, nproc

type :: level_set_halo_requirements_t
  integer :: phi_top = 0
  integer :: phi_bottom = 0
  integer :: velocity_top = 0
  integer :: velocity_bottom = 0
  integer :: stress_top = 0
  integer :: stress_bottom = 0
  integer :: fmm_top = 0
  integer :: fmm_bottom = 0
end type level_set_halo_requirements_t

!private
!public :: phi

!logical, parameter :: global_CD_calc = .true. ! Compute global CD based on inflow velocity
logical :: global_CA_calc = .false. ! Compute global CA based on inflow velocity
integer :: global_CA_nskip = 10     ! Number of time steps to skip between global CA writes

!logical, parameter :: vel_BC = .false.
logical :: vel_BC = .false. !--means we are forcing velocity for
                            !  level set BC
                            !  (default = .false.)
! logical, parameter :: use_log_profile = .false.       !  (default = .false.)
! logical, parameter :: use_enforce_un = .false.        !  (default = .false.)
! logical, parameter :: physBC = .true.                 !  (default = .true.)
! logical, parameter :: use_smooth_tau = .true.         !  (default = .true.)
! logical, parameter :: use_extrap_tau_log = .false.    !  (default = .false.)
! logical, parameter :: use_extrap_tau_simple = .true.  !  (default = .true.)
! logical, parameter :: use_modify_dutdn = .false.  !--only works w/interp_tau; not MPI compliant
!                                                   !--wont work w/extra_tau_log
!                                                   !  (default = .false.)
logical :: use_log_profile = .false.       !  (default = .false.)
logical :: use_enforce_un = .false.        !  (default = .false.)
logical :: physBC = .true.                 !  (default = .true.)
logical :: use_smooth_tau = .true.         !  (default = .true.)
logical :: use_extrap_tau_log = .false.    !  (default = .false.)
logical :: use_extrap_tau_simple = .true.  !  (default = .true.)
logical :: use_modify_dutdn = .false.  !--only works w/interp_tau; not MPI compliant
                                                  !--wont work w/extra_tau_log
                                                  !  (default = .false.)

! ! Enables scale dependent Cs evaluations (not dynamic evaluation)
! ! Used when sgs_model=4 in param module
! logical, parameter :: lag_dyn_modify_beta = .true.

! Enables scale dependent Cs evaluations (not dynamic evaluation)
! Used when sgs_model=4 in param module
logical :: lag_dyn_modify_beta = .true.

! ! Configures the mode in which SOR smoothing is applied in the IB
! ! 'xy' may be safely used in most cases (must be used for MPI cases)
! ! '3d' not MPI compliant
! character (*), parameter :: smooth_mode = 'xy'  !--'xy', '3d'

! Configures the mode in which SOR smoothing is applied in the IB
! 'xy' may be safely used in most cases (must be used for MPI cases)
! '3d' not MPI compliant
character(25) :: smooth_mode = 'xy'  !--'xy', '3d'

!real (rp), parameter :: zo_level_set = 0.0001_rp !--nondimensional roughness length of surface
real (rp) :: zo_level_set = 0.0001_rp !--nondimensional roughness length of surface

logical :: phi_cutoff_is_set = .false.
logical :: phi_0_is_set = .false.


!real (rp) :: phi(ld, ny, lbz:nz)
real(rp), allocatable, target, dimension(:,:,:) :: phi

!--Extended vertical overlap fields used by immersed-surface interpolation.
!--They live with the Level Set data model so both the CPU implementation and
!--the GPU interpolation kernels share one owner and one allocation.
real(rp), allocatable, dimension(:,:,:) :: phitop, phibot
real(rp), allocatable, dimension(:,:,:) :: utop, vtop, wtop
real(rp), allocatable, dimension(:,:,:) :: ubot, vbot, wbot
real(rp), allocatable, dimension(:,:,:) :: txxtop, txytop, txztop,          &
                                           tyytop, tyztop, tzztop
real(rp), allocatable, dimension(:,:,:) :: txxbot, txybot, txzbot,          &
                                           tyybot, tyzbot, tzzbot
real(rp), allocatable, dimension(:,:,:) :: FMMbot, FMMtop

logical :: use_trees

#ifdef PPMPI
  ! Make sure all values (top and bottom) are less than Nz
  integer :: nphitop = 3
  integer :: nphibot = 2
  integer :: nveltop = 1
  integer :: nvelbot = 1
  integer :: ntautop = 3
  integer :: ntaubot = 2
  integer :: nFMMtop = 1
  integer :: nFMMbot = 1
#else
  integer, parameter :: nphitop = 0
  integer, parameter :: nphibot = 0
  integer, parameter :: nveltop = 0
  integer, parameter :: nvelbot = 0
  integer, parameter :: ntautop = 0
  integer, parameter :: ntaubot = 0
  integer, parameter :: nFMMtop = 0
  integer, parameter :: nFMMbot = 0
#endif


contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_required_halos(sgs_enabled, sgs_model, required)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical, intent(in) :: sgs_enabled
integer, intent(in) :: sgs_model
type(level_set_halo_requirements_t), intent(out) :: required

required = level_set_halo_requirements_t()

! Direct log-profile stress treatment does not interpolate stress across rank
! boundaries. The optional velocity boundary treatment still interpolates the
! resolved velocity, including when it enforces a log profile.
if (.not. use_log_profile .or. vel_BC) then
  required%velocity_top = 1
  required%velocity_bottom = 1
end if

! Simple stress extrapolation reflects points across the interface. The
! historical default depths are the minimum supported stencil contract.
if (.not. use_log_profile .and. .not. use_extrap_tau_log .and.              &
    use_extrap_tau_simple) then
  required%phi_top = 3
  required%phi_bottom = 2
  required%stress_top = 3
  required%stress_bottom = 2
end if

if (sgs_enabled .and. (sgs_model == SGS_MODEL_LAGRANGE_SIMILARITY .or.      &
    sgs_model == SGS_MODEL_LAGRANGE_SCALE_DEP)) then
  required%fmm_top = 1
  required%fmm_bottom = 1
end if
end subroutine level_set_required_halos

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_validate_config(sgs_enabled, sgs_model)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logical, intent(in) :: sgs_enabled
integer, intent(in) :: sgs_model
type(level_set_halo_requirements_t) :: required
character(*), parameter :: sub_name = 'level_set_base.level_set_validate_config'
logical :: active_log_law

if (trim(smooth_mode) /= 'xy' .and. trim(smooth_mode) /= '3d') then
  call error(sub_name, 'smooth_mode must be exactly "xy" or "3d": ' //       &
      trim(smooth_mode))
end if

#ifdef PPMPI
if (trim(smooth_mode) == '3d') then
  call error(sub_name, 'smooth_mode="3d" requires a non-MPI build')
end if
#endif

if (nproc > 1) then
  if (use_extrap_tau_log) then
    call error(sub_name, 'use_extrap_tau_log requires a single-rank run')
  end if
  if (.not. use_extrap_tau_log .and. .not. use_extrap_tau_simple) then
    call error(sub_name, 'legacy stress extrapolation requires a single-rank run')
  end if
  if (use_modify_dutdn) then
    call error(sub_name, 'use_modify_dutdn requires a single-rank run')
  end if
end if

! use_log_profile selects the direct stress treatment even when optional
! velocity forcing is disabled, so it is independently a log-law path.
active_log_law = use_extrap_tau_log .or. use_log_profile
if (active_log_law .and. zo_level_set <= 0._rprec) then
  call error(sub_name, 'active Level Set log-law paths require zo_level_set > 0')
end if

if (global_CA_calc .and. global_CA_nskip <= 0) then
  call error(sub_name, 'global_CA_nskip must be positive when global_CA_calc is enabled')
end if

#ifdef PPMPI
call level_set_required_halos(sgs_enabled, sgs_model, required)
call validate_halo('nphitop', nphitop, required%phi_top)
call validate_halo('nphibot', nphibot, required%phi_bottom)
call validate_halo('nveltop', nveltop, required%velocity_top)
call validate_halo('nvelbot', nvelbot, required%velocity_bottom)
call validate_halo('ntautop', ntautop, required%stress_top)
call validate_halo('ntaubot', ntaubot, required%stress_bottom)
call validate_halo('nFMMtop', nFMMtop, required%fmm_top)
call validate_halo('nFMMbot', nFMMbot, required%fmm_bottom)
#endif

contains

subroutine validate_halo(name, actual, minimum)
character(*), intent(in) :: name
integer, intent(in) :: actual, minimum

if (actual < 0 .or. actual >= nz) then
  call error(sub_name, trim(name) // ' must satisfy 0 <= depth < Nz', actual)
end if
if (actual < minimum) then
  call error(sub_name, trim(name) // ' is smaller than the active stencil minimum', &
      actual, ' required=', minimum)
end if
end subroutine validate_halo

end subroutine level_set_validate_config

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine level_set_base_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine initializes all arrays defined in level_set_base
!
implicit none

allocate( phi( ld, ny, lbz:nz ) )

return
end subroutine level_set_base_init


end module level_set_base
