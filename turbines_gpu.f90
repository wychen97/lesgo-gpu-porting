!*******************************************************************************
module turbines_gpu
!*******************************************************************************
! GPU helper module for the actuator-disk turbine path in turbines.f90.
!
! Ownership map:
!   - turbines_cuda_enabled: compile-time/runtime availability query
!   - turbines_interp_w_to_uv_gpu: interpolate w-grid velocity onto uv nodes
!     for turbine sampling
!
! This is separate from the actuator-line/ATM GPU path, which lives in
! actuator_turbine_model.f90 and atm_lesgo_interface.f90.
use types, only : rprec
use param, only : ld, ny, nz, lbz
#ifdef PPMPI
use param, only : coord, nproc
#endif
implicit none

private
public :: turbines_cuda_enabled, turbines_interp_w_to_uv_gpu

contains

!*******************************************************************************
logical function turbines_cuda_enabled()
!*******************************************************************************
implicit none

turbines_cuda_enabled = .false.

end function turbines_cuda_enabled

!*******************************************************************************
subroutine turbines_interp_w_to_uv_gpu(w, w_uv)
!*******************************************************************************
real(rprec), intent(in) :: w(ld,ny,lbz:nz)
real(rprec), intent(out) :: w_uv(ld,ny,lbz:nz)


end subroutine turbines_interp_w_to_uv_gpu

end module turbines_gpu
