!!
!!  Copyright (C) 2016  Johns Hopkins University
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
module iwmles
!*******************************************************************************
!
! This module contains the procedures for using the integral wall model.
! See Yang et al. 2015 for more details
!
! Navigation map:
!   - lifecycle: iwm_init, iwm_finalize, checkpoint/read_checkpoint
!   - timestep entry point: iwm_wallstress
!   - filtered wall-model inputs: iwm_calc_lhs
!   - nonlinear wall solve: iwm_slv and iwm_calc_wallstress
!   - diagnostics: iwm_monitor
!
! Wall-model arrays are point-local in (i,j).  Updates to filtered friction
! velocity and filter timescale must use the same (i,j) location on both sides
! of the assignment.
!
! Ownership map:
!   - this module owns wall-surface state such as iwm_tauw*, iwm_flt_us,
!     iwm_tR, iwm_Dz, and integrated wall-model profiles.
!   - LES velocity/pressure fields are sampled inputs; this module should not
!     make full LES fields host-authoritative during the GPU timestep.
!   - ENABLE_CUDA keeps the legacy managed-memory path, while PPSGS_GPU uses
!     OpenACC-present arrays in the optimized non-LVLSET branch.

use types, only : rprec
#ifdef ENABLE_CUDA
use cudafor
#endif

implicit none

private
public iwm_wallstress, iwm_init, iwm_finalize,                                 &
    iwm_checkpoint, iwm_read_checkpoint, iwm_lhs_update_due

! u_tau,x  u_tau,y
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:), allocatable :: iwm_utx, iwm_uty
#else
real(rprec), dimension(:,:), allocatable :: iwm_utx, iwm_uty
#endif
! tau_wall,x tau_wall,y
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:), allocatable :: iwm_tauwx, iwm_tauwy
#else
real(rprec), dimension(:,:), allocatable :: iwm_tauwx, iwm_tauwy
#endif
! filtered tangential velocity, current and previous
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,:), allocatable :: iwm_flt_tagvel,         &
    iwm_flt_tagvel_m
! filtered pressure
real(rprec), managed, dimension(:,:), allocatable :: iwm_flt_p
#else
real(rprec), dimension(:,:,:), allocatable :: iwm_flt_tagvel, iwm_flt_tagvel_m
! filtered pressure
real(rprec), dimension(:,:), allocatable :: iwm_flt_p
#endif
! direction x
integer :: iwm_dirx = 1
! direction y
integer :: iwm_diry = 2
! dimension of a surface, wall model always deal with 2D surfaces
! (because the world is 3D)
integer :: iwm_DN  = 2
! integrated profiles, current and previous
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,:), allocatable :: iwm_inte, iwm_inte_m
#else
real(rprec), dimension(:,:,:), allocatable :: iwm_inte, iwm_inte_m
#endif
integer :: iwm_Lu = 1   ! index for integral of u
integer :: iwm_Luu = 2  ! index for integral of uu
integer :: iwm_Lv = 3   ! etc.
integer :: iwm_Lvv = 4
integer :: iwm_Luv = 5
! the total number of integrals that need to be calculated
integer :: iwm_LN  = 5
! unsteady, convective, pressure gradient
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,:), allocatable :: iwm_unsdy, iwm_conv,    &
    iwm_PrsGrad
! turbulent diffusion, LHS
real(rprec), managed, dimension(:,:,:), allocatable :: iwm_diff, iwm_LHS
#else
real(rprec), dimension(:,:,:), allocatable :: iwm_unsdy, iwm_conv, iwm_PrsGrad
! turbulent diffusion, LHS
real(rprec), dimension(:,:,:), allocatable :: iwm_diff, iwm_LHS
#endif
! dudz at z=dz/2, dudz at z=zo
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:,:), allocatable :: iwm_dudzT
real(rprec), managed, dimension(:,:,:), allocatable :: iwm_dudzB
#else
real(rprec), dimension(:,:,:), allocatable :: iwm_dudzT, iwm_dudzB
#endif
! filtered friction velocity, filtering time scale
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:), allocatable :: iwm_flt_us, iwm_tR
#else
real(rprec), dimension(:,:), allocatable :: iwm_flt_us, iwm_tR
#endif
! HALF cell height, zo, linear correction in x, y directions
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:,:), allocatable :: iwm_Dz, iwm_z0, iwm_Ax,   &
    iwm_Ay
real(rprec), managed, dimension(:,:), allocatable :: iwm_u_inst, iwm_v_inst,   &
    iwm_w_inst, iwm_p_inst
#else
real(rprec), dimension(:,:), allocatable :: iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay
#if defined(PPSGS_GPU)
real(rprec), dimension(:,:), allocatable :: iwm_u_inst, iwm_v_inst,            &
    iwm_w_inst, iwm_p_inst
#endif
#endif

! number of time steps to skip between wall stress calculations
integer :: iwm_ntime_skip = 5
! time step size seen by the wall model
real(rprec) :: iwm_dt

contains

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function iwm_cuda_enabled()
!*******************************************************************************
implicit none

iwm_cuda_enabled = .true.

end function iwm_cuda_enabled

!*******************************************************************************
subroutine iwm_cuda_sync(where)
!*******************************************************************************
implicit none

character(*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= cudaSuccess) then
    print *, 'IWM CUDA sync failure at ', trim(where), ': ', istat
    stop 1
end if
istat = cudaGetLastError()
if (istat /= cudaSuccess) then
    print *, 'IWM CUDA kernel failure at ', trim(where), ': ', istat
    stop 1
end if

end subroutine iwm_cuda_sync

#endif

!*******************************************************************************
logical function iwm_lhs_update_due()
!*******************************************************************************
! The expensive host-bridge input planes are only consumed when iwm_calc_lhs()
! runs.  On intervening steps iwm_wallstress() just reapplies the stored wall
! stress values, so callers can skip staging u/v/w/p from the device.
use param, only : jt
implicit none

iwm_lhs_update_due = (mod(jt, iwm_ntime_skip) == 0)

end function iwm_lhs_update_due

!*******************************************************************************
subroutine iwm_wallstress
!*******************************************************************************
use param, only : jt, nx, ny, dt
use sim_param , only : dudz, dvdz, txz, tyz
implicit none

integer :: iwm_i, iwm_j

! Calculate the time step used in the integral wall model
!! DO NOT USE iwm_ntime_skip=1 !! !! this number is hard coded to prevent any
!! mis-use...
if (mod(jt,iwm_ntime_skip)==1) then
    iwm_dt = dt
else
    iwm_dt = iwm_dt+dt
end if

! Compute the wall stress
if(mod(jt,iwm_ntime_skip)==0) then
    ! gather flow status, update the integrated unsteady term, convective term,
    ! turbulent diffusion term etc.
    call iwm_calc_lhs()
    ! the subroutine to calculate wall stress
    call iwm_calc_wallstress()
    ! this is to monitor any quantity from the iwm, useful debugging tool
    call iwm_monitor()
end if

! Imposing txz, tyz, dudz, dvdz every time step even iwm_* are not computed
! every time step.
#ifdef ENABLE_CUDA
if (iwm_cuda_enabled()) then
    !$cuf kernel do(2) <<<*,*>>>
    do iwm_j = 1, ny
    do iwm_i = 1, nx
        txz(iwm_i,iwm_j,1) = -iwm_tauwx(iwm_i,iwm_j)
        tyz(iwm_i,iwm_j,1) = -iwm_tauwy(iwm_i,iwm_j)
        dudz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
        dvdz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_diry)
    end do
    end do
    call iwm_cuda_sync('iwm_wallstress apply')
else
#endif
#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
!$acc parallel loop collapse(2) default(present) async(1)
do iwm_j = 1, ny
do iwm_i = 1, nx
    txz(iwm_i,iwm_j,1) = -iwm_tauwx(iwm_i,iwm_j)
    tyz(iwm_i,iwm_j,1) = -iwm_tauwy(iwm_i,iwm_j)
    dudz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
    dvdz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_diry)
end do
end do
return
#endif
do iwm_i = 1, nx
do iwm_j = 1, ny
    ! wall stress, use the value calculated in iwm, note the negative sign
    txz(iwm_i,iwm_j,1) = -iwm_tauwx(iwm_i,iwm_j)
    tyz(iwm_i,iwm_j,1) = -iwm_tauwy(iwm_i,iwm_j)

    ! Use wall gradients computed by the integral wall model rather than
    ! equilibrium estimates.  The positive sign follows the iwm_dudzT
    ! convention used by the wall-stress solve.
    dudz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
    dvdz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_diry)
end do
end do
#ifdef ENABLE_CUDA
end if
#endif

end subroutine iwm_wallstress

!*******************************************************************************
subroutine iwm_init
!*******************************************************************************
!
! This subroutine allocates memory and initializes everything with plug flow
! conditions
!
use types, only : rprec
use param, only : nx, ny, ld, dz, vonk, zo, cfl, L_x

implicit none

real(rprec) :: usinit, uinit, vinit, Dzp

! initial value for us (the friction velocity)
usinit= 1._rprec
! initial value for the x-velocity at first grid point
uinit = usinit/vonk*log(dz/2._rprec/zo)
! initial value for the y-velocity at first grid point
vinit = 0._rprec
! at the height of the first grid point.
Dzp=dz/2._rprec

! us in x, y directions
allocate(iwm_utx(nx,ny))
allocate(iwm_uty(nx,ny))
iwm_utx = usinit
iwm_uty = 0._rprec

! wall stress in x, y directions
allocate(iwm_tauwx(nx,ny))
allocate(iwm_tauwy(nx,ny))
iwm_tauwx = usinit*usinit
iwm_tauwy = 0._rprec

! filitered velocity at the first grid point in x, y directions
allocate(iwm_flt_tagvel  (nx,ny,iwm_DN))
allocate(iwm_flt_tagvel_m(nx,ny,iwm_DN))
iwm_flt_tagvel  (:,:,iwm_dirx) = uinit
iwm_flt_tagvel  (:,:,iwm_diry) = vinit
iwm_flt_tagvel_m(:,:,iwm_dirx) = uinit
iwm_flt_tagvel_m(:,:,iwm_diry) = vinit

! pressure at first grid point
allocate(iwm_flt_p(nx,ny))
iwm_flt_p = 0._rprec

! integrals of Lu, Lv, etc.
allocate(iwm_inte  (nx,ny,iwm_LN))
allocate(iwm_inte_m(nx,ny,iwm_LN))
iwm_inte(:,:,iwm_Lu) = uinit*Dzp
iwm_inte(:,:,iwm_Lv) = 0._rprec
iwm_inte(:,:,iwm_Luu) = uinit*uinit*Dzp
iwm_inte(:,:,iwm_Lvv) = 0._rprec
iwm_inte(:,:,iwm_Luv) = 0._rprec
iwm_inte_m(:,:,iwm_Lu) = uinit*Dzp
iwm_inte_m(:,:,iwm_Lv) = 0._rprec
iwm_inte_m(:,:,iwm_Luu) = uinit*uinit*Dzp
iwm_inte_m(:,:,iwm_Lvv) = 0._rprec
iwm_inte_m(:,:,iwm_Luv) = 0._rprec

! each term in the integral equation and top/bottom derivatives
allocate(iwm_unsdy  (nx,ny,iWM_DN))
allocate(iwm_conv   (nx,ny,iWM_DN))
allocate(iwm_PrsGrad(nx,ny,iWM_DN))
allocate(iwm_diff   (nx,ny,iwm_DN))
allocate(iwm_LHS    (nx,ny,iWM_DN))
allocate(iwm_dudzT  (nx,ny,iwm_DN))
allocate(iwm_dudzB  (nx,ny,iwm_DN))
iWM_unsdy   = 0._rprec
iWM_conv    = 0._rprec
iWM_PrsGrad = 0._rprec
iwm_diff    = 0._rprec
iWM_LHS     = -uinit*Dzp
iwm_dudzT(:,:,iwm_dirx) = usinit/vonk/Dzp
iwm_dudzT(:,:,iwm_diry) = 0._rprec
iwm_dudzB(:,:,iwm_dirx) = usinit/vonk/zo
iwm_dudzB(:,:,iwm_diry) = 0._rprec

! filtered friction velocity and the filtering time scale, tR<1
allocate(iwm_flt_us(nx,ny))
allocate(iwm_tR    (nx,ny))
iwm_flt_us = usinit
iwm_tR = (cfl*L_x/nx/uinit)/(dz/2._rprec/vonk/usinit)

! cell height and imposed roughness length
allocate(iwm_Dz(nx,ny))
allocate(iwm_z0(nx,ny))
iWM_Dz = dz/2._rprec
iWM_z0 = zo !we leave the possibility of zo as a function of x-y

! linear correction to the log profile
allocate(iwm_Ax(nx,ny))
allocate(iwm_Ay(nx,ny))
iwm_Ax = 0._rprec
iwm_Ay = 0._rprec

#if defined(ENABLE_CUDA) || (defined(PPSGS_GPU))
allocate(iwm_u_inst(ld,ny))
allocate(iwm_v_inst(ld,ny))
allocate(iwm_w_inst(ld,ny))
allocate(iwm_p_inst(ld,ny))
iwm_u_inst = 0._rprec
iwm_v_inst = 0._rprec
iwm_w_inst = 0._rprec
iwm_p_inst = 0._rprec
#endif

! time step seen by the iwm
iwm_dt=iwm_ntime_skip*cfl*L_x/nx/uinit

#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
!$acc enter data copyin(iwm_utx, iwm_uty, iwm_tauwx, iwm_tauwy)
!$acc enter data copyin(iwm_flt_tagvel, iwm_flt_tagvel_m, iwm_flt_p)
!$acc enter data copyin(iwm_inte, iwm_inte_m, iwm_unsdy, iwm_conv)
!$acc enter data copyin(iwm_PrsGrad, iwm_diff, iwm_LHS)
!$acc enter data copyin(iwm_dudzT, iwm_dudzB, iwm_flt_us, iwm_tR)
!$acc enter data copyin(iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay)
!$acc enter data copyin(iwm_u_inst, iwm_v_inst, iwm_w_inst, iwm_p_inst)
#endif

end subroutine iwm_init

!*******************************************************************************
subroutine iwm_finalize
!*******************************************************************************
!
! This subroutine deallocates memory used for iwm
!
implicit none

#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
!$acc exit data delete(iwm_utx, iwm_uty, iwm_tauwx, iwm_tauwy)
!$acc exit data delete(iwm_flt_tagvel, iwm_flt_tagvel_m, iwm_flt_p)
!$acc exit data delete(iwm_inte, iwm_inte_m, iwm_unsdy, iwm_conv)
!$acc exit data delete(iwm_PrsGrad, iwm_diff, iwm_LHS)
!$acc exit data delete(iwm_dudzT, iwm_dudzB, iwm_flt_us, iwm_tR)
!$acc exit data delete(iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay)
!$acc exit data delete(iwm_u_inst, iwm_v_inst, iwm_w_inst, iwm_p_inst)
#endif

deallocate(iwm_utx)
deallocate(iwm_uty)

deallocate(iwm_tauwx)
deallocate(iwm_tauwy)

deallocate(iwm_flt_tagvel  )
deallocate(iwm_flt_tagvel_m)

deallocate(iwm_inte  )
deallocate(iwm_inte_m)

deallocate(iwm_unsdy  )
deallocate(iwm_conv   )
deallocate(iwm_PrsGrad)
deallocate(iwm_diff   )
deallocate(iwm_LHS    )
deallocate(iwm_dudzT  )
deallocate(iwm_dudzB  )

deallocate(iwm_flt_us)
deallocate(iwm_tR    )

deallocate(iwm_Dz)
deallocate(iwm_z0)
deallocate(iwm_Ax)
deallocate(iwm_Ay)
#if defined(ENABLE_CUDA) || (defined(PPSGS_GPU))
if (allocated(iwm_u_inst)) deallocate(iwm_u_inst)
if (allocated(iwm_v_inst)) deallocate(iwm_v_inst)
if (allocated(iwm_w_inst)) deallocate(iwm_w_inst)
if (allocated(iwm_p_inst)) deallocate(iwm_p_inst)
#endif

end subroutine iwm_finalize


!*******************************************************************************
subroutine iwm_calc_lhs()
!*******************************************************************************
!
! Ths subroutine calculates the left hand side of the iwm system.
!
use grid_m, only : grid
use types,only : rprec
use param,only : nx,ny,dx,dy,ld
use sim_param,only : u,v,w,p
use test_filtermodule
implicit none

! Wrapped horizontal-neighbor indices from the grid module.
integer, pointer, dimension(:) :: autowrap_i, autowrap_j
integer :: iwm_i,iwm_j
#if defined(ENABLE_CUDA) || (defined(PPSGS_GPU))
integer :: ip, im, jp, jm
#endif
! the instantaneous field
real(rprec), dimension(ld,ny) :: u_inst, v_inst, w_inst, p_inst
! the mean pressure at first grid point
real(rprec) :: p_bar
! Scratch derivatives of integrated wall-model moments such as dLu/dx, dLv/dx.
real(rprec) :: Luux, Luvx, Luvy, Lvvy, Lux, Lvy
real(rprec) :: phip, phim
#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
real(rprec) :: iwm_dt_l
#endif

nullify(autowrap_i, autowrap_j)
autowrap_i => grid % autowrap_i
autowrap_j => grid % autowrap_j

#ifdef ENABLE_CUDA
if (iwm_cuda_enabled()) then
    p_bar = 0._rprec

    !$cuf kernel do(2) <<<*,*>>> reduction(+:p_bar)
    do iwm_j = 1, ny
    do iwm_i = 1, nx
        iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx) =                              &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
        iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry) =                              &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
        iwm_u_inst(iwm_i,iwm_j) = u(iwm_i,iwm_j,1)
        iwm_v_inst(iwm_i,iwm_j) = v(iwm_i,iwm_j,1)
        iwm_w_inst(iwm_i,iwm_j) = w(iwm_i,iwm_j,2)*0.25_rprec
        iwm_p_inst(iwm_i,iwm_j) = p(iwm_i,iwm_j,1)                            &
            - 0.5_rprec*(iwm_u_inst(iwm_i,iwm_j)*iwm_u_inst(iwm_i,iwm_j)      &
            + iwm_v_inst(iwm_i,iwm_j)*iwm_v_inst(iwm_i,iwm_j)                 &
            + iwm_w_inst(iwm_i,iwm_j)*iwm_w_inst(iwm_i,iwm_j))
        p_bar = p_bar + iwm_p_inst(iwm_i,iwm_j)
    end do
    end do
    call iwm_cuda_sync('iwm_calc_lhs pack')

    p_bar = p_bar / real(nx*ny, rprec)
    !$cuf kernel do(2) <<<*,*>>>
    do iwm_j = 1, ny
    do iwm_i = 1, nx
        iwm_p_inst(iwm_i,iwm_j) = iwm_p_inst(iwm_i,iwm_j) - p_bar
    end do
    end do
    call iwm_cuda_sync('iwm_calc_lhs pressure mean')

    call test_filter_3(iwm_u_inst, iwm_v_inst, iwm_w_inst)
    call test_filter(iwm_p_inst)

    !$cuf kernel do(2) <<<*,*>>>
    do iwm_j = 1, ny
    do iwm_i = 1, nx
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx) =                                &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                              &
            * (1._rprec-iwm_tR(iwm_i,iwm_j))                                  &
            + iwm_u_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) =                                &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                              &
            * (1._rprec-iwm_tR(iwm_i,iwm_j))                                  &
            + iwm_v_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
        iwm_flt_p(iwm_i,iwm_j) = iwm_flt_p(iwm_i,iwm_j)                       &
            * (1._rprec-iwm_tR(iwm_i,iwm_j))                                  &
            + iwm_p_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)

        iwm_unsdy(iwm_i,iwm_j,iwm_dirx) =                                     &
            (iwm_inte(iwm_i,iwm_j,iwm_Lu)-iwm_inte_m(iwm_i,iwm_j,iwm_Lu))     &
            / iwm_dt
        iwm_unsdy(iwm_i,iwm_j,iwm_diry) =                                     &
            (iwm_inte(iwm_i,iwm_j,iwm_Lv)-iwm_inte_m(iwm_i,iwm_j,iwm_Lv))     &
            / iwm_dt

        ip = iwm_i + 1
        im = iwm_i - 1
        jp = iwm_j + 1
        jm = iwm_j - 1
        if (ip > nx) ip = 1
        if (im < 1) im = nx
        if (jp > ny) jp = 1
        if (jm < 1) jm = ny

        Luux = (iwm_inte(ip,iwm_j,iwm_Luu) - iwm_inte(im,iwm_j,iwm_Luu))      &
            / dx / 2._rprec
        Luvy = (iwm_inte(iwm_i,jp,iwm_Luv) - iwm_inte(iwm_i,jm,iwm_Luv))      &
            / dy / 2._rprec
        Luvx = (iwm_inte(ip,iwm_j,iwm_Luv) - iwm_inte(im,iwm_j,iwm_Luv))      &
            / dx / 2._rprec
        Lvvy = (iwm_inte(iwm_i,jp,iwm_Lvv) - iwm_inte(iwm_i,jm,iwm_Lvv))      &
            / dy / 2._rprec
        Lux = (iwm_inte(ip,iwm_j,iwm_Lu) - iwm_inte(im,iwm_j,iwm_Lu))         &
            / dx / 2._rprec
        Lvy = (iwm_inte(iwm_i,jp,iwm_Lv) - iwm_inte(iwm_i,jm,iwm_Lv))         &
            / dy / 2._rprec

        iwm_conv(iwm_i,iwm_j,iwm_dirx) = Luux + Luvy                          &
            - iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx)*(Lux+Lvy)
        iwm_conv(iwm_i,iwm_j,iwm_diry) = Luvx + Lvvy                          &
            - iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry)*(Lux+Lvy)

        iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx) =                                   &
            (iwm_flt_p(ip,iwm_j)-iwm_flt_p(im,iwm_j))/dx/2._rprec             &
            * iwm_Dz(iwm_i,iwm_j) - iwm_Dz(iwm_i,iwm_j)
        iwm_PrsGrad(iwm_i,iwm_j,iwm_diry) =                                   &
            (iwm_flt_p(iwm_i,jp)-iwm_flt_p(iwm_i,jm))/dy/2._rprec             &
            * iwm_Dz(iwm_i,iwm_j)

        iwm_lhs(iwm_i,iwm_j,iwm_dirx) = -iwm_inte(iwm_i,iwm_j,iwm_Lu)         &
            + iwm_dt*(iwm_conv(iwm_i,iwm_j,iwm_dirx)                          &
            + iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)                               &
            - iwm_diff(iwm_i,iwm_j,iwm_dirx))
        iwm_lhs(iwm_i,iwm_j,iwm_diry) = -iwm_inte(iwm_i,iwm_j,iwm_Lv)         &
            + iwm_dt*(iwm_conv(iwm_i,iwm_j,iwm_diry)                          &
            + iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)                               &
            - iwm_diff(iwm_i,iwm_j,iwm_diry))
    end do
    end do
    call iwm_cuda_sync('iwm_calc_lhs final')
    nullify(autowrap_i, autowrap_j)
    return
end if
#endif

#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
iwm_dt_l = iwm_dt
p_bar = 0._rprec

! Match the host plane assignment before filtering: the FFT plane has ld rows,
! even though only the physical 1:nx rows contribute to the IWM mean terms.
!$acc parallel loop collapse(2) default(present) async(1)
do iwm_j = 1, ny
do iwm_i = 1, ld
    iwm_u_inst(iwm_i,iwm_j) = u(iwm_i,iwm_j,1)
    iwm_v_inst(iwm_i,iwm_j) = v(iwm_i,iwm_j,1)
    iwm_w_inst(iwm_i,iwm_j) = w(iwm_i,iwm_j,2)*0.25_rprec
    iwm_p_inst(iwm_i,iwm_j) = p(iwm_i,iwm_j,1)
end do
end do

!$acc parallel loop collapse(2) default(present) reduction(+:p_bar) async(1)
do iwm_j = 1, ny
do iwm_i = 1, nx
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx) =                                  &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry) =                                  &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
    iwm_p_inst(iwm_i,iwm_j) = iwm_p_inst(iwm_i,iwm_j)                          &
        - 0.5_rprec*(iwm_u_inst(iwm_i,iwm_j)*iwm_u_inst(iwm_i,iwm_j)          &
        + iwm_v_inst(iwm_i,iwm_j)*iwm_v_inst(iwm_i,iwm_j)                     &
        + iwm_w_inst(iwm_i,iwm_j)*iwm_w_inst(iwm_i,iwm_j))
    p_bar = p_bar + iwm_p_inst(iwm_i,iwm_j)
end do
end do
!$acc wait(1)

p_bar = p_bar / real(nx*ny, rprec)
!$acc parallel loop collapse(2) default(present) async(1)
do iwm_j = 1, ny
do iwm_i = 1, nx
    iwm_p_inst(iwm_i,iwm_j) = iwm_p_inst(iwm_i,iwm_j) - p_bar
end do
end do

call test_filter_plane_gpu(iwm_u_inst)
call test_filter_plane_gpu(iwm_v_inst)
call test_filter_plane_gpu(iwm_w_inst)
call test_filter_plane_gpu(iwm_p_inst)

!$acc parallel loop collapse(2) default(present) async(1)                     &
!$acc private(ip, im, jp, jm, Luux, Luvx, Luvy, Lvvy, Lux, Lvy)
do iwm_j = 1, ny
do iwm_i = 1, nx
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx) =                                    &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                                  &
        * (1._rprec-iwm_tR(iwm_i,iwm_j))                                      &
        + iwm_u_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) =                                    &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                                  &
        * (1._rprec-iwm_tR(iwm_i,iwm_j))                                      &
        + iwm_v_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
    iwm_flt_p(iwm_i,iwm_j) = iwm_flt_p(iwm_i,iwm_j)                           &
        * (1._rprec-iwm_tR(iwm_i,iwm_j))                                      &
        + iwm_p_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)

    iwm_unsdy(iwm_i,iwm_j,iwm_dirx) =                                         &
        (iwm_inte(iwm_i,iwm_j,iwm_Lu)-iwm_inte_m(iwm_i,iwm_j,iwm_Lu))         &
        / iwm_dt_l
    iwm_unsdy(iwm_i,iwm_j,iwm_diry) =                                         &
        (iwm_inte(iwm_i,iwm_j,iwm_Lv)-iwm_inte_m(iwm_i,iwm_j,iwm_Lv))         &
        / iwm_dt_l

    ip = iwm_i + 1
    im = iwm_i - 1
    jp = iwm_j + 1
    jm = iwm_j - 1
    if (ip > nx) ip = 1
    if (im < 1) im = nx
    if (jp > ny) jp = 1
    if (jm < 1) jm = ny

    Luux = (iwm_inte(ip,iwm_j,iwm_Luu) - iwm_inte(im,iwm_j,iwm_Luu))          &
        / dx / 2._rprec
    Luvy = (iwm_inte(iwm_i,jp,iwm_Luv) - iwm_inte(iwm_i,jm,iwm_Luv))          &
        / dy / 2._rprec
    Luvx = (iwm_inte(ip,iwm_j,iwm_Luv) - iwm_inte(im,iwm_j,iwm_Luv))          &
        / dx / 2._rprec
    Lvvy = (iwm_inte(iwm_i,jp,iwm_Lvv) - iwm_inte(iwm_i,jm,iwm_Lvv))          &
        / dy / 2._rprec
    Lux = (iwm_inte(ip,iwm_j,iwm_Lu) - iwm_inte(im,iwm_j,iwm_Lu))             &
        / dx / 2._rprec
    Lvy = (iwm_inte(iwm_i,jp,iwm_Lv) - iwm_inte(iwm_i,jm,iwm_Lv))             &
        / dy / 2._rprec

    iwm_conv(iwm_i,iwm_j,iwm_dirx) = Luux + Luvy                              &
        - iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx)*(Lux+Lvy)
    iwm_conv(iwm_i,iwm_j,iwm_diry) = Luvx + Lvvy                              &
        - iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry)*(Lux+Lvy)

    iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx) =                                       &
        (iwm_flt_p(ip,iwm_j)-iwm_flt_p(im,iwm_j))/dx/2._rprec                 &
        * iwm_Dz(iwm_i,iwm_j) - iwm_Dz(iwm_i,iwm_j)
    iwm_PrsGrad(iwm_i,iwm_j,iwm_diry) =                                       &
        (iwm_flt_p(iwm_i,jp)-iwm_flt_p(iwm_i,jm))/dy/2._rprec                 &
        * iwm_Dz(iwm_i,iwm_j)

    iwm_lhs(iwm_i,iwm_j,iwm_dirx) = -iwm_inte(iwm_i,iwm_j,iwm_Lu)             &
        + iwm_dt_l*(iwm_conv(iwm_i,iwm_j,iwm_dirx)                            &
        + iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)                                   &
        - iwm_diff(iwm_i,iwm_j,iwm_dirx))
    iwm_lhs(iwm_i,iwm_j,iwm_diry) = -iwm_inte(iwm_i,iwm_j,iwm_Lv)             &
        + iwm_dt_l*(iwm_conv(iwm_i,iwm_j,iwm_diry)                            &
        + iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)                                   &
        - iwm_diff(iwm_i,iwm_j,iwm_diry))
end do
end do
nullify(autowrap_i, autowrap_j)
return
#endif

! update the u, v for previous time step
do iwm_i = 1, nx
do iwm_j = 1, ny
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx) =                                   &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
    iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry) =                                   &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
end do
end do

! get the instantaneous field
u_inst = u(:,:,1)
! let us do not worry about the half cell displacement in v
v_inst = v(:,:,1)
! w is quadrantic near the wall
w_inst = w(:,:,2)*0.25_rprec

! the real pressure is needed, this step is CODE SPECIFIC!
p_inst = p(:,:,1)
do iwm_i = 1, nx
do iwm_j = 1, ny
    p_inst(iwm_i,iwm_j)= p_inst(iwm_i,iwm_j)                                   &
        -0.5_rprec*(u_inst(iwm_i,iwm_j)*u_inst(iwm_i,iwm_j)                   &
        + v_inst(iwm_i,iwm_j)*v_inst(iwm_i,iwm_j)                              &
        + w_inst(iwm_i,iwm_j)*w_inst(iwm_i,iwm_j))
end do
end do

! obtain the pressure fluctuations
p_bar = 0._rprec
do iwm_i = 1, nx
do iwm_j = 1, ny
    p_bar = p_bar+p_inst(iwm_i,iwm_j)
end do
end do
p_bar=p_bar/nx/ny

do iwm_i = 1, nx
do iwm_j = 1, ny
    p_inst(iwm_i,iwm_j) = p_inst(iwm_i,iwm_j)-p_bar
end do
end do

! all the data enters must be filtered (see Anderson & Meneveau 2011 JFM)
call test_filter ( u_inst )
call test_filter ( v_inst )
call test_filter ( w_inst )
call test_filter ( p_inst )

!temporal filtering
do iwm_i=1,nx
do iwm_j=1,ny
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx) =                                     &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)*(1._rprec-iwm_tR(iwm_i,iwm_j))    &
        + u_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
    iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) =                                     &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)*(1._rprec-iwm_tR(iwm_i,iwm_j))    &
        + v_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
    iwm_flt_p(iwm_i,iwm_j) = iwm_flt_p(iwm_i,iwm_j)                            &
        * (1._rprec-iwm_tR(iwm_i,iwm_j))                                       &
        + p_inst(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j)
end do
end do

! calculate LHS, calculation of the integrals is done from the last time step
! in the subroutine iwm_calc_wallstress, so is iwm_diff
do iwm_i = 1, nx
do iwm_j = 1, ny
    ! the unsteady term
    iwm_unsdy(iwm_i,iwm_j,iwm_dirx) =                                          &
        (iwm_inte(iwm_i,iwm_j,iwm_Lu)-iwm_inte_m(iwm_i,iwm_j,iwm_Lu))/iwm_dt
    iwm_unsdy(iwm_i,iwm_j,iwm_diry) =                                          &
        (iwm_inte(iwm_i,iwm_j,iwm_Lv)-iwm_inte_m(iwm_i,iwm_j,iwm_Lv))/iwm_dt

    ! the convective term
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luu)
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luu)
    Luux = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Luv)
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Luv)
    Luvy = (phip-phim)/dy/2._rprec
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Luv)
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Luv)
    Luvx = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lvv)
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lvv)
    Lvvy = (phip-phim)/dy/2._rprec
    phip = iwm_inte(autowrap_i(iwm_i+1),iwm_j,iwm_Lu )
    phim = iwm_inte(autowrap_i(iwm_i-1),iwm_j,iwm_Lu )
    Lux = (phip-phim)/dx/2._rprec
    phip = iwm_inte(iwm_i,autowrap_j(iwm_j+1),iwm_Lv )
    phim = iwm_inte(iwm_i,autowrap_j(iwm_j-1),iwm_Lv )
    Lvy = (phip-phim)/dy/2._rprec
    iwm_conv(iwm_i,iwm_j,iwm_dirx) = Luux + Luvy                               &
        - iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_dirx)*(Lux+Lvy)
    iwm_conv(iwm_i,iwm_j,iwm_diry) = Luvx + Lvvy                               &
        - iwm_flt_tagvel_m(iwm_i,iwm_j,iWM_diry)*(Lux+Lvy)

    ! the pressure gradient term
    phip = iwm_flt_p(autowrap_i(iwm_i+1),iwm_j)
    phim = iwm_flt_p(autowrap_i(iwm_i-1),iwm_j)
    ! including the mean unit pressure gradient
    iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx) = (phip-phim)/dx/2._rprec                &
        * iwm_Dz(iwm_i,iwm_j) - 1._rprec*iwm_Dz(iwm_i,iwm_j)
    phip = iwm_flt_p(iwm_i,autowrap_j(iwm_j+1))
    phim = iwm_flt_p(iwm_i,autowrap_j(iwm_j-1))
    iwm_PrsGrad(iwm_i,iwm_j,iwm_diry) = (phip-phim)/dy/2._rprec                &
        * iwm_Dz(iwm_i,iwm_j)

    ! the left hand side
    ! this is the integrated momentum equation, except for the Lu term
    iwm_lhs(iwm_i,iwm_j,iwm_dirx) = -iwm_inte(iwm_i,iwm_j,iwm_Lu)              &
        + iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_dirx)                              &
        + iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)                                    &
        - iwm_diff(iwm_i,iwm_j,iwm_dirx) )
    ! this is the integrated momentum equation, except for the Lv term
    iwm_lhs(iwm_i,iwm_j,iwm_diry) = -iwm_inte(iwm_i,iwm_j,iwm_Lv)              &
        + iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_diry)                              &
        + iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)                                    &
        - iwm_diff(iwm_i,iwm_j,iwm_diry) )
end do
end do

nullify(autowrap_i, autowrap_j)

end subroutine iwm_calc_lhs

!*******************************************************************************
subroutine iwm_slv(lhsx,lhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy)
!*******************************************************************************
use types, only : rprec
use param, only : vonk
implicit none

real(rprec), intent(in)  :: lhsx, lhsy, Ux, Uy, Dz, z0, utx, uty
real(rprec), intent(out) :: fx,fy
real(rprec) :: Ax, Ay, Vel, inteLu, inteLv
real(rprec) :: one_minus_z0dz

one_minus_z0dz = 1._rprec - z0/Dz
Ax = (Ux - utx/vonk*log(Dz/z0)) / one_minus_z0dz
Ay = (Uy - uty/vonk*log(Dz/z0)) / one_minus_z0dz
Vel = sqrt(Ux*Ux + Uy*Uy)
inteLu = 0.5_rprec*Dz*Ax*one_minus_z0dz*one_minus_z0dz                       &
    + 1._rprec/vonk*utx*Dz*(z0/Dz - 1._rprec + log(Dz/z0))
inteLv = 0.5_rprec*Dz*Ay*one_minus_z0dz*one_minus_z0dz                       &
    + 1._rprec/vonk*uty*Dz*(z0/Dz - 1._rprec + log(Dz/z0))
fx = inteLu+lhsx
fy = inteLv+lhsy

end subroutine iwm_slv

!*******************************************************************************
subroutine iwm_calc_wallstress
!*******************************************************************************
use types, only : rprec
use param, only : vonk, nx, ny
use test_filtermodule

implicit none

integer :: iwm_i, iwm_j
real(rprec) :: fx, fy, fxp, fyp
real(rprec) :: iwm_tol, iwm_eps
real(rprec) :: a11, a12, a21, a22
real(rprec) :: iwmutxP,iwmutyP
integer :: iter, MaxIter, equil_flag, div_flag
real(rprec) :: equilWMpara,equilutx,equiluty
real(rprec) :: iwmpAx, iwmpAy, iwmputx,iwmputy,iwmpz0,iwmpDz
real(rprec) :: utaup
real(rprec) :: dVelzT, dVelzB, Vel
real(rprec) :: z0_Dz, one_minus_z0_Dz, log_Dz_z0, vonk_sq
real(rprec) :: lhsx_l, lhsy_l, Ux_l, Uy_l, Dz_l, z0_l, utx_l, uty_l
real(rprec) :: Ax_l, Ay_l, inteLu_l, inteLv_l
#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
real(rprec) :: iwm_dt_l
#endif

MaxIter=1500

iwm_tol = 0.000001_rprec
iwm_eps = 0.000000001_rprec

#ifdef ENABLE_CUDA
if (iwm_cuda_enabled()) then
    !$cuf kernel do(2) <<<*,*>>>
    do iwm_j = 1, ny
    do iwm_i = 1, nx
        iwm_utx(iwm_i,iwm_j) = 1._rprec                                      &
            * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx))
        iwm_uty(iwm_i,iwm_j) = 0.1_rprec                                     &
            * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))

        lhsx_l = iwm_lhs(iwm_i,iwm_j,iwm_dirx)
        lhsy_l = iwm_lhs(iwm_i,iwm_j,iwm_diry)
        Ux_l = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
        Uy_l = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
        Dz_l = iwm_Dz(iwm_i,iwm_j)
        z0_l = iwm_z0(iwm_i,iwm_j)
        utx_l = iwm_utx(iwm_i,iwm_j)
        uty_l = iwm_uty(iwm_i,iwm_j)
        z0_Dz = z0_l/Dz_l
        one_minus_z0_Dz = 1._rprec - z0_Dz
        log_Dz_z0 = log(Dz_l/z0_l)
        Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        fx = inteLu_l + lhsx_l
        fy = inteLv_l + lhsy_l

        iter = 0
        equil_flag = 0
        div_flag = 0
        do while (max(abs(fx),abs(fy)) > iwm_tol)
            iwmutxP = iwm_utx(iwm_i,iwm_j) + iwm_eps
            iwmutyP = iwm_uty(iwm_i,iwm_j)
            utx_l = iwmutxP
            uty_l = iwmutyP
            Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
            Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
            inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz   &
                + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
            inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz   &
                + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
            fxp = inteLu_l + lhsx_l
            fyp = inteLv_l + lhsy_l
            a11 = (fxp-fx)/iwm_eps
            a21 = (fyp-fy)/iwm_eps

            iwmutxP = iwm_utx(iwm_i,iwm_j)
            iwmutyP = iwm_uty(iwm_i,iwm_j) + iwm_eps
            utx_l = iwmutxP
            uty_l = iwmutyP
            Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
            Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
            inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz   &
                + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
            inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz   &
                + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
            fxp = inteLu_l + lhsx_l
            fyp = inteLv_l + lhsy_l
            a12 = (fxp-fx)/iwm_eps
            a22 = (fyp-fy)/iwm_eps

            iwm_utx(iwm_i,iwm_j) = iwm_utx(iwm_i,iwm_j)                      &
                - 0.50_rprec*(a22*fx-a12*fy)/(a11*a22-a12*a21)
            iwm_uty(iwm_i,iwm_j) = iwm_uty(iwm_i,iwm_j)                      &
                - 0.50_rprec*(-a21*fx+a11*fy)/(a11*a22-a12*a21)

            utx_l = iwm_utx(iwm_i,iwm_j)
            uty_l = iwm_uty(iwm_i,iwm_j)
            Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
            Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
            inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz   &
                + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
            inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz   &
                + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
            fx = inteLu_l + lhsx_l
            fy = inteLv_l + lhsy_l
            iter = iter + 1
            if (iter > MaxIter) then
                equil_flag = 1
                div_flag = 1
                exit
            end if
        end do

        if (iwm_utx(iwm_i,iwm_j)-1.0 == iwm_utx(iwm_i,iwm_j) .or.             &
            iwm_uty(iwm_i,iwm_j)-1.0 == iwm_uty(iwm_i,iwm_j)) then
            equil_flag = 1
            div_flag = 1
        end if

        equilutx = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                 &
            / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
        equiluty = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                 &
            / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
        if (equil_flag == 1) then
            iwm_utx(iwm_i,iwm_j) = equilutx
            iwm_uty(iwm_i,iwm_j) = equiluty
        end if

        if (equil_flag == 1) then
            iwmpAx = 0._rprec
            iwmpAy = 0._rprec
        else
            iwmpAx = (iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                   &
                - iwm_utx(iwm_i,iwm_j)/vonk                                  &
                * log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))              &
                / (1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j))
            iwmpAy = (iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                   &
                - iwm_uty(iwm_i,iwm_j)/vonk                                  &
                * log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))              &
                / (1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j))
        end if

        if (abs(iwmpAx) > 1._rprec .or. abs(iwmpAy) > 1._rprec) then
            equil_flag = 1
            iwm_utx(iwm_i,iwm_j) = equilutx
            iwm_uty(iwm_i,iwm_j) = equiluty
            iwmpAx = 0._rprec
            iwmpAy = 0._rprec
        end if

        iwm_Ax(iwm_i,iwm_j) = iwmpAx
        iwm_Ay(iwm_i,iwm_j) = iwmpAy

        iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu )
        iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv )
        iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv)
        iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu)
        iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv)

        iwmputx = iwm_utx(iwm_i,iwm_j)
        iwmputy = iwm_uty(iwm_i,iwm_j)
        iwmpDz = iwm_Dz(iwm_i,iwm_j)
        iwmpz0 = iwm_z0(iwm_i,iwm_j)
        z0_Dz = iwmpz0/iwmpDz
        one_minus_z0_Dz = 1._rprec - z0_Dz
        log_Dz_z0 = log(iwmpDz/iwmpz0)
        vonk_sq = vonk*vonk

        iwm_inte(iwm_i,iwm_j,iwm_Lu) = 0.5_rprec*iwmpDz*iwmpAx               &
            * one_minus_z0_Dz*one_minus_z0_Dz                                &
            + 1._rprec/vonk*iwmputx*iwmpDz*(z0_Dz - 1._rprec + log_Dz_z0)
        iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0.5_rprec*iwmpDz*iwmpAy               &
            * one_minus_z0_Dz*one_minus_z0_Dz                                &
            + 1._rprec/vonk*iwmputy*iwmpDz*(z0_Dz - 1._rprec + log_Dz_z0)
        iwm_inte(iwm_i,iwm_j,iwm_Luv) =                                      &
            1._rprec/vonk_sq*iwmputx*iwmputy*iwmpDz                          &
            * (1._rprec - 2*z0_Dz + (1._rprec - log_Dz_z0)                   &
            * (1._rprec - log_Dz_z0))                                        &
            + 1._rprec/3._rprec*iwmpAx*iwmpAy*iwmpDz                         &
            * one_minus_z0_Dz*one_minus_z0_Dz*one_minus_z0_Dz                &
            - 0.25_rprec/vonk*(iwmpAx*iwmputy + iwmpAy*iwmputx)*iwmpDz       &
            * (1._rprec - 4._rprec*z0_Dz + 3._rprec*z0_Dz*z0_Dz              &
            - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)
        iwm_inte(iwm_i,iwm_j,iwm_Luu) =                                      &
            1._rprec/vonk_sq*iwmputx*iwmputx*iwmpDz                          &
            * ((log_Dz_z0 - 1._rprec)*(log_Dz_z0 - 1._rprec)                 &
            - 2._rprec*z0_Dz + 1._rprec)                                     &
            + 1._rprec/3._rprec*iwmpAx*iwmpAx*iwmpDz                         &
            * one_minus_z0_Dz*one_minus_z0_Dz                                &
            - 0.5_rprec/vonk*iwmputx*iwmpAx*iwmpDz                           &
            * (1._rprec - 4._rprec*z0_Dz + 3._rprec*z0_Dz*z0_Dz              &
            - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)
        iwm_inte(iwm_i,iwm_j,iwm_Lvv) =                                      &
            1._rprec/vonk_sq*iwmputy*iwmputy*iwmpDz                          &
            * ((log_Dz_z0 - 1._rprec)*(log_Dz_z0 - 1._rprec)                 &
            - 2._rprec*z0_Dz + 1._rprec)                                     &
            + 1._rprec/3._rprec*iwmpAy*iwmpAy*iwmpDz                         &
            * one_minus_z0_Dz*one_minus_z0_Dz                                &
            - 0.5_rprec/vonk*iwmputy*iwmpAy*iwmpDz                           &
            * (1._rprec - 4._rprec*z0_Dz - 3._rprec*z0_Dz*z0_Dz              &
            - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)

        iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*(iwmpAx+iwmputx/vonk)
        iwm_dudzT(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*(iwmpAy+iwmputy/vonk)
        iwm_dudzB(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*iwmpAx                  &
            + iwmputx/vonk/iwmpz0
        iwm_dudzB(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*iwmpAy                  &
            + iwmputy/vonk/iwmpz0

        Vel = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                      &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                           &
            + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                           &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))
        dVelzT = abs(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel                &
            * iwm_dudzT(iwm_i,iwm_j,iwm_dirx)                                &
            + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Vel                       &
            * iwm_dudzT(iwm_i,iwm_j,iwm_diry))
        dVelzB = sqrt(iwm_dudzB(iwm_i,iwm_j,iwm_dirx)                        &
            * iwm_dudzB(iwm_i,iwm_j,iwm_dirx)                                &
            + iwm_dudzB(iwm_i,iwm_j,iwm_diry)                                &
            * iwm_dudzB(iwm_i,iwm_j,iwm_diry))

        iwm_diff(iwm_i,iwm_j,iwm_dirx) = vonk_sq*iwmpDz*iwmpDz*dVelzT        &
            * iwm_dudzT(iwm_i,iwm_j,iwm_dirx)                                &
            - vonk_sq*iwmpz0*iwmpz0*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
        iwm_diff(iwm_i,iwm_j,iwm_diry) = vonk_sq*iwmpDz*iwmpDz*dVelzT        &
            * iwm_dudzT(iwm_i,iwm_j,iwm_diry)                                &
            - vonk_sq*iwmpz0*iwmpz0*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_diry)

        if (equil_flag == 1) then
            equilWMpara = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)          &
                * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                       &
                + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                       &
                * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))                      &
                * vonk_sq/(log_Dz_z0*log_Dz_z0)
            iwm_tauwx(iwm_i,iwm_j) = equilWMpara                             &
                * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
            iwm_tauwy(iwm_i,iwm_j) = equilWMpara                             &
                * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
        else
            iwm_tauwx(iwm_i,iwm_j) = vonk_sq*iwmpz0*iwmpz0*dVelzB            &
                * iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
            iwm_tauwy(iwm_i,iwm_j) = vonk_sq*iwmpz0*iwmpz0*dVelzB            &
                * iwm_dudzB(iwm_i,iwm_j,iwm_diry)
        end if

        utaup = sqrt(sqrt(iwm_tauwx(iwm_i,iwm_j)*iwm_tauwx(iwm_i,iwm_j)      &
            + iwm_tauwy(iwm_i,iwm_j)*iwm_tauwy(iwm_i,iwm_j)))
        iwm_flt_us(iwm_i,iwm_j) = iwm_flt_us(iwm_i,iwm_j)                    &
            * (1._rprec-iwm_tR(iwm_i,iwm_j)) + utaup*iwm_tR(iwm_i,iwm_j)
        iwm_tR(iwm_i,iwm_j) = iwm_dt/(iwm_Dz(iwm_i,iwm_j)                    &
            / iwm_flt_us(iwm_i,iwm_j) / vonk)
        if (iwm_tR(iwm_i,iwm_j) > 1._rprec) then
            iwm_tR(iwm_i,iwm_j) = 1._rprec
        end if
    end do
    end do
    call iwm_cuda_sync('iwm_calc_wallstress')
    return
end if
#endif

#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
iwm_dt_l = iwm_dt
! The filtered friction velocity update is point-local in (i,j); keep this
! loop on the same async queue as the surrounding IWM kernels.
! The (i,j) update is independent after the filtered-friction indexing fix, so
! the OpenACC path can expose the full wall plane instead of serializing j.
!$acc parallel loop collapse(2) default(present) async(1)                     &
!$acc private(fx, fy, fxp, fyp, a11, a12, a21, a22, iwmutxP, iwmutyP, iter,   &
!$acc         equil_flag, div_flag, equilWMpara, equilutx, equiluty, iwmpAx,  &
!$acc         iwmpAy, iwmputx, iwmputy, iwmpz0, iwmpDz, utaup, dVelzT,        &
!$acc         dVelzB, Vel, z0_Dz, one_minus_z0_Dz, log_Dz_z0, vonk_sq,        &
!$acc         lhsx_l, lhsy_l, Ux_l, Uy_l, Dz_l, z0_l, utx_l, uty_l, Ax_l,     &
!$acc         Ay_l, inteLu_l, inteLv_l)
do iwm_i = 1, nx
do iwm_j = 1, ny
    iwm_utx(iwm_i,iwm_j) = 1._rprec                                          &
        * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx))
    iwm_uty(iwm_i,iwm_j) = 0.1_rprec                                         &
        * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))

    lhsx_l = iwm_lhs(iwm_i,iwm_j,iwm_dirx)
    lhsy_l = iwm_lhs(iwm_i,iwm_j,iwm_diry)
    Ux_l = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
    Uy_l = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
    Dz_l = iwm_Dz(iwm_i,iwm_j)
    z0_l = iwm_z0(iwm_i,iwm_j)
    utx_l = iwm_utx(iwm_i,iwm_j)
    uty_l = iwm_uty(iwm_i,iwm_j)
    z0_Dz = z0_l/Dz_l
    one_minus_z0_Dz = 1._rprec - z0_Dz
    log_Dz_z0 = log(Dz_l/z0_l)
    Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
    Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
    inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz           &
        + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
    inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz           &
        + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
    fx = inteLu_l + lhsx_l
    fy = inteLv_l + lhsy_l

    iter = 0
    equil_flag = 0
    div_flag = 0
    do while (max(abs(fx),abs(fy)) > iwm_tol)
        iwmutxP = iwm_utx(iwm_i,iwm_j) + iwm_eps
        iwmutyP = iwm_uty(iwm_i,iwm_j)
        utx_l = iwmutxP
        uty_l = iwmutyP
        Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        fxp = inteLu_l + lhsx_l
        fyp = inteLv_l + lhsy_l
        a11 = (fxp-fx)/iwm_eps
        a21 = (fyp-fy)/iwm_eps

        iwmutxP = iwm_utx(iwm_i,iwm_j)
        iwmutyP = iwm_uty(iwm_i,iwm_j) + iwm_eps
        utx_l = iwmutxP
        uty_l = iwmutyP
        Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        fxp = inteLu_l + lhsx_l
        fyp = inteLv_l + lhsy_l
        a12 = (fxp-fx)/iwm_eps
        a22 = (fyp-fy)/iwm_eps

        iwm_utx(iwm_i,iwm_j) = iwm_utx(iwm_i,iwm_j)                          &
            - 0.50_rprec*(a22*fx-a12*fy)/(a11*a22-a12*a21)
        iwm_uty(iwm_i,iwm_j) = iwm_uty(iwm_i,iwm_j)                          &
            - 0.50_rprec*(-a21*fx+a11*fy)/(a11*a22-a12*a21)

        utx_l = iwm_utx(iwm_i,iwm_j)
        uty_l = iwm_uty(iwm_i,iwm_j)
        Ax_l = (Ux_l - utx_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        Ay_l = (Uy_l - uty_l/vonk*log_Dz_z0) / one_minus_z0_Dz
        inteLu_l = 0.5_rprec*Dz_l*Ax_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*utx_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        inteLv_l = 0.5_rprec*Dz_l*Ay_l*one_minus_z0_Dz*one_minus_z0_Dz       &
            + 1._rprec/vonk*uty_l*Dz_l*(z0_Dz - 1._rprec + log_Dz_z0)
        fx = inteLu_l + lhsx_l
        fy = inteLv_l + lhsy_l
        iter = iter + 1
        if (iter > MaxIter) then
            equil_flag = 1
            div_flag = 1
            exit
        end if
    end do

    if (iwm_utx(iwm_i,iwm_j)-1.0 == iwm_utx(iwm_i,iwm_j) .or.                 &
        iwm_uty(iwm_i,iwm_j)-1.0 == iwm_uty(iwm_i,iwm_j)) then
        equil_flag = 1
        div_flag = 1
    end if

    equilutx = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                     &
        / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
    equiluty = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                     &
        / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
    if (equil_flag == 1) then
        iwm_utx(iwm_i,iwm_j) = equilutx
        iwm_uty(iwm_i,iwm_j) = equiluty
    end if

    if (equil_flag == 1) then
        iwmpAx = 0._rprec
        iwmpAy = 0._rprec
    else
        iwmpAx = (iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                       &
            - iwm_utx(iwm_i,iwm_j)/vonk                                      &
            * log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))                  &
            / (1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j))
        iwmpAy = (iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                       &
            - iwm_uty(iwm_i,iwm_j)/vonk                                      &
            * log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j)))                  &
            / (1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j))
    end if

    if (abs(iwmpAx) > 1._rprec .or. abs(iwmpAy) > 1._rprec) then
        equil_flag = 1
        iwm_utx(iwm_i,iwm_j) = equilutx
        iwm_uty(iwm_i,iwm_j) = equiluty
        iwmpAx = 0._rprec
        iwmpAy = 0._rprec
    end if

    iwm_Ax(iwm_i,iwm_j) = iwmpAx
    iwm_Ay(iwm_i,iwm_j) = iwmpAy

    iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu )
    iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv )
    iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv)
    iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu)
    iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv)

    iwmputx = iwm_utx(iwm_i,iwm_j)
    iwmputy = iwm_uty(iwm_i,iwm_j)
    iwmpDz = iwm_Dz(iwm_i,iwm_j)
    iwmpz0 = iwm_z0(iwm_i,iwm_j)
    z0_Dz = iwmpz0/iwmpDz
    one_minus_z0_Dz = 1._rprec - z0_Dz
    log_Dz_z0 = log(iwmpDz/iwmpz0)
    vonk_sq = vonk*vonk

    iwm_inte(iwm_i,iwm_j,iwm_Lu) = 0.5_rprec*iwmpDz*iwmpAx                   &
        * one_minus_z0_Dz*one_minus_z0_Dz                                    &
        + 1._rprec/vonk*iwmputx*iwmpDz*(z0_Dz - 1._rprec + log_Dz_z0)
    iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0.5_rprec*iwmpDz*iwmpAy                   &
        * one_minus_z0_Dz*one_minus_z0_Dz                                    &
        + 1._rprec/vonk*iwmputy*iwmpDz*(z0_Dz - 1._rprec + log_Dz_z0)
    iwm_inte(iwm_i,iwm_j,iwm_Luv) =                                          &
        1._rprec/vonk_sq*iwmputx*iwmputy*iwmpDz                              &
        * (1._rprec - 2*z0_Dz + (1._rprec - log_Dz_z0)                       &
        * (1._rprec - log_Dz_z0))                                            &
        + 1._rprec/3._rprec*iwmpAx*iwmpAy*iwmpDz                             &
        * one_minus_z0_Dz*one_minus_z0_Dz*one_minus_z0_Dz                    &
        - 0.25_rprec/vonk*(iwmpAx*iwmputy + iwmpAy*iwmputx)*iwmpDz           &
        * (1._rprec - 4._rprec*z0_Dz + 3._rprec*z0_Dz*z0_Dz                  &
        - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)
    iwm_inte(iwm_i,iwm_j,iwm_Luu) =                                          &
        1._rprec/vonk_sq*iwmputx*iwmputx*iwmpDz                              &
        * ((log_Dz_z0 - 1._rprec)*(log_Dz_z0 - 1._rprec)                     &
        - 2._rprec*z0_Dz + 1._rprec)                                         &
        + 1._rprec/3._rprec*iwmpAx*iwmpAx*iwmpDz                             &
        * one_minus_z0_Dz*one_minus_z0_Dz                                    &
        - 0.5_rprec/vonk*iwmputx*iwmpAx*iwmpDz                               &
        * (1._rprec - 4._rprec*z0_Dz + 3._rprec*z0_Dz*z0_Dz                  &
        - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)
    iwm_inte(iwm_i,iwm_j,iwm_Lvv) =                                          &
        1._rprec/vonk_sq*iwmputy*iwmputy*iwmpDz                              &
        * ((log_Dz_z0 - 1._rprec)*(log_Dz_z0 - 1._rprec)                     &
        - 2._rprec*z0_Dz + 1._rprec)                                         &
        + 1._rprec/3._rprec*iwmpAy*iwmpAy*iwmpDz                             &
        * one_minus_z0_Dz*one_minus_z0_Dz                                    &
        - 0.5_rprec/vonk*iwmputy*iwmpAy*iwmpDz                               &
        * (1._rprec - 4._rprec*z0_Dz - 3._rprec*z0_Dz*z0_Dz                  &
        - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)

    iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*(iwmpAx+iwmputx/vonk)
    iwm_dudzT(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*(iwmpAy+iwmputy/vonk)
    iwm_dudzB(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*iwmpAx                      &
        + iwmputx/vonk/iwmpz0
    iwm_dudzB(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*iwmpAy                      &
        + iwmputy/vonk/iwmpz0

    Vel = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                          &
        * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                               &
        + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                               &
        * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))
    dVelzT = abs(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel                    &
        * iwm_dudzT(iwm_i,iwm_j,iwm_dirx)                                    &
        + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Vel                           &
        * iwm_dudzT(iwm_i,iwm_j,iwm_diry))
    dVelzB = sqrt(iwm_dudzB(iwm_i,iwm_j,iwm_dirx)                            &
        * iwm_dudzB(iwm_i,iwm_j,iwm_dirx)                                    &
        + iwm_dudzB(iwm_i,iwm_j,iwm_diry)                                    &
        * iwm_dudzB(iwm_i,iwm_j,iwm_diry))

    iwm_diff(iwm_i,iwm_j,iwm_dirx) = vonk_sq*iwmpDz*iwmpDz*dVelzT            &
        * iwm_dudzT(iwm_i,iwm_j,iwm_dirx)                                    &
        - vonk_sq*iwmpz0*iwmpz0*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
    iwm_diff(iwm_i,iwm_j,iwm_diry) = vonk_sq*iwmpDz*iwmpDz*dVelzT            &
        * iwm_dudzT(iwm_i,iwm_j,iwm_diry)                                    &
        - vonk_sq*iwmpz0*iwmpz0*dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_diry)

    if (equil_flag == 1) then
        equilWMpara = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)              &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                           &
            + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                           &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))                          &
            * vonk_sq/(log_Dz_z0*log_Dz_z0)
        iwm_tauwx(iwm_i,iwm_j) = equilWMpara                                 &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
        iwm_tauwy(iwm_i,iwm_j) = equilWMpara                                 &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
    else
        iwm_tauwx(iwm_i,iwm_j) = vonk_sq*iwmpz0*iwmpz0*dVelzB                &
            * iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
        iwm_tauwy(iwm_i,iwm_j) = vonk_sq*iwmpz0*iwmpz0*dVelzB                &
            * iwm_dudzB(iwm_i,iwm_j,iwm_diry)
    end if

    utaup = sqrt(sqrt(iwm_tauwx(iwm_i,iwm_j)*iwm_tauwx(iwm_i,iwm_j)          &
        + iwm_tauwy(iwm_i,iwm_j)*iwm_tauwy(iwm_i,iwm_j)))
    iwm_flt_us(iwm_i,iwm_j) = iwm_flt_us(iwm_i,iwm_j)                        &
        * (1._rprec-iwm_tR(iwm_i,iwm_j)) + utaup*iwm_tR(iwm_i,iwm_j)
    iwm_tR(iwm_i,iwm_j) = iwm_dt_l/(iwm_Dz(iwm_i,iwm_j)                      &
        / iwm_flt_us(iwm_i,iwm_j) / vonk)
    if (iwm_tR(iwm_i,iwm_j) > 1._rprec) then
        iwm_tR(iwm_i,iwm_j) = 1._rprec
    end if
end do
end do
return
#endif

do iwm_i=1,nx
do iwm_j=1,ny

    ! use Newton method to solve the system
    iwm_utx(iwm_i,iwm_j) = 1._rprec                                            &
        * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx))
    iwm_uty(iwm_i,iwm_j) = 0.1_rprec                                           &
        * sign(1._rprec,iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))

    call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx), iwm_lhs(iwm_i,iwm_j,iwm_diry), &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),                                  &
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),                                  &
        iwm_Dz(iwm_i,iwm_j), iwm_z0(iwm_i,iwm_j),                              &
        iwm_utx(iwm_i,iwm_j), iwm_uty(iwm_i,iwm_j),fx,fy )

    iter = 0
    equil_flag = 0
    div_flag = 0
    do while (max(abs(fx),abs(fy))>iwm_tol)
        iwmutxP=iwm_utx(iwm_i,iwm_j)+iWM_eps
        iwmutyP=iwm_uty(iwm_i,iwm_j)
        call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),                            &
            iwm_lhs(iwm_i,iwm_j,iwm_diry),iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),&
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),iwm_Dz(iwm_i,iwm_j),          &
            iwm_z0(iwm_i,iwm_j), iwmutxP, iwmutyP, fxp, fyp )
        a11 = (fxp-fx)/iwm_eps
        a21 = (fyp-fy)/iwm_eps
        iwmutxP = iwm_utx(iwm_i,iwm_j)
        iwmutyP = iwm_uty(iwm_i,iwm_j)+iwm_eps
        call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),                            &
            iwm_lhs(iwm_i,iwm_j,iwm_diry),                                     &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),                              &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),                              &
            iwm_Dz(iwm_i,iwm_j), iwm_z0(iwm_i,iwm_j),                          &
            iwmutxP, iwmutyP, fxp, fyp)
        a12 = (fxp-fx)/iwm_eps
        a22 = (fyp-fy)/iwm_eps
        iwm_utx(iwm_i,iwm_j) = iwm_utx(iwm_i,iwm_j)                            &
            - 0.50*( a22*fx-a12*fy)/(a11*a22-a12*a21)
        iwm_uty(iwm_i,iwm_j) = iwm_uty(iwm_i,iwm_j)                            &
            - 0.50*(-a21*fx+a11*fy)/(a11*a22-a12*a21)
        call iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),                            &
            iwm_lhs(iwm_i,iwm_j,iwm_diry),                                     &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),                              &
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),                              &
            iwm_Dz(iwm_i,iwm_j), iwm_z0(iwm_i,iwm_j),                          &
            iwm_utx(iwm_i,iwm_j), iwm_uty(iwm_i,iwm_j), fx, fy)
        iter = iter+1

        ! maximum iteration reached
        if (iter>MaxIter) then
            equil_flag = 1
            div_flag = 1;
            exit
        end if
    end do

    ! infinity check
    if (iwm_utx(iwm_i,iwm_j)-1.0==iwm_utx(iwm_i,iwm_j) .or.                    &
        iwm_uty(iwm_i,iwm_j)-1.0==iwm_uty(iwm_i,iwm_j)) then
        equil_flag=1
        div_flag  =1
    end if

    ! calculate equilibrium us for equil_flag=1 use
    equilutx = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                       &
    / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
    equiluty = vonk*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                       &
        / log(iwm_Dz(iwm_i,iwm_j)/iwm_z0(iwm_i,iwm_j))
    if (equil_flag==1) then
        iwm_utx(iwm_i,iwm_j) = equilutx
        iwm_uty(iwm_i,iwm_j) = equiluty
    end if

    !calculate Ax, Ay
    if(equil_flag==1)then
        iwmpAx = 0._rprec
        iwmpAy=0._rprec
    else
        ! eq. D2 in Yang et al. 2015
        iwmpAx = ( iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                        &
            -iwm_utx(iwm_i,iwm_j)/vonk*log(iwm_Dz(iwm_i,iwm_j)                 &
            /iwm_z0(iwm_i,iwm_j)))                                             &
            / ((1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j)))
        iwmpAy = ( iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                        &
            -iwm_uty(iwm_i,iwm_j)/vonk*log(iwm_Dz(iwm_i,iwm_j)                 &
            /iwm_z0(iwm_i,iwm_j)))                                             &
            /((1._rprec-iwm_z0(iwm_i,iwm_j)/iwm_Dz(iwm_i,iwm_j)))
    end if

    ! check for excessive linear term correction
    ! after first 100 step this check is rarely invoked
    if (abs(iwmpAx)>1._rprec .or. abs(iwmpAy)>1._rprec) then
        equil_flag = 1
        iwm_utx(iwm_i,iwm_j) = equilutx
        iwm_uty(iwm_i,iwm_j) = equiluty
        iwmpAx = 0._rprec
        iwmpAy = 0._rprec
    end if

    ! store the linear correction
    iwm_Ax(iwm_i,iwm_j) = iwmpAx
    iwm_Ay(iwm_i,iwm_j) = iwmpAy

    ! update integral for last time step
    iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu )
    iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv )
    iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv)
    iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu)
    iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv)

    ! Cache point-local wall-model values used repeatedly in the nonlinear solve.
    iwmputx = iwm_utx(iwm_i,iwm_j)
    iwmputy = iwm_uty(iwm_i,iwm_j)
    iwmpDz  = iwm_Dz (iwm_i,iwm_j)
    iwmpz0  = iwm_z0 (iwm_i,iwm_j)
    z0_Dz = iwmpz0/iwmpDz
    one_minus_z0_Dz = 1._rprec - z0_Dz
    log_Dz_z0 = log(iwmpDz/iwmpz0)
    vonk_sq = vonk*vonk

    ! calculate the needed integrals

    ! Eq. D7 in Yang et al. 2015
    iwm_inte(iwm_i,iwm_j,iwm_Lu) = 0.5_rprec*iwmpDz*iwmpAx                    &
        * one_minus_z0_Dz*one_minus_z0_Dz                                      &
        + 1._rprec/vonk*iwmputx*iwmpDz*(z0_Dz - 1._rprec + log_Dz_z0)
    iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0.5_rprec*iwmpDz*iwmpAy                    &
        * one_minus_z0_Dz*one_minus_z0_Dz                                      &
        + 1._rprec/vonk*iwmputy*iwmpDz*(z0_Dz - 1._rprec + log_Dz_z0)

    ! Eq. D8 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1._rprec/vonk_sq*iwmputx*iwmputy*iwmpDz   &
        * (1._rprec - 2*z0_Dz + (1._rprec - log_Dz_z0)                        &
        * (1._rprec - log_Dz_z0))                                              &
        + 1._rprec/3._rprec*iwmpAx*iwmpAy*iwmpDz*one_minus_z0_Dz              &
        * one_minus_z0_Dz*one_minus_z0_Dz                                      &
        - 0.25_rprec/vonk*(iwmpAx*iwmputy + iwmpAy*iwmputx)*iwmpDz            &
        * (1._rprec - 4._rprec*z0_Dz + 3._rprec*z0_Dz*z0_Dz                   &
        - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)
    iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1._rprec/vonk_sq*iwmputx*iwmputx*iwmpDz   &
        * ((log_Dz_z0 - 1._rprec)*(log_Dz_z0 - 1._rprec)                      &
        - 2._rprec*z0_Dz + 1._rprec)                                           &
        + 1._rprec/3._rprec*iwmpAx*iwmpAx*iwmpDz*one_minus_z0_Dz              &
        * one_minus_z0_Dz*one_minus_z0_Dz                                      &
        - 0.5_rprec/vonk*iwmputx*iwmpAx*iwmpDz                                &
        * (1._rprec - 4._rprec*z0_Dz + 3._rprec*z0_Dz*z0_Dz                   &
        - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)

    ! Eq. D9 in Yang et al 2015
    iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1._rprec/vonk_sq*iwmputy*iwmputy*iwmpDz   &
        * ((log_Dz_z0 - 1._rprec)*(log_Dz_z0 - 1._rprec)                      &
        - 2._rprec*z0_Dz + 1._rprec)                                           &
        + 1._rprec/3._rprec*iwmpAy*iwmpAy*iwmpDz*one_minus_z0_Dz              &
        * one_minus_z0_Dz*one_minus_z0_Dz                                      &
        - 0.5_rprec/vonk*iwmputy*iwmpAy*iwmpDz                                &
        * (1._rprec - 4._rprec*z0_Dz - 3._rprec*z0_Dz*z0_Dz                   &
        - 2._rprec*log_Dz_z0 + 4._rprec*z0_Dz*log_Dz_z0)

    ! calculate top and bottom derivatives
    ! Eq. D5 (a)
    iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*(iwmpAx+iwmputx/vonk)
    iwm_dudzT(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*(iwmpAy+iwmputy/vonk)
    ! Eq. D5 (b)
    iwm_dudzB(iwm_i,iwm_j,iwm_dirx) = 1.0/iwmpDz*iwmpAx+iwmputx/vonk/iwmpz0
    iwm_dudzB(iwm_i,iwm_j,iwm_diry) = 1.0/iwmpDz*iwmpAy+iwmputy/vonk/iwmpz0

    ! calculte the turbulent diffusion term
    !total velocity
    Vel = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                           &
        * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                                 &
        + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                                 &
        * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))
    ! Eq. D6
    dVelzT=abs(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel                        &
        *iwm_dudzT(iwm_i,iwm_j,iwm_dirx)+iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)  &
        /Vel*iwm_dudzT(iwm_i,iwm_j,iwm_diry))
    dvelzB = sqrt(iwm_dudzB(iwm_i,iwm_j,iwm_dirx)                              &
        * iwm_dudzB(iwm_i,iwm_j,iwm_dirx)                                      &
        + iwm_dudzB(iwm_i,iwm_j,iwm_diry)                                      &
        * iwm_dudzB(iwm_i,iwm_j,iwm_diry))

    ! Eq. D4, the eddy viscosity is nu_T=(vonk*y)^2*dudy, hence the formula
    iwm_diff(iwm_i,iwm_j,iwm_dirx) = vonk_sq*iwmpDz*iwmpDz*dVelzT             &
        *iwm_dudzT(iwm_i,iwm_j,iwm_dirx) - vonk_sq*iwmpz0*iwmpz0               &
        *dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
    iwm_diff(iwm_i,iwm_j,iwm_diry) = vonk_sq*iwmpDz*iwmpDz*dVelzT             &
        *iwm_dudzT(iwm_i,iwm_j,iwm_diry) - vonk_sq*iwmpz0*iwmpz0               &
        *dVelzB*iwm_dudzB(iwm_i,iwm_j,iwm_diry)

    ! calculate the wall stress
    if (equil_flag==1) then
        equilWMpara = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)               &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)                             &
            + iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)                             &
            * iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry))                            &
            * vonk_sq/(log_Dz_z0*log_Dz_z0)
        iwm_tauwx(iwm_i,iwm_j) = equilWMpara                                   &
            *iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)
        iwm_tauwy(iwm_i,iwm_j) = equilWMpara                                   &
            *iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)
    else
        ! Eq. D4
        iwm_tauwx(iwm_i,iwm_j) = vonk_sq*iwmpz0*iwmpz0*dVelzB                 &
            *iwm_dudzB(iwm_i,iwm_j,iwm_dirx)
        iwm_tauwy(iwm_i,iwm_j) = vonk_sq*iwmpz0*iwmpz0*dVelzB                 &
            *iwm_dudzB(iwm_i,iwm_j,iwm_diry)
    end if

    ! calculate the friciton velocity
    ! definition of friction velocity
    utaup = sqrt(sqrt(iwm_tauwx(iwm_i,iwm_j)*iwm_tauwx(iwm_i,iwm_j)            &
        + iwm_tauwy(iwm_i,iwm_j)*iwm_tauwy(iwm_i,iwm_j)))
    ! the filtered friction velocity used for filtering time scale
    iwm_flt_us(iwm_i,iwm_j) = iwm_flt_us(iwm_i,iwm_j)                          &
        *(1._rprec-iwm_tR(iwm_i,iwm_j))+utaup*iwm_tR(iwm_i,iwm_j)

    ! update the filtering time scale
    ! Eq. 26
    iwm_tR(iwm_i,iwm_j) = iwm_dt/(iwm_Dz(iwm_i,iwm_j)                          &
        /iwm_flt_us(iwm_i,iwm_j)/vonk)
    ! filtering time scale can only be larger than the time step,
    ! if not, then just use the instantaneous flow field to do the model
    if (iwm_tR(iwm_i,iwm_j)>1._rprec) then
        iwm_tR(iwm_i,iwm_j) = 1._rprec
    end if

end do
end do

end subroutine iwm_calc_wallstress


!*******************************************************************************
subroutine iwm_monitor
!*******************************************************************************
!
! This subroutine is to monitor the parameters at one point, do not call this
! subroutine if you are not interested in how the model works
!
use param, only : nx,ny,jt_total
implicit none

integer :: iwm_i,iwm_j,dmpPrd,fid
character*50 :: fname

dmpPrd = iwm_ntime_skip
iwm_i = int(nx/2._rprec)
iwm_j = int(ny/2._rprec)

if( mod(jt_total,dmpPrd)==0)then
#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
!$acc wait(1)
!$acc update self(iwm_flt_tagvel(iwm_i:iwm_i,iwm_j:iwm_j,1:iwm_DN),             &
!$acc&            iwm_utx(iwm_i:iwm_i,iwm_j:iwm_j),                            &
!$acc&            iwm_uty(iwm_i:iwm_i,iwm_j:iwm_j),                            &
!$acc&            iwm_Ax(iwm_i:iwm_i,iwm_j:iwm_j),                              &
!$acc&            iwm_Ay(iwm_i:iwm_i,iwm_j:iwm_j),                              &
!$acc&            iwm_tR(iwm_i:iwm_i,iwm_j:iwm_j))
#endif
endif

write(fname,'(A,i5.5,A)') 'iwm_track.dat'
open(newunit=fid, file=fname, status='unknown', form='formatted',              &
    position='append')
if( mod(jt_total,dmpPrd)==0)then
write(fid,*) iwm_flt_tagvel(iwm_i,iwm_j,:), iwm_utx(iwm_i,iwm_j),              &
    iwm_uty(iwm_i,iwm_j),  iwm_Ax(iwm_i,iwm_j),                                &
    iwm_Ay(iwm_i,iwm_j), iwm_tR(iwm_i,iwm_j)
end if
close(fid)

end subroutine iwm_monitor


!*******************************************************************************
subroutine iwm_checkPoint()
!*******************************************************************************
!
! This subroutine checkpoints the integral wall model. It is called after making
! sure lbc_mom=3
!
implicit none

integer :: fid

#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
!$acc wait(1)
!$acc update self(iwm_utx, iwm_uty, iwm_tauwx, iwm_tauwy)
!$acc update self(iwm_flt_tagvel, iwm_flt_tagvel_m, iwm_flt_p)
!$acc update self(iwm_inte, iwm_inte_m, iwm_unsdy, iwm_conv)
!$acc update self(iwm_PrsGrad, iwm_diff, iwm_LHS)
!$acc update self(iwm_dudzT, iwm_dudzB, iwm_flt_us, iwm_tR)
!$acc update self(iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay)
#endif

open(newunit=fid, file='iwm_checkPoint.dat', status='unknown',                 &
    form='unformatted', position='rewind')
write(fid) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:),         &
    iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN),              &
    iwm_flt_p(:,:), iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),          &
    iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN),&
    iwm_diff(:,:,1:iwm_DN), iwm_LHS(:,:,1:iwm_DN), iwm_dudzT(:,:,1:iwm_DN),    &
    iwm_dudzB(:,:,1:iwm_DN), iwm_flt_us(:,:), iwm_tR(:,:), iwm_Dz(:,:),        &
    iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), iwm_dt
close(fid)

end subroutine iwm_checkPoint

!*******************************************************************************
subroutine iwm_read_checkPoint()
!*******************************************************************************
!
! Read checkpoint data for the integral wall model.  Call only after confirming
! lbc_mom=3.
!
implicit none

integer :: fid

open(newunit=fid, file='iwm_checkPoint.dat', status='unknown',                 &
    form='unformatted', position='rewind')
read(fid) iwm_utx(:,:), iwm_uty(:,:), iwm_tauwx(:,:), iwm_tauwy(:,:),          &
    iwm_flt_tagvel(:,:,1:iwm_DN), iwm_flt_tagvel_m(:,:,1:iwm_DN),              &
    iwm_flt_p(:,:), iwm_inte(:,:,1:iwm_LN), iWM_inte_m(:,:,1:iwm_LN),          &
    iwm_unsdy(:,:,1:iwm_DN), iwm_conv(:,:,1:iwm_DN), iwm_PrsGrad(:,:,1:iwm_DN),&
    iwm_diff(:,:,1:iwm_DN), iwm_LHS(:,:,1:iwm_DN), iwm_dudzT(:,:,1:iwm_DN),    &
    iwm_dudzB(:,:,1:iwm_DN), iwm_flt_us(:,:), iwm_tR(:,:), iwm_Dz(:,:),        &
    iwm_z0(:,:), iwm_Ax(:,:), iwm_Ay(:,:), iwm_dt
close(fid)

#if defined(PPSGS_GPU) && !defined(ENABLE_CUDA)
!$acc update device(iwm_utx, iwm_uty, iwm_tauwx, iwm_tauwy)
!$acc update device(iwm_flt_tagvel, iwm_flt_tagvel_m, iwm_flt_p)
!$acc update device(iwm_inte, iwm_inte_m, iwm_unsdy, iwm_conv)
!$acc update device(iwm_PrsGrad, iwm_diff, iwm_LHS)
!$acc update device(iwm_dudzT, iwm_dudzB, iwm_flt_us, iwm_tR)
!$acc update device(iwm_Dz, iwm_z0, iwm_Ax, iwm_Ay)
#endif

end subroutine iwm_read_checkPoint

end module iwmles
