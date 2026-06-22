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
module turbine_indicator
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lh
#ifdef ENABLE_CUDA
use cudafor
use cufft
#endif

private
public :: turb_ind_func_t

! Indicator function calculator
type turb_ind_func_t
    real(rprec), dimension(:), allocatable :: r
    real(rprec), dimension(:), allocatable :: R2, R2_t
    real(rprec) :: sqrt6overdelta1, ell, delta1, delta2, thk
    real(rprec) :: M
contains
    procedure, public :: init
    procedure, public :: val
end type turb_ind_func_t

contains

!*******************************************************************************
logical function turbine_indicator_cuda_enabled()
!*******************************************************************************
implicit none

turbine_indicator_cuda_enabled = .true.

end function turbine_indicator_cuda_enabled

!*******************************************************************************
subroutine val(this, r, x, Rval, Rval_t)
!*******************************************************************************
use functions, only : linear_interp
implicit none
class(turb_ind_func_t), intent(in) :: this
real(rprec), intent(in) :: r, x
real(rprec), intent(out) :: Rval, Rval_t
real(rprec) :: R1, R2, R2_t, rr, xx

! Make dimensionless to lookup value
rr = r / this%ell
xx = x / this%ell

! Calculate as a dimensionless value, then return as dimensional
R2 = linear_interp(this%r, this%R2, rr)
R2_t = linear_interp(this%r, this%R2_t, rr)
R1 = 0.5_rprec / this%thk * ( erf(this%sqrt6overdelta1*(xx + 0.5*this%thk))    &
    - erf(this%sqrt6overdelta1*(xx - 0.5*this%thk)) )
Rval = R1 * R2 / this%ell**3
Rval_t = R1 * R2_t / this%ell**3

end subroutine val

!*******************************************************************************
subroutine init(this, delta1, delta2, thk, dia)
!*******************************************************************************
use param, only : write_endian, path, pi
use functions, only : bilinear_interp
implicit none
include'fftw3.f'

class(turb_ind_func_t), intent(inout) :: this
real(rprec), intent(in) :: delta1, delta2, thk, dia

real(rprec) :: L, d, dr
integer :: i, j, N
real(rprec), dimension(:), allocatable :: yz
real(rprec), dimension(:,:), allocatable :: g, fx, ft, h, t

integer*8 plan
complex(rprec), dimension(:,:), allocatable :: ghat, fxhat, fthat, hhat, that
#ifdef ENABLE_CUDA
real(rprec), managed, dimension(:), allocatable :: rtab, r2tab, r2ttab
integer :: plan_fw, plan_bk, istat
integer :: ci, cj
real(rprec) :: norm, msum
#endif

! Units for all values
this%ell = 0.5_rprec * dia
this%delta1 = delta1 / this%ell
this%delta2 = delta2 / this%ell
this%thk = thk / this%ell
this%sqrt6overdelta1 = sqrt(6._rprec) / this%delta1

! Size of radial domain
L = 1 + 10*this%delta2/sqrt(12._rprec)
N = 2 * 1024
d = 2 * L / N

#ifdef ENABLE_CUDA
if (turbine_indicator_cuda_enabled()) then
    allocate(yz(N))
    !$cuf kernel do(1) <<<*, *>>>
    do i = 1, N
        yz(i) = -L + d * (real(i, rprec) - 0.5_rprec)
    end do

    allocate(h(N, N))
    allocate(t(N, N))
    allocate(g(N, N))
    !$cuf kernel do(2) <<<*, *>>>
    do j = 1, N
    do i = 1, N
        g(i,j) = 6._rprec/(pi*this%delta2*this%delta2)                       &
            * exp(-6._rprec*(yz(i)*yz(i)+yz(j)*yz(j))                        &
            /(this%delta2*this%delta2))
        if (sqrt(yz(i)*yz(i) + yz(j)*yz(j)) < 1._rprec) then
            h(i,j) = 1.0_rprec / pi
            t(i,j) = yz(j) / (pi * max(yz(i)*yz(i) + yz(j)*yz(j),             &
                0.000001_rprec))
        else
            h(i,j) = 0._rprec
            t(i,j) = 0._rprec
        end if
    end do
    end do
    istat = cudaDeviceSynchronize()
    if (istat /= cudaSuccess) then
        print *, 'turbine_indicator init kernel failed: ', istat
        stop
    end if

    allocate(ghat(N/2+1, N))
    allocate(hhat(N/2+1, N))
    allocate(fxhat(N/2+1, N))
    allocate(fthat(N/2+1, N))
    allocate(that(N/2+1, N))
    allocate(fx(N, N))
    allocate(ft(N, N))

    istat = cufftPlan2d(plan_fw, N, N, CUFFT_D2Z)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator forward cuFFT plan failed: ', istat
        stop
    end if
    istat = cufftPlan2d(plan_bk, N, N, CUFFT_Z2D)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator inverse cuFFT plan failed: ', istat
        stop
    end if

    istat = cufftExecD2Z(plan_fw, g, ghat)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator g cuFFT failed: ', istat
        stop
    end if
    istat = cufftExecD2Z(plan_fw, h, hhat)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator h cuFFT failed: ', istat
        stop
    end if
    istat = cufftExecD2Z(plan_fw, t, that)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator t cuFFT failed: ', istat
        stop
    end if

    !$cuf kernel do(2) <<<*, *>>>
    do cj = 1, N
    do ci = 1, N/2+1
        fxhat(ci,cj) = ghat(ci,cj) * hhat(ci,cj) * d*d
        fthat(ci,cj) = ghat(ci,cj) * that(ci,cj) * d*d
    end do
    end do

    istat = cufftExecZ2D(plan_bk, fxhat, fx)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator fx inverse cuFFT failed: ', istat
        stop
    end if
    istat = cufftExecZ2D(plan_bk, fthat, ft)
    if (istat /= CUFFT_SUCCESS) then
        print *, 'turbine_indicator ft inverse cuFFT failed: ', istat
        stop
    end if

    allocate(rtab(N/2))
    allocate(r2tab(N/2))
    allocate(r2ttab(N/2))
    !$cuf kernel do(1) <<<*, *>>>
    do i = 1, N/2
        rtab(i) = yz(i) + L
        r2tab(i) = fx(1,i) / real(N*N, rprec)
        r2ttab(i) = ft(1,i) / real(N*N, rprec)
    end do
    istat = cudaDeviceSynchronize()
    if (istat /= cudaSuccess) then
        print *, 'turbine_indicator lookup kernel failed: ', istat
        stop
    end if

    dr = rtab(2) - rtab(1)
    norm = 0._rprec
    !$cuf kernel do(1) <<<*, *>>> reduction(+:norm)
    do i = 1, N/2
        norm = norm + 2._rprec*pi*rtab(i)*r2tab(i)*dr
    end do

    !$cuf kernel do(1) <<<*, *>>>
    do i = 1, N/2
        r2ttab(i) = r2ttab(i) / norm
        r2tab(i) = r2tab(i) / norm
    end do

    msum = 0._rprec
    !$cuf kernel do(1) <<<*, *>>> reduction(+:msum)
    do i = 1, N/2
        msum = msum + 2._rprec*pi*rtab(i)*r2tab(i)*r2tab(i)*dr
    end do
    this%M = msum * pi

    allocate(this%r(N/2))
    allocate(this%R2(N/2))
    allocate(this%R2_t(N/2))
    this%r = rtab
    this%R2 = r2tab
    this%R2_t = r2ttab

    istat = cufftDestroy(plan_fw)
    istat = cufftDestroy(plan_bk)
    deallocate(yz, h, t, g, ghat, hhat, fxhat, fthat, that, fx, ft)
    deallocate(rtab, r2tab, r2ttab)
    return
end if
#endif

! Create yz grid
allocate(yz(N))
do i = 1, N
    yz(i) = -L + d * (i - 0.5_rprec)
end do

! Create Circle and Gaussian
allocate(h(N, N))
allocate(t(N, N))
allocate(g(N, N))
do i = 1, N
    do j = 1, N
        g(i,j) = 6._rprec/(pi*this%delta2**2)                                  &
            * exp(-6*(yz(i)**2+yz(j)**2)/this%delta2**2)
        if (sqrt(yz(i)**2 + yz(j)**2) < 1) then
            h(i,j) = 1.0_rprec / pi
            t(i,j) = 1.0_rprec / pi * yz(j) / max(yz(i)**2 + yz(j)**2,0.000001_rprec)
        else
            h(i,j) = 0._rprec
            t(i,j) = 0._rprec
        end if
    end do
end do

! Do the convolution f = g*h in fourier space
allocate(ghat(N/2+1, N))
allocate(hhat(N/2+1, N))
allocate(fxhat(N/2+1, N))
allocate(fthat(N/2+1, N))
allocate(that(N/2+1, N))
allocate(fx(N, N))
allocate(ft(N, N))

call dfftw_plan_dft_r2c_2d(plan, N, N, g, ghat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, g, ghat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, h, hhat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, h, hhat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, t, that, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, t, that)
call dfftw_destroy_plan(plan)

fxhat = ghat * hhat * d**2!
fthat = ghat * that * d**2!

call dfftw_plan_dft_c2r_2d(plan, N, N, fxhat, fx, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fxhat, fx)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_2d(plan, N, N, fthat, ft, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fthat, ft)
call dfftw_destroy_plan(plan)

! Place into the lookup table
allocate(this%r(N/2))
allocate(this%R2(N/2))
allocate(this%R2_t(N/2))
this%r = yz(1:N/2) + L
this%R2 = fx(1,1:N/2) / N**2
this%R2_t = ft(1,1:N/2) / N**2

! Normalize to integrate to unity
dr = this%r(2) - this%r(1)
this%R2_t = this%R2_t/sum(2*pi*this%r*this%R2*dr)
this%R2 = this%R2/sum(2*pi*this%r*this%R2*dr)

! Calculate correction factor
this%M = sum(2*pi*this%r*this%R2**2*dr)*pi

end subroutine init

end module turbine_indicator
