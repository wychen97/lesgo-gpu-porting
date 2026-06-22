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
subroutine scaledep_dynamic(Cs_1D)
!*******************************************************************************
!
! Subroutine uses the scale dependent dynamic model to calculate the Smagorinsky
! coefficient Cs_1D and |S|. This is done layer-by-layer to save memory.
!
use types, only : rprec
use param, only : ld, nx, ny, nz, coord
use sim_param, only : u, v, w
use test_filtermodule
use sgs_stag_util, only : rtnewt
use sgs_param, only : S11, S12, S13, S22, S23, S33, delta, S, ee_now,          &
    u_bar, v_bar, w_bar, L11, L12, L13, L22, L23, L33, S_bar,                  &
    S11_bar, S12_bar, S13_bar, S22_bar, S23_bar, S33_bar,                      &
    S_S11_bar, S_S12_bar, S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar,          &
    u_hat, v_hat, w_hat,                                                       &
    S_hat, S11_hat, S12_hat, S13_hat, S22_hat, S23_hat, S33_hat,               &
    S_S11_hat, S_S12_hat, S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat
#ifdef ENABLE_CUDA
use cudafor
#endif
implicit none

integer :: jx, jy, jz
real(rprec), dimension(nz), intent (inout) :: Cs_1D
#ifdef ENABLE_CUDA
real(rprec), managed, save, allocatable, target, dimension(:,:) :: Q11, Q12, Q13, Q22,  &
    Q23, Q33
real(rprec), managed, save, dimension(:), allocatable :: beta
#else
real(rprec), save, allocatable, target, dimension(:,:) :: Q11, Q12, Q13, Q22,  &
    Q23, Q33
real(rprec), save, dimension(:), allocatable :: beta
#endif
real(rprec), pointer, dimension(:,:) :: M11, M12, M13, M22, M23, M33
logical, save :: arrays_allocated = .false.
real(rprec) :: const
real(rprec), dimension(0:5) :: A
real(rprec) :: a1, b1, c1, d1, e1, a2, b2, c2, d2, e2
real(rprec) :: lm_sum, mm_sum, cs_val, beta_val
real(rprec), dimension(ld,ny) :: LM,MM
#ifdef ENABLE_CUDA
logical :: scaledep_dynamic_cuda_enabled
logical :: use_scaledep_cuda

use_scaledep_cuda = scaledep_dynamic_cuda_enabled()
#endif

! Allocate arrays
if( .not. arrays_allocated ) then
    allocate ( Q11(ld,ny), Q12(ld,ny), Q13(ld,ny), &
        Q22(ld,ny), Q23(ld,ny), Q33(ld,ny) )
    allocate ( beta(nz) )
    arrays_allocated = .true.
endif

! Associate pointers
M11 => Q11
M12 => Q12
M13 => Q13
M22 => Q22
M23 => Q23
M33 => Q33

do jz = 1, nz
#ifdef ENABLE_CUDA
    if (use_scaledep_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            if ((coord == 0) .and. (jz == 1)) then
                L11(jx,jy) = u(jx,jy,1)*u(jx,jy,1)
                L12(jx,jy) = u(jx,jy,1)*v(jx,jy,1)
                L13(jx,jy) = u(jx,jy,1)*0.25_rprec*w(jx,jy,2)
                L22(jx,jy) = v(jx,jy,1)*v(jx,jy,1)
                L23(jx,jy) = v(jx,jy,jz)*0.25_rprec*w(jx,jy,2)
                L33(jx,jy) = (0.25_rprec*w(jx,jy,2))**2
                u_bar(jx,jy) = u(jx,jy,1)
                v_bar(jx,jy) = v(jx,jy,1)
                w_bar(jx,jy) = 0.25_rprec*w(jx,jy,2)
            else
                u_bar(jx,jy) = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
                v_bar(jx,jy) = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
                w_bar(jx,jy) = w(jx,jy,jz)
                L11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
                L12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
                L13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
                L22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
                L23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
                L33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
            end if
            u_hat(jx,jy) = u_bar(jx,jy)
            v_hat(jx,jy) = v_bar(jx,jy)
            w_hat(jx,jy) = w_bar(jx,jy)
        end do
        end do
        call scaledep_dynamic_cuda_sync('velocity and L setup')

        call test_filter_3 ( u_bar, v_bar, w_bar )
        call test_filter_6 ( L11, L12, L13, L22, L23, L33 )
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            L11(jx,jy) = L11(jx,jy) - u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = L12(jx,jy) - u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = L13(jx,jy) - u_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = L22(jx,jy) - v_bar(jx,jy)*v_bar(jx,jy)
            L23(jx,jy) = L23(jx,jy) - v_bar(jx,jy)*w_bar(jx,jy)
            L33(jx,jy) = L33(jx,jy) - w_bar(jx,jy)*w_bar(jx,jy)
            Q11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
            Q12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
            Q13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
            Q22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
            Q23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
            Q33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
        end do
        end do
        call scaledep_dynamic_cuda_sync('filtered Lij correction and Q setup')

        call test_test_filter_3 ( u_hat, v_hat, w_hat )
        call test_test_filter_6 ( Q11, Q12, Q13, Q22, Q23, Q33 )
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            Q11(jx,jy) = Q11(jx,jy) - u_hat(jx,jy)*u_hat(jx,jy)
            Q12(jx,jy) = Q12(jx,jy) - u_hat(jx,jy)*v_hat(jx,jy)
            Q13(jx,jy) = Q13(jx,jy) - u_hat(jx,jy)*w_hat(jx,jy)
            Q22(jx,jy) = Q22(jx,jy) - v_hat(jx,jy)*v_hat(jx,jy)
            Q23(jx,jy) = Q23(jx,jy) - v_hat(jx,jy)*w_hat(jx,jy)
            Q33(jx,jy) = Q33(jx,jy) - w_hat(jx,jy)*w_hat(jx,jy)
        end do
        end do
        call scaledep_dynamic_cuda_sync('filtered Qij correction')

!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +      &
                S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +                 &
                S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
            S11_bar(jx,jy) = S11(jx,jy,jz)
            S12_bar(jx,jy) = S12(jx,jy,jz)
            S13_bar(jx,jy) = S13(jx,jy,jz)
            S22_bar(jx,jy) = S22(jx,jy,jz)
            S23_bar(jx,jy) = S23(jx,jy,jz)
            S33_bar(jx,jy) = S33(jx,jy,jz)
            S11_hat(jx,jy) = S11_bar(jx,jy)
            S12_hat(jx,jy) = S12_bar(jx,jy)
            S13_hat(jx,jy) = S13_bar(jx,jy)
            S22_hat(jx,jy) = S22_bar(jx,jy)
            S23_hat(jx,jy) = S23_bar(jx,jy)
            S33_hat(jx,jy) = S33_bar(jx,jy)
            S_S11_bar(jx,jy) = S(jx,jy)*S11(jx,jy,jz)
            S_S12_bar(jx,jy) = S(jx,jy)*S12(jx,jy,jz)
            S_S13_bar(jx,jy) = S(jx,jy)*S13(jx,jy,jz)
            S_S22_bar(jx,jy) = S(jx,jy)*S22(jx,jy,jz)
            S_S23_bar(jx,jy) = S(jx,jy)*S23(jx,jy,jz)
            S_S33_bar(jx,jy) = S(jx,jy)*S33(jx,jy,jz)
            S_S11_hat(jx,jy) = S_S11_bar(jx,jy)
            S_S12_hat(jx,jy) = S_S12_bar(jx,jy)
            S_S13_hat(jx,jy) = S_S13_bar(jx,jy)
            S_S22_hat(jx,jy) = S_S22_bar(jx,jy)
            S_S23_hat(jx,jy) = S_S23_bar(jx,jy)
            S_S33_hat(jx,jy) = S_S33_bar(jx,jy)
        end do
        end do
        call scaledep_dynamic_cuda_sync('S setup')

        call test_filter_6 ( S11_bar, S12_bar, S13_bar, S22_bar, S23_bar,      &
            S33_bar )
        call test_test_filter_6 ( S11_hat, S12_hat, S13_hat, S22_hat,          &
            S23_hat, S33_hat )
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            S_bar(jx,jy) = sqrt(2._rprec*(S11_bar(jx,jy)**2 +                 &
                S22_bar(jx,jy)**2 + S33_bar(jx,jy)**2 + 2._rprec*(            &
                S12_bar(jx,jy)**2 + S13_bar(jx,jy)**2 + S23_bar(jx,jy)**2)))
            S_hat(jx,jy) = sqrt(2._rprec*(S11_hat(jx,jy)**2 +                 &
                S22_hat(jx,jy)**2 + S33_hat(jx,jy)**2 + 2._rprec*(            &
                S12_hat(jx,jy)**2 + S13_hat(jx,jy)**2 + S23_hat(jx,jy)**2)))
        end do
        end do
        call scaledep_dynamic_cuda_sync('filtered S magnitudes')

        call test_filter_6 ( S_S11_bar, S_S12_bar, S_S13_bar, S_S22_bar,       &
            S_S23_bar, S_S33_bar )
        call test_test_filter_6 ( S_S11_hat, S_S12_hat, S_S13_hat, S_S22_hat,  &
            S_S23_hat, S_S33_hat )

        a1 = 0._rprec
        b1 = 0._rprec
        c1 = 0._rprec
        d1 = 0._rprec
        e1 = 0._rprec
        a2 = 0._rprec
        b2 = 0._rprec
        c2 = 0._rprec
        d2 = 0._rprec
        e2 = 0._rprec
!$cuf kernel do(2) <<<*,*>>> reduction(+:a1,b1,c1,d1,e1,a2,b2,c2,d2,e2)
        do jy = 1, ny
        do jx = 1, ld
            a1 = a1 + S_bar(jx,jy)*(S11_bar(jx,jy)*L11(jx,jy) +                &
                S22_bar(jx,jy)*L22(jx,jy) + S33_bar(jx,jy)*L33(jx,jy) +        &
                2._rprec*(S12_bar(jx,jy)*L12(jx,jy) +                         &
                S13_bar(jx,jy)*L13(jx,jy) + S23_bar(jx,jy)*L23(jx,jy)))
            b1 = b1 + S_S11_bar(jx,jy)*L11(jx,jy) +                            &
                S_S22_bar(jx,jy)*L22(jx,jy) + S_S33_bar(jx,jy)*L33(jx,jy) +    &
                2._rprec*(S_S12_bar(jx,jy)*L12(jx,jy) +                       &
                S_S13_bar(jx,jy)*L13(jx,jy) + S_S23_bar(jx,jy)*L23(jx,jy))
            c1 = c1 + S_S11_bar(jx,jy)**2 + S_S22_bar(jx,jy)**2 +              &
                S_S33_bar(jx,jy)**2 + 2._rprec*(S_S12_bar(jx,jy)**2 +          &
                S_S13_bar(jx,jy)**2 + S_S23_bar(jx,jy)**2)
            d1 = d1 + 0.5_rprec*S_bar(jx,jy)*S_bar(jx,jy)*                    &
                S_bar(jx,jy)*S_bar(jx,jy)
            e1 = e1 + S_bar(jx,jy)*(S11_bar(jx,jy)*S_S11_bar(jx,jy) +          &
                S22_bar(jx,jy)*S_S22_bar(jx,jy) +                              &
                S33_bar(jx,jy)*S_S33_bar(jx,jy) + 2._rprec*(                  &
                S12_bar(jx,jy)*S_S12_bar(jx,jy) +                              &
                S13_bar(jx,jy)*S_S13_bar(jx,jy) +                              &
                S23_bar(jx,jy)*S_S23_bar(jx,jy)))
            a2 = a2 + S_hat(jx,jy)*(S11_hat(jx,jy)*Q11(jx,jy) +                &
                S22_hat(jx,jy)*Q22(jx,jy) + S33_hat(jx,jy)*Q33(jx,jy) +        &
                2._rprec*(S12_hat(jx,jy)*Q12(jx,jy) +                         &
                S13_hat(jx,jy)*Q13(jx,jy) + S23_hat(jx,jy)*Q23(jx,jy)))
            b2 = b2 + S_S11_hat(jx,jy)*Q11(jx,jy) +                            &
                S_S22_hat(jx,jy)*Q22(jx,jy) + S_S33_hat(jx,jy)*Q33(jx,jy) +    &
                2._rprec*(S_S12_hat(jx,jy)*Q12(jx,jy) +                       &
                S_S13_hat(jx,jy)*Q13(jx,jy) + S_S23_hat(jx,jy)*Q23(jx,jy))
            c2 = c2 + S_S11_hat(jx,jy)**2 + S_S22_hat(jx,jy)**2 +              &
                S_S33_hat(jx,jy)**2 + 2._rprec*(S_S12_hat(jx,jy)**2 +          &
                S_S13_hat(jx,jy)**2 + S_S23_hat(jx,jy)**2)
            d2 = d2 + 0.5_rprec*S_hat(jx,jy)*S_hat(jx,jy)*                    &
                S_hat(jx,jy)*S_hat(jx,jy)
            e2 = e2 + S_hat(jx,jy)*(S11_hat(jx,jy)*S_S11_hat(jx,jy) +          &
                S22_hat(jx,jy)*S_S22_hat(jx,jy) +                              &
                S33_hat(jx,jy)*S_S33_hat(jx,jy) + 2._rprec*(                  &
                S12_hat(jx,jy)*S_S12_hat(jx,jy) +                              &
                S13_hat(jx,jy)*S_S13_hat(jx,jy) +                              &
                S23_hat(jx,jy)*S_S23_hat(jx,jy)))
        end do
        end do
        call scaledep_dynamic_cuda_sync('A-coefficient reductions')

        a1 = -2._rprec*(delta**2)*4._rprec*a1/(nx*ny)
        b1 = -2._rprec*(delta**2)*b1/(nx*ny)
        c1 = (2._rprec*delta**2)**2*c1/(nx*ny)
        d1 = (2._rprec*delta**2)**2*16._rprec*d1/(nx*ny)
        e1 = 2._rprec*(2._rprec*delta**2)**2*4._rprec*e1/(nx*ny)
        a2 = -2._rprec*(delta**2)*16._rprec*a2/(nx*ny)
        b2 = -2._rprec*(delta**2)*b2/(nx*ny)
        c2 = (2._rprec*delta**2)**2*c2/(nx*ny)
        d2 = (2._rprec*delta**2)**2*256._rprec*d2/(nx*ny)
        e2 = 2._rprec*(2._rprec*delta**2)**2*16._rprec*e2/(nx*ny)

        A(0) = b2*c1 - b1*c2
        A(1) = a1*c2 - b2*e1
        A(2) = b2*d1 + b1*e2 - a2*c1
        A(3) = a2*e1 - a1*e2
        A(4) = -a2*d1 - b1*d2
        A(5) = a1*d2

        beta_val = rtnewt(A,jz)
        beta(jz) = beta_val

        const = 2._rprec*delta**2
        lm_sum = 0._rprec
        mm_sum = 0._rprec
!$cuf kernel do(2) <<<*,*>>> reduction(+:lm_sum,mm_sum)
        do jy = 1, ny
        do jx = 1, ld
            Q11(jx,jy) = const*(S_S11_bar(jx,jy) -                              &
                4._rprec*beta_val*S_bar(jx,jy)*S11_bar(jx,jy))
            Q12(jx,jy) = const*(S_S12_bar(jx,jy) -                              &
                4._rprec*beta_val*S_bar(jx,jy)*S12_bar(jx,jy))
            Q13(jx,jy) = const*(S_S13_bar(jx,jy) -                              &
                4._rprec*beta_val*S_bar(jx,jy)*S13_bar(jx,jy))
            Q22(jx,jy) = const*(S_S22_bar(jx,jy) -                              &
                4._rprec*beta_val*S_bar(jx,jy)*S22_bar(jx,jy))
            Q23(jx,jy) = const*(S_S23_bar(jx,jy) -                              &
                4._rprec*beta_val*S_bar(jx,jy)*S23_bar(jx,jy))
            Q33(jx,jy) = const*(S_S33_bar(jx,jy) -                              &
                4._rprec*beta_val*S_bar(jx,jy)*S33_bar(jx,jy))
            lm_sum = lm_sum + L11(jx,jy)*Q11(jx,jy) + L22(jx,jy)*Q22(jx,jy) +    &
                L33(jx,jy)*Q33(jx,jy) + 2._rprec*(L12(jx,jy)*Q12(jx,jy) +        &
                L13(jx,jy)*Q13(jx,jy) + L23(jx,jy)*Q23(jx,jy))
            mm_sum = mm_sum + Q11(jx,jy)**2 + Q22(jx,jy)**2 + Q33(jx,jy)**2 +    &
                2._rprec*(Q12(jx,jy)**2 + Q13(jx,jy)**2 + Q23(jx,jy)**2)
        end do
        end do
        call scaledep_dynamic_cuda_sync('Mij reductions')

        cs_val = max(0._rprec, real(lm_sum/mm_sum,rprec))
        Cs_1D(jz) = cs_val
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            ee_now(jx,jy,jz) = L11(jx,jy)**2 + L22(jx,jy)**2 +                  &
                L33(jx,jy)**2 + 2._rprec*(L12(jx,jy)**2 + L13(jx,jy)**2 +        &
                L23(jx,jy)**2) - 2._rprec*(L11(jx,jy)*Q11(jx,jy) +              &
                L22(jx,jy)*Q22(jx,jy) + L33(jx,jy)*Q33(jx,jy) +                 &
                2._rprec*(L12(jx,jy)*Q12(jx,jy) + L13(jx,jy)*Q13(jx,jy) +        &
                L23(jx,jy)*Q23(jx,jy)))*cs_val + (Q11(jx,jy)**2 +               &
                Q22(jx,jy)**2 + Q33(jx,jy)**2 + 2._rprec*(Q12(jx,jy)**2 +       &
                Q13(jx,jy)**2 + Q23(jx,jy)**2))*cs_val**2
        end do
        end do
        call scaledep_dynamic_cuda_sync('ee_now')
    else
#endif
    ! using L_ij as temp storage here
    if ( (coord == 0) .and. (jz == 1) ) then
        !!! watch the 0.25's: w = c*z**z near wall, so get 0.25
        ! put on uvp-nodes
        L11(:,:) = u(:,:,1)*u(:,:,1)  ! uv-node
        L12(:,:) = u(:,:,1)*v(:,:,1)  ! uv-node
        L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2)  ! assume parabolic near wall
        L22(:,:) = v(:,:,1)*v(:,:,1)  ! uv-node
        L23(:,:) = v(:,:,jz)*0.25_rprec*w(:,:,2)  ! uv-node
        L33(:,:) = (0.25_rprec*w(:,:,2))**2  ! uv-node
        u_bar(:,:) = u(:,:,1)
        v_bar(:,:) = v(:,:,1)
        w_bar(:,:) = 0.25_rprec*w(:,:,2)
    else  ! w-nodes
        L11(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*                        &
                0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        L12(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*                        &
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        L13(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*w(:,:,jz)
        L22(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*                        &
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        L23(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*w(:,:,jz)
        L33(:,:) = w(:,:,jz)*w(:,:,jz)
        u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        w_bar(:,:) = w(:,:,jz)
    end if
    u_hat = u_bar
    v_hat = v_bar
    w_hat = w_bar

    ! Filter first term and add the second term to get the final value
    ! in-place filtering
    call test_filter ( u_bar )
    call test_filter ( v_bar )
    call test_filter ( w_bar )
    call test_filter ( L11 )
    L11 = L11 - u_bar*u_bar
    call test_filter ( L12 )
    L12 = L12 - u_bar*v_bar
    call test_filter ( L13 )
    L13 = L13 - u_bar*w_bar
    call test_filter ( L22 )
    L22 = L22 - v_bar*v_bar
    call test_filter ( L23 )
    L23 = L23 - v_bar*w_bar
    call test_filter ( L33 )
    L33 = L33 - w_bar*w_bar

    Q11 = u_bar*u_bar
    Q12 = u_bar*v_bar
    Q13 = u_bar*w_bar
    Q22 = v_bar*v_bar
    Q23 = v_bar*w_bar
    Q33 = w_bar*w_bar

    call test_test_filter ( u_hat )
    call test_test_filter ( v_hat )
    call test_test_filter ( w_hat )
    call test_test_filter ( Q11 )
    Q11 = Q11 - u_hat*u_hat
    call test_test_filter ( Q12 )
    Q12 = Q12 - u_hat*v_hat
    call test_test_filter ( Q13 )
    Q13 = Q13 - u_hat*w_hat
    call test_test_filter ( Q22 )
    Q22 = Q22 - v_hat*v_hat
    call test_test_filter ( Q23 )
    Q23 = Q23 - v_hat*w_hat
    call test_test_filter ( Q33 )
    Q33 = Q33 - w_hat*w_hat

    ! calculate |S|
    S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +&
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

    ! already on w-nodes
    S11_bar(:,:) = S11(:,:,jz)
    S12_bar(:,:) = S12(:,:,jz)
    S13_bar(:,:) = S13(:,:,jz)
    S22_bar(:,:) = S22(:,:,jz)
    S23_bar(:,:) = S23(:,:,jz)
    S33_bar(:,:) = S33(:,:,jz)

    S11_hat = S11_bar
    S12_hat = S12_bar
    S13_hat = S13_bar
    S22_hat = S22_bar
    S23_hat = S23_bar
    S33_hat = S33_bar

    call test_filter ( S11_bar )
    call test_filter ( S12_bar )
    call test_filter ( S13_bar )
    call test_filter ( S22_bar )
    call test_filter ( S23_bar )
    call test_filter ( S33_bar )

    call test_test_filter ( S11_hat )
    call test_test_filter ( S12_hat )
    call test_test_filter ( S13_hat )
    call test_test_filter ( S22_hat )
    call test_test_filter ( S23_hat )
    call test_test_filter ( S33_hat )

    S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

    S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
        2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

    S_S11_bar = S*S11(:,:,jz)
    S_S12_bar = S*S12(:,:,jz)
    S_S13_bar = S*S13(:,:,jz)
    S_S22_bar = S*S22(:,:,jz)
    S_S23_bar = S*S23(:,:,jz)
    S_S33_bar = S*S33(:,:,jz)

    S_S11_hat = S_S11_bar
    S_S12_hat = S_S12_bar
    S_S13_hat = S_S13_bar
    S_S22_hat = S_S22_bar
    S_S23_hat = S_S23_bar
    S_S33_hat = S_S33_bar

    call test_filter ( S_S11_bar )
    call test_filter ( S_S12_bar )
    call test_filter ( S_S13_bar )
    call test_filter ( S_S22_bar )
    call test_filter ( S_S23_bar )
    call test_filter ( S_S33_bar )

    call test_test_filter ( S_S11_hat )
    call test_test_filter ( S_S12_hat )
    call test_test_filter ( S_S13_hat )
    call test_test_filter ( S_S22_hat )
    call test_test_filter ( S_S23_hat )
    call test_test_filter ( S_S33_hat )

    ! note: check that the Nyquist guys are zero!
    ! the 1./(nx*ny) is not really neccessary, but in practice it does
    ! affect the results.
    a1 = -2._rprec*(delta**2)*4._rprec*sum( S_bar*(S11_bar*L11 + S22_bar*L22 + &
        S33_bar*L33 + 2._rprec*( S12_bar*L12 + S13_bar*L13 + S23_bar*L23)))/   &
        (nx*ny)
    b1 = -2._rprec*(delta**2)*sum( S_S11_bar*L11 + S_S22_bar*L22 +             &
        S_S33_bar*L33 + 2._rprec*( S_S12_bar*L12 + S_S13_bar*L13 +             &
        S_S23_bar*L23))/(nx*ny)
    c1 = (2._rprec*delta**2)**2 * sum( S_S11_bar**2 + S_S22_bar**2 +           &
        S_S33_bar**2 + &
        2._rprec*( S_S12_bar**2 + S_S13_bar**2 + S_S23_bar**2))/(nx*ny)
    d1 = (2._rprec*delta**2)**2 * 16._rprec*sum( 0.5_rprec*S_bar**4 )/(nx*ny)
    e1 = 2._rprec*(2._rprec*delta**2)**2*4._rprec*sum( S_bar*( S11_bar*        &
        S_S11_bar +S22_bar*S_S22_bar + S33_bar*S_S33_bar + 2._rprec*(          &
        S12_bar*S_S12_bar + S13_bar*S_S13_bar + S23_bar*S_S23_bar)))/          &
        (nx*ny)

    a2 = -2._rprec*(delta**2)*16._rprec*sum( S_hat*(S11_hat*Q11 + S22_hat*Q22 +&
        S33_hat*Q33 + 2._rprec*( S12_hat*Q12 + S13_hat*Q13 + S23_hat*Q23)))/&
        (nx*ny)
    b2 = -2._rprec*(delta**2)*sum( S_S11_hat*Q11 + S_S22_hat*Q22 +&
        S_S33_hat*Q33 + 2._rprec*( S_S12_hat*Q12 + S_S13_hat*Q13 + &
        S_S23_hat*Q23))/(nx*ny)
    c2 = (2._rprec*delta**2)**2 * sum( S_S11_hat**2 + S_S22_hat**2 +&
        S_S33_hat**2 + &
        2._rprec*( S_S12_hat**2 + S_S13_hat**2 + S_S23_hat**2))/(nx*ny)
    d2 = (2._rprec*delta**2)**2 *256._rprec*sum( 0.5_rprec*S_hat**4 )/(nx*ny)
    e2 = 2._rprec*(2._rprec*delta**2)**2*16._rprec*sum( S_hat*( S11_hat*       &
        S_S11_hat + S22_hat*S_S22_hat + S33_hat*S_S33_hat + 2.*(               &
        S12_hat*S_S12_hat + S13_hat*S_S13_hat + S23_hat*S_S23_hat)))/          &
        (nx*ny)

    A(0) = b2*c1 - b1*c2
    A(1) = a1*c2 - b2*e1
    A(2) = b2*d1 + b1*e2 - a2*c1
    A(3) = a2*e1 - a1*e2
    A(4) = -a2*d1 - b1*d2
    A(5) = a1*d2

    beta(jz) = rtnewt(A,jz)

    ! now put beta back into M_ij: using Q_ij as storage
    const = 2._rprec*delta**2
    M11 = const*(S_S11_bar - 4._rprec*beta(jz)*S_bar*S11_bar)
    M12 = const*(S_S12_bar - 4._rprec*beta(jz)*S_bar*S12_bar)
    M13 = const*(S_S13_bar - 4._rprec*beta(jz)*S_bar*S13_bar)
    M22 = const*(S_S22_bar - 4._rprec*beta(jz)*S_bar*S22_bar)
    M23 = const*(S_S23_bar - 4._rprec*beta(jz)*S_bar*S23_bar)
    M33 = const*(S_S33_bar - 4._rprec*beta(jz)*S_bar*S33_bar)

    Cs_1D(jz) = sum(L11*M11 + L22*M22 + L33*M33 + 2._rprec*(L12*M12 +          &
        L13*M13 + L23*M23)) /                                                  &
        sum(M11**2 + M22**2 + M33**2 + 2._rprec*(M12**2 + M13**2 + M23**2))
    Cs_1D(jz) = max(0._rprec, real(Cs_1D(jz),rprec))

    ! Calculate ee_now (the current value of eij*eij)
    LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
    MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)
    ee_now(:,:,jz) = L11**2+L22**2+L33**2+2._rprec*(L12**2+L13**2+L23**2)      &
                    -2._rprec*LM*Cs_1D(jz) + MM*Cs_1D(jz)**2
#ifdef ENABLE_CUDA
    end if
#endif
end do

! Nullify pointers
nullify( M11, M12, M13, M22, M23, M33 )

end subroutine scaledep_dynamic

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function scaledep_dynamic_cuda_enabled()
!*******************************************************************************
implicit none

scaledep_dynamic_cuda_enabled = .true.

end function scaledep_dynamic_cuda_enabled

!*******************************************************************************
subroutine scaledep_dynamic_cuda_sync(where)
!*******************************************************************************
use cudafor
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'scaledep_dynamic CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'scaledep_dynamic CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine scaledep_dynamic_cuda_sync
#endif
